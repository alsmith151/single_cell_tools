use rayon::prelude::*;
// use rust_htslib::{bam, bam::record::Aux, bam::Header, bam::Read};
use anyhow::{Result};

use noodles::{bam, bgzf, sam};
use sam::alignment::record::data::field::value::Value;
use std::fs::File;

use std::path::PathBuf;
use std::str::from_utf8;
use std::str::{FromStr};
use std::{collections::HashMap};
use tempdir::TempDir;



pub struct SampleBarcodes {
    sample: String,
    barcodes: Vec<String>,
}

impl SampleBarcodes {
    pub fn new(sample: String, barcodes: Vec<String>) -> SampleBarcodes {
        SampleBarcodes {
            sample: sample,
            barcodes: barcodes,
        }
    }
    
}


pub struct SampleBarcodeCollection {
    samples: Vec<String>,
    barcodes: HashMap<String, Vec<String>>,
}

impl SampleBarcodeCollection {
    pub fn new(
        samples: Vec<String>,
        barcodes: HashMap<String, Vec<String>>,
    ) -> SampleBarcodeCollection {
        SampleBarcodeCollection {
            samples: samples,
            barcodes: barcodes,
        }
    }

    pub fn from_sample_barcodes(
        sample_barcodes: Vec<SampleBarcodes>,
    ) -> Result<SampleBarcodeCollection> {
        let samples = sample_barcodes
            .iter()
            .map(|s| s.sample.to_owned())
            .collect();
        let barcodes = sample_barcodes
            .iter()
            .map(|s| (s.sample.to_owned(), s.barcodes.to_owned()))
            .collect();
        Ok(SampleBarcodeCollection {
            samples: samples,
            barcodes: barcodes,
        })
    }

    fn barcodes_inversed(&self) -> HashMap<String, String> {
        self.barcodes
            .iter()
            .flat_map(|(sample, barcodes)| {
                barcodes
                    .iter()
                    .map(move |barcode| (barcode.to_owned(), sample.to_owned()))
            })
            .collect::<HashMap<String, String>>()
    }
}

pub struct SplitBam {
    path: PathBuf,
    barcodes: SampleBarcodeCollection,
    temp_dir: TempDir,
}

impl SplitBam {
    pub fn new<P>(path: P, barcodes: SampleBarcodeCollection) -> SplitBam
    where
        P: Into<PathBuf>,
    {   
        let path = path.into();
        let bam_name = path
            .file_name()
            .expect("Cannot read bam file name")
            .to_str()
            .expect("Cannot convert bam file name to string");

        SplitBam {
            path: path.clone(),
            barcodes: barcodes,
            temp_dir: TempDir::new(&format!("split_bam_{}", bam_name)).unwrap(),
        }
    }

    fn reader(&self) -> Result<bam::io::IndexedReader<bgzf::reader::Reader<File>>> {
        let reader = bam::io::indexed_reader::Builder::default().build_from_path(&self.path)?;
        Ok(reader)
    }

    fn chromosomes(&self) -> Result<Vec<String>> {
        let mut bam = self.reader()?;
        let header = bam.read_header()?;

        let chromosomes = header
            .reference_sequences()
            .iter()
            .map(|(chrom, _map)| chrom.to_string())
            .collect();

        Ok(chromosomes)
    }

    pub fn split_bam(&self, outdir: PathBuf) -> Result<()> {
        let chromosomes = self.chromosomes()?;
        let bam_files: HashMap<String, Vec<PathBuf>> = chromosomes
            .par_iter()
            .map(|chromosome| self._split_bam_by_chromosome(chromosome.to_owned()))
            .filter(|x| x.is_ok())
            .fold(
                || HashMap::new(),
                |mut acc: HashMap<String, Vec<PathBuf>>, x| {
                    let x = x.unwrap();
                    for (k, v) in x {
                        acc.entry(k).or_insert(vec![]).push(v);
                    }
                    acc
                },
            )
            .reduce(
                || HashMap::new(),
                |mut acc, x| {
                    for (k, v) in x {
                        acc.entry(k).or_insert(vec![]).extend(v);
                    }
                    acc
                },
            );

        for (sample, bam_paths) in bam_files {
            let output = self.temp_dir.path().join(format!("{}.bam", sample));
            let merged_bams = self._merge_bam(bam_paths, output)?;
            std::fs::rename(merged_bams, outdir.join(format!("{}.bam", sample)))?;
        }
        Ok(())
    }

    fn _merge_bam<P>(&self, bam_files: Vec<P>, output: P) -> Result<PathBuf>
    where
        P: Into<PathBuf> + Ord + Clone + FromStr,
    {
        let output = output.into();

        // Sort the bam file names
        let mut bam_files = bam_files;
        bam_files.sort();

        // Set up writer
        let header = self.reader()?.read_header()?;
        let mut writer = bam::io::Writer::new(std::io::BufWriter::new(File::create(output.clone())?));

        for bam_in in bam_files {
            let mut bam = bam::io::Reader::new(std::io::BufReader::new(File::open(bam_in.into())?));
            for result in bam.records() {
                let record = result?;
                writer.write_record(&header, &record)?;
            }
        }

        Ok(output)
    }

    fn _split_bam_by_chromosome(&self, chromosome: String) -> Result<HashMap<String, PathBuf>> {
        // Open bam file and extract required chromosome
        let mut bam_file = self.reader()?;
        let header = bam_file.read_header()?;

        // Create a bam file for each sample
        // Paths
        let bam_files: HashMap<String, PathBuf> = self
            .barcodes
            .samples
            .iter()
            .map(|sample| {
                let outfile = self
                    .temp_dir
                    .path()
                    .join(format!("{}_{}.bam", sample, chromosome));
                (sample.to_owned(), outfile)
            })
            .collect();

        // Writers
        let mut bam_writers: HashMap<
            String,
            bam::io::Writer<bgzf::Writer<std::io::BufWriter<File>>>,
        > = bam_files
            .iter()
            .map(|(sample, path)| {
                let _header = header.to_owned();
                let file = File::create(path).expect("Cannot create file");
                let writer = std::io::BufWriter::new(file);
                let bam_writer = bam::io::Writer::new(writer);
                (sample.to_owned(), bam_writer)
            })
            .collect();

        // Iterate through all reads in the bam file on the specified chromosome
        // If the barcodes provided are present, write to the specified bam file
        let barcodes_inversed = self.barcodes.barcodes_inversed();
        let region = &chromosome.parse()?;
        let query_region = bam_file.query(&header, region)?;

        for result in query_region {
            let record = result?;
            let record_data = record.data();
            let cell_barcode = record_data.get(b"CB");

            match cell_barcode {
                Some(Ok(cb)) => match cb {
                    Value::String(s) => {
                        let barcode = from_utf8(s).unwrap();
                        if barcodes_inversed.contains_key(barcode) {
                            let sample = barcodes_inversed.get(barcode).unwrap();
                            let writer = bam_writers.get_mut(sample).unwrap();
                            writer.write_record(&header, &record)?;
                        }
                    }
                    _ => {}
                },
                _ => {}
            };
        }

        Ok(bam_files)
    }
}

// fn split_region_reads_on_barcodes(
//     bam: String,
//     barcodes: HashMap<String, Vec<String>>,
//     chromosome: String,
//     outdir: String,
// ) -> Result<Vec<String>, rust_htslib::errors::Error> {
//     // Open bam file and extract required chromosome
//     let mut bam = bam::IndexedReader::from_path(bam)?;
//     let header = Header::from_template(bam.header());
//     bam.fetch(&chromosome)?;

//     // Sample: Writer
//     let mut bam_writers = barcodes
//         .keys()
//         .map(|sample| {
//             // Create bam file for each sample
//             let path = format!("{}/{}/{}.bam", outdir, sample, chromosome);
//             std::fs::create_dir_all(format!("{}/{}", outdir, sample)).unwrap();

//             let mut bam_writer = bam::Writer::from_path(path, &header, bam::Format::Bam)
//                 .expect(&format!("Failed to create {} file", &sample));
//             bam_writer.set_threads(1).unwrap();

//             (sample.to_owned(), bam_writer)
//         })
//         .collect::<HashMap<String, bam::Writer>>();

//     // Barcode: Sample
//     let barcodes_inversed = barcodes
//         .iter()
//         .flat_map(|(sample, barcodes)| {
//             barcodes
//                 .iter()
//                 .map(move |barcode| (barcode.to_owned(), sample.to_owned()))
//         })
//         .collect::<HashMap<String, String>>();

//     // Iterate through all reads in the bam file on the specified chromosome
//     // If the barcodes provided are present, write to the specified bam file
//     for result in bam.records() {
//         let record = result?;

//         match record.aux(b"CB") {
//             Ok(Aux::String(bc)) => {
//                 if barcodes_inversed.contains_key(bc) {
//                     let sample = barcodes_inversed.get(bc).unwrap();
//                     bam_writers
//                         .get_mut(sample)
//                         .unwrap()
//                         .write(&record.to_owned())?;
//                 }
//             }
//             _ => {}
//         }
//     }

//     let bam_files = barcodes.keys().map(|k| k.to_owned()).collect();

//     Ok(bam_files)
// }

// fn merge_bam_files(
//     bam_files: Vec<String>,
//     output: String,
//     n_threads: Option<usize>,
// ) -> Result<(), Error> {
//     let header = Header::from_template(bam::Reader::from_path(&bam_files[0]).unwrap().header());

//     let mut bam_writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;
//     let n_threads = n_threads.unwrap_or(1);
//     bam_writer.set_threads(n_threads).unwrap();

//     for bam_file in bam_files {
//         let mut bam = bam::Reader::from_path(&bam_file).unwrap();
//         for result in bam.records() {
//             let record = result.unwrap();
//             bam_writer.write(&record.to_owned()).unwrap();
//         }

//         // remove bam file
//         std::fs::remove_file(&bam_file).unwrap();
//     }

//     Ok(())
// }

// // Run through each chromosome of a bam file and split all reads that match a barcode in the hashmap into a new bam file
// pub fn split_bam(
//     bam_file: String,
//     barcodes: HashMap<String, Vec<String>>,
//     n_threads: usize,
//     output_prefix: String,
// ) -> Result<(), Error> {
//     // Split the bam file into separate bam files for each chromosome
//     // Merge the bam files for each sample together

//     let tmp_dir_prefix = format!("split_{}", &bam_file);
//     let tmp_dir = TempDir::new(&tmp_dir_prefix)?;

//     let bam = bam::IndexedReader::from_path(&bam_file).unwrap();
//     let header = bam.header();
//     let chromosomes: Vec<String> = header
//         .target_names()
//         .iter()
//         .map(|c| String::from_utf8(c.to_vec()).expect("Cannot read chromosome name"))
//         .collect();

//     chromosomes
//         .par_iter()
//         .for_each(|chromosome| {
//             split_region_reads_on_barcodes(
//                 bam_file.to_owned(),
//                 barcodes.to_owned(),
//                 chromosome.to_owned(),
//                 tmp_dir.path().to_str().unwrap().to_owned(),
//             );
//         })
//         .filter(|x| x.is_err())
//         .for_each(|x| x.unwrap());

//     // Merge the bam files for each chromosome together
//     let samples = barcodes
//         .keys()
//         .map(|k| k.to_owned())
//         .collect::<Vec<String>>();

//     let n_threads = match n_threads / samples.len() as usize {
//         0 => 1,
//         n => n,
//     };

//     samples.par_iter().for_each(|sample| {
//         let mut bam_files = vec![];

//         for chromosome in chromosomes.iter() {
//             bam_files.push(format!("{}/{}/{}.bam", output_prefix, sample, chromosome));
//         }

//         merge_bam_files(
//             bam_files,
//             format!("{}/{}.bam", output_prefix, sample),
//             Some(n_threads),
//         )
//         .unwrap();
//     });

//     Ok(())
// }
