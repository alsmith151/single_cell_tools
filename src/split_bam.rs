use rayon::prelude::*;
use rust_htslib::{bam, bam::record::Aux, bam::Header, bam::Read, errors::Error};
use std::fs::File;
use std::io::prelude::*;
use std::str::from_utf8;
use std::str::{Bytes, FromStr};
use std::{collections::HashMap, result};

fn split_region_reads_on_barcodes(
    bam: String,
    barcodes: HashMap<String, Vec<String>>,
    chromosome: String,
    outdir: String,
) -> Result<Vec<String>, rust_htslib::errors::Error> {
    // Open bam file and extract required chromosome
    let mut bam = bam::IndexedReader::from_path(bam)?;
    let header = Header::from_template(bam.header());
    bam.fetch(&chromosome)?;

    // Sample: Writer
    let mut bam_writers = barcodes
        .keys()
        .map(|sample| {
            // Create bam file for each sample
            let path = format!("{}/{}/{}.bam", outdir, sample, chromosome);
            std::fs::create_dir_all(format!("{}/{}", outdir, sample)).unwrap();

            let mut bam_writer = bam::Writer::from_path(path, &header, bam::Format::Bam)
                .expect(&format!("Failed to create {} file", &sample));
            bam_writer.set_threads(1).unwrap();

            (sample.to_owned(), bam_writer)
        })
        .collect::<HashMap<String, bam::Writer>>();

    // Barcode: Sample
    let barcodes_inversed = barcodes
        .iter()
        .flat_map(|(sample, barcodes)| {
            barcodes
                .iter()
                .map(move |barcode| (barcode.to_owned(), sample.to_owned()))
        })
        .collect::<HashMap<String, String>>();

    // Iterate through all reads in the bam file on the specified chromosome
    // If the barcodes provided are present, write to the specified bam file
    for result in bam.records() {
        let record = result?;

        match record.aux(b"CB") {
            Ok(Aux::String(bc)) => {
                if barcodes_inversed.contains_key(bc) {
                    let sample = barcodes_inversed.get(bc).unwrap();
                    bam_writers
                        .get_mut(sample)
                        .unwrap()
                        .write(&record.to_owned())?;
                }
            }
            _ => {}
        }
    }

    let bam_files = barcodes.keys().map(|k| k.to_owned()).collect();

    Ok(bam_files)
}


fn merge_bam_files(
    bam_files: Vec<String>,
    output: String,
    n_threads: Option<usize>,
) -> Result<(), Error> {
    let header = Header::from_template(bam::Reader::from_path(&bam_files[0]).unwrap().header());

    let mut bam_writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;
    let n_threads = n_threads.unwrap_or(1);
    bam_writer.set_threads(n_threads).unwrap();

    for bam_file in bam_files {
        let mut bam = bam::Reader::from_path(&bam_file).unwrap();
        for result in bam.records() {
            let record = result.unwrap();
            bam_writer.write(&record.to_owned()).unwrap();
        }

        // remove bam file
        std::fs::remove_file(&bam_file).unwrap();


    }

    Ok(())
}

// Run through each chromosome of a bam file and split all reads that match a barcode in the hashmap into a new bam file
pub fn split_bam(
    bam_file: String,
    barcodes: HashMap<String, Vec<String>>,
    n_threads: usize,
    output_prefix: String,
) -> Result<(), Error> {

    // Split the bam file into separate bam files for each chromosome
    // Merge the bam files for each sample together


    let bam = bam::IndexedReader::from_path(&bam_file).unwrap();
    let header = bam.header();
    let chromosomes: Vec<String> = header
        .target_names()
        .iter()
        .map(|c| String::from_utf8(c.to_vec()).expect("Cannot read chromosome name"))
        .collect();

    chromosomes.par_iter().for_each(|chromosome| {
        split_region_reads_on_barcodes(
            bam_file.to_owned(),
            barcodes.to_owned(),
            chromosome.to_owned(),
            output_prefix.to_owned(),
        )
        .unwrap();
    });

    // Merge the bam files for each chromosome together
    let samples = barcodes
        .keys()
        .map(|k| k.to_owned())
        .collect::<Vec<String>>();

    let n_threads = n_threads / samples.len() as usize;

    samples.par_iter().for_each(|sample| {
        let mut bam_files = vec![];

        for chromosome in chromosomes.iter() {
            bam_files.push(format!("{}/{}/{}.bam", output_prefix, sample, chromosome));
        }

        merge_bam_files(
            bam_files,
            format!("{}/{}.bam", output_prefix, sample),
            Some(n_threads),
        )
        .unwrap();
    });

    // remove the sample directories
    // for sample in samples {
    //     std::fs::remove_dir_all(format!("{}/{}", output_prefix, sample)).unwrap();
    // }
    

    Ok(())
}
