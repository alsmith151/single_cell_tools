use pyo3::prelude::*;
use std::collections::HashMap;
use std::prelude::rust_2021::*;

pub mod split_bam;

#[pyclass]
struct BamSplitter {
    bam_file: String,
    barcodes: HashMap<String, Vec<String>>,
    output_prefix: String,
}

#[pymethods]
impl BamSplitter {
    #[new]
    fn new(
        bam_file: String,
        barcodes: HashMap<String, Vec<String>>,
        output_prefix: String,
    ) -> Self {
        BamSplitter {
            bam_file,
            barcodes,
            output_prefix,
        }
    }

    fn split(&self, _n_threads: usize) -> PyResult<()> {
        let barcodes: Vec<split_bam::SampleBarcodes> = self
            .barcodes
            .iter()
            .map(|(sample, bc)| split_bam::SampleBarcodes::new(sample.clone(), bc.clone()))
            .collect();

        let barcode_collection = split_bam::SampleBarcodeCollection::from_sample_barcodes(barcodes)
            .expect("Error creating barcode collection");

        let bam_splitter = split_bam::SplitBam::new(self.bam_file.clone(), barcode_collection);

        bam_splitter
            .split_bam(self.output_prefix.clone().into())
            .expect("Error splitting bam file");

        Ok(())
    }
}

// /// Splits a bam file by provided barcodes
// /// Barcodes must be in the format {SAMPLE: [BARCODE_1...BARCODE_N]}
// #[pyfunction]
// #[pyo3(text_signature = "(bam_file, barcodes, n_threads, output_prefix)")]
// fn split_bam_file(
//     bam_file: String,
//     barcodes: HashMap<String, Vec<String>>,
//     n_threads: usize,
//     output_prefix: String,
// ) -> PyResult<()> {
//     split_bam::split_bam(bam_file, barcodes, n_threads, output_prefix)
//         .expect("Error splitting bam file");
//     Ok(())
// }

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn scnado(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(split_bam_file, m)?)?;

    m.add_class::<BamSplitter>()?;
    Ok(())
}
