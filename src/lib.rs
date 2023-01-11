use pyo3::prelude::*;
use std::collections::HashMap;
use std::prelude::rust_2021::*;

pub mod split_bam;

/// Splits a bam file by provided barcodes
/// Barcodes must be in the format {SAMPLE: [BARCODE_1...BARCODE_N]}
#[pyfunction]
#[pyo3(text_signature = "(bam_file, barcodes, n_threads, output_prefix)")]
fn split_bam_file(
    bam_file: String,
    barcodes: HashMap<String, Vec<String>>,
    n_threads: usize,
    output_prefix: String,
) -> PyResult<()> {
    split_bam::split_bam(bam_file, barcodes, n_threads, output_prefix)
        .expect("Error splitting bam file");
    Ok(())
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn single_cell_tools(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(split_bam_file, m)?)?;
    Ok(())
}
