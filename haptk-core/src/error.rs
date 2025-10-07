use std::path::PathBuf;

use thiserror::Error as ThisError;

#[rustfmt::skip]
#[derive(ThisError, Debug)]
pub enum Error {
    #[error("Failed to parse coords: {coord}")]
    CoordParse { coord: String },

    #[error("Position: {value:?} is not an integer in coords {coord}")]
    PosParse { coord: String, value: String },

    #[error("Htslib error: {0}")]
    HtsLib(#[from] rust_htslib::errors::Error),

    #[error("Io error: {0} {1}")]
    Io(PathBuf, std::io::Error),

    #[error("None of the samples are found in the vcf")]
    SamplesNotFound,

    #[error("Io error: {0}")]
    StdIo(#[from] std::io::Error),

    #[error("Input file was too short")]
    FileTooShort,
    
    #[error("Io error: {0} {1}")]
    Path(PathBuf, String),

    #[error("The VCF file is not sorted: {prev_pos} > {pos} at {coord}")]
    Order { prev_pos: u64, pos: u64, coord: String, },

    #[error("Unknown file extension: {0:?}")]
    UnknownExtension(PathBuf),

    #[error("No file type extension in path: {path}")]
    NoFileType { path: PathBuf },

    #[error("File contains zero rows: {path}")]
    EmptyFile { path: PathBuf },

    #[error("File type: {ext} is not supported")]
    FileNotSupported { ext: String },

    #[error("VCF header has no contig length for {contig}")]
    NoContigLength { contig: String },

    #[error("Error parsing VCF header: {0}")]
    HeaderParse(#[from] std::num::ParseIntError),

    #[error("Something failed in reading variants from VCF {0}")]
    ShapeError(#[from] ndarray::ShapeError),

    #[error("At pos {pos} allele count != 2. Normalize alleles using bcftools norm")]
    Normalize { pos: u64 },

    #[error("{0}")]
    New(String),

    #[cfg(feature = "noodles")]
    #[error("parse error: {0}")]
    Parse(#[from] noodles::sam::header::record::value::map::tag::ParseError),

}
