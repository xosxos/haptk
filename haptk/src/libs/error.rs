use std::path::PathBuf;

use thiserror::Error as ThisError;

#[rustfmt::skip]
#[derive(ThisError, Debug)]
pub enum Error {
    #[error("Failed to open file: {path:?}")]
    Io { path: PathBuf },

    #[error("Failed to parse coords: {coord}")]
    CoordParse { coord: String },

    #[error("Position {value:?} is not an integer in coords {coord}")]
    PosParse { coord: String, value: String },

    #[error("No variants were found at given coordinate {contig}:{variant_pos}")]
    NonExistantVariant { contig: String, variant_pos: u64 },

    #[error("At pos {pos} allele count != 2. Normalize alleles using bcftools norm")]
    Normalize { pos: u64 },

    #[error("The VCF file is not sorted: {prev_pos} > {pos} at {coord}")]
    Order { prev_pos: u64, pos: u64, coord: String, },

    #[error("The file is not sorted: {prev_pos} > {pos}, file: {path:?}")]
    Sort { prev_pos: f32, pos: f32, path: PathBuf },

    #[error("None of the wanted samples was found in the vcf.")]
    SamplesNotFound,

    #[error("File type: {ext} is not supported")]
    FileNotSupported { ext: String },

    #[error("VCF header has no contig length for {contig}")]
    NoContigLength { contig: String },

    #[error("VCF header has no contig {contig}")]
    NoContig { contig: String },

    #[error("No file type extension in path: {path}")]
    NoFileType { path: PathBuf },

    #[error("File contains zero rows: {path}")]
    EmptyFile { path: PathBuf },

    #[error("Make sure no headers are present and that the recombination file is in order chr,pos,rate,cm. If the issue is not fixed, you have an invalid field in the file.")]
    RecombinationDeserialization { path: PathBuf},
    
    #[error("VCF has less haplotypes than the required minimum node size ({n_haplotypes} < {min_size})" )]
    MinNodes { n_haplotypes: usize, min_size: usize },

    #[error("The HST has less than 3 nodes, use other methods to study haplotype sharing")]
    HstTooSmall,

    #[error("Missing genotypes at position: {pos}")]
    MissingGenotypes { pos: u64},

    #[error("None of the variants in the case file are present in the controls file.")]
    NoMatchingVariants,

    #[error("HST construction ended prematurely due to windowed file-read")]
    HstEnd,

    #[error("HST construction ended prematurely on both sides due to windowed file-read")]
    HstBothEnd,

    #[error("HST construction ended prematurely to the right due to windowed file-read")]
    HstRightEnd,

    #[error("HST construction ended prematurely to the left due to windowed file-read")]
    HstLeftEnd,

    #[error("Sample does not have the amount of genotypes required by ploidy at position: {pos}. Missing genotypes are not allowed. {num}" )]
    Ploidy { num: usize, pos: usize },

    #[error("The given coordinate {variant_pos} is larger than the largest found position. Total records read: {records_n}. Check your coordinates and vcf file. Comparing to a haplotype file also automatically limits min and max coordinates to the haplotype coordinates.")]
    VariantPosNotFound { variant_pos: u64, records_n: usize },
}
