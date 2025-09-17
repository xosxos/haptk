#![allow(dead_code)]
use std::path::PathBuf;

use color_eyre::Result;

use haptk::{
    args::{Selection, StandardArgs},
    clap::LogAndVerbosity,
    read_vcf::read_vcf_to_matrix,
    structs::PhasedMatrix,
};

pub const TEST_VCF: &str = "tests/data/test.vcf.gz";
pub const TEST_HAPLOTYPE: &str = "tests/data/test_haplotype.csv";
pub const TEST_FASTA: &str = "tests/data/unc13a_exons.fa";
pub const OUTDIR: &str = "tests/results";
pub const REC_RATES: &str = "tests/data/recombination_rates.tsv";
pub const TEST_CONTRA: &str = "tests/data/test_contrahomozygote.vcf.gz";
pub const COORDS: &str = "chr9:32";
pub const COORD_RANGE: &str = "chr9:0-70";

pub fn create_test_matrix() -> Result<PhasedMatrix> {
    let args = StandardArgs {
        file: PathBuf::from(TEST_VCF),
        ..Default::default()
    };
    read_vcf_to_matrix(&args, "chr9", 32, None, None, None, false)
}

pub fn standard_args(selection: Selection) -> StandardArgs {
    StandardArgs {
        file: PathBuf::from(TEST_VCF),
        output: PathBuf::from("tests/results"),
        coords: String::from(COORDS),
        selection,
        ..Default::default()
    }
}

#[cfg(feature = "clap")]
pub fn clap_standard_args(select: Selection) -> haptk::args::StandardArgs {
    haptk::args::StandardArgs {
        file: PathBuf::from(TEST_VCF),
        output: PathBuf::from("tests/results"),
        coords: String::from(COORDS),
        selection: select,
        samples: None,
        prefix: None,
        no_alt: false,
        list: None,
    }
}

#[cfg(feature = "clap")]
pub fn silent_verbosity() -> LogAndVerbosity {
    LogAndVerbosity {
        verbosity: 1,
        log_file: None,
        silent: false,
    }
}
