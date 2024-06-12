#![allow(dead_code)]
use std::path::PathBuf;

use haptk::{
    args::{Selection, StandardArgs},
    clap::LogAndVerbosity,
    read_vcf::read_vcf_to_matrix,
    structs::PhasedMatrix,
};

pub const TEST_VCF: &str = "tests/data/test.vcf.gz";
pub const TEST_HAPLOTYPE: &str = "tests/data/test_haplotype.csv";
pub const OUTDIR: &str = "tests/results";
pub const REC_RATES: &str = "tests/data/recombination_rates.tsv";
pub const TEST_CONTRA: &str = "tests/data/test_contrahomozygote.vcf.gz";
pub const COORDS: &str = "chr9:32";

pub fn create_test_matrix() -> PhasedMatrix {
    let args = StandardArgs {
        file: PathBuf::from(TEST_VCF),
        ..Default::default()
    };
    read_vcf_to_matrix(&args, "chr9", 32, None, None).unwrap()
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
pub fn clap_standard_args(select: Selection) -> haptk::clap::ClapStandardArgs {
    haptk::clap::ClapStandardArgs {
        file: PathBuf::from(TEST_VCF),
        outdir: PathBuf::from("tests/results"),
        coords: String::from(COORDS),
        alleles: select.into(),
        samples: None,
        // info_limit: None,
        prefix: None,
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
