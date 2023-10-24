mod common;

use std::path::PathBuf;

use haptk::subcommands::read_sample_names::{self, get_sample_names};

use crate::common::TEST_VCF;

#[test]
#[cfg(feature = "clap")]
fn samples() {
    let cmd = haptk::clap::SubCommand::Samples {
        file: PathBuf::from(TEST_VCF),
        log_and_verbosity: crate::common::silent_verbosity(),
    };
    haptk::clap::run_cmd(cmd).unwrap();
}

#[test]
fn sample_names_subcommand() {
    let path = PathBuf::from(TEST_VCF);
    let ids = get_sample_names(path).unwrap();
    let vec: Vec<_> = (1..15).map(|v| format!("SAMPLE{v}")).collect();
    assert_eq!(vec, ids);

    let path = PathBuf::from("tests/data/samples.fam");
    let ids = get_sample_names(path).unwrap();
    let vec = vec!["SAMPLE1", "SAMPLE2", "SAMPLE3", "SAMPLE4"];
    assert_eq!(vec, ids);

    let path = PathBuf::from("tests/data/samples.fam");
    let res = read_sample_names::run(path);
    assert!(res.is_ok());
}
