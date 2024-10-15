mod common;

use std::path::PathBuf;

use color_eyre::Result;

use haptk::{
    args::Selection,
    subcommands::{
        list_markers::read_hst_coords,
        list_samples::{self, get_sample_names},
    },
};

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
    let res = list_samples::run(path);
    assert!(res.is_ok());
}

#[test]
fn list_markers_subcommand() -> Result<()> {
    let mut args = common::clap_standard_args(Selection::All);

    args.prefix = Some(String::from("list_markers"));
    args.coords = String::from("chr9:32000000");

    let cmd = haptk::clap::SubCommand::Bhst {
        args: args.clone(),
        log_and_verbosity: crate::common::silent_verbosity(),
        threads: 1,
        min_size: 1,
        publish: false,
        window: None,
    };
    haptk::clap::run_cmd(cmd)?;

    let coords = read_hst_coords(PathBuf::from("tests/results/bhst.hst.gz"))?;
    assert_eq!(coords.len(), 63);

    Ok(())
}
