mod common;
use common::{standard_args, OUTDIR};

use std::path::PathBuf;

use haptk::{args::Selection, subcommands::uhst};

#[test]
#[cfg(feature = "clap")]
fn compare_haplotypes() {
    let args = standard_args(Selection::All);
    uhst::run(args, 1, false, None).unwrap();

    let args = standard_args(Selection::OnlyLongest);
    uhst::run(args, 1, false, None).unwrap();

    let files = vec![
        PathBuf::from("tests/results/uhst_mbah.csv"),
        PathBuf::from("tests/results/uhst_shared_core_haplotype_only_longest.csv"),
    ];
    let output = PathBuf::from(OUTDIR);
    let cmd = haptk::clap::SubCommand::CompareHaplotypes {
        haplotypes: files,
        output,
        prefix: None,
        csv: true,
        hide_missing: false,
        tag_rows: false,
        log_and_verbosity: crate::common::silent_verbosity(),
        nucleotides: false,
    };

    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/ht_comparison.csv").unwrap();
    insta::assert_yaml_snapshot!(res);
}
