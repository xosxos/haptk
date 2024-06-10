mod common;
use common::{standard_args, OUTDIR};

use std::path::PathBuf;

use haptk::{
    args::{GraphArgs, Selection},
    subcommands::uhst,
};

#[test]
#[cfg(feature = "clap")]
fn compare_haplotypes() {
    let args = standard_args(Selection::All);
    uhst::run(args, GraphArgs::default(), None, None, None, 1, false, true).unwrap();

    let args = standard_args(Selection::OnlyLongest);
    uhst::run(args, GraphArgs::default(), None, None, None, 1, false, true).unwrap();

    let files = vec![
        PathBuf::from("tests/results/uhst_mbah.csv"),
        PathBuf::from("tests/results/uhst_shared_core_haplotype_only_longest.csv"),
    ];
    let outdir = PathBuf::from(OUTDIR);
    let cmd = haptk::clap::SubCommand::CompareHaplotypes {
        haplotypes: files,
        outdir,
        prefix: None,
        csv: true,
        hide_missing: false,
        tag_rows: false,
        log_and_verbosity: crate::common::silent_verbosity(),
    };

    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/ht_comparison.csv").unwrap();
    insta::assert_yaml_snapshot!(res);
}
