mod common;
use std::path::PathBuf;

use haptk::args::Selection;
#[test]
#[cfg(feature = "clap")]
fn mrca() {
    let args = common::clap_standard_args(Selection::All);

    let cmd = haptk::clap::SubCommand::Mrca {
        args,
        recombination_rates: PathBuf::from(common::REC_RATES),
        log_and_verbosity: crate::common::silent_verbosity(),
        start: None,
        stop: None,
    };
    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/mrca_gamma_method.txt").unwrap();
    insta::assert_yaml_snapshot!(res);
}
