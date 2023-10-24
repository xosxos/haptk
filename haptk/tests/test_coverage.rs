mod common;
use std::path::PathBuf;

#[test]
#[cfg(feature = "clap")]
fn coverage() {
    let cmd = haptk::clap::SubCommand::Coverage {
        file: PathBuf::from(common::TEST_VCF),
        npipes: 10,
        bp_per_snp: 2,
        threads: 1,
        log_and_verbosity: crate::common::silent_verbosity(),
    };
    let res = haptk::clap::run_cmd(cmd);
    // Contig has under 100 records
    assert!(res.is_err());
}
