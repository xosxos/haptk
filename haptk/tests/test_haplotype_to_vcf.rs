mod common;

use crate::common::TEST_HAPLOTYPE;
use std::path::PathBuf;

#[test]
#[cfg(feature = "clap")]
fn haplotype_to_vcf() {
    let cmd = haptk::clap::SubCommand::HaplotypeToVcf {
        file: PathBuf::from(TEST_HAPLOTYPE),
        sample_name: String::from("TEST1"),
        log_and_verbosity: crate::common::silent_verbosity(),
        output: PathBuf::from("tests/results/TEST_haplotype.vcf"),
    };

    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/TEST_haplotype.vcf").unwrap();
    let res = res.split("ID=chr9").nth(1).unwrap();
    insta::assert_yaml_snapshot!(res);
}
