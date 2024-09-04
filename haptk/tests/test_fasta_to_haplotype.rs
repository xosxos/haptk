mod common;

use crate::common::TEST_FASTA;
use std::path::PathBuf;

#[test]
#[cfg(feature = "clap")]
fn fasta_to_haplotype() {
    let cmd = haptk::clap::SubCommand::FastaToHaplotype {
        file: PathBuf::from(TEST_FASTA),
        seq_name: vec![String::from(
            "chr19:17676011-17676041_ENST00000519716.7_exon_1_0_chr19_17676012_r",
        )],
        log_and_verbosity: crate::common::silent_verbosity(),
        output: PathBuf::from("tests/results/exon_haplotype.csv"),
    };

    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/exon_haplotype.csv").unwrap();
    insta::assert_yaml_snapshot!(res);
}
