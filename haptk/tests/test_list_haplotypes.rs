mod common;

use crate::common::TEST_VCF;
use std::path::PathBuf;

#[test]
#[cfg(feature = "clap")]
fn read_haplotypes() {
    use common::COORD_RANGE;
    use haptk::args::{Selection, StandardArgs};

    let cmd = haptk::clap::SubCommand::Haplotypes {
        args: StandardArgs {
            file: PathBuf::from(TEST_VCF),
            output: PathBuf::from("tests/results"),
            coords: String::from(COORD_RANGE),
            selection: Selection::All,
            info_limit: None,
            prefix: None,
            samples: Some(vec![PathBuf::from("tests/data/SAMPLE1.ids")]),
            no_alt: false,
        },
        selection_variant: None,
        log_and_verbosity: crate::common::silent_verbosity(),
        nucleotides: false,
    };
    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/haplotype_file.csv").unwrap();
    insta::assert_yaml_snapshot!(res);

    let res = std::fs::read_to_string("tests/results/genotype_file.csv").unwrap();
    insta::assert_yaml_snapshot!(res);
}

#[test]
#[cfg(feature = "clap")]
fn read_haplotypes_only_alt() {
    use common::COORD_RANGE;
    use haptk::args::{Selection, StandardArgs};

    let cmd = haptk::clap::SubCommand::Haplotypes {
        args: StandardArgs {
            file: PathBuf::from(TEST_VCF),
            output: PathBuf::from("tests/results"),
            coords: String::from(COORD_RANGE),
            selection: Selection::OnlyAlts,
            info_limit: None,
            prefix: None,
            samples: None,
            no_alt: false,
        },
        selection_variant: Some(String::from("chr9:32")),
        log_and_verbosity: crate::common::silent_verbosity(),
        nucleotides: false,
    };
    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/haplotype_file_only_alts.csv").unwrap();
    insta::assert_yaml_snapshot!(res);

    let res = std::fs::read_to_string("tests/results/genotype_file_only_alts.csv").unwrap();
    insta::assert_yaml_snapshot!(res);
}

#[test]
#[cfg(feature = "clap")]
fn read_haplotypes_only_ref() {
    use common::COORD_RANGE;
    use haptk::args::{Selection, StandardArgs};

    let cmd = haptk::clap::SubCommand::Haplotypes {
        args: StandardArgs {
            file: PathBuf::from(TEST_VCF),
            output: PathBuf::from("tests/results"),
            coords: String::from(COORD_RANGE),
            selection: Selection::OnlyRefs,
            info_limit: None,
            prefix: None,
            samples: None,
            no_alt: false,
        },
        selection_variant: Some(String::from("chr9:32")),
        log_and_verbosity: crate::common::silent_verbosity(),
        nucleotides: false,
    };
    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/haplotype_file_only_refs.csv").unwrap();
    insta::assert_yaml_snapshot!(res);

    let res = std::fs::read_to_string("tests/results/genotype_file_only_refs.csv").unwrap();
    insta::assert_yaml_snapshot!(res);
}

#[test]
#[cfg(feature = "clap")]
fn read_haplotypes_only_longest() {
    use common::COORD_RANGE;
    use haptk::args::{Selection, StandardArgs};

    let cmd = haptk::clap::SubCommand::Haplotypes {
        args: StandardArgs {
            file: PathBuf::from(TEST_VCF),
            output: PathBuf::from("tests/results"),
            coords: String::from(COORD_RANGE),
            selection: Selection::OnlyLongest,
            info_limit: None,
            prefix: None,
            samples: None,
            no_alt: false,
        },
        selection_variant: Some(String::from("chr9:32")),
        log_and_verbosity: crate::common::silent_verbosity(),
        nucleotides: false,
    };
    haptk::clap::run_cmd(cmd).unwrap();

    let res = std::fs::read_to_string("tests/results/haplotype_file_only_longest.csv").unwrap();
    insta::assert_yaml_snapshot!(res);

    let res = std::fs::read_to_string("tests/results/genotype_file_only_longest.csv").unwrap();
    insta::assert_yaml_snapshot!(res);
}
