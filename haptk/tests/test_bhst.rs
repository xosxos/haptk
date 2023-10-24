mod common;

#[cfg(test)]
#[cfg(feature = "clap")]
mod test_bhst {
    use super::*;
    use haptk::args::Selection;
    use std::path::PathBuf;

    //--- Binary tests

    #[test]
    fn bhst_all() {
        run_bhst(Selection::All);

        let res = std::fs::read_to_string("tests/results/bhst.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_mbah.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_shared_core_haplotype.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn bhst_only_ref() {
        run_bhst(Selection::OnlyRefs);

        let res = std::fs::read_to_string("tests/results/bhst_only_refs.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_mbah_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_shared_core_haplotype_only_refs.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn bhst_only_alt() {
        run_bhst(Selection::OnlyAlts);

        let res = std::fs::read_to_string("tests/results/bhst_only_alts.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_mbah_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_shared_core_haplotype_only_alts.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn bhst_only_longest() {
        run_bhst(Selection::OnlyLongest);

        let res = std::fs::read_to_string("tests/results/bhst_only_longest.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_mbah_only_longest.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res =
            std::fs::read_to_string("tests/results/bhst_shared_core_haplotype_only_longest.csv")
                .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn test_bhst_core_ht() {
        let args = haptk::clap::ClapStandardArgs {
            file: PathBuf::from("tests/data/test_bhst_core_ht.vcf.gz"),
            outdir: PathBuf::from("tests/results"),
            coords: String::from("chr9:32"),
            select: haptk::clap::ClapSelection::All,
            samples: None,
            // info_limit: None,
            prefix: Some("test_core".to_string()),
        };

        let cmd = haptk::clap::SubCommand::Bhst {
            args,
            graph_args: haptk::clap::ClapGraphArgs::default(),
            mark_samples: None,
            variable_data: None,
            variables: None,
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
        };
        haptk::clap::run_cmd(cmd).unwrap();

        let res = std::fs::read_to_string("tests/results/test_core_bhst_shared_core_haplotype.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    fn run_bhst(selection: Selection) {
        let args = common::clap_standard_args(selection);

        let cmd = haptk::clap::SubCommand::Bhst {
            args,
            graph_args: haptk::clap::ClapGraphArgs::default(),
            mark_samples: None,
            variable_data: Some(PathBuf::from("tests/data/clinical_data.csv")),
            variables: Some(vec!["aoo".to_string()]),
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
        };
        haptk::clap::run_cmd(cmd).unwrap();
    }
}
