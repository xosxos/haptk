mod common;

#[cfg(test)]
#[cfg(feature = "clap")]
mod test_uhst {
    use super::*;

    use std::io::Write;
    use std::path::PathBuf;

    use haptk::{
        args::{Selection, StandardArgs},
        read_vcf::read_vcf_to_matrix,
        subcommands::uhst::{self, LocDirection},
    };

    #[test]
    fn matrix_to_graph() {
        let args = StandardArgs {
            file: PathBuf::from(common::TEST_VCF),
            ..Default::default()
        };
        let vcf = read_vcf_to_matrix(&args, "chr9", 32, None, None).unwrap();

        let g = uhst::construct_uhst(&vcf, &LocDirection::Left, 32, 1, false);
        let mut f = std::fs::File::create("tests/results/test.dot").unwrap();
        f.write_all(format!("{}", petgraph::dot::Dot::new(&g)).as_bytes())
            .unwrap();
    }

    //--- Binary tests

    #[test]
    fn uhst_all() {
        run_uhst(Selection::All);

        let res = std::fs::read_to_string("tests/results/uhst_left.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_right.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_mbah.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_shared_core_haplotype.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn uhst_only_ref() {
        run_uhst(Selection::OnlyRefs);

        let res = std::fs::read_to_string("tests/results/uhst_left_only_refs.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_right_only_refs.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_mbah_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_shared_core_haplotype_only_refs.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn uhst_only_alt() {
        run_uhst(Selection::OnlyAlts);

        let res = std::fs::read_to_string("tests/results/uhst_left_only_alts.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_right_only_alts.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_mbah_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_shared_core_haplotype_only_alts.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn uhst_only_longest() {
        run_uhst(Selection::OnlyLongest);

        let res = std::fs::read_to_string("tests/results/uhst_left_only_longest.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_right_only_longest.svg").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_mbah_only_longest.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res =
            std::fs::read_to_string("tests/results/uhst_shared_core_haplotype_only_longest.csv")
                .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    fn run_uhst(selection: Selection) {
        let args = common::clap_standard_args(selection);

        let cmd = haptk::clap::SubCommand::Uhst {
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
