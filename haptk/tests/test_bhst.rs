mod common;

#[cfg(test)]
#[cfg(feature = "clap")]
mod test_bhst {
    use super::*;
    use haptk::{args::Selection, subcommands::bhst_shard::read_hst_file};
    use std::path::PathBuf;

    use color_eyre::Result;

    //--- Binary tests

    #[test]
    fn bhst_all() {
        run_bhst(Selection::All);

        let res = std::fs::read_to_string("tests/results/bhst_mbah.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_shared_core_haplotype.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn bhst_only_ref() {
        run_bhst(Selection::OnlyRefs);

        let res = std::fs::read_to_string("tests/results/bhst_mbah_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_shared_core_haplotype_only_refs.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn bhst_only_alt() {
        run_bhst(Selection::OnlyAlts);

        let res = std::fs::read_to_string("tests/results/bhst_mbah_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/bhst_shared_core_haplotype_only_alts.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn bhst_only_longest() {
        run_bhst(Selection::OnlyLongest);

        let res = std::fs::read_to_string("tests/results/bhst_mbah_only_longest.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res =
            std::fs::read_to_string("tests/results/bhst_shared_core_haplotype_only_longest.csv")
                .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn test_bhst_core_ht() {
        let args = haptk::args::StandardArgs {
            file: PathBuf::from("tests/data/test_bhst_core_ht.vcf.gz"),
            output: PathBuf::from("tests/results"),
            coords: String::from("chr9:32"),
            selection: haptk::args::Selection::All,
            samples: None,
            info_limit: None,
            prefix: Some("test_core".to_string()),
        };

        let cmd = haptk::clap::SubCommand::Bhst {
            args,
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
            sharded: false,
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
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
            sharded: false,
        };
        haptk::clap::run_cmd(cmd).unwrap();
    }

    #[test]
    fn bhst_sharded() -> Result<()> {
        let mut args = common::clap_standard_args(Selection::All);

        args.file = PathBuf::from("tests/data/test_sharded.vcf.gz");
        args.prefix = Some(String::from("sharded_one_run"));
        args.coords = String::from("chr9:32000000");

        let cmd = haptk::clap::SubCommand::Bhst {
            args: args.clone(),
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
            sharded: false,
        };
        haptk::clap::run_cmd(cmd)?;

        println!("Constructed normal tree");

        args.prefix = Some(String::from("sharded_multi_run"));

        let cmd = haptk::clap::SubCommand::Bhst {
            args,
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
            sharded: true,
        };
        haptk::clap::run_cmd(cmd)?;

        let hst1 = read_hst_file(PathBuf::from("tests/results/sharded_multi_run_bhst.hst.gz"))?;
        let hst2 = read_hst_file(PathBuf::from("tests/results/sharded_one_run_bhst.hst.gz"))?;

        assert_eq!(hst1.hst.node_count(), hst2.hst.node_count());
        assert_eq!(hst1.hst.edge_count(), hst2.hst.edge_count());

        for (idx1, idx2) in hst1.hst.node_indices().zip(hst2.hst.node_indices()) {
            assert_eq!(idx1, idx2);

            let data1 = hst1.hst.node_weight(idx1).unwrap();
            let data2 = hst2.hst.node_weight(idx2).unwrap();
            assert_eq!(data1, data2);
        }

        Ok(())
    }
}
