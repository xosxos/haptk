mod common;

#[cfg(test)]
#[cfg(feature = "clap")]
mod test_uhst {
    use std::path::PathBuf;

    use color_eyre::Result;

    use super::*;
    use haptk::{args::Selection, subcommands::bhst::read_hst_file};

    #[test]
    fn uhst_all() -> Result<()> {
        run_uhst(Selection::All)?;

        let res = std::fs::read_to_string("tests/results/ancestral_haplotype.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/core_haplotype.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
        Ok(())
    }

    #[test]
    fn uhst_only_ref() -> Result<()> {
        run_uhst(Selection::OnlyRefs)?;

        let res =
            std::fs::read_to_string("tests/results/ancestral_haplotype_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/core_haplotype_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
        Ok(())
    }

    #[test]
    fn uhst_only_alt() -> Result<()> {
        run_uhst(Selection::OnlyAlts)?;

        let res =
            std::fs::read_to_string("tests/results/ancestral_haplotype_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/core_haplotype_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
        Ok(())
    }

    #[test]
    fn uhst_only_longest() -> Result<()> {
        run_uhst(Selection::OnlyLongest)?;

        let res =
            std::fs::read_to_string("tests/results/ancestral_haplotype_only_longest.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/core_haplotype_only_longest.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
        Ok(())
    }

    fn run_uhst(selection: Selection) -> Result<()> {
        let args = common::clap_standard_args(selection);

        let cmd = haptk::clap::SubCommand::Hst {
            args,
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
            window: 20_000,
        };
        haptk::clap::run_cmd(cmd)
    }

    #[test]
    fn uhst_sharded() -> Result<()> {
        let mut args = common::clap_standard_args(Selection::All);

        args.file = PathBuf::from("tests/data/test_sharded.vcf.gz");
        args.prefix = Some(String::from("sharded_one_run"));
        args.coords = String::from("chr9:32000000");

        let cmd = haptk::clap::SubCommand::Hst {
            args: args.clone(),
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 1,
            min_size: 1,
            publish: false,
            window: 20_000,
        };
        haptk::clap::run_cmd(cmd)?;

        println!("Constructed normal tree");

        args.prefix = Some(String::from("sharded_multi_run"));

        let cmd = haptk::clap::SubCommand::Hst {
            args,
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 1,
            min_size: 1,
            publish: false,
            window: 20_000,
        };
        haptk::clap::run_cmd(cmd)?;

        let hst1 = read_hst_file(PathBuf::from("tests/results/sharded_multi_run_left.hst.gz"))?;
        let hst2 = read_hst_file(PathBuf::from("tests/results/sharded_one_run_left.hst.gz"))?;

        assert_eq!(hst1.hst.node_count(), hst2.hst.node_count());
        assert_eq!(hst1.hst.edge_count(), hst2.hst.edge_count());

        for (idx1, idx2) in hst1.hst.node_indices().zip(hst2.hst.node_indices()) {
            assert_eq!(idx1, idx2);
        }

        Ok(())
    }
}
