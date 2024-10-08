mod common;

#[cfg(test)]
#[cfg(feature = "clap")]
mod test_uhst {
    use super::*;

    use color_eyre::Result;

    use std::io::Write;
    use std::path::PathBuf;

    use haptk::{
        args::{Selection, StandardArgs},
        read_vcf::read_vcf_to_matrix,
        structs::Coord,
        subcommands::uhst_shard::{self, LocDirection},
    };

    #[test]
    fn matrix_to_graph() -> Result<()> {
        let args = StandardArgs {
            file: PathBuf::from(common::TEST_VCF),
            ..Default::default()
        };
        let mut vcf = read_vcf_to_matrix(&args, "chr9", 32, None, None, false).unwrap();

        let coord = Coord {
            contig: String::from("chr9"),
            reference: String::from("G"),
            alt: String::from("T"),
            pos: 32,
        };

        let g = uhst_shard::construct_uhst(&mut vcf, &LocDirection::Left, &coord, 1, false)?;
        let mut f = std::fs::File::create("tests/results/test.dot").unwrap();
        f.write_all(format!("{}", petgraph::dot::Dot::new(&g)).as_bytes())
            .unwrap();
        Ok(())
    }

    //--- Binary tests

    #[test]
    fn uhst_all() {
        run_uhst(Selection::All);

        let res = std::fs::read_to_string("tests/results/uhst_mbah.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_shared_core_haplotype.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn uhst_only_ref() {
        run_uhst(Selection::OnlyRefs);

        let res = std::fs::read_to_string("tests/results/uhst_mbah_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_shared_core_haplotype_only_refs.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn uhst_only_alt() {
        run_uhst(Selection::OnlyAlts);

        let res = std::fs::read_to_string("tests/results/uhst_mbah_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);

        let res = std::fs::read_to_string("tests/results/uhst_shared_core_haplotype_only_alts.csv")
            .unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn uhst_only_longest() {
        run_uhst(Selection::OnlyLongest);

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
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
            sharded: false,
        };
        haptk::clap::run_cmd(cmd).unwrap();
    }
}
