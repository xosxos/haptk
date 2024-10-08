mod common;
use common::{standard_args, TEST_HAPLOTYPE, TEST_VCF};

use std::path::PathBuf;

use haptk::{
    args::{GraphArgs, Selection, StandardArgs},
    io::read_haplotype_file,
    read_vcf::read_vcf_to_matrix,
    subcommands::{compare_to_haplotype::transform_gt_matrix_to_match_matrix, uhst_shard},
};
#[cfg(test)]
#[cfg(feature = "clap")]
mod test_compare_to_haplotype {
    use haptk::{args::SortOption, structs::MatrixSlice};

    use super::*;

    #[test]
    fn gt_matrix_to_match_matrix() {
        let file = PathBuf::from(TEST_HAPLOTYPE);
        let ht = read_haplotype_file(file).unwrap();
        let start = ht.first().unwrap();
        let end = ht.last().unwrap();

        let args = StandardArgs {
            file: PathBuf::from(TEST_VCF),
            ..Default::default()
        };
        let vcf = read_vcf_to_matrix(
            &args,
            "chr9",
            32,
            Some((Some(start.pos), Some(end.pos))),
            None,
            false,
        )
        .unwrap();

        // let ht = remove_unused_ht(&vcf, ht.clone()).unwrap();
        let vcf = transform_gt_matrix_to_match_matrix(vcf, &ht, 32).unwrap();

        for (idx, row) in vcf.matrix_axis_iter(0).enumerate() {
            println!(
                "{} {}",
                vcf.get_sample_name(idx),
                row.to_vec()
                    .iter()
                    .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"))
            );
        }

        let homozygous_variants = vec![
            vcf.matrix_slice(MatrixSlice::All, MatrixSlice::Point(13))
                .to_vec(),
            vcf.matrix_slice(MatrixSlice::All, MatrixSlice::Point(15))
                .to_vec(),
        ];
        for col in homozygous_variants {
            assert_eq!(vec![1; 28], col)
        }

        // let mut vcf = read_vcf_to_matrix(
        //     &args,
        //     "chr9",
        //     32,
        //     Some((Some(start.pos), Some(end.pos))),
        //     None,
        // )
        // .unwrap();
        // vcf.select_carriers(32, &Selection::OnlyAlts).unwrap();
        // let vcf = transform_gt_matrix_to_match_matrix(vcf, &ht, 32).unwrap();
        // for row in vcf.matrix_axis_iter(0) {
        // assert_eq!((17..46).map(|_| 1).collect::<Vec<u8>>(), row.to_vec())
        // }
    }

    // #[test]
    // fn matrix_graph_png() {
    //     let args = standard_args(Selection::All);
    //     let vcf = read_vcf_to_matrix(&args, "chr9", 32, None, None, false).unwrap();
    //     let indexes: Vec<_> = (0..vcf.samples().len()).collect();
    //     haptk::graphs::matrix_graph::matrix_graph_png(
    //         &vcf,
    //         GraphArgs::default(),
    //         false,
    //         None,
    //         &indexes,
    //     );
    // }

    // --- Integration tests

    #[test]
    fn compare_to_haplotype_all() {
        compare_to_haplotype(Selection::All, true, false);
        let res = std::fs::read_to_string("tests/results/ht_shared_segments.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn compare_to_haplotype_only_refs() {
        compare_to_haplotype(Selection::OnlyRefs, false, false);
        let res =
            std::fs::read_to_string("tests/results/ht_shared_segments_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn compare_to_haplotype_only_alts() {
        compare_to_haplotype(Selection::OnlyAlts, false, false);
        let res =
            std::fs::read_to_string("tests/results/ht_shared_segments_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn compare_to_haplotype_only_longest() {
        compare_to_haplotype(Selection::OnlyLongest, false, false);
        let res =
            std::fs::read_to_string("tests/results/ht_shared_segments_only_longest.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    fn compare_to_haplotype(selection: Selection, mark: bool, png: bool) {
        let args = standard_args(Selection::OnlyLongest);
        uhst_shard::run(args, 1, false, false).unwrap();

        let args = common::clap_standard_args(selection);

        let cmd = haptk::clap::SubCommand::CompareToHaplotype {
            args,
            haplotype: PathBuf::from("tests/results/uhst_shared_core_haplotype_only_longest.csv"),
            threads: 8,
            mark_samples: None,
            mark_shorter_alleles: mark,
            png,
            graph_args: GraphArgs::default(),
            log_and_verbosity: crate::common::silent_verbosity(),
            npy: false,
            sort_option: SortOption::Left,
        };
        assert_eq!(8, cmd.threads());
        haptk::clap::run_cmd(cmd).unwrap();
    }

    // #[test]
    // fn compare_to_haplotype_png() {
    //     compare_to_haplotype(Selection::All, false, true);
    //     let res = std::fs::read_to_string("tests/results/differences.csv").unwrap();
    //     insta::assert_yaml_snapshot!(res);
    // }

    #[test]
    fn compare_to_haplotype_is_error() {
        let args = common::clap_standard_args(Selection::Unphased);

        let cmd = haptk::clap::SubCommand::CompareToHaplotype {
            args,
            haplotype: PathBuf::from("tests/results/shared_haplotype_only_longest.csv"),
            threads: 8,
            mark_samples: None,
            mark_shorter_alleles: false,
            png: false,
            graph_args: GraphArgs::default(),
            log_and_verbosity: crate::common::silent_verbosity(),
            npy: false,
            sort_option: SortOption::Left,
        };
        assert!(haptk::clap::run_cmd(cmd).is_err());
    }
}
