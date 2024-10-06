mod common;
use common::{create_test_matrix, standard_args};

#[cfg(test)]
#[cfg(feature = "clap")]
mod test_check_for_haplotype {
    use super::*;
    use std::fs::read_to_string;
    use std::path::PathBuf;

    use common::clap_standard_args;

    use haptk::{args::Selection, subcommands::uhst_shard};
    use haptk::{structs::HapVariant, subcommands::check_for_haplotype};

    #[rustfmt::skip]
    #[test]
    fn test_check_for_haplotype() {
        let vcf = create_test_matrix();

        let ht = vec![
            HapVariant{ pos: 32, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 0 },
            HapVariant{ pos: 33, contig: "chr9".into(), reference: "A".into(), alt: "C".into(), gt: 0 }
        ];

        let matches = check_for_haplotype::identical_haplotype_count(&vcf, &ht);
        assert_eq!(vec![0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24], matches);

        let ht = vec![
            HapVariant{ pos: 32, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 0 },
            HapVariant{ pos: 33, contig: "chr9".into(), reference: "A".into(), alt: "C".into(), gt: 0 },
            HapVariant{ pos: 34, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 1 },
        ];

        let matches = check_for_haplotype::identical_haplotype_count(&vcf, &ht);
        assert_eq!(matches.len(), 1);

        let ht = vec![
            HapVariant{ pos: 32, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 0 },
            HapVariant{ pos: 34, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 1 },
        ];

        let matches = check_for_haplotype::identical_haplotype_count(&vcf, &ht);
        assert_eq!(matches.len(), 1);

        let ht = vec![
            HapVariant{ pos: 32, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 0 },
            HapVariant{ pos: 34, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 1 },
            // This does not exist in the test matrix so it will not be taken into account
            HapVariant{ pos: 35, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 1 }
        ];

        let matches = check_for_haplotype::identical_haplotype_count(&vcf, &ht);
        assert_eq!(matches.len(), 1);

        let ht = vec![
            HapVariant{ pos: 32, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 0 },
            HapVariant{ pos: 34, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 1 },
            HapVariant{ pos: 35, contig: "chr9".into(), reference: "A".into(), alt: "C".into(), gt: 1 }
        ];

        let matches = check_for_haplotype::identical_haplotype_count(&vcf, &ht);
        assert_eq!(matches.len(), 0);

        let ht = vec![
            HapVariant{ pos: 9999, contig: "chr9".into(), reference: "G".into(), alt: "T".into(), gt: 0 },
        ];

        let matches = check_for_haplotype::identical_haplotype_count(&vcf, &ht);

        assert!(matches.is_empty());
    }

    #[test]
    fn check_for_haplotype_all() {
        check_for_haplotype(Selection::All);

        let res = read_to_string("tests/results/haplotype_check.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn check_for_haplotype_only_refs() {
        check_for_haplotype(Selection::OnlyRefs);

        let res = read_to_string("tests/results/haplotype_check_only_refs.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn check_for_haplotype_only_alts() {
        check_for_haplotype(Selection::OnlyAlts);

        let res = read_to_string("tests/results/haplotype_check_only_alts.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    #[test]
    fn check_for_haplotype_only_longest() {
        check_for_haplotype(Selection::OnlyLongest);

        let res = read_to_string("tests/results/haplotype_check_only_longest.csv").unwrap();
        insta::assert_yaml_snapshot!(res);
    }

    fn check_for_haplotype(selection: Selection) {
        let args = standard_args(Selection::OnlyLongest);
        uhst_shard::run(args, 1, false).unwrap();

        let args = clap_standard_args(selection);
        let cmd = haptk::clap::SubCommand::CheckForHaplotype {
            args,
            haplotype: PathBuf::from("tests/results/uhst_shared_core_haplotype_only_longest.csv"),
            log_and_verbosity: haptk::clap::LogAndVerbosity {
                verbosity: 1,
                log_file: None,
                silent: false,
            },
        };

        haptk::clap::run_cmd(cmd).unwrap();
    }
}
