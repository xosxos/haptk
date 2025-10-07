mod common;
use common::{REC_RATES, TEST_HAPLOTYPE};

#[cfg(test)]
mod io {
    use super::*;

    use std::path::PathBuf;

    #[test]
    fn read_haplotype_file() {
        let file = PathBuf::from(TEST_HAPLOTYPE);
        let ht = haptk::io::read_haplotype_file(file).unwrap();
        let pos: Vec<u64> = (18..=46).collect();

        assert_eq!(pos, ht.iter().map(|h| h.pos).collect::<Vec<u64>>());

        let file = PathBuf::from("tests/data/test_haplotype_unsorted.csv");
        let result = haptk::io::read_haplotype_file(file);
        assert!(result.is_err());
    }

    #[test]
    fn read_recombination_rates() {
        let rates = PathBuf::from(REC_RATES);
        let btree = haptk::io::read_recombination_file(rates).unwrap();
        let pos: Vec<_> = (1..64).collect();
        for (i, (p, _)) in btree.iter().enumerate() {
            assert_eq!(&pos[i], p);
        }
        let rates = PathBuf::from("tests/data/recombination_rates_with_header.tsv");
        let res = haptk::io::read_recombination_file(rates);
        assert!(res.is_err());
    }

    #[test]
    fn read_samples_file() {
        let rates = PathBuf::from("tests/data/samples.ids");
        let wanted = haptk::io::read_multiple_sample_ids(&Some(vec![rates]))
            .unwrap()
            .unwrap();
        let vec = vec![
            "SAMPLE1",
            "SAMPLE2",
            "SAMPLE3",
            "SAMPLE10",
            "#*@#_)i%*#",
            "foo_foosa_23e41_das",
            "!!!!!!!",
            "#######",
        ];
        assert_eq!(vec, wanted);

        let path = PathBuf::from("tests/data/does_not_exist.ids");
        let res = haptk::io::read_multiple_sample_ids(&Some(vec![path]));
        assert!(res.is_err());
    }
}
