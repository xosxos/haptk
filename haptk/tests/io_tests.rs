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
        let wanted = haptk::io::read_sample_ids(&Some(rates)).unwrap().unwrap();
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
        let res = haptk::io::read_sample_ids(&Some(path));
        assert!(res.is_err());
    }

    #[test]
    fn read_variable_data() {
        let path = PathBuf::from("tests/data/clinical_data.csv");
        let df = haptk::io::read_variable_data_file(path).unwrap();

        let array = df["id"].utf8().unwrap();
        let vec: Vec<_> = array.into_iter().flatten().collect();
        assert_eq!(vec[0], "SAMPLE1");

        let array = df["aoo"].i64().unwrap();
        let vec: Vec<_> = array.into_iter().flatten().collect();
        assert_eq!(vec[0], 88);
        assert_eq!(vec[2], 58);

        let path = PathBuf::from("tests/data/does_not_exist.csv");
        let res = haptk::io::read_variable_data_file(path);
        assert!(res.is_err());
    }
}
