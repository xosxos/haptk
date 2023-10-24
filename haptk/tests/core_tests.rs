mod common;
use common::{create_test_matrix, TEST_VCF};

#[cfg(test)]
mod core {
    use super::*;
    use std::path::PathBuf;

    use haptk::args::{Selection, StandardArgs};
    use haptk::io::read_variable_data_file;
    use haptk::read_vcf::read_vcf_to_matrix;
    use haptk::structs::Coord;

    fn create_samples() -> Vec<String> {
        vec![
            "SAMPLE1".to_string(),
            "SAMPLE2".to_string(),
            "SAMPLE3".to_string(),
            "SAMPLE4".to_string(),
            "SAMPLE5".to_string(),
            "SAMPLE6".to_string(),
            "SAMPLE7".to_string(),
            "SAMPLE8".to_string(),
            "SAMPLE9".to_string(),
            "SAMPLE10".to_string(),
            "SAMPLE11".to_string(),
            "SAMPLE12".to_string(),
            "SAMPLE13".to_string(),
            "SAMPLE14".to_string(),
        ]
    }

    #[test]
    fn read_vcf_to_phased_matrix() {
        let vcf = create_test_matrix();

        for row in vcf.matrix.rows() {
            println!(
                "{}",
                row.to_vec()
                    .iter()
                    .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"))
            );
        }

        let row = vcf.matrix.rows().into_iter().next().unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000000000000000001000100000000000000000000000000000")
        );

        let row = vcf.matrix.rows().into_iter().nth(1).unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000100000000000000010000000000000001000000000000000")
        );

        let samples = create_samples();

        let coords = (1..=63)
            .map(|p| Coord {
                pos: p,
                contig: "chr9".to_string(),
                reference: if p % 2 == 0 {
                    "G".to_string()
                } else {
                    "A".to_string()
                },
                alt: if p % 2 == 0 {
                    "T".to_string()
                } else {
                    "C".to_string()
                },
            })
            .collect::<Vec<Coord>>();

        assert_eq!(vcf.samples(), &samples);
        assert_eq!(vcf.variant_idx(), 31);
        assert_eq!(vcf.coords(), &coords);
    }

    #[test]
    fn select_rows() {
        let mut vcf = create_test_matrix();
        vcf.select_rows(vec![0, 3, 20]);

        let row = vcf.matrix.rows().into_iter().next().unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000000000000000001000100000000000000000000000000000")
        );

        assert_eq!(
            vcf.samples(),
            &vec![
                "SAMPLE1".to_string(),
                "SAMPLE2".to_string(),
                "SAMPLE11".to_string()
            ]
        );
    }

    #[test]
    fn select_columns_by_range() {
        let mut vcf = create_test_matrix();
        let (pos, idx) = (vcf.variant_idx_pos(), vcf.variant_idx());

        vcf.select_columns_by_range(10..50);

        assert_ne!(idx, vcf.variant_idx());
        assert_eq!(pos, vcf.variant_idx_pos());
    }

    #[test]
    fn select_columns() {
        let mut vcf = create_test_matrix();
        vcf.select_columns_by_idx(&mut [0, 3, 20]);
        assert_eq!(
            vcf.coords().iter().map(|c| c.pos).collect::<Vec<u64>>(),
            vec![1, 4, 21]
        );
    }

    #[test]
    fn select_alt_carriers() {
        let mut vcf = create_test_matrix();

        vcf.select_carriers(vcf.variant_idx_pos(), &Selection::OnlyAlts)
            .unwrap();

        let mut samples = create_samples();
        samples.push("SAMPLE14".to_string());

        let row = vcf.matrix.rows().into_iter().next().unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000100000000000000010000000000000001000000000000000")
        );

        let row = vcf.matrix.rows().into_iter().nth(1).unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000001000000000000000010000000000000000100000000000000")
        );

        assert_eq!(&samples, vcf.samples());
    }

    #[test]
    fn select_ref_carriers() {
        let mut vcf = create_test_matrix();
        vcf.select_carriers(vcf.variant_idx_pos(), &Selection::OnlyRefs)
            .unwrap();

        let mut samples = create_samples();
        samples.truncate(samples.len() - 1);

        let row = vcf.matrix.rows().into_iter().next().unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000000000000000001000100000000000000000000000000000")
        );

        assert_eq!(&samples, vcf.samples());
    }

    #[test]
    fn select_only_longest() {
        let mut vcf = create_test_matrix();
        vcf.select_only_longest();

        let samples = create_samples();

        let row = vcf.matrix.rows().into_iter().next().unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000100000000000000010000000000000001000000000000000")
        );

        let row = vcf.matrix.rows().into_iter().nth(1).unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000001000000000000000010000000000000000100000000000000")
        );

        assert_eq!(&samples, vcf.samples());
    }

    #[test]
    fn only_longest_indexes() {
        let vcf = create_test_matrix();
        let indexes = vcf.only_longest_indexes();
        assert_eq!(
            vec![1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27],
            indexes
        );
    }

    #[test]
    fn get_sample_names_from_phased_matrix() {
        let vcf = create_test_matrix();
        let names = vcf.get_sample_names(&[0, 3, 20]);
        assert_eq!(
            names,
            vec![
                "SAMPLE1".to_string(),
                "SAMPLE2".to_string(),
                "SAMPLE11".to_string()
            ]
        );

        let idxs = vcf.get_sample_idxs(&names).unwrap();
        assert_eq!(idxs, vec![0, 1, 2, 3, 20, 21]);

        let name = vcf.get_sample_name(0);
        assert_eq!("SAMPLE1".to_string(), name);
    }

    #[test]
    fn set_vcf_coords_and_variant_idx() {
        let mut vcf = create_test_matrix();
        vcf.set_coords(vec![]);
        assert_eq!(&Vec::<Coord>::new(), vcf.coords());

        vcf.set_variant_idx(99);
        assert_eq!(99, vcf.variant_idx());

        vcf.set_variant_idx(0);
        assert_eq!(0, vcf.variant_idx());

        let coords = vcf.coords_mut();
        coords.push(Coord {
            contig: "".to_string(),
            pos: 10,
            reference: "".to_string(),
            alt: "".to_string(),
        });

        assert_eq!(10, vcf.variant_idx_pos());
        assert_eq!(1, vcf.ncoords());
        assert_eq!(14, vcf.nsamples());
        assert_eq!(14 * 2, vcf.nrows());
    }

    #[test]
    fn nearest_idx_by_pos() {
        let mut vcf = create_test_matrix();
        let coords = vcf.coords_mut();

        coords.iter_mut().for_each(|c| c.pos = c.pos.pow(2));

        assert_eq!(4, vcf.get_nearest_idx_by_pos(26));
        assert_eq!(4, vcf.get_nearest_idx_by_pos(23));
        assert_eq!(5, vcf.get_nearest_idx_by_pos(31));
        assert_eq!(vcf.ncoords(), vcf.get_nearest_idx_by_pos(9999999));

        let coord = Coord {
            contig: "chr9".to_string(),
            pos: 25,
            reference: "A".to_string(),
            alt: "C".to_string(),
        };

        assert_eq!(4, vcf.idx_by_coord(&coord).unwrap());
    }

    // #[test]
    // fn vcf_info_filtering() {
    //     let args = StandardArgs {
    //         file: PathBuf::from(TEST_VCF),
    //         info_limit: Some(0.9),
    //         ..Default::default()
    //     };
    //     let vcf = read_vcf_to_matrix(&args, "chr9", 32, None, None).unwrap();

    //     assert_eq!(vcf.ncoords(), 31);
    // }

    // #[test]
    // fn read_vcf_to_phased_matrix_filter_info() {
    //     let args = standardargs {
    //         file: pathbuf::from(test_vcf),
    //         info_limit: some(0.9),
    //         ..default::default()
    //     };
    //     let vcf = read_vcf_to_matrix(&args, "chr9", 32, none, none).unwrap();

    //     assert_eq!(vcf.ncoords(), 31);
    // }

    #[test]
    fn read_vcf_to_phased_matrix_by_coords() {
        let vcf1 = create_test_matrix();
        let args = StandardArgs {
            file: PathBuf::from(TEST_VCF),
            ..Default::default()
        };
        let vcf2 = read_vcf_to_matrix(&args, "chr9", 16, Some((1, 32)), None).unwrap();

        assert_eq!(
            &vcf1
                .coords()
                .iter()
                .take(32)
                .cloned()
                .collect::<Vec<Coord>>(),
            vcf2.coords()
        );
    }

    #[test]
    fn error_on_not_normalized_vcf() {
        let args = StandardArgs {
            file: PathBuf::from("tests/data/test_not_normalized.vcf.gz"),
            ..Default::default()
        };

        let res = read_vcf_to_matrix(&args, "chr9", 32, None, None);
        assert!(res.is_err());
    }

    #[test]
    fn variable_data() {
        let path = PathBuf::from("tests/data/clinical_data.csv");
        let df = read_variable_data_file(path).unwrap();

        let args = StandardArgs {
            file: PathBuf::from(TEST_VCF),
            ..Default::default()
        };
        let mut vcf = read_vcf_to_matrix(&args, "chr9", 32, None, None).unwrap();
        vcf.set_variable_data(df).unwrap();

        // Make sure NA is calculated right in the average
        let means = vcf
            .get_variable_data_mean(&[0, 2], &["aoo".to_string(), "dur".to_string()])
            .unwrap()
            .unwrap();
        assert_eq!(68.50, means[0]);
        assert_eq!(6.0, means[1]);

        // Request variable data on a sample that has no variable data
        let result = vcf.get_variable_data_mean(&[27], &["aoo".to_string()]);
        println!("{result:?}");
        assert!(result.is_err());

        let vecs = vcf
            .get_variable_data_vecs(&[0, 2], &["aoo".to_string(), "dur".to_string()])
            .unwrap()
            .unwrap();
        assert_eq!(vec![88.0, 49.0], vecs[0]);
        assert_eq!(vec![2.0, 10.0], vecs[1]);
    }
}
