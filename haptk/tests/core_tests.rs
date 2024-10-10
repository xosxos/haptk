mod common;
use common::{create_test_matrix, TEST_VCF};

#[cfg(test)]
mod core {
    use super::*;
    use std::collections::BTreeSet;
    use std::path::PathBuf;

    use color_eyre::Result;

    use haptk::args::StandardArgs;
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
    fn read_vcf_to_phased_matrix() -> Result<()> {
        let vcf = create_test_matrix()?;

        for row in vcf.matrix_axis_iter(0) {
            println!(
                "{}",
                row.to_vec()
                    .iter()
                    .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"))
            );
        }

        let row = vcf.matrix_axis_iter(0).next().unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000000000000000001000100000000000000000000000000000")
        );

        let row = vcf.matrix_axis_iter(0).nth(1).unwrap().to_vec();
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
            .collect::<BTreeSet<Coord>>();

        assert_eq!(vcf.samples(), &samples);
        assert_eq!(vcf.variant_idx(), 31);
        assert_eq!(vcf.coords(), &coords);
        Ok(())
    }

    #[test]
    fn select_rows() -> Result<()> {
        let mut vcf = create_test_matrix()?;
        vcf.select_rows(vec![0, 3, 20]);

        let row = vcf.matrix_axis_iter(0).next().unwrap().to_vec();
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
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_nearest_coord_by_pos() {
        assert_eq!(true, false)
    }

    #[test]
    fn select_columns_by_range() -> Result<()> {
        let mut vcf = create_test_matrix()?;
        let (pos, idx) = (vcf.variant_idx_pos(), vcf.variant_idx());

        println!("coords len pre {}", vcf.ncoords());
        vcf.select_columns_by_range_idx(10..50);
        println!("coords len post {}", vcf.ncoords());
        println!("coords len post {:?}", vcf.coords());

        println!("variant_idx pre {idx} and post {}", vcf.variant_idx());
        assert_ne!(idx, vcf.variant_idx());
        assert_eq!(pos, vcf.variant_idx_pos());
        Ok(())
    }

    // #[test]
    // fn select_alt_carriers() {
    //     let mut vcf = create_test_matrix()?;

    //     vcf.select_carriers(vcf.variant_idx_pos(), &Selection::OnlyAlts)
    //         .unwrap();

    //     let mut samples = create_samples();
    //     samples.push("SAMPLE14".to_string());

    //     let row = vcf.matrix_axis_iter(0).next().unwrap().to_vec();
    //     let row = row
    //         .to_vec()
    //         .iter()
    //         .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
    //     assert_eq!(
    //         row,
    //         String::from("000000000000000100000000000000010000000000000001000000000000000")
    //     );

    //     let row = vcf.matrix_axis_iter(0).nth(1).unwrap().to_vec();
    //     let row = row
    //         .to_vec()
    //         .iter()
    //         .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
    //     assert_eq!(
    //         row,
    //         String::from("000000000000001000000000000000010000000000000000100000000000000")
    //     );

    //     assert_eq!(&samples, vcf.samples());
    // }

    // #[test]
    // fn select_ref_carriers() {
    //     let mut vcf = create_test_matrix()?;
    //     vcf.select_carriers(vcf.variant_idx_pos(), &Selection::OnlyRefs)
    //         .unwrap();

    //     let mut samples = create_samples();
    //     samples.truncate(samples.len() - 1);

    //     let row = vcf.matrix_axis_iter(0).next().unwrap().to_vec();
    //     let row = row
    //         .to_vec()
    //         .iter()
    //         .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
    //     assert_eq!(
    //         row,
    //         String::from("000000000000000000000000000001000100000000000000000000000000000")
    //     );

    //     assert_eq!(&samples, vcf.samples());
    // }

    #[test]
    fn select_only_longest() -> Result<()> {
        let mut vcf = create_test_matrix()?;
        vcf.select_only_longest()?;

        let samples = create_samples();

        let row = vcf.matrix_axis_iter(0).next().unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000000100000000000000010000000000000001000000000000000")
        );

        let row = vcf.matrix_axis_iter(0).nth(1).unwrap().to_vec();
        let row = row
            .to_vec()
            .iter()
            .fold(String::new(), |cur, nxt| format!("{cur}{nxt}"));
        assert_eq!(
            row,
            String::from("000000000000001000000000000000010000000000000000100000000000000")
        );

        assert_eq!(&samples, vcf.samples());
        Ok(())
    }

    #[test]
    fn only_longest_indexes() -> Result<()> {
        let mut vcf = create_test_matrix()?;
        let indexes = vcf.only_longest_indexes()?;
        assert_eq!(
            vec![1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27],
            indexes
        );
        Ok(())
    }

    #[test]
    fn get_sample_names_from_phased_matrix() -> Result<()> {
        let vcf = create_test_matrix()?;
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
        Ok(())
    }

    #[test]
    fn set_vcf_coords_and_variant_idx() -> Result<()> {
        let mut vcf = create_test_matrix()?;
        vcf.set_coords(BTreeSet::new());
        assert_eq!(&BTreeSet::<Coord>::new(), vcf.coords());

        let coords = vcf.coords_mut();
        coords.insert(Coord {
            contig: "chr9".to_string(),
            pos: 10,
            reference: "T".to_string(),
            alt: "C".to_string(),
        });
        coords.insert(Coord {
            contig: "chr9".to_string(),
            pos: 11,
            reference: "A".to_string(),
            alt: "T".to_string(),
        });

        vcf.set_variant_idx(1);
        assert_eq!(1, vcf.variant_idx());

        assert_eq!(11, vcf.variant_idx_pos());
        assert_eq!(2, vcf.ncoords());
        assert_eq!(14, vcf.nsamples());
        assert_eq!(14 * 2, vcf.nhaplotypes());
        Ok(())
    }

    #[test]
    fn nearest_idx_by_pos() -> Result<()> {
        let mut vcf = create_test_matrix()?;
        let coords = vcf.coords();

        let mut vec = Vec::from_iter(coords.iter().cloned());
        vec.iter_mut().for_each(|c| c.pos = c.pos.pow(2));
        vcf.set_coords(vec.into_iter().collect());

        let one = vcf.get_nearest_coord_by_pos(26);
        let two = vcf.get_nearest_coord_by_pos(23);
        let three = vcf.get_nearest_coord_by_pos(31);
        let four = vcf.get_nearest_coord_by_pos(9999999);
        // println!("{one:?}");

        assert_eq!(4, vcf.get_coord_idx(one));
        println!("{one:?}");
        assert_eq!(4, vcf.get_coord_idx(two));
        println!("{two:?}");
        assert_eq!(5, vcf.get_coord_idx(three));
        assert_eq!(vcf.ncoords() - 1, vcf.get_coord_idx(four));

        let coord = Coord {
            contig: "chr9".to_string(),
            pos: 25,
            reference: "A".to_string(),
            alt: "C".to_string(),
        };

        let idx = vcf.coords().iter().position(|c| c == &coord).unwrap();
        assert_eq!(4, idx);
        Ok(())
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
    fn read_vcf_to_phased_matrix_by_coords() -> Result<()> {
        let vcf1 = create_test_matrix()?;
        let args = StandardArgs {
            file: PathBuf::from(TEST_VCF),
            ..Default::default()
        };
        let vcf2 =
            read_vcf_to_matrix(&args, "chr9", 16, Some((Some(1), Some(32))), None, None).unwrap();

        assert_eq!(
            &vcf1
                .coords()
                .iter()
                .take(32)
                .cloned()
                .collect::<BTreeSet<Coord>>(),
            vcf2.coords()
        );
        Ok(())
    }

    #[test]
    fn error_on_not_normalized_vcf() -> Result<()> {
        let args = StandardArgs {
            file: PathBuf::from("tests/data/test_not_normalized.vcf.gz"),
            ..Default::default()
        };

        let res = read_vcf_to_matrix(&args, "chr9", 32, None, None, None);
        assert!(res.is_err());
        Ok(())
    }
}
