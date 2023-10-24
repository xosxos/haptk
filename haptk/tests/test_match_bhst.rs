mod common;

#[cfg(test)]
#[cfg(feature = "clap")]
mod test_bhst {
    use haptk::args::Selection;
    use haptk::subcommands::bhst::Hst;
    use std::path::PathBuf;

    //--- Binary tests

    #[test]
    fn compare_to_bhst_all() {
        run_compare_to_hst(Selection::All);

        let path = "tests/results/match_hst.hst.gz";

        let file = std::fs::File::open(path).expect("Error opening {path:?}");
        let reader = bgzip::BGZFReader::new(file).unwrap();
        let hst: Hst = serde_json::from_reader(reader).unwrap();
        let hst_string = serde_json::to_string(&hst).unwrap();
        insta::assert_yaml_snapshot!(hst_string);
    }

    #[test]
    fn compare_to_bhst_only_ref() {
        run_compare_to_hst(Selection::OnlyRefs);

        let path = "tests/results/match_hst_only_refs.hst.gz";

        let file = std::fs::File::open(path).expect("Error opening {path:?}");
        let reader = bgzip::BGZFReader::new(file).unwrap();
        let hst: Hst = serde_json::from_reader(reader).unwrap();
        let hst_string = serde_json::to_string(&hst).unwrap();
        insta::assert_yaml_snapshot!(hst_string);
    }

    #[test]
    fn compare_to_bhst_only_alt() {
        run_compare_to_hst(Selection::OnlyAlts);

        let path = "tests/results/match_hst_only_alts.hst.gz";

        let file = std::fs::File::open(path).expect("Error opening {path:?}");
        let reader = bgzip::BGZFReader::new(file).unwrap();
        let hst: Hst = serde_json::from_reader(reader).unwrap();
        let hst_string = serde_json::to_string(&hst).unwrap();
        insta::assert_yaml_snapshot!(hst_string);
    }

    #[test]
    fn compare_to_bhst_only_longest() {
        run_compare_to_hst(Selection::OnlyLongest);

        let path = "tests/results/match_hst_only_longest.hst.gz";

        let file = std::fs::File::open(path).expect("Error opening {path:?}");
        let reader = bgzip::BGZFReader::new(file).unwrap();
        let hst: Hst = serde_json::from_reader(reader).unwrap();
        let hst_string = serde_json::to_string(&hst).unwrap();
        insta::assert_yaml_snapshot!(hst_string);
    }

    fn run_compare_to_hst(selection: Selection) {
        // Build match HST
        let args = haptk::clap::ClapStandardArgs {
            file: PathBuf::from("tests/data/test.vcf.gz"),
            outdir: PathBuf::from("tests/results"),
            coords: String::from("chr9:32"),
            select: selection.clone().into(),
            samples: None,
            // info_limit: None,
            prefix: None,
        };

        let cmd = haptk::clap::SubCommand::Bhst {
            args,
            graph_args: haptk::clap::ClapGraphArgs::default(),
            mark_samples: None,
            variable_data: None,
            variables: None,
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            min_size: 1,
            publish: false,
        };

        // Build match HST
        haptk::clap::run_cmd(cmd).unwrap();

        // Compare vcf to the built HST
        let args = haptk::clap::ClapStandardArgs {
            file: PathBuf::from("tests/data/test_match.vcf.gz"),
            outdir: PathBuf::from("tests/results"),
            coords: String::from("chr9:32"),
            select: selection.into(),
            samples: None,
            // info_limit: None,
            prefix: None,
        };

        let cmd = haptk::clap::SubCommand::CompareToHst {
            args,
            hst: PathBuf::from("tests/results/bhst_only_longest.hst.gz"),
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
        };

        haptk::clap::run_cmd(cmd).unwrap();
    }
}
