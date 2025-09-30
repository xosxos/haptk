mod common;

#[cfg(test)]
#[cfg(feature = "clap")]
mod compare_hst {
    use haptk::args::Selection;
    use haptk::subcommands::bhst::Hst;
    use std::path::PathBuf;
    // use std::thread::sleep;
    // use std::time::Duration;

    //--- Binary tests

    #[ignore]
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

    #[ignore]
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

    #[ignore]
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

    #[ignore]
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
        // Compare vcf to the built HST
        let args = haptk::args::StandardArgs {
            file: PathBuf::from("tests/data/test.vcf.gz"),
            output: PathBuf::from("tests/results"),
            coords: String::from("chr9:32"),
            selection: selection.clone(),
            samples: None,
            prefix: None,
            no_alt: false,
            list: None,
            only_snv: true,
        };

        let cmd = haptk::clap::SubCommand::CompareToHst {
            args: args.clone(),
            hst: PathBuf::from("tests/results/bhst.hst.gz"),
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            only_longest_leafs: false,
        };

        haptk::clap::run_cmd(cmd).unwrap();

        let args = haptk::args::StandardArgs {
            file: PathBuf::from("tests/data/test.vcf.gz"),
            output: PathBuf::from("tests/results"),
            coords: String::from("chr9:32"),
            selection,
            samples: None,
            prefix: Some(String::from("left")),
            no_alt: false,
            list: None,
            only_snv: true,
        };

        let cmd = haptk::clap::SubCommand::CompareToHst {
            args,
            hst: PathBuf::from("tests/results/uhst_left.hst.gz"),
            log_and_verbosity: crate::common::silent_verbosity(),
            threads: 8,
            only_longest_leafs: false,
        };

        haptk::clap::run_cmd(cmd).unwrap();
    }
}
