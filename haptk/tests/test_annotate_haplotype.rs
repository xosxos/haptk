mod common;

// use std::path::PathBuf;

// #[test]
// #[cfg(feature = "clap")]
// fn annotate_haplotype() {
//     use haptk::clap::LogAndVerbosity;

//     let cmd = haptk::clap::SubCommand::AnnotateHaplotype {
//         file: PathBuf::from(common::TEST_HAPLOTYPE),
//         ann: PathBuf::from("tests/data/test.gtf.gz"),
//         outdir: PathBuf::from("tests/results/annotated_haplotype.csv"),
//         log_and_verbosity: LogAndVerbosity {
//             verbosity: 1,
//             log_file: None,
//         },
//     };
//     haptk::clap::run_cmd(cmd).unwrap();

//     let res = std::fs::read_to_string("tests/results/annotated_haplotype.csv").unwrap();
//     insta::assert_yaml_snapshot!(res);
// }
