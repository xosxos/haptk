use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::args::{Selection, StandardArgs};

pub fn push_to_output(args: &StandardArgs, output: &mut PathBuf, name: &str, suffix: &str) {
    if let Some(prefix) = &args.prefix {
        match args.selection {
            Selection::All => output.push(format!("{prefix}_{name}.{suffix}")),
            Selection::OnlyAlts => output.push(format!("{prefix}_{name}_only_alts.{suffix}")),
            Selection::OnlyRefs => output.push(format!("{prefix}_{name}_only_refs.{suffix}")),
            Selection::OnlyLongest => output.push(format!("{prefix}_{name}_only_longest.{suffix}")),
            Selection::Unphased => output.push(format!("{prefix}_{name}_rwc.{suffix}")),
            Selection::Haploid => output.push(format!("{prefix}_{name}_haploid.{suffix}")),
        }
    } else {
        match args.selection {
            Selection::All => output.push(format!("{name}.{suffix}")),
            Selection::OnlyAlts => output.push(format!("{name}_only_alts.{suffix}")),
            Selection::OnlyRefs => output.push(format!("{name}_only_refs.{suffix}")),
            Selection::OnlyLongest => output.push(format!("{name}_only_longest.{suffix}")),
            Selection::Unphased => output.push(format!("{name}_rwc.{suffix}")),
            Selection::Haploid => output.push(format!("{name}_haploid.{suffix}")),
        }
    }
}

pub fn filter_samples(samples: &Vec<String>, wanted: Option<Vec<String>>) -> Vec<usize> {
    if let Some(wanted) = wanted {
        for i in &wanted {
            if !samples.contains(i) {
                tracing::warn!("Wanted sample {i} is not in the VCF");
            }
        }

        samples
            .iter()
            .enumerate()
            .filter(|(_, s)| wanted.contains(s))
            .map(|(i, _)| i)
            .collect()
    } else {
        (0..samples.len()).collect()
    }
}

pub fn _to_centimorgans(rates: &BTreeMap<u64, f32>, start: u64, end: u64) -> f32 {
    let start_cm = match rates.range(start..).next() {
        Some(cm) => cm,
        None => rates.range(..start).next_back().unwrap(),
    };

    let end_cm = match rates.range(end..).next() {
        Some(cm) => cm,
        None => rates.range(..end).next_back().unwrap(),
    };
    (start_cm.1 - end_cm.1).abs()
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn sample_filtering() {
        let samples = vec!["foo".to_string(), "fii".to_string()];
        let wanted = vec!["foo".to_string(), "fii".to_string()];
        let sample_indexes = filter_samples(&samples, Some(wanted));
        assert_eq!(sample_indexes, vec![0,1]);

        let samples = vec!["foo".to_string(), "fii".to_string()];
        let wanted = vec!["fee".to_string(), "foo".to_string(), "fii".to_string()];
        let sample_indexes = filter_samples(&samples, Some(wanted));
        assert_eq!(sample_indexes, vec![0,1]);

        let samples = vec!["faa".to_string(), "foo".to_string(), "fii".to_string()];
        let wanted = vec!["fee".to_string(), "foo".to_string(), "fii".to_string()];
        let sample_indexes = filter_samples(&samples, Some(wanted));
        assert_eq!(sample_indexes, vec![1,2]);

        let sample_indexes = filter_samples(&samples, None);
        assert_eq!(sample_indexes, vec![0,1,2]);
    }

    #[test]
    fn test_push_to_output() {
        let mut output = std::path::PathBuf::new();
        let args = crate::args::StandardArgs::default();
        push_to_output(&args, &mut output, "picture", "png");
        assert_eq!(output, std::path::PathBuf::from("picture.png"));

        let mut output = std::path::PathBuf::from("./foo");
        let args = crate::args::StandardArgs::default();
        push_to_output(&args, &mut output, "picture", "png");
        assert_eq!(output, std::path::PathBuf::from("./foo/picture.png"));

        let mut output = std::path::PathBuf::from("./foo");
        let args = crate::args::StandardArgs { 
            prefix: Some("nice".to_string()), 
            ..Default::default()
        };
        push_to_output(&args, &mut output, "picture", "png");
        assert_eq!(output, std::path::PathBuf::from("./foo/nice_picture.png"));
    }
}
