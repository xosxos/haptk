use std::path::Path;
use std::path::PathBuf;

use rust_htslib::bcf::IndexedReader;
use rust_htslib::bcf::Read;

use crate::io::read_multiple_sample_ids;
use crate::Error;

pub mod header;

pub use header::Header;

pub fn contig_len_from_vcf(path: &Path, contig: &str) -> Result<u64, Error> {
    header::Header::try_get(path)?
        .reference_sequences()?
        .iter()
        .find(|(ctg, _len)| ctg == contig)
        .map(|(_, len)| *len)
        .ok_or_else(|| Error::NoContigLength {
            contig: contig.to_string(),
        })
}

pub fn filter_out_samples_not_in_vcf(
    path: &Path,
    wanted_samples: Vec<String>,
) -> Result<Vec<String>, Error> {
    let samples = header::Header::try_get(path)?.samples()?;

    for i in &wanted_samples {
        if !samples.contains(i) {
            tracing::warn!("Wanted sample {i} is not in the VCF");
        }
    }

    Ok(samples
        .into_iter()
        .filter(|s| wanted_samples.contains(s))
        .collect())
}

pub fn get_indexes_and_sample_ids_from_vcf(
    path: &Path,
    path_samples: &Option<Vec<PathBuf>>,
    mut wanted_samples: Option<Vec<String>>,
) -> Result<(Vec<usize>, Vec<String>), Error> {
    let vcf_samples = IndexedReader::from_path(path)?
        .header()
        .samples()
        .into_iter()
        .map(|sample| Ok(String::from_utf8_lossy(sample).to_string()))
        .collect::<Result<Vec<String>, Error>>()?;

    if wanted_samples.is_none() {
        wanted_samples = read_multiple_sample_ids(path_samples)?;
    }

    let sample_indexes = filter_samples(&vcf_samples, wanted_samples);

    if sample_indexes.is_empty() {
        return Err(Error::SamplesNotFound);
    }

    let samples: Vec<String> = sample_indexes
        .iter()
        .map(|s| &vcf_samples[*s])
        .cloned()
        .collect();

    Ok((sample_indexes, samples))
}

fn filter_samples(vcf_samples: &[String], wanted_samples: Option<Vec<String>>) -> Vec<usize> {
    if let Some(wanted_samples) = wanted_samples {
        for sample in &wanted_samples {
            if !vcf_samples.contains(sample) {
                tracing::warn!("Wanted sample {sample} is not in the VCF");
            }
        }

        vcf_samples
            .iter()
            .enumerate()
            .filter(|(_, s)| wanted_samples.contains(s))
            .map(|(i, _)| i)
            .collect()
    } else {
        (0..vcf_samples.len()).collect()
    }
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
}
