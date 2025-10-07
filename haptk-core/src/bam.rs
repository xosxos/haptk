mod header;
mod reader;
mod record;
mod writer;

pub use header::Header;
pub use reader::Reader;
pub use record::Record;
pub use writer::Writer;

use crate::Error;
use rayon::prelude::*;
use std::path::PathBuf;

pub fn zip_bam_paths_to_sample_name(paths: Vec<PathBuf>) -> Result<Vec<(String, PathBuf)>, Error> {
    let mut names_and_paths: Vec<(String, PathBuf)> = paths
        .into_par_iter()
        // .into_iter()
        .enumerate()
        .map(|(i, path)| {
            let header = Header::try_get(&path)?;

            let samples = header.read_group_samples()?;

            if samples.is_empty() {
                tracing::warn!(
                    "No sample name in the .bam header for {path:?}. Naming the sample unnamed_{i}"
                );
                let sample_name = format!("unnamed_{i}");

                return Ok((sample_name, path));
            }

            if samples.len() > 1 {
                tracing::warn!(
                    "Multiple read groups, selecting the first one: {}",
                    &samples[0]
                )
            }

            Ok((samples[0].to_string(), path))
        })
        .collect::<Result<Vec<(String, PathBuf)>, Error>>()?;

    names_and_paths.sort_by(|a, b| a.0.cmp(&b.0));

    let name_vec: Vec<String> = names_and_paths
        .iter()
        .map(|(name, _path)| name.clone())
        .collect();

    names_and_paths.iter_mut().enumerate().for_each(|(i, v)| {
        if name_vec[i + 1..].contains(&v.0) {
            v.0.push_str(&i.to_string());
        }
    });

    Ok(names_and_paths)
}
