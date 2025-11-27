use std::fs;
use std::path::PathBuf;

use crate::error::Error;
use crate::io::get_extension;
use crate::io::read_lines;

use haptk_core::vcf;

use super::list_markers::HstMetadata;

#[doc(hidden)]
pub fn run(path: PathBuf) -> Result<(), Error> {
    for id in get_sample_names(path)? {
        println!("{id}");
    }

    Ok(())
}

pub fn get_sample_names(path: PathBuf) -> Result<Vec<String>, Error> {
    let ext = get_extension(&path)?;

    let ids = match ext.as_str() {
        "hst.gz" | "hst" => read_hst_samples(path)?,
        "vcf.gz" | "vcf" | "bcf" | "bcf.gz" => read_vcf_samples(path)?,
        "fam" => read_fam_samples(path)?,
        _ => return Err(Error::FileNotSupported { ext }),
    };

    Ok(ids)
}

pub fn read_fam_samples(path: PathBuf) -> Result<Vec<String>, Error> {
    read_lines(path)?
        .map(|line| {
            let line = line?;
            let mut split = line.split('\t');

            let id = split
                .nth(1)
                .expect("error splitting by column, is the file tab delimited?");

            Ok(id.to_string())
        })
        .collect()
}

pub fn read_vcf_samples(path: PathBuf) -> Result<Vec<String>, Error> {
    Ok(vcf::Header::try_get(&path)?.samples()?)
}

pub fn read_hst_samples(path: PathBuf) -> Result<Vec<String>, Error> {
    let file = fs::File::open(&path).map_err(|e| Error::Io { e, path })?;

    let reader = bgzip::BGZFReader::new(file)?;

    let hst: HstMetadata = serde_json::from_reader(reader)?;

    Ok(hst.metadata.samples)
}
