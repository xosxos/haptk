use std::path::PathBuf;

use color_eyre::{
    eyre::{ensure, eyre, Context},
    Result,
};
use rust_htslib::faidx::Reader;
// use serde::{Deserialize, Serialize};

use crate::error::HatkError::PosParseError;
use crate::{io::open_csv_writer, structs::HapVariant};

pub fn parse_seq_name(coords: &str) -> Result<(&str, u64, u64, &str)> {
    let error_msg = format!("Error parsing {coords:?}, make sure the sequence naming is as follows: [contig]:[start]-[stop];[annotation]");

    let mut coord_split = coords.split(':');
    let (contig, rest) = (coord_split.next(), coord_split.next());
    let contig = contig.expect(&error_msg);
    let rest = rest.expect(&error_msg);

    let mut rest_split = rest.split('-');
    let (start, rest) = (rest_split.next(), rest_split.next());
    let start = start.expect(&error_msg);
    let rest = rest.expect(&error_msg);

    let (stop, ann) = rest.split_once('_').expect(&error_msg);

    let start = start
        .parse::<u64>()
        .wrap_err(eyre!(PosParseError((coords.into(), start.into()))))?;

    let stop = stop
        .parse::<u64>()
        .wrap_err(eyre!(PosParseError((coords.into(), stop.into()))))?;

    return Ok((contig, start, stop, ann));
}

pub fn run(path: PathBuf, seq_names: Vec<String>, output: PathBuf) -> Result<()> {
    ensure!(
        !seq_names.is_empty(),
        "Please enter one or more sequence names from the fasta using the --seq-name parameter"
    );

    let fasta_reader =
        Reader::from_path(&path).wrap_err(eyre!("Could not read the path: {:?}", path))?;

    let mut fasta: Vec<HapVariant> = vec![];

    for seq_name in seq_names {
        let sequence = fasta_reader.fetch_seq_string(&seq_name, 0, 20000)?;
        let (contig, start, _stop, annotation) = parse_seq_name(&seq_name)?;
        tracing::info!(
            " Reading sequence: {seq_name:?}, length: {}",
            sequence.len()
        );
        for (n, gt) in sequence.chars().enumerate() {
            let hap_variant = HapVariant {
                contig: contig.to_string(),
                pos: start + n as u64,
                reference: gt.to_string(),
                alt: "-".to_string(),
                annotation: Some(annotation.to_string()),
                gt: 0,
            };
            fasta.push(hap_variant);
        }
    }

    let mut writer = open_csv_writer(output)?;
    println!("contig,pos,ref,alt,gt,ann");
    fasta
        .into_iter()
        .try_for_each(|hap_variant| -> Result<()> {
            writer.write_record(vec![
                hap_variant.contig,
                hap_variant.pos.to_string(),
                hap_variant.reference,
                hap_variant.alt,
                hap_variant.gt.to_string(),
                format!("{}", hap_variant.annotation.unwrap()),
            ])?;
            Ok(())
        })?;

    tracing::debug!("Finished writing the fasta to haplotype");

    Ok(())
}
