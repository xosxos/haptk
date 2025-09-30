use std::path::PathBuf;

use color_eyre::{
    eyre::{ensure, WrapErr},
    Result,
};
use rust_htslib::faidx::Reader;

use crate::error::Error;
use crate::io::open_csv_writer;

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

    let start = start.parse::<u64>().wrap_err(Error::PosParse {
        coord: coords.into(),
        value: start.into(),
    })?;

    let stop = stop.parse::<u64>().wrap_err(Error::PosParse {
        coord: coords.into(),
        value: stop.into(),
    })?;

    Ok((contig, start, stop, ann))
}

pub fn run(path: PathBuf, seq_names: Vec<String>, output: PathBuf) -> Result<()> {
    ensure!(
        !seq_names.is_empty(),
        "Please enter one or more sequence names from the fasta using the --seq-name parameter"
    );

    let fasta_reader = Reader::from_path(&path).wrap_err(Error::Io { path })?;

    let mut fasta: Vec<Vec<String>> = vec![];

    for seq_name in seq_names {
        let sequence = fasta_reader.fetch_seq_string(&seq_name, 0, 20000)?;
        let (contig, start, _stop, annotation) = parse_seq_name(&seq_name)?;
        tracing::info!(
            " Reading sequence: {seq_name:?}, length: {}",
            sequence.len()
        );
        for (n, gt) in sequence.chars().enumerate() {
            let csv_row = vec![
                contig.to_string(),
                (start + n as u64).to_string(),
                gt.to_string(),
                "-".to_string(),
                format!("0"),
                annotation.to_string(),
            ];
            fasta.push(csv_row);
        }
    }

    let mut writer = open_csv_writer(output)?;
    println!("contig,pos,ref,alt,gt,ann");
    fasta.into_iter().try_for_each(|row| -> Result<()> {
        writer.write_record(row)?;
        Ok(())
    })?;

    tracing::debug!("Finished writing the fasta to haplotype");

    Ok(())
}
