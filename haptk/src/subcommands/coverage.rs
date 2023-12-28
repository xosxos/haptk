use std::path::PathBuf;

use color_eyre::{eyre::eyre, Result};
use rayon::prelude::*;
use rust_htslib::bcf::{IndexedReader, Read};
use termion::color;

use crate::io::get_htslib_bcf_contigs;

#[doc(hidden)]
#[tracing::instrument]
pub fn run(path: PathBuf, bp_per_snp: i64, npipes: i64) -> Result<()> {
    let contigs = get_htslib_bcf_contigs(&path)?;

    let lines = contigs
        .par_iter()
        .map_with(path, |path: &mut PathBuf, (contig, contig_size)| {
            // Fetch readers
            let mut reader = IndexedReader::from_path(path)?;
            let header_view = reader.header();

            let rid = header_view.name2rid(contig.as_bytes())?;
            reader.fetch(rid, 0, None)?;

            let records = reader.records();

            // Program logic
            let counts_per_window = counts_per_window(records, *contig_size, npipes)?;
            let line = create_line(counts_per_window, contig, bp_per_snp);

            Ok(line)
        })
        .collect::<Result<Vec<String>>>()?;

    for line in lines {
        println!("{line}");
    }

    Ok(())
}

fn counts_per_window<R: Read>(
    records: rust_htslib::bcf::Records<R>,
    contig_size: i64,
    npipes: i64,
) -> Result<Vec<i64>> {
    let window_size = contig_size / npipes;

    // Count hits per window in a single iteration
    let mut window = 1;
    let mut counter = 1;
    let mut counts_per_window = vec![];
    let mut record_counter = 0;
    for record in records {
        record_counter += 1;
        let record = record?;

        if record.pos() < window * window_size {
            counter += 1;
        } else {
            counts_per_window.push(window_size / counter);
            window += 1;
            counter = 1;
        }
    }

    if counts_per_window.len() != npipes as usize {
        counts_per_window.push(window_size / counter);
    }
    while counts_per_window.len() != npipes as usize {
        counts_per_window.push(window_size)
    }

    if record_counter < 99 {
        Err(eyre!("Contig has under 100 records: {record_counter}"))
    } else {
        Ok(counts_per_window)
    }
}

fn create_line(window_counts: Vec<i64>, contig: &str, bp_per_snp: i64) -> String {
    let start = format!("\n{}{contig}\t", color::Fg(color::Reset));

    window_counts.iter().fold(start, |acc, count| {
        if count <= &bp_per_snp {
            format!("{acc}{}|", color::Fg(color::AnsiValue::rgb(0, 5, 1)))
        } else {
            format!("{acc}{}-", color::Fg(color::Magenta))
        }
    })
}
