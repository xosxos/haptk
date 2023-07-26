use std::collections::HashMap;
use std::ffi::{OsStr, OsString};
use std::io;
use std::path::PathBuf;
use std::thread::{self, JoinHandle};

use bgzip::tabix::Tabix;
use color_eyre::{eyre::eyre, Result};
use csv::{QuoteStyle, Reader, ReaderBuilder, Writer, WriterBuilder};
use eyre::Context;
use rust_htslib::bcf::{header::HeaderRecord, IndexedReader, Read};
use std::sync::mpsc::Receiver;

use crate::error::HatkError::{CoordsParseError, PosParseError};

pub fn get_tsv_reader<R: io::Read>(input: R, has_headers: bool) -> Reader<R> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(has_headers)
        .flexible(false)
        .from_reader(input)
}

pub fn get_csv_reader<R: io::Read>(input: R) -> Reader<R> {
    ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .flexible(false)
        .from_reader(input)
}

pub fn get_csv_writer<W: io::Write>(output: W) -> Writer<W> {
    WriterBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .flexible(true)
        .from_writer(output)
}

pub fn get_strict_tsv_writer<W: io::Write>(output: W) -> Writer<W> {
    WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_writer(output)
}

pub fn get_vcf_writer<W: io::Write>(output: W) -> Writer<W> {
    WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .double_quote(false)
        .quote_style(QuoteStyle::Never)
        .from_writer(output)
}

pub fn get_input(filename: Option<PathBuf>) -> Result<Box<dyn io::Read>> {
    let input: Box<dyn io::Read> = match filename {
        Some(name) => match name.to_str() {
            Some("-") => Box::new(io::stdin()),
            Some(name) => {
                let r = match niffler::from_path(name) {
                    Ok(x) => x.0,
                    Err(err) => {
                        let msg = format!("failed to open \"{name}\": {err}");
                        return Err(eyre!(msg))?;
                    }
                };
                Box::new(r)
            }
            None => return Err(eyre!("Unknown I/O error")),
        },
        None => Box::new(io::stdin()),
    };
    Ok(input)
}

pub fn get_output(filename: Option<PathBuf>) -> Result<Box<dyn io::Write>> {
    let output: Box<dyn io::Write> = match filename {
        Some(name) => match name.to_str() {
            Some("-") => Box::new(io::stdout()),
            Some(name) => Box::new(
                match std::fs::File::options()
                    .create(true)
                    .write(true)
                    .truncate(true)
                    .open(name)
                {
                    Ok(x) => x,
                    Err(err) => return Err(eyre!("failed to open \"{name}\": {err}"))?,
                },
            ),
            None => return Err(eyre!("Unknown I/O error")),
        },
        None => Box::new(io::stdout()),
    };
    Ok(output)
}

pub fn open_csv_writer(name: PathBuf) -> Result<Writer<Box<dyn io::Write>>> {
    Ok(get_csv_writer(get_output(Some(name))?))
}

pub fn open_strict_tsv_writer(name: PathBuf) -> Result<Writer<Box<dyn io::Write>>> {
    Ok(get_strict_tsv_writer(get_output(Some(name))?))
}

pub fn spawn_csv_collector(
    rx: Receiver<Vec<String>>,
    header: Vec<String>,
    outname: Option<PathBuf>,
    sort: bool,
) -> JoinHandle<Result<()>> {
    thread::spawn(move || -> Result<()> {
        let output = get_output(outname)?;
        let mut output = get_csv_writer(output);
        output.write_record(header)?;

        match sort {
            true => {
                let mut lines = vec![];

                while let Ok(line) = rx.recv() {
                    lines.push(line);
                }

                lines.sort_by(|a, b| alphanumeric_sort::compare_str(&a[0], &b[0]));

                for line in lines {
                    output.write_record(&line)?;
                }
            }
            false => {
                while let Ok(line) = rx.recv() {
                    output.write_record(&line)?;
                }
            }
        }

        Ok(())
    })
}

pub fn append_ext(ext: impl AsRef<OsStr>, path: &PathBuf) -> PathBuf {
    let mut os_string: OsString = path.into();
    os_string.push(".");
    os_string.push(ext.as_ref());
    os_string.into()
}

pub fn read_tabix(path: &PathBuf) -> Result<HashMap<String, (u64, usize)>> {
    let path = append_ext("tbi", path);
    let mut file = match std::fs::File::open(&path) {
        Ok(x) => x,
        Err(err) => {
            let msg = format!("failed to open \"{}\": {err}", path.display());
            return Err(eyre!(msg))?;
        }
    };
    let tabix = Tabix::from_reader(&mut file)?;

    let mut contigs: HashMap<String, (u64, usize)> = HashMap::new();

    let names = tabix
        .names
        .iter()
        .map(|n| Ok(std::str::from_utf8(n)?))
        .collect::<Result<Vec<&str>>>()?;

    let start = tabix
        .sequences
        .iter()
        .map(|n| {
            Ok((
                *n.intervals.first().unwrap(),
                n.number_of_distinct_bin as usize,
            ))
        })
        .collect::<Result<Vec<(u64, usize)>>>()?;

    names.iter().zip(start).for_each(|(name, (start, end))| {
        let _ = contigs.insert(name.trim_matches(char::from(0)).to_string(), (start, end));
    });

    Ok(contigs)
}

pub fn parse_snp_coord(coords: &str) -> Result<(&str, u64)> {
    let mut coord_split = coords.split(':');

    match (coord_split.next(), coord_split.next()) {
        (Some(contig), Some(value)) => {
            let value = value
                .parse::<u64>()
                .wrap_err(eyre!(PosParseError((coords.into(), value.into()))))?;
            Ok((contig, value))
        }
        _ => Err(eyre!(CoordsParseError(coords.into()))),
    }
}

pub fn parse_coords(coords: &str) -> Result<(&str, Option<u64>, Option<u64>)> {
    let mut coord_split = coords.split(':');

    match (coord_split.next(), coord_split.next()) {
        (Some(contig), Some(value)) => {
            let mut pos_split = value.split('-');
            match (pos_split.next(), pos_split.next()) {
                (None, None) => Ok((contig, None, None)),
                (Some(start), None) => {
                    let start = start
                        .parse::<u64>()
                        .wrap_err(eyre!(PosParseError((coords.into(), start.into()))))?;
                    Ok((contig, Some(start), None))
                }
                (Some(start), Some(stop)) => {
                    let start = start
                        .parse::<u64>()
                        .wrap_err(eyre!(PosParseError((coords.into(), start.into()))))?;
                    let stop = stop
                        .parse::<u64>()
                        .wrap_err(eyre!(PosParseError((coords.into(), stop.into()))))?;
                    Ok((contig, Some(start), Some(stop)))
                }
                _ => panic!(),
            }
        }
        (Some(contig), None) => Ok((contig, None, None)),
        _ => Err(eyre!(CoordsParseError(coords.into()))),
    }
}

pub fn get_htslib_contig_len(path: &PathBuf, contig: &str) -> Result<u64> {
    for r in IndexedReader::from_path(path)?
        .header()
        .header_records()
        .iter()
        .filter(|r| matches!(r, HeaderRecord::Contig { .. }))
    {
        match r {
            HeaderRecord::Contig { values, .. } => {
                if let Some(id) = values.get("ID") {
                    if id == contig {
                        return Ok(values
                            .get("length")
                            .ok_or_else(|| eyre!("VCF header has no contig length for {id}"))?
                            .parse::<u64>()?);
                    }
                }
            }
            _ => {
                panic!("Header contig filtering failed for some reason. Create an issue at Github")
            }
        }
    }
    Err(eyre!(
        "Cannot get length for contig {contig}. Check the vcf header."
    ))
}

pub fn get_htslib_bcf_contigs(path: &PathBuf) -> Result<Vec<(String, i64)>> {
    IndexedReader::from_path(path)?
        .header()
        .header_records()
        .iter()
        .filter(|r| matches!(r, HeaderRecord::Contig { .. }))
        .map(|r| match r {
            HeaderRecord::Contig { values, .. } => {
                let id = values
                    .get("ID")
                    .expect("The input VCF has no ID for some contig record in the header");
                Ok((
                    id.clone(),
                    values
                        .get("length")
                        .expect(&format!("VCF header has no contig length for {id}"))
                        .parse::<i64>()?,
                ))
            }
            _ => {
                panic!("Header contig filtering failed for some reason. Create an issue at Github")
            }
        })
        .collect()
}

pub fn centromeres_hg38(chr: &str) -> (u64, u64) {
    match chr {
        "chr1" => (121700000, 125100000),
        "chr2" => (91800000, 96000000),
        "chr3" => (87800000, 94000000),
        "chr4" => (48200000, 51800000),
        "chr5" => (46100000, 51400000),
        "chr6" => (58500000, 62600000),
        "chr7" => (58100000, 62100000),
        "chr8" => (43200000, 47200000),
        "chr9" => (42200000, 45500000),
        "chr10" => (38000000, 41600000),
        "chr11" => (51000000, 55800000),
        "chr12" => (33200000, 37800000),
        "chr13" => (16500000, 18900000),
        "chr14" => (16100000, 18200000),
        "chr15" => (17500000, 20500000),
        "chr16" => (35300000, 38400000),
        "chr17" => (22700000, 27400000),
        "chr18" => (15400000, 21500000),
        "chr19" => (24200000, 28100000),
        "chr20" => (25700000, 30400000),
        "chr21" => (10900000, 13000000),
        "chr22" => (13700000, 17400000),
        "chrX" => (58100000, 63800000),
        "chrY" => (10300000, 10600000),
        _ => panic!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_snp_coord() {
        let (chr, pos) = parse_snp_coord("chr9:1920").unwrap();
        assert_eq!(chr, "chr9");
        assert_eq!(pos, 1920);
    }

    #[test]
    fn test_parse_coord() {
        let (contig, start, stop) = parse_coords("chr9").unwrap();
        assert_eq!(contig, "chr9");
        assert_eq!(start, None);
        assert_eq!(stop, None);
        let (contig, start, stop) = parse_coords("chr9:1920").unwrap();
        assert_eq!(contig, "chr9");
        assert_eq!(start, Some(1920));
        assert_eq!(stop, None);
        let (contig, start, stop) = parse_coords("chr9:1920-2500").unwrap();
        assert_eq!(contig, "chr9");
        assert_eq!(start, Some(1920));
        assert_eq!(stop, Some(2500));
        let result = parse_coords("chr9:1920--2500");
        assert!(result.is_err())
    }

    #[test]
    #[ignore]
    #[should_panic(
        expected = "called `Result::unwrap()` on an `Err` value: Position \"foo\" is not an integer in coords \"chr9:foo\"\n\nCaused by:\n    invalid digit found in string\n\nLocation:\n    src/libs/core.rs:175:38"
    )]
    fn test_parse_snp_coord_int_parser_error() {
        parse_snp_coord("chr9:foo").unwrap();
    }
}
