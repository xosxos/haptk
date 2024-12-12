use std::collections::BTreeMap;

use std::collections::HashMap;
use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::sync::mpsc::Receiver;
use std::thread::{self, JoinHandle};

use bgzip::tabix::Tabix;
use color_eyre::eyre::{eyre, OptionExt, WrapErr};
use color_eyre::Result;
use csv::{QuoteStyle, Reader, ReaderBuilder, Writer, WriterBuilder};
use rust_htslib::bcf::{header::HeaderRecord, IndexedReader, Read};

use crate::args::{Selection, StandardArgs};
use crate::read_vcf::get_reader;
use crate::structs::Coord;
use crate::structs::HapVariant;
use crate::utils::strip_prefix;

#[cfg(feature = "experimental")]
use crate::subcommands::scan::scan_segregate::CoordSamples;

#[allow(clippy::upper_case_acronyms)]
pub enum FileType {
    BED,
    VCF,
    CSV,
    TSV,
    HST,
}

impl FileType {
    pub fn from_path(path: &PathBuf) -> Result<Self> {
        let extension: &str = Path::new(&path)
            .extension()
            .and_then(OsStr::to_str)
            .ok_or_else(|| eyre!("No filetype in path"))?;

        let extension = match extension {
            "gz" | "bgz" => return_double_extension_filetype(path, extension)?,
            _ => extension.to_string(),
        };

        Ok(match extension.as_str() {
            "vcf.gz" | "vcf" | "bcf" => Self::VCF,
            "hst.gz" | "hst" => Self::HST,
            "bed.gz" | "bed" => Self::BED,
            "csv" | "ids" => Self::CSV,
            "tsv" => Self::TSV,
            _ => return Err(eyre!("File extension: {extension} is not supported")),
        })
    }
}

pub fn return_double_extension_filetype(path: &Path, e1: &str) -> Result<String> {
    let stem = path
        .file_stem()
        .and_then(OsStr::to_str)
        .ok_or_else(|| eyre!("file has no stem"))?;
    let e2 = Path::new(&stem)
        .extension()
        .and_then(OsStr::to_str)
        .ok_or_else(|| eyre!("file has no other filetype"))?;
    Ok(format!("{e2}.{e1}"))
}

pub fn read_coords_file(path: &PathBuf) -> Result<Vec<Coord>> {
    let input = get_input(Some(path.clone()))?;

    let mut rdr = match FileType::from_path(path)? {
        FileType::CSV => get_csv_reader(input),
        FileType::TSV => get_tsv_reader(input, false),
        FileType::BED => get_tsv_reader(input, false),
        _ => return Err(eyre!("Filetype from {path:?} is not supported")),
    };

    let mut variants: Vec<Coord> = vec![];

    for line in rdr.records() {
        let record = line?;
        let variant: Coord = record.deserialize(None)?;
        variants.push(variant);
    }

    Ok(variants)
}

#[cfg(feature = "experimental")]
pub fn read_coords_sample_file(path: &PathBuf) -> Result<Vec<CoordSamples>> {
    let input = get_input(Some(path.clone()))?;

    let mut rdr = match FileType::from_path(path)? {
        FileType::CSV => get_csv_reader(input),
        FileType::TSV => get_tsv_reader(input, false),
        FileType::BED => get_tsv_reader(input, false),
        _ => return Err(eyre!("Filetype from {path:?} is not supported")),
    };

    let mut variants: Vec<CoordSamples> = vec![];

    for line in rdr.records() {
        let record = line?;
        let variant: CoordSamples = record.deserialize(None)?;
        variants.push(variant);
    }

    Ok(variants)
}

pub fn filter_out_non_vcf_sample_names(
    path: &PathBuf,
    contig: &str,
    wanted_samples: Vec<String>,
) -> Result<Vec<String>> {
    let reader = get_reader(path, contig, None)?;

    let samples = reader
        .header()
        .samples()
        .into_iter()
        .map(|sample| Ok(String::from_utf8_lossy(sample).to_string()))
        .collect::<Result<Vec<String>>>()?;

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

#[derive(Debug, Clone, serde::Deserialize)]
struct RecombinationRow<'a> {
    chr: &'a str,
    pos: u64,
    _rate: f32,
    cm: f32,
}

pub fn check_recombination_file_matches_contig(path: PathBuf, contig: &str) -> Result<bool> {
    let mut rdr = get_tsv_reader(get_input(Some(path))?, false);
    let mut lines = rdr.records();

    let line = lines
        .next()
        .ok_or_eyre("Recombination rates file {path} contains 0 rows.")?;

    let record = line?;
    let row: RecombinationRow = record.deserialize(None).
            wrap_err(eyre!("Make sure no headers are present and that the recombination file is in order chr,pos,rate,cm. If the issue is not fixed, you have an invalid field in the file."))?;
    Ok(row.chr == contig)
}

pub fn read_recombination_file(path: PathBuf) -> Result<BTreeMap<u64, f32>> {
    let mut rdr = get_tsv_reader(get_input(Some(path))?, false);
    let mut rates = BTreeMap::new();

    let mut last_cm = 0.0;
    let mut last_pos = 0;
    for line in rdr.records() {
        let record = line?;
        let row: RecombinationRow = record.deserialize(None).
            wrap_err(eyre!("Make sure no headers are present and that the recombination file is in order chr,pos,rate,cm. If the issue is not fixed, you have an invalid field in the file."))?;
        if last_cm > row.cm {
            return Err(eyre!(
                "Recombination rates file is not sorted for centimorgans"
            ));
        }
        if last_pos > row.pos {
            return Err(eyre!(
                "Recombination rates file is not sorted for positions"
            ));
        }
        last_cm = row.cm;
        last_pos = row.pos;
        rates.insert(row.pos, row.cm);
    }

    Ok(rates)
}

pub fn read_haplotype_file(ht_path: PathBuf) -> Result<Vec<HapVariant>> {
    let input = get_input(Some(ht_path))?;
    let mut rdr = get_csv_reader(input);

    let mut variants: Vec<HapVariant> = vec![];

    for line in rdr.records() {
        let record = line?;
        let variant: HapVariant = record.deserialize(None)?;
        if let Some(latest) = variants.last() {
            if latest.pos > variant.pos {
                return Err(eyre!(
                    "The haplotype file is not sorted by position. {} is larger than {}",
                    latest,
                    variant
                ));
            }
        }
        variants.push(variant);
    }

    Ok(variants)
}

pub fn read_lines<P>(filename: P) -> Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    // let file = File::open(filename)?;
    let name = filename.as_ref().display();
    let file = match File::open(&filename) {
        Ok(x) => x,
        Err(err) => {
            let msg = format!("failed to open {name}: {err}");
            return Err(std::io::Error::new(std::io::ErrorKind::NotFound, msg))?;
        }
    };
    Ok(io::BufReader::new(file).lines())
}

pub fn read_multiple_sample_ids(path: &Option<Vec<PathBuf>>) -> Result<Option<Vec<String>>> {
    match path {
        Some(paths) => {
            let mut samples = vec![];
            for path in paths {
                for line in read_lines(path)?.map_while(Result::ok) {
                    let line = line.trim();
                    samples.push(line.to_string());
                }
            }
            Ok(Some(samples))
        }
        None => Ok(None),
    }
}

pub fn read_sample_ids(path: &Option<PathBuf>) -> Result<Option<Vec<String>>> {
    match path {
        Some(path) => {
            let mut samples = vec![];

            for line in read_lines(path)?.map_while(Result::ok) {
                let line = line.trim();
                samples.push(line.to_string());
            }
            Ok(Some(samples))
        }
        None => Ok(None),
    }
}

pub fn push_to_output(args: &StandardArgs, output: &mut PathBuf, name: &str, suffix: &str) {
    if let Some(prefix) = &strip_prefix(args.prefix.clone()) {
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
                    .ok_or_eyre("The input VCF has no ID for some contig record in the header")?;
                Ok((
                    id.clone(),
                    values
                        .get("length")
                        .ok_or_eyre("VCF header has no contig length for {id}")?
                        .parse::<i64>()?,
                ))
            }
            _ => {
                panic!("Header contig filtering failed for some reason. Create an issue at Github")
            }
        })
        .collect()
}

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

pub fn write_haplotype(
    haplotype: Vec<HapVariant>,
    mut writer: csv::Writer<Box<dyn std::io::Write>>,
) -> Result<()> {
    writer.write_record(vec!["contig", "pos", "ref", "alt", "gt"])?;

    for coord in haplotype {
        writer.write_record(vec![
            coord.contig,
            coord.pos.to_string(),
            coord.reference,
            coord.alt,
            coord.gt.to_string(),
        ])?;
    }
    Ok(())
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

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
