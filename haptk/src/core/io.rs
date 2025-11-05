use std::collections::BTreeMap;

use std::collections::HashMap;
use std::ffi::{OsStr, OsString};
use std::io;
use std::path::PathBuf;
use std::sync::mpsc::Receiver;
use std::thread::{self, JoinHandle};

use bgzip::tabix::Tabix;
use color_eyre::eyre::{eyre, OptionExt, WrapErr};
use color_eyre::Result;
use csv::{Reader, ReaderBuilder, Writer, WriterBuilder};
use serde::Deserialize;
use serde::Serialize;

use crate::args::{Selection, StandardArgs};
use crate::core::Coord;
use crate::core::HapVariant;
use crate::error::Error;
use crate::utils::strip_prefix;

#[cfg(feature = "experimental")]
use crate::subcommands::scan::scan_segregate::CoordSamples;

pub use haptk_core::io::get_csv_reader;
pub use haptk_core::io::get_csv_writer;
pub use haptk_core::io::get_extension;
pub use haptk_core::io::get_input;
pub use haptk_core::io::get_output;
pub use haptk_core::io::get_vcf_writer;
pub use haptk_core::io::read_lines;
pub use haptk_core::io::read_multiple_sample_ids;
pub use haptk_core::io::FileType;
pub use haptk_core::vcf::contig_len_from_vcf;
pub use haptk_core::vcf::filter_out_samples_not_in_vcf;
pub use haptk_core::vcf::get_indexes_and_sample_ids_from_vcf;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct SampleHaplotypeList {
    pub sample: String,
    pub haplotypes: String,
}

pub fn read_sample_ht_list_file(path: &PathBuf) -> Result<HashMap<String, [bool; 2]>> {
    let input = get_input(Some(path.clone()))?;

    let mut rdr = match FileType::from_path(path)? {
        FileType::CSV => get_csv_reader(input, false),
        FileType::TSV => get_tsv_reader(input, false),
        FileType::BED => get_tsv_reader(input, false),
        _ => {
            return Err(eyre!(Error::FileNotSupported {
                ext: format!("{path:?}")
            }))
        }
    };

    let mut variants = HashMap::new();

    for line in rdr.records() {
        let record = line?;
        let variant: SampleHaplotypeList = record.deserialize(None)?;

        let haplotypes: Vec<_> = variant.haplotypes.split(":").collect();

        if haplotypes.len() > 2 {
            return Err(eyre!(
                "Detected more than 2 chromosomes  for {}. Only max. diploid organisms are supported.", variant.sample
            ));
        }

        if haplotypes.is_empty() {
            return Err(eyre!("No haplotypes for sample {}.", variant.sample));
        }

        if haplotypes.len() == 2 {
            variants.insert(variant.sample.clone(), [true, true]);
        }

        if haplotypes.len() == 1 {
            if haplotypes[0] == "1" {
                variants.insert(variant.sample, [false, true]);
            } else {
                variants.insert(variant.sample, [true, false]);
            }
        }
    }

    Ok(variants)
}

pub fn read_coords_file(path: &PathBuf) -> Result<Vec<Coord>> {
    let input = get_input(Some(path.clone()))?;

    let mut rdr = match FileType::from_path(path)? {
        FileType::CSV => get_csv_reader(input, false),
        FileType::TSV => get_tsv_reader(input, false),
        FileType::BED => get_tsv_reader(input, false),
        _ => {
            return Err(eyre!(Error::FileNotSupported {
                ext: format!("{path:?}")
            }))
        }
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
        FileType::CSV => get_csv_reader(input, false),
        FileType::TSV => get_tsv_reader(input, false),
        FileType::BED => get_tsv_reader(input, false),
        _ => {
            return Err(eyre!(Error::FileNotSupported {
                ext: format!("{path:?}")
            }))
        }
    };

    let mut variants: Vec<CoordSamples> = vec![];

    for line in rdr.records() {
        let record = line?;
        let variant: CoordSamples = record.deserialize(None)?;
        variants.push(variant);
    }

    Ok(variants)
}

#[derive(Debug, Clone, serde::Deserialize)]
struct RecombinationRow<'a> {
    chr: &'a str,
    pos: u64,
    _rate: f32,
    cm: f32,
}

pub fn check_recombination_file_matches_contig(path: PathBuf, contig: &str) -> Result<bool> {
    let mut rdr = get_tsv_reader(get_input(Some(path.clone()))?, false);
    let mut lines = rdr.records();

    let record = lines
        .next()
        .ok_or_eyre(Error::EmptyFile { path: path.clone() })??;

    let row: RecombinationRow = record
        .deserialize(None)
        .wrap_err(Error::RecombinationDeserialization { path })?;

    Ok(row.chr == contig)
}

pub fn read_recombination_file(path: PathBuf) -> Result<BTreeMap<u64, f32>> {
    let mut rdr = get_tsv_reader(get_input(Some(path.clone()))?, false);
    let mut rates = BTreeMap::new();

    let mut last_cm = 0.0;
    let mut last_pos = 0;

    for line in rdr.records() {
        let record = line?;

        let row: RecombinationRow = record
            .deserialize(None)
            .wrap_err(Error::RecombinationDeserialization { path: path.clone() })?;

        // Check cm sorting
        if last_cm > row.cm {
            return Err(eyre!(Error::Sort {
                prev_pos: last_cm,
                pos: row.cm,
                path,
            }));
        }

        // Check position sorting
        if last_pos > row.pos {
            return Err(eyre!(Error::Sort {
                prev_pos: last_pos as f32,
                pos: row.pos as f32,
                path,
            }));
        }

        last_cm = row.cm;
        last_pos = row.pos;
        rates.insert(row.pos, row.cm);
    }

    Ok(rates)
}

pub fn read_haplotype_file(ht_path: PathBuf) -> Result<Vec<HapVariant>> {
    let input = get_input(Some(ht_path.clone()))?;
    let mut rdr = get_csv_reader(input, true);

    let mut variants: Vec<HapVariant> = vec![];

    for line in rdr.records() {
        let record = line?;
        let variant: HapVariant = record.deserialize(None)?;

        if let Some(latest) = variants.last() {
            // Check sorting
            if latest.pos > variant.pos {
                return Err(eyre!(Error::Sort {
                    prev_pos: latest.pos as f32,
                    pos: variant.pos as f32,
                    path: ht_path,
                }));
            }
        }
        variants.push(variant);
    }

    Ok(variants)
}

pub fn push_to_output(args: &StandardArgs, output: &mut PathBuf, name: &str, suffix: &str) {
    if let Some(prefix) = &strip_prefix(args.prefix.clone()) {
        match args.selection {
            Selection::All => output.push(format!("{prefix}_{name}.{suffix}")),
            Selection::OnlyAlts => output.push(format!("{prefix}_{name}_only_alts.{suffix}")),
            Selection::OnlyRefs => output.push(format!("{prefix}_{name}_only_refs.{suffix}")),
            Selection::OnlyLongest => output.push(format!("{prefix}_{name}_only_longest.{suffix}")),
            Selection::List => output.push(format!("{prefix}_{name}_list.{suffix}")),
            Selection::Unphased => output.push(format!("{prefix}_{name}_rwc.{suffix}")),
            Selection::Haploid => output.push(format!("{prefix}_{name}_haploid.{suffix}")),
        }
    } else {
        match args.selection {
            Selection::All => output.push(format!("{name}.{suffix}")),
            Selection::OnlyAlts => output.push(format!("{name}_only_alts.{suffix}")),
            Selection::OnlyRefs => output.push(format!("{name}_only_refs.{suffix}")),
            Selection::OnlyLongest => output.push(format!("{name}_only_longest.{suffix}")),
            Selection::List => output.push(format!("{name}_list.{suffix}")),
            Selection::Unphased => output.push(format!("{name}_rwc.{suffix}")),
            Selection::Haploid => output.push(format!("{name}_haploid.{suffix}")),
        }
    }
}

pub fn get_tsv_reader<R: io::Read>(input: R, has_headers: bool) -> Reader<R> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(has_headers)
        .flexible(false)
        .from_reader(input)
}

pub fn get_strict_tsv_writer<W: io::Write>(output: W) -> Writer<W> {
    WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_writer(output)
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

pub fn read_tabix(path: &PathBuf) -> Result<HashMap<String, (u64, usize)>> {
    fn append_ext(ext: impl AsRef<OsStr>, path: &PathBuf) -> PathBuf {
        let mut os_string: OsString = path.into();
        os_string.push(".");
        os_string.push(ext.as_ref());
        os_string.into()
    }

    let path = append_ext("tbi", path);
    let mut file = std::fs::File::open(&path).wrap_err(Error::Io { path: path.clone() })?;

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

    #[test]
    fn test_extension_filetype() {
        let path = std::path::PathBuf::from("test.vcf.gz");
        let ftype = get_extension(&path).unwrap();
        assert_eq!(String::from("vcf.gz"), ftype);
    }
}
