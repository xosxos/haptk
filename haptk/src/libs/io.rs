use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};

use color_eyre::{
    eyre::{eyre, WrapErr},
    Result,
};
use polars::prelude::*;

use crate::core::{get_csv_reader, get_input, get_tsv_reader};
use crate::structs::HapVariant;

pub fn read_variable_data_file(path: PathBuf) -> Result<DataFrame> {
    let df = CsvReader::from_path(path)?
        .with_null_values(Some(NullValues::AllColumnsSingle("NA".to_string())))
        .infer_schema(Some(500))
        .has_header(true)
        .finish()?;

    Ok(df)
}

#[derive(Debug, Clone, serde::Deserialize)]
struct RecombinationRow<'a> {
    _chr: &'a str,
    pos: u64,
    _rate: f32,
    cm: f32,
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

pub fn read_sample_ids(path: &Option<PathBuf>) -> Result<Option<Vec<String>>> {
    match path {
        Some(path) => {
            let mut samples = vec![];

            for line in read_lines(path)?.flatten() {
                let line = line.trim();
                samples.push(line.to_string());
            }
            Ok(Some(samples))
        }
        None => Ok(None),
    }
}
