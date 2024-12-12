use std::ffi::OsStr;
use std::io::Read;
use std::path::{Path, PathBuf};

use bio::data_structures::interval_tree::IntervalTree;
use color_eyre::{eyre::eyre, Result};
use csv::Reader;
use polars::frame::DataFrame;
use polars::io::SerReader;
use polars::prelude::{CsvReader, NullValues};

use crate::io::{
    get_csv_writer, get_input, get_output, get_tsv_reader, return_double_extension_filetype,
    FileType,
};
use crate::subcommands::compare_haplotypes::get_csv_reader;

pub fn read_variable_data_file(path: PathBuf) -> Result<DataFrame> {
    let df = CsvReader::from_path(path)?
        .with_null_values(Some(NullValues::AllColumnsSingle("NA".to_string())))
        .infer_schema(Some(500))
        .has_header(true)
        .finish()?;

    Ok(df)
}

pub fn read_delimited_file(file: &PathBuf) -> Result<Reader<Box<dyn Read>>> {
    let input = get_input(Some(file.clone()))?;

    let extension: &str = Path::new(&file)
        .extension()
        .and_then(OsStr::to_str)
        .ok_or_else(|| eyre!("No filetype in path"))?;

    let extension = match extension {
        "gz" | "bgz" => return_double_extension_filetype(file, extension)?,
        _ => extension.to_string(),
    };

    let file_type = match extension.as_str() {
        "vcf.gz" | "vcf" | "bcf" => FileType::VCF,
        "hst.gz" | "hst" => FileType::HST,
        "bed.gz" | "bed" => FileType::BED,
        "csv" | "ids" => FileType::CSV,
        "tsv" => FileType::TSV,
        _ => FileType::CSV,
    };

    let rdr = match file_type {
        FileType::CSV => get_csv_reader(input),
        FileType::TSV => get_tsv_reader(input, false),
        FileType::BED => get_tsv_reader(input, false),
        _ => return Err(eyre!("Filetype from {file:?} is not supported")),
    };

    Ok(rdr)
}
#[doc(hidden)]
pub fn run(file: PathBuf, annotate_file: PathBuf) -> Result<()> {
    let mut annotation: IntervalTree<i64, String> = IntervalTree::new();

    let start = 1;
    let stop = 2;
    let data = 3;

    for line in read_delimited_file(&annotate_file)?.records() {
        let record = line?;
        let mut row: Vec<String> = record.deserialize(None)?;

        let start = row[start].parse::<i64>().unwrap();
        let stop = row[stop].parse::<i64>().unwrap();

        annotation.insert(start..stop, row.swap_remove(data));
    }

    let file_start = 7;
    let file_stop = 8;

    let mut wrtr = get_csv_writer(get_output(None).unwrap());

    for line in read_delimited_file(&file)?.records() {
        let record = line?;
        let mut row: Vec<String> = record.deserialize(None)?;

        let start = row[file_start].parse::<i64>().unwrap();
        let stop = row[file_stop].parse::<i64>().unwrap();

        let anns = annotation
            .find(start..stop)
            .map(|data| data.data())
            .cloned()
            .collect::<Vec<String>>();

        match anns.is_empty() {
            true => {
                row.push("".to_string());
                wrtr.write_record(row)?;
            }
            false => {
                row.push(anns.join(",").to_string());
                wrtr.write_record(row)?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn test_extension_filetype() {
        let path = std::path::PathBuf::from("test.vcf.gz");
        let ftype = return_double_extension_filetype(&path, "gz").unwrap();
        assert_eq!(String::from("vcf.gz"), ftype);
    }
}
