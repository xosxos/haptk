use std::path::PathBuf;
use std::{collections::HashMap, fs::File};

use bgzip::BGZFReader;
use bio::data_structures::interval_tree::IntervalTree;
use color_eyre::{
    eyre::{eyre, WrapErr},
    Result,
};
use csv::{QuoteStyle, ReaderBuilder, WriterBuilder};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::core::get_output;
use crate::io::read_haplotype_file;

#[derive(Deserialize, Serialize, Default, Debug, Clone)]
pub struct GtfRow {
    pub seqname: String,
    pub source: String,
    pub feature: String,
    pub start: u64,
    pub stop: u64,
    pub score: String,
    pub strand: String,
    pub frame: String,
    pub attribute: String,
}

pub fn run(path: PathBuf, gtf_path: PathBuf, output: PathBuf) -> Result<()> {
    let ht = read_haplotype_file(path)?;
    let gtf_hm = read_gtf_file(gtf_path)?;
    println!("{gtf_hm:?}");

    let output = get_output(Some(output))?;

    let mut wrtr = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .double_quote(false)
        .quote_style(QuoteStyle::Never)
        .from_writer(output);

    ht.into_iter().try_for_each(|row| -> Result<()> {
        if let Some(gtf_it) = gtf_hm.get(&row.contig) {
            let gtf_row = gtf_it.find(row.pos..row.pos + 1);
            let record = match gtf_row.last() {
                Some(gtf_row) => vec![
                    row.contig,
                    row.pos.to_string(),
                    row.reference,
                    row.alt,
                    row.gt.to_string(),
                    gtf_row.data().attribute.clone(),
                ],

                None => {
                    vec![
                        row.contig,
                        row.pos.to_string(),
                        row.reference,
                        row.alt,
                        row.gt.to_string(),
                        "".to_string(),
                    ]
                }
            };
            wrtr.write_record(record)?;
        } else {
            tracing::warn!("Annotation file has no seqid {:?}", row.contig);
        }
        Ok(())
    })
}

pub fn read_gtf_file(path: PathBuf) -> Result<HashMap<String, IntervalTree<u64, GtfRow>>> {
    let file = File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;

    let reader = BGZFReader::new(file)?;

    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(reader);

    let mut hm: HashMap<String, Vec<GtfRow>> = HashMap::new();
    for line in rdr.records() {
        let record = line?;
        let row: GtfRow = record
            .deserialize(None)
            .wrap_err(eyre!("Error record {record:?} from {path:?}."))?;

        hm.entry(row.seqname.clone())
            .or_insert(Vec::new())
            .push(row);
    }

    let keys: Vec<String> = hm.keys().cloned().collect();
    let vec = keys
        .par_iter()
        .map(|key| {
            let mut it = IntervalTree::<u64, GtfRow>::new();
            for row in hm.get(key).unwrap() {
                it.insert(row.start..row.stop + 1, row.clone());
            }
            (key.clone(), it)
        })
        .collect::<Vec<(String, IntervalTree<u64, GtfRow>)>>();
    let mut it_hm: HashMap<String, IntervalTree<u64, GtfRow>> = HashMap::new();

    for (key, it) in vec.into_iter() {
        it_hm.insert(key, it);
    }

    Ok(it_hm)
}

pub fn read_gtf_file_to_vec(path: PathBuf, contig: &String) -> Result<Vec<GtfRow>> {
    let file = File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;

    let reader = BGZFReader::new(file)?;

    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(reader);

    let mut hm: HashMap<String, Vec<GtfRow>> = HashMap::new();
    for line in rdr.records() {
        let record = line?;
        let row: GtfRow = record
            .deserialize(None)
            .wrap_err(eyre!("Error record {record:?} from {path:?}."))?;

        hm.entry(row.seqname.clone())
            .or_insert(Vec::new())
            .push(row);
    }

    Ok(hm.get(contig).unwrap().clone())
}
