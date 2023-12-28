use std::path::PathBuf;

use color_eyre::Result;
use serde::{Deserialize, Serialize};

use crate::{
    io::{get_output, get_vcf_writer, read_haplotype_file},
    structs::{HapVariant, Ploidy},
};

pub fn run(path: PathBuf, sample_name: String, output: PathBuf) -> Result<()> {
    let ht = read_haplotype_file(path)?;
    let header = header(vec![sample_name], &ht.first().unwrap().contig);

    let output = get_output(Some(output))?;
    let mut wrtr = get_vcf_writer(output);

    for line in header {
        wrtr.write_record(vec![line])?;
    }

    ht.into_iter()
        .map(Into::<VcfRow>::into)
        .try_for_each(|row| -> Result<()> {
            wrtr.write_record(Into::<Vec<String>>::into(row))?;
            Ok(())
        })
}

#[derive(Deserialize, Serialize, Default, Debug, Clone)]
pub struct VcfRow {
    pub seqid: String,
    pub pos: u64,
    pub varid: String,
    pub reference: String,
    pub alt: String,
    pub quality: String,
    pub filter: String,
    pub info: String,
    pub format: String,
    pub samples: Vec<String>,
}

impl VcfRow {
    pub fn new(
        seqid: String,
        pos: u64,
        reference: String,
        alt: String,
        samples: Vec<(u8, u8)>,
        ploidy: Ploidy,
    ) -> Self {
        let samples: Vec<String> = samples
            .into_iter()
            .map(|(a, b)| match ploidy {
                Ploidy::Haploid => format!("{a}"),
                Ploidy::Diploid => format!("{a}|{b}"),
                Ploidy::Mixed => format!("{a}|{b}"),
            })
            .collect();

        Self {
            varid: format!("{seqid}_{pos}_{reference}_{alt}"),
            pos,
            seqid,
            reference,
            alt,
            quality: String::from("."),
            filter: String::from("PASS"),
            info: String::from("."),
            format: String::from("GT"),
            samples,
        }
    }
}

impl From<HapVariant> for VcfRow {
    fn from(value: HapVariant) -> Self {
        Self {
            pos: value.pos,
            varid: format!(
                "{}_{}_{}_{}",
                &value.contig, value.pos, value.reference, value.alt
            ),
            seqid: value.contig,
            reference: value.reference,
            alt: value.alt,
            quality: String::from("."),
            filter: String::from("PASS"),
            info: String::from("."),
            format: String::from("GT"),
            samples: vec![format!("{}", value.gt)],
        }
    }
}

impl From<VcfRow> for Vec<String> {
    fn from(row: VcfRow) -> Vec<String> {
        let mut record = vec![
            row.seqid,
            row.pos.to_string(),
            row.varid,
            row.reference,
            row.alt,
            row.quality,
            row.filter,
            row.info,
            row.format,
        ];
        record.extend(row.samples);
        record
    }
}

pub fn header<S: AsRef<str>, T: AsRef<str>>(sample_names: Vec<S>, contig: T) -> Vec<String> {
    let time = chrono::offset::Utc::now();
    let mut header = vec![];

    let samples_string =
        sample_names
            .iter()
            .enumerate()
            .fold(String::new(), |mut acc, (i, name)| {
                match i == sample_names.len() - 1 {
                    true => acc.push_str(name.as_ref()),
                    false => acc.push_str(&format!("{}\t", name.as_ref())),
                };
                acc
            });

    header.push("##fileformat=VCFv4.2".into());
    header.push(format!("##fileDate={time}"));
    header.push("##source=HAPTK".to_string());
    header.push("##FILTER=<ID=PASS,Description=\"All filters passed\">".into());
    header.push("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".into());
    header.push(format!("##contig=<ID={}>", contig.as_ref()));
    header.push(format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples_string}",
    ));

    header
}
