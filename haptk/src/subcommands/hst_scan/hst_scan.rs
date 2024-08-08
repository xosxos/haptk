use std::{borrow::Borrow, collections::HashMap, hash::DefaultHasher, path::PathBuf, sync::Arc};

use std::hash::Hasher;
use std::sync::mpsc::Sender;

use petgraph::graph::NodeIndex;
use petgraph::Graph;
use rayon::prelude::*;
use rust_htslib::bam::{
    header::HeaderRecord,
    record::{Aux, AuxArray, Cigar, CigarString},
    Format, Header, Record, Writer,
};

use crate::args::{Selection, StandardArgs};
use crate::io::push_to_output;
use crate::structs::{Coord, HapVariant, PhasedMatrix};
use crate::utils::parse_coords;

use color_eyre::{
    eyre::{ensure, eyre, WrapErr},
    Result,
};
use serde::{Deserialize, Serialize};

use crate::io::get_output;
use crate::read_vcf::read_vcf_to_matrix;
use crate::subcommands::bhst::{construct_bhst, Node};

#[doc(hidden)]
pub fn run(args: StandardArgs, step_size: usize) -> Result<()> {
    ensure!(
        args.selection == Selection::All || args.selection == Selection::Haploid,
        "Running only with phased data and all chromsomes is supported."
    );

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "trees", "json.gz");

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = match (start, stop) {
        (Some(start), Some(stop)) => {
            read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None)?
        }
        _ => read_vcf_to_matrix(&args, contig, 0, None, None)?,
    };

    // Create a Vec because par_iter cannot be used with pure ranges and par_bridge does not
    // return in ordered fashion with .collect()
    let range: Vec<_> = (0..vcf.ncoords()).collect();

    let hsts = range
        .par_iter()
        .filter(|idx| *idx % step_size == 0)
        .map(|idx| {
            if idx % (vcf.ncoords() / 10) == 0 {
                tracing::info!("Finished constructing 10% of the trees...");
            }

            HstScanRow {
                idx: *idx,
                hst: construct_bhst(&vcf, *idx, 4),
            }
        })
        .collect();

    let metadata = Metadata::new(&vcf, args);

    let trees = HstScan { hsts, metadata };

    write_hsts(trees, output)?;

    Ok(())
}

pub type Limits = (usize, usize, usize, usize);

#[derive(Serialize, Deserialize, Clone)]
pub struct HstScanRow {
    pub idx: usize,
    pub hst: Graph<Node, u8>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Metadata {
    pub coords: Vec<Coord>,
    pub contig: String,
    pub samples: Vec<String>,
    pub selection: Selection,
    pub vcf_name: PathBuf,
    pub info_limit: Option<f32>,
}

impl Metadata {
    fn new(vcf: &PhasedMatrix, args: StandardArgs) -> Self {
        let samples = vcf.samples().to_vec();

        Self {
            coords: vcf.coords().clone(),
            contig: vcf.get_contig().clone(),
            samples,
            selection: args.selection,
            vcf_name: args.file,
            info_limit: args.info_limit,
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct HstScan {
    pub metadata: Metadata,
    pub hsts: Vec<HstScanRow>,
}

impl HstScan {
    pub fn get_pos(&self, idx: usize) -> u64 {
        self.metadata.coords[idx].pos
    }

    pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
        let mut names = vec![];
        for i in indexes {
            names.push(self.metadata.samples[*i].clone());
        }
        names
    }
}

fn write_hsts(trees: HstScan, path: PathBuf) -> Result<()> {
    // use std::io::{BufWriter, Write};

    tracing::info!("Writing trees to file.");
    let now = std::time::Instant::now();
    let mut output = get_output(Some(path))?;

    let mut writer = bgzip::BGZFWriter::new(&mut output, bgzip::Compression::default());
    // let mut writer = BufWriter::new(output);

    serde_json::to_writer(&mut writer, &trees)?;

    writer.close()?;
    // writer.flush()?;

    tracing::info!("Finished writing trees to file.");
    tracing::debug!("Wrote trees file in {:?}", now.elapsed());
    Ok(())
}

pub fn read_tree_file(path: PathBuf) -> Result<HstScan> {
    let file = std::fs::File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;
    let reader = bgzip::BGZFReader::new(file)?;
    let hst_scan: HstScan = serde_json::from_reader(reader)?;
    tracing::info!(
        "Read HSTs with the following metadata: \nselection: {:?},\nnsamples:{},\ninfo_limit: {:?}",
        hst_scan.metadata.selection,
        hst_scan.metadata.samples.len(),
        hst_scan.metadata.info_limit
    );

    Ok(hst_scan)
}

// impl std::fmt::Display for Node {
//     fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
//         let nmarkers = match self.stop_idx.cmp(&self.start_idx) {
//             std::cmp::Ordering::Greater => self.stop_idx - self.start_idx + 1,
//             std::cmp::Ordering::Less => self.start_idx - self.stop_idx + 1,
//             std::cmp::Ordering::Equal => 0,
//         };
//         let line = format!(
//             "n: {}\nstart_idx: {}\nstop_idx: {}\nnmarkers: {}\nblock length (bp): {}\n",
//             self.indexes.len(),
//             self.start_idx,
//             self.stop_idx,
//             nmarkers,
//             self.haplotype.len(),
//         );
//         write!(f, "{line}")
//     }
// }

pub fn get_sender(
    write_bam: bool,
    args: &StandardArgs,
    hsts: Arc<HstScan>,
) -> Option<Sender<(usize, NodeIndex, f64)>> {
    match write_bam {
        true => {
            let mut bam_output = args.output.clone();
            push_to_output(args, &mut bam_output, "hst_scan", "bam");
            let header = return_header(&hsts.metadata.contig);
            let mut out = Writer::from_path(bam_output, &header, Format::Bam).unwrap();
            let (tx, rx) = std::sync::mpsc::channel();

            std::thread::spawn(move || {
                while let Ok((tree_idx, top_node_idx, optimized_value)) = rx.recv() {
                    let record =
                        create_bam_record(hsts.clone(), tree_idx, top_node_idx, optimized_value);
                    out.write(&record).unwrap();
                }
            });
            Some(tx)
        }
        false => None,
    }
}

pub fn return_assoc<F, U>(
    hsts: Arc<HstScan>,
    args: &StandardArgs,
    limits: Limits,
    write_bam: bool,
    optimizer: F,
    rower: U,
) -> Vec<AssocRow>
where
    F: Fn(Arc<HstScan>, usize, Limits) -> Option<(usize, NodeIndex, f64)> + Sync + Send + Clone,
    U: Fn(Arc<HstScan>, usize, NodeIndex, f64) -> AssocRow + Sync + Send,
{
    let tx = get_sender(write_bam, &args, hsts.clone());

    hsts.hsts
        .par_iter()
        .filter_map(|hst_row| optimizer(hsts.clone(), hst_row.idx, limits))
        .map_with(tx, |tx, optimizer_tuple| {
            if let Some(tx) = tx {
                tx.send(optimizer_tuple).unwrap();
            }
            optimizer_tuple
        })
        .map(|(tree_idx, top_node_idx, optimized_value)| {
            rower(hsts.clone(), tree_idx, top_node_idx, optimized_value)
        })
        .collect()
}

pub fn write_assoc(
    header: &[&str],
    rows: Vec<AssocRow>,
    mut writer: csv::Writer<Box<dyn std::io::Write>>,
) -> Result<()> {
    writer.write_record(header)?;

    rows.into_iter().try_for_each(|row| -> Result<()> {
        let record = row.to_csv_row();
        writer.write_record(record)?;
        Ok(())
    })
}

#[derive(Clone, Debug)]
pub struct AssocRow {
    pub contig: String,
    pub pos: u64,
    pub marker_id: String,
    pub bp_len: f32,
    pub marker_len: usize,
    pub start: u64,
    pub stop: u64,
    pub start_idx: usize,
    pub stop_idx: usize,
    pub opt_value: f64,
    pub node_idx: NodeIndex,
    pub first_sample_idx: usize,
    pub hashmap: HashMap<String, String>,
}

impl AssocRow {
    pub fn new(
        hsts: Arc<HstScan>,
        tree_idx: usize,
        top_node_idx: NodeIndex,
        opt_value: f64,
        hashmap: HashMap<String, String>,
    ) -> AssocRow {
        let hst = &hsts.hsts[tree_idx].hst;
        let top_node = hst.node_weight(top_node_idx).unwrap();

        AssocRow {
            contig: hsts.metadata.contig.clone(),
            pos: hsts.get_pos(tree_idx),
            marker_id: get_marker_id(top_node, hsts.clone()),
            bp_len: (hsts.get_pos(top_node.stop_idx) - hsts.get_pos(top_node.start_idx)) as f32,
            marker_len: top_node.stop_idx - top_node.start_idx + 1,
            start: hsts.get_pos(top_node.start_idx),
            stop: hsts.get_pos(top_node.stop_idx),
            start_idx: top_node.start_idx,
            stop_idx: top_node.stop_idx,
            opt_value,
            hashmap,
            node_idx: top_node_idx,
            first_sample_idx: top_node.indexes[0],
        }
    }

    fn pos(&self) -> u64 {
        self.pos
    }
    fn bp_len(&self) -> f32 {
        self.bp_len
    }
    fn marker_len(&self) -> usize {
        self.marker_len
    }
    fn node_idx(&self) -> NodeIndex {
        self.node_idx
    }
    fn start_idx(&self) -> usize {
        self.start_idx
    }
    fn stop_idx(&self) -> usize {
        self.stop_idx
    }
    fn first_sample_idx(&self) -> usize {
        self.first_sample_idx
    }
    fn opt_value(&self) -> f64 {
        self.opt_value
    }
    fn hashmap(&self) -> &std::collections::HashMap<String, String> {
        &self.hashmap
    }
    fn show(&self) -> String {
        self.to_string()
    }
    fn to_csv_row(&self) -> Vec<String> {
        self.into()
    }
}

impl From<&AssocRow> for Vec<String> {
    fn from(row: &AssocRow) -> Vec<String> {
        vec![
            row.contig.to_string(),
            row.pos.to_string(),
            row.marker_id.to_string(),
            row.bp_len.to_string(),
            row.marker_len.to_string(),
            row.start.to_string(),
            row.stop.to_string(),
            row.opt_value.to_string(),
        ]
    }
}

impl std::fmt::Display for AssocRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let vector: Vec<String> = self.into();
        let line = format!("{:?}", vector);
        write!(f, "{line}")
    }
}

pub fn get_marker_id<T: Borrow<HstScan>>(top_node: &Node, hsts: T) -> String {
    let hsts = hsts.borrow();
    let start = hsts.get_pos(top_node.start_idx);
    let stop = hsts.get_pos(top_node.stop_idx);

    let ht = top_node.haplotype.iter().fold(
        format!("{}_{}_{}_", hsts.metadata.contig, start, stop),
        |acc, e| format!("{acc}{e}"),
    );

    let mut hasher = DefaultHasher::new();
    hasher.write(ht.as_bytes());

    format!("{:02x}", hasher.finish())
}

/// BAM creator
pub fn return_header(contig: &str) -> Header {
    let mut header = Header::new();
    let mut header_record = HeaderRecord::new(b"SQ");
    header_record.push_tag(b"SN", contig);
    header_record.push_tag(b"LN", 0);
    header.push_record(&header_record);
    header
}

pub fn u8_to_hapvariant(hsts: Arc<HstScan>, ht: &[u8]) -> Vec<HapVariant> {
    panic!("unimplented")
}

pub fn create_bam_record(
    hsts: Arc<HstScan>,
    tree_idx: usize,
    top_node_idx: NodeIndex,
    lowest_pvalue: f64,
) -> Record {
    let hst = &hsts.hsts[tree_idx].hst;
    let top_node = hst.node_weight(top_node_idx).unwrap();
    let mut record = record_from_ht(u8_to_hapvariant(hsts.clone(), &top_node.haplotype));
    add_aux_data(&mut record, lowest_pvalue, top_node);
    record
}

fn add_aux_data(record: &mut Record, lowest_pvalue: f64, top_node: &Node) {
    // Add an integer field
    let aux_integer_field = Aux::Float(lowest_pvalue as f32);
    record.push_aux(b"PV", aux_integer_field).unwrap();

    let aux_array: &Vec<u32> = &top_node.indexes.iter().map(|v| *v as u32).collect();
    let aux_array: AuxArray<u32> = aux_array.into();
    let aux_array_field = Aux::ArrayU32(aux_array);
    record.push_aux(b"ID", aux_array_field).unwrap();
}

fn record_from_ht(ht: Vec<HapVariant>) -> Record {
    let start = ht.first().unwrap().pos;
    let stop = ht.last().unwrap().pos;

    let hash_str = ht
        .iter()
        .map(|ht| ht.genotype())
        .fold(format!("{}_{}_{}_", "chr9", start, stop), |acc, e| {
            format!("{acc}{e}")
        });

    let mut hasher = DefaultHasher::new();
    hasher.write(hash_str.as_bytes());

    let qname = format!("{:02x}", hasher.finish());

    let mut cigar = vec![];
    let mut seq = vec![];

    for (i, item) in ht.iter().enumerate() {
        let curr_pos = ht[i].pos;
        let last_pos = ht.get(i - 1).map_or(0, |last| last.pos);

        if last_pos != 0 {
            let distance = curr_pos as u32 - last_pos as u32 - 1;
            cigar.push(Cigar::RefSkip(distance));
        }

        seq.push(item.genotype().as_bytes()[0]);

        match item.gt {
            0 => cigar.push(Cigar::Equal(1)),
            1 => cigar.push(Cigar::Diff(1)),
            _ => unreachable!(),
        }
    }

    let mut record = Record::new();
    record.set_pos(start as i64);
    record.set(
        qname.as_bytes(),
        Some(&CigarString(cigar)),
        &seq,
        &vec![255; seq.len()],
    );
    record
}
