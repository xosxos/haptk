use std::{collections::BTreeMap, path::PathBuf, sync::Arc};
use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;

use petgraph::graph::NodeIndex;
use petgraph::Graph;
use rayon::prelude::*;

use crate::args::{Selection, StandardArgs};
use crate::io::push_to_output;
use crate::structs::Coord;
use crate::subcommands::immutable_hst::construct_bhst_no_mut;
use crate::utils::parse_coords;

use color_eyre::{
    eyre::{ensure, eyre, WrapErr},
    Result,
};
use serde::{Deserialize, Serialize};

use crate::io::get_output;
use crate::read_vcf::read_vcf_to_matrix;
use crate::subcommands::bhst::{HstType, Metadata, Node};

#[doc(hidden)]
pub fn run(args: StandardArgs, step_size: usize) -> Result<()> {
    ensure!(
        args.selection == Selection::All || args.selection == Selection::Haploid,
        "Running only with phased data and all chromsomes is supported."
    );

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "trees", "json.gz");

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None, None)?;

    let hsts = Vec::from_iter(vcf.coords().clone())
        .par_iter()
        // .into_iter()
        .enumerate()
        .filter(|(n, _)| *n % step_size == 0)
        .map(|(_, coord)| {
            (
                coord.clone(),
                construct_bhst_no_mut(&vcf, coord, 4).unwrap(),
            )
        })
        .collect();

    let metadata = Metadata::new(&vcf, &args, vcf.samples().clone(), HstType::Bhst);

    let trees = HstScan { hsts, metadata };

    write_hsts(trees, output)?;

    Ok(())
}

pub type Limits = (usize, usize, usize, usize);
pub type Hst = Graph<Node, ()>;
pub type HstMap = BTreeMap<Coord, Hst>;

#[derive(Serialize, Deserialize, Clone)]
pub struct HstScan {
    pub metadata: Metadata,
    pub hsts: HstMap,
}

impl HstScan {
    pub fn nhaplotypes(&self) -> usize {
        self.metadata.samples.len() * *self.metadata.ploidy
    }

    pub fn get_pos(&self, idx: usize) -> u64 {
        self.metadata.coords.iter().nth(idx).unwrap().pos
    }

    pub fn get_coord_idx(&self, coord: &Coord) -> usize {
        self.metadata
            .coords
            .iter()
            .position(|c| c == coord)
            .unwrap()
    }

    pub fn get_sample_name(&self, index: usize) -> String {
        self.metadata
            .samples
            .get(index / *self.metadata.ploidy)
            .unwrap()
            .clone()
    }

    pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
        indexes.iter().map(|i| self.get_sample_name(*i)).collect()
    }

    pub fn get_sample_idxs(&self, samples: &[String]) -> Result<Vec<usize>> {
        let idxs: Vec<_> = self
            .metadata
            .samples
            .iter()
            .enumerate()
            .filter(|(_, s)| samples.contains(s))
            .flat_map(|(i, _)| {
                ((i * *self.metadata.ploidy)..(i * *self.metadata.ploidy) + *self.metadata.ploidy)
                    .collect::<Vec<usize>>()
            })
            .collect();

        ensure!(
            !idxs.is_empty(),
            "None of the control samples are found in the vcf."
        );
        Ok(idxs)
    }

    pub fn get_sample_idx(&self, sample: &str) -> Result<usize> {
        self.metadata
            .samples
            .iter()
            .position(|s| s == sample)
            .ok_or_else(|| eyre!("sample {sample} not found in the vcf"))
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

    let hst_scan: HstScan = serde_json::from_reader(reader).wrap_err(eyre!(
        "Failed deserializing HSTs from the json.gz. Are you sure the input file is correct?"
    ))?;

    tracing::info!(
        "Read HSTs with the following metadata: \nselection: {:?},\nnsamples:{},",
        hst_scan.metadata.selection,
        hst_scan.metadata.samples.len(),
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

// #[allow(unused_variables)]
// pub fn get_sender<'a>(
//     write_bam: bool,
//     args: &StandardArgs,
//     hsts: Arc<HstScan>,
// ) -> Option<Sender<(&'a Coord, NodeIndex, f64)>> {
//     match write_bam {
//         true => {
//             let (tx, rx) = std::sync::mpsc::channel();

//             std::thread::spawn(move || {
//                 while let Ok((tree_idx, top_node_idx, optimized_value)) = rx.recv() {
//                     todo!()
//                     // let record =
//                     // create_bam_record(hsts.clone(), tree_idx, top_node_idx, optimized_value);
//                     // out.write(&record).unwrap();
//                 }
//             });
//             Some(tx)
//         }
//         false => None,
//     }
// }

pub fn return_assoc<F, U>(
    hsts: &Arc<HstScan>,
    _args: &StandardArgs,
    limits: Limits,
    _write_bam: bool,
    optimizer: F,
    rower: U,
) -> Vec<AssocRow>
where
    F: Fn(Arc<HstScan>, &Hst, &Coord, Limits) -> Option<(Coord, NodeIndex, f64)>
        + Sync
        + Send
        + Clone,
    U: Fn(Arc<HstScan>, &Coord, NodeIndex, f64) -> AssocRow + Sync + Send,
{
    // let tx = get_sender(write_bam, args, hsts.clone());

    hsts.hsts
        .par_iter()
        .filter_map(|(coord, hst)| optimizer(hsts.clone(), hst, coord, limits))
        // .map_with(tx, |tx, optimizer_tuple| {
        //     if let Some(tx) = tx {
        //         tx.send(optimizer_tuple).unwrap();
        //     }
        //     optimizer_tuple
        // })
        .map(|(coord, top_node_idx, optimized_value)| {
            rower(hsts.clone(), &coord, top_node_idx, optimized_value)
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
    pub start_coord: Coord,
    pub stop_coord: Coord,
    pub opt_value: f64,
    pub node_idx: NodeIndex,
    pub first_sample_idx: usize,
    pub hashmap: BTreeMap<String, String>,
}

#[allow(dead_code)]
impl AssocRow {
    pub fn new(
        hsts: Arc<HstScan>,
        coord: &Coord,
        top_node_idx: NodeIndex,
        opt_value: f64,
        hashmap: BTreeMap<String, String>,
    ) -> AssocRow {
        let hst = hsts.hsts.get(coord).unwrap();
        let top_node = hst.node_weight(top_node_idx).unwrap();

        AssocRow {
            contig: coord.contig.clone(),
            pos: coord.pos,
            marker_id: get_marker_id(top_node),
            bp_len: (top_node.stop.pos - top_node.start.pos) as f32,
            marker_len: hsts.get_coord_idx(&top_node.stop) - hsts.get_coord_idx(&top_node.start)
                + 1,
            start: top_node.start.pos,
            stop: top_node.stop.pos,
            start_coord: top_node.start.clone(),
            stop_coord: top_node.stop.clone(),
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
    fn start_coord(&self) -> &Coord {
        &self.start_coord
    }
    fn stop_coord(&self) -> &Coord {
        &self.stop_coord
    }
    fn first_sample_idx(&self) -> usize {
        self.first_sample_idx
    }
    fn opt_value(&self) -> f64 {
        self.opt_value
    }
    fn hashmap(&self) -> &BTreeMap<String, String> {
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
        let mut csv_row = vec![
            row.contig.to_string(),
            row.pos.to_string(),
            row.marker_id.to_string(),
            row.bp_len.to_string(),
            row.marker_len.to_string(),
            row.start.to_string(),
            row.stop.to_string(),
            row.opt_value.to_string(),
        ];

        for v in row.hashmap.values() {
            csv_row.push(v.to_string());
        }

        csv_row
    }
}

impl std::fmt::Display for AssocRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let vector: Vec<String> = self.into();
        let line = format!("{:?}", vector);
        write!(f, "{line}")
    }
}

pub fn get_marker_id(top_node: &Node) -> String {
    let start = top_node.start.pos;
    let stop = top_node.stop.pos;

    let ht = top_node.haplotype.iter().fold(
        format!("{}_{}_{}_", top_node.start.contig, start, stop),
        |acc, e| format!("{acc}{e}"),
    );

    let mut hasher = DefaultHasher::new();
    hasher.write(ht.as_bytes());

    format!("{:02x}", hasher.finish())
}

pub fn top_node_from_hsts<'a>(
    hsts: &'a HstMap,
    coord: &Coord,
    top_node_idx: NodeIndex,
) -> &'a Node {
    let hst = hsts.get(coord).unwrap();
    hst.node_weight(top_node_idx).unwrap()
}

pub fn zygosity_from_node(node: &Node) -> (usize, usize) {
    let (nhet_cases, nhom_cases) = calculate_zygote_n(node.indexes.clone());
    (nhet_cases, nhom_cases)
}

pub fn case_ctrl_zygosity_from_node(
    node: &Node,
    case_list: &[usize],
) -> (usize, usize, usize, usize) {
    let controls = node
        .indexes
        .iter()
        .filter(|i| !case_list.contains(i))
        .copied()
        .collect::<Vec<usize>>();

    let cases = node
        .indexes
        .iter()
        .filter(|i| case_list.contains(i))
        .copied()
        .collect::<Vec<usize>>();

    let (nhet_ctrls, nhom_ctrls) = calculate_zygote_n(controls);
    let (nhet_cases, nhom_cases) = calculate_zygote_n(cases);

    (nhet_cases, nhom_cases, nhet_ctrls, nhom_ctrls)
}

fn calculate_zygote_n(mut indexes: Vec<usize>) -> (usize, usize) {
    let prior_len = indexes.len();
    indexes.sort();
    indexes.dedup();
    let post_len = indexes.len();

    let nhomozygotes = prior_len - post_len;
    let nheterozygotes = post_len - nhomozygotes;

    (nheterozygotes, nhomozygotes)
}
