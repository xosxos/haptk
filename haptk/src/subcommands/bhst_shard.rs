use std::{borrow::Cow, collections::BTreeSet, path::PathBuf};

use color_eyre::{
    eyre::{ensure, eyre, Context},
    Result,
};
use petgraph::prelude::NodeIndex;
use petgraph::{Direction, Graph};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{
    args::{Selection, StandardArgs},
    io::{get_output, open_csv_writer, push_to_output, write_haplotype},
    libs::read_vcf::{get_sample_names, read_vcf_to_matrix},
    libs::structs::{CoordDataSlot, PhasedMatrix},
    structs::{Coord, HapVariant, Ploidy},
    utils::parse_snp_coord,
};

pub fn read_vcf_with_selections(args: &StandardArgs) -> Result<PhasedMatrix> {
    let (contig, pos) = parse_snp_coord(&args.coords)?;

    let (indexes, _) = get_sample_names(args, contig, None)?;
    ensure!(
        indexes.len() > 1,
        "Cannot build a tree with less than 2 samples."
    );
    let mut vcf = read_vcf_to_matrix(args, contig, pos, None, None)?;

    match &args.selection {
        Selection::OnlyAlts | Selection::OnlyRefs => vcf.select_carriers(pos, &args.selection)?,
        Selection::OnlyLongest => vcf.select_only_longest(),
        _ => (),
    }

    Ok(vcf)
}

#[doc(hidden)]
pub fn run(args: StandardArgs, min_size: usize, publish: bool) -> Result<()> {
    ensure!(
        args.selection != Selection::Unphased,
        "Running with unphased data is not supported."
    );

    let vcf = read_vcf_with_selections(&args)?;

    ensure!(
        vcf.nhaplotypes() >= min_size,
        "VCF has less haplotypes than the required minimum node size ({} < {min_size})",
        vcf.nhaplotypes()
    );

    // Construct the bilateral HST
    let bhst = construct_bhst(&vcf, vcf.start_coord(), min_size);
    tracing::info!("Finished HST construction.");

    // Majority based haplotype
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "bhst_mbah", "csv");
    let writer = open_csv_writer(sh_output)?;
    let mbah = find_mbah(&bhst, &vcf)?;
    write_haplotype(mbah, writer)?;

    // Shared haplotype
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "bhst_shared_core_haplotype", "csv");
    let writer = open_csv_writer(sh_output)?;
    let ht = find_shared_haplotype(&bhst, &vcf);
    write_haplotype(ht, writer)?;

    // Write to .hst
    let mut hst_output = args.output.clone();
    push_to_output(&args, &mut hst_output, "bhst", "hst.gz");
    write_hst_file(bhst, &vcf, hst_output, publish, args, HstType::Bhst)?;

    Ok(())
}

#[derive(Serialize, Deserialize, Default, Clone, Debug, PartialEq, PartialOrd)]
pub struct Node<'matrix> {
    pub indexes: Vec<usize>,
    pub start: Cow<'matrix, Coord>,
    pub stop: Cow<'matrix, Coord>,
    pub haplotype: Vec<u8>,
}

impl<'matrix> std::fmt::Display for Node<'matrix> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // let nmarkers = match self.stop_idx.cmp(&self.start_idx) {
        // std::cmp::Ordering::Greater => self.stop_idx - self.start_idx + 1,
        // std::cmp::Ordering::Less => self.start_idx - self.stop_idx + 1,
        // std::cmp::Ordering::Equal => 0,
        // };
        // let line = format!(
        //     "n: {}\nstart_idx: {}\nstop_idx: {}\nnmarkers: {}\nblock length (bp): {}\n",
        //     self.indexes.len(),
        //     self.start_idx,
        //     self.stop_idx,
        //     nmarkers,
        //     self.block_len,
        // );
        let line = format!("n: {}\n", self.indexes.len(),);
        write!(f, "{line}")
    }
}

#[doc(hidden)]
pub fn construct_bhst<'matrix>(
    vcf: &'matrix PhasedMatrix,
    start_coord: &'matrix Coord,
    min_size: usize,
) -> Graph<Node<'matrix>, u8> {
    let mut hst = initiate_hst(vcf, start_coord);

    let mut min_size_blacklist = vec![];
    loop {
        // From all nodes find all leaf nodes with more indexes than min_size
        let indices = find_leaf_nodes(&hst, min_size, &min_size_blacklist);

        // Iterate through the filtered leaf nodes in parallel to find new nodes to add
        let nodes = indices
            .into_par_iter()
            .filter_map(|node_idx| {
                find_contradictory_gt(vcf, &hst, node_idx).map(|nodes| (node_idx, nodes))
            })
            .collect::<Vec<(NodeIndex, Vec<Node>)>>();

        if nodes.is_empty() {
            break;
        }

        // Add new nodes to the tree
        for (node_idx, new_nodes) in nodes {
            let hst_node_count_before = hst.node_count();

            for new_node in new_nodes.into_iter() {
                if new_node.indexes.len() >= min_size {
                    let new_node_idx = hst.add_node(new_node);

                    hst.add_edge(node_idx, new_node_idx, 0);
                }
            }

            if hst.node_count() == hst_node_count_before {
                min_size_blacklist.push(node_idx);
            }
        }
    }
    hst
}

pub fn initiate_hst<'matrix>(
    vcf: &'matrix PhasedMatrix,
    start_coord: &'matrix Coord,
) -> Graph<Node<'matrix>, u8> {
    let mut hst = Graph::<_, _>::new();

    let indexes: Vec<usize> = (0..vcf.nhaplotypes()).collect();

    hst.add_node(Node {
        start: Cow::Borrowed(start_coord),
        stop: Cow::Borrowed(start_coord),
        indexes: indexes.clone(),
        haplotype: vec![],
    });

    // If the first variant is already contradictory, split into two nodes
    if vcf.is_contradictory(start_coord, &indexes) {
        let genotypes = vcf.get_slot(start_coord);
        let (mut ones, mut zeroes) = (vec![], vec![]);

        for i in indexes.iter() {
            match genotypes[*i] == 1 {
                true => ones.push(*i),
                false => zeroes.push(*i),
            }
        }

        let node1 = Node {
            start: Cow::Borrowed(start_coord),
            stop: Cow::Borrowed(start_coord),
            haplotype: vcf.find_u8_haplotype_for_sample(start_coord..=start_coord, ones[0]),
            indexes: ones,
        };

        let node2 = Node {
            start: Cow::Borrowed(start_coord),
            stop: Cow::Borrowed(start_coord),
            haplotype: vcf.find_u8_haplotype_for_sample(start_coord..=start_coord, zeroes[0]),
            indexes: zeroes,
        };

        let node1 = hst.add_node(node1);
        let node2 = hst.add_node(node2);

        hst.add_edge(NodeIndex::new(0), node1, 0);
        hst.add_edge(NodeIndex::new(0), node2, 0);
    }

    hst
}

pub fn find_leaf_nodes(
    hst: &Graph<Node, u8>,
    min_size: usize,
    min_size_blacklist: &[NodeIndex],
) -> Vec<NodeIndex> {
    hst.node_indices()
        .filter(|node_idx| !min_size_blacklist.contains(node_idx))
        .filter(|node_idx| {
            let node = hst.node_weight(*node_idx).unwrap();

            let count = hst
                .neighbors_directed(*node_idx, Direction::Outgoing)
                .count();

            // Children count is 0, but there are still over min_size samples left
            // Collect such nodes into the `indices` vector
            // min_size.max(2) means use min_size if it's higher than 2, if lower, use 2
            count == 0 && node.indexes.len() >= min_size.max(2)
        })
        .collect()
}

#[doc(hidden)]
fn find_contradictory_gt<'matrix>(
    vcf: &'matrix PhasedMatrix,
    bhst: &Graph<Node<'matrix>, u8>,
    node_idx: NodeIndex,
) -> Option<Vec<Node<'matrix>>> {
    let node = bhst.node_weight(node_idx).unwrap();
    let (left_coord, right_coord) = (node.start.clone(), node.stop.clone());

    // Use the CoordDataSlot trait from structs.rs to access the genotype matrix
    // Accessing the matrix has to be behind an interface to allow for the back end logic
    // to be switched out in the future
    let prev = vcf.prev_contradictory(&left_coord, &node.indexes);
    let next = vcf.next_contradictory(&right_coord, &node.indexes);

    // Allocate Vecs with capacity so no reallocation is required
    // NOTE: Benchmarking required
    let (mut zo, mut zz, mut oz, mut oo) = (
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
    );

    // Fill the genotype buckets and then create HST nodes from the buckets
    match (prev, next) {
        (Some(left), Some(right)) => {
            let left_vec = vcf.get_slot(left);
            let right_vec = vcf.get_slot(right);

            for i in node.indexes.iter() {
                let left_bit = left_vec[*i] == 1;
                let right_bit = right_vec[*i] == 1;

                match (left_bit, right_bit) {
                    (false, false) => zz.push(*i),
                    (false, true) => zo.push(*i),
                    (true, false) => oz.push(*i),
                    (true, true) => oo.push(*i),
                }
            }
            // Create a new node for each bucket if it is not empty
            let nodes = create_nodes_from_buckets(vcf, left, right, oo, oz, zo, zz);
            Some(nodes)
        }
        (Some(left), None) => {
            tracing::warn!(
                "Genotyping data ran out on the right with samples {:?}",
                vcf.get_sample_names(&node.indexes)
            );

            let left_vec = vcf.get_slot(left);

            for i in node.indexes.iter() {
                match left_vec[*i] == 1 {
                    true => oz.push(*i),
                    false => zz.push(*i),
                }
            }
            let nodes =
                create_nodes_from_buckets(vcf, left, vcf.coords().last().unwrap(), oo, oz, zo, zz);
            Some(nodes)
        }
        (None, Some(right)) => {
            tracing::warn!(
                "Genotyping data ran out on the left with samples {:?}",
                vcf.get_sample_names(&node.indexes)
            );

            let right_vec = vcf.get_slot(right);

            for i in node.indexes.iter() {
                match right_vec[*i] == 1 {
                    true => zo.push(*i),
                    false => zz.push(*i),
                }
            }
            let nodes = create_nodes_from_buckets(
                vcf,
                vcf.coords().first().unwrap(),
                right,
                oo,
                oz,
                zo,
                zz,
            );
            Some(nodes)
        }
        (None, None) => None,
    }
}

fn create_nodes_from_buckets<'matrix>(
    vcf: &'matrix PhasedMatrix,
    left: &'matrix Coord,
    right: &'matrix Coord,
    oo: Vec<usize>,
    oz: Vec<usize>,
    zo: Vec<usize>,
    zz: Vec<usize>,
) -> Vec<Node<'matrix>> {
    let mut nodes = Vec::with_capacity(4);
    if !zz.is_empty() {
        let node = Node {
            start: Cow::Borrowed(left),
            stop: Cow::Borrowed(right),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, zz[0]),
            indexes: zz,
        };
        nodes.push(node);
    }

    if !zo.is_empty() {
        let node = Node {
            start: Cow::Borrowed(left),
            stop: Cow::Borrowed(right),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, zo[0]),
            indexes: zo,
        };
        nodes.push(node);
    }

    if !oz.is_empty() {
        let node = Node {
            start: Cow::Borrowed(left),
            stop: Cow::Borrowed(right),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, oz[0]),
            indexes: oz,
        };
        nodes.push(node);
    }

    if !oo.is_empty() {
        let node = Node {
            start: Cow::Borrowed(left),
            stop: Cow::Borrowed(right),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, oo[0]),
            indexes: oo,
        };
        nodes.push(node);
    }
    nodes
}

// HST UTILS

pub fn find_mbah(g: &Graph<Node, u8>, vcf: &PhasedMatrix) -> Result<Vec<HapVariant>> {
    let start_idx = NodeIndex::new(0);
    let nodes = find_majority_nodes(g, start_idx);

    ensure!(
        nodes.len() > 2,
        "The majority branch has less than 3 nodes."
    );

    let mut iter = nodes.iter().rev();
    let (last_node, _last_node_idx) = iter.next().unwrap();
    let (second_last_node, _) = iter.next().unwrap();

    tracing::info!(
        "Majority-based ancestral samples:{}",
        second_last_node
            .indexes
            .iter()
            .fold(String::new(), |cur, acc| {
                let acc = vcf.get_sample_name(*acc);
                format!("{cur} {acc}")
            }),
    );

    let idx = vcf.get_coord_idx(&last_node.start);

    // Start and stop idxs are the first deviating genotype and it cannot be a part of the haplotype
    // last_node.start_idx..last_node.stop_idx + 1,

    let start = vcf.get_coord(idx + 1);

    let start = match start > &last_node.stop {
        true => last_node.start.as_ref(),
        false => start,
    };

    Ok(vcf.find_haplotype_for_sample(start..last_node.stop.as_ref(), last_node.indexes[0]))
}

pub fn find_shared_haplotype(g: &Graph<Node, u8>, vcf: &PhasedMatrix) -> Vec<HapVariant> {
    let start_idx = NodeIndex::new(0);
    let nodes = find_majority_nodes(g, start_idx);
    let (second_node, _second_node_idx) = nodes[1];

    let idx = vcf.get_coord_idx(&second_node.start);
    let start = vcf.get_coord(idx + 1);

    let start = match start > &second_node.stop {
        true => second_node.start.as_ref(),
        false => start,
    };

    vcf.find_haplotype_for_sample(start..second_node.stop.as_ref(), second_node.indexes[0])
}

pub fn calculate_block_len(node: &Node) -> u64 {
    match node.stop.cmp(&node.start) {
        std::cmp::Ordering::Greater => node.stop.pos - node.start.pos,
        std::cmp::Ordering::Less => node.start.pos - node.stop.pos,
        std::cmp::Ordering::Equal => 0,
    }
}

pub fn find_majority_nodes<'a, 'b>(
    g: &'b Graph<Node<'a>, u8>,
    start_idx: NodeIndex,
) -> Vec<(&'b Node<'a>, NodeIndex)> {
    let mut idx = start_idx;
    let start_node = g.node_weight(idx).unwrap();
    let mut majority_nodes = vec![(start_node, idx)];

    loop {
        let nodes = g
            .neighbors_directed(idx, Direction::Outgoing)
            .map(|idx| (g.node_weight(idx).unwrap(), idx))
            .collect::<Vec<(&_, NodeIndex)>>();

        if nodes.is_empty() {
            break;
        }

        if let Some(max) = nodes
            .iter()
            .max_by(|x, y| x.0.indexes.len().cmp(&y.0.indexes.len()))
        {
            majority_nodes.push(*max);
            idx = max.1;
        } else {
            // Max by returns None if the iterator is empty, thus when a leaf node is reached
            break;
        }
    }

    majority_nodes
}

pub fn find_lowest_start(g: &Graph<Node, u8>) -> Coord {
    let mut start = &Coord {
        pos: u64::MAX,
        contig: String::from(""),
        reference: String::from(""),
        alt: String::from(""),
    };

    for node in g.node_weights() {
        if node.start.pos < start.pos {
            start = &node.start;
        }
    }
    assert_ne!(start.pos, u64::MAX);
    start.clone()
}

pub fn find_highest_stop(g: &Graph<Node, u8>) -> Coord {
    let mut stop = &Coord {
        pos: u64::MIN,
        contig: String::from(""),
        reference: String::from(""),
        alt: String::from(""),
    };

    for node in g.node_weights() {
        if node.stop.pos > stop.pos {
            stop = &node.stop;
        }
    }
    assert_ne!(stop.pos, u64::MIN);
    stop.clone()
}

pub fn find_coord_list(g: &Graph<Node, u8>, vcf: &PhasedMatrix) -> Vec<Coord> {
    let start_idx = find_lowest_start(g);
    let stop_idx = find_highest_stop(g);

    vcf.coords().range(start_idx..=stop_idx).cloned().collect()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hst<'a> {
    // pub coords: Vec<Coord>,
    pub hst: Graph<Node<'a>, u8>,
    // pub samples: Vec<String>,
    // #[serde(default)]
    pub metadata: Metadata,
}

impl<'a> Hst<'a> {
    pub fn get_haplotype(&self, node_idx: NodeIndex) -> Vec<HapVariant> {
        let node = self.hst.node_weight(node_idx).unwrap();

        if node.haplotype.is_empty() && node.start.pos == node.stop.pos {
            return vec![];
        }

        self.metadata
            .coords
            .range(node.start.as_ref()..=node.stop.as_ref())
            .enumerate()
            .map(|(hap_index, coord)| HapVariant {
                contig: coord.contig.to_string(),
                pos: coord.pos,
                alt: coord.alt.clone(),
                reference: coord.reference.clone(),
                gt: node.haplotype[hap_index],
            })
            .collect()
    }

    // pub fn get_pos(&self, idx: usize) -> u64 {
    // self.metadata.coords[idx].pos
    // }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum HstType {
    UhstLeft,
    UhstRight,
    Bhst,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Metadata {
    pub start_coord: Coord,
    pub hst_type: HstType,
    pub coords: BTreeSet<Coord>,
    pub contig: String,
    pub samples: Vec<String>,
    pub selection: Selection,
    pub ploidy: Ploidy,
    pub vcf_name: PathBuf,
}

impl Metadata {
    pub fn new(
        vcf: &PhasedMatrix,
        args: &StandardArgs,
        samples: Vec<String>,
        hst_type: HstType,
    ) -> Self {
        Self {
            start_coord: vcf.start_coord().clone(),
            coords: vcf.coords().clone(),
            contig: vcf.get_contig().clone(),
            samples,
            hst_type,
            selection: args.selection.clone(),
            ploidy: vcf.ploidy.clone(),
            vcf_name: args.file.clone(),
        }
    }
}

pub fn write_hst_file(
    mut hst: Graph<Node, u8>,
    vcf: &PhasedMatrix,
    path: PathBuf,
    publish: bool,
    args: StandardArgs,
    hst_type: HstType,
) -> Result<()> {
    let samples = match publish {
        true => vec!["r".to_string()],
        false => vcf.samples().clone(),
    };

    if publish {
        for node in hst.node_weights_mut() {
            node.indexes = vec![0; node.indexes.len()];
        }
    }

    let hst_export = Hst {
        // coords: vcf.coords().clone(),
        hst,
        // samples,
        metadata: Metadata::new(vcf, &args, samples, hst_type),
    };

    tracing::info!("HST output: {path:?}.");
    let now = std::time::Instant::now();
    let mut output = get_output(Some(path))?;

    let mut writer = bgzip::BGZFWriter::new(&mut output, bgzip::Compression::default());

    serde_json::to_writer(&mut writer, &hst_export)?;

    writer.close()?;

    tracing::info!("Wrote HSTs in {:?}", now.elapsed());
    Ok(())
}

pub fn read_hst_file<'a>(path: PathBuf) -> Result<Hst<'a>> {
    let file = std::fs::File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;
    let reader = bgzip::BGZFReader::new(file)?;
    let hst: Hst = serde_json::from_reader(reader)?;

    Ok(hst)
}
