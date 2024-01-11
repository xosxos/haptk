use std::path::PathBuf;

use color_eyre::{
    eyre::{ensure, eyre, Context},
    Result,
};
use ndarray::s;
use petgraph::prelude::NodeIndex;
use petgraph::{Direction, Graph};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{
    args::{GraphArgs, Selection, StandardArgs},
    graphs::HstGraph,
    io::{get_output, open_csv_writer, push_to_output, read_sample_ids, read_variable_data_file},
    read_vcf::{get_sample_names, read_vcf_to_matrix},
    structs::{Coord, HapVariant, PhasedMatrix, CoordDataSlot},
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
pub fn run(
    args: StandardArgs,
    graph_args: GraphArgs,
    decoy_samples: Option<PathBuf>,
    variable_data: Option<PathBuf>,
    variable_names: Option<Vec<String>>,
    min_size: usize,
    publish: bool,
) -> Result<()> {
    ensure!(
        args.selection != Selection::Unphased,
        "Running with unphased data is not supported."
    );

    let decoy_samples = read_sample_ids(&decoy_samples)?;

    let mut vcf = read_vcf_with_selections(&args)?;

    // Set clinical data
    if let Some(path) = variable_data {
        vcf.set_variable_data(read_variable_data_file(path)?)?;
    }

    ensure!(
        vcf.nrows() >= min_size,
        "VCF has less haplotypes than the required minimum node size ({} < {min_size})",
        vcf.nrows()
    );

    // Construct the bilateral HST
    let bhst = construct_bhst(&vcf, vcf.variant_idx(), min_size);
    tracing::info!("Finished HST construction.");

    // Write to svg
    let mut dg = HstGraph::new(&bhst, &vcf, graph_args, decoy_samples, min_size);
    dg.draw_graph()?;
    let mut decay_graph = args.output.clone();
    push_to_output(&args, &mut decay_graph, "bhst", "svg");
    svg::save(decay_graph, &dg.document)?;

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

    // Calculate variable data
    let mut vars_output = args.output.clone();
    push_to_output(&args, &mut vars_output, "bhst_vars", "csv");
    calculate_and_write_var_data(&bhst, &vcf, variable_names, vars_output)?;

    // Write to .hst
    let mut hst_output = args.output.clone();
    push_to_output(&args, &mut hst_output, "bhst", "hst.gz");
    write_hst_file(bhst, &vcf, hst_output, publish, args)?;

    Ok(())
}

#[derive(Serialize, Deserialize, Default, Clone, Debug, PartialEq, PartialOrd)]
pub struct Node {
    pub indexes: Vec<usize>,
    pub start_idx: usize,
    pub stop_idx: usize,
    pub haplotype: Vec<u8>,
}

impl std::fmt::Display for Node {
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
pub fn construct_bhst(vcf: &PhasedMatrix, idx: usize, min_size: usize) -> Graph<Node, u8> {
    let mut bhst = Graph::<_, _>::new();

    bhst.add_node(Node {
        start_idx: idx,
        stop_idx: idx,
        indexes: (0..vcf.matrix.nrows()).collect(),
        haplotype: vec![],
    });

    let mut min_size_blacklist = vec![];
    loop {
        // Filter out nodes with children and nodes with less indexes than min_size
        let indices = bhst
            .node_indices()
            .filter(|node_idx| !min_size_blacklist.contains(node_idx))
            .filter(|node_idx| {
                let node = bhst.node_weight(*node_idx).unwrap();
                let count = bhst
                    .neighbors_directed(*node_idx, Direction::Outgoing)
                    .count();

                // Children count is 0, but there are still over min_size samples left
                // Collect such nodes into the `indices` vector
                count == 0 && node.indexes.len() >= min_size.max(2)
            })
            .collect::<Vec<NodeIndex>>();

        // Multithread horizontally all childless nodes
        let nodes = indices
            .into_par_iter()
            .filter_map(|node_idx| {
                find_contradictory_gt(vcf, &bhst, node_idx).map(|nodes| (node_idx, nodes))
            })
            .collect::<Vec<(NodeIndex, Vec<Node>)>>();

        // Terminate if no new nodes will be added to the tree
        if nodes.is_empty() {
            break;
        }

        // Add nodes to the tree
        for (node_idx, new_nodes) in nodes {
            let bhst_node_count_before = bhst.node_count();

            for new_node in new_nodes {
                if new_node.indexes.len() >= min_size {
                    let new_node_idx = bhst.add_node(new_node.clone());

                    bhst.add_edge(node_idx, new_node_idx, 0);
                }
            }

            if bhst.node_count() == bhst_node_count_before {
                min_size_blacklist.push(node_idx);
            }
        }
    }
    bhst
}

#[doc(hidden)]
fn find_contradictory_gt(
    vcf: &PhasedMatrix,
    bhst: &Graph<Node, u8>,
    node_idx: NodeIndex,
) -> Option<Vec<Node>> {
    let node = bhst.node_weight(node_idx).unwrap();
    let (left_idx, mut right_idx) = (node.start_idx, node.stop_idx);
    // Minus 1 to account for the starting variant itself as well
    if node_idx == NodeIndex::new(0) {
        right_idx = right_idx.saturating_sub(1);
    }

    // Allocate Vecs with capacity so no reallocation is required
    let prev = vcf.prev_contradictory(left_idx, &node.indexes);
    let next = vcf.next_contradictory(right_idx, &node.indexes);

    let (mut zo, mut zz, mut oz, mut oo) = (
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
    );
    match (prev, next) {
        (Some(left), Some(right)) => {
            let left_vec =vcf.get_slot(left);
            let right_vec =vcf.get_slot(right);
            for i in node.indexes.iter() {
                let left_bit = left_vec[*i] == 1;
                let right_bit = right_vec[*i] == 1;
                match (left_bit, right_bit) {
                    (false, false) => zz.push(*i),
                    (false, true)  => zo.push(*i),
                    (true,  false) => oz.push(*i),
                    (true,  true)  => oo.push(*i),
                }
            }
            // Create a new node for each bucket if it is not empty
            let nodes = create_nodes_from_buckets(vcf, left, right, oo, oz, zo, zz);
            Some(nodes)
        }
        (_, _) => None
    }
}

fn create_nodes_from_buckets(
    vcf: &PhasedMatrix,
    left: usize,
    right: usize,
    oo: Vec<usize>,
    oz: Vec<usize>,
    zo: Vec<usize>,
    zz: Vec<usize>,
) -> Vec<Node> {
    let mut nodes = Vec::with_capacity(4);
    if !zz.is_empty() {
        let node = Node {
            start_idx: left,
            stop_idx: right,
            haplotype: vcf.find_u8_haplotype_for_sample(left..right + 1, zz[0]),
            indexes: zz,
        };
        nodes.push(node);
    }

    if !zo.is_empty() {
        let node = Node {
            start_idx: left,
            stop_idx: right,
            haplotype: vcf.find_u8_haplotype_for_sample(left..right + 1, zo[0]),
            indexes: zo,
        };
        nodes.push(node);
    }

    if !oz.is_empty() {
        let node = Node {
            start_idx: left,
            stop_idx: right,
            haplotype: vcf.find_u8_haplotype_for_sample(left..right + 1, oz[0]),
            indexes: oz,
        };
        nodes.push(node);
    }

    if !oo.is_empty() {
        let node = Node {
            start_idx: left,
            stop_idx: right,
            haplotype: vcf.find_u8_haplotype_for_sample(left..right + 1, oo[0]),
            indexes: oo,
        };
        nodes.push(node);
    }
    nodes
}

pub fn find_mbah(g: &Graph<Node, u8>, vcf: &PhasedMatrix) -> Result<Vec<HapVariant>> {
    let nodes = find_majority_nodes(g);

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

    Ok(vcf.find_haplotype_for_sample(
        vcf.get_contig(),
        // Start and stop idxs are the first deviating genotype and it cannot be a part of the haplotype
        // last_node.start_idx..last_node.stop_idx + 1,
        last_node.start_idx + 1..last_node.stop_idx,
        last_node.indexes[0],
    ))
}

pub fn find_shared_haplotype(g: &Graph<Node, u8>, vcf: &PhasedMatrix) -> Vec<HapVariant> {
    let nodes = find_majority_nodes(g);
    let (second_node, _second_node_idx) = nodes[1];

    vcf.find_haplotype_for_sample(
        vcf.get_contig(),
        second_node.start_idx + 1..second_node.stop_idx,
        second_node.indexes[0],
    )
}

// HST UTILS

pub fn calculate_block_len(vcf: &PhasedMatrix, node: &Node) -> u64 {
    match node.stop_idx.cmp(&node.start_idx) {
        std::cmp::Ordering::Greater => vcf.get_pos(node.stop_idx) - vcf.get_pos(node.start_idx),
        std::cmp::Ordering::Less => vcf.get_pos(node.start_idx) - vcf.get_pos(node.stop_idx),
        std::cmp::Ordering::Equal => 0,
    }
}

pub fn find_majority_nodes(g: &Graph<Node, u8>) -> Vec<(&Node, NodeIndex)> {
    let mut idx = NodeIndex::new(0);
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

pub fn calculate_and_write_var_data(
    g: &Graph<Node, u8>,
    vcf: &PhasedMatrix,
    variable_names: Option<Vec<String>>,
    vars_output: PathBuf,
) -> Result<()> {
    // Calculate variable data
    if let Some(vars) = &variable_names {
        let mut writer = open_csv_writer(vars_output)?;

        let mut header: Vec<String> = ["side", "pos", "1st_n", "2nd_n"]
            .iter()
            .map(|s| s.to_string())
            .collect();

        for var in vars {
            header.push(format!("{var}_1st"));
            header.push(format!("{var}_2nd"));
            header.push(format!("{var}_t_test"));
        }

        writer.write_record(header)?;

        let mnodes = find_majority_nodes(g);
        for (_node, idx) in mnodes {
            let children: Vec<NodeIndex> = g.neighbors_directed(idx, Direction::Outgoing).collect();

            let mut children: Vec<(&_, NodeIndex)> = children
                .iter()
                .map(|c| (g.node_weight(*c).unwrap(), *c))
                .collect();

            children.sort_by(|a, b| b.0.indexes.len().cmp(&a.0.indexes.len()));

            for (i, (child_node, _child_idx)) in children.iter().enumerate() {
                if i > 0 && (children[0].0.indexes.len() > 5 && child_node.indexes.len() > 5) {
                    let majority = children[0].0;
                    let sibling_node = child_node;

                    let v1 = vcf
                        .get_variable_data_vecs(&majority.indexes, vars)?
                        .ok_or(eyre!("Error in variable data."))?;
                    let v2 = vcf
                        .get_variable_data_vecs(&sibling_node.indexes, vars)?
                        .ok_or(eyre!("Error in variable data."))?;

                    let n1 = vcf
                        .get_variable_data_mean(&majority.indexes, vars)?
                        .ok_or(eyre!("Cannot calculate variable data means."))?;
                    let n2 = vcf
                        .get_variable_data_mean(&sibling_node.indexes, vars)?
                        .ok_or(eyre!("Cannot calculate variable data means."))?;

                    let block_len = calculate_block_len(vcf, majority);

                    let mut csv_row = vec![
                        block_len.to_string(),
                        majority.indexes.len().to_string(),
                        sibling_node.indexes.len().to_string(),
                    ];
                    for i in 0..n1.len() {
                        csv_row.push(format!("{:.2}", n1[i]));
                        csv_row.push(format!("{:.2}", n2[i]));
                        let t_test = crate::stats::two_tail_welch_t_test(&v1[i], &v2[i]);
                        csv_row.push(format!("{t_test:.5}"));
                    }
                    writer.write_record(csv_row)?;
                }
            }
        }
    }
    Ok(())
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

pub fn find_lowest_start_idx(g: &Graph<Node, u8>) -> usize {
    let mut start_idx = usize::MAX;

    for node in g.node_weights() {
        if node.start_idx < start_idx {
            start_idx = node.start_idx;
        }
    }
    assert_ne!(start_idx, usize::MAX);
    start_idx
}

pub fn find_highest_stop_idx(g: &Graph<Node, u8>) -> usize {
    let mut stop_idx = usize::MIN;

    for node in g.node_weights() {
        if node.stop_idx > stop_idx {
            stop_idx = node.stop_idx;
        }
    }
    assert_ne!(stop_idx, usize::MIN);
    stop_idx
}

pub fn find_coord_list(g: &Graph<Node, u8>, vcf: &PhasedMatrix) -> Vec<Coord> {
    let start_idx = find_lowest_start_idx(g);
    let stop_idx = find_highest_stop_idx(g);
    let mut coords = vec![];

    for i in start_idx..=stop_idx {
        coords.push(vcf.get_coord(i).clone());
    }
    coords
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    pub start_coord: String,
    pub selection: Selection,
    pub vcf_name: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hst {
    pub coords: Vec<Coord>,
    pub hst: Graph<Node, u8>,
    pub samples: Vec<String>,
    #[serde(default)]
    pub metadata: Metadata,
}

impl Hst {
    pub fn get_haplotype(&self, node_idx: NodeIndex) -> Vec<HapVariant> {
        let node = self.hst.node_weight(node_idx).unwrap();

        if node.start_idx == node.stop_idx {
            return vec![];
        }

        (node.start_idx..node.stop_idx + 1)
            .enumerate()
            .map(|(hap_index, coord_index)| HapVariant {
                contig: self.coords[0].contig.to_string(),
                pos: self.coords[coord_index].pos,
                alt: self.coords[coord_index].alt.clone(),
                reference: self.coords[coord_index].reference.clone(),
                gt: node.haplotype[hap_index],
            })
            .collect()
    }

    pub fn get_pos(&self, idx: usize) -> u64 {
        self.coords[idx].pos
    }
}

pub fn write_hst_file(
    mut hst: Graph<Node, u8>,
    vcf: &PhasedMatrix,
    path: PathBuf,
    publish: bool,
    args: StandardArgs,
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
        coords: vcf.coords().clone(),
        hst,
        samples,
        metadata: Metadata {
            start_coord: args.coords,
            selection: args.selection,
            vcf_name: args.file,
        },
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

pub fn read_hst_file(path: PathBuf) -> Result<Hst> {
    let file = std::fs::File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;
    let reader = bgzip::BGZFReader::new(file)?;
    let hst: Hst = serde_json::from_reader(reader)?;

    Ok(hst)
}
