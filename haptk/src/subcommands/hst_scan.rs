use std::path::PathBuf;

use color_eyre::{
    eyre::{ensure, eyre, WrapErr},
    Result,
};
use ndarray::s;
use petgraph::prelude::NodeIndex;
use petgraph::{Direction, Graph};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::args::{Selection, StandardArgs};
use crate::core::get_output;
use crate::core::parse_coords;
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::PhasedMatrix;
use crate::subcommands::bhst;
use crate::utils::push_to_output;

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

    let trees = range
        .par_iter()
        .filter(|idx| *idx % step_size == 0)
        .map(|idx| {
            if idx % (vcf.ncoords() / 10) == 0 {
                tracing::info!("Finished constructing 10% of the trees...");
            }

            TreeRow {
                idx: *idx,
                tree: construct_bhst(&vcf, *idx, 4),
            }
        })
        .collect();

    let metadata = Metadata::new(&vcf, args);

    let trees = Trees { trees, metadata };

    write_trees(trees, output)?;

    Ok(())
}
#[derive(Serialize, Deserialize, Default, Clone, Debug, PartialEq, PartialOrd)]
pub struct Node {
    pub indexes: Vec<usize>,
    pub start_idx: usize,
    pub stop_idx: usize,
    pub block_len: f32,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct TreeRow {
    pub idx: usize,
    pub tree: Graph<Node, u8>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Metadata {
    pub coords: String,
    pub samples: Vec<String>,
    pub selection: Selection,
    pub vcf_name: PathBuf,
    pub info_limit: Option<f32>,
}

impl Metadata {
    fn new(vcf: &PhasedMatrix, args: StandardArgs) -> Self {
        let samples = vcf.samples().to_vec();

        Self {
            coords: args.coords,
            samples,
            selection: args.selection,
            vcf_name: args.file,
            info_limit: args.info_limit,
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Trees {
    pub metadata: Metadata,
    pub trees: Vec<TreeRow>,
}

fn write_trees(trees: Trees, path: PathBuf) -> Result<()> {
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

pub fn read_tree_file(path: PathBuf) -> Result<Trees> {
    let file = std::fs::File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;
    let reader = bgzip::BGZFReader::new(file)?;
    let trees: Trees = serde_json::from_reader(reader)?;
    tracing::info!("Read HSTs with the following metadata: \ncoords: {},\nselection: {:?},\nnsamples:{},\ninfo_limit: {:?}", trees.metadata.coords, trees.metadata.selection, trees.metadata.samples.len(), trees.metadata.info_limit);

    Ok(trees)
}

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let nmarkers = match self.stop_idx.cmp(&self.start_idx) {
            std::cmp::Ordering::Greater => self.stop_idx - self.start_idx + 1,
            std::cmp::Ordering::Less => self.start_idx - self.stop_idx + 1,
            std::cmp::Ordering::Equal => 0,
        };
        let line = format!(
            "n: {}\nstart_idx: {}\nstop_idx: {}\nnmarkers: {}\nblock length (bp): {}\n",
            self.indexes.len(),
            self.start_idx,
            self.stop_idx,
            nmarkers,
            self.block_len,
        );
        write!(f, "{line}")
    }
}

impl bhst::Indexes for Node {
    fn indexes(&self) -> &Vec<usize> {
        &self.indexes
    }
    fn block_len(&self) -> f32 {
        self.block_len
    }

    fn show(&self) -> String {
        self.to_string()
    }
}

#[doc(hidden)]
pub fn construct_bhst(vcf: &PhasedMatrix, idx: usize, min_size: usize) -> Graph<Node, u8> {
    let mut bhst = Graph::<_, _>::new();

    bhst.add_node(Node {
        block_len: 0.0,
        start_idx: idx,
        stop_idx: idx,
        indexes: (0..vcf.matrix.nrows()).collect(),
    });

    loop {
        // Filter out nodes with children and nodes with less indexes than min_size
        let indices = bhst
            .node_indices()
            .filter(|node_idx| {
                let node = bhst.node_weight(*node_idx).unwrap();
                let count = bhst
                    .neighbors_directed(*node_idx, Direction::Outgoing)
                    .count();

                count == 0 && node.indexes.len() > min_size
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
            for new_node in new_nodes {
                let new_node_idx = bhst.add_node(new_node.clone());

                // FOR TESTING
                // let ht = vcf.find_haplotype_for_sample(
                //     vcf.get_contig(),
                //     new_node.start_idx..new_node.stop_idx + 1,
                //     new_node.indexes[0],
                // );
                // let matching_idxs = crate::subcommands::check_for_haplotype::identical_haplotype_count(&vcf, &ht);
                // let len = new_node.indexes.len();

                // if matching_idxs.len() != len && new_node.start_idx != new_node.stop_idx {
                //     for h in &ht {
                //         println!("{:?} {}", h, vcf.get_nearest_idx_by_pos(h.pos));
                //     }
                //     println!(
                //         "{new_node_idx:?} {} {len} {}",
                //         matching_idxs.len(),
                //         ht.len()
                //     );
                //     panic!();
                // }

                bhst.add_edge(node_idx, new_node_idx, 0);
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
    let (mut left, mut right) = (node.start_idx, node.stop_idx);

    // Minus 1 to account for the starting variant itself as well
    if node_idx == NodeIndex::new(0) {
        right = right.saturating_sub(1);
    }

    let (mut rz, mut ro);
    // Expand to the right until a contradictory genotype is found or the end of sequencing data
    loop {
        (rz, ro) = (vec![], vec![]);
        if right < vcf.matrix.ncols() - 1 {
            right += 1;
        }

        for i in node.indexes.iter() {
            let right_gt = vcf.matrix.slice(s![*i, right]);

            match right_gt.into_scalar() {
                0 => rz.push(*i),
                1 => ro.push(*i),
                _ => unreachable!("Other genotypes than 0 and 1 are present in the matrix"),
            }
        }

        if (!rz.is_empty() && !ro.is_empty()) || right == vcf.matrix.ncols() - 1 {
            break;
        }
    }

    let (mut lz, mut lo);
    // Expand to the left until a contradictory genotype is found or the end of sequencing data
    loop {
        (lz, lo) = (vec![], vec![]);
        left = left.saturating_sub(1);

        for i in node.indexes.iter() {
            let left_gt = vcf.matrix.slice(s![*i, left]);

            match left_gt.into_scalar() {
                0 => lz.push(*i),
                1 => lo.push(*i),
                _ => unreachable!("Other genotypes than 0 and 1 are present in the matrix"),
            }
        }

        if (!lz.is_empty() && !lo.is_empty()) || left == 0 {
            break;
        }
    }

    // Divide left and right side contradictory genotypes to buckets:
    // [0,1] [0,0] [1, 0] [0, 0]
    let (mut zo, mut zz, mut oz, mut oo) = (vec![], vec![], vec![], vec![]);

    for i in &lz {
        if rz.contains(i) {
            zz.push(*i);
        } else if ro.contains(i) {
            zo.push(*i)
        }
    }

    for i in &lo {
        if rz.contains(i) {
            oz.push(*i);
        } else if ro.contains(i) {
            oo.push(*i)
        }
    }

    // Create a new node for each bucket if it is not empty
    let mut nodes = Vec::with_capacity(4);
    let block_len = vcf.get_pos(right) as f32 - vcf.get_pos(left) as f32;

    if !zz.is_empty() {
        let node = Node {
            block_len,
            start_idx: left,
            stop_idx: right,
            indexes: zz,
        };
        nodes.push(node);
    }

    if !zo.is_empty() {
        let node = Node {
            block_len,
            start_idx: left,
            stop_idx: right,
            indexes: zo,
        };
        nodes.push(node);
    }

    if !oz.is_empty() {
        let node = Node {
            block_len,
            start_idx: left,
            stop_idx: right,
            indexes: oz,
        };
        nodes.push(node);
    }

    if !oo.is_empty() {
        let node = Node {
            block_len,
            start_idx: left,
            stop_idx: right,
            indexes: oo,
        };
        nodes.push(node);
    }

    // Return buckets if more than one bucket exists
    match nodes.len() > 1 {
        true => Some(nodes),
        false => None,
    }
}
