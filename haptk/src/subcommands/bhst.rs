use std::{
    collections::{BTreeSet, HashSet},
    hash::{DefaultHasher, Hasher},
    path::PathBuf,
    sync::mpsc::sync_channel,
};

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
    error::HaptkError,
    io::{get_output, open_csv_writer, push_to_output, write_haplotype},
    libs::{
        read_vcf::{get_sample_names, read_vcf_to_matrix},
        structs::{CoordDataSlot, PhasedMatrix},
    },
    structs::{Coord, HapVariant, Ploidy},
    utils::{centromeres_hg38, parse_snp_coord},
};

#[derive(Serialize, Deserialize, Default, Clone, Debug, PartialEq, PartialOrd, Hash, Eq)]
pub struct Node {
    pub indexes: Vec<usize>,
    pub start: Coord,
    pub stop: Coord,
    pub haplotype: Vec<u8>,
}

impl Node {
    pub fn identifier(&self) -> String {
        let ht = self.haplotype.iter().fold(
            format!(
                "{}_{}_{}_",
                self.start.contig, self.start.pos, self.stop.pos
            ),
            |acc, e| format!("{acc}{e}"),
        );

        let mut hasher = DefaultHasher::new();
        hasher.write(ht.as_bytes());

        format!("{:02x}", hasher.finish())
    }

    pub fn zygosity(&self, samples: &[String], ploidy: usize) -> (usize, usize) {
        let prior = self.indexes.len();
        let post: HashSet<_> = self.indexes.iter().map(|i| &samples[i / ploidy]).collect();

        let nhomozygotes = prior - post.len();
        let nheterozygotes = (2 * post.len()) - prior;

        (nheterozygotes, nhomozygotes)
    }

    pub fn sample_name_list(&self, samples: &[String], ploidy: usize) -> String {
        self.indexes
            .iter()
            .map(|i| &*samples[i / ploidy])
            .collect::<Vec<&str>>()
            .join(";")
    }

    pub fn check_for_centromere_hg38(&self) -> bool {
        let start = self.start.pos;
        let stop = self.stop.pos;

        let (cen_start, cen_stop) = centromeres_hg38(&self.start.contig);

        // Start pos is inside centromere
        let c1 = start > cen_start && start < cen_stop;

        // Stop pos is inside centromere
        let c2 = stop > cen_start && stop < cen_stop;

        // Starts before centromere and ends after centromere
        let c3 = start < cen_start && stop > cen_stop;

        c1 || c2 || c3
    }
}

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "n: {}\nstart: {}\nstop: {}\nblock length (bp): {}\n",
            self.indexes.len(),
            self.start.pos,
            self.stop.pos,
            self.stop.pos.saturating_sub(self.start.pos),
        );
        write!(f, "{line}")
    }
}

pub fn read_vcf_with_selections(args: &StandardArgs, window: Option<u64>) -> Result<PhasedMatrix> {
    let (contig, pos) = parse_snp_coord(&args.coords)?;

    let (indexes, _) = get_sample_names(args, contig, None)?;
    ensure!(
        indexes.len() > 1,
        "Cannot build a tree with less than 2 samples."
    );
    let mut vcf = read_vcf_to_matrix(args, contig, pos, None, None, window, false)?;

    if args.selection == Selection::OnlyLongest {
        vcf.select_only_longest()?
    }

    Ok(vcf)
}

#[doc(hidden)]
pub fn run(args: StandardArgs, min_size: usize, publish: bool, window: Option<u64>) -> Result<()> {
    ensure!(
        args.selection != Selection::Unphased,
        "Running with unphased data is not supported."
    );

    let mut vcf = read_vcf_with_selections(&args, window)?;

    ensure!(
        vcf.nhaplotypes() >= min_size,
        "VCF has less haplotypes than the required minimum node size ({} < {min_size})",
        vcf.nhaplotypes()
    );

    let start_coord = vcf.start_coord().clone();

    // Construct the bidirectional HST
    let bhst = construct_bhst(&mut vcf, &start_coord, min_size)?;
    tracing::info!("Finished HST construction.");

    // Find the majority based haplotype
    let mbah = find_mbah(&bhst, &vcf)?;

    // Write it to file
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "bhst_ancestral_haplotype", "csv");
    write_haplotype(mbah, open_csv_writer(sh_output)?)?;

    // Find the shared core haplotype
    let ht = find_shared_haplotype(&bhst, &vcf);

    // Write it to file
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "bhst_core_haplotype", "csv");
    write_haplotype(ht, open_csv_writer(sh_output)?)?;

    // Write the HST to file
    let mut hst_output = args.output.clone();
    push_to_output(&args, &mut hst_output, "bhst", "hst.gz");
    write_hst_file(bhst, &vcf, hst_output, publish, args, HstType::Bhst)?;

    Ok(())
}

#[doc(hidden)]
pub fn construct_bhst(
    vcf: &mut PhasedMatrix,
    start_coord: &Coord,
    min_size: usize,
) -> Result<Graph<Node, ()>> {
    // Initate the HST by inserting the root node, and two child nodes if there is a contradictory genotype right at the starting variant
    let mut hst = initiate_hst(vcf, start_coord, None);

    let mut blacklist_nodes: Vec<NodeIndex> = vec![];

    loop {
        // Find contradictory genotypes for samples in leaf nodes
        // Insert new nodes into a the HST using a parallel iterator and  channels
        // Return a blacklist of nodes to not exclude from iteration each round
        let blacklist = insert_nodes_to_bhst(vcf, &mut hst, &blacklist_nodes, min_size);

        match blacklist {
            Err(e) => {
                // If data runs, and the file has not ended, read more from disk
                vcf.read_more(e)?;
                continue;
            }
            Ok(blacklist) => {
                if blacklist.is_empty() {
                    break;
                }
                blacklist_nodes.extend(blacklist.iter().flatten());
            }
        }
    }

    Ok(hst)
}

pub fn initiate_hst(
    vcf: &PhasedMatrix,
    start_coord: &Coord,
    indexes: Option<Vec<usize>>,
) -> Graph<Node, ()> {
    let mut hst = Graph::new();

    let indexes = indexes.unwrap_or_else(|| (0..vcf.nhaplotypes()).collect());

    let root = Node {
        // start: Cow::Borrowed(start_coord),
        // stop: Cow::Borrowed(start_coord),
        start: start_coord.clone(),
        stop: start_coord.clone(),
        indexes: indexes.clone(),
        haplotype: vec![],
    };

    let _ = hst.add_node(root.clone());

    // If the first variant is already contradictory, split into two nodes
    if vcf.is_contradictory(start_coord, &indexes) {
        let genotypes = vcf.matrix_column(start_coord);
        let (mut ones, mut zeroes) = (vec![], vec![]);

        for i in indexes.iter() {
            match genotypes[*i] == 1 {
                true => ones.push(*i),
                false => zeroes.push(*i),
            }
        }

        let node1 = Node {
            start: start_coord.clone(),
            stop: start_coord.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(start_coord..=start_coord, ones[0]),
            indexes: ones,
        };

        let node2 = Node {
            start: start_coord.clone(),
            stop: start_coord.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(start_coord..=start_coord, zeroes[0]),
            indexes: zeroes,
        };

        let node_idx1 = hst.add_node(node1.clone());
        let node_idx2 = hst.add_node(node2.clone());

        hst.add_edge(NodeIndex::new(0), node_idx1, ());
        hst.add_edge(NodeIndex::new(0), node_idx2, ());
    };

    hst
}

pub fn insert_nodes_to_bhst(
    vcf: &PhasedMatrix,
    hst: &mut Graph<Node, ()>,
    blacklist_nodes: &[NodeIndex],
    min_size: usize,
) -> std::result::Result<Vec<Vec<NodeIndex>>, HaptkError> {
    let (tx, rx) = sync_channel(2024);

    let blacklist = hst
        .node_indices()
        .par_bridge()
        // Filter out nodes that were blacklisted (genotyping data ran out or node haplotypes min_size was reached )
        .filter(|n| !blacklist_nodes.contains(n))
        // Filter in leaf nodes with 2 or more haplotypes
        .filter_map(|node_idx| {
            let node = hst.node_weight(node_idx).unwrap();

            let count = hst
                .neighbors_directed(node_idx, Direction::Outgoing)
                .count();

            match count == 0 && node.indexes.len() >= min_size.max(2) {
                true => Some((node_idx, node.clone())),
                false => None,
            }
        })
        // Search for contradictory genotypes by extending the node haplotype
        .map(
            |(parent_idx, node)| match find_contradictory_gt_bhst(vcf, &node) {
                Err(e) => Err(e),
                Ok(Some(nodes)) => {
                    let black_list_nodes: Vec<NodeIndex> = nodes
                        .into_iter()
                        .filter_map(|insert_node| {
                            if insert_node.indexes.len() >= min_size {
                                tx.send((parent_idx, insert_node)).unwrap();
                                None
                            } else {
                                Some(parent_idx)
                            }
                        })
                        .collect();
                    Ok(black_list_nodes)
                }
                _ => Ok(vec![parent_idx]),
            },
        )
        .collect::<std::result::Result<Vec<Vec<NodeIndex>>, HaptkError>>();

    drop(tx);

    // Receive nodes for the iterator and add them to the tree
    // We use a channel because if genotyping data needes to be dynamically extended, we dont want to lose the already computed nodes
    while let Ok((parent_idx, insert_node)) = rx.recv() {
        let idx = hst.add_node(insert_node.clone());
        hst.add_edge(parent_idx, idx, ());
    }

    blacklist
}

#[doc(hidden)]
pub fn find_contradictory_gt_bhst(
    vcf: &PhasedMatrix,
    node: &Node,
) -> std::result::Result<Option<Vec<Node>>, HaptkError> {
    let (left_coord, right_coord) = (&node.start, &node.stop);

    // NOTE: Remains to be seen what the overhead here would be
    // let (prev, next) = thread::scope(|s| {
    //     let prev_handle = s.spawn(|| vcf.prev_contradictory(left_coord, &node.indexes));
    //     let next = vcf.next_contradictory(right_coord, &node.indexes);
    //     let prev = prev_handle.join().unwrap();
    //     (prev, next)
    // });

    let prev = vcf.prev_contradictory(left_coord, &node.indexes);
    let next = vcf.next_contradictory(right_coord, &node.indexes);

    let (prev, next) = match (prev, next) {
        (Err(_), Err(_)) => return Err(HaptkError::HstBothEndError),
        (_, Err(_)) => return Err(HaptkError::HstRightEndError),
        (Err(_), _) => return Err(HaptkError::HstLeftEndError),
        (Ok(prev), Ok(next)) => (prev, next),
    };

    // NOTE: Benchmarking required
    // Allocate Vecs with capacity so no reallocation is required
    let (mut zo, mut zz, mut oz, mut oo) = (
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
    );

    // Fill the genotype buckets and then create HST nodes from the buckets
    match (prev, next) {
        (Some((left, left_vec)), Some((right, right_vec))) => {
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
            Ok(Some(nodes))
        }
        // Handle if data ran out on the right side
        (Some((left, left_vec)), None) => {
            if !vcf.is_genome_wide() {
                tracing::warn!(
                    "Genotyping data ran out on the right with samples {:?}",
                    vcf.get_sample_names(&node.indexes)
                );
            }

            for i in node.indexes.iter() {
                match left_vec[*i] == 1 {
                    true => oz.push(*i),
                    false => zz.push(*i),
                }
            }
            let nodes =
                create_nodes_from_buckets(vcf, left, vcf.coords().last().unwrap(), oo, oz, zo, zz);
            Ok(Some(nodes))
        }
        // Handle if data ran out on the left side
        (None, Some((right, right_vec))) => {
            if !vcf.is_genome_wide() {
                tracing::warn!(
                    "Genotyping data ran out on the left with samples {:?}",
                    vcf.get_sample_names(&node.indexes)
                );
            }

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
            Ok(Some(nodes))
        }
        // Handle if data ran out on both sides
        (None, None) => Ok(None),
    }
}

fn create_nodes_from_buckets(
    vcf: &PhasedMatrix,
    left: &Coord,
    right: &Coord,
    oo: Vec<usize>,
    oz: Vec<usize>,
    zo: Vec<usize>,
    zz: Vec<usize>,
) -> Vec<Node> {
    let mut nodes = Vec::with_capacity(4);
    if !zz.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, zz[0]),
            indexes: zz,
        };
        nodes.push(node);
    }

    if !zo.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, zo[0]),
            indexes: zo,
        };
        nodes.push(node);
    }

    if !oz.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, oz[0]),
            indexes: oz,
        };
        nodes.push(node);
    }

    if !oo.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, oo[0]),
            indexes: oo,
        };
        nodes.push(node);
    }
    nodes
}

// HST UTILS

pub fn find_mbah(g: &Graph<Node, ()>, vcf: &PhasedMatrix) -> Result<Vec<HapVariant>> {
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

    // Start and stop idxs are the first deviating genotype and it cannot be a part of the haplotype
    let start = vcf.coords().range(&last_node.start..).nth(1).unwrap();

    let start = match start > &last_node.stop {
        true => &last_node.start,
        false => start,
    };

    Ok(vcf.find_haplotype_for_sample(start..&last_node.stop, last_node.indexes[0]))
}

pub fn find_shared_haplotype(g: &Graph<Node, ()>, vcf: &PhasedMatrix) -> Vec<HapVariant> {
    let start_idx = NodeIndex::new(0);
    let nodes = find_majority_nodes(g, start_idx);
    let (second_node, _second_node_idx) = nodes[1];

    let start = vcf.coords().range(&second_node.start..).nth(1).unwrap();

    let start = match start > &second_node.stop {
        true => &second_node.start,
        false => start,
    };

    vcf.find_haplotype_for_sample(start..&second_node.stop, second_node.indexes[0])
}

pub fn calculate_block_len(node: &Node) -> u64 {
    match node.stop.cmp(&node.start) {
        std::cmp::Ordering::Greater => node.stop.pos - node.start.pos,
        std::cmp::Ordering::Less => node.start.pos - node.stop.pos,
        std::cmp::Ordering::Equal => 0,
    }
}

pub fn find_majority_nodes(g: &Graph<Node, ()>, start_idx: NodeIndex) -> Vec<(&Node, NodeIndex)> {
    let mut idx = start_idx;
    let start_node = g.node_weight(idx).unwrap();
    let mut majority_nodes = vec![(start_node, idx)];

    // Max by returns None if the iterator is empty, thus when a leaf node is reached
    while let Some(max) = g
        .neighbors_directed(idx, Direction::Outgoing)
        .map(|idx| (g.node_weight(idx).unwrap(), idx))
        .max_by(|x, y| x.0.indexes.len().cmp(&y.0.indexes.len()))
    {
        majority_nodes.push(max);
        idx = max.1;
    }

    majority_nodes
}

pub fn find_lowest_start(g: &Graph<Node, ()>) -> Coord {
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

pub fn find_highest_stop(g: &Graph<Node, ()>) -> Coord {
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

pub fn find_coord_list(g: &Graph<Node, ()>, vcf: &PhasedMatrix) -> Vec<Coord> {
    let start_idx = find_lowest_start(g);
    let stop_idx = find_highest_stop(g);

    vcf.coords().range(start_idx..=stop_idx).cloned().collect()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hst {
    pub hst: Graph<Node, ()>,
    // #[serde(default)]
    pub metadata: Metadata,
}

impl Hst {
    pub fn get_haplotype(&self, node: &Node) -> Vec<HapVariant> {
        if node.haplotype.is_empty() {
            return vec![];
        }

        self.metadata
            .coords
            .range(&node.start..=&node.stop)
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
}

#[derive(Default, Debug, Serialize, Deserialize, Clone)]
pub enum HstType {
    HstLeft,
    HstRight,
    #[default]
    Bhst,
}

#[derive(Default, Debug, Serialize, Deserialize, Clone)]
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
    mut hst: Graph<Node, ()>,
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
        hst,
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

pub fn read_hst_file(path: PathBuf) -> Result<Hst> {
    let file = std::fs::File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;
    let reader = bgzip::BGZFReader::new(file)?;
    let hst: Hst = serde_json::from_reader(reader)?;

    Ok(hst)
}
