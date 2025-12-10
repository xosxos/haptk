use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::DefaultHasher;
use std::hash::Hasher;
use std::path::Path;
use std::path::PathBuf;

use petgraph::graph::EdgeIndex;
use petgraph::graph::NodeIndices;
use petgraph::graph::NodeReferences;
use petgraph::graph::NodeWeightsMut;
use petgraph::prelude::NodeIndex;
use petgraph::visit::IntoNodeReferences;
use petgraph::Direction;
use petgraph::Graph;
use serde::{Deserialize, Serialize};

use crate::args::ConciseArgs;
use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::Coord;
use crate::core::HapVariant;
use crate::core::PhasedMatrix;
use crate::core::Ploidy;
use crate::core::SelectedHaplotypes;
use crate::error::Error;
use crate::io::get_output;
use crate::utils::centromeres_hg38;

use super::initiate_hst;

pub trait HstMetadata {
    fn selection(&self) -> Selection;
    fn input_file(&self) -> PathBuf;
}

impl HstMetadata for &ConciseArgs {
    fn selection(&self) -> Selection {
        self.selection.clone()
    }

    fn input_file(&self) -> PathBuf {
        self.file.clone()
    }
}

impl HstMetadata for &StandardArgs {
    fn selection(&self) -> Selection {
        self.selection.clone()
    }

    fn input_file(&self) -> PathBuf {
        self.file.clone()
    }
}

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

    pub fn sample_names(
        &self,
        samples: &[String],
        ploidy: usize,
        selection: &SelectedHaplotypes,
    ) -> String {
        self.indexes
            .iter()
            .map(|i| format!("{}_{}", &*samples[i / ploidy], selection.find_ht_num(*i)))
            .collect::<Vec<String>>()
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hst {
    pub hst: Graph<Node, ()>,
    // #[serde(default)]
    pub metadata: Metadata,
}

impl Hst {
    pub fn new<S: HstMetadata>(vcf: &PhasedMatrix, args: S, hst_type: HstType) -> Self {
        let start_coord = vcf.start_coord();

        Hst {
            hst: initiate_hst(vcf, start_coord, None),
            metadata: Metadata::new(
                vcf,
                vcf.samples().clone(),
                hst_type,
                args.selection(),
                args.input_file(),
            ),
        }
    }

    pub fn new_with_indexes<S: HstMetadata>(
        vcf: &PhasedMatrix,
        args: &S,
        hst_type: HstType,
        root_indexes: Vec<usize>,
    ) -> Self {
        let start_coord = vcf.start_coord();

        Hst {
            hst: initiate_hst(vcf, start_coord, Some(root_indexes)),
            metadata: Metadata::new(
                vcf,
                vcf.samples().clone(),
                hst_type,
                args.selection(),
                args.input_file(),
            ),
        }
    }
}

// Inner methods
impl Hst {
    pub fn add_node(&mut self, node: Node) -> NodeIndex {
        self.hst.add_node(node)
    }

    pub fn add_edge(&mut self, parent: NodeIndex, child: NodeIndex) -> EdgeIndex {
        self.hst.add_edge(parent, child, ())
    }

    pub fn node_indices(&self) -> NodeIndices<u32> {
        self.hst.node_indices()
    }

    pub fn node_weights_mut(&mut self) -> NodeWeightsMut<'_, Node> {
        self.hst.node_weights_mut()
    }

    pub fn nodes(&self) -> NodeReferences<'_, Node> {
        self.hst.node_references()
    }

    pub fn node_weight(&self, node_idx: NodeIndex) -> Option<&Node> {
        self.hst.node_weight(node_idx)
    }

    pub fn n_children(&self, node_idx: NodeIndex) -> usize {
        self.hst
            .neighbors_directed(node_idx, Direction::Outgoing)
            .count()
    }

    pub fn get_children(&self, node_idx: NodeIndex) -> petgraph::graph::Neighbors<'_, ()> {
        self.hst.neighbors_directed(node_idx, Direction::Outgoing)
    }

    pub fn get_parent(&self, node_idx: NodeIndex) -> Option<NodeIndex> {
        self.hst
            .neighbors_directed(node_idx, Direction::Incoming)
            .next()
    }

    pub fn get_siblings(&self, node_idx: NodeIndex) -> Vec<NodeIndex> {
        self.hst
            // Get parent
            .neighbors_directed(node_idx, Direction::Incoming)
            // Get the children of the parent, i.e. siblings
            .flat_map(|parent| self.get_children(parent))
            // Filter self out
            .filter(|idx| *idx != node_idx)
            .collect()
    }
}

// Hst helpers
impl Hst {
    pub fn coords(&self) -> &BTreeSet<Coord> {
        &self.metadata.coords
    }

    pub fn get_sample_name(&self, index: usize) -> String {
        self.metadata
            .samples
            .get(index / *self.metadata.ploidy)
            .unwrap()
            .clone()
    }

    pub fn get_ht_num(&self, index: usize) -> usize {
        self.metadata.selected_haplotypes.find_ht_num(index)
    }

    pub fn nhaplotypes(&self) -> usize {
        let idx = NodeIndex::new(0);
        let root = self.node_weight(idx).unwrap();
        root.indexes.len()
    }

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

// Io
impl Hst {
    pub fn write_to_file(&mut self, path: PathBuf, publish: bool) -> Result<(), Error> {
        if publish {
            self.metadata.samples = vec!["r".to_string()];

            for node in self.hst.node_weights_mut() {
                node.indexes = vec![0; node.indexes.len()];
            }
        }

        tracing::info!("HST output: {path:?}.");
        let now = std::time::Instant::now();
        let mut output = get_output(Some(path))?;

        let mut writer = bgzip::BGZFWriter::new(&mut output, bgzip::Compression::default());

        serde_json::to_writer(&mut writer, &self)?;

        writer.close()?;

        tracing::info!("Wrote HSTs in {:?}", now.elapsed());
        Ok(())
    }

    pub fn from_file(path: &Path) -> Result<Hst, Error> {
        let file = std::fs::File::open(path).map_err(|e| Error::Io {
            e,
            path: path.to_path_buf(),
        })?;

        let reader = bgzip::BGZFReader::new(file)?;
        let hst: Hst = serde_json::from_reader(reader)?;

        Ok(hst)
    }
}

impl Hst {
    pub fn majority_branch(&self) -> Vec<NodeIndex> {
        let mut idx = NodeIndex::new(0);
        let mut majority_nodes = vec![idx];

        while let Some((_data, node_idx)) = self
            .get_children(idx)
            .map(|idx| (self.node_weight(idx).unwrap(), idx))
            //
            // Max by returns None if the iterator is empty, thus when a leaf node is reached
            .max_by(|x, y| x.0.indexes.len().cmp(&y.0.indexes.len()))
        {
            majority_nodes.push(node_idx);
            idx = node_idx;
        }

        majority_nodes
    }

    pub fn subtree(&mut self, root: NodeIndex) {
        let mut nodes = vec![NodeIndex::from(0)];
        let mut edges = vec![];
        let mut count = 0;
        let mut mapping = HashMap::new();
        mapping.insert(NodeIndex::from(0), root);

        self.recurse_down(
            root,
            NodeIndex::from(0),
            &mut nodes,
            &mut edges,
            &mut mapping,
            &mut count,
        );

        // println!("{} {}", nodes.len(), edges.len());
        let mut graph = Graph::<Node, ()>::from_edges(edges);
        // println!("{:#?}", graph);

        for node in nodes {
            // println!("{node:?}");
            let weight = graph.node_weight_mut(node).unwrap();
            *weight = self
                .node_weight(mapping.get(&node).unwrap().clone())
                .unwrap()
                .clone();
        }

        self.hst = graph;
        // self.hst.retain_nodes(|_g, idx| subtree.contains(&idx));
    }

    fn recurse_down<'a>(
        &'a self,
        parent: NodeIndex,
        parent_new: NodeIndex,
        nodes: &mut Vec<NodeIndex>,
        edges: &'a mut Vec<(NodeIndex, NodeIndex)>,
        mapping: &mut HashMap<NodeIndex, NodeIndex>,
        count: &mut u32,
    ) {
        for orig_child_idx in self.get_children(parent) {
            *count += 1;
            let child = NodeIndex::from(*count);
            nodes.push(child);
            edges.push((parent_new, child));
            mapping.insert(child, orig_child_idx);
            self.recurse_down(orig_child_idx, child, nodes, edges, mapping, count)
        }
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
    pub selected_haplotypes: SelectedHaplotypes,
}

impl Metadata {
    pub fn new(
        vcf: &PhasedMatrix,
        samples: Vec<String>,
        hst_type: HstType,
        selection: Selection,
        input_vcf_name: PathBuf,
    ) -> Self {
        Self {
            start_coord: vcf.start_coord().clone(),
            coords: vcf.coords().clone(),
            contig: vcf.get_contig().clone(),
            samples,
            hst_type,
            selection,
            ploidy: vcf.ploidy.clone(),
            vcf_name: input_vcf_name,
            selected_haplotypes: vcf.metadata.get_selected_haplotypes().clone(),
        }
    }
}
