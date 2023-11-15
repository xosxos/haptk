use std::collections::BTreeSet;
use std::rc::Rc;

use egui::{Color32, RichText};
use egui_plot::PlotPoint;
use petgraph::{graph::NodeIndex, Direction, Graph};

use haptk::{structs::PhasedMatrix, subcommands::bhst::Node};

#[derive(Default, Clone, Debug)]
pub struct TreePoint {
    pub text: String,
    pub x: f64,
    pub y: f64,
    pub idx: NodeIndex,
    pub is_decoy_leaf_node: bool,
}

impl TreePoint {
    fn new(x: f64, y: f64, idx: NodeIndex, text: String, is_decoy_leaf_node: bool) -> Self {
        Self {
            x,
            y,
            idx,
            text,
            is_decoy_leaf_node,
        }
    }
}

impl Ord for TreePoint {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.x.total_cmp(&other.x)
    }
}

impl PartialOrd for TreePoint {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TreePoint {
    fn eq(&self, other: &Self) -> bool {
        let other_x = f64::trunc(other.x * 1.0) / 1.0;
        let other_y = f64::trunc(other.y * 1.0) / 1.0;
        let self_x = f64::trunc(self.x * 1.0) / 1.0;
        let self_y = f64::trunc(self.y * 1.0) / 1.0;
        self_x == other_x && self_y == other_y
    }
}
impl Eq for TreePoint {}

impl From<TreePoint> for egui_plot::Text {
    fn from(value: TreePoint) -> Self {
        let color = match value.is_decoy_leaf_node {
            true => Color32::from_rgb(255, 0, 0),
            false => Color32::from_rgb(0, 0, 0),
        };

        let text = RichText::new(value.text)
            .size(13.0)
            .strong()
            .color(color)
            .background_color(Color32::from_rgb(255, 255, 255));

        Self::new(PlotPoint::new(value.x, value.y), text).color(color)
    }
}

pub fn return_tree_text_and_lines(
    tree: Rc<Graph<Node, u8>>,
    nmin_samples: usize,
    vcf: Rc<PhasedMatrix>,
    decoy_samples: Rc<Vec<String>>,
) -> (BTreeSet<TreePoint>, Vec<Vec<[f64; 2]>>) {
    let width = 1000.0;
    let root_idx = NodeIndex::new(0);
    let root_node = tree.node_weight(root_idx).unwrap();

    let root = TreePoint::new(
        width / 2.0,
        0.0,
        root_idx,
        root_node.indexes.len().to_string(),
        false,
    );

    let mut plot_text: BTreeSet<TreePoint> = BTreeSet::new();
    let mut plot_lines: Vec<_> = vec![];

    plot_text.insert(root.clone());

    recursive_draw(
        tree.clone(),
        vcf,
        decoy_samples,
        root_idx,
        root,
        (0.0, width),
        0.0,
        &mut plot_lines,
        &mut plot_text,
        nmin_samples,
    );
    (plot_text, plot_lines)
}

fn recursive_draw(
    tree: Rc<Graph<Node, u8>>,
    vcf: Rc<PhasedMatrix>,
    decoy_samples: Rc<Vec<String>>,
    parent_idx: NodeIndex,
    parent_point: TreePoint,
    width: (f64, f64),
    parent_y: f64,
    plot_lines: &mut Vec<Vec<[f64; 2]>>,
    plot_text: &mut BTreeSet<TreePoint>,
    nmin_samples: usize,
) {
    let (start, stop) = width;

    let children_idxs: Vec<NodeIndex> = tree
        .neighbors_directed(parent_idx, Direction::Outgoing)
        .collect();

    let parent_node = tree.node_weight(parent_idx).unwrap();

    if children_idxs.len() < 2 || parent_node.indexes.len() < nmin_samples {
        return;
    }

    let mut children: Vec<(&_, NodeIndex)> = children_idxs
        .iter()
        .map(|c| (tree.node_weight(*c).unwrap(), *c))
        .collect();

    children.sort_by(|a, b| b.0.indexes.len().cmp(&a.0.indexes.len()));

    let mut split_positions = vec![];
    let mut prev_split_pos = start;

    let sum_of_idxs: usize = children.iter().map(|c| c.0.indexes.len()).sum();

    for (i, node_idx) in children.into_iter() {
        let ratio = i.indexes.len() as f64 / sum_of_idxs as f64;
        let curr_split_pos = prev_split_pos + ((stop - start) * ratio);
        split_positions.push((node_idx, i, prev_split_pos, curr_split_pos));
        prev_split_pos = curr_split_pos;
    }

    let mut points = vec![];
    for (i, (node_idx, node, startn, stopn)) in split_positions.into_iter().enumerate() {
        let center = startn + (stopn - startn) / 2.0;

        let y = parent_y - 1.0;

        let decoys = vcf
            .get_sample_idxs(decoy_samples.as_ref())
            .unwrap_or(vec![]);
        let is_decoy_node = node.indexes.iter().any(|i| decoys.contains(i));

        let point = TreePoint {
            text: node.indexes.len().to_string(),
            x: center,
            y,
            idx: node_idx,
            is_decoy_leaf_node: is_decoy_node,
        };

        points.push(point.clone());

        let ps = if i >= 1 {
            create_path(&points[i - 1], &point, false)
        } else {
            create_path(&parent_point, &point, true)
        };
        plot_lines.push(ps);
        plot_text.insert(point.clone());

        recursive_draw(
            tree.clone(),
            vcf.clone(),
            decoy_samples.clone(),
            node_idx,
            point,
            (startn, stopn),
            y,
            plot_lines,
            plot_text,
            nmin_samples,
        );
    }
}

pub fn create_path(parent: &TreePoint, child: &TreePoint, is_parent: bool) -> Vec<[f64; 2]> {
    match is_parent {
        true => vec![
            [parent.x, parent.y],
            [parent.x, child.y],
            [child.x, child.y],
        ]
        .into_iter()
        .collect(),
        false => vec![[parent.x, parent.y], [child.x, child.y]]
            .into_iter()
            .collect(),
    }
}
