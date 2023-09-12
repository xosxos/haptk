use std::collections::{BTreeSet, HashMap};

use eframe::egui::plot::PlotPoint;
use egui::{Color32, RichText};
use ndarray::s;
use petgraph::{graph::NodeIndex, Direction, Graph};

use fishers_exact::fishers_exact;
use haptk::structs::PhasedMatrix;
use haptk::subcommands::bhst::Node;
use haptk::subcommands::hst_scan::TreeRow;

//// Tree graph
#[derive(Default, Clone, Debug)]
pub struct TreePoint {
    pub text: String,
    pub x: f64,
    pub y: f64,
    pub p: f64,
    pub idx: NodeIndex,
    pub ctrls_true: usize,
}

impl TreePoint {
    fn new(x: f64, y: f64, p: f64, ctrls_true: usize, idx: NodeIndex, text: String) -> Self {
        Self {
            x,
            y,
            p,
            ctrls_true,
            idx,
            text,
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
        self.x == other.x
    }
}
impl Eq for TreePoint {}

impl From<TreePoint> for egui::plot::Text {
    fn from(value: TreePoint) -> Self {
        let logp = -(f64::log10(value.p)) as i64;
        let v = match logp * 20 > 255 {
            true => 255,
            false => (logp * 20) as u8,
        };
        let color = Color32::from_rgb(v, 0, 0);

        let text = RichText::new(value.text)
            .size(13.0)
            .strong()
            .color(color)
            .background_color(Color32::from_rgb(255, 255, 255));

        if logp >= 8 {
            Self::new(PlotPoint::new(value.x, value.y), text)
                .color(color)
                .highlight(true)
        } else {
            Self::new(PlotPoint::new(value.x, value.y), text).color(color)
        }
    }
}

pub fn return_tree_text_and_lines(
    vcf: &PhasedMatrix,
    ctrl_vcf: &PhasedMatrix,
    trees: &[TreeRow],
    idx: usize,
) -> (BTreeSet<TreePoint>, Vec<Vec<[f64; 2]>>) {
    let TreeRow { tree, .. } = &trees[idx];
    let pvalues = get_pvalues(tree, vcf, ctrl_vcf);

    let width = 1000.0;
    let root_idx = NodeIndex::new(0);
    let root_node = tree.node_weight(root_idx).unwrap();

    let (p, ctrls_true) = pvalues.get(&root_idx).unwrap();
    let root = TreePoint::new(
        width / 2.0,
        0.0,
        *p,
        *ctrls_true,
        root_idx,
        root_node.indexes.len().to_string(),
    );

    let mut plot_text: BTreeSet<TreePoint> = BTreeSet::new();
    let mut plot_lines: Vec<_> = vec![];

    plot_text.insert(root.clone());

    recursive_draw(
        tree,
        root_idx,
        root,
        (0.0, width),
        0.0,
        &pvalues,
        &mut plot_lines,
        &mut plot_text,
    );
    (plot_text, plot_lines)
}

fn get_pvalues(
    tree: &Graph<Node, u8>,
    vcf: &PhasedMatrix,
    ctrl_vcf: &PhasedMatrix,
) -> HashMap<NodeIndex, (f64, usize)> {
    let mut hm = HashMap::new();
    for node_idx in tree.node_indices() {
        let node = tree.node_weight(node_idx).unwrap();
        let index = node.indexes[0];
        let ht = vcf.matrix.slice(s![index, node.start_idx..node.stop_idx]);

        let mut ctrls_true = 0;
        for j in 0..ctrl_vcf.matrix.nrows() {
            let ctrl_slice = ctrl_vcf.matrix.slice(s![j, node.start_idx..node.stop_idx]);
            if ht == ctrl_slice {
                ctrls_true += 1;
            }
        }

        let cases_true = node.indexes.len();
        let cases_false = vcf.matrix.nrows() - cases_true;
        let controls_false = ctrl_vcf.matrix.nrows() - ctrls_true;

        let p = fishers_exact(&[
            cases_true as u32,
            ctrls_true as u32,
            cases_false as u32,
            controls_false as u32,
        ])
        .unwrap();
        hm.insert(node_idx, (p.greater_pvalue, ctrls_true));
    }
    hm
}

fn recursive_draw(
    tree: &Graph<Node, u8>,
    parent_idx: NodeIndex,
    parent_point: TreePoint,
    width: (f64, f64),
    parent_y: f64,
    pvalues: &HashMap<NodeIndex, (f64, usize)>,
    plot_lines: &mut Vec<Vec<[f64; 2]>>,
    plot_text: &mut BTreeSet<TreePoint>,
) {
    let (start, stop) = width;

    let parent_node = tree.node_weight(parent_idx).unwrap();
    let children_idxs: Vec<NodeIndex> = tree
        .neighbors_directed(parent_idx, Direction::Outgoing)
        .collect();

    if children_idxs.len() < 2 {
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

        let additional_y = f64::max((node.block_len - parent_node.block_len) as f64, 1.0 * 2.0);
        let additional_y = f64::min(additional_y, 1.0 * 14.0);

        let y = parent_y - additional_y;

        let (p, ctrls_true) = pvalues.get(&node_idx).unwrap();
        let point = TreePoint {
            text: node.indexes.len().to_string(),
            x: center,
            y,
            p: *p,
            ctrls_true: *ctrls_true,
            idx: node_idx,
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
            tree,
            node_idx,
            point,
            (startn, stopn),
            y,
            pvalues,
            plot_lines,
            plot_text,
        );
    }
}

pub fn create_path(parent: &TreePoint, child: &TreePoint, is_parent: bool) -> Vec<[f64; 2]> {
    let font_size = 13.0;

    match is_parent {
        true => vec![
            [parent.x, parent.y - font_size * 0.18],
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
