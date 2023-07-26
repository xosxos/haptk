use color_eyre::Result;
use petgraph::{graph::NodeIndex, Direction, Graph};
use svg::node::element::path::Data;
use svg::node::element::Element;
use svg::node::element::Path;
use svg::node::Text;
use svg::Document;
use svg::Node;

use crate::args::GraphArgs;
use crate::structs::PhasedMatrix;
use crate::subcommands::bhst;

#[derive(Clone, Debug)]
pub struct Point {
    text: String,
    y: f32,
    x: f32,
    color: String,
    opacity: f32,
}

impl std::fmt::Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!("Point {{ x: {}, y: {} }}", self.x, self.y);
        write!(f, "{line}")
    }
}

#[derive(Debug)]
pub struct HstGraph<'a, T: bhst::Indexes + Clone> {
    pub g: &'a Graph<T, u8>,
    pub vcf: &'a PhasedMatrix,
    pub s: GraphArgs,
    pub document: Document,
    pub tree_height: f32,
    // Basepairs per pixel
    pub scale: f32,
    pub nsamples: usize,
    pub top_padding: f32,
    pub side_padding: f32,
    pub largest: usize,
    pub min_size: usize,
    pub variables: Option<Vec<String>>,
    pub decoys: Option<Vec<String>>,
}

impl<'a, T: bhst::Indexes + Clone> HstGraph<'a, T> {
    pub fn new(
        g: &'a Graph<T, u8>,
        vcf: &'a PhasedMatrix,
        mut s: GraphArgs,
        decoys: Option<Vec<String>>,
        min_size: usize,
    ) -> Self {
        let mut largest = 0.0;
        for idx in g.node_indices() {
            let node = g.node_weight(idx).unwrap();
            if node.block_len() > largest {
                largest = node.block_len();
            }
        }
        let nsamples = g.node_weight(NodeIndex::new(0)).unwrap().indexes().len();

        let min_width = find_width(g, &s, min_size);
        if min_width > s.width {
            tracing::warn!("Width has been incremented to {min_width}");
            s.width = min_width;
        }

        let min_height = find_height(g, &s, min_size);
        if min_height > s.height {
            tracing::warn!("Height has been incremented to {min_height}");
            s.height = min_height;
        }

        let top_padding = s.height * 0.02;
        let btm_padding = s.height * 0.05;
        let side_padding = 0.02;

        let document = Document::new()
            .set("viewBox", (0, 0, s.width, s.height))
            .set("style", format!("background-color:{}", s.background_color));

        Self {
            document,
            tree_height: s.height - top_padding - btm_padding,
            g,
            vcf,
            // Basepairs per pixel
            scale: largest / (s.height - top_padding - btm_padding),
            s,
            nsamples,
            top_padding,
            side_padding,
            largest: largest as usize,
            variables: None,
            min_size,
            decoys,
        }
    }

    pub fn draw_graph(&mut self) -> Result<()> {
        let root_idx = NodeIndex::new(0);
        let width = (0.0, self.s.width as f32);

        // Rescale by subtracting the pixels lost to font scaling from the tree height
        // let extra_y = 0.0;
        // self.scale = (self.scale * self.tree_height) / (self.tree_height - extra_y);

        let root_node = self.g.node_weight(root_idx).unwrap();

        let root = Point {
            text: root_node.indexes().len().to_string(),
            x: self.s.width / 2.,
            y: self.top_padding,
            color: self.s.color.clone(),
            opacity: 1.0,
        };

        self.document.append(draw_point(&root, &self.s));

        self.recursive_draw(root_idx, root, width, self.top_padding, self.min_size)?;

        Ok(())
    }

    fn recursive_draw(
        &mut self,
        parent_idx: NodeIndex,
        mut parent_point: Point,
        width: (f32, f32),
        parent_y: f32,
        min_size: usize,
    ) -> Result<()> {
        let (start, stop) = width;

        let parent_node = self.g.node_weight(parent_idx).unwrap();
        let children_idxs: Vec<NodeIndex> = self
            .g
            .neighbors_directed(parent_idx, Direction::Outgoing)
            .collect();

        if children_idxs.len() < 2 || parent_node.indexes().len() < min_size {
            // self.draw_branch_data(parent_point.clone(), parent_idx);

            //Redraw parent point for possible decoys
            parent_point.color = return_color(parent_node, &self.decoys, &self.s, self.vcf, true)?;
            self.document.append(draw_point(&parent_point, &self.s));

            return Ok(());
        }

        let mut children: Vec<(&_, NodeIndex)> = children_idxs
            .iter()
            .map(|c| (self.g.node_weight(*c).unwrap(), *c))
            .collect();

        children.sort_by(|a, b| b.0.indexes().len().cmp(&a.0.indexes().len()));

        let mut split_positions = vec![];
        let mut prev_split_pos = start;
        let sum_of_idxs: usize = children.iter().map(|c| c.0.indexes().len()).sum();
        let children_len = children.len();

        for (i, node_idx) in children.into_iter() {
            let ratio = i.indexes().len() as f32 / sum_of_idxs as f32;
            let curr_split_pos = prev_split_pos + ((stop - start) * ratio);
            split_positions.push((node_idx, i, prev_split_pos, curr_split_pos));
            prev_split_pos = curr_split_pos;
        }

        let mut points = vec![];
        for (i, (node_idx, node, startn, stopn)) in split_positions.into_iter().enumerate() {
            let mut center = startn + (stopn - startn) / 2.0;
            if i == children_len - 2 && children_len == 3 && node.indexes().len() < 4 {
                center -= self.s.font_size * 0.5;
            }

            if i == children_len.saturating_sub(2) && children_len == 4 && node.indexes().len() < 4
            {
                center -= self.s.font_size * 0.5;
            }

            if i == children_len.saturating_sub(3) && children_len == 4 && node.indexes().len() < 6
            {
                center -= self.s.font_size * 0.5;
            }

            let additional_y = f32::max(
                (node.block_len() - parent_node.block_len()) / self.scale,
                self.s.font_size * 1.2,
            );
            let additional_y = f32::min(additional_y, self.s.font_size * 1.2);

            let y = parent_y + additional_y;

            let point = Point {
                text: node.indexes().len().to_string(),
                x: center,
                y,
                color: return_color(node, &self.decoys, &self.s, self.vcf, false)?,
                opacity: return_opacity(node, &self.decoys, &self.s, self.vcf)?,
            };
            points.push(point.clone());

            self.document.append(draw_point(&point, &self.s));
            if i > 1 {
                self.document
                    .append(create_path(&points[i - 1], &point, &self.s, 1.0, false));
            } else {
                self.document
                    .append(create_path(&parent_point, &point, &self.s, 1.0, true));
            }

            self.recursive_draw(node_idx, point, (startn, stopn), y, min_size)?;
        }

        Ok(())
    }
}

fn find_width<'a, T: bhst::Indexes>(g: &Graph<T, u8>, s: &GraphArgs, min_size: usize) -> f32 {
    let leaf_node_n = g
        .node_indices()
        .filter(|idx| {
            let parent_idx = g.neighbors_directed(*idx, Direction::Incoming).next();
            let node = g.node_weight(*idx).unwrap();
            if let Some(parent_idx) = parent_idx {
                let parent = g.node_weight(parent_idx).unwrap();
                parent.indexes().len() > min_size && node.indexes().len() <= min_size
            } else {
                false
            }
        })
        .count();
    leaf_node_n as f32 * s.font_size * 1.1
}

fn find_height<'a, T: bhst::Indexes>(g: &Graph<T, u8>, s: &GraphArgs, min_size: usize) -> f32 {
    if let Some(height) = g
        .node_indices()
        .filter(|idx| {
            let parent_idx = g.neighbors_directed(*idx, Direction::Incoming).next();
            let node = g.node_weight(*idx).unwrap();
            if let Some(parent_idx) = parent_idx {
                let parent = g.node_weight(parent_idx).unwrap();
                parent.indexes().len() > min_size && node.indexes().len() <= min_size
            } else {
                false
            }
        })
        .map(|idx| {
            petgraph::algo::all_simple_paths::<Vec<_>, _>(g, NodeIndex::new(0), idx, 0, None)
                .map(|path| path.iter().count())
        })
        .flatten()
        .max()
    {
        tracing::debug!("Height of the tree is {height}");
        height as f32 * s.font_size * 1.2
    } else {
        0.0
    }
}

fn return_color<'a, T: bhst::Indexes>(
    node: &'a T,
    decoys: &Option<Vec<String>>,
    s: &'a GraphArgs,
    vcf: &'a PhasedMatrix,
    last_node: bool,
) -> Result<String> {
    if let Some(decoys) = decoys {
        let decoys = vcf.get_sample_idxs(decoys)?;
        if !last_node {
            return match node.indexes().iter().any(|i| decoys.contains(i))
                && node.indexes().len() == 1
            {
                true => Ok("red".to_string()),
                false => Ok(s.color.to_string()),
            };
        } else {
            return match node.indexes().iter().any(|i| decoys.contains(i)) {
                true => Ok("red".to_string()),
                false => Ok(s.color.to_string()),
            };
        }
    }
    Ok(s.color.to_string())
}

fn return_opacity<'a, T: bhst::Indexes>(
    node: &'a T,
    decoys: &Option<Vec<String>>,
    _s: &'a GraphArgs,
    vcf: &PhasedMatrix,
) -> Result<f32> {
    if let Some(decoys) = decoys {
        let decoys = vcf.get_sample_idxs(decoys)?;
        return match node.indexes().iter().any(|i| decoys.contains(i)) && node.indexes().len() == 1
        {
            true => Ok(1.0),
            false => Ok(1.0),
        };
    }
    Ok(1.0)
}

pub fn draw_point(p: &Point, s: &GraphArgs) -> Element {
    let mut element = Element::new("text");
    element.assign("x", p.x);
    element.assign("y", p.y);
    element.assign("opacity", p.opacity);
    element.assign("font-weight", "bold");
    element.assign("fill", p.color.clone());
    element.assign("font-size", format!("{}px", s.font_size));
    element.append(Text::new(p.text.clone()));
    element
}

pub fn branch_data(p: &Point, s: &GraphArgs, o: f32, color: &str) -> Element {
    let mut element = Element::new("text");
    element.assign("x", p.x);
    element.assign("y", p.y);
    element.assign("opacity", o);
    element.assign("font-weight", "bold");
    element.assign("fill", color);
    element.assign("font-size", format!("{}px", s.font_size * 0.7));
    element.append(Text::new(p.text.clone()));
    element
}

pub fn create_path(parent: &Point, child: &Point, s: &GraphArgs, o: f32, is_parent: bool) -> Path {
    let font_size = s.font_size;
    let x_start = match parent.text.chars().count() {
        4 => font_size * 1.2,
        3 => font_size * 1.4,
        2 => font_size * 1.0,
        1 => font_size * 0.3,
        _ => font_size * 3.0,
    };

    let shorten_path = if parent.x > child.x {
        match child.text.chars().count() {
            4 => font_size * 2.5,
            3 => font_size * 1.9,
            2 => font_size * 1.25,
            1 => font_size * 0.7,
            _ => font_size * 3.0,
        }
    } else {
        match child.text.chars().count() {
            4 => font_size * -0.74,
            3 => font_size * -0.54,
            2 => font_size * -0.32,
            1 => font_size * -0.10,
            _ => font_size * -1.0,
        }
    };

    let data = match is_parent {
        true => Data::new()
            .move_to((parent.x + x_start, parent.y + font_size * 0.18))
            .line_to((parent.x + x_start, child.y - font_size * 0.2))
            .line_to((child.x + shorten_path, child.y - font_size * 0.2)),
        false => Data::new()
            .move_to((parent.x + x_start * 1.5, parent.y - font_size * 0.2))
            .line_to((child.x + shorten_path, child.y - font_size * 0.2)),
    };

    Path::new()
        .set("fill", "none")
        .set("stroke", s.color.as_str())
        .set("opacity", o)
        .set("stroke-width", s.stroke_width)
        .set("d", data)
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn point_display() {
        let point = Point { x: 5.0, y: 6.0, text: "".to_string(), color: "".to_string(), opacity: 1.0 };
        assert_eq!("Point { x: 5, y: 6 }".to_string(), format!("{point}"));
    }
}
