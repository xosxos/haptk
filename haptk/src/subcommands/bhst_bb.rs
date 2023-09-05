use std::fmt;

use color_eyre::eyre::ensure;
use ndarray::s;
use petgraph::dot::{Config, Dot};
use petgraph::prelude::NodeIndex;
use petgraph::{Direction, Graph};
use rayon::prelude::*;
use rust_htslib::bcf::Read;
use serde::{Deserialize, Serialize};

use crate::args::{Selection, StandardArgs};
use crate::core::{get_htslib_contig_len, parse_snp_coord};
use crate::read_vcf::{get_reader, read_vcf_to_matrix};
use crate::structs::{Coord, HapVariant, PhasedMatrix};

#[derive(Serialize, Deserialize, Default, Clone, Debug, PartialEq, PartialOrd)]
pub struct Node {
    pub indexes: Vec<usize>,
    pub haplotype: Vec<HapVariant>,
}

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "n: {}\nnmarkers: {}\nblock length (bp): {}\n",
            self.indexes.len(),
            self.haplotype.len(),
            (self.haplotype.last().unwrap().pos - self.haplotype.first().unwrap().pos) as f32
        );
        write!(f, "{line}")
    }
}

pub trait Indexes {
    fn indexes(&self) -> &Vec<usize>;
    fn block_len(&self) -> f32;
    fn show(&self) -> String;
}

impl Indexes for Node {
    fn indexes(&self) -> &Vec<usize> {
        &self.indexes
    }
    fn block_len(&self) -> f32 {
        (self.haplotype.last().unwrap().pos - self.haplotype.first().unwrap().pos) as f32
    }

    fn show(&self) -> String {
        self.to_string()
    }
}

type Result<T> = std::result::Result<T, HstError>;

#[derive(Debug)]
enum HstError {
    PrematureEndOfSequencingData,
    EndOfSequencingData,
    EndOfLeafNodes,
}

impl fmt::Display for HstError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::PrematureEndOfSequencingData => write!(f, "Premature end of sequencing data"),
            Self::EndOfSequencingData => write!(f, "End of sequencing data"),
            Self::EndOfLeafNodes => write!(f, "End of leaf nodes"),
        }
    }
}

#[doc(hidden)]
pub fn run(args: StandardArgs) -> color_eyre::Result<()> {
    ensure!(
        args.selection == Selection::All,
        "For the biobank version of hst, only selecting all is supported for now."
    );

    let (contig, pos) = parse_snp_coord(&args.coords)?;

    let (vcf_first_pos, vcf_last_pos) = get_vcf_first_and_last_pos(&args, contig)?;

    let stop = if vcf_last_pos > pos + 1_000_000 {
        pos + 1_000_000
    } else {
        vcf_last_pos + 1
    };

    let (mut last_start, mut last_stop) = (pos.saturating_sub(1_000_000), stop);

    let mut vcf = read_vcf_to_matrix(&args, contig, pos, Some((last_start, last_stop)), None)?;
    let min_size = 1;

    let mut bhst = Graph::<_, _>::new();

    bhst.add_node(Node {
        indexes: (0..vcf.matrix.nrows()).collect(),
        haplotype: vec![],
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
            .map(|node_idx| {
                find_contradictory_gt(&vcf, &bhst, node_idx, vcf_first_pos, vcf_last_pos)
                    .map(|nodes| (node_idx, nodes))
            })
            .collect::<Result<Vec<(NodeIndex, Vec<Node>)>>>();

        let nodes: Vec<(NodeIndex, Vec<Node>)> = match nodes {
            Ok(nodes) => nodes,
            Err(HstError::PrematureEndOfSequencingData) => {
                (vcf, last_start, last_stop) =
                    match read_more_vcf(&args, &vcf, last_start, last_stop, vcf_last_pos) {
                        Err(HstError::EndOfSequencingData) => break,
                        Ok((vcf, last_start, last_stop)) => (vcf, last_start, last_stop),
                        _ => unreachable!(),
                    };
                continue;
            }
            Err(HstError::EndOfLeafNodes) => break,
            Err(HstError::EndOfSequencingData) => unreachable!(),
        };

        // Add nodes to the tree
        for (node_idx, new_nodes) in nodes {
            for new_node in new_nodes {
                let new_node_idx = bhst.add_node(new_node.clone());

                bhst.add_edge(node_idx, new_node_idx, 0);
            }
        }
    }

    // tracing::debug!("Finished HST construction: {bhst:?}");

    Ok(())
}

fn find_coord(vcf: &PhasedMatrix, idx: usize) -> Coord {
    vcf.coords()[idx].clone()
}

fn coord_to_hapvariant(coord: Coord, gt: u8) -> HapVariant {
    HapVariant {
        gt,
        contig: coord.contig,
        pos: coord.pos,
        reference: coord.reference,
        alt: coord.alt,
    }
}

fn find_previous_idx_from_hapvariant(
    vcf: &PhasedMatrix,
    hap: &HapVariant,
    vcf_first_pos: u64,
) -> Result<usize> {
    match vcf.coords().iter().rev().position(|c| c.pos < hap.pos) {
        Some(idx) => Ok(vcf.ncoords() - idx - 1),
        None => {
            if vcf.coords().first().unwrap().pos == vcf_first_pos {
                return Err(HstError::EndOfSequencingData);
            } else {
                return Err(HstError::PrematureEndOfSequencingData);
            }
        }
    }
}

fn find_next_idx_from_hapvariant(
    vcf: &PhasedMatrix,
    hap: &HapVariant,
    vcf_last_pos: u64,
) -> Result<usize> {
    match vcf.coords().iter().position(|c| c.pos > hap.pos) {
        Some(idx) => Ok(idx),
        None => {
            if vcf.coords().last().unwrap().pos == vcf_last_pos {
                return Err(HstError::EndOfSequencingData);
            } else {
                return Err(HstError::PrematureEndOfSequencingData);
            }
        }
    }
}

fn find_haplotype(
    vcf: &PhasedMatrix,
    haplotype: &Vec<HapVariant>,
    left_coord: &Coord,
    right_coord: &Coord,
    left_gt: u8,
    right_gt: u8,
) -> Vec<HapVariant> {
    let right_idx = vcf.idx_by_coord(right_coord).unwrap();
    let left_idx = vcf.idx_by_coord(left_coord).unwrap();

    let left_coord = vcf.coords()[left_idx].clone();
    let right_coord = vcf.coords()[right_idx].clone();
    let left_hap = coord_to_hapvariant(left_coord, left_gt);
    let right_hap = coord_to_hapvariant(right_coord, right_gt);
    let haplotype = vec![vec![left_hap], haplotype.clone(), vec![right_hap]];
    haplotype.into_iter().flatten().collect()
}

pub fn get_vcf_first_and_last_pos(
    args: &StandardArgs,
    contig: &str,
) -> color_eyre::Result<(u64, u64)> {
    let contig_len = get_htslib_contig_len(&args.file, contig)?;
    let mut reader = get_reader(&args.file, contig, None)?;

    let mut start_pos = 0;
    for record in reader.records() {
        let record = record?;
        // RUST-HTSLIB is 0-based so add 1
        start_pos = (record.pos() + 1) as u64;
        break;
    }

    let mut end_pos = 0;
    let mut x = 0;
    loop {
        let mut reader = get_reader(
            &args.file,
            contig,
            Some((contig_len - 10_000 - x, contig_len - x)),
        )?;

        let mut found_record = false;
        for record in reader.records() {
            let record = record?;
            // RUST-HTSLIB is 0-based so add 1
            end_pos = (record.pos() + 1) as u64;
            found_record = true;
        }
        if found_record {
            break;
        } else {
            x += 10_000;
        }
    }
    Ok((start_pos, end_pos))
}

fn read_more_vcf(
    args: &StandardArgs,
    vcf: &PhasedMatrix,
    old_start: u64,
    old_stop: u64,
    contig_len: u64,
) -> Result<(PhasedMatrix, u64, u64)> {
    let mut left_start = old_start.saturating_sub(1_000_000);
    let mut right_start = old_stop + 1_000_000;

    let mut right = None;
    let mut have_run_right_last = false;
    loop {
        if right_start > contig_len && !have_run_right_last {
            tracing::debug!("setting run right last to true {right_start} > {contig_len}");
            have_run_right_last = true;
        } else if have_run_right_last {
            break;
        }
        let right_vcf = read_vcf_to_matrix(
            &args,
            vcf.get_contig(),
            0,
            Some((old_stop, right_start)),
            None,
        )
        .expect("vcf read error in bhst bb right function read_more_vcf");

        if !right_vcf.coords().is_empty() {
            right = Some(right_vcf);
            break;
        }
        right_start = right_start + 1_000_000;
    }

    let mut left = None;
    let mut have_run_left_last = false;
    loop {
        if left_start == 0 && !have_run_left_last {
            tracing::debug!("setting run left last to true {left_start}");
            have_run_left_last = true;
        } else if have_run_left_last {
            break;
        }
        let left_vcf = read_vcf_to_matrix(
            &args,
            vcf.get_contig(),
            0,
            Some((left_start, old_start)),
            None,
        )
        .expect("vcf read error in bhst bb left function read_more_vcf");
        if !left_vcf.coords().is_empty() {
            left = Some(left_vcf);
            break;
        }
        left_start = left_start.saturating_sub(1_000_000);
    }

    let (matrix, coords, samples) = match (left, right) {
        (Some(l), Some(r)) => {
            let coords = l
                .coords()
                .iter()
                .chain(r.coords().iter())
                .cloned()
                .collect();
            let samples = l.samples().clone();
            let matrix = ndarray::concatenate!(ndarray::Axis(1), l.matrix, r.matrix);
            (matrix, coords, samples)
        }
        (Some(l), None) => {
            let coords = l.coords().clone();
            let samples = l.samples().clone();
            (l.matrix, coords, samples)
        }
        (None, Some(r)) => {
            let coords = r.coords().clone();
            let samples = r.samples().clone();
            (r.matrix, coords, samples)
        }
        (None, None) => return Err(HstError::EndOfSequencingData),
    };

    Ok((
        PhasedMatrix::new(
            //Variant idx
            0,
            matrix,
            samples,
            coords,
            &args.selection,
        ),
        left_start,
        right_start,
    ))
}

fn find_contradictory_gt(
    vcf: &PhasedMatrix,
    bhst: &Graph<Node, u8>,
    node_idx: NodeIndex,
    vcf_first_pos: u64,
    vcf_last_pos: u64,
) -> Result<Vec<Node>> {
    let node = bhst.node_weight(node_idx).unwrap();

    let (mut left, mut right) = if node_idx == NodeIndex::new(0) {
        let left = find_coord(vcf, vcf.variant_idx());
        // Minus 1 to account for the starting variant itself as well
        let right = find_coord(vcf, vcf.variant_idx() - 1);
        (coord_to_hapvariant(left, 1), coord_to_hapvariant(right, 1))
    } else {
        (
            node.haplotype.first().unwrap().clone(),
            node.haplotype.last().unwrap().clone(),
        )
    };
    tracing::debug!("Starting node: {node_idx:?} from left: {left} and right: {right}");

    let (mut l_haplotype, mut r_haplotype) = (vec![], vec![]);
    let (mut left_coord, mut right_coord) = (left.clone().into(), right.clone().into());

    let (mut rz, mut ro);

    // Expand to the right until a contradictory genotype is found or the end of sequencing data
    loop {
        (rz, ro) = (vec![], vec![]);
        let idx = find_next_idx_from_hapvariant(vcf, &right, vcf_last_pos);

        let idx = match idx {
            Ok(idx) => idx,
            Err(e) => match e {
                HstError::EndOfSequencingData => break,
                HstError::PrematureEndOfSequencingData => {
                    return Err(HstError::PrematureEndOfSequencingData)
                }
                HstError::EndOfLeafNodes => unreachable!(),
            },
        };

        for i in node.indexes.iter() {
            let right_gt = vcf.matrix.slice(s![*i, idx]);

            match right_gt.into_scalar() {
                0 => rz.push(*i),
                1 => ro.push(*i),
                _ => unreachable!("Other genotypes than 0 and 1 are present in the matrix"),
            }
        }

        let coord = find_coord(vcf, idx);
        right_coord = coord.clone();

        if !rz.is_empty() && !ro.is_empty() {
            tracing::debug!("{node_idx:?}: Found a contradictory gt on the right");
            break;
        }

        let hapvariant = if rz.is_empty() {
            coord_to_hapvariant(coord, 1)
        } else {
            coord_to_hapvariant(coord, 0)
        };
        r_haplotype.push(hapvariant.clone());
        right = hapvariant;
    }

    let (mut lz, mut lo);
    // Expand to the left until a contradictory genotype is found or the end of sequencing data
    loop {
        (lz, lo) = (vec![], vec![]);
        let temp_idx = vcf.idx_by_coord(&left.clone().into());
        // tracing::debug!("{node_idx:?}: old idx: {temp_idx:?}");

        let idx = find_previous_idx_from_hapvariant(vcf, &left, vcf_first_pos);
        // tracing::debug!("{node_idx:?}: new idx: {idx:?}");

        let idx = match idx {
            Ok(idx) => idx,
            Err(e) => match e {
                HstError::EndOfSequencingData => break,
                HstError::PrematureEndOfSequencingData => {
                    return Err(HstError::PrematureEndOfSequencingData)
                }
                HstError::EndOfLeafNodes => unreachable!(),
            },
        };

        for i in node.indexes.iter() {
            let left_gt = vcf.matrix.slice(s![*i, idx]);

            match left_gt.into_scalar() {
                0 => lz.push(*i),
                1 => lo.push(*i),
                _ => unreachable!("Other genotypes than 0 and 1 are present in the matrix"),
            }
        }

        let coord = find_coord(vcf, idx);
        left_coord = coord.clone();

        if !lz.is_empty() && !lo.is_empty() {
            tracing::debug!("{node_idx:?}: Found a contradictory gt on the left");
            break;
        }

        let hapvariant = if lz.is_empty() {
            coord_to_hapvariant(coord, 1)
        } else {
            coord_to_hapvariant(coord, 0)
        };

        l_haplotype.push(hapvariant.clone());
        left = hapvariant;
    }

    l_haplotype.extend(r_haplotype);

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

    if !zz.is_empty() {
        let node = Node {
            haplotype: find_haplotype(vcf, &l_haplotype, &left_coord, &right_coord, 0, 0),
            indexes: zz,
        };
        nodes.push(node);
    }

    if !zo.is_empty() {
        let node = Node {
            haplotype: find_haplotype(vcf, &l_haplotype, &left_coord, &right_coord, 0, 1),
            indexes: zo,
        };
        nodes.push(node);
    }

    if !oz.is_empty() {
        let node = Node {
            haplotype: find_haplotype(vcf, &l_haplotype, &left_coord, &right_coord, 1, 0),
            indexes: oz,
        };
        nodes.push(node);
    }

    if !oo.is_empty() {
        let node = Node {
            haplotype: find_haplotype(vcf, &l_haplotype, &left_coord, &right_coord, 1, 1),
            indexes: oo,
        };
        nodes.push(node);
    }

    tracing::debug!("Returning nodes");
    // Return buckets if more than one bucket exists
    match nodes.len() > 1 {
        true => Ok(nodes),
        false => Err(HstError::EndOfLeafNodes),
    }
}
