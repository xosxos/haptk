use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;
use std::sync::Arc;

use ndarray::s;

use haptk::structs::PhasedMatrix;
use haptk::subcommands::hst_gwas::Assoc;

//// Alignment graph

pub struct AlignmentPoint {
    pub x: usize,
    pub y: usize,
    pub opt_value: f64,
    pub gt: u8,
}

#[derive(Debug, Clone)]
struct Level<'a> {
    level: usize,
    start_idx: usize,
    stop_idx: usize,
    ht: &'a [u8],
    opt_value: f64,
}

impl<'a> Level<'a> {
    fn new(level: usize, start_idx: usize, stop_idx: usize, ht: &'a [u8], opt_value: f64) -> Self {
        Self {
            level,
            start_idx,
            stop_idx,
            ht,
            opt_value,
        }
    }
}

impl<'a> Ord for Level<'a> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.level.cmp(&other.level)
    }
}

impl<'a> PartialOrd for Level<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> PartialEq for Level<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.level == other.level
    }
}
impl<'a> Eq for Level<'a> {}

pub fn return_alignment(
    assoc_rows: Rc<BTreeMap<u64, Box<dyn Assoc + Send>>>,
    vcf: Arc<PhasedMatrix>,
    p_limit: usize,
) -> (Vec<AlignmentPoint>, Vec<[f64; 2]>) {
    let mut levels = BTreeSet::<Level>::new();

    let mut weights: Vec<[f64; 2]> = vec![];
    let mut points = vec![];
    for pos in 0..vcf.ncoords() {
        let _ = levels.extract_if(|level| pos > level.stop_idx);

        for (_, node) in assoc_rows
            .iter()
            .filter(|(_, n)| pos == n.start_idx() && -f64::log10(n.opt_value()) > p_limit as f64)
        {
            let ht = vcf
                .matrix
                .slice(s![
                    node.first_sample_idx(),
                    node.start_idx()..node.stop_idx() + 1
                ])
                .to_slice()
                .unwrap();

            let nlevel = find_new_level(&levels);
            let lvl = Level::new(nlevel, pos, node.stop_idx(), ht, node.opt_value());
            levels.insert(lvl);
        }

        if !levels.is_empty() {
            weights.push([vcf.get_pos(pos) as f64, levels.len() as f64]);
        }
        for level in &levels {
            let idx = pos.saturating_sub(level.start_idx + 1);
            match level.ht[idx] {
                0 => points.push(AlignmentPoint {
                    x: pos,
                    y: level.level,
                    opt_value: -f64::log10(level.opt_value),
                    gt: 0,
                }),

                1 => points.push(AlignmentPoint {
                    x: pos,
                    y: level.level,
                    opt_value: -f64::log10(level.opt_value),
                    gt: 1,
                }),
                _ => panic!(),
            }
        }
    }
    (points, weights)
}

fn find_new_level(levels: &BTreeSet<Level>) -> usize {
    let mut prev_lev = 0;

    for (i, level) in levels.iter().enumerate() {
        if level.level != 0 && i == 0 {
            return 0;
        } else if level.level - prev_lev > 1 {
            return prev_lev + 1;
        } else {
            prev_lev = level.level;
        }
    }

    levels.len()
}
