use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;
use std::sync::Arc;

use eframe::egui;
use eframe::egui::plot::Plot;

use haptk::structs::PhasedMatrix;
use haptk::subcommands::hst_gwas::Assoc;
use haptk::subcommands::hst_scan::TreeRow;

use crate::app::{AlignmentPlot, Ps, State};
use crate::tree::TreePoint;
use crate::utils::{write_indexes_to_file, write_tp_ht_to_file};

// Init plot
pub fn init_plot(
    state: State,
    assoc: Rc<BTreeMap<u64, Box<dyn Assoc + Send>>>,
    vcf: Arc<PhasedMatrix>,
    tree_cache: Rc<RefCell<BTreeMap<usize, (BTreeSet<TreePoint>, Vec<Ps>)>>>,
    trees: Rc<Vec<TreeRow>>,
) -> egui::plot::Plot {
    let statec = state.clone();
    let plot = Plot::new("measurements")
        .auto_bounds_x()
        .auto_bounds_y()
        .label_formatter(move |_name, value| {
            let print_row = |row: Option<&Box<dyn Assoc + Send>>| match row {
                Some(row) => format!("{}", row.show()),
                None => format!("x: {}, y: {}", value.x, value.y),
            };

            match statec {
                State::Alignment(AlignmentPlot::Pvalue | AlignmentPlot::RefAlt) => {
                    let x = core::cmp::min(value.x as usize, vcf.ncoords() - 1);
                    let x = vcf.get_pos(x);
                    print_row(assoc.get(&x))
                }

                State::HaplotypeStateTree(tree_idx) => {
                    if let Some((points, _)) = &tree_cache.borrow().get(&tree_idx) {
                        match points.get(&TreePoint {
                            x: value.x,
                            ..Default::default()
                        }) {
                            Some(tp) => {
                                let TreeRow { tree, .. } = &trees[tree_idx];
                                let node = tree.node_weight(tp.idx).unwrap();

                                write_tp_ht_to_file(node, vcf.clone(), tree_idx, tp.ctrls_true);
                                write_indexes_to_file(node, vcf.clone(), tree_idx, tp.ctrls_true);

                                format!(
                                    "n: {}\nstart: {}\nstop:{}\nnmarkers: {}\nlength:{}\nctrls_true: {}\np: {}\n",
                                    node.indexes.len(),
                                    vcf.get_pos(node.start_idx),
                                    vcf.get_pos(node.stop_idx),
                                    node.stop_idx - node.start_idx + 1,
                                    node.block_len,
                                    tp.ctrls_true,
                                    tp.p
                                )
                            }
                            None => format!("x: {}, y: {}", value.x, value.y),
                        }
                    } else {
                        format!("x: {}, y: {}", value.x, value.y)
                    }
                }
                _ => print_row(assoc.get(&(value.x as u64))),
            }
        });

    match state {
        State::HaplotypeStateTree(_) => plot.show_axes([false, false]),
        _ => plot,
    }
}
