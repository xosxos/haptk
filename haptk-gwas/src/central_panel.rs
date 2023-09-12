use std::collections::BTreeMap;
use std::rc::Rc;
use std::sync::Arc;

use eframe::egui;
use eframe::egui::plot::{PlotPoints, Points};
use egui::plot::{Line, PlotBounds};
use egui::Color32;
use itertools::Itertools;

use haptk::structs::PhasedMatrix;
use haptk::subcommands::annotate_haplotype::GtfRow;
use haptk::subcommands::bhst::Node;
use haptk::subcommands::hst_gwas::return_binary_assoc;
use haptk::subcommands::hst_mrca_gwas::calculate_mrca;
use haptk::subcommands::hst_scan::TreeRow;

use crate::alignment::return_alignment;
use crate::app::{
    AlignmentCache, AlignmentPlot, AssocCache, Key, MbahPlot, State, TopNodeMrcaCache, TreeCache,
};
use crate::plot::init_plot;
use crate::tree::return_tree_text_and_lines;
use crate::utils::assoc_vec_to_bm;

pub fn central_panel(
    ctx: &egui::Context,
    min_samples: usize,
    max_samples: usize,
    min_ht_len: usize,
    max_ht_len: usize,
    pvalue: usize,
    pvalue_limit: f64,
    state: &mut State,
    alignment_bounds: &mut PlotBounds,
    alignment_bounds_set: &mut bool,
    local_bounds_set: &mut bool,
    vcf: Arc<PhasedMatrix>,
    ctrl_vcf: &PhasedMatrix,
    trees: &Rc<Vec<TreeRow>>,
    tree_mbah_lengths: &[(usize, Node)],
    assoc_cache: &mut AssocCache,
    alignment_cache: &mut AlignmentCache,
    tree_cache: &mut TreeCache,
    top_node_mrca_cache: &mut TopNodeMrcaCache,
    rec_rates: &Option<BTreeMap<u64, f32>>,
    gff3: &Option<Vec<GtfRow>>,
) -> egui::InnerResponse<()> {
    egui::CentralPanel::default().show(ctx, |ui| {
        let key = Key {
            min_samples,
            max_samples,
            min_ht_len,
            max_ht_len,
        };

        let assoc = match assoc_cache.get(&key) {
            Some(assoc) => assoc,
            None => {
                let assoc = return_binary_assoc(
                    trees.as_ref(),
                    vcf.clone(),
                    (min_samples, max_samples, min_ht_len, min_ht_len),
                    ctrl_vcf,
                    None,
                );
                let assoc = assoc_vec_to_bm(assoc);

                assoc_cache.insert(key.clone(), Rc::new(assoc));
                assoc_cache.get(&key).unwrap()
            }
        };

        let (points, weights) = match alignment_cache.get(&(key.clone(), pvalue)) {
            Some(v) => v,
            None => {
                let alignment = return_alignment(assoc.clone(), vcf.clone(), pvalue);
                alignment_cache.insert((key.clone(), pvalue), alignment);
                alignment_cache.get(&(key.clone(), pvalue)).unwrap()
            }
        };

        let mrcas = if let Some(rec_rates) = rec_rates {
            match top_node_mrca_cache.get(&key) {
                Some(mrcas) => Some(mrcas),
                None => {
                    let mut mrcas = vec![];
                    for (pos, assoc) in assoc.iter() {
                        let tree = trees.get(vcf.get_nearest_idx_by_pos(*pos)).unwrap();
                        let mrca =
                            calculate_mrca(assoc.node_idx(), &tree.tree, rec_rates, vcf.clone());
                        mrcas.push((*pos, mrca));
                    }
                    top_node_mrca_cache.insert(key.clone(), mrcas);
                    top_node_mrca_cache.get(&key)
                }
            }
        } else {
            None
        };

        if !*alignment_bounds_set {
            let min = [
                points.iter().min_by_key(|p| p.x).unwrap().x as f64 - 50.0,
                0.0,
            ];
            let max = [
                points.iter().max_by_key(|p| p.x).unwrap().x as f64 + 50.0,
                points.iter().max_by_key(|p| p.y).unwrap().y as f64 + 50.0,
            ];
            let plot_bounds = PlotBounds::from_min_max(min, max);
            *alignment_bounds = plot_bounds;
            *alignment_bounds_set = true;
        }

        let plot = init_plot(
            state.clone(),
            assoc.clone(),
            vcf.clone(),
            tree_cache.clone(),
            trees.clone(),
        );

        plot.show(ui, |plot_ui| match state {
            State::AlignmentWeight => {
                if let Some(gtf_rows) = gff3 {
                    let min = weights.iter().max_by_key(|w| w[0] as u64).unwrap();
                    let max = weights.iter().max_by_key(|w| w[0] as u64).unwrap();
                    gtf_rows
                        .iter()
                        .filter(|r| r.stop > min[0] as u64)
                        .filter(|r| r.start < max[0] as u64)
                        .map(|r| vec![[r.start as f64, -10.0], [r.stop as f64, -10.0]])
                        .map(|l| Line::new(l).width(10.0).color(Color32::from_rgb(0, 0, 0)))
                        .for_each(|line| {
                            plot_ui.line(line);
                        });

                    //     gtf_rows
                    //         .iter()
                    //         .filter(|r| r.stop > alignment_bounds.min()[0] as u64)
                    //         .filter(|r| r.start < alignment_bounds.max()[0] as u64)
                    //         .map(|r| {
                    //             let text = RichText::new(r.attribute.clone());
                    //             plot::Text::new(
                    //                 [(r.stop - (r.stop - r.start) / 2) as f64, 0.0].into(),
                    //                 text,
                    //             )
                    //         })
                    //         .for_each(|text| {
                    //             plot_ui.text(text);
                    //         });
                }
                plot_ui.line(Line::new(weights.iter().cloned().collect::<PlotPoints>()))
            }
            State::MbahLengths(state) => {
                let points: PlotPoints = match state {
                    MbahPlot::NumberOfMarkers => tree_mbah_lengths
                        .iter()
                        .map(|(idx, node)| {
                            [
                                vcf.get_pos(*idx) as f64,
                                (node.stop_idx - node.start_idx) as f64,
                            ]
                        })
                        .collect(),

                    MbahPlot::BasePairs => tree_mbah_lengths
                        .iter()
                        .map(|(idx, node)| [vcf.get_pos(*idx) as f64, node.block_len as f64])
                        .collect(),
                };
                plot_ui.line(egui::plot::Line::new(points));
            }
            State::TopNodeMrca => {
                if let Some(mrcas) = mrcas {
                    // Last 10 elements moving average
                    let mut moving_average = vec![];
                    for (a, b, c, d, e, f, g, h, i, j) in mrcas.iter().tuples() {
                        let sum = a.1 + b.1 + c.1 + d.1 + e.1 + f.1 + g.1 + h.1 + i.1 + j.1;
                        let mean = sum / 10.0;
                        moving_average.push((j.0, mean));
                    }

                    plot_ui.line(Line::new(
                        moving_average
                            .iter()
                            .map(|(idx, node)| [*idx as f64, f64::log10(*node)])
                            .collect::<PlotPoints>(),
                    ));
                }
            }
            State::NumberOfMarkers => {
                plot_ui.line(Line::new(
                    assoc
                        .iter()
                        .map(|(key, row)| [*key as f64, row.marker_len() as f64])
                        .collect::<PlotPoints>(),
                ));
            }
            State::BlockLength => {
                plot_ui.line(Line::new(
                    assoc
                        .iter()
                        .map(|(key, row)| [*key as f64, row.bp_len() as f64])
                        .collect::<PlotPoints>(),
                ));
            }
            State::Association => {
                if plot_ui.plot_clicked() {
                    if let Some(point) = plot_ui.pointer_coordinate() {
                        let idx = vcf.get_nearest_idx_by_pos(point.x as u64);
                        *state = State::HaplotypeStateTree(idx);
                    }
                }

                plot_ui.points(
                    Points::new(
                        assoc
                            .iter()
                            .map(|(key, row)| [*key as f64, -(f64::log10(row.opt_value()))])
                            .collect::<PlotPoints>(),
                    )
                    .radius(3.0),
                );
            }

            State::HaplotypeStateTree(idx) => {
                let exists = tree_cache.borrow().get(idx).is_some();

                match exists {
                    true => {
                        let cache = tree_cache.borrow();
                        let (text, points) = cache.get(idx).unwrap();
                        for p in points {
                            let p: PlotPoints = p.iter().cloned().collect();
                            let line = Line::new(p).color(Color32::from_rgb(0, 0, 0));
                            plot_ui.line(line);
                        }
                        for t in text.iter() {
                            plot_ui.text(t.clone().into());
                        }
                    }
                    false => {
                        let (text, points) =
                            return_tree_text_and_lines(&vcf, ctrl_vcf, trees.as_ref(), *idx);

                        for p in &points {
                            let p: PlotPoints = p.iter().cloned().collect();
                            let line = Line::new(p).color(Color32::from_rgb(0, 0, 0));
                            plot_ui.line(line);
                        }
                        for t in &text {
                            plot_ui.text(t.clone().into());
                        }
                        tree_cache.borrow_mut().insert(*idx, (text, points));
                    }
                }
            }

            State::Alignment(state) => {
                let (rpoints, gpoints) = match state {
                    AlignmentPlot::Pvalue => {
                        if *alignment_bounds_set && !*local_bounds_set {
                            plot_ui.set_plot_bounds(*alignment_bounds);
                            *local_bounds_set = true;
                        }
                        *alignment_bounds = plot_ui.plot_bounds();

                        let (mut high_points, mut low_points) = (vec![], vec![]);

                        points
                            .iter()
                            .filter(|p| p.x % 2 == 0)
                            .filter(|p| p.x > alignment_bounds.min()[0] as usize)
                            .filter(|p| p.x < alignment_bounds.max()[0] as usize)
                            .map(|p| ([p.x as f64, p.y as f64], p.opt_value))
                            .for_each(|(p, pv)| match pv >= pvalue_limit {
                                true => high_points.push(p),
                                false => low_points.push(p),
                            });

                        (
                            Points::new(high_points)
                                .radius(1.0)
                                .color(egui::Color32::from_black_alpha(250))
                                .shape(egui::widgets::plot::MarkerShape::Square),
                            Points::new(low_points)
                                .radius(1.0)
                                .color(egui::Color32::from_black_alpha(50))
                                .shape(egui::widgets::plot::MarkerShape::Square),
                        )
                    }
                    AlignmentPlot::RefAlt => {
                        if *alignment_bounds_set && !*local_bounds_set {
                            plot_ui.set_plot_bounds(*alignment_bounds);
                            *local_bounds_set = true;
                        }

                        *alignment_bounds = plot_ui.plot_bounds();

                        let (mut ref_points, mut alt_points) = (vec![], vec![]);

                        points
                            .iter()
                            .filter(|p| p.x % 2 == 0)
                            .filter(|p| p.x > alignment_bounds.min()[0] as usize)
                            .filter(|p| p.x < alignment_bounds.max()[0] as usize)
                            .map(|p| ([p.x as f64, p.y as f64], p.gt))
                            .for_each(|(p, gt)| match gt == 0 {
                                true => ref_points.push(p),
                                false => alt_points.push(p),
                            });

                        (
                            Points::new(ref_points)
                                .radius(1.0)
                                .color(egui::Color32::from_rgb(0, 255, 0))
                                .shape(egui::widgets::plot::MarkerShape::Square),
                            Points::new(alt_points)
                                .radius(1.0)
                                .color(egui::Color32::from_rgb(255, 0, 0))
                                .shape(egui::widgets::plot::MarkerShape::Square),
                        )
                    }
                };

                plot_ui.points(rpoints);
                plot_ui.points(gpoints);
            }
        });
    })
}
