use std::collections::BTreeSet;
use std::rc::Rc;

use eframe::egui;
use egui::Color32;
use egui_plot::PlotPoints;
use egui_plot::{Line, PlotPoint, PlotUi};

use haptk::structs::PhasedMatrix;

use crate::app::{Hst, State};
use crate::plot::init_plot;
use crate::tree::{return_tree_text_and_lines, TreePoint};
use crate::utils::{write_indexes_to_file, write_tp_ht_to_csv};

pub fn central_panel(
    ctx: &egui::Context,
    state: &mut State,
    vcf: Rc<PhasedMatrix>,
    bhst: Rc<Hst>,
    uhst_left: Rc<Hst>,
    uhst_right: Rc<Hst>,
    nmin_samples: usize,
    decoy_samples: Rc<Vec<String>>,
) -> egui::InnerResponse<()> {
    egui::CentralPanel::default().show(ctx, |ui| {
        let plot = init_plot(
            state.clone(),
            vcf.clone(),
            bhst.clone(),
            uhst_left.clone(),
            uhst_right.clone(),
            nmin_samples,
            decoy_samples.clone(),
        );

        plot.show(ui, |plot_ui| match state {
            State::Bhst => draw_hst(plot_ui, bhst, "bhst", nmin_samples, vcf, decoy_samples),
            State::UhstLeft => draw_hst(
                plot_ui,
                uhst_left,
                "uhst_left",
                nmin_samples,
                vcf,
                decoy_samples,
            ),
            State::UhstRight => draw_hst(
                plot_ui,
                uhst_right,
                "uhst_right",
                nmin_samples,
                vcf,
                decoy_samples,
            ),
        });
    })
}

fn draw_hst(
    plot_ui: &mut PlotUi,
    hst: Rc<Hst>,
    name: &str,
    nmin_samples: usize,
    vcf: Rc<PhasedMatrix>,
    decoy_samples: Rc<Vec<String>>,
) {
    let (text, points) = return_tree_text_and_lines(
        hst.clone(),
        nmin_samples,
        vcf.clone(),
        decoy_samples.clone(),
    );

    for p in &points {
        let p: PlotPoints = p.iter().cloned().collect();
        let line = Line::new(p).color(Color32::from_rgb(0, 0, 0));
        plot_ui.line(line);
    }
    for t in &text {
        plot_ui.text(t.clone().into());
    }

    if plot_ui.response().clicked() {
        if let Some(pp) = plot_ui.pointer_coordinate() {
            if let Some(p) = find_closest(&text, pp) {
                if let Some(node) = text.get(p).map(|p| hst.node_weight(p.idx).unwrap()) {
                    write_indexes_to_file(name, node, vcf.clone());
                    write_tp_ht_to_csv(name, node, vcf.clone());
                };
            }
        }
    }
}

fn find_closest(texts: &BTreeSet<TreePoint>, point: PlotPoint) -> Option<&TreePoint> {
    let interact_radius_sq = (16.0f64).powi(2);

    texts
        .iter()
        .map(|item| {
            let (a, b) = (item.x, item.y);
            let (x, y) = (point.x, point.y);
            ((a - x).abs().powi(2) + (b - y).abs().powi(2), item)
        })
        .min_by(|(a, _), (b, _)| a.total_cmp(b))
        .filter(|(dist, _)| *dist <= interact_radius_sq)
        .map(|(_, elem)| elem)
}
