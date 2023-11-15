use std::rc::Rc;

use egui_plot::{Plot, PlotPoint};
use haptk::structs::PhasedMatrix;

use crate::app::{Hst, State};
use crate::tree::{return_tree_text_and_lines, TreePoint};

// Init plot
pub fn init_plot(
    state: State,
    vcf: Rc<PhasedMatrix>,
    bhst: Rc<Hst>,
    uhst_left: Rc<Hst>,
    uhst_right: Rc<Hst>,
    nmin_samples: usize,
    decoy_samples: Rc<Vec<String>>,
) -> Plot {
    let decoy_samples = decoy_samples.clone();

    let plot = Plot::new("measurements")
        .auto_bounds_x()
        .auto_bounds_y()
        .label_formatter(move |_name, value| {
            let text = match state {
                State::Bhst => show_label(
                    bhst.clone(),
                    nmin_samples,
                    vcf.clone(),
                    decoy_samples.clone(),
                    value,
                ),
                State::UhstLeft => show_label(
                    uhst_left.clone(),
                    nmin_samples,
                    vcf.clone(),
                    decoy_samples.clone(),
                    value,
                ),
                State::UhstRight => show_label(
                    uhst_right.clone(),
                    nmin_samples,
                    vcf.clone(),
                    decoy_samples.clone(),
                    value,
                ),
            };

            text.unwrap_or(format!("x: {}, y: {}", value.x, value.y))
        });

    plot.show_axes([false, false])
}

fn show_label(
    hst: Rc<Hst>,
    nmin_samples: usize,
    vcf: Rc<PhasedMatrix>,
    decoy_samples: Rc<Vec<String>>,
    value: &PlotPoint,
) -> Option<String> {
    let (points, _) = return_tree_text_and_lines(hst.clone(), nmin_samples, vcf, decoy_samples);
    points
        .get(&TreePoint {
            x: value.x,
            ..Default::default()
        })
        .map(|p| hst.node_weight(p.idx).unwrap())
        .map(|node| format!("{node}"))
}
