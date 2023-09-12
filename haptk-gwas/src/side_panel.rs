use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

use eframe::egui;
use egui::Margin;
use egui_file::FileDialog;

use haptk::io::read_recombination_file;
use haptk::structs::PhasedMatrix;
use haptk::subcommands::annotate_haplotype::{read_gtf_file_to_vec, GtfRow};

use crate::app::{AlignmentPlot, FiletypeToOpen, MbahPlot, State};

pub fn side_panel(
    ctx: &egui::Context,
    vcf: Arc<PhasedMatrix>,
    min_samples: &mut usize,
    max_samples: &mut usize,
    min_ht_len: &mut usize,
    max_ht_len: &mut usize,
    pvalue_limit: &mut f64,
    pvalue: &mut usize,
    rec_rates: &mut Option<BTreeMap<u64, f32>>,
    gff3: &mut Option<Vec<GtfRow>>,
    state: &mut State,
    open_file_dialog: &mut Option<FileDialog>,
    opened_file: &mut Option<PathBuf>,
    filetype_to_open: &mut FiletypeToOpen,
) -> egui::InnerResponse<()> {
    egui::SidePanel::left("sidebar")
        .default_width(360.0)
        .frame(egui::Frame {
            fill: egui::Color32::WHITE,
            inner_margin: Margin::symmetric(5.0, 5.0),
            ..Default::default()
        })
        .show(ctx, |ui| {
            ui.heading("HST filtering parameters");

            let slider = egui::Slider::new(min_samples, 2..=200).text("min_samples");
            ui.add(slider);

            let slider = egui::Slider::new(max_samples, 2..=1000).text("max_samples");
            ui.add(slider);

            let slider = egui::Slider::new(min_ht_len, 1..=2000).text("min_ht_len");
            ui.add(slider);

            let slider = egui::Slider::new(max_ht_len, 2..=20000).text("max_ht_len");
            ui.add(slider);

            let slider = egui::Slider::new(pvalue, 1..=400).text("pvalue");
            ui.add(slider);

            ui.add_space(5.0);
            ui.heading("Views");
            if ui.button("MBAH lengths").clicked() {
                *state = State::MbahLengths(MbahPlot::NumberOfMarkers);
            }

            if state == &mut State::MbahLengths(MbahPlot::BasePairs) {
                if ui.button("MBAH: number of markers").clicked() {
                    *state = State::MbahLengths(MbahPlot::NumberOfMarkers);
                }
            } else if state == &mut State::MbahLengths(MbahPlot::NumberOfMarkers)
                && ui.button("MBAH: base pairs").clicked()
            {
                *state = State::MbahLengths(MbahPlot::BasePairs);
            }

            if ui.button("Number of markers").clicked() {
                *state = State::NumberOfMarkers;
            }

            if ui.button("Block length").clicked() {
                *state = State::BlockLength;
            }

            if ui.button("P-value").clicked() {
                *state = State::Association;
            }

            if rec_rates.is_some() && ui.button("Top node MRCA").clicked() {
                *state = State::TopNodeMrca;
            }

            if ui.button("Alignment").clicked() {
                *state = State::Alignment(AlignmentPlot::RefAlt);
            }

            if state == &mut State::Alignment(AlignmentPlot::Pvalue) {
                let slider = egui::Slider::new(pvalue_limit, 1.0..=400.0).text("pvalue_limit");
                ui.add(slider);

                if ui.button("RefAlt").clicked() {
                    *state = State::Alignment(AlignmentPlot::RefAlt);
                }
            } else if state == &mut State::Alignment(AlignmentPlot::RefAlt)
                && ui.button("P-value").clicked()
            {
                *state = State::Alignment(AlignmentPlot::Pvalue);
            }

            if ui.button("Alignment weight").clicked() {
                *state = State::AlignmentWeight;
            }

            ui.add_space(5.0);
            ui.heading("Files");
            if (ui.button("Recombination rates")).clicked() {
                let mut dialog = FileDialog::open_file(opened_file.clone());
                dialog.open();
                *open_file_dialog = Some(dialog);
                *filetype_to_open = FiletypeToOpen::RecombinationRates;
            }

            if (ui.button("GFF3")).clicked() {
                let mut dialog = FileDialog::open_file(opened_file.clone());
                dialog.open();
                *open_file_dialog = Some(dialog);
                *filetype_to_open = FiletypeToOpen::Gff3;
            }

            if let Some(dialog) = open_file_dialog {
                if dialog.show(ctx).selected() {
                    if let Some(file) = dialog.path() {
                        match filetype_to_open {
                            FiletypeToOpen::RecombinationRates => {
                                *rec_rates = Some(read_recombination_file(file).unwrap());
                            }
                            FiletypeToOpen::Gff3 => {
                                let gff = read_gtf_file_to_vec(file, vcf.get_contig()).unwrap();
                                *gff3 = Some(gff);
                            }
                        }
                    }
                }
            }
        })
}
