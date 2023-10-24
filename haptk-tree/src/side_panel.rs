use eframe::egui;
use egui::Margin;

use crate::app::State;

pub fn side_panel(ctx: &egui::Context, state: &mut State) -> egui::InnerResponse<()> {
    egui::SidePanel::left("sidebar")
        .default_width(200.0)
        .frame(egui::Frame {
            fill: egui::Color32::WHITE,
            inner_margin: Margin::symmetric(5.0, 5.0),
            ..Default::default()
        })
        .show(ctx, |ui| {
            ui.heading("HSTs");

            if ui.button("BHST").clicked() {
                *state = State::Bhst;
            }
            if ui.button("UHST left").clicked() {
                *state = State::UhstLeft;
            }
            if ui.button("UHST right").clicked() {
                *state = State::UhstRight;
            }
        })
}
