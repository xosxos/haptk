use std::rc::Rc;

use eframe::egui;
use petgraph::Graph;

use haptk::structs::PhasedMatrix;
use haptk::subcommands::bhst_shard::Node;

use crate::central_panel::central_panel;
use crate::side_panel::side_panel;

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum State {
    Bhst,
    UhstLeft,
    UhstRight,
}

pub type Hst<'a> = Graph<Node<'a>, u8>;

pub struct TreeApp<'a> {
    vcf: Rc<PhasedMatrix>,
    state: State,
    bhst: Rc<Hst<'a>>,
    uhst_left: Rc<Hst<'a>>,
    uhst_right: Rc<Hst<'a>>,
    nmin_samples: usize,
    decoy_samples: Rc<Vec<String>>,
}

impl<'a> TreeApp<'a> {
    pub fn new(
        vcf: PhasedMatrix,
        bhst: Hst<'a>,
        uhst_left: Hst<'a>,
        uhst_right: Hst<'a>,
        nmin_samples: usize,
        decoy_samples: Option<Vec<String>>,
    ) -> Self {
        let decoy_samples = decoy_samples.map(Rc::new).unwrap_or(Rc::new(vec![]));
        Self {
            vcf: Rc::new(vcf),
            bhst: Rc::new(bhst),
            uhst_left: Rc::new(uhst_left),
            uhst_right: Rc::new(uhst_right),
            state: State::Bhst,
            nmin_samples,
            decoy_samples,
        }
    }
}

impl<'a> eframe::App for TreeApp<'a> {
    /// Called by the frame work to save state before shutdown.
    /// Note that you must enable the `persistence` feature for this to work.
    #[cfg(feature = "persistence")]
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(storage, eframe::APP_KEY, self);
    }

    /// Called each time the UI needs repainting, which may be many times per second.
    /// Put your widgets into a `SidePanel`, `TopPanel`, `CentralPanel`, `Window` or `Area`.
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let TreeApp {
            state,
            vcf,
            bhst,
            uhst_left,
            uhst_right,
            nmin_samples,
            decoy_samples,
            ..
        } = self;

        side_panel(ctx, state);

        central_panel(
            ctx,
            state,
            vcf.clone(),
            bhst.clone(),
            uhst_left.clone(),
            uhst_right.clone(),
            *nmin_samples,
            decoy_samples.clone(),
        );

        // ctx.request_repaint();
    }
}
