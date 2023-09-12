use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::Arc;

use color_eyre::Result;
use eframe::egui;
use egui::plot::PlotBounds;
use egui_file::FileDialog;

use haptk::args::StandardArgs;
use haptk::structs::PhasedMatrix;
use haptk::subcommands::annotate_haplotype::GtfRow;
use haptk::subcommands::bhst::Node;
use haptk::subcommands::hst_gwas::{
    read_tree_file, read_vcfs, return_binary_assoc, return_mbah_lengths, Assoc,
};
use haptk::subcommands::hst_scan::{TreeRow, Trees};

use crate::alignment::{return_alignment, AlignmentPoint};
use crate::central_panel::central_panel;
use crate::side_panel::side_panel;
use crate::tree::TreePoint;
use crate::utils::assoc_vec_to_bm;
use crate::Args;

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum AlignmentPlot {
    Pvalue,
    RefAlt,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum MbahPlot {
    NumberOfMarkers,
    BasePairs,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum FiletypeToOpen {
    RecombinationRates,
    Gff3,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum State {
    MbahLengths(MbahPlot),
    NumberOfMarkers,
    Association,
    BlockLength,
    TopNodeMrca,
    Alignment(AlignmentPlot),
    HaplotypeStateTree(usize),
    AlignmentWeight,
}

pub struct Parameters {
    min_samples: usize,
    max_samples: usize,
    min_ht_len: usize,
    max_ht_len: usize,
}

#[derive(Clone, Default, Eq, PartialEq, Hash)]
pub struct Key {
    pub min_samples: usize,
    pub max_samples: usize,
    pub min_ht_len: usize,
    pub max_ht_len: usize,
}

pub type Ps = Vec<[f64; 2]>;
pub type TreeCache = Rc<RefCell<BTreeMap<usize, (BTreeSet<TreePoint>, Vec<Ps>)>>>;
pub type AlignmentCache = HashMap<(Key, usize), (Vec<AlignmentPoint>, Ps)>;
pub type AssocCache = HashMap<Key, Rc<BTreeMap<u64, Box<dyn Assoc + Send>>>>;
pub type TopNodeMrcaCache = HashMap<Key, Vec<(u64, f64)>>;

pub struct GwasApp {
    vcf: Arc<PhasedMatrix>,
    ctrl_vcf: PhasedMatrix,
    trees: Rc<Vec<TreeRow>>,
    assoc_cache: AssocCache,
    alignment_cache: AlignmentCache,
    tree_cache: TreeCache,
    top_node_mrca_cache: TopNodeMrcaCache,
    state: State,
    parameters: Parameters,
    pvalue_limit: f64,
    pvalue: usize,
    alignment_bounds_set: bool,
    local_bounds_set: bool,
    alignment_bounds: PlotBounds,
    tree_mbah_lengths: Vec<(usize, Node)>,
    opened_file: Option<PathBuf>,
    open_file_dialog: Option<FileDialog>,
    rec_rates: Option<BTreeMap<u64, f32>>,
    gff3: Option<Vec<GtfRow>>,
    filetype_to_open: FiletypeToOpen,
}

impl GwasApp {
    pub fn new(args: Args) -> Result<Self> {
        let now = std::time::Instant::now();
        let Trees { metadata, trees } = read_tree_file(args.trees)?;
        let trees = Rc::new(trees);

        let mut st_args = StandardArgs {
            file: args.args.file,
            output: args.args.outdir,
            info_limit: metadata.info_limit,
            coords: metadata.coords,
            selection: metadata.selection.clone(),
            prefix: args.args.prefix,
            samples: None,
        };

        tracing::debug!("Loaded trees at {:?}", now.elapsed());

        let (vcf, ctrl_vcf) = read_vcfs(
            &mut st_args,
            args.controls,
            args.control_vcf,
            metadata.samples,
        )?;

        tracing::debug!("Loaded vcfs at {:?}", now.elapsed());

        let min_samples = 10;
        let max_samples = 200;
        let min_ht_len = 1;
        let max_ht_len = 20000;
        let pvalue = 6;

        let vcf = Arc::new(vcf);

        let assoc = return_binary_assoc(
            trees.as_ref(),
            vcf.clone(),
            (min_samples, max_samples, min_ht_len, min_ht_len),
            &ctrl_vcf,
            None,
        );

        let assoc = Rc::new(assoc_vec_to_bm(assoc));

        tracing::info!("Read assoc at {:?}", now.elapsed());

        let (points, weights) = return_alignment(assoc.clone(), vcf.clone(), pvalue);

        tracing::info!("Created alignment at {:?}", now.elapsed());

        let mbahs = return_mbah_lengths(&trees);
        tracing::info!("Calculated MBAH lengths at {:?}", now.elapsed());

        let mut assoc_cache = HashMap::new();
        let mut alignment_cache = HashMap::new();

        let key = Key {
            min_samples,
            max_samples,
            min_ht_len,
            max_ht_len,
        };

        assoc_cache.insert(key.clone(), assoc);
        alignment_cache.insert((key, pvalue), (points, weights));

        Ok(Self {
            vcf,
            ctrl_vcf,
            trees,
            parameters: Parameters {
                min_samples,
                max_samples,
                min_ht_len,
                max_ht_len,
            },
            assoc_cache,
            alignment_cache,
            state: State::Alignment(AlignmentPlot::Pvalue),
            pvalue_limit: 8.0,
            pvalue: 6,
            alignment_bounds_set: false,
            local_bounds_set: false,
            alignment_bounds: PlotBounds::from_min_max([0.0, 0.0], [0.0, 0.0]),
            tree_mbah_lengths: mbahs,
            tree_cache: Rc::new(RefCell::new(BTreeMap::new())),
            opened_file: None,
            open_file_dialog: None,
            rec_rates: None,
            top_node_mrca_cache: HashMap::new(),
            gff3: None,
            filetype_to_open: FiletypeToOpen::RecombinationRates,
        })
    }
}

impl eframe::App for GwasApp {
    /// Called by the frame work to save state before shutdown.
    /// Note that you must enable the `persistence` feature for this to work.
    #[cfg(feature = "persistence")]
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(storage, eframe::APP_KEY, self);
    }

    /// Called each time the UI needs repainting, which may be many times per second.
    /// Put your widgets into a `SidePanel`, `TopPanel`, `CentralPanel`, `Window` or `Area`.
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let GwasApp {
            state,
            pvalue_limit,
            pvalue,
            opened_file,
            open_file_dialog,
            parameters,
            alignment_bounds,
            alignment_bounds_set,
            local_bounds_set,
            vcf,
            ctrl_vcf,
            trees,
            tree_mbah_lengths,
            assoc_cache,
            alignment_cache,
            tree_cache,
            top_node_mrca_cache,
            rec_rates,
            gff3,
            filetype_to_open,
            ..
        } = self;

        let Parameters {
            min_samples,
            max_samples,
            min_ht_len,
            max_ht_len,
        } = parameters;

        side_panel(
            ctx,
            vcf.clone(),
            min_samples,
            max_samples,
            min_ht_len,
            max_ht_len,
            pvalue_limit,
            pvalue,
            rec_rates,
            gff3,
            state,
            open_file_dialog,
            opened_file,
            filetype_to_open,
        );

        central_panel(
            ctx,
            *min_samples,
            *max_samples,
            *min_ht_len,
            *max_ht_len,
            *pvalue,
            *pvalue_limit,
            state,
            alignment_bounds,
            alignment_bounds_set,
            local_bounds_set,
            vcf.clone(),
            ctrl_vcf,
            trees,
            tree_mbah_lengths,
            assoc_cache,
            alignment_cache,
            tree_cache,
            top_node_mrca_cache,
            rec_rates,
            gff3,
        );

        // ctx.request_repaint();
    }
}
