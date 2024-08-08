// Genome-wide HSTs
pub mod hst_gwas;
pub mod hst_mrca_gwas;
pub mod hst_scan;
// pub mod hst_qt_gwas;
pub mod hst_segregation;
pub use hst_scan::{
    get_marker_id, get_sender, read_tree_file, return_assoc, run, write_assoc, AssocRow, HstScan,
    HstScanRow, Metadata,
};
