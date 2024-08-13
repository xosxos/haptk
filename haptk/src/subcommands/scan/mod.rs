// Genome-wide HSTs
pub mod hst_scan;
pub mod scan_branch_mrca;
pub mod scan_nodes;
pub mod scan_quantitative;
pub mod scan_segregate;
pub use hst_scan::{
    get_marker_id, get_sender, read_tree_file, return_assoc, run, top_node_from_hsts, write_assoc,
    zygosity_from_node, AssocRow, HstScan, HstScanRow, Limits, Metadata,
};
