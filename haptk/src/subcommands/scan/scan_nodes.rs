use std::{collections::BTreeMap, path::PathBuf, sync::Arc};

use color_eyre::Result;
use petgraph::graph::NodeIndex;

use crate::{
    args::{ConciseArgs, StandardArgs},
    io::{open_csv_writer, push_to_output},
    structs::Coord,
    subcommands::scan::{read_tree_file, return_assoc, write_assoc, AssocRow, HstScan, Limits},
    utils::{centromeres_hg38, parse_coords},
};

use super::Hst;

const HEADER: &[&str] = &[
    "contig",
    "pos",
    "marker_id",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "opt_value",
    "centromere",
];

fn check_for_centromere(contig: &str, start: u64, stop: u64) -> bool {
    let (c_start, c_stop) = centromeres_hg38(contig);

    let c1 = start > c_start && start < c_stop;
    let c2 = stop > c_start && stop < c_stop;
    let c3 = start < c_start && stop > c_stop;

    c1 || c2 || c3
}

#[doc(hidden)]
pub fn run(args: ConciseArgs, limits: Limits) -> Result<()> {
    let hsts = read_tree_file(args.file)?;
    let hsts = Arc::new(hsts);

    let args = StandardArgs {
        file: PathBuf::new(),
        output: args.output,
        info_limit: None,
        coords: hsts.metadata.contig.clone(),
        selection: hsts.metadata.selection.clone(),
        prefix: args.prefix,
        samples: None,
    };

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "node_scan", "csv");
    let writer = open_csv_writer(output)?;

    let (_contig, _start, _stop) = parse_coords(&args.coords)?;

    // Define the minimizer/maximiser function for the HST
    let optimizer =
        |hsts: Arc<HstScan>,
         hst: &Hst,
         coord: &Coord,
         limits: Limits|
         -> Option<(Coord, NodeIndex, f64)> { optimizer_inner(hsts, hst, coord, limits) };

    // Define what information is wanted for each CSV row
    let rower = |hsts: Arc<HstScan>, coord: &Coord, top_node_idx, optimized_value| -> AssocRow {
        rower_inner(hsts, coord, top_node_idx, optimized_value)
    };

    let write_bam = false;
    let assoc = return_assoc(&hsts, &args, limits, write_bam, optimizer, rower);

    write_assoc(HEADER, assoc, writer)?;
    Ok(())
}

fn rower_inner(
    hsts: Arc<HstScan>,
    coord: &Coord,
    top_node_idx: NodeIndex,
    optimized_value: f64,
) -> AssocRow {
    let hst = hsts.hsts.get(coord).unwrap();
    let top_node = hst.node_weight(top_node_idx).unwrap();
    let start = top_node.start.pos;
    let stop = top_node.stop.pos;

    let is_centromere = check_for_centromere(&coord.contig, start, stop);

    let mut assoc_hm = BTreeMap::new();
    assoc_hm.insert("centromere".to_string(), is_centromere.to_string());

    AssocRow::new(hsts, coord, top_node_idx, optimized_value, assoc_hm)
}

fn optimizer_inner(
    _hsts: Arc<HstScan>,
    hst: &Hst,
    coord: &Coord,
    limits: Limits,
) -> Option<(Coord, NodeIndex, f64)> {
    let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;
    let mut top_node_idx = NodeIndex::new(0);

    // START VALUE
    let start_value = f64::MIN;

    let mut optimized_value = start_value;

    let mut iterator = hst.node_indices();

    // Skip first node because it contains no haplotypes
    iterator.next();

    for idx in iterator {
        let node = hst.node_weight(idx).unwrap();

        if node.indexes.len() >= nmin_samples
            && node.indexes.len() <= nmax_samples
            && node.haplotype.len() >= nmin_variants
            && node.haplotype.len() <= nmax_variants
        {
            // OPTIMIZER CODE
            let node = hst.node_weight(idx).unwrap();

            // Longest haplotype in bp
            let start = node.start.pos;
            let stop = node.stop.pos;
            let value = stop
                .checked_sub(start)
                .expect("Overflowed subtraction, report to HAPTK github")
                as f64;

            // Longest haplotype in markers
            // let value = node.haplotype.len() as f64;

            // Clause
            let clause = value > optimized_value;

            // OPTIMIZER CODE END

            if clause {
                optimized_value = value;
                top_node_idx = idx;
            }
        }
    }

    if optimized_value == start_value {
        tracing::warn!(
            "No optimized value for the HST in contig {} at pos {}",
            coord.contig,
            coord.pos,
        );
        return None;
    }

    Some((coord.clone(), top_node_idx, optimized_value))
}
