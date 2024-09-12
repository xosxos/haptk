use std::collections::BTreeMap;
use std::path::PathBuf;

use color_eyre::{eyre::eyre, Result};
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Gamma;
use statrs::statistics::Statistics;

use crate::{
    args::{Selection, StandardArgs},
    io::{push_to_output, read_recombination_file},
    read_vcf::read_vcf_to_matrix,
    structs::PhasedMatrix,
    subcommands::bhst::Node,
    utils::parse_snp_coord,
};

pub type Age = (f64, f64, f64);

#[doc(hidden)]
pub fn run(args: StandardArgs, rec_rates: PathBuf) -> Result<()> {
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    let (contig, variant_pos) = parse_snp_coord(&args.coords)?;

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "mrca_gamma_method", "txt");

    let rates = read_recombination_file(rec_rates)?;

    let mut vcf = read_vcf_to_matrix(&args, contig, variant_pos, None, None)?;

    if args.selection == Selection::OnlyAlts || args.selection == Selection::OnlyRefs {
        vcf.select_carriers(variant_pos, &args.selection)?;
    } else if args.selection == Selection::OnlyLongest {
        vcf.select_only_longest();
    };

    let lengths = vcf.get_lengths_from_uhst(vcf.variant_idx());
    let ((i_tau_hat, i_l, i_u), (c_tau_hat, c_l, c_u)) =
        mrca_gamma_method(&vcf, lengths, variant_pos, &rates)?;

    let data = format!("Independent genealogy:\nage: {i_tau_hat:.3} CI ({i_l:.3}, {i_u:.3})\nCorrelated genealogy:\nage: {c_tau_hat:.3} CI ({c_l:.3}, {c_u:.3})");

    std::fs::write(&output, data).unwrap_or_else(|_| panic!("Unable to write to {output:?}"));

    Ok(())
}

pub fn get_sums(
    vcf: &PhasedMatrix,
    shared_lengths: Vec<(Node, Node)>,
    variant_pos: u64,
    rec_rates: &BTreeMap<u64, f32>,
) -> (Vec<f64>, Vec<f64>) {
    // let l_sum = shared_lengths
    //     .iter()
    //     .map(|(lnode, _)| vcf.get_pos(lnode.stop_idx) - vcf.get_pos(lnode.start_idx))
    //     .sum::<u64>();
    // let r_sum = shared_lengths
    //     .iter()
    //     .map(|(_, rnode)| vcf.get_pos(rnode.stop_idx) - vcf.get_pos(rnode.start_idx))
    //     .sum::<u64>();
    // tracing::debug!("Sums bp: left {l_sum:.3}, right {r_sum:.3})",);

    let variant_cm = match rec_rates.range(variant_pos..).next() {
        None => {
            let last_rate = rec_rates.range(..variant_pos).next_back().unwrap();
            tracing::warn!("Recombination rate data ends at {}, but a marker has position {}. No centimorgans are added after the end of recombination rate data.", last_rate.0, variant_pos);
            last_rate
        }
        Some(variant_cm) => variant_cm,
    };

    let mut l_lengths = shared_lengths
        .iter()
        .map(|(lnode, _)| vcf.get_pos(lnode.start_idx))
        .map(|start| {
            // Transform from pos to centimorgans by selecting nearest value to the left in the
            // BTreeMap
            let break_cm = if let Some(cm) = rec_rates.range(..start).next_back() {
                cm
            } else {
                rec_rates.range(start..).next().unwrap()
            };

            (variant_cm.1 - break_cm.1) as f64 / 100.
        })
        .collect::<Vec<f64>>();

    let r_lengths = shared_lengths
        .iter()
        .map(|(_, rnode)| {
            if rnode.stop_idx == vcf.ncoords() {
                vcf.get_pos(rnode.stop_idx - 1)
            } else {
                vcf.get_pos(rnode.stop_idx)
            }
        })
        .map(|stop| {
            let break_cm = if let Some(cm) = rec_rates.range(stop..).next() {
                cm
            } else {
                rec_rates.range(..stop).next_back().unwrap()
            };

            (break_cm.1 - variant_cm.1) as f64 / 100.
        })
        .collect::<Vec<f64>>();

    (l_lengths, r_lengths)
}

pub fn independent(
    n: f64,
    l_lengths: &[f64],
    r_lengths: &[f64],
    cc: f64,
    cs_corr: f64,
) -> Result<Age> {
    let l_sum: f64 = l_lengths.iter().sum();
    let r_sum: f64 = r_lengths.iter().sum();

    // Independent genealogy
    let length_corr = (l_sum + r_sum - 2.0 * (n - 1.0) * cs_corr) / (2.0 * n);

    // Sum of ancestral segment lengths with corrections
    let sum: f64 = l_sum + r_sum + 2.0 * length_corr - 2.0 * (n - 1.0) * cs_corr;

    // Gamma function MLE bias correction factor
    let b_c = (2.0 * n - 1.0) / (2.0 * n);

    // Minimum variance unbiased estimate of T
    // Compare it with the pure MLE i.e. the mean of the gamma distribution = α / λ
    let i_tau_hat = (b_c * 2.0 * n) / sum;

    let gamma = Gamma::new(2.0 * n, 2.0 * n * b_c)?;
    let g_l = gamma.inverse_cdf((1.0 - cc) / 2.0);
    let g_u = gamma.inverse_cdf(cc + (1.0 - cc) / 2.0);

    let i_tau_hat_l = i_tau_hat * g_l;
    let i_tau_hat_u = i_tau_hat * g_u;

    Ok((i_tau_hat, i_tau_hat_l, i_tau_hat_u))
}

// Correlated genealogy
pub fn correlated(
    mut l_lengths: Vec<f64>,
    mut r_lengths: Vec<f64>,
    cc: f64,
    cs_corr: f64,
) -> Result<Age> {
    let l_sum: f64 = l_lengths.iter().sum::<f64>();
    let r_sum: f64 = r_lengths.iter().sum::<f64>();
    tracing::debug!("Sums: left {l_sum:.3}, right {r_sum:.3})",);
    let n = l_lengths.len() as f64;

    let length_corr: f64 = (l_sum + r_sum - 2.0 * (n - 1.0) * cs_corr) / (2.0 * n);

    let highest_l = l_lengths.iter_mut().max_by(|a, b| a.total_cmp(b)).unwrap();
    *highest_l = *highest_l + length_corr + cs_corr;

    let highest_r = r_lengths.iter_mut().max_by(|a, b| a.total_cmp(b)).unwrap();
    *highest_r = *highest_r + length_corr + cs_corr;

    let lengths = &l_lengths
        .iter()
        .zip(r_lengths.iter())
        .map(|(a, b)| a + b - 2.0 * cs_corr)
        .collect::<Vec<f64>>();

    let sum: f64 = lengths.iter().sum();

    let term1 = n * lengths.mean().powi(2) + lengths.variance() * (n - 1.0);
    let term2 = n * lengths.mean().powi(2) - lengths.variance() * (1.0 + 2.0 * n);

    let rho_hat = term2 / term1;

    let mut n_star = n / (1.0 + (n - 1.0) * rho_hat);

    if n_star > n {
        n_star = n;
    }
    if n_star < -n {
        n_star = -n;
    }

    let b_c = (2.0 * n_star - 1.0) / (2.0 * n_star);
    let c_tau_hat = (b_c * 2.0 * n) / sum;

    if -2.0 / (n - 1.0) <= rho_hat && rho_hat < -1.0 / (n - 1.0) {
        n_star = n;
    } else if rho_hat < -2.0 / (n - 1.0) {
        n_star = n / (1.0 + (n - 1.0) * rho_hat.abs());
    };

    let gamma = Gamma::new(2.0 * n_star, 2.0 * n_star * b_c)?;
    let c_l = gamma.inverse_cdf((1.0 - cc) / 2.0);
    let c_u = gamma.inverse_cdf(cc + (1.0 - cc) / 2.0);

    let c_tau_hat_l = c_tau_hat * c_l;
    let c_tau_hat_u = c_tau_hat * c_u;

    Ok((c_tau_hat, c_tau_hat_l, c_tau_hat_u))
}

///
/// The original R algorithm by Gandolfo et al translated to Rust.
/// <https://github.com/bahlolab/DatingRareMutations>
///
pub fn mrca_gamma_method(
    vcf: &PhasedMatrix,
    shared_lengths: Vec<(Node, Node)>,
    variant_pos: u64,
    rec_rates: &BTreeMap<u64, f32>,
) -> Result<(Age, Age)> {
    let cc = 0.95;
    let cs_corr = 0.0;
    let n = shared_lengths.len() as f64;

    let (l_lengths, r_lengths) = get_sums(vcf, shared_lengths, variant_pos, rec_rates);

    let (left_tau_hat, _, _) = independent(n, &l_lengths, &[0.1], cc, cs_corr)?;
    let (right_tau_hat, _, _) = independent(n, &[0.1], &r_lengths, cc, cs_corr)?;

    let (i_tau_hat, i_tau_hat_l, i_tau_hat_u) =
        independent(n, &l_lengths, &r_lengths, cc, cs_corr)?;
    let (c_tau_hat, c_tau_hat_l, c_tau_hat_u) = correlated(l_lengths, r_lengths, cc, cs_corr)?;

    tracing::debug!("Independent. Left age: {left_tau_hat:.3}, right age: {right_tau_hat:.3})",);

    tracing::debug!(
        "Independent genealogy:\nage: {i_tau_hat:.3} CI ({i_tau_hat_l:.3}, {i_tau_hat_u:.3})",
    );

    tracing::debug!(
        "Correlated genealogy:\nage: {c_tau_hat:.3} CI ({c_tau_hat_l:.3}, {c_tau_hat_u:.3})",
    );

    let independent = (i_tau_hat, i_tau_hat_l, i_tau_hat_u);
    let correlated = (c_tau_hat, c_tau_hat_l, c_tau_hat_u);

    Ok((independent, correlated))
}
