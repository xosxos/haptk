use std::collections::BTreeMap;
use std::path::PathBuf;

use color_eyre::{eyre::eyre, Result};
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Gamma;
use statrs::statistics::Statistics;

use crate::{
    args::{Selection, StandardArgs},
    io::{push_to_output, read_recombination_file},
    libs::read_vcf::read_vcf_to_matrix,
    subcommands::bhst::Node,
    utils::parse_snp_coord,
};

pub type Age = (f64, f64, f64);

#[doc(hidden)]
pub fn run(args: StandardArgs, rec_rates: PathBuf, window: Option<u64>) -> Result<()> {
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    let (contig, variant_pos) = parse_snp_coord(&args.coords)?;

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "mrca_gamma_method", "txt");

    let rates = read_recombination_file(rec_rates)?;

    let mut vcf = read_vcf_to_matrix(&args, contig, variant_pos, None, None, window, false)?;

    if args.selection == Selection::OnlyLongest {
        vcf.select_only_longest()?;
    };

    let start = vcf.start_coord().clone();
    let lengths = vcf.get_lengths_from_uhst(&start)?;

    let ((i_tau_hat, i_l, i_u), (c_tau_hat, c_l, c_u)) =
        mrca_gamma_method(lengths, variant_pos, &rates)?;

    let data = format!("Independent genealogy:\nage: {i_tau_hat:.3} CI ({i_l:.3}, {i_u:.3})\nCorrelated genealogy:\nage: {c_tau_hat:.3} CI ({c_l:.3}, {c_u:.3})");

    std::fs::write(&output, data).unwrap_or_else(|_| panic!("Unable to write to {output:?}"));

    Ok(())
}

pub fn get_sums(
    shared_lengths: Vec<(Node, Node)>,
    variant_pos: u64,
    rec_rates: &BTreeMap<u64, f32>,
) -> (Vec<f64>, Vec<f64>) {
    let l_sum = shared_lengths
        .iter()
        .map(|(lnode, _)| lnode.stop.pos - lnode.start.pos)
        .sum::<u64>();

    let r_sum = shared_lengths
        .iter()
        .map(|(_, rnode)| rnode.stop.pos - rnode.start.pos)
        .sum::<u64>();

    tracing::debug!("Sums bp: left {l_sum:.3}, right {r_sum:.3})",);

    let (_var_pos, var_cm) = rec_rates.range(variant_pos..).next().unwrap_or_else(|| {
            let rate = rec_rates.range(..variant_pos).next_back().unwrap();
            tracing::warn!("Recombination rate data ends at {}, but a marker has position {}. No centimorgans are added after the end of recombination rate data.", rate.0, variant_pos);
            rate
    });

    // Transform from positions to centimorgans by selecting the nearest values from the recombination map
    shared_lengths
        .iter()
        .map(|(lnode, rnode)| (lnode.start.pos, rnode.stop.pos))
        .fold((vec![], vec![]), |mut acc, (start, stop)| {
            let (_, start_cm) = rec_rates
                .range(..=start)
                .next_back()
                .unwrap_or_else(|| rec_rates.range(start..).nth(1).unwrap());

            let (_, stop_cm) = rec_rates
                .range(stop..)
                .next()
                .unwrap_or_else(|| rec_rates.range(..stop).next_back().unwrap());

            let left_len = (var_cm - start_cm) as f64 / 100.;
            let right_len = (stop_cm - var_cm) as f64 / 100.;

            acc.0.push(left_len);
            acc.1.push(right_len);
            acc
        })
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

    match Gamma::new(2.0 * n, 2.0 * n * b_c) {
        Ok(gamma) => {
            let g_l = gamma.inverse_cdf((1.0 - cc) / 2.0);
            let g_u = gamma.inverse_cdf(cc + (1.0 - cc) / 2.0);

            let i_tau_hat_l = i_tau_hat * g_l;
            let i_tau_hat_u = i_tau_hat * g_u;

            Ok((i_tau_hat, i_tau_hat_l, i_tau_hat_u))
        }
        Err(e) => {
            tracing::error!(
                "Could not estimate the mrca for {} samples: {e:?}",
                l_lengths.len()
            );
            tracing::error!("Returning all zeros for the MRCA estimate");
            Ok((0.0, 0.0, 0.0))
        }
    }
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

    match Gamma::new(2.0 * n_star, 2.0 * n_star * b_c) {
        Ok(gamma) => {
            let c_l = gamma.inverse_cdf((1.0 - cc) / 2.0);
            let c_u = gamma.inverse_cdf(cc + (1.0 - cc) / 2.0);

            let c_tau_hat_l = c_tau_hat * c_l;
            let c_tau_hat_u = c_tau_hat * c_u;

            Ok((c_tau_hat, c_tau_hat_l, c_tau_hat_u))
        }
        Err(e) => {
            tracing::error!(
                "Could not estimate the mrca for {} samples: {e:?}",
                l_lengths.len()
            );
            tracing::error!("Returning all zeros for the MRCA estimate");
            Ok((0.0, 0.0, 0.0))
        }
    }
}

///
/// The original R algorithm by Gandolfo et al translated to Rust.
/// <https://github.com/bahlolab/DatingRareMutations>
///
pub fn mrca_gamma_method(
    shared_lengths: Vec<(Node, Node)>,
    variant_pos: u64,
    rec_rates: &BTreeMap<u64, f32>,
) -> Result<(Age, Age)> {
    let cc = 0.95;
    let cs_corr = 0.0;
    let n = shared_lengths.len() as f64;

    let (l_lengths, r_lengths) = get_sums(shared_lengths, variant_pos, rec_rates);

    let (print_left_tau_hat, _, _) = independent(n, &l_lengths, &[0.1], cc, cs_corr)?;
    let (print_right_tau_hat, _, _) = independent(n, &[0.1], &r_lengths, cc, cs_corr)?;

    tracing::debug!(
        "Independent. Left age: {print_left_tau_hat:.3}, right age: {print_right_tau_hat:.3})",
    );

    let (i_tau_hat, i_tau_hat_l, i_tau_hat_u) =
        independent(n, &l_lengths, &r_lengths, cc, cs_corr)?;

    let (c_tau_hat, c_tau_hat_l, c_tau_hat_u) = correlated(l_lengths, r_lengths, cc, cs_corr)?;

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
