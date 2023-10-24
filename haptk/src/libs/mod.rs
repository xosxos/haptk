// HATK - Haplotype analysis toolkit
// Copyright (C) 2022  Osma S. Rautila
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
//

// #![doc(
//     html_logo_url = "https://raw.githubusercontent.com/rust-bio/rust-bio/master/img/bioferris.svg",
//     html_favicon_url = "https://raw.githubusercontent.com/rust-bio/rust-bio/master/img/bioferris.svg"
// )]

//! HATK - Haplotype analysis toolkit
//!
//! This library and program provides a multitude of tools designed for the analysis of phased and
//! imputed genomes.
//! All provided tools are tested via continuous integration.
//!
//! For navigating the documentation of the available modules, see [the `Modules` section below](#modules).
//! If you want to contribute to `hatk`, see [the `Contribute` section in the repo]()
//!
//! HATK toolkit commands
//!
//! * Check for haplotype
//! * Find the longest common haplotype at a locus
//! * Contrahomozygosity
//! * Random sampled contrahomozygosity
//! * Decay of shared ancestral haplotype
//! * Recursive decay of shared ancestral haplotypes
//! * Most recent common ancestor estimation using the Gamma method
//! * Pairwise Identity-By-Descent at a specific locus
//! * Sample homozygosity
//! * Utils to check VCF coverage and sample names
//!
//! # Getting started
//!
//! Users who already know `HATK` might want to jump right into the [modules docs](https://docs.rs/bio/#modules)
//!
//! ## Installing HATK
//!
//! Rust and its package manager cargo can be installed following the instruction for [rustup](https://rustup.rs/).
//!
//! After installing cargo, run the following command
//!
//! ```bash
//! cargo install hatk
//! ```
//!
//! ## Running HATK
//!
//! To print the available commands use (more about the commands in the subcommnds section):
//! ```bash
//! hatk --help
//! ```
//! To run for example the longest common haplotype at a specific locus tool use:
//! ```bash
//! hatk common-haplotype $file --coords $coords > ${outdir}/longest_haplotype.tsv
//! ```
//!
//! For a full analysis of a given locus using only the longest shared haplotypes run the following:
//! ```bash
//!hatk common-haplotype $file --coords $coords --only-longest > ${outdir}/longest_haplotype.tsv
//!
//!hatk decay $file --coords $coords --only-longest > ${outdir}/decay.csv
//!
//!hatk multi-decay $file --coords $coords --only-longest -o "${outdir}/decay_graph"
//!
//!hatk mrca $file --coords $coords --only-longest > ${outdir}/age.txt
//!
//!hatk ibs $file --coords $coords --only-longest -o "${outdir}/ibs.png"
//!
//!hatk rand-homozygosity $file -t 8 --iters 300 > ${outdir}/rand_homozygosity.csv
//!
//!hatk sample-homozygosity \
//!  $file --coords $coords -t 8 --iters 100 > ${outdir}/sampled-homozygosity.csv
//!
//!hatk homozygosity $file -t 8 > ${outdir}/homozygosity.csv
//!
//!hatk check-for-haplotype $file --coords $coords --haplotype ${outdir}/longest_haplotype.tsv
//!
//!
//!```
//!

#[doc(hidden)]
pub mod args;

#[doc(hidden)]
pub mod io;

/// Functions for reading vcfs into matrices
pub mod read_vcf;

/// HATK structs
pub mod structs;

#[doc(hidden)]
pub mod utils;

#[doc(hidden)]
pub mod core;

#[doc(hidden)]
pub mod stats;

#[doc(hidden)]
pub mod error;

#[cfg(feature = "clap")]
pub mod clap;
