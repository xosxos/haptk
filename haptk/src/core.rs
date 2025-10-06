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

/// Functions for reading vcfs into matrices
pub mod read_vcf;

/// Input args
pub mod args;

/// Errors
pub mod error;

/// IO functions
pub mod io;

/// Statistics, currently not used
pub mod stats;

/// Structs
pub mod structs;

/// Utils
// pub mod utils;
pub mod clap;

pub use haptk_core::bam;

pub use haptk_core::utils;

// #[cfg(feature = "bam")]
// pub mod bam;

pub mod phased_matrix;

pub mod cigar_iterator;

pub use phased_matrix::PhasedMatrix;
pub use structs::CigarVariant;
pub use structs::Coord;
pub use structs::HapVariant;
pub use structs::Ploidy;
