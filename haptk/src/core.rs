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

/// Clap API for HAPTK
pub mod clap;

// Re-export core modules
pub use haptk_core::bam;
pub use haptk_core::cigar_iterator;
pub use haptk_core::phased_matrix::Matrix;
pub use haptk_core::phased_matrix::PhasedMatrix;
pub use haptk_core::phased_matrix::ReadMetadata;
pub use haptk_core::ploidy::Ploidy;
pub use haptk_core::utils;
pub use haptk_core::variant::CigarVariant;
pub use haptk_core::variant::Coord;
pub use haptk_core::variant::HapVariant;
