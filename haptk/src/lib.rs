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

pub mod libs;
pub use libs::{args, error, io, read_vcf, stats, structs, utils};

#[cfg(feature = "clap")]
pub use libs::clap;

/// Graphs implemented in Rust, HST graph will be replaced by the Python implementation
pub mod graphs;

/// HAPTK commands
pub mod subcommands;
