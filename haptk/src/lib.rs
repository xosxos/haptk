#![allow(
    clippy::large_enum_variant,
    clippy::too_many_arguments,
    clippy::new_without_default,
    clippy::uninlined_format_args,
    clippy::missing_errors_doc,
    clippy::unreadable_literal,
    clippy::too_many_lines,
    clippy::must_use_candidate,
    clippy::enum_glob_use,
    clippy::missing_panics_doc,
    clippy::module_name_repetitions,
    clippy::match_bool,
    clippy::single_match_else,
    clippy::cast_possible_wrap,
    clippy::return_self_not_must_use,
    clippy::used_underscore_binding,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::wildcard_in_or_patterns,
    clippy::needless_pass_by_value,
    clippy::default_trait_access,
    clippy::struct_field_names,
    clippy::unused_self
)]

// HATK - Haplotype analysis toolkit
// Copyright (C) 2024  Osma S. Rautila
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
