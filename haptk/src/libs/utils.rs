use color_eyre::eyre::{eyre, OptionExt, WrapErr};
use color_eyre::Result;

use crate::error::HatkError::{CoordsParseError, PosParseError};

// Round to n significant digits
// https://stackoverflow.com/questions/28655362/how-does-one-round-a-floating-point-number-to-a-specified-number-of-digits
pub fn precision_f64(x: f64, decimals: u32) -> f64 {
    if x == 0. || decimals == 0 {
        0.
    } else {
        let shift = decimals as i32 - x.abs().log10().ceil() as i32;
        let shift_factor = 10_f64.powi(shift);

        (x * shift_factor).round() / shift_factor
    }
}

//NOTE: This should be parsed by clap automatically, but Option<String> parsing is not supported out of the box as of now
pub fn strip_prefix(prefix: Option<String>) -> Option<String> {
    if let Some(prefix) = prefix {
        match prefix.as_ref() {
            "" => None,
            "\\0" => None,
            v => Some(v.to_string()),
        }
    } else {
        None
    }
}

// Coords are in the format [contig]:[position]
pub fn parse_snp_coord(coords: &str) -> Result<(&str, u64)> {
    let mut coord_split = coords.split(':');

    match (coord_split.next(), coord_split.next()) {
        (Some(contig), Some(value)) => {
            let value = value
                .parse::<u64>()
                .wrap_err(eyre!(PosParseError((coords.into(), value.into()))))?;
            Ok((contig, value))
        }
        _ => Err(eyre!(CoordsParseError(coords.into()))),
    }
}

// Coords are in the format [contig] or [contig]:[start]-[stop]
pub fn parse_coords(coords: &str) -> Result<(&str, Option<u64>, Option<u64>)> {
    let mut coord_split = coords.split(':');

    let contig = coord_split
        .next()
        .ok_or_eyre(CoordsParseError(coords.into()))?;

    let positions = coord_split.next();

    if positions.is_none() {
        return Ok((contig, None, None));
    }

    let mut pos_split = positions.unwrap().split('-');
    let (start, stop) = (pos_split.next(), pos_split.next());

    if let (Some(start), Some(stop)) = (start, stop) {
        let start = start
            .parse::<u64>()
            .wrap_err(eyre!(PosParseError((coords.into(), start.into()))))?;
        let stop = stop
            .parse::<u64>()
            .wrap_err(eyre!(PosParseError((coords.into(), stop.into()))))?;

        Ok((contig, Some(start), Some(stop)))
    } else {
        Err(eyre!(CoordsParseError(coords.into())))
    }
}

pub fn centromeres_hg38(chr: &str) -> (u64, u64) {
    match chr {
        "chr1" => (121700000, 125100000),
        "chr2" => (91800000, 96000000),
        "chr3" => (87800000, 94000000),
        "chr4" => (48200000, 51800000),
        "chr5" => (46100000, 51400000),
        "chr6" => (58500000, 62600000),
        "chr7" => (58100000, 62100000),
        "chr8" => (43200000, 47200000),
        "chr9" => (42200000, 45500000),
        "chr10" => (38000000, 41600000),
        "chr11" => (51000000, 55800000),
        "chr12" => (33200000, 37800000),
        "chr13" => (16500000, 18900000),
        "chr14" => (16100000, 18200000),
        "chr15" => (17500000, 20500000),
        "chr16" => (35300000, 38400000),
        "chr17" => (22700000, 27400000),
        "chr18" => (15400000, 21500000),
        "chr19" => (24200000, 28100000),
        "chr20" => (25700000, 30400000),
        "chr21" => (10900000, 13000000),
        "chr22" => (13700000, 17400000),
        "chrX" => (58100000, 63800000),
        "chrY" => (10300000, 10600000),
        _ => panic!(),
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;


    #[test]
    fn test_parse_snp_coord() {
        let (chr, pos) = parse_snp_coord("chr9:1920").unwrap();
        assert_eq!(chr, "chr9");
        assert_eq!(pos, 1920);
    }

    #[test]
    fn test_parse_coord() {
        let (contig, start, stop) = parse_coords("chr9").unwrap();
        assert_eq!(contig, "chr9");
        assert_eq!(start, None);
        assert_eq!(stop, None);

        let res = parse_coords("chr9:1920");
        assert!(res.is_err());

        let (contig, start, stop) = parse_coords("chr9:1920-2500").unwrap();
        assert_eq!(contig, "chr9");
        assert_eq!(start, Some(1920));
        assert_eq!(stop, Some(2500));

        let result = parse_coords("chr9:1920--2500");
        assert!(result.is_err())
    }

    #[test]
    #[ignore]
    #[should_panic(
        expected = "called `Result::unwrap()` on an `Err` value: Position \"foo\" is not an integer in coords \"chr9:foo\"\n\nCaused by:\n    invalid digit found in string\n\nLocation:\n    src/libs/core.rs:175:38"
    )]
    fn test_parse_snp_coord_int_parser_error() {
        parse_snp_coord("chr9:foo").unwrap();
    }
}
