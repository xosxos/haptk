#[derive(Debug)]
pub enum HaptkError {
    CoordsParseError(String),
    PosParseError((String, String)),
    NormalizeError(u64),
    PloidyError((u64, usize)),
    SamplesNotFoundError,
    VariantPosNotFoundError((u64, usize)),
    MissingGenotypesError(u64),
    NoMatchingVariantsError,
    HstEndError,
    HstBothEndError,
    HstRightEndError,
    HstLeftEndError,
}

impl std::fmt::Display for HaptkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CoordsParseError(coords) => write!(f, "Failed to parse coords: {coords:?}"),
            Self::PosParseError((coords, value)) => write!(
                f,
                "Position {value:?} is not an integer in coords {coords:?}"
            ),
            Self::NormalizeError(pos) => write!(
                f,
                "At pos {pos} allele count != 2. Normalize alleles using bcftools norm"
            ),
            Self::PloidyError((pos, num)) => write!(
                f,
                "Ploidy error: A sample does not have the amount of genotypes required by ploidy at position: {pos}. Remember: missing genotypes are not allowed. {num}"
            ),
            Self::MissingGenotypesError(pos) => write!(f, "Missing genotypes at position: {pos}"),
            Self::SamplesNotFoundError => {
                write!(f, "None of the wanted samples was found in the vcf.")
            }
            Self::NoMatchingVariantsError => {
                write!(
                    f,
                    "None of the variants in the case file are present in the controls file."
                )
            }
            Self::HstEndError => {
                write!(
                    f,
                    "HST construction ended prematurely due to windowed file-read"
                )
            }
            Self::HstBothEndError => {
                write!(
                    f,
                    "HST construction ended prematurely due to windowed file-read"
                )
            }
            Self::HstLeftEndError => {
                write!(
                    f,
                    "HST construction ended prematurely due to windowed file-read"
                )
            }
            Self::HstRightEndError => {
                write!(
                    f,
                    "HST construction ended prematurely due to windowed file-read"
                )
            }
            Self::VariantPosNotFoundError((variant_pos, records_n)) => {
                write!(f, "The given coordinate {variant_pos} is larger than the largest found position. Total records read: {records_n}. Check your coordinates and vcf file. Comparing to a haplotype file also automatically limits min and max coordinates to the haplotype coordinates.")
            }
        }
    }
}
