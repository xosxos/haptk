/// Identify whether hapotype is present in samples
pub mod check_for_haplotype;

/// Find the most recent common ancestor. The original R algorithm by Gandolfo et al translated to Rust. <https://github.com/bahlolab/DatingRareMutations>
pub mod mrca;

/// MRCA scan
pub mod mrca_scan;

/// Shortcut to read vcf sample namples
pub mod list_samples;

/// Output haplotype files for samples or HSTs
pub mod list_haplotypes;

/// Output marker lists for HSTs
pub mod list_markers;

/// Bidirectional haplotype sharing tree algorithm
pub mod bhst;

/// Unidirectional haplotype sharing tree algorithm
pub mod uhst;

/// Haplotype to VCF
pub mod haplotype_to_vcf;

/// Fasta to haplotype
pub mod fasta_to_haplotype;

/// Identify differences to a haplotype
pub mod compare_to_haplotype;

/// Compare samples to a HST
pub mod compare_to_hst;

/// Compare haplotypes
pub mod compare_haplotypes;

// / Annotate haplotype
// pub mod annotate_haplotype;

pub mod scan;
