/// Identify whether hapotype is present in samples
pub mod check_for_haplotype;

/// Check for vcf coverage
pub mod coverage;

/// Find the most recent common ancestor. The original R algorithm by Gandolfo et al translated to Rust. <https://github.com/bahlolab/DatingRareMutations>
pub mod mrca;

/// Shortcut to read vcf sample namples
pub mod read_sample_names;

/// Print out haplotype files for samples
pub mod read_sample_haplotype;

/// Bidirectional haplotype sharing tree algorithm
pub mod bhst;

/// Unidirectional haplotype sharing tree algorithm
pub mod uhst;

/// Haplotype to VCF
pub mod haplotype_to_vcf;

/// Identify differences to a haplotype
pub mod compare_to_haplotype;

/// Compare samples to a HST
pub mod compare_to_hst;

/// Compare haplotypes
pub mod compare_haplotypes;

// / Annotate haplotype
// pub mod annotate_haplotype;
