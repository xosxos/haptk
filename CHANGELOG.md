# CHANGELOG

## haptk version 0.3.2
- Regression: Bhst and uhst got stuck in a loop if genotyping data ran out on both sides
- Python: Recombination rates not mandatory for ancestral segment lengths python plots
- Python: Fixed HST construction for --all, published new version py-haptk 0.0.8

## haptk version 0.3.1

Features:
- Added sharded VCF reading to mrca
- Refactored a much faster haplotype comparison
- Faster longest-haplotype reading for compare-to-haplotype and compare-to-hst.
- Added sharded longest-haplotype selection to list-haplotypes
- Compare-to-haplotype speed up optimizations and sharded only-longest reading
- Added a circular visualization for match HSTs. Also, added a longest-leafs-only parameter to compare to HSTs
- Support reading ids from a hst scan
- --no-alt paramater now diregards all non-contradictory genotypes to save memory and compute

Enhancements:
- Changed tracing::error into a warning when duplicate genotypes are present
- Better error message for mrca gamma distribution failure
- Better logging

Fixes:
- Fixed marker listing from HSTs
- Fixed a hst comparison regression
- Fixed a regression on longest-haplotype selection for comparing to haplotypes


## haptk version 0.3
- HST construction by reading VCFs in parts instead of the whole file
- The --window parameter defines the starting basepair window when constructing HSTs
- Node start_idx and stop_idx changed to coordinates: `Coord`
- Bidirectional HST now always constructs the first nodes as biallelic if the starting variant is contradictory
- Haplotype comparison tools updated to support multihaplotype files and calculate frequencies
- Fasta2haplotype conversion added
- `Coverage` subcommand deprecated
- Rust based graphs deprecated except for `MatrixGraph`
- Experimental genome-wide methods outlined

