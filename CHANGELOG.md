# CHANGELOG

haptk version 0.3
- HST construction by reading VCFs in parts instead of the whole file
- The --window parameter defines the starting basepair window when constructing HSTs
- Node start_idx and stop_idx changed to coordinates: `Coord`
- Bidirectional HST now always constructs the first nodes as biallelic if the starting variant is contradictory
- Haplotype comparison tools updated to support multihaplotype files and calculate frequencies
- Fasta2haplotype conversion added
- `Coverage` subcommand deprecated
- Rust based graphs deprecated except for `MatrixGraph`
- Experimental genome-wide methods outlined
