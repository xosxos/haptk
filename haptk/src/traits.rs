pub use helper::MrcaHelper;

pub mod helper {
    use color_eyre::Result;
    use petgraph::graph::NodeIndex;
    use petgraph::Graph;

    use crate::core::Coord;
    use crate::core::PhasedMatrix;
    use crate::subcommands::hst::find_majority_nodes;
    use crate::subcommands::hst::initiate_hst;
    use crate::subcommands::hst::Node;
    use crate::subcommands::immutable_hst;
    use crate::subcommands::uhst;
    use crate::subcommands::uhst::LocDirection;

    pub trait MrcaHelper {
        fn get_lengths_from_uhst_no_mut(&self, start_coord: &Coord) -> Result<Vec<(Node, Node)>>;
        fn get_lengths_from_uhst(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>>;
        fn get_lengths(
            &self,
            uhst_left: Graph<Node, ()>,
            uhst_right: Graph<Node, ()>,
        ) -> Vec<(Node, Node)>;
    }

    impl MrcaHelper for PhasedMatrix {
        fn get_lengths_from_uhst(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
            let mut uhst_left = initiate_hst(self, start_coord, None);
            let mut uhst_right = initiate_hst(self, start_coord, None);

            uhst::populate_uhst(
                self,
                &mut uhst_left,
                &LocDirection::Left,
                start_coord,
                1,
                true,
            )?;
            uhst::populate_uhst(
                self,
                &mut uhst_right,
                &LocDirection::Right,
                start_coord,
                1,
                true,
            )?;

            Ok(self.get_lengths(uhst_left, uhst_right))
        }

        fn get_lengths_from_uhst_no_mut(&self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
            let uhst_right = immutable_hst::construct_uhst_no_mut(
                self,
                &LocDirection::Right,
                start_coord,
                1,
                true,
            )?;
            let uhst_left = immutable_hst::construct_uhst_no_mut(
                self,
                &LocDirection::Left,
                start_coord,
                1,
                true,
            )?;

            Ok(self.get_lengths(uhst_left, uhst_right))
        }

        fn get_lengths(
            &self,
            uhst_left: Graph<Node, ()>,
            uhst_right: Graph<Node, ()>,
        ) -> Vec<(Node, Node)> {
            let start_idx = NodeIndex::new(0);
            let lmaj_branch = find_majority_nodes(&uhst_left, start_idx);
            let rmaj_branch = find_majority_nodes(&uhst_right, start_idx);

            (0..self.nhaplotypes())
                .map(|idx| {
                    let mut lnode = lmaj_branch.len()
                        - lmaj_branch
                            .iter()
                            .rev()
                            .position(|(node, _)| node.indexes.contains(&idx))
                            .unwrap();
                    let mut rnode = rmaj_branch.len()
                        - rmaj_branch
                            .iter()
                            .rev()
                            .position(|(node, _)| node.indexes.contains(&idx))
                            .unwrap();
                    if rnode == rmaj_branch.len() {
                        rnode -= 1;
                    }

                    if lnode == lmaj_branch.len() {
                        lnode -= 1;
                    }

                    (
                        idx,
                        lmaj_branch.get(lnode).unwrap(),
                        rmaj_branch.get(rnode).unwrap(),
                    )
                })
                .map(|(_idx, (lnode, _), (rnode, _))| ((*lnode).clone(), (*rnode).clone()))
                .collect()
        }
    }
}

pub use only_longest::OnlyLongest;
mod only_longest {
    use color_eyre::Result;

    use crate::core::Coord;
    use crate::core::PhasedMatrix;
    use crate::core::Ploidy;
    use crate::subcommands::hst::Node;

    use crate::traits::MrcaHelper;

    pub trait OnlyLongest {
        fn select_only_longest(&mut self) -> Result<()>;
        fn select_only_longest_no_shard(&mut self) -> Result<()>;
        fn only_longest_indexes(&mut self) -> Result<Vec<usize>>;
        fn only_longest_indexes_no_shard(&self, start: &Coord) -> Result<Vec<usize>>;
        fn get_only_longest_lookups(&mut self) -> Result<Vec<[bool; 2]>>;
        fn get_only_longest_lookups_no_shard(&self) -> Result<Vec<[bool; 2]>>;
        fn only_longest_lengths(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>>;
        fn only_longest_lengths_no_shard(&self, start_coord: &Coord) -> Result<Vec<(Node, Node)>>;
        fn lengths_to_indexes(&self, lengths: Vec<(Node, Node)>) -> Vec<usize>;
        fn nodes_to_lengths(&self, lengths: Vec<(Node, Node)>) -> Vec<(Node, Node)>;
    }

    impl OnlyLongest for PhasedMatrix {
        fn select_only_longest(&mut self) -> Result<()> {
            let longest_indexes = self.only_longest_indexes()?;

            // Update lookups
            let mut lookups = vec![];

            for idx in &longest_indexes {
                let idxs = self.get_idxs_for_samples(&[self.get_sample_name(*idx)])?;
                let pos = idxs.iter().position(|i| idx == i).unwrap();
                let lookup = match pos {
                    0 => [true, false],
                    1 => [false, true],
                    _ => unreachable!("Only diploid genotypes are supported"),
                };
                lookups.push(lookup);
            }

            self.metadata.lookups = lookups;

            self.select_rows(longest_indexes);
            self.ploidy = Ploidy::Haploid;

            tracing::info!("Finished only-longest selection.");
            Ok(())
        }

        fn select_only_longest_no_shard(&mut self) -> Result<()> {
            let longest_indexes = self.only_longest_indexes_no_shard(self.start_coord())?;

            self.select_rows(longest_indexes);
            self.ploidy = Ploidy::Haploid;

            tracing::info!("Finished only-longest selection.");
            Ok(())
        }

        fn only_longest_indexes(&mut self) -> Result<Vec<usize>> {
            let start = self.start_coord().clone();
            let lengths = self.get_lengths_from_uhst(&start)?;

            Ok(self.lengths_to_indexes(lengths))
        }

        fn get_only_longest_lookups(&mut self) -> Result<Vec<[bool; 2]>> {
            let indexes = self.only_longest_indexes()?;

            let mut lookups = vec![];

            for idx in indexes {
                let idxs = self.get_idxs_for_samples(&[self.get_sample_name(idx)])?;
                let pos = idxs.iter().position(|i| idx == *i).unwrap();
                let lookup = match pos {
                    0 => [true, false],
                    1 => [false, true],
                    _ => unreachable!("Only diploid genotypes are supported"),
                };
                lookups.push(lookup);
            }

            Ok(lookups)
        }

        fn get_only_longest_lookups_no_shard(&self) -> Result<Vec<[bool; 2]>> {
            let indexes = self.only_longest_indexes_no_shard(self.start_coord())?;

            let mut lookups = vec![];

            for idx in indexes {
                let idxs = self.get_idxs_for_samples(&[self.get_sample_name(idx)])?;
                let pos = idxs.iter().position(|i| idx == *i).unwrap();
                let lookup = match pos {
                    0 => [true, false],
                    1 => [false, true],
                    _ => unreachable!("Only diploid genotypes are supported"),
                };
                lookups.push(lookup);
            }

            Ok(lookups)
        }

        fn only_longest_indexes_no_shard(&self, start: &Coord) -> Result<Vec<usize>> {
            let lengths = self.get_lengths_from_uhst_no_mut(start)?;

            Ok(self.lengths_to_indexes(lengths))
        }

        fn lengths_to_indexes(&self, lengths: Vec<(Node, Node)>) -> Vec<usize> {
            let lengths = lengths
                .iter()
                .map(|(lnode, rnode)| {
                    if lnode.start == rnode.stop {
                        0
                    } else {
                        let stop = self.coords().range(..=&rnode.stop).next_back().unwrap();

                        let start = self.coords().range(&lnode.start..).nth(1).unwrap();
                        stop.pos - start.pos + 1
                    }
                })
                .collect::<Vec<u64>>();

            (0..self.nsamples())
                .map(|i| {
                    let lengths = lengths
                        .iter()
                        .enumerate()
                        .skip(i * *self.ploidy)
                        .take(*self.ploidy);

                    let (max_idx, max_len) = lengths.clone().max_by_key(|(_, l)| *l).unwrap();

                    if !self.is_genome_wide() && lengths.filter(|(_, l)| *l == max_len).count() > 1
                    {
                        tracing::warn!(
                            "Sample {} has two equally long haplotypes in only-longest selection.",
                            self.get_sample_name(i)
                        );
                    }

                    max_idx
                })
                .collect()
        }

        fn only_longest_lengths(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
            let lengths = self.get_lengths_from_uhst(start_coord)?;

            Ok(self.nodes_to_lengths(lengths))
        }

        fn only_longest_lengths_no_shard(&self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
            let lengths = self.get_lengths_from_uhst_no_mut(start_coord)?;

            Ok(self.nodes_to_lengths(lengths))
        }

        fn nodes_to_lengths(&self, lengths: Vec<(Node, Node)>) -> Vec<(Node, Node)> {
            let calculate_len = |(lnode, rnode): &(Node, Node)| {
                if lnode.start == rnode.stop {
                    0
                } else {
                    let stop = self.coords().range(..&rnode.stop).next_back().unwrap();
                    let start = self.coords().range(&lnode.start..).nth(1).unwrap();
                    stop.pos.saturating_sub(start.pos + 1)
                }
            };

            (0..self.nsamples())
                .map(|i| {
                    let lengths = lengths
                        .iter()
                        .enumerate()
                        .skip(i * *self.ploidy)
                        .take(*self.ploidy);

                    let (_, max_nodes) = lengths
                        .clone()
                        .max_by_key(|(_, l)| calculate_len(l))
                        .unwrap();

                    let max_len = calculate_len(max_nodes);

                    if !self.is_genome_wide()
                        && lengths.filter(|(_, l)| calculate_len(l) == max_len).count() > 1
                    {
                        tracing::warn!(
                            "Sample {} has two equally long haplotypes.",
                            self.get_sample_name(i)
                        );
                    }

                    max_nodes.clone()
                })
                .collect()
        }
    }
}
