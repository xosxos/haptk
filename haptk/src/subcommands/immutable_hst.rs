use color_eyre::Result;
use petgraph::Graph;

use crate::{
    libs::structs::PhasedMatrix,
    structs::Coord,
    subcommands::{
        bhst::{initiate_hst, insert_nodes_to_bhst, Node},
        uhst::{insert_nodes_to_uhst, LocDirection},
    },
};

#[doc(hidden)]
pub fn construct_bhst_no_mut(
    vcf: &PhasedMatrix,
    start_coord: &Coord,
    min_size: usize,
    start_indexes: Option<Vec<usize>>,
) -> Result<Graph<Node, ()>> {
    let mut hst = initiate_hst(vcf, start_coord, start_indexes);
    let mut blacklist_nodes = vec![];

    loop {
        let blacklist = insert_nodes_to_bhst(vcf, &mut hst, &blacklist_nodes, min_size);

        match blacklist {
            Err(_) => {
                unreachable!()
            }
            Ok(blacklist) => {
                if blacklist.is_empty() {
                    break;
                }
                blacklist_nodes.extend(blacklist.iter().flatten());
            }
        }
    }

    Ok(hst)
}

#[doc(hidden)]
pub fn construct_uhst_no_mut(
    vcf: &PhasedMatrix,
    direction: &LocDirection,
    start_coord: &Coord,
    min_size: usize,
    only_majority: bool,
) -> Result<Graph<Node, ()>> {
    let mut hst = initiate_hst(vcf, start_coord, None);
    let mut blacklist_nodes = vec![];

    loop {
        let blacklist = insert_nodes_to_uhst(
            vcf,
            &mut hst,
            &blacklist_nodes,
            min_size,
            direction,
            only_majority,
            start_coord,
        );

        match blacklist {
            Err(_) => unreachable!(),
            Ok(blacklist) => {
                if blacklist.is_empty() {
                    break;
                }
                if only_majority && blacklist[0].is_empty() {
                    break;
                }
                blacklist_nodes.extend(blacklist.iter().flatten());
            }
        }
    }
    Ok(hst)
}
