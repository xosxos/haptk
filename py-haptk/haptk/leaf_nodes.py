import itertools

def find_leaf_neighbors(hst):
    # A set of sibling nodes from the tree
    neighbor_sets = []

    for idx in hst.node_indices():
        # A node is leaf if it has no outgoing edges
        if len(hst.out_edges(idx)) == 0:

            # Checks if all siblings have exactly 1 sample
            is_paired_leaf = True

            # Add the first sibling
            neighbors = [idx]

            # Get parents by moving up to adj_direction
            parents = hst.adj_direction(idx, True)

            # Iterate parent for siblings (I think this loops only once in HSTs)
            for parent in parents:
                for sibling in hst.adj_direction(parent, False):
                    # Dont re-add the first sibling
                    if sibling != idx:
                        data = hst.get_node_data(sibling)

                        # Dont add if there is more than 1 sample in the node
                        if len(data['indexes']) == 1:
                            neighbors.append(sibling)
                            pass
                        else:
                            # A sibling has more than 1 sample
                            is_paired_leaf = False

            if is_paired_leaf:
                neighbors.sort()
                neighbor_sets.append(neighbors)

                
    # Sort
    neighbor_sets.sort()

    # Dedup because looping all leaf nodes produces duplicate pairs
    neighbor_sets = list(k for k,_ in itertools.groupby(neighbor_sets))

    return neighbor_sets

   
