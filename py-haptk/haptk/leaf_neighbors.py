import itertools

def find_leaf_neighbors(hst):
    neighbor_sets = []
    for idx in hst.node_indices():
        if len(hst.out_edges(idx)) == 0:

            is_paired_leaf = True
            neighbors = [idx]

            parent = hst.adj_direction(idx, True)

            for k in parent:
                for neighbor in hst.adj_direction(k, False):
                    if neighbor not in neighbors:
                        data = hst.get_node_data(neighbor)
                        # not a leaf node with only 1 sample neighbors
                        if len(data['indexes']) == 1:
                            neighbors.append(neighbor)
                            pass
                        else:
                            is_paired_leaf = False

            if is_paired_leaf:
                neighbor_sets.append(neighbors)

            neighbors.sort()
                
    neighbor_sets.sort()
    neighbor_sets = list(k for k,_ in itertools.groupby(neighbor_sets))
    return neighbor_sets
   
