from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip

def read_path_to_hst(path):
    with gzip.open(path, 'rt') as f:
        json_bytes = f.read()
        hst_struct = json.loads(json_bytes)
        hst = hst_struct["hst"]
        coords = hst_struct["coords"]
        samples = hst_struct["samples"]
        metadata = hst_struct["metadata"]

        G = rx.PyDiGraph()

        edges_tuple_vec = []
        for edges in hst["edges"]:
            edges_tuple_vec.append((edges[0], edges[1], edges[2]))

        G.add_nodes_from(hst["nodes"])
        G.add_edges_from(edges_tuple_vec)

        if metadata["selection"] == "All":
            samples = [val for val in samples for _ in (0, 1)]

        return G, coords, samples, metadata

def hst_to_graph_dict(hst):
    graph_dict = {}
    for node_idx in hst.node_indices():
        neighbors = hst.neighbors(node_idx)
        if not len(neighbors) == 0:
            children_dict = {}

            for nbr_idx in neighbors:
                children_dict[nbr_idx] = hst.get_node_data(nbr_idx)

            graph_dict[node_idx] = children_dict
    return graph_dict

# Originally from https://stackoverflow.com/questions/50003007/how-to-convert-python-dictionary-to-newick-form-format
def newickify(hst, graph_dict, min_samples, hard_cut) -> str:
    visited_nodes = set()

    # Recursion
    def newick_render_node(node_idx, node_data) -> str:
        assert node_idx not in visited_nodes, "Error: The tree may not be circular!"

        if node_idx not in graph_dict or len(node_data["indexes"]) < min_samples:
            # Leaf node string
            return F'{node_idx}:0'
        else:
            # Add node to visited nodes
            visited_nodes.add(node_idx)

            # Get children
            children = graph_dict[node_idx]

            # Sort keys by size
            keys_and_size = [(child, len(hst.get_node_data(child)["indexes"])) for child in children.keys()]
            keys_and_size.sort(key = lambda x: x[1], reverse=False)

            if hard_cut:
                keys_and_size = list(filter(lambda x: x[1] >= min_samples, keys_and_size))

                if len(keys_and_size) == 0:
                    return F'{node_idx}:0'

            # Children strings
            children_strings = [newick_render_node(child, children[child]) for (child, _) in keys_and_size]
            children_strings = ",".join(children_strings)

            # Parent node string
            return F'({children_strings}){node_idx}:0'

    # Root node is 0
    root_node = 0
    node_data = hst.get_node_data(root_node)
    newick_string = newick_render_node(root_node, node_data) + ';'

    # assert visited_nodes == set(graph_dict.keys()), "Error: some nodes aren't in the tree"

    return newick_string

