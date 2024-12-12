from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip
import matplotlib as mpl
import numpy as np

from bisect import bisect_left, bisect_right

def find_gt(a, x):
    'Find leftmost value greater than x'
    i = bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_eq(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError


def color_fader(c1,c2,mix): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def return_node_style(color, size):
    style = NodeStyle()
    style["fgcolor"] = color
    style["size"] = size
    style["hz_line_type"] = 0
    style["vt_line_type"] = 0
    style["shape"] = "square"
    return style


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
def newickify(hst, graph_dict, min_samples, hard_cut, min_start, max_stop, right_up) -> str:
    visited_nodes = set()

    # Recursion
    def newick_render_node(node_idx, node_data) -> str:
        assert node_idx not in visited_nodes, "Error: The tree may not be circular!"

        start_stop_clause = False

        if min_start or max_stop:
            # parents = hst.G.predecessors(node_idx)
            # if len(parents) > 0:
                # parent_node = parents[0]
                # start = hst.coords[parent_node['start_idx']]['pos']
                # stop = hst.coords[parent_node['stop_idx']]['pos']
                # print(start, stop, len(parent_node['haplotype']), min_start, max_stop)
            start = node_data['start']['pos']
            stop = node_data['stop']['pos']

            # print(start, stop, len(node_data['haplotype']), min_start, max_stop)

            if min_start and start <= min_start:
                start_stop_clause = True

            if max_stop and stop >= max_stop:
                start_stop_clause = True              

        if node_idx not in graph_dict or len(node_data["indexes"]) < min_samples or start_stop_clause:
            # Leaf node string
            return F'{node_idx}:0'
        else:
            # Add node to visited nodes
            visited_nodes.add(node_idx)

            # Get children
            children = graph_dict[node_idx]

            # Sort keys by size
            keys_and_size = [(child, len(hst.G.get_node_data(child)["indexes"])) for child in children.keys()]
            if right_up:
                keys_and_size.sort(key = lambda x: x[1], reverse=True)
            else:
                keys_and_size.sort(key = lambda x: x[1], reverse=False)

            if hard_cut:
                keys_and_size = list(filter(lambda x: x[1] >= min_samples, keys_and_size))

                if len(keys_and_size) == 0:
                    return F'{node_idx}:0'

            # Children strings recursion
            children_strings = [newick_render_node(child, children[child]) for (child, _) in keys_and_size]
            children_strings = ",".join(children_strings)

            # Parent node string
            return F'({children_strings}){node_idx}:0'

    # Root node is 0
    root_node = 0
    node_data = hst.G.get_node_data(root_node)
    newick_string = newick_render_node(root_node, node_data) + ';'

    # assert visited_nodes == set(graph_dict.keys()), "Error: some nodes aren't in the tree"

    return newick_string

