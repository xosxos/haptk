from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip

from haptk import utils
from haptk.circle_tree import draw_circle_tree
from haptk.index_tree import draw_index_tree
from haptk.normal_tree import draw_normal_tree
from haptk.match_tree import draw_match_tree
from haptk.leaf_neighbors import find_leaf_neighbors

class HST:
    def __init__(self, G, coords, samples, metadata):
        self.G = G
        self.coords = coords
        self.samples = samples
        self.metadata = metadata

def read_hst(path):
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

        hst = HST(G, coords, samples, metadata)

        return hst

def create_ete3_tree(hst, min_size, hard_cut):
    graph_dict = utils.hst_to_graph_dict(hst.G)
    newick = utils.newickify(hst.G, graph_dict, min_size, hard_cut)
    return Tree(newick, format=1)

def leaf_neighbors(hst): 
    return find_leaf_neighbors(hst.G)

def match_tree(match_hst, other_hst, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, proportions=False): 
    other_t = create_ete3_tree(other_hst, min_size, hard_cut)
    draw_match_tree(match_hst.G, other_t, output, match_hst.samples, to_tag, colors, proportions)

def index_tree(hst, output, min_size=1, hard_cut=False): 
    t = create_ete3_tree(hst, min_size, hard_cut)
    draw_index_tree(t, output)

def circle_tree(hst, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False): 
    if len(to_tag) > len(colors):
        raise ValueError("more samples to tag than available colors")

    t = create_ete3_tree(hst, min_size, hard_cut)

    draw_circle_tree(hst.G, t, output, hst.samples, to_tag, colors)

def normal_tree(hst, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, proportions=False): 
    if len(to_tag) > len(colors):
        raise ValueError("more samples to tag than available colors")

    ete3_tree = create_ete3_tree(hst, min_size, hard_cut)

    draw_normal_tree(hst.G, ete3_tree, output, hst.samples, to_tag, colors, proportions)


