from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip

from haptk import utils
from haptk.circle_tree import draw_circle_tree, circle_tree_style
from haptk.index_tree import draw_index_tree, index_tree_style
from haptk.normal_tree import draw_normal_tree, normal_tree_style
from haptk.match_tree import draw_match_tree, match_tree_style
from haptk.leaf_neighbors import find_leaf_neighbors

class HST:
    def __init__(self, G, coords, samples, metadata):
        self.G = G
        self.coords = coords
        self.samples = samples
        self.metadata = metadata

    def get_node_data(self, idx):
        return self.G.get_node_data(idx)

    
    def leaf_neighbors(self): 
        return find_leaf_neighbors(self.G)

    def leaf_nodes(self):
        leaf_nodes = []

        for idx in self.G.node_indices():
            if len(self.G.out_edges(idx)) == 0:
                leaf_nodes.append(idx)

        return leaf_nodes

    def node_samples(self, nodes):
        list_of_samples = []

        for node in nodes:
            data = self.G.get_node_data(node)
            node_samples = []

            for index in data['indexes']:
                node_samples.append(self.samples[index])

            list_of_samples.append(node_samples)
        return list_of_samples

    def get_nodes_data(self, nodes):
        data = []

        for node in nodes:
            data.append(self.G.get_node_data(node))

        return data

    def match_tree(self, other_hst, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, proportions=False, dpi=600, tree_style=match_tree_style()): 
        other_t = create_ete3_tree(other_hst, min_size, hard_cut)
        t = draw_match_tree(self.G, other_t, output, self.samples, to_tag, colors, proportions)

        t.render(output, units="px", tree_style=tree_style, dpi=dpi)

    def index_tree(self, output, min_size=1, hard_cut=False, dpi=600, tree_style=index_tree_style()): 
        t = create_ete3_tree(self, min_size, hard_cut)
        t = draw_index_tree(t, output)

        t.render(output, units="px", tree_style=tree_style, dpi=dpi)

    def circle_tree(self, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, w=1624, h=1624, dpi=600, tree_style=circle_tree_style()): 
        if len(to_tag) > len(colors):
            raise ValueError("more samples to tag than available colors")

        t = create_ete3_tree(self, min_size, hard_cut)
        t = draw_circle_tree(self.G, t, output, self.samples, to_tag, colors)

        t.render(output, w=w, h=h, units="px", tree_style=tree_style, dpi=dpi)

    def normal_tree(self, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, proportions=False, dpi=600, tree_style=normal_tree_style()): 
        if len(to_tag) > len(colors):
            raise ValueError("more samples to tag than available colors")

        t = create_ete3_tree(self, min_size, hard_cut)
        draw_normal_tree(self.G, t, output, self.samples, to_tag, colors, proportions)

        t.render(output, units="px", tree_style=tree_style, dpi=dpi)

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



        
    


