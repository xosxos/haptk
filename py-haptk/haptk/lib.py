from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip
import math

from haptk import utils
from haptk.circle_tree import draw_circle_tree, circle_tree_style
from haptk.index_tree import draw_index_tree, index_tree_style
from haptk.normal_tree import draw_normal_tree, normal_tree_style
from haptk.match_tree import circle_match_tree_style, draw_match_tree, match_tree_style
from haptk.leaf_neighbors import find_leaf_neighbors
from haptk.iterate_tree import iterate_tree_inner, iterate_tree_style

class HST:
    def __init__(self, G, coords, samples, metadata):
        self.G = G
        self.coords = coords
        self.samples = samples
        self.metadata = metadata

    def get_node_data(self, idx):
        idx = int(idx)
        return self.G.get_node_data(idx)

    def get_node_haplotype(self, idx):
        idx = int(idx)
        data = self.G.get_node_data(idx)
        return data["haplotype"]

    def get_node_indexes(self, idx):
        idx = int(idx)
        data = self.G.get_node_data(idx)
        return data["indexes"]

    def get_node_start(self, idx):
        idx = int(idx)
        data = self.G.get_node_data(idx)
        return data["start"]

    def get_node_stop(self, idx):
        idx = int(idx)
        data = self.G.get_node_data(idx)
        return data["stop"]
    
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
                node_samples.append(self.get_sample_name(index))

            list_of_samples.append(node_samples)
        return list_of_samples

    def get_sample_name(self, index):
        # This division always results in whole numbers on the Rust side
        return self.samples[int(index / self.ploidy_int())]

    def ploidy_int(self):
        match self.metadata['ploidy']:
            case "Mixed":
                return 1
            case "Haploid":
                return 1
            case "Diploid":
                return 2
        
    
    def get_nodes_data(self, nodes):
        data = []

        for node in nodes:
            data.append(self.G.get_node_data(node))

        return data

    def haplotype_string(self, node_data, size):
        start = Coord.new(node_data['start'])
        start_idx = utils.find_eq(self.coords, start)

        haplotype = [self.coords[start_idx + i].get_genotype(gt) for i, gt in enumerate(node_data['haplotype'])]

        if len(haplotype) > size:
            between_len = len(haplotype) - size - 1
            if between_len < 2:
                ht_string = '-'.join(str(x) for x in haplotype)

                return f'{ht_string}'
            
            if self.metadata['hst_type'] == 'UhstLeft':
                ending = haplotype[-size:]
                ending = '-'.join(str(x) for x in ending)

                return f'{haplotype[0]}..({between_len})..{ending}'

            elif self.metadata['hst_type'] == 'UhstRight':
                beginning = haplotype[0:size]
                beginning = '-'.join(str(x) for x in beginning)

                return f'{beginning}..({between_len})..{haplotype[-1]}'

            elif self.metadata['hst_type'] == 'Bhst':
                left = int(math.floor(size/2))
                right = int(math.floor(size/2))
                between_len = len(haplotype) - left - right

                beginning = haplotype[0:left]
                ending = haplotype[-right:]
                beginning = '-'.join(str(x) for x in beginning)
                ending = '-'.join(str(x) for x in ending)

                return f'{beginning}..({between_len})..{ending}'

            else:
                print(f"hst_type {self.metadata['hst_type']} is not supported")
        else:
            haplotype_string = '-'.join(str(x) for x in haplotype)

            return f'{haplotype_string}'


    def is_maj_up_to(self, node, value):
        iter_node = node
        # print(f"checking {iter_node.name} for majority")
        for _ in range(0, value):
            if iter_node.name == "0":
                return True

            if not iter_node.up.name == "0":
                # print("input node:", iter_node.up)
                is_big = self.is_majority_node(iter_node.up)
                if is_big == False:
                    # print(f"run {_}: {iter_node.up.name} is not big")
                    return False

            iter_node = iter_node.up
            
        return True

    def is_majority_node(self, node):
        parent = node.up
        input_node_data = self.get_node_data(int(node.name))
        value = len(input_node_data["indexes"])
        # print(f"input node {node.name} has {value} indexes")

        siblings = {}
        for child in parent.children:
            data = self.get_node_data(int(child.name))
            if not child.name == node.name:
                siblings[child.name] = len(data["indexes"])


        # print(value, siblings.values())
        for v in siblings.values():
            if value < v: 
                return False
  

    def get_children_idx_amounts(self, node):
        return [len(self.get_node_data(c.name)["indexes"]) for c in node.children]

    def iterate_tree(self, df, optimizer, output, to_tag=[], min_size=1, hard_cut=False, min_start=None, max_stop=None, right_up=False, tree_style = iterate_tree_style(), dpi=600, w=1624, h=1624):
        t = create_ete3_tree(self, min_size, hard_cut, min_start, max_stop, right_up)
        t = iterate_tree_inner(self, t, to_tag, df, optimizer)
        t.render(output, units="px", w=w, h=h, tree_style=tree_style, dpi=dpi)

    def match_tree(self, other_hst, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, min_start=None, max_stop=None, proportions=False, dpi=600, tree_style=match_tree_style(), circle=False): 
        other_t = create_ete3_tree(other_hst, min_size, hard_cut, min_start, max_stop, right_up=False)
        t = draw_match_tree(self, other_t, output, self.samples, to_tag, colors, proportions, circle)

        if circle:
            tree_style=circle_match_tree_style()

        t.render(output, units="px", tree_style=tree_style, dpi=dpi)

    def index_tree(self, output, min_size=1, hard_cut=False, min_start=None, max_stop=None, dpi=600, tree_style=index_tree_style()): 
        t = create_ete3_tree(self, min_size, hard_cut, min_start, max_stop, right_up=False)
        t = draw_index_tree(t, output)

        t.render(output, units="px", tree_style=tree_style, dpi=dpi)

    def circle_tree(self, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, min_start=None, max_stop=None, w=1624, h=1624, dpi=600, tree_style=circle_tree_style(), branch_length_as_majority = 999999999, branch_point_size = 999999999): 
        if len(to_tag) > len(colors):
            raise ValueError("more samples to tag than available colors")

        t = create_ete3_tree(self, min_size, hard_cut, min_start, max_stop, right_up=False)
        t = draw_circle_tree(self, t, to_tag, colors, branch_point_size, branch_length_as_majority)

        t.render(output, w=w, h=h, units="px", tree_style=tree_style, dpi=dpi)

    def normal_tree(self, output, to_tag=[], colors = ["#ff0000", "#FF69B4", "#4cfe92", "#4ccbfe", "#c9efff", "orange", "yellow"], min_size=1, hard_cut=False, min_start=None, max_stop=None, proportions=False, dpi=600, tree_style=normal_tree_style(), show_haplotype=False, n_markers=3, show_pos=False, branch_length_as_majority = 999999999, branch_point_size = 999999999): 
        if len(to_tag) > len(colors):
            raise ValueError("more samples to tag than available colors")

        t = create_ete3_tree(self, min_size, hard_cut, min_start, max_stop, right_up=False)
        t = draw_normal_tree(self, t, to_tag, colors, proportions, show_haplotype, n_markers, show_pos, branch_point_size, branch_length_as_majority)

        t.render(output, units="px", tree_style=tree_style, dpi=dpi)

class Coord:
    def __init__(self, contig, pos, ref, alt):
        self.contig = contig
        self.ref = ref
        self.alt = alt
        self.pos = pos

    def new(dct):
        return Coord(dct['contig'], dct['pos'], dct['ref'], dct['alt'])

    def get_genotype(self, gt):
        if gt == 0:
            return self.ref
        if gt == 1:
            return self.alt
        else:
            raise ValueError

    def __repr__(self):
        return f'{self.contig}_{self.pos}_{self.ref}_{self.alt}'

    def __lt__(self, other):
        # return self.pos < other.pos
        if self.pos < other.pos:
            return True

        if self.pos == other.pos:
            if self.ref < other.ref:
                return True
            elif self.ref == other.ref:
                return self.alt < other.alt

        return False


    def __gt__(self, other):
        # return self.pos > other.pos
        if self.pos > other.pos:
            return True

        if self.pos == other.pos:
            if self.ref > other.ref:
                return True
            elif self.ref == other.ref:
                return self.alt > other.alt

        return False

    def __eq__(self, other):
        return self.pos == other.pos and self.ref == other.ref and self.alt == other.alt
      
def read_hst(path):
    with gzip.open(path, 'rt') as f:
        json_bytes = f.read()
        hst_struct = json.loads(json_bytes)
        hst = hst_struct["hst"]
        metadata = hst_struct["metadata"]
        coords = [Coord(dct['contig'], dct['pos'], dct['ref'], dct['alt']) for dct in metadata["coords"]]
        samples = metadata["samples"]
        del metadata["coords"]
        del metadata["samples"]

        G = rx.PyDiGraph()

        edges_tuple_vec = []
        for edges in hst["edges"]:
            edges_tuple_vec.append((edges[0], edges[1], edges[2]))

        G.add_nodes_from(hst["nodes"])
        G.add_edges_from(edges_tuple_vec)

        # if metadata["selection"] == "All":
            # samples = [val for val in samples for _ in (0, 1)]

        hst = HST(G, coords, samples, metadata)

        return hst

def create_ete3_tree(hst, min_size, hard_cut, min_start, max_stop, right_up):
    graph_dict = utils.hst_to_graph_dict(hst.G)
    newick = utils.newickify(hst, graph_dict, min_size, hard_cut, min_start, max_stop, right_up)
    return Tree(newick, format=1)



        
    


