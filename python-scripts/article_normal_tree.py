from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import PyQt5
import pygraphviz as pgv
import networkx as nx
import json
import argparse

# Originally from https://stackoverflow.com/questions/50003007/how-to-convert-python-dictionary-to-newick-form-format
def newickify(graph_dict, min_samples) -> str:
    visited_nodes = set()

    # Recursion
    def newick_render_node(node_idx, node_data) -> str:
        assert node_idx not in visited_nodes, "Error: The tree may not be circular!"

        if node_idx not in graph_dict or node_data["n"] < min_samples:
            # Leaf node string
            return F'{node_idx}:0'
        else:
            # Add node to visited nodes
            visited_nodes.add(node_idx)

            # Get children
            children = graph_dict[node_idx]

            # Sort keys by size
            keys_and_size = [(child, json.loads(G.nodes[child]["label"])["n"]) for child in children.keys()]
            keys_and_size.sort(key = lambda x: x[1], reverse=False)

            # Children strings
            children_strings = [newick_render_node(child, children[child]) for (child, _) in keys_and_size]
            children_strings = ",".join(children_strings)

            # Parent node string
            return F'({children_strings}){node_idx}:0'

    root_node = '0'
    node_label = G.nodes[root_node]["label"]
    node_data = json.loads(node_label)
    newick_string = newick_render_node(root_node, node_data) + ';'

    # assert visited_nodes == set(graph_dict.keys()), "Error: some nodes aren't in the tree"

    return newick_string

def draw_tree(t, output, samples_to_tag):
    ts = TreeStyle()
    style = NodeStyle()
    style["size"] = 0

    nnode_color = "#669e00"
    narea_color = "#DAF7A6"

    lnode_color = "#9F2B68"
    larea_color = "#FFB6C1"

    lstyle = NodeStyle()
    lstyle["fgcolor"] = "#000"
    lstyle["shape"] = "square"
    lstyle["size"] = 1
    lstyle["hz_line_type"] = 0
    lstyle["vt_line_type"] = 0

    nstyle = NodeStyle()
    nstyle["fgcolor"] = nnode_color
    nstyle["size"] = 1
    nstyle["hz_line_type"] = 0
    nstyle["vt_line_type"] = 0
    nstyle["shape"] = "square"

    for n in t.traverse():
        n.set_style(style)
        node_label = G.nodes[n.name]["label"]
        node_data = json.loads(node_label)

        F = TextFace(node_data["n"], tight_text=True, penwidth=30)
        # F = TextFace(n.name, tight_text=True)
        F.rotation = 270
        n.add_face(F, column=0)

        result = [i for i in samples_to_tag if i in node_data["samples"]]

        if n.is_leaf():
            if result:
                n.set_style(nstyle)
            else:
                n.set_style(lstyle)

    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False
    # ts.orientation = 1

    # ts.mode = "c"
    # ts.arc_start = 0 # 0 degrees = 15:00 o'clock
    # ts.arc_span = 345
    # ts.root_opening_factor = 1
    ts.min_leaf_separation = 0
    ts.branch_vertical_margin = 5
    ts.optimal_scale_level = "full"
    ts.allow_face_overlap = True

    t.render(output, units="px", tree_style=ts, dpi=600)

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('file', type=str)
parser.add_argument('--min-size', type=int, default=1)
parser.add_argument('--ids', type=str)
parser.add_argument('--output', type=str)
args = parser.parse_args()

G = nx.DiGraph(pgv.AGraph(args.file, strict=True, directed=True))

graph_dict = {}
for node_idx, neighbors_dict in G.adjacency():
    if not len(neighbors_dict) == 0:

        children_dict = {}
        for nbr_idx in neighbors_dict.keys():
            node_label = G.nodes[nbr_idx]["label"]
            children_dict[nbr_idx] = json.loads(node_label)

        graph_dict[node_idx] = children_dict

# print(graph_dict)


newick = newickify(graph_dict, args.min_size)
t = Tree(newick, format=1)

samples_to_tag = []
if args.ids:
    file = open(args.ids, 'r')
    for line in file.readlines():
        samples_to_tag.append(line.strip())
    
if args.output:
    draw_tree(t, args.output, samples_to_tag)
else:
    draw_tree(t, "mytree.png", samples_to_tag)
    

