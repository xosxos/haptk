from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip

import hst_utils

def draw_tree(hst, t, output):
    ts = TreeStyle()
    style = NodeStyle()
    style["size"] = 0

    lnode_color = "#9F2B68"
    larea_color = "#FFB6C1"

    lstyle = NodeStyle()
    lstyle["fgcolor"] = "#000"
    lstyle["shape"] = "square"
    lstyle["size"] = 1
    lstyle["hz_line_type"] = 0
    lstyle["vt_line_type"] = 0

    for n in t.traverse():
        n.set_style(style)

        label = n.name

        F = TextFace(label, tight_text=True, penwidth=30)
        F.rotation = 270
        n.add_face(F, column=0)
        n.set_style(lstyle)

    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.min_leaf_separation = 0
    ts.branch_vertical_margin = 5
    ts.optimal_scale_level = "full"
    ts.allow_face_overlap = True

    t.render(output, units="px", tree_style=ts, dpi=600)

### Parameters
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()


# Read .hst.gz and wrangle it to an ete3 tree
HST, coords, samples, metadata = hst_utils.read_path_to_hst(args.hst)

graph_dict = hst_utils.hst_to_graph_dict(HST)

newick = hst_utils.newickify(HST, graph_dict, args.min_size, args.hard_cut)

ete3_tree = Tree(newick, format=1)

if args.output:
    draw_tree(HST, ete3_tree, args.output)
else:
    draw_tree(HST, ete3_tree,  "mytree.png")



