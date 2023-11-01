from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip

import hst_utils

def draw_tree(hst, t, output, samples, samples_to_tag):
    ts = TreeStyle()
    style = NodeStyle()
    style["size"] = 0

    nnode_color = "#669e00"
    narea_color = "#DAF7A6"

    lnode_color = "#9F2B68"
    larea_color = "#FFB6C1"

    lstyle = NodeStyle()
    lstyle["fgcolor"] = lnode_color
    lstyle["shape"] = "square"
    lstyle["size"] = 5
    lstyle["hz_line_type"] = 0
    lstyle["vt_line_type"] = 0

    nstyle = NodeStyle()
    nstyle["fgcolor"] = nnode_color
    nstyle["size"] = 5
    nstyle["hz_line_type"] = 0
    nstyle["vt_line_type"] = 0
    nstyle["shape"] = "square"

    for n in t.traverse():
        n.set_style(style)

        node_data = hst.get_node_data(int(n.name))
        label = len(node_data["indexes"])

        if not n.is_leaf():
            F = TextFace(label, tight_text=True, penwidth=30, fsize=8)
            F.rotation = 270
            n.add_face(F, column=0)

        if n.is_leaf():
            results = []
            if len(samples_to_tag) > 0:
                results = [i for i in node_data["indexes"] if samples[i] in samples_to_tag]

            if results:
                n.set_style(nstyle)
                n.img_style["bgcolor"] = narea_color
            else:
                n.set_style(lstyle)
                n.img_style["bgcolor"] = larea_color

    ts.branch_vertical_margin = 10
    ts.show_leaf_name = False
    ts.show_scale = False

    ts.mode = "c"
    ts.arc_start = 0 # 0 degrees = 15:00 o'clock
    ts.arc_span = 350
    ts.min_leaf_separation = 0
    ts.root_opening_factor = 1
    ts.optimal_scale_level = "full"
    ts.allow_face_overlap = True

    t.render(output, w=2024, h=2024, units="px", tree_style=ts, dpi=600)


### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--ids', type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()


# Read .hst.gz and wrangle it to an ete3 tree
HST, coords, samples = hst_utils.read_path_to_hst(args.hst)
graph_dict = hst_utils.hst_to_graph_dict(HST)
newick = hst_utils.newickify(HST, graph_dict, args.min_size, args.hard_cut)

ete3_tree = Tree(newick, format=1)

# Read in potential sample names to tag
samples_to_tag = []
if args.ids:
    file = open(args.ids, 'r')
    for line in file.readlines():
        samples_to_tag.append(line.strip())

if args.output:
    draw_tree(HST, ete3_tree, args.output, samples, samples_to_tag)
else:
    draw_tree(HST, ete3_tree,  "mytree.png", samples, samples_to_tag)



