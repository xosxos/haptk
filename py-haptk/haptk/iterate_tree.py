from ete3 import TreeStyle, TextFace, CircleFace,  NodeStyle, add_face_to_node
from scipy.stats import fisher_exact

from haptk import utils

def iterate_tree_style():
    ts = TreeStyle()

    ts.branch_vertical_margin = 10
    ts.show_leaf_name = False
    ts.show_scale = False

    ts.mode = "c"
    ts.arc_start = 0 # 0 degrees = 15:00 o'clock
    ts.arc_span = 355
    ts.min_leaf_separation = 0
    ts.root_opening_factor = 0
    ts.optimal_scale_level = "mid"
    ts.allow_face_overlap = False
    # ts.allow_face_overlap = False

    return ts

def iterate_tree_inner(hst, t, samples_to_tag, df, optimizer):
    style = NodeStyle()
    style["size"] = 0

    for n in t.traverse():
        n.set_style(style)
        
        F = optimizer(hst, n, df, samples_to_tag)      

        n.add_face(F, column=0)
        

    return t

