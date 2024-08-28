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
    ts.allow_face_overlap = True
    # ts.allow_face_overlap = False

    return ts

def iterate_tree_inner(hst, t, samples_to_tag, df, optimizer):
    style = NodeStyle()
    style["size"] = 0

    for n in t.traverse():
        n.set_style(style)
        
        node_name = n.name
        indexes = hst.get_node_indexes(node_name)
        optimized_value, label, color, text_color = optimizer(hst, n, indexes, df, samples_to_tag)
        
        # if optimized_value < 0.0005:
        #     label = len(indexes)
        #     F = TextFace(f"{optimized_value:.4}", tight_text=True, penwidth=30, fsize=8, fgcolor='red')
        #     F.rotation = 270
        #     F.background.color = "#FFF"
        #     F.margin_right = 10
        #     n.add_face(F, column=0)

        F = TextFace(label, tight_text=True, penwidth=30, fsize=8, fgcolor="black")
        F.rotation = 270

        # mix = optimized_value * 100
        # if mix > 1:
        #     mix = 1
        # color = utils.color_fader('red', 'blue', mix)

        F.background.color = color
        F.fgcolor = text_color
        n.add_face(F, column=0)
        

    return t

