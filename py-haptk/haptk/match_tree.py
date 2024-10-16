from ete3 import TreeStyle, TextFace, CircleFace,  NodeStyle, add_face_to_node

from haptk import utils

def match_tree_style():
    ts = TreeStyle()
    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.min_leaf_separation = 0
    ts.branch_vertical_margin = 10
    ts.optimal_scale_level = "full"
    ts.allow_face_overlap = True
    return ts

def circle_match_tree_style():
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

    return ts

def draw_match_tree(hst, t, output, samples, samples_to_tag, colors, proportions, circle):
    style = NodeStyle()
    style["size"] = 0

    if circle:
        for n in t.traverse():
            n.set_style(style)
        
            indexes = hst.get_node_indexes(n.name)
            if not n.is_leaf():

                if n.name == '0':
                    label = "ROOT"
                else:
                    label = len(indexes)

                F = TextFace(label, tight_text=True, penwidth=30, fsize=8)
                F.rotation = 270
                n.add_face(F, column=0)

            if n.is_leaf():
                if samples_to_tag == []:
                    n.set_style(utils.return_node_style("#000", 0))
                    n.img_style["bgcolor"] = "#FFF"

                for (color_i, sample_list) in enumerate(samples_to_tag):
                    results = [i for i in indexes if hst.samples[i] in sample_list]
                    if results:
                        n.set_style(utils.return_node_style(colors[color_i], 0))
                        n.img_style["bgcolor"] = colors[color_i]

    else:
        for n in t.traverse():
            n.set_style(style)

            indexes = hst.get_node_indexes(n.name)
            if n.name == '0':
                label = "ROOT"
            else:
                label = len(indexes)

                if proportions:
                    prc = label / len(samples)
                    label = f"{label} ({prc:.2f})"

            F = TextFace(label, tight_text=True, penwidth=30)
            F.rotation = 270
            n.add_face(F, column=0)

            if n.is_leaf():

                if samples_to_tag == []:
                    n.set_style(utils.return_node_style("#000", 0))
                    n.img_style["bgcolor"] = "#FFF"

                for (color_i, sample_list) in enumerate(samples_to_tag):
                    results = [i for i in indexes if hst.samples[i] in sample_list]
                    if results:
                        n.set_style(utils.return_node_style(colors[color_i], 0))
                        n.img_style["bgcolor"] = colors[color_i]


    return t


