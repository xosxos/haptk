from ete3 import TreeStyle, TextFace, CircleFace,  NodeStyle, add_face_to_node

from haptk import utils

def circle_tree_style():
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

def draw_circle_tree(hst, t, samples_to_tag, colors, tag_indexes):
    style = NodeStyle()
    style["size"] = 0

    narea_color = "#DAF7A6"
    larea_color = "#FFB6C1"

    for n in t.traverse():
        n.set_style(style)
        
        indexes = hst.get_node_indexes(n.name)
        if not n.is_leaf():
            label = len(indexes)
            F = TextFace(label, tight_text=True, penwidth=30, fsize=8)
            F.rotation = 270
            n.add_face(F, column=0)

        if n.is_leaf():
            if samples_to_tag == []:
                n.set_style(utils.return_node_style("#000", 0))
                n.img_style["bgcolor"] = "#FFF"

            for (color_i, sample_list) in enumerate(samples_to_tag):
                if tag_indexes:
                    results = [i for i in indexes if str(i) in sample_list]
                else:
                    results = [i for i in indexes if hst.get_sample_name(i) in sample_list]

                if results:
                    n.set_style(utils.return_node_style(colors[color_i], 0))
                    n.img_style["bgcolor"] = colors[color_i]


    return t

