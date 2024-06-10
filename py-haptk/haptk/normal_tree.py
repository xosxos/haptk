from ete3 import TreeStyle, TextFace, NodeStyle

from haptk import utils

def normal_tree_style():
    ts = TreeStyle()
    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.min_leaf_separation = 0
    ts.branch_vertical_margin = 5
    ts.optimal_scale_level = "full"
    ts.allow_face_overlap = True
    return ts

def draw_normal_tree(hst, t, output, samples, samples_to_tag, colors, proportions):
    style = NodeStyle()
    style["size"] = 0

    for n in t.traverse():
        n.set_style(style)

        node_data = hst.get_node_data(int(n.name))
        label = len(node_data["indexes"])

        if proportions:
            top_data = hst.get_node_data(int('0'))
            samplen = len(top_data["indexes"])
            prc = label / samplen
            label = f"{label} ({prc:.2f})"

        F = TextFace(label, tight_text=True, penwidth=30)
        F.rotation = 270
        n.add_face(F, column=0)

        if n.is_leaf():
            if samples_to_tag == []:
                n.set_style(utils.return_node_style("#000", 1))

            for (color_i, sample_list) in enumerate(samples_to_tag):
                results = [i for i in node_data["indexes"] if samples[i] in sample_list]
                if results:
                    n.set_style(utils.return_node_style(colors[color_i], 7 - color_i))


    return t

