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


def draw_normal_tree(hst, t, samples_to_tag, colors, proportions, show_haplotype, n_markers, show_pos, branch_point_size, branch_length_as_majority):
    style = NodeStyle()
    style["size"] = 0

    # spent_parents = []
    for n in t.traverse():
        n.set_style(style)

        node_data = hst.get_node_data(int(n.name))

        label = len(node_data["indexes"])

        if proportions:
            top_data = hst.get_node_data(int('0'))
            samplen = len(top_data["indexes"])
            prc = label / samplen
            label = f"N={label} ({prc:.2f})"
        else:
            label = f"{label}"

        if show_haplotype:
            label = f"N={len(node_data["indexes"])}\n{hst.haplotype_string(node_data, n_markers)}"

        if show_pos:
            label = f"N={len(node_data["indexes"])}"
            haplotype = hst.haplotype_string(node_data, n_markers)
            start = node_data['start']['pos']
            stop = node_data['stop']['pos']
            label = f"{label}\n{haplotype}\nstart: {start}\nstop: {stop}"

        F = TextFace(label, tight_text=True, penwidth=30)
        F.rotation = 270
        n.add_face(F, column=0)

        # text_face = TextFace("", fsize=30, tight_text=True, fgcolor="red")
        # text_face.rotation = 270
        # text_face.margin_right = 10
        # utils.tag_branching_point_to_node(hst, n, text_face, branch_point_size)

        if n.is_leaf():
            if samples_to_tag == []:
                n.set_style(utils.return_node_style("#000", 1))

            for (color_i, sample_list) in enumerate(samples_to_tag):
                results = [i for i in node_data["indexes"] if hst.samples[i] in sample_list]
                if results:
                    n.set_style(utils.return_node_style(colors[color_i], 7 - color_i))

            # text_face = TextFace("", fsize=30, fgcolor="blue")
            # spent_parents = utils.tag_big_branch_to_node(hst, n, spent_parents, text_face, branch_length_as_majority)


    return t

