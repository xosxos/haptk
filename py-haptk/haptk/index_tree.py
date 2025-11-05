from ete3 import TreeStyle, TextFace, NodeStyle

def index_tree_style():
    ts = TreeStyle()
    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.min_leaf_separation = 0
    ts.branch_vertical_margin = 5
    ts.optimal_scale_level = "full"
    ts.allow_face_overlap = True
    return ts
    

def draw_index_tree(hst, t, output, leaf_as_name):
    style = NodeStyle()
    style["size"] = 0

    lstyle = NodeStyle()
    lstyle["fgcolor"] = "#000"
    lstyle["shape"] = "square"
    lstyle["size"] = 1
    lstyle["hz_line_type"] = 0
    lstyle["vt_line_type"] = 0

    for n in t.traverse():
        n.set_style(style)

        node_data = hst.get_node_data(int(n.name))

        if leaf_as_name and len(node_data["indexes"]) == 1:
            sample_index = node_data["indexes"][0]
            label = hst.get_sample_name(sample_index)
        else:
            label = n.name

        node_data = hst.get_node_data(int(n.name))

        F = TextFace(label, tight_text=True, penwidth=30)
        F.rotation = 270
        n.add_face(F, column=0)
        n.set_style(lstyle)

    return t
