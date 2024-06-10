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
    

def draw_index_tree(t, output):
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

        label = n.name

        F = TextFace(label, tight_text=True, penwidth=30)
        F.rotation = 270
        n.add_face(F, column=0)
        n.set_style(lstyle)

    return t
