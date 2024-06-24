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

def is_maj_up_to(hst, node, value):
    iter_node = node
    # print(f"checking {iter_node.name} for majority")
    for _ in range(0, value):
        if iter_node.name == "0":
            return True

        if not iter_node.up.name == "0":
            # print("input node:", iter_node.up)
            is_big = is_majority_node(hst, iter_node.up)
            if is_big == False:
                # print(f"run {_}: {iter_node.up.name} is not big")
                return False

        iter_node = iter_node.up
            
    return True

def is_majority_node(hst, node):
    parent = node.up
    input_node_data = hst.get_node_data(int(node.name))
    value = len(input_node_data["indexes"])
    # print(f"input node {node.name} has {value} indexes")

    siblings = {}
    for child in parent.children:
        data = hst.get_node_data(int(child.name))
        if not child.name == node.name:
            siblings[child.name] = len(data["indexes"])


    # print(value, siblings.values())
    for v in siblings.values():
        if value < v: 
            return False
  

def get_children_idx_amounts(hst, node):
    return [len(hst.get_node_data(c.name)["indexes"]) for c in node.children]

def draw_circle_tree(hst, t, output, samples, samples_to_tag, colors):
    style = NodeStyle()
    style["size"] = 0

    narea_color = "#DAF7A6"
    larea_color = "#FFB6C1"

    spent_parents = []
    spent_grand_parents = []

    for n in t.traverse():
        n.set_style(style)
        
        if not n.is_leaf():
            indexes = hst.get_node_indexes(n.name)
            label = len(indexes)
            F = TextFace(label, tight_text=True, penwidth=30, fsize=8)
            F.rotation = 270
            n.add_face(F, column=0)

            count = 0
            for i in get_children_idx_amounts(hst, n):
                if i > 80:
                    count += 1

            if count > 1:
                haplotype = hst.get_node_haplotype(n.name)
                F = TextFace(f"{len(haplotype)}", fsize=30, tight_text=True, fgcolor="red")
                # F.rotation = 270
                n.add_face(F, column=0)
                
                

        if n.is_leaf():
            if samples_to_tag == []:
                n.set_style(utils.return_node_style("#000", 0))
                n.img_style["bgcolor"] = "#FFF"

            for (color_i, sample_list) in enumerate(samples_to_tag):
                results = [i for i in indexes if samples[i] in sample_list]
                if results:
                    n.set_style(utils.return_node_style(colors[color_i], 0))
                    n.img_style["bgcolor"] = colors[color_i]

            if n.up:
                parent = n.up
                parent_indexes = hst.get_node_indexes(parent.name)

                if len(parent_indexes) == 2:
                    if parent in spent_parents:
                        if is_maj_up_to(hst, n, 14):
                            # print(f"{n.name} is_maj")

                            haplotype = hst.get_node_haplotype(n.name)
                            start_idx = hst.get_node_start_idx(n.name)
                            stop_idx = hst.get_node_stop_idx(n.name)
                            start = hst.coords[start_idx]
                            stop = hst.coords[stop_idx]
                            # print(start, stop)
                            bp = stop['pos'] - start['pos']
                            F = TextFace(f"{len(haplotype)}, {bp}", fsize=30, fgcolor="white")
                            # F.margin_left = 50
                            n.add_face(F, column=0)

                            # n.set_style(utils.return_node_style("white", 0))
                            # n.img_style["bgcolor"] = "white"
                    else:
                        spent_parents.append(parent)


    return t

