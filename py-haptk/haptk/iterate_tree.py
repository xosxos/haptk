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

    return ts

def iterate_tree_inner(hst, t, samples_to_tag):
    style = NodeStyle()
    style["size"] = 0

    for n in t.traverse():
        n.set_style(style)
        
        indexes = hst.get_node_indexes(n.name)

        tot_cases, tot_ctrls = 0, 0
        n_cases, n_ctrls = [], []

        for (i, sample_list) in enumerate(samples_to_tag):
            if i == 0:
                tot_cases = len(sample_list)
                n_cases = [i for i in indexes if hst.samples[i] in sample_list]
            if i == 1:
                tot_ctrls = len(sample_list)
                n_ctrls = [i for i in indexes if hst.samples[i] in sample_list]

        n_cases = len(n_cases)
        n_ctrls = len(n_ctrls)
        # res = fisher_exact([[n_cases, n_ctrls], [tot_cases - n_cases, tot_ctrls - n_ctrls]], alternative="greater")
        res = fisher_exact([[n_cases, n_ctrls], [tot_cases - n_cases, tot_ctrls - n_ctrls]])
        
        print(f"{n.name},{res.pvalue},{n_cases},{n_ctrls},{tot_cases},{tot_ctrls}")
        if res.pvalue < 0.0005:
            label = len(indexes)
            F = TextFace(f"{res.pvalue:.4}", tight_text=True, penwidth=30, fsize=8, fgcolor='red')
            F.rotation = 270
            F.background.color = "#FFF"
            F.margin_right = 10
            n.add_face(F, column=0)

        # if n.is_leaf() and res.pvalue < 0.005:
        label = len(indexes)
        F = TextFace(label, tight_text=True, penwidth=30, fsize=8)
        F.rotation = 270

        mix = res.pvalue * 1000
        if mix > 1:
            mix = 1
        color = utils.color_fader('red', 'blue', mix)

        F.background.color = color
        n.add_face(F, column=0)

    return t
