import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from ete3 import TreeStyle, TextFace

# Import the HAPTK python library
import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--df', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--min-start', type=int, default=None)    
parser.add_argument('--max-stop', type=int, default=None)    
parser.add_argument('--adjust-for-root', action="store_true")    
parser.add_argument('--rm', action="store_true")    
parser.add_argument('--ids', nargs="+", type=str)
parser.add_argument('--w', type=int)
parser.add_argument('--h', type=int)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

ADJUST_FOR_ROOT = args.adjust_for_root
RECESSIVE = args.rm

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)

df = pd.read_csv(args.df, sep=",")

orig_nsamples = df.shape[0]

# print(hst.metadata)

mask = df['id'].isin(hst.samples)
df = df[mask]

print(f'removed {orig_nsamples - df.shape[0]} samples not in the HST')
print(f'Samples left: {df.shape[0]}')
# print(df)

calculate_dur = True
# Calculate disease duration for those not deceased
if calculate_dur:
    # Up to date survival data is not available on 83 samples
    df["status"] = df.apply(lambda row: 0 if np.isnan(row["dur"]) else 1, axis=1)

    # Calculate duration from the year of onset (data has been last checked in 2021, thus subtract from 2021)
    check_yoo = lambda row: np.nan if np.isnan(row['yoo']) else 2021 - row['yoo']

    df["dur"] = df.apply(lambda row: check_yoo(row) if np.isnan(row["dur"]) else row["dur"], axis=1)
    df = df.drop("yoo", axis = 1)
    df.dropna(inplace=True)
else:
    df["status"] = 1
    df = df.drop("yoo", axis = 1)
    df.dropna(inplace=True)

print(f'Dropped NAs, samples left: {df.shape[0]}')
# for id in df['id']:
    # print(id)

df["site"] = df["site"].apply(lambda x: 0 if x == "SPINAL" else 1)
    
# Code sex as male=0 and female=1
df["sex"] = df["sex"].apply(lambda x: 0 if x == "M" else 1)

# print(df)

def tree_style(hst):
    ts = TreeStyle()

    if hst.metadata['hst_type'] == 'UhstLeft':
        ts.rotation = 180
    elif hst.metadata['hst_type'] == 'UhstRight':
        ts.rotation = 0
    elif hst.metadata['hst_type'] == 'Bhst':
        ts.rotation = 270

    ts.show_leaf_name = False
    ts.show_scale = False
    # ts.min_leaf_separation = 10
    ts.branch_vertical_margin = 25
    ts.optimal_scale_level = "full"
    ts.allow_face_overlap = False

    # ts.legend.add_face(TextFace("p < 0.05"), column=0)

    # ts.branch_vertical_margin = 10
    # ts.show_leaf_name = False
    # ts.show_scale = False

    # ts.mode = "c"
    # ts.arc_start = 0 # 0 degrees = 15:00 o'clock
    # ts.arc_span = 355
    # ts.min_leaf_separation = 0
    # ts.root_opening_factor = 0
    # ts.optimal_scale_level = "mid"
    # ts.allow_face_overlap = True
    # # ts.allow_face_overlap = False

    return ts


# Create list of samples to tag
samples_to_tag = []
if args.ids:
    for file in args.ids:
        ids = []
        file = open(file, 'r')
        for line in file.readlines():
            ids.append(line.strip())

        samples_to_tag.append(ids)

def optimizer_fisher(hst, _node_name, indexes, _df, samples_to_tag):
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

    res = fisher_exact([[n_cases, n_ctrls], [tot_cases - n_cases, tot_ctrls - n_ctrls]])

    #print(f"{_node_name},{res.pvalue},{n_cases},{n_ctrls},{tot_cases},{tot_ctrls}")

    return res.pvalue


def optimizer(hst, n, df, _samples_to_tag):
    node_data = hst.get_node_data(n.name)
    indexes = node_data['indexes']

    label = "o"
    if not n.is_leaf():
        label = len(indexes)
    else:
        label = "o"

    color = "#fafafa"
    # color = "#ffffff"
    text_color = "black"

    # print(n.name)
    if n.name != "0":

        names = [hst.get_sample_name(x) for x in indexes]

        df["gt"] = df["id"].apply(lambda id: names.count(id))

         # Additive model
        hets = list(df["gt"]).count(1)
        alt_homs = list(df["gt"]).count(2)
        sum = hets + alt_homs * 2
        freq = sum
        rows = df.loc[(df['gt'] == 1) | (df['gt'] ==2)]
        avg = f'{rows['dur'].mean():.2g}'
        # median = f'{rows['dur'].median():.2g}'
        clause = (hets + alt_homs) < 5

        # Recessive model
        if RECESSIVE:
            df["gt"] = df["gt"].apply(lambda x: 0 if x < 2 else 1)
            alt_homs = list(df["gt"]).count(1)
            freq = alt_homs
            rows = df.loc[df['gt'] == 1]
            avg = f'{rows['dur'].mean():.2g}'
            # median = f'{rows['dur'].median():.2g}'
            clause = alt_homs < 5

        if clause:
            label = f" N={freq:.2g} AVG={avg}\n{hst.haplotype_string(node_data, 5)}"
            return create_text_face(label, hst, color, text_color)


        # Adjust for root variant
        if ADJUST_FOR_ROOT:
            if n.name != "2" and n.name != "1":
                root_indexes = hst.get_node_indexes("2")
                root_names = [hst.get_sample_name(x) for x in root_indexes]
                df["root_gt"] = df["id"].apply(lambda id: root_names.count(id))

        # Drop IDs before fitting the model
        df = df.drop("id", axis = 1)
        df = df.drop("sod1", axis = 1)

        # Fit the model
        cph = CoxPHFitter()

        cph.fit(df, duration_col = 'dur', event_col = 'status')

        coef = round(cph.summary['coef'].iloc[4], 6)
        e_coef = round(cph.summary['exp(coef)'].iloc[4], 6)
        p = cph.summary['p'].iloc[4]
        label = f"N={freq:.3g} AVG={avg}\nHR={e_coef:.3g} P={p:.2g}\n{hst.haplotype_string(node_data, 5)}"
        # print(n.name, len(indexes), hets, alt_homs, coef, e_coef, p)

        if p < 0.05:
            color = "#ffbda6"
            color = '#ffb8b8'
            color = '#ffc4c4'
            color = '#ffcccc'

        if p < 0.0036:
            color = "#ff6f6f"
            color = '#ffc9c9'
            color = '#ff7c7c'
            color = "#ff8585"
            # color = "red"


        return create_text_face(label, hst, color, text_color)
    else:

        return create_root_face(hst, "#ffd3f2")

def create_text_face(label, hst, bgcolor, fgcolor):
        F = TextFace(label, penwidth=30, fsize=14, fgcolor="black")
        F.rotation = get_rotation(hst)
        F.inner_background.color = bgcolor
        F.fgcolor = fgcolor
        F.margin_right = 15
        return F

def create_root_face(hst, bgcolor):
        if hst.metadata['hst_type'] == 'UhstLeft':
            root_label = f"<- {hst.metadata['start_coord']}"
        elif hst.metadata['hst_type'] == 'UhstRight':
            root_label = f"{hst.metadata['start_coord']} ->"
        elif hst.metadata['hst_type'] == 'Bhst':
            root_label = f"<- {hst.metadata['start_coord']} ->"
        else:
            root_label = f"{hst.metadata.start_coord}"

        F = TextFace(root_label, penwidth=30, fsize=15, fgcolor="black")
        F.rotation = get_rotation(hst)

        F.background.color = bgcolor
        F.margin_right = 15
        F.margin_bottom = 15

        return F

def get_rotation(hst):
    if hst.metadata['hst_type'] == 'UhstLeft':
        rotation = 540
    elif hst.metadata['hst_type'] == 'UhstRight':
        rotation = 360
    elif hst.metadata['hst_type'] == 'Bhst':
        rotation = 270
    else:
        rotation = 0

    return rotation


# # Render the tree
if hst.metadata['hst_type'] == 'UhstRight':
    hst.iterate_tree(df, optimizer, args.output, w=args.w, h=args.h, to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut, min_start=args.min_start, max_stop=args.max_stop, tree_style=tree_style(hst), right_up=True)
else:
    hst.iterate_tree(df, optimizer, args.output, w=args.w, h=args.h, to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut, min_start=args.min_start, max_stop=args.max_stop, tree_style=tree_style(hst))



