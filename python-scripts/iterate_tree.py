import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from ete3 import TreeStyle

# Import the HAPTK python library
import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--df', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--ids', nargs="+", type=str)
parser.add_argument('--w', type=int)
parser.add_argument('--h', type=int)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)

df = pd.read_csv(args.df, sep=",")

# print(hst.metadata)
samples = hst.metadata["samples"]

mask = df['id'].isin(samples)
df = df[mask]


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


df["site"] = df["site"].apply(lambda x: 0 if x == "SPINAL" else 1)
    
# Code sex as male=0 and female=1
df["sex"] = df["sex"].apply(lambda x: 0 if x == "M" else 1)

# print(df)

def tree_style():
    ts = TreeStyle()

    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.min_leaf_separation = 35
    # ts.branch_vertical_margin = 10
    ts.optimal_scale_level = "mid"
    ts.allow_face_overlap = False

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
            n_cases = [i for i in indexes if hst.metadata["samples"][i] in sample_list]
        if i == 1:
            tot_ctrls = len(sample_list)
            n_ctrls = [i for i in indexes if hst.metadata["samples"][i] in sample_list]

    n_cases = len(n_cases)
    n_ctrls = len(n_ctrls)

    res = fisher_exact([[n_cases, n_ctrls], [tot_cases - n_cases, tot_ctrls - n_ctrls]])

    #print(f"{_node_name},{res.pvalue},{n_cases},{n_ctrls},{tot_cases},{tot_ctrls}")

    return res.pvalue

def optimizer(hst, n, indexes, df, _samples_to_tag):
    label = "o"
    if not n.is_leaf():
        label = len(indexes)
    else:
        label = "o"

    color = "blue"
    text_color = "white"

    # print(node_name)
    if n.name != "0":
        names = [hst.samples[x] for x in indexes]

        df["gt"] = df["id"].apply(lambda id: names.count(id))

         # # Additive model
        hets = list(df["gt"]).count(1)
        alt_homs = list(df["gt"]).count(2)

        if hets + alt_homs < 5:
            return (1.0, "NA", color, text_color)

        # Recessive model
        # df["gt"] = df["gt"].apply(lambda x: 0 if x < 2 else 1)
        # hets = list(df["gt"]).count(0)
        # alt_homs = list(df["gt"]).count(1)

        # if alt_homs < 5:
        #     return (1.0, "NA", color, text_color)

        # Drop IDs before fitting the model
        df = df.drop("id", axis = 1)

        # Fit the model
        cph = CoxPHFitter()

        cph.fit(df, duration_col = 'dur', event_col = 'status')

        coef = round(cph.summary['coef'].iloc[5], 6)
        e_coef = round(cph.summary['exp(coef)'].iloc[5], 6)
        p = cph.summary['p'].iloc[5]
        label = f"{e_coef:.3g}"
        print(n.name, len(indexes), hets, alt_homs, coef, e_coef, p)

        if p < 0.05:
            text_color = "white"
            color = "darkviolet"

        if p < 0.01:
            text_color = "black"
            color = "pink"

        if p < 0.005:
            text_color = "black"
            color = "red"

        return (float(p), label, color, text_color)
    else:
        return (1.0, "ROOT", "blue", "white")



# # Render the tree
hst.iterate_tree(df, optimizer, args.output, w=args.w, h=args.h, to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut, tree_style=tree_style())




