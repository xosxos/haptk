from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip
# import pandas as pd
import numpy as np
# import seaborn as sns
from matplotlib import pyplot as plt

import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--cases', type=str)
parser.add_argument('--ctrls', type=str)

args = parser.parse_args()

# Create list of samples to tag
def read_ids(filename):
    file = open(filename, 'r')
    ids = []
    for line in file.readlines():
        ids.append(line.strip())
    return ids

def get_status(cases, controls, sample):
    if sample in cases:
        return 'CASE'
    elif sample in controls:
        return 'CTRL'
    else:
        return 'ERROR'

      
cases = read_ids(args.cases)
ctrls = read_ids(args.ctrls)

hst = haptk.read_hst(args.hst)
neighbors = hst.leaf_neighbors()

nneighbors = len(neighbors)
ncases = len(cases)
nctrls = len(ctrls)
total = ncases + nctrls

print(f"ncases: {ncases}, nctrls: {nctrls}, % cases: {(ncases / (nctrls+ncases)):.3} % ctrls: {(nctrls / (nctrls+ncases)):.3}")
pair_cases = ((ncases * (ncases - 1)) / 2 )
pair_ctrls = ((nctrls * (nctrls - 1)) / 2 )
pair_total = ((total * (total - 1)) / 2)
pair_mix = pair_total - pair_cases - pair_ctrls

print(f"pair_cases: {pair_cases/pair_total}, pair_ctrls: {pair_ctrls/pair_total}, pair_mix: {pair_mix/pair_total}")

print(f"n_neighbors: {nneighbors}")

full_case = 0
full_ctrl = 0
mix = 0
for pair in neighbors:
    ncases = 0
    nctrls = 0
    for n in pair:
        data = hst.G.get_node_data(n)
        sample_idx = data['indexes'][0]
        sample = hst.samples[sample_idx]

        match get_status(cases, ctrls, sample):
            case 'CASE':
                ncases += 1
            case 'CTRL':
                nctrls += 1
            case 'ERROR':
                print(f'error: sample {sample} not in cases nor controls')
                quit()

    if ncases > 0 and nctrls > 0:
        mix += 1
    elif ncases == 0 and nctrls > 0:
        full_ctrl += 1 
    elif nctrls == 0 and ncases > 0:
        full_case += 1
    else:
        print('error: logic not handled')
        quit()

print(f"full_case: {full_case}, full_ctrl: {full_ctrl}, mix: {mix}")
print(f"full_case: {full_case/nneighbors}, full_ctrl: {full_ctrl/nneighbors}, mix: {mix/nneighbors}")
