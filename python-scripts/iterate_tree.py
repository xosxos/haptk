import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

# Import the HAPTK python library
import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--df', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--ids', nargs="+", type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)


df = pd.read_csv(args.df, sep=",")

print(hst.metadata)
# samples = hst.metadata

# df = df['id'].isin(samples)
# print(df)

# # Create list of samples to tag
# samples_to_tag = []
# if args.ids:
#     for file in args.ids:
#         ids = []
#         file = open(file, 'r')
#         for line in file.readlines():
#             ids.append(line.strip())

#         samples_to_tag.append(ids)

# # Render the tree
# hst.iterate_tree(df, args.output, w=4600, h=4600, to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut)




