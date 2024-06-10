from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--data', type=str)    
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

df = pd.read_csv(args.data)
df = df.dropna()
print(df)

def get_repeat_n(data, samples, df):
    sample_idx = data['indexes'][0]
    sample_name = samples[sample_idx]
    n = df.loc[df['id'] == sample_name]['n']
    return n


hst = haptk.read_hst(args.hst)
neighbors = hst.leaf_neighbors()

node_n = []
neighbor_n = []
for pair in neighbors:
    repeat_n = []
    for n in pair:
        data = hst.G.get_node_data(n)
        node_repeat_n = get_repeat_n(data, hst.samples, df)
        repeat_n.append(node_repeat_n)

    first = repeat_n[0]
    rest = repeat_n[1:]
    node_n.append(node_n)
    neighbor_n.append(sum(rest)/len(rest))
        

df = pd.DataFrame(
    {'repeat_n': node_n,
     'neighbor_avg': neighbor_n,
    })

df["status"] = np.select(
    [
        df['repeat_n'] >= 200,
        df['repeat_n'].between(150, 199),
        df['repeat_n'].between(100, 149),
        df['repeat_n'].between(50, 99),
        df['repeat_n'] < 50
    ],
    ['+200', '150-199', '100-149', '50-99', '0-49'],
    np.nan
)

# print(over_200['neighbor_avg'].mean())
# print(btw_150_200['neighbor_avg'].mean())

print(df.loc[df['status'] == "+200"]['neighbor_avg'].mean())
print(df.loc[df['status'] == "150-199"]['neighbor_avg'].mean())
print(df.loc[df['status'] == "100-149"]['neighbor_avg'].mean())
print(df.loc[df['status'] == "50-99"]['neighbor_avg'].mean())
print(df.loc[df['status'] == "0-49"]['neighbor_avg'].mean())

sns.regplot(x=df['repeat_n'], y=df['neighbor_avg'])
plt.savefig("regression.png")

         
