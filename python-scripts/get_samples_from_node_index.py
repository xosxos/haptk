from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx

import hst_utils

### Parameters
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('-i','--idx', type=int)    

args = parser.parse_args()

HST, coords, samples, metadata = hst_utils.read_path_to_hst(args.hst)

node_data = HST.get_node_data(args.idx)

for i in node_data["indexes"]:
    print(samples[i])
