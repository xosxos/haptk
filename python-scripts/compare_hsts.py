from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import json
import argparse
import rustworkx as rx
import gzip

from haptk import graph

### Parameters
parser = argparse.ArgumentParser()
parser.add_argument('--hst', type=str)
parser.add_argument('--match-hst', type=str)
parser.add_argument('--min-size', type=int)
parser.add_argument('--ids', type=str)
parser.add_argument('--proportions', action="store_true")
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

hst = graph.read_hst(args.hst)
match_hst = graph.read_hst(args.match_hst)

samples_to_tag = []
if args.ids:
    for file in args.ids:
        ids = []
        file = open(file, 'r')
        for line in file.readlines():
            ids.append(line.strip())

        samples_to_tag.append(ids)

graph.match_tree(match_hst, hst, args.output, to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut)
