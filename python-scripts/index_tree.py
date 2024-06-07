import argparse
from haptk import graph


### Parameters
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

# Read .hst.gz and wrangle it to an ete3 tree
hst = graph.read_hst(args.hst)

graph.index_tree(hst, args.output, min_size=args.min_size, hard_cut=args.hard_cut)



