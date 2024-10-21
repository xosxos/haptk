import argparse

# Import the HAPTK python library
import haptk

parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)

# Render the tree
hst.index_tree(args.output, min_size=args.min_size, hard_cut=args.hard_cut)



