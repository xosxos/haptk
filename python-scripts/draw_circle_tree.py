import argparse

# Import the HAPTK python library
import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--branch-point-size', type=int, default=9999999)    
parser.add_argument('--branch-length', type=int, default=99999999)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--ids', nargs="+", type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)

# Create list of samples to tag
samples_to_tag = []
if args.ids:
    for file in args.ids:
        ids = []
        file = open(file, 'r')
        for line in file.readlines():
            ids.append(line.strip())

        samples_to_tag.append(ids)

# Render the tree
hst.circle_tree(args.output, colors=["red", "blue"], to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut, branch_point_size=args.branch_point_size, branch_length_as_majority=args.branch_length)




