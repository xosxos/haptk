import argparse
import haptk
from haptk import graph

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--ids', type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

# Read .hst.gz and wrangle it to an ete3 tree
hst = graph.read_hst(args.hst)

samples_to_tag = []
if args.ids:
    for file in args.ids:
        ids = []
        file = open(file, 'r')
        for line in file.readlines():
            ids.append(line.strip())

        samples_to_tag.append(ids)

graph.circle_tree(hst, args.output, to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut)




