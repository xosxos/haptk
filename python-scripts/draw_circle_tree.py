import argparse

# Import the HAPTK python library
import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--ids', nargs="+", type=str)
parser.add_argument('--tag-indexes', action="store_true")
parser.add_argument('--colors', nargs="+", type=str)
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

if args.colors:
    colors = args.colors
else:
    # kyme, kanta-hame, karjala, paijat-hame, pirkanmaa
    colors=["#866adb", "#f7708b", "#6ec9b8", "#ffe1bd", "#6ed4ff"]
    # colors=["#c34757", "#9a69b1", "#4071ab", "#90a720", "#d75bae"],


# Render the tree
hst.circle_tree(args.output,  to_tag=samples_to_tag, colors=colors, min_size=args.min_size, hard_cut=args.hard_cut, tag_indexes=args.tag_indexes)




