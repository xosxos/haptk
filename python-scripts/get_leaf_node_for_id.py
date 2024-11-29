import argparse

# Import the HAPTK python library
import haptk

parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--id', type=str)

args = parser.parse_args()

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)

idx = hst.find_leaf_node_for_id(args.id)

if idx != -1:
    data = hst.get_node_data(idx)
    print(data['stop']['pos'] - data['start']['pos'])

