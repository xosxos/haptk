import argparse

# Import the HAPTK python library
import haptk

parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('-i','--idx', type=int)    

args = parser.parse_args()

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)

# Get all samples from a node index
samples = hst.node_samples([args.idx])

for sample_list in samples:
    for sample in sample_list:
        print(sample)
        
