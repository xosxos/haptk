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
haplotype = hst.get_node_haplotype(args.idx)
start_idx = hst.get_node_start_idx(args.idx)
stop_idx = hst.get_node_stop_idx(args.idx)

coords = hst.coords[start_idx:stop_idx+1]

print("contig,pos,ref,alt,gt")
for coord, ht in zip(coords, haplotype):
    k = f"{coord['contig']},{coord['pos']},{coord['reference']},{coord['alt']},{ht}"
    print(k)
    


