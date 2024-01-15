import argparse
import rustworkx as rx

def read_path_to_hst(path):
    with gzip.open(path, 'rt') as f:
        json_bytes = f.read()
        hst_struct = json.loads(json_bytes)
        hst = hst_struct["hst"]
        coords = hst_struct["coords"]
        samples = hst_struct["samples"]
        metadata = hst_struct["metadata"]

        G = rx.PyDiGraph()

        edges_tuple_vec = []
        for edges in hst["edges"]:
            edges_tuple_vec.append((edges[0], edges[1], edges[2]))

        G.add_nodes_from(hst["nodes"])
        G.add_edges_from(edges_tuple_vec)

        if metadata["selection"] == "All":
            samples = [val for val in samples for _ in (0, 1)]

        return G, coords, samples, metadata

### Parameters
parser = argparse.ArgumentParser()
parser.add_argument('--hst', type=str)
parser.add_argument('--match-hst', type=str)

args = parser.parse_args()

# Read in both .hst.gz
HST, coords1, samples1, metadata1 = read_path_to_hst(args.hst)
HST2, coords2, samples2, metadata2 = read_path_to_hst(args.match_hst)

if coords1 != coords2:
    print("coords error")

if samples1 != samples2:
    print("samples error")

if metadata1 != metadata2:
    print("metadata error")

e1 = HST.num_edges()
e2 = HST2.num_edges()

n1 = HST.num_nodes()
n2 = HST2.num_nodes()

print(f"e1: {e1}, e2: {e2}")
print(f"n1: {n1}, n2: {n2}")

el1 = HST.edge_list()
el2 = HST2.edge_list()
nl1 = HST.nodes()
nl2 = HST2.nodes()

emin = min(e1, e2)

for i in range(0, emin):
    if el1[i] != el2[i] or nl1[i] != nl2 [i]:
        print(f"failed equality {i}")
        # print(f"edge  {i}: {el1[i]} != {el2[i]}")
        # print(f"node  {i}: {nl1[i+1]} != {nl2[i+1]}")
        # break
