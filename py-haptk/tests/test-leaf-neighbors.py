from haptk import graph
from haptk.graph import leaf_neighbors, read_hst

path = "data/bhst.hst.gz"
hst = read_hst(path)

neighbors = leaf_neighbors(hst)

print(neighbors)
