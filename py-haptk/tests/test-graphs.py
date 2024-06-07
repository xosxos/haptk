from haptk import graph

path = "data/bhst.hst.gz"
hst = graph.read_hst(path)

match_hst = graph.read_hst(path)

graph.index_tree(hst, "results/index-tree.png")
graph.circle_tree(hst, "results/circle-tree.png")
graph.normal_tree(hst, "results/normal-tree.png")
graph.match_tree(match_hst, hst, "results/match-tree.png")

