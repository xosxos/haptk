import unittest
import haptk

class testGraphs(unittest.TestCase):
    def init(self):
        path = "tests/data/bhst.hst.gz"
        hst = haptk.read_hst(path)
        match_hst = haptk.read_hst(path)
        return hst, match_hst

    def test_index_tree(self):
        hst, _ = self.init()
        hst.index_tree("tests/results/index-tree.png")

    def test_match_tree(self):
        match_hst, hst = self.init()
        match_hst.match_tree(hst, "tests/results/match-tree.png")

    def test_circle_tree(self):
        hst, _ = self.init()
        hst.circle_tree("tests/results/circle-tree.png")

    def test_normal_tree(self):
        hst, _ = self.init()
        hst.normal_tree("tests/results/normal-tree.png")
        

if __name__ == '__main__':
    unittest.main()


