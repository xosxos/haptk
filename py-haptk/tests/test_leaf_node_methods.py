import unittest
import haptk

class testLeafMethods(unittest.TestCase):
    def init(self):
        path = "tests/data/bhst.hst.gz"
        hst = haptk.read_hst(path)
        return hst

    def test_leaf_neighbors(self):
        hst = self.init()
        neighbors = hst.leaf_neighbors()

        self.assertEqual(neighbors, [[44, 45], [52, 53]])

    def test_leaf_nodes(self):
        hst = self.init()
        leaf_nodes = hst.leaf_nodes()

        self.assertEqual(leaf_nodes, [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 44, 45, 47, 49, 51, 52, 53])
       
    def test_nodes_data_on_leaf_nodes(self):
        hst = self.init()
        leaf_nodes = hst.leaf_nodes()
        data = hst.get_nodes_data(leaf_nodes)
        # print(data)
        # for d in data:
            # print("nmarkers", len(d['haplotype']))

    def test_node_samples(self):
        hst = self.init()
        leaf_nodes = hst.leaf_nodes()
        samples = hst.node_samples(leaf_nodes)
        correct_samples = [
            ['SAMPLE1'], ['SAMPLE2'], ['SAMPLE14'], ['SAMPLE3'], ['SAMPLE1'], ['SAMPLE4'],
            ['SAMPLE2'], ['SAMPLE5'], ['SAMPLE3'], ['SAMPLE6'], ['SAMPLE4'], ['SAMPLE7'],
            ['SAMPLE5'], ['SAMPLE8'], ['SAMPLE6'], ['SAMPLE9'], ['SAMPLE7'], ['SAMPLE10'],
            ['SAMPLE8'], ['SAMPLE11'], ['SAMPLE9'], ['SAMPLE13'], ['SAMPLE12'], ['SAMPLE10'],
            ['SAMPLE11'], ['SAMPLE12'], ['SAMPLE14'], ['SAMPLE13']
        ]
        self.assertEqual(samples, correct_samples)
        


if __name__ == '__main__':
    unittest.main()

