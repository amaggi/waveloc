import unittest
from waveloc.clustering import compute_nbsta, do_clustering
import numpy as np


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(ClusteringTest))
    return suite


class ClusteringTest(unittest.TestCase):

    def setUp(self):

        self.event = 5
        self.threshold = 0.8
        self.coeff = {'A': [[1, 1, 0.3, 0.7, 0.8], [0, 1, 0.5, 0.2, 1],
                            [0, 0, 1, 0.6, 0.8], [0, 0, 0, 1, 0.4],
                            [0, 0, 0, 0, 1]],
                      'B': [[1, 1, 0.3, 0.7, 0.8], [0, 1, 0.5, 0.2, 1],
                            [0, 0, 1, 0.6, 0.8], [0, 0, 0, 1, 0.4],
                            [0, 0, 0, 0, 1]],
                      'C': [[1, 1, 0.3, 0.7, 0.8], [0, 1, 0.5, 0.2, 1],
                            [0, 0, 1, 0.6, 0.8], [0, 0, 0, 1, 0.4],
                            [0, 0, 0, 0, 1]]}
        self.nbsta = np.matrix([[3, 3, 0, 0, 3], [0, 3, 0, 0, 3],
                                [0, 0, 3, 0, 3], [0, 0, 0, 3, 0],
                                [0, 0, 0, 0, 3]])
        self.nbmin = 2
        self.cluster = {1: [1, 2, 3, 5]}

    def test_nbsta(self):
        exp_nbsta = compute_nbsta(self.event, self.coeff, self.threshold)
        self.assertEqual(exp_nbsta.all(), self.nbsta.all())

    def test_clustering(self):
        exp_cluster = do_clustering(self.event, self.nbsta, self.nbmin)
        self.assertEqual(exp_cluster, self.cluster)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
