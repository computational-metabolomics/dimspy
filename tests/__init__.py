import os
import unittest
from collections import OrderedDict

import dimspy

import pprint


pp = pprint.PrettyPrinter(indent=4)


class test_dimspy(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test(self):
        pass


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMSnPyFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

