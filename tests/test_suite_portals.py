#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_portals

author(s): Albert
origin: 05-14-2017

"""


import unittest, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))
from . import test_txt_portal, test_hdf5_portal#, test_paths_portal


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_txt_portal))
    suite.addTest(unittest.findTestCases(test_hdf5_portal))
    #suite.addTest(unittest.findTestCases(test_paths_portal))

    runner = unittest.TextTestRunner()
    runner.run(suite)

