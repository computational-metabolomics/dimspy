#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_portals

author(s): Albert
origin: 05-14-2017

"""


import unittest, os
import test_txt_portal, test_hdf5_portal, test_paths_portal
from test_suite_runner import runTestSuite


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_txt_portal))
    suite.addTest(unittest.findTestCases(test_hdf5_portal))
    suite.addTest(unittest.findTestCases(test_paths_portal))

    report = os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), 'test_results', 'results_test_suite_portals')
    runTestSuite(suite, report, title = 'Portals Test Suite Report', verbosity = 2)
