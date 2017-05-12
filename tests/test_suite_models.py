#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_models

author(s): Albert
origin: 04-29-2017

"""


import unittest, os
import test_peaklist_metadata, test_peaklist_tags, test_peaklist, test_peak_matrix
from test_suite_runner import runTestSuite


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_peaklist_metadata))
    suite.addTest(unittest.findTestCases(test_peaklist_tags))
    suite.addTest(unittest.findTestCases(test_peaklist))
    suite.addTest(unittest.findTestCases(test_peak_matrix))

    report = os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), 'test_results', 'results_test_suite_models')
    runTestSuite(suite, report, title = 'Models Test Suite Report', verbosity = 2)
