#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_models

author(s): Albert
origin: 04-29-2017

"""


import unittest, os
import test_peak_alignment, test_peak_filters
from test_suite_runner import runTestSuite


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_peak_alignment))
    suite.addTest(unittest.findTestCases(test_peak_filters))

    report = os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), 'test_results', 'results_test_suite_process')
    runTestSuite(suite, report, title = 'Process Test Suite Report', verbosity = 2)
