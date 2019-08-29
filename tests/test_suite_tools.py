#!/usr/bin/env python
#  -*- coding: utf-8 -*-


import unittest, sys
from . import test_tools
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_tools))

    report = os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), 'test_results', 'results_test_suite_tools')
    runTestSuite(suite, report, title = 'Process Test Suite Report', verbosity = 2)
