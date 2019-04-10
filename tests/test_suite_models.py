#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_models

author(s): Albert
origin: 04-29-2017

"""


import unittest, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))
import test_peaklist_metadata, test_peaklist_tags, test_peaklist, test_peak_matrix


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_peaklist_metadata))
    suite.addTest(unittest.findTestCases(test_peaklist_tags))
    suite.addTest(unittest.findTestCases(test_peaklist))
    suite.addTest(unittest.findTestCases(test_peak_matrix))

    runner = unittest.TextTestRunner()
    runner.run(suite)
