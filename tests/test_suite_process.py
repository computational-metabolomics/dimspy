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
from . import test_peak_alignment, test_peak_filters


if __name__ == '__main__':
    suite = unittest.TestSuite()

    # suite.addTest(unittest.findTestCases(test_peak_alignment))
    suite.addTest(unittest.findTestCases(test_peak_filters))

    runner = unittest.TextTestRunner()
    runner.run(suite)
