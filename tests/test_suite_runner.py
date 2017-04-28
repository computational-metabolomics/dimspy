#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_runner: provide functions to run test suites

author(s): Albert
origin: 04-29-2017

"""


import unittest, logging


def runTextSuite(suite, file_name, mode = 'w', **kwargs):
    with open(file_name, mode) as f:
        runner = unittest.TextTestRunner(stream = f, **kwargs)
        runner.run(suite)

def runHtmlSuite(suite, file_name, mode = 'w', **kwargs):
    try: from HTMLTestRunner import HTMLTestRunner
    except ImportError: logging.error('required module [HTMLTestRunner] not found, terminated'); return

    with open(file_name, mode) as f:
        runner = HTMLTestRunner(stream = f, **kwargs)
        runner.run(suite)

