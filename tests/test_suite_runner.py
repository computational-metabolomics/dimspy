#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_runner: provide functions to run test suites

author(s): Albert
origin: 04-29-2017

"""


import logging
from inspect import getargspec


def runTestSuite(suite, file_name, mode = 'w', **kwargs):
    try:
        from HTMLTestRunner import HTMLTestRunner as runner
        file_name += '.html'
    except ImportError:
        logging.warning('Module [HTMLTestRunner] not found, use default test runner.\n'
                        'Available: http://tungwaiyip.info/software/HTMLTestRunner_0_8_2/HTMLTestRunner.py')
        from unittest import TextTestRunner as runner
        file_name += '.txt'

    with open(file_name, mode) as f:
        kwargs = dict((k, v) for k,v in kwargs.items() if k in getargspec(runner.__init__)[0] and k != 'self')
        runner = runner(stream = f, **kwargs)
        runner.run(suite)

