#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_suite_runner: provide functions to run test suites

author(s): Albert Zhou, Ralf Weber
origin: 04-29-2017

"""


import logging
from inspect import getargspec


def runTestSuite(suite, filename, mode = 'w', **kwargs):
    try:
        from HTMLTestRunner import HTMLTestRunner as runner
        filename += '.html'
    except ImportError:
        logging.warning('Module [HTMLTestRunner] not found, use default test runner.\n'
                        'Available: http://tungwaiyip.info/software/HTMLTestRunner_0_8_2/HTMLTestRunner.py')
        from unittest import TextTestRunner as runner
        filename += '.txt'

    with open(filename, mode) as f:
        kwargs = dict((k, v) for k,v in kwargs.items() if k in getargspec(runner.__init__)[0] and k != 'self')
        runner = runner(stream = f, **kwargs)
        runner.run(suite)

