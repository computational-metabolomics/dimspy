#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017-2019 Ralf Weber, Albert Zhou.
#
# This file is part of DIMSpy.
#
# DIMSpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIMSpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DIMSpy.  If not, see <https://www.gnu.org/licenses/>.
#


import sys
import unittest
from pathlib import Path

from . import test_tools

sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_tools))

    report = os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), 'results', 'results_test_suite_tools')
    runTestSuite(suite, report, title = 'Process Test Suite Report', verbosity = 2)
