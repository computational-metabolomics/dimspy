#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017-2020 Ralf Weber, Albert Zhou.
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

sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))
from . import test_peak_filters, test_peak_alignment


if __name__ == '__main__':
    suite = unittest.TestSuite()

    suite.addTest(unittest.findTestCases(test_peak_alignment))
    suite.addTest(unittest.findTestCases(test_peak_filters))

    runner = unittest.TextTestRunner()
    runner.run(suite)
