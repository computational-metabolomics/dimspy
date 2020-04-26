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


import os
import unittest
import zipfile
from dimspy.process.replicate_processing import read_scans


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset", *args)

def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class ReplicateProcessingTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        zip_ref = zipfile.ZipFile(to_test_data("MTBLS79_mzml_single.zip"), 'r')
        zip_ref.extractall(to_test_results("zip_data"))
        zip_ref.close()

    def test_read_scans(self):

        scans = read_scans(to_test_data("mzml", "batch04_QC17_rep01_262.mzML"), function_noise="median",
                           min_scans=1, filter_scan_events={"exclude": [["70.0", "170.0", "sim"]]})
        self.assertListEqual(list(scans.keys()), ['FTMS + p ESI w SIM ms [140.00-240.00]',
                                                  'FTMS + p ESI w SIM ms [210.00-310.00]',
                                                  'FTMS + p ESI w SIM ms [280.00-380.00]',
                                                  'FTMS + p ESI w SIM ms [350.00-450.00]',
                                                  'FTMS + p ESI w SIM ms [420.00-520.00]',
                                                  'FTMS + p ESI w SIM ms [490.00-590.00]'])

        scans = read_scans(to_test_data("mzml", "batch04_QC17_rep01_262.mzML"), function_noise="median",
                           min_scans=1, filter_scan_events={"include": [["70.0", "170.0", "sim"]]})
        self.assertListEqual(list(scans.keys()), ['FTMS + p ESI w SIM ms [70.00-170.00]'])

        scans = read_scans(to_test_data("mzml", "batch04_QC17_rep01_262.mzML"), function_noise="median",
                           min_scans=1, filter_scan_events={"exclude": ["FTMS + p ESI w SIM ms [70.00-170.00]"]})
        self.assertListEqual(list(scans.keys()), ['FTMS + p ESI w SIM ms [140.00-240.00]',
                                                  'FTMS + p ESI w SIM ms [210.00-310.00]',
                                                  'FTMS + p ESI w SIM ms [280.00-380.00]',
                                                  'FTMS + p ESI w SIM ms [350.00-450.00]',
                                                  'FTMS + p ESI w SIM ms [420.00-520.00]',
                                                  'FTMS + p ESI w SIM ms [490.00-590.00]'])

        scans = read_scans(to_test_data("mzml", "batch04_QC17_rep01_262.mzML"), function_noise="median",
                           min_scans=1, filter_scan_events={"include": ["FTMS + p ESI w SIM ms [70.00-170.00]"]})
        self.assertListEqual(list(scans.keys()), ['FTMS + p ESI w SIM ms [70.00-170.00]'])
