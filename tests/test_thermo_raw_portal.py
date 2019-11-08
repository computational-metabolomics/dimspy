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


import os
import unittest
import platform

from dimspy.portals.thermo_raw_portal import ThermoRaw


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset", *args)

def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class ThermoRawPortalsTestCase(unittest.TestCase):

    def test_thermo_raw_portal(self):

        run = ThermoRaw(to_test_data("raw", "batch04_QC17_rep01_262.RAW"))

        if platform.system() == "Darwin":
            self.assertEqual(str(run.timestamp), "02/04/2011 03:28:02")
        else:
            self.assertEqual(str(run.timestamp), "4/2/2011 3:28:02 AM")

        self.assertListEqual(list(run.headers().keys()), ['FTMS + p ESI w SIM ms [70.00-170.00]',
                                                          'FTMS + p ESI w SIM ms [140.00-240.00]',
                                                          'FTMS + p ESI w SIM ms [210.00-310.00]',
                                                          'FTMS + p ESI w SIM ms [280.00-380.00]',
                                                          'FTMS + p ESI w SIM ms [350.00-450.00]',
                                                          'FTMS + p ESI w SIM ms [420.00-520.00]',
                                                          'FTMS + p ESI w SIM ms [490.00-590.00]'])
        self.assertListEqual(list(run.scan_ids().keys()), list(range(1,89)))
        self.assertListEqual(list(run.tics().values())[0:2], [39800032.0, 38217892.0])
        self.assertEqual(len(run.tics()), 88)
        self.assertListEqual(list(run.ion_injection_times().values())[0:2], [40.434, 40.095])
        self.assertEqual(len(run.ion_injection_times()), 88)
        self.assertListEqual(run.scan_dependents(), [])
        pl = run.peaklist(1)
        self.assertEqual(pl.ID, 1)
        self.assertEqual(pl.metadata["header"], "FTMS + p ESI w SIM ms [70.00-170.00]")
        self.assertEqual(pl.metadata["ms_level"], 1.0)
        self.assertEqual(pl.metadata["ion_injection_time"], 40.434)
        self.assertEqual(pl.metadata["scan_time"], 0.5010899999999999)
        self.assertEqual(pl.metadata["elapsed_scan_time"], 1.05)
        self.assertEqual(pl.metadata["tic"], 39800032.0)
        self.assertEqual(pl.metadata["function_noise"], "noise_packets")
        self.assertEqual(pl.metadata["mz_range"], [70.0, 170.0])
        run.close()


if __name__ == '__main__':
    unittest.main()
