#!/usr/bin/env python
#  -*- coding: utf-8 -*-


import unittest
import os
from dimspy.portals.mzml_portal import Mzml


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset", *args)

def to_test_result(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)


class MzmlPortalsTestCase(unittest.TestCase):

    def test_mzml_portal(self):

        run = Mzml(to_test_data("mzml", "batch04_QC17_rep01_262.mzML"))
        self.assertEqual((run.run.get_spectrum_count(), run.run.get_spectrum_count()), (88, 88))
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
        self.assertListEqual(list(run.ion_injection_times().values())[0:2], [40.433891296387, 40.094646453857])
        self.assertEqual(len(run.ion_injection_times()), 88)
        self.assertListEqual(run.scan_dependents(), [])
        run.close()

if __name__ == '__main__':
    unittest.main()
