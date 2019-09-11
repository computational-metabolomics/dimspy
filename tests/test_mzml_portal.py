#!/usr/bin/env python
#  -*- coding: utf-8 -*-


import unittest
import os
import zipfile

from dimspy.portals.mzml_portal import Mzml


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)

def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class MzmlPortalsTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        zip_ref = zipfile.ZipFile(to_test_data("mzml_DIMSn.zip"), 'r')
        zip_ref.extractall(to_test_results("zip_data", "mzml"))
        zip_ref.close()

    def test_mzml_portal(self):

        run = Mzml(to_test_data("MTBLS79_subset", "mzml", "batch04_QC17_rep01_262.mzML"))
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

        pl = run.peaklist(1)
        self.assertEqual(pl.ID, 1)
        self.assertEqual(pl.metadata["header"], "FTMS + p ESI w SIM ms [70.00-170.00]")
        self.assertEqual(pl.metadata["ms_level"], 1.0)
        self.assertEqual(pl.metadata["ion_injection_time"], 40.433891296387)
        self.assertEqual(pl.metadata["scan_time"], 0.50109)
        self.assertEqual(pl.metadata["tic"], 39800032.0)
        self.assertEqual(pl.metadata["function_noise"], "median")
        self.assertEqual(pl.metadata["mz_range"], [70.0, 170.0])
        run.close()

        run = Mzml(to_test_results("zip_data", "mzml", "A08_Apolar_Daph_AMP1_C30_LCMS_Pos_DIMSn_subset.mzML"))
        sd = run.scan_dependents()
        self.assertListEqual(list(run.tics().values())[0:2], [120293696.0, 13602.5234375])
        self.assertEqual(len(run.tics()), 36)
        self.assertListEqual(sd[0], [1, 3])
        self.assertListEqual(sd[-1], [511, 512])
        self.assertEqual(len(sd), 30)
        run.close()

    @classmethod
    def tearDownClass(cls):

        import shutil
        shutil.rmtree(to_test_results(""))
        os.makedirs(to_test_results(""))


if __name__ == '__main__':
    unittest.main()
