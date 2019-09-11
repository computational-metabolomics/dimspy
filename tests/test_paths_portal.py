#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_txt_portal

author(s): Albert
origin: 05-14-2017

"""


import unittest, os
from dimspy.portals import paths


class PathsPortalsTestCase(unittest.TestCase):
    def test_paths_portal(self):

        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset")

        source_raw = os.path.join(path, "raw")
        fn_filelist_raw = os.path.join(path, "filelist_raw_triplicates.txt")
        fns = paths.check_paths(fn_filelist_raw, source_raw)
        fns_c = [os.path.join(source_raw, 'batch04_QC17_rep01_262.RAW'),
                 os.path.join(source_raw, 'batch04_QC17_rep02_263.RAW'),
                 os.path.join(source_raw, 'batch04_QC17_rep03_264.RAW')]
        self.assertListEqual(fns, fns_c)

        fns = [os.path.join(source_raw, "batch04_QC17_rep01_262.RAW")]
        fns_out = paths.check_paths(None, fns)
        self.assertListEqual(fns, fns_out)

        fns = [os.path.join(source_raw, "batch04_QC17_rep01_262.RAW"),
               os.path.join(source_raw, "batch04_QC17_rep02_263.RAW"),
               os.path.join(source_raw, "batch04_QC17_rep03_264.RAW")]
        fns_out = paths.check_paths(fn_filelist_raw, fns)
        self.assertListEqual(fns, fns_out)

        source_mzml = os.path.join(path, "mzml")
        fns = [os.path.join(source_mzml, 'batch04_QC17_rep01_262.mzML')]
        fns_out = paths.check_paths(None, fns)
        self.assertListEqual(fns, fns_out)

        fn_filelist_mzml = os.path.join(path, "filelist_mzml_triplicates.txt")
        source_mzml_fns = [os.path.join(source_mzml, "batch04_QC17_rep01_262.mzML"),
                           os.path.join(source_mzml, "batch04_QC17_rep02_263.mzML"),
                           os.path.join(source_mzml, "batch04_QC17_rep03_264.mzML")]

        with self.assertRaises(IOError):
            paths.check_paths(fn_filelist_mzml, source_mzml_fns)

        with self.assertRaises(IOError):
            paths.check_paths(fn_filelist_mzml, source_mzml)


if __name__ == '__main__':
    unittest.main()
