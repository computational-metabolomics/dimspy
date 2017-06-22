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

        source = os.path.join("data", "MTBLS79_subset", "MTBLS79_mzml_triplicates.zip")
        fns = paths.check_paths(None, source)
        fns_c = ['batch04_QC17_rep01_262.mzML', 'batch04_QC17_rep02_263.mzML', 'batch04_QC17_rep03_264.mzML',
                    'batch04_B02_rep01_301.mzML', 'batch04_B02_rep02_302.mzML', 'batch04_B02_rep03_303.mzML']
        self.assertListEqual(fns, fns_c)

        fn_filelist = os.path.join("data", "MTBLS79_subset", "filelist_mzml_triplicates.txt")
        fns = paths.check_paths(fn_filelist, source)
        fns_c = ['batch04_QC17_rep01_262.mzML', 'batch04_QC17_rep02_263.mzML', 'batch04_QC17_rep03_264.mzML',
                    'batch04_B02_rep01_301.mzML', 'batch04_B02_rep02_302.mzML', 'batch04_B02_rep03_303.mzML']
        self.assertListEqual(fns, fns_c)

        source_raw = os.path.join("data", "MTBLS79_subset", "raw")
        fn_filelist_raw = os.path.join("data", "MTBLS79_subset", "filelist_raw_triplicates.txt")
        fns = paths.check_paths(fn_filelist_raw, source_raw)
        fns_c = ['data/MTBLS79_subset/raw/batch04_QC17_rep01_262.RAW',
                    'data/MTBLS79_subset/raw/batch04_QC17_rep02_263.RAW',
                    'data/MTBLS79_subset/raw/batch04_QC17_rep03_264.RAW']
        self.assertListEqual(fns, fns_c)

        fns = paths.check_paths(None, [os.path.join(source_raw, "batch04_QC17_rep01_262.RAW")])
        fns_c = ['data/MTBLS79_subset/raw/batch04_QC17_rep01_262.RAW']
        self.assertListEqual(fns, fns_c)


        source_raw_fns = [os.path.join(source_raw, "batch04_QC17_rep01_262.RAW"),
                os.path.join(source_raw, "batch04_QC17_rep02_263.RAW"),
                os.path.join(source_raw, "batch04_QC17_rep03_264.RAW")]
        fns = paths.check_paths(fn_filelist_raw, source_raw_fns)
        fns_c = ['data/MTBLS79_subset/raw/batch04_QC17_rep01_262.RAW',
                'data/MTBLS79_subset/raw/batch04_QC17_rep02_263.RAW',
                'data/MTBLS79_subset/raw/batch04_QC17_rep03_264.RAW']
        self.assertListEqual(fns, fns_c)

        source_mzml = os.path.join("data", "MTBLS79_subset", "mzml")
        fns = paths.check_paths(None, os.path.join(source_mzml, "batch04_QC17_rep01_262.mzML"))
        fns_c = ['data/MTBLS79_subset/mzml/batch04_QC17_rep01_262.mzML']
        self.assertListEqual(fns, fns_c)

        fn_filelist_mzml = os.path.join("data", "MTBLS79_subset", "filelist_mzml_triplicates.txt")
        source_mzml_fns = [os.path.join(source_mzml, "batch04_QC17_rep01_262.mzML"),
                os.path.join(source_mzml, "batch04_QC17_rep02_263.mzML"),
                os.path.join(source_mzml, "batch04_QC17_rep03_264.mzML")]

        with self.assertRaises(IOError):
            paths.check_paths(fn_filelist_mzml, source_mzml_fns)

        with self.assertRaises(IOError):
            paths.check_paths(fn_filelist_mzml, source_mzml)

if __name__ == '__main__':
    unittest.main()
