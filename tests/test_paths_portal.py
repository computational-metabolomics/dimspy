#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_txt_portal

author(s): Albert
origin: 05-14-2017

"""


import os
import unittest
import platform

from dimspy.portals import paths


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", *args)


def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "results", *args)


class PathsPortalsTestCase(unittest.TestCase):
    def test_paths_portal(self):

        files_correct = [to_test_data("MTBLS79_subset", "raw", "batch04_QC17_rep01_262.RAW"),
                          to_test_data("MTBLS79_subset", "raw", "batch04_QC17_rep02_263.RAW"),
                          to_test_data("MTBLS79_subset", "raw", "batch04_QC17_rep03_264.RAW")]
        tsv = to_test_data("MTBLS79_subset", "filelist_raw_triplicates.txt")

        source = to_test_data("MTBLS79_subset", "raw")
        files = paths.validate_and_sort_paths(source, tsv)
        self.assertListEqual(files, files_correct)

        source = to_test_data("MTBLS79_subset", "raw")
        files = paths.validate_and_sort_paths(source, tsv)
        self.assertListEqual(files, files_correct)

        source = [to_test_data("MTBLS79_subset", "raw", "batch04_QC17_rep03_264.RAW"),
                  to_test_data("MTBLS79_subset", "raw", "batch04_QC17_rep02_263.RAW"),
                  to_test_data("MTBLS79_subset", "raw", "batch04_QC17_rep01_262.RAW")]
        files = paths.validate_and_sort_paths(source, tsv)
        self.assertListEqual(files, files_correct)

        files = paths.validate_and_sort_paths(tsv=None, source=source)
        self.assertListEqual(files, source)

        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset")

        source_raw = os.path.join(path, "raw")
        fn_filelist_raw = os.path.join(path, "filelist_raw_triplicates.txt")
        fns = paths.validate_and_sort_paths(source_raw, fn_filelist_raw)
        fns_c = [os.path.join(source_raw, 'batch04_QC17_rep01_262.RAW'),
                 os.path.join(source_raw, 'batch04_QC17_rep02_263.RAW'),
                 os.path.join(source_raw, 'batch04_QC17_rep03_264.RAW')]
        self.assertListEqual(fns, fns_c)

        fns = [os.path.join(source_raw, "batch04_QC17_rep01_262.RAW")]
        fns_out = paths.validate_and_sort_paths(fns, None)
        self.assertListEqual(fns, fns_out)

        fns = [os.path.join(source_raw, "batch04_QC17_rep01_262.RAW"),
               os.path.join(source_raw, "batch04_QC17_rep02_263.RAW"),
               os.path.join(source_raw, "batch04_QC17_rep03_264.RAW")]
        fns_out = paths.validate_and_sort_paths(fns, fn_filelist_raw)
        self.assertListEqual(fns, fns_out)

        source_mzml = os.path.join(path, "mzml")
        fns = [os.path.join(source_mzml, 'batch04_QC17_rep01_262.mzML')]
        fns_out = paths.validate_and_sort_paths(fns, None)
        self.assertListEqual(fns, fns_out)

        fn_filelist_mzml = os.path.join(path, "filelist_mzml_triplicates.txt")
        source_mzml_fns = [os.path.join(source_mzml, "batch04_QC17_rep01_262.mzML"),
                           os.path.join(source_mzml, "batch04_QC17_rep02_263.mzML"),
                           os.path.join(source_mzml, "batch04_QC17_rep03_264.mzML")]

        with self.assertRaises(IOError):
            paths.validate_and_sort_paths(source_mzml_fns, fn_filelist_mzml)

        with self.assertRaises(IOError):
            paths.validate_and_sort_paths(source_mzml, fn_filelist_mzml)

    def test_sort_ms_files_by_timestamp(self):
        p = to_test_data("MTBLS79_subset", "mzml")
        ps = [os.path.join(p, fn) for fn in os.listdir(p)]
        files_sorted = paths.sort_ms_files_by_timestamp(ps)
        self.assertEqual(files_sorted[0], (os.path.join(p, "batch04_QC17_rep01_262.mzML"), '2011-04-02T03:28:02Z'))
        self.assertEqual(files_sorted[1], (os.path.join(p, "batch04_QC17_rep02_263.mzML"), '2011-04-02T03:31:04Z'))
        self.assertEqual(files_sorted[2], (os.path.join(p, "batch04_QC17_rep03_264.mzML"), '2011-04-02T03:34:08Z'))

        ps.reverse()
        files_sorted = paths.sort_ms_files_by_timestamp(ps)
        self.assertEqual(files_sorted[0], (os.path.join(p, "batch04_QC17_rep01_262.mzML"), '2011-04-02T03:28:02Z'))
        self.assertEqual(files_sorted[1], (os.path.join(p, "batch04_QC17_rep02_263.mzML"), '2011-04-02T03:31:04Z'))
        self.assertEqual(files_sorted[2], (os.path.join(p, "batch04_QC17_rep03_264.mzML"), '2011-04-02T03:34:08Z'))

        p = to_test_data("MTBLS79_subset", "raw")
        ps = [os.path.join(p, fn) for fn in os.listdir(p)]
        files_sorted = paths.sort_ms_files_by_timestamp(ps)

        if platform.system() == "Darwin":
            self.assertEqual(files_sorted[0], (os.path.join(p, "batch04_QC17_rep01_262.RAW"), '02/04/2011 03:28:02'))
            self.assertEqual(files_sorted[1], (os.path.join(p, "batch04_QC17_rep02_263.RAW"), '02/04/2011 03:31:05'))
            self.assertEqual(files_sorted[2], (os.path.join(p, "batch04_QC17_rep03_264.RAW"), '02/04/2011 03:34:09'))
        else:
            self.assertEqual(files_sorted[0], (os.path.join(p, "batch04_QC17_rep01_262.RAW"), '4/2/2011 3:28:02 AM'))
            self.assertEqual(files_sorted[1], (os.path.join(p, "batch04_QC17_rep02_263.RAW"), '4/2/2011 3:31:05 AM'))
            self.assertEqual(files_sorted[2], (os.path.join(p, "batch04_QC17_rep03_264.RAW"), '4/2/2011 3:34:09 AM'))

        ps.reverse()
        files_sorted = paths.sort_ms_files_by_timestamp(ps)
        if platform.system() == "Darwin":
            self.assertEqual(files_sorted[0], (os.path.join(p, "batch04_QC17_rep01_262.RAW"), '02/04/2011 03:28:02'))
            self.assertEqual(files_sorted[1], (os.path.join(p, "batch04_QC17_rep02_263.RAW"), '02/04/2011 03:31:05'))
            self.assertEqual(files_sorted[2], (os.path.join(p, "batch04_QC17_rep03_264.RAW"), '02/04/2011 03:34:09'))
        else:
            self.assertEqual(files_sorted[0], (os.path.join(p, "batch04_QC17_rep01_262.RAW"), '4/2/2011 3:28:02 AM'))
            self.assertEqual(files_sorted[1], (os.path.join(p, "batch04_QC17_rep02_263.RAW"), '4/2/2011 3:31:05 AM'))
            self.assertEqual(files_sorted[2], (os.path.join(p, "batch04_QC17_rep03_264.RAW"), '4/2/2011 3:34:09 AM'))


if __name__ == '__main__':
    unittest.main()
