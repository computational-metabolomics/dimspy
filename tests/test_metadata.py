#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_metadata

author(s): Albert, Ralf, Tom
origin: 29-10-2019

"""

import unittest
import os

from dimspy.metadata import validate_metadata


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset", *args)


class ValidateMetadataTestCase(unittest.TestCase):

    def test_filelist_standard(self):
        # filename	replicate	batch	injectionOrder	classLabel
        fm_dict = validate_metadata(to_test_data("filelist_csl_MTBLS79_mzml_triplicates.txt"))
        self.assertEqual(fm_dict['filename'], ['batch04_B02_rep01_301.mzML', 'batch04_B02_rep02_302.mzML',
                                               'batch04_B02_rep03_303.mzML', 'batch04_QC17_rep01_262.mzML',
                                               'batch04_QC17_rep02_263.mzML', 'batch04_QC17_rep03_264.mzML',
                                               'batch04_S01_rep01_247.mzML', 'batch04_S01_rep02_248.mzML',
                                               'batch04_S01_rep03_249.mzML'])
        self.assertEqual(fm_dict['replicate'], [1, 2, 3, 1, 2, 3, 1, 2, 3])
        self.assertEqual(fm_dict['batch'], [1] * 9)
        self.assertEqual(fm_dict['injectionOrder'], [1, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(fm_dict['classLabel'], ['blank', 'blank', 'blank',
                                                 'QC', 'QC', 'QC', 'sample', 'sample', 'sample'])

    def test_filelist_multi(self):
        fm_dict = validate_metadata(to_test_data("filelist_multi.txt"))
        self.assertEqual(fm_dict['multilist'], [1, 1, 2])

    def test_filename_error(self):
        with self.assertRaises(Exception) as context:
            validate_metadata(to_test_data("filelist_filename_error.txt"))
        self.assertTrue("Duplicate filename in list" in str(context.exception))

    def test_filelist_multilist_error(self):
        with self.assertRaises(Exception) as context:
            validate_metadata(to_test_data("filelist_multi_error.txt"))
        self.assertTrue("Column 'multilist' values should be integers" in str(context.exception))

    def test_filelist_injection_order_error(self):
        with self.assertRaises(Exception) as context:
            validate_metadata(to_test_data("filelist_injection_order_error.txt"))
        self.assertTrue("samples not in order" in str(context.exception))

    def test_filelist_class_label_error(self):
        with self.assertRaises(Exception) as context:
            validate_metadata(to_test_data("filelist_class_label_error.txt"))
        self.assertTrue("class names do not match with number of replicates" in str(context.exception))

    def test_filelist_replicate_error_zero_value(self):
        with self.assertRaises(Exception) as context:
            validate_metadata(to_test_data("filelist_replicate_error_1.txt"))
        self.assertTrue("Incorrect replicate number in list" in str(context.exception))

    def test_filelist_replicate_error_zero_value(self):
        with self.assertRaises(Exception) as context:
            validate_metadata(to_test_data("filelist_replicate_error_2.txt"))
        self.assertTrue("Incorrect numbering for replicates" in str(context.exception))



if __name__ == '__main__':
    unittest.main()
