#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_txt_portal

author(s): Albert
origin: 05-14-2017

"""


import unittest, os
import numpy as np
from dimspy.models.peaklist import PeakList
from dimspy.process.peak_alignment import align_peaks
from dimspy.portals.txt_portal import save_peaklist_as_txt, load_peaklist_from_txt
from dimspy.portals.txt_portal import save_peak_matrix_as_txt, load_peak_matrix_from_txt


class TxtPortalsTestCase(unittest.TestCase):
    def test_peaklist_portal(self):
        pkl = PeakList('peaklist', np.sort(np.random.uniform(100, 1200, size = 100)), np.random.normal(100, 10, size = 100))
        pkl.add_attribute('odd_flag', [0, 1] * 50, is_flag = True)

        save_peaklist_as_txt(pkl, 'peaklist.txt')
        npkl = load_peaklist_from_txt('peaklist.txt', 'peaklist')

        self.assertEqual(npkl.size, 50)
        self.assertEqual(npkl.full_size, 100)
        self.assertTrue(np.allclose(pkl.mz_all, npkl.mz_all))
        self.assertTrue(np.allclose(pkl.intensity, npkl.intensity))

    def test_peak_matrix_portal(self):
        _mzs = lambda: sorted(np.random.uniform(100, 1200, size = 100))
        _ints = lambda: np.abs(np.random.normal(100, 10, size = 100))

        pkls = [
            PeakList('sample_1_1', _mzs(), _ints()),
            PeakList('sample_1_2', _mzs(), _ints()),
            PeakList('QC_1', _mzs(), _ints()),
            PeakList('sample_2_1', _mzs(), _ints()),
            PeakList('sample_2_2', _mzs(), _ints()),
            PeakList('QC_2', _mzs(), _ints()),
        ]
        pkls[0].tags.add_tags('sample', treatment = 'compound_1', time_point = '1hr', plate = 1)
        pkls[1].tags.add_tags('sample', treatment = 'compound_1', time_point = '6hr', plate = 1)
        pkls[2].tags.add_tags('qc', plate = 1)
        pkls[3].tags.add_tags('sample', treatment = 'compound_2', time_point = '1hr', plate = 2)
        pkls[4].tags.add_tags('sample', treatment = 'compound_2', time_point = '6hr', plate = 2)
        pkls[5].tags.add_tags('qc', plate = 2)

        pm = align_peaks(pkls, ppm = 2.0, block_size = 10, ncpus = 2)

        save_peak_matrix_as_txt(pm, 'peak_matrix.txt', transpose = True, extend = True)
        npm = load_peak_matrix_from_txt('peak_matrix.txt', transposed = True, extended = True)

        self.assertEqual(pm.shape, npm.shape)
        self.assertTrue(np.allclose(pm.intensity_matrix, npm.intensity_matrix))
        self.assertTupleEqual(pm.peaklist_tag_types, npm.peaklist_tag_types)
        self.assertTupleEqual(pm.peaklist_tag_values, npm.peaklist_tag_values)

    def tearDown(self):
        if os.path.isfile('peaklist.txt'): os.remove('peaklist.txt')
        if os.path.isfile('peak_matrix.txt'): os.remove('peak_matrix.txt')


if __name__ == '__main__':
    unittest.main()
