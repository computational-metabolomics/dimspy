#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_peak_filters

author(s): Albert
origin: 05-14-2017

"""


import unittest
import numpy as np
from dimspy.models.peaklist_tags import PeakList_Tags
from dimspy.models.peaklist import PeakList
from dimspy.models.peak_matrix import PeakMatrix, unmask_peakmatrix
from dimspy.process.peak_filters import filter_attr, filter_rsd, filter_fraction, filter_blank_peaks, filter_sparsity


class PeakFiltersTestCase(unittest.TestCase):
    def test_peaklist_attr_filter(self):
        pkl = PeakList('peaklist', np.arange(10, dtype = float), np.ones(10, dtype = float))
        pkl.add_attribute('snr', (np.arange(10, dtype = float) + 1) / 10)

        try:
            filter_attr(pkl, 'snr', 0.5, flag_index = 2)
        except Exception, e:
            self.fail('filter peaklist attribute failed: ' + str(e))
        self.assertListEqual(pkl['snr'].tolist(), [0.1, 0.2, 0.3, 0.4, 0.5])
        self.assertTupleEqual(pkl.attributes, ('mz', 'intensity', 'snr_flag', 'snr'))

        self.assertRaises(AttributeError, lambda: filter_attr(pkl, 'not_exists', 0.5))
        self.assertRaises(AttributeError, lambda: filter_attr(pkl, 'snr', 0.6))
        self.assertRaises(ValueError, lambda: filter_attr(pkl, 'snr'))

        filter_attr(pkl, 'snr', min_threshold = 0.4, max_threshold = 0.4, flag_name = 'new_snr_flag')
        self.assertListEqual(pkl.mz.tolist(), [3])

    @staticmethod
    def _createPeakMatrix():
        pids, tags = zip(*[
            ('sample_1_1', PeakList_Tags('sample', treatment = 'compound_1', time_point = '1hr', plate = 1)),
            ('sample_1_2', PeakList_Tags('sample', treatment = 'compound_1', time_point = '6hr', plate = 1)),
            ('QC_1',       PeakList_Tags('qc', plate = 1)),
            ('Blank_1',    PeakList_Tags('blank', plate = 1)),
            ('sample_2_1', PeakList_Tags('sample', treatment = 'compound_2', time_point = '1hr', plate = 2)),
            ('sample_2_2', PeakList_Tags('sample', treatment = 'compound_2', time_point = '6hr', plate = 2)),
            ('QC_2',       PeakList_Tags('qc', plate = 2)),
            ('Blank_2',    PeakList_Tags('blank', plate = 2)),
        ])

        mzs = np.tile(np.arange(0, 1000, step = 100, dtype = float) + 1, (8, 1))
        ints = np.arange(80, dtype = float).reshape((8, 10)) / 20.
        ics = np.array([[1, 2] * 5] * 8)

        return PeakMatrix(pids, tags, mz = mzs, intensity = ints, intra_count = ics)

    def test_peak_matrix_rsd_filter(self):
        pm = self._createPeakMatrix()
        pm = filter_rsd(pm, 60)
        self.assertTrue(np.allclose(pm.rsd,
            [50., 58.75097045, 57.28219619, 55.88506945, 54.55447256, 53.28576389, 52.07472381]))

        pm = self._createPeakMatrix()
        pm = filter_rsd(pm, qc_label = 'qc')
        self.assertTrue(pm.is_empty())

        pm = self._createPeakMatrix()
        self.assertRaises(AttributeError, lambda: filter_rsd(pm, qc_label = 'not_QC'))
        def _maskedcall():
            with unmask_peakmatrix('qc'): filter_rsd(pm, qc_label = 'qc')
        self.assertRaises(AttributeError,  _maskedcall)
        self.assertRaises(ValueError, lambda: filter_rsd(pm))

    def test_peak_matrix_fraction_filter(self):
        pm = self._createPeakMatrix()
        pm = filter_fraction(pm, 1)
        self.assertEqual(pm.shape[1], 9)

        pm = self._createPeakMatrix()
        pm = filter_fraction(pm, 1, within_classes = True, class_tag_type = 'plate')
        self.assertEqual(pm.shape[1], 9)
        self.assertRaises(AttributeError, lambda: filter_fraction(pm, 1, within_classes = True, class_tag_type = 'time_point'))
        self.assertRaises(AttributeError, lambda: filter_fraction(pm, 1, within_classes = True))

        pm = self._createPeakMatrix()
        pm = filter_fraction(pm, 2)
        self.assertTrue(pm.is_empty())

    def test_peak_matrix_blank_filter(self):
        pm = self._createPeakMatrix()
        pm = filter_blank_peaks(pm, 'blank', 0.3)
        self.assertTupleEqual(pm.shape, (6, 10))

        pm = self._createPeakMatrix()
        pm = filter_blank_peaks(pm, 'blank', 0.4)
        self.assertTrue(pm.is_empty())

        pm = self._createPeakMatrix()
        pm = filter_blank_peaks(pm, 'blank', 0.3, method = 'max')
        self.assertTrue(pm.is_empty())

        pm = self._createPeakMatrix()
        self.assertRaises(ValueError, lambda: filter_blank_peaks(pm, 'Not_blank', 0.3))

    def test_peak_matrix_sparsity_filter(self):
        pm = self._createPeakMatrix()
        pm = filter_sparsity(pm, 1e+1)
        self.assertTupleEqual(pm.shape, (8, 10))
        pm = filter_sparsity(pm, 5e+5)
        self.assertTupleEqual(pm.shape, (8, 2))


if __name__ == '__main__':
    unittest.main()
