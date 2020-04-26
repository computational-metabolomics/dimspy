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


import unittest

from dimspy.models.peaklist_tags import PeakList_Tags
from dimspy.process.peak_filters import *


class PeakFiltersTestCase(unittest.TestCase):
    @staticmethod
    def _createPeakList():
        pkl = PeakList('peaklist', np.arange(10, dtype = float), np.arange(10, dtype = float) + 1)
        pkl.add_attribute('snr', (np.arange(10, dtype = float) + 1) / 10)
        return pkl

    @staticmethod
    def _createPeakMatrix():
        pids, tags = list(zip(*[
            ('sample_1_1', PeakList_Tags('sample', treatment = 'compound_1', time_point = '1hr', plate = 1, order = 1)),
            ('sample_1_2', PeakList_Tags('sample', treatment = 'compound_1', time_point = '6hr', plate = 1, order = 2)),
            ('QC_1',       PeakList_Tags('qc', plate = 1, order = 3)),
            ('Blank_1',    PeakList_Tags('blank', plate = 1, order = 4)),
            ('sample_2_1', PeakList_Tags('sample', treatment = 'compound_2', time_point = '1hr', plate = 2, order = 1)),
            ('sample_2_2', PeakList_Tags('sample', treatment = 'compound_2', time_point = '6hr', plate = 2, order = 2)),
            ('QC_2',       PeakList_Tags('qc', plate = 2, order = 3)),
            ('Blank_2',    PeakList_Tags('blank', plate = 2, order = 4)),
        ]))

        mzs = np.tile(np.arange(0, 1000, step = 100, dtype = float), (8, 1))
        ints = np.arange(80, dtype = float).reshape((8, 10)) / 20.
        ints[3, 1] = ints[7, 1] = ints[7, 3] = 0 # test blank filter
        ics = np.array([[1, 2] * 5] * 8)

        return PeakMatrix(pids, tags, (('mz', mzs), ('intensity', ints), ('intra_count', ics)))

    # peaklist filters
    def test_peaklist_attr_filter(self):
        pkl = self._createPeakList()

        try:
            filter_attr(pkl, 'snr', 0.5, flag_index = 2)
        except Exception as e:
            self.fail('filter peaklist attribute failed: ' + str(e))
        self.assertListEqual(pkl.snr.tolist(), [0.1, 0.2, 0.3, 0.4, 0.5])
        self.assertTupleEqual(pkl.attributes, ('mz', 'intensity', 'snr_flag', 'snr'))

        self.assertRaises(AttributeError, lambda: filter_attr(pkl, 'not_exists', 0.5))
        self.assertRaises(AttributeError, lambda: filter_attr(pkl, 'snr', 0.6))
        self.assertRaises(ValueError, lambda: filter_attr(pkl, 'snr'))

        filter_attr(pkl, 'snr', min_threshold = 0.4, max_threshold = 0.4, flag_name = 'new_snr_flag')
        self.assertListEqual(pkl.mz.tolist(), [3])

    def test_peaklist_ringing_filter(self):
        pkl = self._createPeakList()

        try:
            filter_ringing(pkl, threshold = 0.9, bin_size = 3.0)
        except Exception as e:
            self.fail('filter peaklist ringing failed: ' + str(e))
        self.assertListEqual(pkl.mz.tolist(), [2., 5., 8., 9.])

    def test_peaklist_mz_ranges(self):
        pkl = self._createPeakList()

        try:
            filter_mz_ranges(pkl, [(1.,3.), (5.,8.)])
        except Exception as e:
            self.fail('filter peaklist mz ranges failed: ' + str(e))
        self.assertListEqual(pkl.mz.tolist(), [0., 4., 9.])

    # peakmatrix filters
    def test_peak_matrix_rsd_filter(self):
        pm = self._createPeakMatrix()

        try:
            pm = filter_rsd(pm, 62, 'qc')
        except Exception as e:
            self.fail('filter peak_matrix rsd failed: ' + str(e))
        self.assertTrue(np.allclose(pm.rsd('qc'),
            [61.48754619, 60.17930052, 58.92556509, 57.72300254]))

        self.assertRaises(AttributeError, lambda: filter_rsd(pm, 45, 'not_QC'))

    def test_peak_matrix_fraction_filter(self):
        pm = self._createPeakMatrix()
        for attr in ('mz', 'intensity', 'intra_count'): pm._attr_dict[attr][:,1] = 0

        try:
            pm = filter_fraction(pm, 1)
        except Exception as e:
            self.fail('filter peak_matrix fraction failed: ' + str(e))
        self.assertEqual(pm.shape[1], 9)

        pm = self._createPeakMatrix()
        for attr in ('mz', 'intensity', 'intra_count'):
            pm._attr_dict[attr][:,1] *= [1, 1, 1, 0, 1, 1, 1, 0]
            pm._attr_dict[attr][:,2] *= [1, 1, 1, 1, 1, 1, 0, 0]

        pm = filter_fraction(pm, 0.6, within_classes = True, class_tag_type = 'plate')
        self.assertEqual(pm.shape[1], 10)
        self.assertRaises(AttributeError, lambda: filter_fraction(pm, 1, within_classes = True, class_tag_type = 'time_point'))
        self.assertRaises(KeyError, lambda: filter_fraction(pm, 1, within_classes = True))

    def test_peak_matrix_blank_filter(self):
        pm = self._createPeakMatrix()
        pm = filter_blank_peaks(pm, 'blank', 0.3)
        self.assertTupleEqual(pm.shape, (6, 10))

        pm = self._createPeakMatrix()
        pm = filter_blank_peaks(pm, 'blank', 0.4)
        self.assertTupleEqual(pm.shape, (6, 2))

        pm = self._createPeakMatrix()
        pm = filter_blank_peaks(pm, 'blank', 0.3, method = 'max')
        self.assertTupleEqual(pm.shape, (6, 2))

        pm = self._createPeakMatrix()
        pm = filter_blank_peaks(pm, 'blank', 0.3, fold_threshold = 2)
        self.assertTupleEqual(pm.shape, (6, 1))

        pm = self._createPeakMatrix()
        self.assertRaises(ValueError, lambda: filter_blank_peaks(pm, 'Not_blank', 0.3))


if __name__ == '__main__':
    unittest.main()
