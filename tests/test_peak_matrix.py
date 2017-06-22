#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_peak_matrix: 

author(s): Albert
origin: 05-10-2017

"""


import unittest
import numpy as np
import cPickle as cp
from dimspy.models.peaklist_tags import PeakList_Tags
from dimspy.models.peak_matrix import PeakMatrix, mask_peakmatrix, unmask_peakmatrix, unmask_all_peakmatrix


class PeakListTestCase(unittest.TestCase):
    @staticmethod
    def _createPeakMatrix():
        pids, tags = zip(*[
            ('sample_1_1', PeakList_Tags('sample', treatment = 'compound_1', time_point = '1hr', plate = 1)),
            ('sample_1_2', PeakList_Tags('sample', treatment = 'compound_1', time_point = '6hr', plate = 1)),
            ('QC_1',       PeakList_Tags('qc', plate = 1)),
            ('sample_2_1', PeakList_Tags('sample', treatment = 'compound_2', time_point = '1hr', plate = 2)),
            ('sample_2_2', PeakList_Tags('sample', treatment = 'compound_2', time_point = '6hr', plate = 2)),
            ('QC_2',       PeakList_Tags('qc', plate = 2)),
        ])

        mzs = np.tile(np.arange(0, 1000, step = 100, dtype = float) + 1, (6, 1))
        ints = np.arange(60, dtype = float).reshape((6, 10)) / 20.
        ics = np.array([[0, 2] * 5] * 6)
        map(lambda x: np.fill_diagonal(x, 0), (mzs, ints, ics)) # simulate missing values

        return PeakMatrix(pids, tags, mz = mzs, intensity = ints, intra_count = ics)

    def test_properties(self):
        pm = self._createPeakMatrix()

        pm.mask = [True, False] * 3
        self.assertTupleEqual(pm.shape, (3, 10))
        self.assertTupleEqual(pm.full_shape, (6, 10))
        self.assertTrue(np.all(pm.mask == [True, False, True, False, True, False]), str(pm.mask))
        pm.mask = None
        self.assertTrue(np.all(pm.mask == [True] * 6))

        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'sample_2_2', 'QC_2'))
        self.assertEqual(len(pm.peaklist_tags), 6)

        self.assertListEqual(sorted(pm.peaklist_tag_types), sorted(('treatment', 'time_point', 'plate')))
        self.assertListEqual(sorted(pm.peaklist_tag_values),
                             sorted(('sample', 'qc', 'compound_1', 'compound_2', '1hr', '6hr', 1, 2)))

        self.assertListEqual(sorted(pm.attributes), sorted(('mz', 'intensity', 'intra_count')))

        self.assertTrue(np.all(pm.present == [5] * 6 + [6] * 4))
        self.assertTrue(np.allclose(pm.fraction, [0.83333333] * 6 + [1] * 4))
        self.assertTrue(np.all(pm.missing_values == [1] * 6))
        self.assertTrue(np.allclose(pm.rsd, [47.14045208, 59.32638115, 66.24013211, 68.69347034, 66.17173282,
                                             56.56854249, 55.09113315, 53.36953524, 51.7522766 , 50.23015081,]))
        self.assertTrue(np.all(pm.occurrence == [0, 10] * 3 + [0, 12] * 2))
        self.assertTrue(np.allclose(pm.impure, [0, 0.83333333] * 3 + [0, 1] * 2))

        self.assertTrue(np.allclose(pm.mz_mean_vector, np.arange(0, 1000, step = 100)+1))
        self.assertTrue(np.allclose(pm.ints_mean_vector * 20, [30, 29, 28, 27, 26, 25, 31, 32, 33, 34]))

    def test_mask(self):
        pm = self._createPeakMatrix()

        self.assertListEqual(sorted(pm.tags_of('plate')), [1, 2])
        self.assertListEqual(sorted(pm.tags_of()), sorted(('sample', 'qc')))
        self.assertRaises(ValueError, lambda: pm.tags_of('treatment'))
        self.assertRaises(ValueError, lambda: pm.tags_of('not_exist'))

        pm.mask_tags('qc', plate = 1) # retain samples with any if the two
        self.assertListEqual(sorted(pm.peaklist_ids), sorted(('sample_1_1', 'sample_1_2', 'QC_1', 'QC_2')))
        pm.mask = None
        pm.mask_tags('qc').mask_tags(plate = 1) # retain samples with both if the two
        self.assertListEqual(sorted(pm.peaklist_ids), sorted(('QC_1',)))
        pm.mask = None
        pm.mask_tags('not_exist')
        self.assertListEqual(sorted(pm.peaklist_ids), [])

        pm.mask = None
        pm.unmask_tags('qc', plate = 1) # remove samples with both if the two
        self.assertListEqual(sorted(pm.peaklist_ids),
                             sorted(('sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2', 'QC_2')))
        pm.mask = None
        pm.unmask_tags('qc').unmask_tags(plate = 1) # remove samples with any if the two
        self.assertListEqual(sorted(pm.peaklist_ids),
                             sorted(('sample_2_1', 'sample_2_2')))
        pm.mask = None
        pm.unmask_tags(treatment = 'not_exist')
        self.assertListEqual(sorted(pm.peaklist_ids),
                             sorted(('sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'sample_2_2', 'QC_2')))

        self.mask = None
        with mask_peakmatrix(pm, plate = 1) as m:
            self.assertListEqual(sorted(m.peaklist_ids),
                                 sorted(('sample_1_1', 'sample_1_2', 'QC_1')))
        self.assertEqual(np.sum(pm.mask), 6)

        with unmask_peakmatrix(pm, plate = 1) as m:
            self.assertListEqual(sorted(m.peaklist_ids),
                                 sorted(('sample_2_1', 'sample_2_2', 'QC_2')))
            with unmask_all_peakmatrix(pm) as mm:
                self.assertListEqual(sorted(mm.peaklist_ids),
                                     sorted(('sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'sample_2_2', 'QC_2')))

    def test_access(self):
        pm = self._createPeakMatrix()

        pkl = pm.to_peaklist('merged_pkl')
        self.assertTrue(np.allclose(pkl.rsd, [47.14045208, 59.32638115, 66.24013211, 68.69347034, 66.17173282,
                                              56.56854249, 55.09113315, 53.36953524, 51.7522766 , 50.23015081,]))

        def _with_remove():
            with unmask_peakmatrix(pm, 'qc') as m: m.remove_samples((0, 1))
        self.assertRaises(ReferenceError, _with_remove)

        pm.unmask_tags('qc')
        self.assertListEqual(sorted(pm.peaklist_ids), sorted(('sample_2_1', 'sample_2_2')))

        pm._attr_dict['intensity'][0][2:] = 0 # QC_1
        pm.remove_peaks((0, 1))
        self.assertListEqual(sorted(pm.peaklist_ids), sorted(('sample_2_1', 'sample_2_2')))
        self.assertEqual(pm.shape, (2, 8))
        pm.mask = None
        self.assertListEqual(sorted(pm.peaklist_ids), sorted(('sample_2_1', 'sample_2_2', 'QC_2')))
        self.assertEqual(pm.shape, (3, 8))

        self.assertFalse(pm.is_empty())

        try:
            pm.to_str(comprehensive = True)
        except Exception, e:
            self.fail('PeakMatrix to_str() method failed: ' + str(e))

    def test_pickle(self):
        pm = self._createPeakMatrix()
        try:
            pstr = cp.dumps(pm)
            pm = cp.loads(pstr)
        except Exception, e:
            self.fail('PeakMatrix pickle failed: ' + str(e))
        self.assertListEqual(sorted(pm.attributes), sorted(('mz', 'intensity', 'intra_count')))

    def test_get_peaklist(self):
        pm = self._createPeakMatrix()
        peaklists = pm.get_peaklists()

        self.assertEqual(len(peaklists), 6)

        mzs = np.array([[101., 201., 301., 401., 501., 601., 701., 801., 901.],
               [1., 201., 301., 401., 501., 601., 701., 801., 901.],
               [1., 101., 301., 401., 501., 601., 701., 801., 901.],
               [1., 101., 201., 401., 501., 601., 701., 801., 901.],
               [1., 101., 201., 301., 501., 601., 701., 801., 901.],
               [1., 101., 201., 301., 401., 601., 701., 801., 901.]])

        ints = np.array([[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45],
               [0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95],
               [1., 1.05, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45],
               [1.5, 1.55, 1.6, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95],
               [2., 2.05, 2.1, 2.15, 2.25, 2.3, 2.35, 2.4, 2.45],
               [2.5, 2.55, 2.6, 2.65, 2.7, 2.8, 2.85, 2.9, 2.95]])

        mz_test = np.array([p.to_list()[0] for p in peaklists])
        ints_test = np.array([p.to_list()[1] for p in peaklists])

        self.assertTrue(np.alltrue(mz_test == mzs))
        self.assertTrue(np.alltrue(ints_test == ints))



if __name__ == '__main__':
    unittest.main()
