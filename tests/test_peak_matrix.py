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
from dimspy.models.peak_matrix import PeakMatrix, peak_matrix_rsd
from dimspy.models.peak_matrix import mask_peakmatrix, unmask_peakmatrix, mask_all_peakmatrix, unmask_all_peakmatrix


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
        ics = np.array([[2] * 10] * 6)
        # simulate missing values
        for m in (mzs, ints, ics):
            np.fill_diagonal(m, 0)
            m[:,2] = 0
        return PeakMatrix(pids, tags, [('mz', mzs), ('intensity', ints), ('intra_count', ics)])

    def test_properties(self):
        pm = self._createPeakMatrix()

        pm.mask = [True, False] * 3
        self.assertTrue(np.all(pm.mask == [True, False, True, False, True, False]))
        pm.mask = None
        self.assertTrue(np.all(pm.mask == [False] * 6))

        self.assertTupleEqual(pm.flag_names, ())
        self.assertTrue(np.all(pm.flags == np.ones(10)))

        self.assertTupleEqual(pm.attributes, ('mz', 'intensity', 'intra_count'))

        self.assertTupleEqual(pm.peaklist_ids,
            ('sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'sample_2_2', 'QC_2'))

        self.assertEqual(len(pm.peaklist_tags), 6)
        self.assertListEqual(sorted(pm.peaklist_tag_types),
                             sorted(('treatment', 'time_point', 'plate')))
        self.assertListEqual(sorted(pm.peaklist_tag_values),
                             sorted(('sample', 'qc', 'compound_1', 'compound_2', '1hr', '6hr', 1, 2)))

        pm.mask = [True, False] * 3
        self.assertTupleEqual(pm.shape, (3, 10))
        self.assertTupleEqual(pm.full_shape, (6, 10))
        pm.mask = None

        self.assertTrue(np.all(pm.present == [5]*2+[0]+[5]*3+[6]*4))
        self.assertTrue(np.allclose(pm.fraction, [0.83333333]*2+[0]+[0.83333333]*3+[1]*4))
        self.assertTrue(np.all(pm.missing_values == [2]*2+[1]+[2]*3))
        self.assertTrue(np.all(pm.occurrence == [10]*2+[0]+[10]*3+[12]*4))
        self.assertTrue(np.allclose(pm.purity[~np.isnan(pm.purity)], [0]*9))

        pm.add_flag('odd_flag', [1, 0] * 5)
        self.assertTrue(np.all(pm.property('present') == [5, 0, 5, 6, 6]))
        self.assertTrue(np.all(pm.property('present', flagged_only = False) == [5]*2+[0]+[5]*3+[6]*4))
        pm.drop_flag('odd_flag')

        mmz = np.arange(0, 1000, step = 100, dtype = float) + 1
        mmz[2] = np.nan
        self.assertTrue(np.allclose(*map(np.nan_to_num, (pm.mz_mean_vector, mmz))))
        mit = [30., 29., np.nan, 27., 26., 25., 31., 32., 33., 34.]
        self.assertTrue(np.allclose(*map(np.nan_to_num, (pm.intensity_mean_vector*20, mit))))

    def test_mask(self):
        pm = self._createPeakMatrix()

        self.assertListEqual(sorted(pm.tags_of('plate')), [1, 2])
        self.assertListEqual(sorted(pm.tags_of()), sorted(('sample', 'qc')))
        self.assertRaises(ValueError, lambda: pm.tags_of('treatment'))
        self.assertRaises(ValueError, lambda: pm.tags_of('not_exist'))

        pm.mask_tags('qc', plate = 1) # mask samples with both of the two
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2', 'QC_2'))
        pm.mask = None
        pm.mask_tags('qc')
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2'))
        pm.mask = None
        pm.mask_tags('qc').mask_tags(plate = 1)
        self.assertTupleEqual(pm.peaklist_ids, ('sample_2_1', 'sample_2_2'))
        pm.mask = None
        pm.mask_tags('not_exist')
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'sample_2_2', 'QC_2'))

        pm.mask = [True] * 6
        pm.unmask_tags('qc', plate = 1) # unmask samples with both of the two
        self.assertTupleEqual(pm.peaklist_ids, ('QC_1',))
        pm.mask = [True] * 6
        pm.unmask_tags('qc')
        self.assertTupleEqual(pm.peaklist_ids, ('QC_1', 'QC_2'))
        pm.mask = [True] * 6
        pm.unmask_tags('qc').unmask_tags(plate = 1)
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'QC_1', 'QC_2'))
        pm.mask = [True] * 6
        pm.unmask_tags('not_exist')
        self.assertTupleEqual(pm.peaklist_ids, ())

        pm.unmask_tags('qc', override = True)
        self.assertTupleEqual(pm.peaklist_ids, ('QC_1', 'QC_2'))
        with mask_all_peakmatrix(pm) as m:
            m.unmask_tags('sample')
            self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2'))

        pm.mask_tags('qc', override = True)
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2'))
        with unmask_all_peakmatrix(pm) as m:
            m.mask_tags('sample')
            self.assertTupleEqual(m.peaklist_ids, ('QC_1', 'QC_2'))
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2'))

        pm.mask = None
        with unmask_peakmatrix(pm, plate = 1) as m:
            self.assertTupleEqual(m.peaklist_ids, ('sample_1_1', 'sample_1_2', 'QC_1'))
            self.assertTupleEqual(m.full_shape, (6, 10))
        self.assertEqual(len(pm.peaklist_ids), 6)

        with mask_peakmatrix(pm, plate = 2) as m:
            self.assertTupleEqual(m.peaklist_ids, ('sample_1_1', 'sample_1_2', 'QC_1'))
            with unmask_all_peakmatrix(pm) as mm:
                self.assertTupleEqual(mm.peaklist_ids,
                                      ('sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'sample_2_2', 'QC_2'))

        with mask_peakmatrix(pm, 'qc') as m:
            m.remove_samples((1, 2))
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'QC_1', 'sample_2_2', 'QC_2'))

    def test_flags(self):
        pm = self._createPeakMatrix()

        self.assertTrue(np.sum(pm.flags) == 10)

        pm.add_flag('qua_flag', [1, 1, 0, 1] * 2 + [1, 1], flagged_only = True)
        pm.add_flag('odd_flag', [1, 0] * 5, flagged_only = False)

        self.assertTupleEqual(pm.flag_names, ('qua_flag', 'odd_flag'))
        self.assertTrue(np.all(pm.flags == [1, 0, 0, 0, 1, 0, 0, 0, 1, 0]))
        self.assertTrue(np.all(pm.flag_values('odd_flag') == [1, 0] * 5))

        with mask_peakmatrix(pm, 'qc') as m:
            self.assertTupleEqual(m.shape, (4, 3))
        self.assertTupleEqual(m.shape, (6, 3))
        self.assertTupleEqual(m.full_shape, (6, 10))

        with mask_peakmatrix(pm, plate = 1) as m:
            mzs = np.array([
                [   1.,    0.,  401.,  601.,  801.],
                [   1.,    0.,    0.,  601.,  801.],
                [   1.,    0.,  401.,  601.,  801.],
            ])
            m.drop_flag('qua_flag')
            self.assertTrue(np.allclose(m.mz_matrix, mzs))
        self.assertTupleEqual(pm.shape, (6, 5))

    def test_access(self):
        pm = self._createPeakMatrix()

        pm.add_flag('even_flag', [0, 1] * 5)
        self.assertTrue(np.allclose(pm.attr_mean_vector('mz'),
                                    [101.0, 301.0, 501.0, 701.0, 901.0]))
        self.assertTrue(np.allclose(*map(np.nan_to_num, (pm.attr_mean_vector('mz', flagged_only = False),
                                    [1.0, 101.0, np.nan, 301.0, 401.0, 501.0, 601.0, 701.0, 801.0, 901.0]))))
        self.assertTrue(np.allclose((lambda x: x[~np.isnan(x)])(peak_matrix_rsd(pm, 'qc')),
                                    [41.666667, 39.473684, 35.714285, 34.090909]))
        self.assertTrue(np.allclose((lambda x: x[~np.isnan(x)])(peak_matrix_rsd(pm)),
                                    [59.32638115, 68.6934703, 56.56854249, 53.36953524, 50.23015081]))

        pm.remove_peaks((0, 1), flagged_only = False)
        self.assertTrue(np.allclose((lambda x: x[~np.isnan(x)])(peak_matrix_rsd(pm, 'qc')),
                                    [39.473684, 35.714285, 34.090909]))
        pm.remove_peaks((0, 1), flagged_only = True)
        self.assertTrue(np.allclose(peak_matrix_rsd(pm, 'qc'),
                                    [35.714285, 34.090909]))

        self.assertRaises(AttributeError, lambda: peak_matrix_rsd(pm, 'no_such_tag'))

        with mask_peakmatrix(pm, 'sample', plate = 1):
            pm.remove_samples((0, 1))
            self.assertTupleEqual(pm.peaklist_ids, ('sample_2_2', 'QC_2'))
            pm.remove_samples((1, 2), masked_only = False)
            self.assertTupleEqual(pm.peaklist_ids, ('QC_2',))
        self.assertTupleEqual(pm.peaklist_ids, ('sample_1_1', 'QC_2'))

    def test_exports(self):
        pm = self._createPeakMatrix()

        pm.add_flag('even_flag', [0, 1] * 5)
        with mask_peakmatrix(pm, plate = 1):
            peaklists = pm.extract_peaklists()
        self.assertListEqual(map(lambda x: x.ID, peaklists), ['sample_2_1', 'sample_2_2', 'QC_2'])

        mzs = [
            [101.0, 501.0, 701.0, 901.0],
            [101.0, 301.0, 501.0, 701.0, 901.0],
            [101.0, 301.0, 701.0, 901.0],
        ]
        self.assertTrue(all(map(lambda x: np.allclose(x[0].mz, x[1]), zip(peaklists, mzs))))

        pm.drop_flag('even_flag')
        pkl = pm.to_peaklist('merged_pkl')
        self.assertTrue(np.allclose(pkl.mz, [1.0, 101.0, 301.0, 401.0, 501.0, 601.0, 701.0, 801.0, 901.0]))

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
        self.assertTupleEqual(pm.attributes, ('mz', 'intensity', 'intra_count'))


if __name__ == '__main__':
    unittest.main()
