#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Albert
origin: 04-28-2017

"""


import unittest
import numpy as np
from dimspy.models.peaklist import PeakList
from dimspy.process.peak_alignment import align_peaks


class PeakAlignmentTestCase(unittest.TestCase):
    mz = [
        [10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        [10,     30,     50, 60, 70, 80, 90, 100],
        [20,     30, 40,         70, 80, 90, 100],
        [10, 20, 30,                 80, 90, 100],
        [10, 20,         50, 60, 70, 80,        ],
        [                50,                    ],
    ]

    ints = [
        [11, 12, 13, 14, 15, 16, 17, 18, 19, 110],
        [21,     23,     25, 26, 27, 28, 29, 210],
        [    32, 33, 34,         37, 38, 39, 310],
        [41, 42, 43,                 48, 49, 410],
        [51, 52,         55, 56, 57, 58,        ],
        [                65,                    ],
    ]

    strs = [
        ['a','b','c','d','e','f','g','h','i','j'],
        ['k',    'l',    'm','n','o','p','q','r'],
        [    's','t','u',        'v','w','x','y'],
        ['z','a','b',                'c','d','e'],
        ['f','g',        'h','i','j','k',       ],
        [                'l',                   ],
    ]

    def _createPeakLists(self):
        mz = [np.array(m) + np.random.normal(0, 1e-5, len(m)) for m in self.mz]
        pkls = []
        for i in range(len(mz)):
            pl = PeakList('peaklist_' + str(i), mz[i], self.ints[i])
            pl.add_attribute('str_attr', self.strs[i])
            pkls += [pl]
        return pkls

    def _checkAlignmentResults(self, pm):
        self.assertTrue(np.allclose(np.unique(np.round(pm.to_peaklist('merged').mz)), np.arange(10, 110, step = 10)))
        self.assertTrue(all(np.allclose(mi[mm != 0], ri) for mi, mm, ri in zip(pm.intensity_matrix, pm.mz_matrix, self.ints)))

    def test_normal_alignment(self):
        pkls = self._createPeakLists()

        try:
            pm = align_peaks(pkls, ppm = 2.0, block_size = 5, fixed_block = True, edge_extend = 10, ncpus = 2)
            # print pm.attr_matrix('str_attr')
            # print pm.attr_mean_vector('str_attr')
        except Exception, e:
            self.fail('alignment failed: ' + str(e))

        self._checkAlignmentResults(pm)

    def test_block_size(self):
        pkls = self._createPeakLists()
        try:
            pm = align_peaks(pkls, ppm = 2.0, block_size = 1, fixed_block = True, edge_extend = 10, ncpus = 2)
        except Exception, e:
            self.fail('alignment failed: ' + str(e))
        self._checkAlignmentResults(pm)

        pkls = self._createPeakLists()
        try:
            pm = align_peaks(pkls, ppm = 2.0, block_size = 20, fixed_block = True, edge_extend = 10, ncpus = 2)
        except Exception, e:
            self.fail('alignment failed: ' + str(e))
        self._checkAlignmentResults(pm)

    def test_ppm(self):
        pkls = self._createPeakLists()

        try:
            pm = align_peaks(pkls, ppm = 1e+10, block_size = 5, fixed_block = True, edge_extend = 10, ncpus = 2)
        except Exception, e:
            self.fail('alignment failed: ' + str(e))

        self.assertTrue(np.allclose(pm.to_peaklist('merged').mz, [np.mean(map(np.mean, self.mz))]))
        self.assertTrue(np.allclose(pm.intensity_matrix.flatten(), map(np.mean, self.ints)))
        self.assertTrue(np.allclose(pm.attr_matrix('intra_count').flatten(), map(len, self.mz)))

        try:
            pm = align_peaks(pkls, ppm = 1e-10, block_size = 5, fixed_block = True, edge_extend = 10, ncpus = 2)
        except Exception, e:
            self.fail('alignment failed: ' + str(e))

        self.assertTrue(np.allclose(pm.to_peaklist('merged').mz, np.sort(reduce(lambda x,y: x+y, map(list, self.mz)))))
        self.assertTrue(np.allclose(np.sort(np.sum(pm.intensity_matrix, axis = 0)), np.sort(reduce(lambda x,y: x+y, self.ints))))
        self.assertTrue(np.allclose(np.sum(pm.attr_matrix('intra_count'), axis = 0), np.ones(pm.shape[1])))

    def test_single_peaklist(self):
        pkls = [PeakList('peaklist_0', np.arange(10, 110, step = 10), np.arange(10) + 11)]

        try:
            pm = align_peaks(pkls, ppm = 2.0, block_size = 5, fixed_block = True, edge_extend = 10, ncpus = 2)
        except Exception, e:
            self.fail('alignment failed: ' + str(e))

        self.assertTrue(np.allclose(pm.to_peaklist('merged').mz, np.arange(10, 110, step = 10)))
        self.assertTrue(np.allclose(pm.intensity_matrix, [np.arange(10) + 11]))

    def test_special_peaklists(self):
        pkls = [PeakList('peaklist_' + str(i), np.ones(10) * 10, np.ones(10)) for i in range(6)]

        try:
            pm = align_peaks(pkls, ppm = 2.0, block_size = 5, fixed_block = False, edge_extend = 10, ncpus = 2)
        except Exception, e:
            self.fail('alignment failed: ' + str(e))

        self.assertTrue(np.allclose(pm.to_peaklist('merged').mz, [10.]))
        self.assertTrue(np.allclose(np.sum(pm.intensity_matrix, axis = 0), [6]))
        self.assertTrue(np.allclose(np.sum(pm.attr_matrix('intra_count'), axis = 0), [60]))

        try:
            pm = align_peaks(pkls, ppm = 1e-10, block_size = 1, fixed_block = True, edge_extend = 1, ncpus = 2)
        except Exception, e:
            self.fail('alignment failed: ' + str(e))

        self.assertTrue(np.allclose(pm.to_peaklist('merged').mz, [10.]))
        self.assertTrue(np.allclose(np.sum(pm.intensity_matrix, axis = 0), [6]))
        self.assertTrue(np.allclose(np.sum(pm.attr_matrix('intra_count'), axis = 0), [60]))

    # may take a while to run
    # def test_large_peaklists(self):
    #     pkls = [PeakList('peaklist_' + str(i),
    #                      np.sort(np.random.uniform(100, 1200, size = 10000)),
    #                      np.random.normal(100, 10, size = 10000))
    #             for i in range(100)]
    #
    #     try:
    #         pm = align_peaks(pkls, ppm = 2.0, block_size = 5000, fixed_block = False, edge_extend = 10, ncpus = 2)
    #     except Exception, e:
    #         self.fail('alignment failed: ' + str(e))


if __name__ == '__main__':
    unittest.main()
