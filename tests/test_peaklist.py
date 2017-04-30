#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Albert
origin: 04-28-2017

"""


import unittest
import numpy as np
import cPickle as cp
from dimspy.models.peaklist import PeakList


class PeakListTestCase(unittest.TestCase):
    @staticmethod
    def _createPeakList():
        mzs = sorted(np.random.uniform(100, 1000, size = 10))
        ints = np.abs(np.random.normal(10, 3, size = 10))
        pl = PeakList('sample_peaklist', mzs, ints, mz_range = (100, 1000), frag_mode = 'slb')

        # snr = np.abs(np.random.normal(1000, 400, size = 10))
        # pl.add_attribute('snr', snr)
        # pl.add_tags('sample', 'passed_qc', treatment = 'high_dose')
        # pl.metadata.type = 'blank'

        return pl

    def test_creation(self):
        try:
            self._createPeakList()
        except Exception, e:
            self.fail('create PeakList object failed: ' + str(e))

    def test_properties(self):
        pl = self._createPeakList()
        self.assertEqual(pl.ID, 'sample_peaklist')

        pl.add_attribute('odd_flag', [1, 0] * 5, is_flag = True)
        self.assertEqual(pl.size, 5)
        self.assertEqual(pl.full_size, 10)
        self.assertTupleEqual(pl.shape, (5, 3))
        self.assertTupleEqual(pl.full_shape, (10, 3))

        try:
            pl.metadata.type = 'blank'
        except Exception, e:
            self.fail('access metadata failed: ' + str(e))
        self.assertListEqual(sorted(pl.metadata.keys()), ['frag_mode', 'mz_range', 'type'])

        try:
            pl.tags.add_tags('sample', 'passed_qc', treatment = 'high_dose')
        except Exception, e:
            self.fail('access tags failed: ' + str(e))
        self.assertListEqual(sorted(pl.tags.untyped_tags), ['passed_qc', 'sample'])
        self.assertListEqual(sorted(pl.tags.typed_tags), [('treatment', 'high_dose')])

        self.assertTupleEqual(pl.attributes, ('mz', 'intensity', 'odd_flag'))
        self.assertTupleEqual(pl.flag_attributes, ('odd_flag',))

        self.assertTrue(np.all(pl.flags == [1, 0] * 5))

        self.assertTupleEqual((len(pl.peaks), len(pl.peaks[0])), (5, 3))
        self.assertTupleEqual((len(pl.dtable), len(pl.dtable[0])), (10, 3))

    def test_pickle(self):
        pl = self._createPeakList()
        try:
            pstr = cp.dumps(pl)
            pl = cp.loads(pstr)
        except Exception, e:
            self.fail('PeakList pickle failed: ' + str(e))


if __name__ == '__main__':
    unittest.main()
