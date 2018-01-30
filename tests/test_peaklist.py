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
        mzs = np.arange(0, 1000, step = 100)
        ints = np.abs(np.random.normal(10, 3, size = 10))
        pl = PeakList('sample_peaklist', mzs, ints, mz_range = (100, 1000), frag_mode = 'slb')
        return pl

    def test_pl_creation(self):
        try:
            self._createPeakList()
        except Exception, e:
            self.fail('create PeakList object failed: ' + str(e))

    def test_pl_properties(self):
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
            pl.tags.add_tag('sample')
            pl.tags.add_tag('passed_qc')
            pl.tags.add_tag('high_dose', tag_type = 'treatment')
        except Exception, e:
            self.fail('access tags failed: ' + str(e))
        self.assertEqual(set(pl.tags.tag_types), {None, 'treatment'})
        self.assertEqual(set(pl.tags.tag_values), {'sample', 'passed_qc', 'high_dose'})

        self.assertTupleEqual(pl.attributes, ('mz', 'intensity', 'odd_flag'))
        self.assertTupleEqual(pl.flag_attributes, ('odd_flag',))

        self.assertTrue(np.all(pl.flags == [1, 0] * 5))

        self.assertTupleEqual((len(pl.peaks), len(pl.peaks[0])), (5, 3))
        self.assertTupleEqual((len(pl.dtable), len(pl.dtable[0])), (10, 3))

    def test_pl_attribute_operations(self):
        pl = self._createPeakList()

        self.assertTrue(pl.has_attribute('mz'))
        self.assertFalse(pl.has_attribute('snr'))
        self.assertFalse(pl.has_attribute('flag')) # flag is not a real attribute

        snr = np.array([20, 0] * 5, dtype = int)
        pl.add_attribute('snr', snr, attr_dtype = float)
        pl.add_attribute('snr_flag', snr > 10, is_flag = True)
        self.assertTrue(np.all(pl.get_attribute('snr') > 10))
        self.assertTrue(np.all(pl.get_attribute('snr', flagged_only = False) == snr))

        pl.add_attribute('values_1', [0, 1] * 5, on_index = 2, flagged_only = False)
        self.assertEqual(pl.attributes[2], 'values_1')
        pl.set_attribute('values_1', [1] * 5) # snr_flag already masked odd peaks
        self.assertTrue(np.all(pl.get_attribute('values_1', flagged_only = False) == np.ones(10)))
        pl.set_attribute('values_1', [0] * 10, flagged_only = False)
        self.assertTrue(np.all(pl.get_attribute('values_1', flagged_only = False) == np.zeros(10)))
        pl.drop_attribute('values_1')
        self.assertFalse(pl.has_attribute('values_1'))

        self.assertRaises(AttributeError, lambda: pl.add_attribute('mz', np.ones(pl.size)))
        self.assertRaises(AttributeError, lambda: pl.add_attribute('snr', np.ones(pl.size)))
        self.assertRaises(AttributeError, lambda: pl.add_attribute('_dtable', np.ones(pl.size)))
        self.assertRaises(ValueError, lambda: pl.add_attribute('flags_1', np.arange(pl.size), is_flag = True))
        self.assertRaises(IndexError, lambda: pl.add_attribute('values_2', np.arange(pl.size), on_index = 0))
        self.assertRaises(IndexError, lambda: pl.add_attribute('values_2', np.arange(pl.size), on_index = -pl.shape[1]))
        self.assertRaises(ValueError, lambda: pl.add_attribute('values_2', np.arange(pl.full_size)))

        self.assertRaises(AttributeError, lambda: pl.set_attribute('flags', np.ones_like(pl.size)))
        self.assertRaises(AttributeError, lambda: pl.set_attribute('values_3', np.arange(pl.size)))
        self.assertRaises(ValueError, lambda: pl.set_attribute('mz', np.arange(10)[::-1], flagged_only = False))

        try:
            pl.set_attribute('mz', np.arange(10)[::-1], flagged_only = False, unsorted_mz = True)
        except Exception, e:
            self.fail('unsorted_mz flag failed: ' + str(e))
        self.assertTrue(np.all(pl.get_attribute('mz') == np.arange(10)[1::2])) # setting mz reversed the snr_flag

        self.assertRaises(AttributeError, lambda: pl.get_attribute('values_4'))
        self.assertRaises(AttributeError, lambda: pl.drop_attribute('values_4'))
        self.assertRaises(AttributeError, lambda: pl.drop_attribute('mz'))

    def test_pl_peaks_operations(self):
        pl = self._createPeakList()
        pl.add_attribute('value_flag', [1, 0] * 5, is_flag = True)

        # mz = 0, (100), 200, (300), 400, (500), 600, (700), 800, (900)
        pl.set_peak(4, (50, 10., True), flagged_only = False)
        self.assertTupleEqual((0, 50, 200, 600, 800), tuple(pl.get_attribute('mz')))

        # mz = 0, 50, (100), 200, (300), (500), 600, (700), 800, (900)
        pl.insert_peak((150, 10., True))
        self.assertTupleEqual((0, 50, 150, 200, 600, 800), tuple(pl.get_attribute('mz')))
        self.assertEqual(pl.full_size, 11)

        # mz = 0, 50, (100), 150, 200, (300), (500), 600, (700), 800, (900)
        pl.remove_peak((1,2))
        self.assertTupleEqual((0, 100, 200, 300, 500, 600, 700, 800, 900), tuple(pl.get_attribute('mz', flagged_only = False)))
        pl.remove_peak(1, flagged_only = False)
        self.assertTupleEqual((0, 200, 300, 500, 600, 700, 800, 900), tuple(pl.get_attribute('mz', flagged_only = False)))
        self.assertEqual(pl.size, 4)
        self.assertEqual(pl.full_size, 8)

        # mz = 0, 200, (300), (500), 600, (700), 800, (900)
        self.assertRaises(AttributeError, lambda: pl.cleanup_unflagged_peaks('mz'))
        self.assertRaises(AttributeError, lambda: pl.cleanup_unflagged_peaks('not_exists'))
        pl.cleanup_unflagged_peaks('value_flag')
        self.assertEqual(pl.full_size, pl.size)
        pl.cleanup_unflagged_peaks()
        self.assertTupleEqual((0, 200, 600, 800), tuple(pl.get_attribute('mz')))

    def test_pl_build_ins(self):
        pl = self._createPeakList()

        try:
            str(pl)
        except Exception, e:
            self.fail('__str__ failed: ' + str(e))
        self.assertEqual(len(pl), 10)

        pl.add_attribute('value_flag', [1, 0] * 5, is_flag = True)
        # mz = 0, (100), 200, (300), 400, (500), 600, (700), 800, (900)
        self.assertEqual(len(pl), 5)

        self.assertListEqual([0, 200, 400, 600, 800], pl.mz.tolist())
        self.assertListEqual(np.arange(0, 1000, step = 100).tolist(), pl.mz_all.tolist())

        self.assertListEqual([0, 200, 400, 600, 800], pl['mz'].tolist())
        self.assertListEqual([0, 200, 400], list(zip(*pl[:3].tolist())[0]))

    def test_pl_exports(self):
        pl = self._createPeakList()

        try:
            lst = pl.to_list()
        except Exception, e:
            self.fail('to_list function failed: ' + str(e))
        self.assertListEqual(np.arange(0, 1000, step = 100).tolist(), list(lst[0]))

        try:
            psr = pl.to_str(',')
        except Exception, e:
            self.fail('to_str function failed: ' + str(e))
        self.assertListEqual(np.arange(0, 1000, step = 100).tolist(),
                             map(float, zip(*map(lambda x: x.split(','), psr.split('\n')[1:]))[0]))

    def test_pl_pickle(self):
        pl = self._createPeakList()
        try:
            pstr = cp.dumps(pl)
            pl = cp.loads(pstr)
        except Exception, e:
            self.fail('PeakList pickle failed: ' + str(e))
        self.assertTupleEqual(pl.attributes, ('mz', 'intensity'))


if __name__ == '__main__':
    unittest.main()
