#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Albert
origin: 04-28-2017

"""


import unittest


class PeakListTestCase(unittest.TestCase):
    def setUp(self):
        size = 100
        mzs = sorted(np.random.uniform(100, 1000, size = size))
        ints = np.abs(np.random.normal(10, 3, size = size))
        snr = np.abs(np.random.normal(1000, 400, size = size))

        pl = PeakList('sample_peaklist', mzs, ints, mz_range = (100, 1000))
        pl.add_tags('sample', 'passed_qc', main_class = 'treatment_1')
        pl.add_attribute('snr', snr)
        pl.metadata.type = 'blank'

    def test_something(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
