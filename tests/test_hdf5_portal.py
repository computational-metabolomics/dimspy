#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_hdf5_portal: 

author(s): Albert
origin: 05-14-2017

"""


import unittest, os
import numpy as np
from dimspy.models.peaklist_tags import Tag
from dimspy.models.peaklist import PeakList
from dimspy.process.peak_alignment import align_peaks
from dimspy.portals.hdf5_portal import save_peaklists_as_hdf5, load_peaklists_from_hdf5
from dimspy.portals.hdf5_portal import save_peak_matrix_as_hdf5, load_peak_matrix_from_hdf5


class HDF5PortalsTestCase(unittest.TestCase):
    @staticmethod
    def _createPeaklists():
        _mzs = lambda: sorted(np.random.uniform(100, 1200, size = 100))
        _ints = lambda: np.abs(np.random.normal(100, 10, size = 100))

        pkls = [
            PeakList('sample_1_1', _mzs(), _ints(), mz_range = (100, 1200)),
            PeakList('sample_1_2', _mzs(), _ints(), mz_range = (100, 1200)),
            PeakList('QC_1', _mzs(), _ints(), mz_range = (100, 1200)),
            PeakList('sample_2_1', _mzs(), _ints(), mz_range = (100, 1200)),
            PeakList('sample_2_2', _mzs(), _ints(), mz_range = (100, 1200)),
            PeakList('QC_2', _mzs(), _ints(), mz_range = (100, 1200)),
        ]

        for t in ('sample', Tag('compound_1', 'treatment'), Tag('1hr', 'time_point'), Tag(1, 'plate')): pkls[0].tags.add_tag(t)
        for t in ('sample', Tag('compound_1', 'treatment'), Tag('6hr', 'time_point'), Tag(1, 'plate')): pkls[1].tags.add_tag(t)
        for t in ('qc', Tag(1, 'plate')): pkls[2].tags.add_tag(t)
        for t in ('sample', Tag('compound_2', 'treatment'), Tag('1hr', 'time_point'), Tag(2, 'plate')): pkls[3].tags.add_tag(t)
        for t in ('sample', Tag('compound_2', 'treatment'), Tag('6hr', 'time_point'), Tag(2, 'plate')): pkls[4].tags.add_tag(t)
        for t in ('qc', Tag(2, 'plate')): pkls[5].tags.add_tag(t)

        for p in pkls: p.add_attribute('snr', np.random.uniform(300, 400, size = 100))
        for p in pkls: p.add_attribute('quad_flag', [0, 1, 1, 1] * 25, is_flag = True)
        for p in pkls: p.add_attribute('lab', [chr(i%26+97) for i in range(100)], flagged_only = False)
        return pkls

    def test_peaklist_portal(self):
        pkls = self._createPeaklists()

        save_peaklists_as_hdf5(pkls, '.test_peaklist.hdf5')
        npkls = load_peaklists_from_hdf5('.test_peaklist.hdf5')

        self.assertListEqual(map(lambda x: x.size, npkls), [75] * 6)
        self.assertListEqual(map(lambda x: x.full_size, npkls), [100] * 6)
        self.assertTrue(all(map(lambda x: np.allclose(x[0].mz_all, x[1].mz_all), zip(pkls, npkls))))
        self.assertTrue(all(map(lambda x: np.allclose(x[0].intensity, x[1].intensity), zip(pkls, npkls))))
        self.assertTrue(all(map(lambda x: np.allclose(x[0].snr, x[1].snr, atol = 1e-30), zip(pkls, npkls))))
        self.assertTrue(all(map(lambda x: np.all(x[0].quad_flag == x[1].quad_flag), zip(pkls, npkls))))
        self.assertTrue(all(map(lambda x: np.all(x[0].lab == x[1].lab), zip(pkls, npkls))))
        self.assertTrue(all(map(lambda x: x[0].metadata.keys() == x[1].metadata.keys(), zip(pkls, npkls))))
        self.assertTrue(all(map(lambda x: x[0].tags.tag_types == x[1].tags.tag_types, zip(pkls, npkls))))
        self.assertTrue(all(map(lambda x: x[0].tags.tag_values == x[1].tags.tag_values, zip(pkls, npkls))))

    def test_peak_matrix_portal(self):
        pkls = self._createPeaklists()
        pm = align_peaks(pkls, ppm = 2.0, block_size = 10, ncpus = 2)

        pm.mask_tags('qc')

        pnum = pm.full_shape[1]
        pm.add_flag('odd_flag', ([0, 1] * int(pnum/2.+1))[:pnum])
        pm.add_flag('qua_flag', ([0, 0, 0, 1] * int(pnum/4.+1))[:pnum], flagged_only = False)

        save_peak_matrix_as_hdf5(pm, '.test_peak_matrix.hdf5')
        npm = load_peak_matrix_from_hdf5('.test_peak_matrix.hdf5')

        self.assertEqual(pm.shape, npm.shape)
        self.assertEqual(pm.full_shape, npm.full_shape)
        self.assertTupleEqual(pm.attributes, npm.attributes)
        self.assertTrue(np.allclose(pm.mz_matrix, npm.mz_matrix))
        self.assertTrue(np.allclose(pm.intensity_matrix, npm.intensity_matrix))
        self.assertTrue(np.allclose(pm.attr_matrix('snr'), npm.attr_matrix('snr')))
        self.assertTrue(np.all(pm.attr_matrix('lab') == npm.attr_matrix('lab')))
        self.assertTrue(np.all( pm.property('present_matrix', flagged_only = False) ==
                               npm.property('present_matrix', flagged_only = False)))
        self.assertEqual(pm.peaklist_tag_types, npm.peaklist_tag_types)
        self.assertEqual(pm.peaklist_tag_values, npm.peaklist_tag_values)
        self.assertTrue(np.all(pm.mask == npm.mask))
        self.assertTrue(np.all(pm.flag_values('odd_flag') == npm.flag_values('odd_flag')))
        self.assertTrue(np.all(pm.flag_values('qua_flag') == npm.flag_values('qua_flag')))
        self.assertTrue(np.all(pm.flags == npm.flags))

    def tearDown(self):
        if os.path.isfile('.test_peaklist.hdf5'): os.remove('.test_peaklist.hdf5')
        if os.path.isfile('.test_peak_matrix.hdf5'): os.remove('.test_peak_matrix.hdf5')


if __name__ == '__main__':
    unittest.main()
