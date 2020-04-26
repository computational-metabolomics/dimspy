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


import os
import unittest

import numpy as np
from dimspy.models.peaklist import PeakList
from dimspy.models.peaklist_tags import Tag
from dimspy.portals.txt_portal import save_peak_matrix_as_txt, load_peak_matrix_from_txt
from dimspy.portals.txt_portal import save_peaklist_as_txt, load_peaklist_from_txt
from dimspy.process.peak_alignment import align_peaks


class TxtPortalsTestCase(unittest.TestCase):
    def test_peaklist_portal(self):
        pkl = PeakList('peaklist', np.sort(np.random.uniform(100, 1200, size = 100)), np.random.normal(100, 10, size = 100))
        pkl.add_attribute('odd_flag', [0, 1] * 50, is_flag = True)

        save_peaklist_as_txt(pkl, '.test_peaklist.txt')
        npkl = load_peaklist_from_txt('.test_peaklist.txt', 'peaklist')

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
        for t in ('sample', Tag('compound_1', 'treatment'), Tag('1hr', 'time_point'), Tag(1, 'plate')): pkls[0].tags.add_tag(t)
        for t in ('sample', Tag('compound_1', 'treatment'), Tag('6hr', 'time_point'), Tag(1, 'plate')): pkls[1].tags.add_tag(t)
        for t in ('qc', Tag(1, 'plate')): pkls[2].tags.add_tag(t)
        for t in ('sample', Tag('compound_2', 'treatment'), Tag('1hr', 'time_point'), Tag(2, 'plate')): pkls[3].tags.add_tag(t)
        for t in ('sample', Tag('compound_2', 'treatment'), Tag('6hr', 'time_point'), Tag(2, 'plate')): pkls[4].tags.add_tag(t)
        for t in ('qc', Tag(2, 'plate')): pkls[5].tags.add_tag(t)

        pm = align_peaks(pkls, ppm = 2e+4, block_size = 10, ncpus = 2)
        pm.add_flag('odd_flag', ([0, 1] * int(pm.shape[1]/2+1))[:pm.shape[1]])
        pm.add_flag('qua_flag', ([0, 0, 1, 1] * int(pm.shape[1]/4+1))[:pm.shape[1]])

        save_peak_matrix_as_txt(pm, '.test_peak_matrix.txt', samples_in_rows = True, comprehensive = True,
                                rsd_tags = ('qc', Tag('compound_1', 'treatment'), Tag('compound_2', 'treatment')))
        npm = load_peak_matrix_from_txt('.test_peak_matrix.txt', samples_in_rows = True, comprehensive = 'auto')

        self.assertEqual(pm.shape, npm.shape)
        self.assertEqual(pm.full_shape, npm.full_shape)
        self.assertTrue(np.all(pm.flags == npm.flags))
        self.assertTrue(np.all(pm.flag_names == npm.flag_names))
        self.assertTrue(np.allclose(pm.intensity_matrix, npm.intensity_matrix))
        self.assertEqual(pm.peaklist_tag_types, npm.peaklist_tag_types)
        self.assertEqual(pm.peaklist_tag_values, npm.peaklist_tag_values)

    def tearDown(self):
        if os.path.isfile('.test_peaklist.txt'): os.remove('.test_peaklist.txt')
        if os.path.isfile('.test_peak_matrix.txt'): os.remove('.test_peak_matrix.txt')


if __name__ == '__main__':
    unittest.main()
