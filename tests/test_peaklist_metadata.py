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


import pickle as cp
import unittest

from dimspy.models.peaklist_metadata import PeakList_Metadata


class PeakListMetadataTestCase(unittest.TestCase):
    @staticmethod
    def _createMetadata():
        return PeakList_Metadata((('a', 1), ('b', 2), ('c', 3)))

    def test_pl_meta_creation(self):
        try:
            self._createMetadata()
        except Exception as e:
            self.fail('create metadata object failed: ' + str(e))

    def test_pl_meta_operations(self):
        meta = self._createMetadata()

        self.assertListEqual(sorted(meta.keys()), ['a', 'b', 'c'])
        self.assertListEqual(sorted(meta.values()), [1, 2, 3])
        self.assertListEqual(sorted(meta.items()), [('a', 1), ('b', 2), ('c', 3)])
        self.assertTrue(meta['a'] == 1 and meta['b'] == 2 and meta['c'] == 3)
        self.assertTrue(('a' in meta) == True and ('d' in meta) == False)
        self.assertTrue(meta.get('a', 4) == 1 and meta.get('d', 4) == 4)

        meta['a'] = 4
        self.assertEqual(meta['a'], 4)
        meta['d'] = 5
        self.assertEqual(meta['d'], 5)
        del meta['b']
        self.assertFalse('b' in meta)

    def test_pl_meta_pickle(self):
        meta = self._createMetadata()
        try:
            mstr = cp.dumps(meta)
            meta = cp.loads(mstr)
        except Exception as e:
            self.fail('metadata pickle failed: ' + str(e))
        self.assertTrue(meta['a'] == 1 and meta['b'] == 2 and meta['c'] == 3)


if __name__ == '__main__':
    unittest.main()
