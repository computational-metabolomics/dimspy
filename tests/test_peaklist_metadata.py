#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_peaklist_metadata

author(s): Albert
origin: 04-28-2017

"""


import unittest
import cPickle as cp
from dimspy.models.peaklist_metadata import PeakList_Metadata


class PeakListMetadataTestCase(unittest.TestCase):
    @staticmethod
    def _createMetadata():
        return PeakList_Metadata((('a', 1), ('b', 2), ('c', 3)))

    def test_pl_meta_creation(self):
        try:
            self._createMetadata()
        except Exception, e:
            self.fail('create metadata object failed: ' + str(e))

    def test_pl_meta_operations(self):
        meta = self._createMetadata()

        self.assertListEqual(sorted(meta.keys()), ['a', 'b', 'c'])
        self.assertListEqual(sorted(meta.values()), [1, 2, 3])
        self.assertListEqual(sorted(meta.items()), [('a', 1), ('b', 2), ('c', 3)])
        self.assertTrue(meta['a'] == 1 and meta['b'] == 2 and meta['c'] == 3)
        self.assertTrue(meta.has_key('a') == True and meta.has_key('d') == False)
        self.assertTrue(meta.get('a', 4) == 1 and meta.get('d', 4) == 4)

        meta['a'] = 4
        self.assertEqual(meta['a'], 4)
        meta['d'] = 5
        self.assertEqual(meta['d'], 5)
        del meta['b']
        self.assertFalse(meta.has_key('b'))

    def test_pl_meta_pickle(self):
        meta = self._createMetadata()
        try:
            mstr = cp.dumps(meta)
            meta = cp.loads(mstr)
        except Exception, e:
            self.fail('metadata pickle failed: ' + str(e))
        self.assertTrue(meta['a'] == 1 and meta['b'] == 2 and meta['c'] == 3)


if __name__ == '__main__':
    unittest.main()
