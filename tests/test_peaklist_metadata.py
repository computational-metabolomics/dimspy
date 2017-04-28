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
    def test_create_metadata(self):
        try:
            PeakList_Metadata((('a', 1), ('b', 2), ('c', 3)))
        except Exception, e:
            self.fail('create metadata object failed: ' + str(e))

    def test_dict_operations(self):
        meta = PeakList_Metadata((('a', 1), ('b', 2), ('c', 3)))
        self.assertSetEqual(set(meta.keys()), {'a', 'b', 'c'}, 'metadata keys not match')
        self.assertSetEqual(set(meta.values()), {1, 2, 3}, 'metadata values not match')
        self.assertTrue(meta['a'] == 1 and meta['b'] == 2 and meta['c'] == 3, 'metadata items not match')
        self.assertTrue(meta.has_key('a') == True and meta.has_key('d') == False, 'metadata has_key() function failed')
        self.assertTrue(meta.get('a', 4) == 1 and meta.get('d', 4) == 4, 'metadata get() function failed')
        meta['a'] = 4
        self.assertEqual(meta['a'], 4, 'metadata setting function failed')
        meta['d'] = 5
        self.assertEqual(meta['d'], 5, 'metadata adding function failed')
        del meta['b']
        self.assertFalse(meta.has_key('b'), 'metadata deleting function failed')

    def test_pickle(self):
        meta = PeakList_Metadata((('a', 1), ('b', 2), ('c', 3)))
        try:
            mstr = cp.dumps(meta)
            meta = cp.loads(mstr)
        except Exception, e:
            self.fail('metadata pickle failed: ' + str(e))
        self.assertTrue(meta['a'] == 1 and meta['b'] == 2 and meta['c'] == 3, 'metadata pickle incorrect')


if __name__ == '__main__':
    unittest.main()
