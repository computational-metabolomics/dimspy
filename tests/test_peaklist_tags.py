#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_peaklist_tags

author(s): Albert
origin: 04-29-2017

"""


import unittest
import cPickle as cp
from dimspy.models.peaklist_tags import PeakList_Tags


class PeakListTagsTestCase(unittest.TestCase):
    @staticmethod
    def _createTags():
        return PeakList_Tags(1, [2, 3], 'str_tag', typed_tag1 = (4, 5), typed_tag2 = u'ustr_tag')

    def test_creation(self):
        try:
            self._createTags()
        except Exception, e:
            self.fail('create tags object failed: ' + str(e))

    def test_properties(self):
        tags = self._createTags()
        self.assertListEqual(sorted(tags.tag_types), ['typed_tag1', 'typed_tag2'])
        self.assertListEqual(sorted(tags.tag_values), [1, [2, 3], 'str_tag', (4, 5), u'ustr_tag'])
        self.assertListEqual(sorted(tags.typed_tags), [('typed_tag1', (4, 5)), ('typed_tag2', u'ustr_tag')])
        self.assertListEqual(sorted(tags.untyped_tags), [1, [2, 3], 'str_tag'])

    def test_checking_methods(self):
        tags = self._createTags()

        self.assertTrue(tags.has_tag_type('typed_tag1') and tags.has_tag_type('typed_tag2'))
        self.assertFalse(tags.has_tag_type('None') or tags.has_tag_type('not_exist'))

        self.assertTrue(tags.has_tag(1) and tags.has_tag([2, 3]) and tags.has_tag('str_tag'))
        self.assertTrue(tags.has_tag(typed_tag1 = (4, 5)) and tags.has_tag(typed_tag2 = u'ustr_tag'))
        self.assertTrue(tags.has_tag((4, 5)) and tags.has_tag(u'ustr_tag'))
        self.assertFalse(tags.has_tag('not_exist') or tags.has_tag(not_such_type = None))
        self.assertFalse(tags.has_tag(typed_tag1 = u'ustr_tag'))
        self.assertRaises(ValueError, lambda: tags.has_tag(1, typed_tag1 = (2, 3)))

        self.assertListEqual(sorted(tags.tag_of()), [1, [2, 3], 'str_tag'])
        self.assertTrue(tags.tag_of('typed_tag1') == (4, 5) and tags.tag_of('typed_tag2') == u'ustr_tag')
        self.assertRaises(KeyError, lambda: tags.tag_of('not_such_type'))

    def test_adding_methods(self):
        tags = self._createTags()

        self.assertRaises(KeyError, lambda: tags.add_tags(**{'None': 'none_tag'}))
        self.assertRaises(ValueError, lambda: tags.add_tags(None))
        self.assertRaises(ValueError, lambda: tags.add_tags(typed_tag3 = None))
        self.assertRaises(ValueError, lambda: tags.add_tags('str_tag3', 'str_tag4', typed_tag3 = 'str_tag3'))
        self.assertRaises(ValueError, lambda: tags.add_tags('str_tag', 'str_tag2'))
        self.assertRaises(ValueError, lambda: tags.add_tags(u'ustr_tag', ))
        self.assertRaises(ValueError, lambda: tags.add_tags(typed_tag2 = u'ustr_tag'))
        self.assertRaises(ValueError, lambda: tags.add_tags(typed_tag3 = u'ustr_tag'))

        try:
            tags.add_tags() # empty input should work
        except Exception, e:
            self.fail('tags add_tags failed: ' + str(e))
        self.assertListEqual(sorted(tags.typed_tags), [('typed_tag1', (4, 5)), ('typed_tag2', u'ustr_tag')])
        self.assertListEqual(sorted(tags.untyped_tags), [1, [2, 3], 'str_tag'])

        try:
            tags.add_tags(2, 'str_tag2', typed_tag3 = u'ustr_tag2')
        except Exception, e:
            self.fail('tags add_tags failed: ' + str(e))
        self.assertListEqual(sorted(tags.typed_tags),
                             [('typed_tag1', (4, 5)), ('typed_tag2', u'ustr_tag'), ('typed_tag3', u'ustr_tag2')])
        self.assertListEqual(sorted(tags.untyped_tags), [1, 2, [2, 3], 'str_tag', 'str_tag2'])

    def test_dropping_methods(self):
        tags = self._createTags()

        tags.drop_tags([2, 3], (4, 5), {6, 7})
        self.assertListEqual(sorted(tags.typed_tags), [('typed_tag2', u'ustr_tag')])
        self.assertListEqual(sorted(tags.untyped_tags), [1, 'str_tag'])

        tags.drop_tag_types('typed_tag1', 'typed_tag2', 'typed_tag3')
        self.assertTupleEqual(tags.typed_tags, ())

        tags.drop_all_tags()
        self.assertTupleEqual(tags.typed_tags, ())
        self.assertTupleEqual(tags.untyped_tags, ())

    def test_portals(self):
        tags = self._createTags()

        self.assertListEqual(tags.to_list(), [1, [2, 3], 'str_tag', ('typed_tag1', (4, 5)), ('typed_tag2', u'ustr_tag')])
        self.assertEqual(tags.to_str(), '1, [2, 3], str_tag, typed_tag1:(4, 5), typed_tag2:ustr_tag')
        self.assertEqual(str(tags), '1, [2, 3], str_tag, typed_tag1:(4, 5), typed_tag2:ustr_tag')

    def test_pickle(self):
        tags = self._createTags()
        try:
            tstr = cp.dumps(tags)
            tags = cp.loads(tstr)
        except Exception, e:
            self.fail('tags pickle failed: ' + str(e))
        self.assertListEqual(tags.to_list(), [1, [2, 3], 'str_tag', ('typed_tag1', (4, 5)), ('typed_tag2', u'ustr_tag')])


if __name__ == '__main__':
    unittest.main()
