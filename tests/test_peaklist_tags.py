#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_peaklist_tags

author(s): Albert
origin: 04-29-2017

"""


import unittest
import cPickle as cp
from dimspy.models.peaklist_tags import Tag, PeakList_Tags


class TagTestCase(unittest.TestCase):
    def test_tag_creation(self):
        try:
            tag1 = Tag('1')
            tag2 = Tag(2, 'batch')
            tag3 = Tag(tag2)
        except Exception, e:
            self.fail('create tag object failed: ' + str(e))

        self.assertTrue(tag1.value == '1' and tag1.ttype is None)
        self.assertTrue(tag2.value == 2 and tag2.ttype == 'batch')
        self.assertTrue(tag3.value == 2 and tag3.ttype == 'batch')

        self.assertRaises(TypeError, lambda: Tag((3,4,5)))
        self.assertRaises(TypeError, lambda: Tag(6, ttype = [7,8]))
        self.assertRaises(KeyError, lambda: Tag(9, ttype = 'None'))

    def test_tag_property(self):
        tag = Tag('value', ttype = 'type')
        self.assertTrue(tag.typed)

        tag.value = 1
        tag.ttype = None
        self.assertTrue(tag.value == 1 and tag.ttype is None)
        self.assertFalse(tag.typed)

        def _assign_val(): tag.value = (1,2,3)
        self.assertRaises(TypeError, _assign_val)
        def _assign_type(t): tag.ttype = t
        self.assertRaises(TypeError, lambda: _assign_type([4,5]))
        self.assertRaises(KeyError, lambda: _assign_type('None'))

    def test_tag_magic(self):
        tag = Tag(1, ttype = 'type')

        self.assertEqual(tag, Tag(1, 'type'))
        self.assertNotEqual(tag, 1)

        tag.ttype = None
        self.assertEqual(tag, 1)
        self.assertTrue(1 == tag)
        self.assertFalse(1 != tag)
        self.assertTrue(2 != tag)
        self.assertTrue(tag in (1, 2, 3))
        self.assertTrue(1 in (tag, 2, 3))

        self.assertEqual(str(tag), '1')
        tag.ttype = 'type'
        self.assertEqual(str(tag), 'type:1')

class PeakListTagsTestCase(unittest.TestCase):
    @staticmethod
    def _createTags():
        return PeakList_Tags(0, 'str_tag', u'ustr_tag', Tag(1, 'typed_tag1'), typed_tag2 = 2)

    def test_pl_tags_creation(self):
        try:
            self._createTags()
        except Exception, e:
            self.fail('create tags object failed: ' + str(e))

    def test_pl_tags_properties(self):
        tags = self._createTags()
        self.assertEqual(tags.tag_types, {None, 'typed_tag1', 'typed_tag2'})
        self.assertEqual(tags.tag_values, {0, 1, 2, 'str_tag', u'ustr_tag'})
        self.assertEqual(len(tags), 5)
        self.assertTrue(all(map(lambda x: x.ttype is not None, tags.typed_tags)))
        self.assertTrue(all(map(lambda x: x.ttype is None, tags.untyped_tags)))

    def test_pl_tags_checking_methods(self):
        tags = self._createTags()

        self.assertTrue(tags.has_tag_type('typed_tag1') and tags.has_tag_type('typed_tag2'))
        self.assertTrue(tags.has_tag_type(None))
        self.assertFalse(tags.has_tag_type('not_exist'))

        self.assertTrue(Tag(2, 'typed_tag2') in tags)
        self.assertTrue(tags.has_tag(0) and tags.has_tag('str_tag') and tags.has_tag(u'ustr_tag'))
        self.assertTrue(tags.has_tag(1, 'typed_tag1') and tags.has_tag(2, 'typed_tag2'))
        self.assertTrue(tags.has_tag(Tag(1, 'typed_tag1')))
        self.assertFalse(tags.has_tag(0, 'typed_tag1'))
        self.assertFalse(tags.has_tag(1) or tags.has_tag(2))
        self.assertFalse(tags.has_tag('not_exist') or tags.has_tag(1, 'wrong_type'))

        self.assertTupleEqual(tags.tag_of(), (0, 'str_tag', u'ustr_tag'))
        self.assertTrue(tags.tag_of('typed_tag1').value == 1 and tags.tag_of('typed_tag2').value == 2)
        self.assertTrue(tags.tag_of('not_such_type') is None)

    def test_pl_tags_adding_methods(self):
        tags = self._createTags()

        self.assertRaises(KeyError, lambda: tags.add_tag(3, 'typed_tag1'))
        self.assertRaises(ValueError, lambda: tags.add_tag(0))
        self.assertRaises(ValueError, lambda: tags.add_tag(u'ustr_tag'))

        tags.add_tag(1)
        tags.add_tag(1, 'typed_tag3')
        tags.add_tag(Tag('new_value', 'typed_tag4'))
        self.assertEqual(tags.tag_types, {None, 'typed_tag1', 'typed_tag2', 'typed_tag3', 'typed_tag4'})
        self.assertEqual(tags.tag_values, {0, 1, 2, 'new_value', 'str_tag', u'ustr_tag'})

    def test_pl_tags_dropping_methods(self):
        tags = self._createTags()

        tags.drop_tag(0)
        tags.drop_tag(1)
        tags.drop_tag(1, 'wrong_type')
        self.assertEqual(tags.tag_types, {None, 'typed_tag1', 'typed_tag2'})
        self.assertEqual(tags.tag_values, {1, 2, 'str_tag', u'ustr_tag'})
        tags.drop_tag('str_tag')
        tags.drop_tag(u'ustr_tag')
        self.assertEqual(tags.tag_types, {'typed_tag1', 'typed_tag2'})
        self.assertEqual(tags.tag_values, {1, 2})

        tags.drop_tag_type('typed_tag1')
        self.assertEqual(tags.tag_types, {'typed_tag2'})
        self.assertEqual(tags.tag_values, {2})

        tags.drop_all_tags()
        self.assertTupleEqual(tags.tags, ())

    def test_pl_tags_portals(self):
        tags = self._createTags()
        self.assertListEqual(tags.to_list(), [(0, None), ('str_tag', None), (u'ustr_tag', None), (1, 'typed_tag1'), (2, 'typed_tag2')])
        self.assertEqual(tags.to_str(), '0, str_tag, ustr_tag, typed_tag1:1, typed_tag2:2')
        self.assertEqual(str(tags), '0, str_tag, ustr_tag, typed_tag1:1, typed_tag2:2')

    def test_pl_tags_pickle(self):
        tags = self._createTags()
        try:
            tstr = cp.dumps(tags)
            tags = cp.loads(tstr)
        except Exception, e:
            self.fail('tags pickle failed: ' + str(e))
        self.assertEqual(tags.tag_types, {None, 'typed_tag1', 'typed_tag2'})
        self.assertEqual(tags.tag_values, {0, 1, 2, 'str_tag', u'ustr_tag'})
        self.assertEqual(len(tags), 5)


if __name__ == '__main__':
    unittest.main()
