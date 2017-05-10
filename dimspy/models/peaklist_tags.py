#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
peaklist_tags: PeakList tags class, for internal use only

author(s): Albert
origin: Sep. 27, 2016

"""


from string import join


# For internal use only.
class PeakList_Tags(object):
    def __init__(self, *args, **kwargs):
        self._untyped_tags = []
        self._typed_tags = {}
        self.add_tags(*args, **kwargs)

    # build-ins
    def __str__(self):
        return self.to_str()

    # properties
    @property
    def tag_types(self):
        return tuple(self._typed_tags.keys())

    @property
    def tag_values(self):
        return tuple(self._untyped_tags + self._typed_tags.values())

    @property
    def typed_tags(self):
        return tuple(self._typed_tags.items())

    @property
    def untyped_tags(self):
        return tuple(self._untyped_tags)

    # methods
    def has_tag_type(self, tag_type):
        return tag_type in self.tag_types

    def has_tag(self, *args, **kwargs):
        if len(args) + len(kwargs) > 1:
            raise ValueError('searching multiple tags is not allowded')
        return (args[0] in self.tag_values) if len(args) > 0 else (kwargs.items()[0] in self.typed_tags)

    def tag_of(self, tag_type = None):
        if not (tag_type is None or self.has_tag_type(tag_type)):
            raise KeyError('unknown tag type [%s]' % tag_type)
        return self.untyped_tags if tag_type is None else self._typed_tags[tag_type]

    def add_tags(self, *args, **kwargs):
        if kwargs.has_key('None'):
            raise KeyError('["None"] is not an acceptable tag type') # reserve for hdf5 protal
        if None in args or None in kwargs.values():
            raise ValueError('None is not an acceptable tag value')
        if any(map(lambda x: x in kwargs.values(), args)) or \
           any(map(lambda x: x in args, kwargs.values())):
            raise ValueError('assigned tags have duplication')
        if any(map(self.has_tag, args)) or any(map(self.has_tag, kwargs.values())):
            raise ValueError('tag(s) already exists')
        self._untyped_tags += list(args)
        self._typed_tags.update(kwargs)

    def drop_tags(self, *args):
        self._untyped_tags = filter(lambda x: x not in args, self._untyped_tags)
        self._typed_tags = dict(filter(lambda x: x[1] not in args, self._typed_tags.items()))

    def drop_tag_types(self, *args):
        self._typed_tags = dict(filter(lambda x: x[0] not in args, self._typed_tags.items()))

    def drop_all_tags(self):
        self._untyped_tags = []
        self._typed_tags = {}

    # portals
    def to_list(self):
        return list(self.untyped_tags + self.typed_tags)

    def to_str(self):
        return join(map(str, self.untyped_tags) + map(lambda x: join(map(str, x), ':'), self.typed_tags), ', ')
