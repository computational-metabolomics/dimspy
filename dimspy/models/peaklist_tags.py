#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
The PeakList tags class.

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 1.0.0

.. warning::
    This class is designed for PeakList and PeakMatrix internal use only.
   
"""


from string import join


class PeakList_Tags(object):
    """
    The PeakList_Tags class.

    Container for both typed and untyped tags. This class is mainly used in PeakList and PeakMatrix classes for sample filtering.

    :param args: list of untyped tags
    :param kwargs: list of typed tags. Only one tag value can be assigned to a specific tag type

    >>> PeakList_Tags('untyped_tag1', 'untyped_tag2')
    >>> PeakList_Tags(tag_type1 = 'tag_value1', tag_type2 = 'tag_value2')

    """

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
        """
        Property of included tag types.

        :getter: returns a tuple containing all the tag types of the typed tags
        :type: tuple

        """
        return tuple(self._typed_tags.keys())

    @property
    def tag_values(self):
        """
        Property of included tag values.

        :getter: returns a tuple containing all the tag values, both typed and untyped tags
        :type: tuple

        """
        return tuple(self._untyped_tags + self._typed_tags.values())

    @property
    def typed_tags(self):
        """
        Property of included typed tags.

        :getter: returns a tuple containing all the typed tags, each in the format of (tag_type, tag_value)
        :type: tuple

        """
        return tuple(self._typed_tags.items())

    @property
    def untyped_tags(self):
        """
        Property of included untyped tags.

        :getter: returns a tuple containing all the untyped tags
        :type: tuple

        """
        return tuple(self._untyped_tags)

    # methods
    def has_tag_type(self, tag_type):
        """
        Checks whether there exists a specific tag type.

        :param tag_type: the tag type for checking, None tag_type (i.e. untyped) always return True
        :rtype: bool

        """
        return tag_type is None or tag_type in self.tag_types

    def has_tag(self, *args, **kwargs):
        """
        Checks whether there exists a specific tag.

        :param args: **one** tag value, either typed or untyped
        :param kwargs: **one** tag_type = tag_value
        :rtype: bool

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> tags.has_tag('untyped_tag1')
        True
        >>> tags.has_tag(tag_type1 = 'untyped_tag1')
        False
        >>> tags.has_tag('untyped_tag1', 'tag_value1')
        ...
        ValueError: searching multiple tags is not allowded

        """
        if len(args) + len(kwargs) > 1:
            raise ValueError('searching multiple tags is not allowded')
        return (args[0] in self.tag_values) if len(args) > 0 else (kwargs.items()[0] in self.typed_tags)

    def tag_of(self, tag_type = None):
        """
        Returns tag value of the given tag type, or tuple of untyped tags if tag_type is None.

        :param tag_type: valid tag type, or None for untyped tags
        :rtype: same as tag_value (tag_type is not None), or tupe (tag_type is None)

        """
        if not (tag_type is None or self.has_tag_type(tag_type)):
            raise KeyError('unknown tag type [%s]' % tag_type)
        return self.untyped_tags if tag_type is None else self._typed_tags[tag_type]

    def add_tags(self, *args, **kwargs):
        """
        Adds multiple typed and untyped tags.

        :param args: list of untyped tags
        :param kwargs: list of typed tags. Only one tag value can be assigned to a specific tag type

        >>> PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')

        """
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
        """
        Drops multiple typed and untyped tags.

        :param args: list of tag values, both typed and untyped

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> tags.drop_tags('tag_value1')
        >>> print tags
        untyped_tag1

        """
        self._untyped_tags = filter(lambda x: x not in args, self._untyped_tags)
        self._typed_tags = dict(filter(lambda x: x[1] not in args, self._typed_tags.items()))

    def drop_tag_types(self, *args):
        """
        Drops multiple tag types and their tag values.

        :param args: list of tag types

        """
        self._typed_tags = dict(filter(lambda x: x[0] not in args, self._typed_tags.items()))

    def drop_all_tags(self):
        """
        Drops all tags, both typed and untyped.

        """
        self._untyped_tags = []
        self._typed_tags = {}

    # portals
    def to_list(self):
        """
        Exports tags to a list.

        :rtype: list

        """
        return list(self.untyped_tags + self.typed_tags)

    def to_str(self):
        """
        Exports tags to a string. It can also be used inexplicitly as

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> print tags
        untyped_tag1, tag_type1:tag_value1

        :rtype: str

        """
        return join(map(str, self.untyped_tags) + map(lambda x: join(map(str, x), ':'), self.typed_tags), ', ')
