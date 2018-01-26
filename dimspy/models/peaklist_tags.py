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
        Property of included tag values. Same typed tag values will be combined

        :getter: returns a tuple containing all the tag values, both typed and untyped tags
        :type: tuple

        """
        return tuple(self._untyped_tags) + \
               reduce(lambda x,y: x + ((y,) if y not in x else ()), self._typed_tags.values(), ())

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

    def has_tag(self, tag, tag_type = None):
        """
        Checks whether there exists a specific tag.

        :param tag: the tag value for checking
        :param tag_type: the tag type for checking, None indicates untyped tags
        :rtype: bool

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> tags.has_tag('untyped_tag1')
        True
        >>> tags.has_tag(tag_type1 = 'untyped_tag1')
        False
        >>> tags.has_tag('tag_value1')
        False

        """
        return (tag in self._untyped_tags) if tag_type is None else \
               (self._typed_tags.has_key(tag_type) and self._typed_tags[tag_type] == tag)

    def tag_of(self, tag_type = None):
        """
        Returns tag value of the given tag type, or tuple of untyped tags if tag_type is None.

        :param tag_type: valid tag type, or None for untyped tags
        :rtype: tuple (tag_type is None), or any (tag_type is not None)

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
        if len(args) != len(reduce(lambda x,y: x+([y] if y not in x else []), args, [])): # set cannot not apply to list
            raise ValueError('assigned untyped tags have duplication')
        if any(map(self.has_tag, args)) or \
           any(map(lambda x: self.has_tag(x[1], tag_type = x[0]), kwargs.items())):
            raise ValueError('tag(s) already exists')
        self._untyped_tags += list(args)
        self._typed_tags.update(kwargs)

    def drop_tags(self, *args, **kwargs):
        """
        Drops multiple typed and untyped tags.

        :param args: list of untyped tags
        :param kwargs: list of typed tags

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> tags.drop_tags(tag_type1 = 'tag_value1')
        >>> print tags
        untyped_tag1

        """
        if not all(map(self.has_tag, args)) or \
           not all(map(lambda x: self.has_tag(x[1], tag_type = x[0]), kwargs.items())):
            raise ValueError('tag(s) not exists')
        self._untyped_tags = list(filter(lambda x: x not in args, self.untyped_tags))
        self._typed_tags = dict(filter(lambda x: x not in kwargs.items(), self.typed_tags))

    def drop_tag_types(self, *args):
        """
        Drops multiple tag types and their tag values.

        :param args: list of tag types

        """
        self._typed_tags = dict(filter(lambda x: x[0] not in args, self.typed_tags))

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
