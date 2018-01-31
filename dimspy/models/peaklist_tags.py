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
from types import NoneType


class Tag(object):
    """
    The Tag class.

    This class is mainly used in PeakList and PeakMatrix classes for sample filtering.

    :param value: tag value, must be number (int, float), string (ascii, unicode), or Tag object (ignore ttype setting)
    :param ttype: tag type, must be string or None (untyped), default = None

    Single value will be treated as untyped tag:

    >>> tag = Tag(1)
    >>> tag == 1
    True
    >>> tag = Tag(1, 'batch')
    >>> tag == 1
    False

    """
    _value_valid_types = (int, float, str, unicode)
    _ttype_valid_types = (NoneType, str, unicode)

    def __init__(self, value, ttype = None):
        self._value, self._type = None, None
        self.value, self.ttype = (value.value, value.ttype) if isinstance(value, Tag) else (value, ttype)

    @property
    def value(self):
        """
        Property of tag value.

        :getter: returns the value of the tag
        :setter: set the tag value, must be number or string
        :type: int, float, str, unicode

        """
        return self._value

    @value.setter
    def value(self, value): # numpy types should be manually converted
        if not isinstance(value, Tag._value_valid_types):
            raise TypeError('Tag value must be string or number')
        self._value = value

    @property
    def ttype(self):
        """
        Property of tag type. None indicates untyped tag.

        :getter: returns the type of the tag
        :setter: set the tag type, must be None or string
        :type: None, str, unicode

        """
        return self._type

    @ttype.setter
    def ttype(self, value):
        if not isinstance(value, Tag._ttype_valid_types):
            raise TypeError('Tag type must be string or None')
        if value in ('None', ''): # reserve for hdf5 protal
            raise KeyError('["%s"] is not an acceptable tag type' % str(value))
        self._type = None if value is None else str(value)

    @property
    def typed(self):
        """
        Property to decide if the tag is typed or untyped.

        :getter: returns typed status of the tag
        :type: bool

        """
        return not self._type is None

    def __eq__(self, other):
        if not isinstance(other, Tag) and not isinstance(other, Tag._value_valid_types):
            raise TypeError('undefined comparison between Tag and %s' % type(other).__name__)
        v, t = (other.value, other.ttype) if isinstance(other, Tag) else (other, None)
        # Tag.ttype can only be str or None. It should be safe to use == for None comparison here
        # But still need to be aware that it (although highly unlikely) may fail
        return v == self.value and t == self.ttype

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return str(self._value) if self._type is None else (self._type + ':' + str(self._value))


class PeakList_Tags(object):
    """
    The PeakList_Tags class.

    Container for both typed and untyped tags. This class is mainly used in PeakList and PeakMatrix classes for sample filtering.
    For a PeakList the tag types must be unique, but not the tag values (unless they are untyped).
    For instance, PeakList can have tags batch = 1 and plate = 1, but not batch = 1 and batch = 2, or (untyped) 1 and (untyped) 1.
    Single value will be treated as untyped tag.

    :param args: list of untyped tags
    :param kwargs: list of typed tags. Only one tag value can be assigned to a specific tag type

    >>> PeakList_Tags('untyped_tag1', Tag('untyped_tag2'), Tag('typed_tag', 'tag_type'))
    >>> PeakList_Tags(tag_type1 = 'tag_value1', tag_type2 = 'tag_value2')

    """

    def __init__(self, *args, **kwargs):
        self._tags = []
        for v in args: self.add_tag(v)
        for k,v in kwargs.items(): self.add_tag(v,k)

    # build-ins
    def __str__(self):
        return self.to_str()

    def __contains__(self, item):
        return item in self._tags

    def __len__(self):
        return len(self._tags)

    # properties
    @property
    def tag_types(self):
        """
        Property of included tag types. None indicates untyped tags included.

        :getter: returns a set containing all the tag types of the typed tags
        :type: set

        """
        return set(map(lambda x: x.ttype, self._tags))

    @property
    def tag_values(self):
        """
        Property of included tag values. Same tag values will be merged

        :getter: returns a set containing all the tag values, both typed and untyped tags
        :type: set

        """
        return set(map(lambda x: x.value, self._tags))

    @property
    def tags(self):
        """
        Property of all included tags.

        :getter: returns a tuple containing all the tags, both typed and untyped
        :type: tuple

        """
        return tuple(self._tags)

    @property
    def typed_tags(self):
        """
        Property of included typed tags.

        :getter: returns a tuple containing all the typed tags
        :type: tuple

        """
        return tuple(filter(lambda x: x.typed, self._tags))

    @property
    def untyped_tags(self):
        """
        Property of included untyped tags.

        :getter: returns a tuple containing all the untyped tags
        :type: tuple

        """
        return tuple(filter(lambda x: not x.typed, self._tags))

    # methods
    def has_tag(self, tag, tag_type = None):
        """
        Checks whether there exists a specific tag.

        :param tag: the tag for checking
        :param tag_type: the type of the tag
        :rtype: bool

        >>> tags = PeakList_Tags('untyped_tag1', Tag('tag_value1', 'tag_type1'))
        >>> tags.has_tag('untyped_tag1')
        True
        >>> tags.has_tag('typed_tag1')
        False
        >>> tags.has_tag(Tag('tag_value1', 'tag_type1'))
        True
        >>> tags.has_tag('tag_value1', 'tag_type1')
        True

        """
        return (tag in self._tags) if isinstance(tag, Tag) or tag_type is None else \
               (Tag(tag, tag_type) in self._tags)

    def has_tag_type(self, tag_type = None):
        """
        Checks whether there exists a specific tag type.

        :param tag_type: the tag type for checking, None indicates untyped tags
        :rtype: bool

        """
        return tag_type in self.tag_types

    def tag_of(self, tag_type = None):
        """
        Returns tag value of the given tag type, or tuple of untyped tags if tag_type is None.

        :param tag_type: valid tag type, None for untyped tags
        :rtype: Tag, or None if tag_type not exists

        """
        t = filter(lambda x: x.ttype == tag_type, self._tags)
        return None if len(t) == 0 else tuple(t) if tag_type is None else t[0]

    def add_tag(self, tag, tag_type = None):
        """
        Adds typed or untyped tag.

        :param tag: tag or tag value to add
        :param tag_type: type of the tag value

        >>> tags = PeakList_Tags()
        >>> tags.add_tag('untyped_tag1')
        >>> tags.add_tag(Tag('typed_tag1', 'tag_type1'))
        >>> tags.add_tag(tag_type2 = 'typed_tag2')

        """
        if tag_type is not None and self.has_tag_type(tag_type):
            raise KeyError('tag type %s already exists' % tag_type)
        tag = Tag(tag, tag_type)
        if self.has_tag(tag):
            raise ValueError('tag already exist')
        self._tags += [tag]

    def drop_tag(self, tag, tag_type = None):
        """
        Drops typed and untyped tag.

        :param tag: tag or tag value to drop
        :param tag_type: type of the tag value

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> tags.drop_tag(Tag('tag_value1', 'tag_type1'))
        >>> print tags
        untyped_tag1

        """
        t = Tag(tag, tag_type)
        self._tags = filter(lambda x: x != t, self._tags)

    def drop_tag_type(self, tag_type = None):
        """
        Drops the tag with the given type.

        :param tag_type: tag type to drop, None (untyped) may drop multiple tags

        """
        self._tags = filter(lambda x: x.ttype != tag_type, self._tags)

    def drop_all_tags(self):
        """
        Drops all tags, both typed and untyped.

        """
        self._tags = []

    # portals
    def to_list(self):
        """
        Exports tags to a list. Each element is a tuple of (tag value, tag type).

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> tags.to_list()
        [('untyped_tag1', None), ('tag_value1', 'tag_type1')]

        :rtype: list

        """
        return [(t.value, t.ttype) for t in self._tags]

    def to_str(self):
        """
        Exports tags to a string. It can also be used inexplicitly as

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> print tags
        untyped_tag1, tag_type1:tag_value1

        :rtype: str

        """
        return join(map(str, self._tags), ', ')
