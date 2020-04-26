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


from __future__ import annotations

from typing import Union


class Tag(object):
    """
    The Tag class.

    This class is mainly used in PeakList and PeakMatrix classes for sample filtering.

    :param value: Tag value, must be number (int, float), string (ascii, unicode), or Tag object (ignore ttype setting)
    :param ttype: Tag type, must be string or None (untyped), default = None

    Single value will be treated as untyped tag:

    >>> tag = Tag(1)
    >>> tag == 1
    True
    >>> tag = Tag(1, 'batch')
    >>> tag == 1
    False

    """

    def __init__(self, value: Union[int, float, str, Tag], ttype: Union[str, None] = None):
        self._value, self._type = None, None
        self.value, self.ttype = (value.value, value.ttype) if isinstance(value, Tag) else (value, ttype)

    @property
    def value(self):
        """
        Property of tag value.

        :getter: Returns the value of the tag
        :setter: Set the tag value, must be number or string
        :type: int, float, str, unicode

        """
        return self._value

    @value.setter
    def value(self, value: Union[int, float, str]):  # numpy types should be manually converted
        self._value = value

    @property
    def ttype(self):
        """
        Property of tag type. None indicates untyped tag.

        :getter: Returns the type of the tag
        :setter: Set the tag type, must be None or string
        :type: None, str, unicode

        """
        return self._type

    @ttype.setter
    def ttype(self, value: Union[str, None]):
        if value in ('None', ''):  # reserve for hdf5 protal
            raise KeyError('["%s"] is not an acceptable tag type' % value)
        self._type = None if value is None else value

    @property
    def typed(self):
        """
        Property to decide if the tag is typed or untyped.

        :getter: Returns typed status of the tag
        :type: bool

        """
        return not self._type is None

    def __eq__(self, other: Union[int, float, str, Tag]):
        v, t = (other.value, other.ttype) if isinstance(other, Tag) else (other, None)
        return v == self.value and ((t is None and self.ttype is None) or (t == self.ttype))

    def __ne__(self, other: Union[int, float, str, Tag]):
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

    :param args: List of untyped tags
    :param kwargs: List of typed tags. Only one tag value can be assigned to a specific tag type

    >>> PeakList_Tags('untyped_tag1', Tag('untyped_tag2'), Tag('typed_tag', 'tag_type'))
    >>> PeakList_Tags(tag_type1 = 'tag_value1', tag_type2 = 'tag_value2')

    """

    def __init__(self, *args, **kwargs):
        self._tags = []
        for v in args: self.add_tag(v)
        for k, v in list(kwargs.items()): self.add_tag(v, k)

    # build-ins
    def __str__(self):
        return self.to_str()

    def __contains__(self, item: Union[int, float, str, Tag]):
        return item in self._tags

    def __len__(self):
        return len(self._tags)

    # properties
    @property
    def tag_types(self):
        """
        Property of included tag types. None indicates untyped tags included.

        :getter: Returns a set containing all the tag types of the typed tags
        :type: set

        """
        return set([x.ttype for x in self._tags])

    @property
    def tag_values(self):
        """
        Property of included tag values. Same tag values will be merged

        :getter: Returns a set containing all the tag values, both typed and untyped tags
        :type: set

        """
        return set([x.value for x in self._tags])

    @property
    def tags(self):
        """
        Property of all included tags.

        :getter: Returns a tuple containing all the tags, both typed and untyped
        :type: tuple

        """
        return tuple(self._tags)

    @property
    def typed_tags(self):
        """
        Property of included typed tags.

        :getter: Returns a tuple containing all the typed tags
        :type: tuple

        """
        return tuple([x for x in self._tags if x.typed])

    @property
    def untyped_tags(self):
        """
        Property of included untyped tags.

        :getter: Returns a tuple containing all the untyped tags
        :type: tuple

        """
        return tuple([x for x in self._tags if not x.typed])

    # methods
    def has_tag(self, tag: Union[int, float, str, Tag], tag_type: Union[str, None] = None):
        """
        Checks whether there exists a specific tag.

        :param tag: The tag for checking
        :param tag_type: The type of the tag
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

    def has_tag_type(self, tag_type: Union[str, None] = None):
        """
        Checks whether there exists a specific tag type.

        :param tag_type: The tag type for checking, None indicates untyped tags
        :rtype: bool

        """
        return tag_type in self.tag_types

    def tag_of(self, tag_type: Union[str, None] = None):
        """
        Returns tag value of the given tag type, or tuple of untyped tags if tag_type is None.

        :param tag_type: Valid tag type, None for untyped tags
        :rtype: Tag, or None if tag_type not exists

        """
        t = [x for x in self._tags if x.ttype == tag_type]
        return None if len(t) == 0 else tuple(t) if tag_type is None else t[0]

    def add_tag(self, tag: Union[int, float, str, Tag], tag_type: Union[str, None] = None):
        """
        Adds typed or untyped tag.

        :param tag: Tag or tag value to add
        :param tag_type: Type of the tag value

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

    def drop_tag(self, tag: Union[int, float, str, Tag], tag_type: Union[str, None] = None):
        """
        Drops typed and untyped tag.

        :param tag: Tag or tag value to drop
        :param tag_type: Type of the tag value

        >>> tags = PeakList_Tags('untyped_tag1', tag_type1 = 'tag_value1')
        >>> tags.drop_tag(Tag('tag_value1', 'tag_type1'))
        >>> print(tags)
        untyped_tag1

        """
        t = Tag(tag, tag_type)
        self._tags = [x for x in self._tags if x != t]

    def drop_tag_type(self, tag_type: Union[str, None] = None):
        """
        Drops the tag with the given type.

        :param tag_type: Tag type to drop, None (untyped) may drop multiple tags

        """
        self._tags = [x for x in self._tags if x.ttype != tag_type]

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
        >>> print(tags)
        untyped_tag1, tag_type1:tag_value1

        :rtype: str

        """
        return str.join(', ', map(str, self._tags))
