#!/usr/bin/env python
#  -*- coding: utf-8 -*-

# DO NOT try metadata.metadata.attr.
# All attribute methods overrided
class PeakList_Metadata(dict):
    """
    The PeakList_Metadata class.

    Dictionary-like container for PeakList metadata storage.

    :param args: iterable object of key-value pairs
    :param kwargs: metadata key-value pairs

    >>> PeakList_Metadata([('name', 'sample_1'), ('qc', False)])
    >>> PeakList_Metadata(name = 'sample_1', qc = False)

    metadata attributes can be accessed in both dictionary-like and property-like manners.

    >>> meta = PeakList_Metadata(name = 'sample_1', qc = False)
    >>> meta['name']
    sample_1
    >>> meta.qc
    False
    >>> del meta.qc
    >>> meta.has_key('qc')
    False

    .. warning::
        The *__getattr__*, *__setattr__*, and *__delattr__* methods are overrided. **DO NOT** assign a metadata object
        to another metadata object, e.g., metadata.metadata.attr = value.

    """

    def __getattr__(self, item):
        return self[item] if item in self else super().__getattribute__(item)

    def __setattr__(self, item, value):
        if item == '__dict__':
            raise ValueError('"__dict__" is not an acceptable metadata key')
        if type(value) == PeakList_Metadata:
            raise ValueError('metadata object is not an acceptable metadata value')

        if item not in self.__dict__:
            self[item] = value
        else:
            super().__setattr__(item, value)

    def __delattr__(self, item):
        if item in self:
            del self[item]
        else:
            super().__delattr__(item)
