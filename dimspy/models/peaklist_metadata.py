#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
peaklist_metadata: PeakList metadata class, for internal use only

author(s): Albert
origin: Sep. 27, 2016

"""


# For internal use only.
# DO NOT try metadata.metadata.attr.
# All attribute methods overrided
class PeakList_Metadata(dict):
    def __init__(self, *args, **kwargs):
        super(PeakList_Metadata, self).__init__(*args, **kwargs)

    def __getattr__(self, item):
        return self[item] if self.has_key(item) else \
               super(PeakList_Metadata, self).__getattribute__(item)

    def __setattr__(self, item, value):
        if item == '__dict__': raise ValueError('"__dict__" is not an acceptable metadata key')
        if type(value) == PeakList_Metadata: raise ValueError('metadata object is not an acceptable metadata value')
        if not self.__dict__.has_key(item):
            self[item] = value
        else:
            super(PeakList_Metadata, self).__setattr__(item, value)

    def __delattr__(self, item):
        if self.has_key(item):
            del self[item]
        else:
            super(PeakList_Metadata, self).__delattr__(item)

