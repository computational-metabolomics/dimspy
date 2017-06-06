#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
peakList: the PeakList data object class.

author(s): Albert Zhou
origin: Sep. 27, 2016

"""


import logging
import numpy as np
import numpy.lib.recfunctions as rfn
from collections import OrderedDict, Iterable
from string import join
from copy import deepcopy
from peaklist_metadata import PeakList_Metadata
from peaklist_tags import PeakList_Tags


class PeakList(object):
    _is_ordered = staticmethod(lambda vals: all(map(lambda x: x[0] - x[1] >= 0, zip(vals[1:], vals[:-1]))))

    def __init__(self, ID, mz, intensity, **metadata):
        if not self._is_ordered(mz):
            raise ValueError('mz values not in ascending order, check input data')
        if len(intensity) != len(mz):
            raise ValueError('mz values and intensities must have the same size')

        self._dtable = np.array(zip(mz, intensity), dtype = [('mz', 'f8'), ('intensity', 'f8')])

        self._id = str(ID)
        self._metadata = PeakList_Metadata(**metadata)

        self._tags = PeakList_Tags()
        self._flags = np.ones_like(mz, dtype = np.bool)
        self._flag_attrs = []

    # privates
    # reference: https://stackoverflow.com/questions/21818140/numpy-how-to-fill-multiple-fields-in-a-structured-array-at-once
    @staticmethod
    def _fields_view(arr, fields):
        dtype2 = np.dtype({name: arr.dtype.fields[name] for name in fields})
        return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

    # built-ins
    def __len__(self):
        return self.size

    def __str__(self):
        return self.to_str(', ')

    def __getitem__(self, item):
        if type(item) in (int, slice, list, np.ndarray):
            return self.get_peak(item)
        else:
            return self.get_attribute(item)

    def __setitem__(self, item, value):
        if type(item) in (int, slice, list, np.ndarray):
            self.set_peak(item, value)
        else:
            self.set_attribute(item, value)

    def __getattr__(self, item):
        if item.endswith('_all'):
            return self.get_attribute(item[:-4], flagged_only = False)
        else:
            return self.get_attribute(item, flagged_only = True)

    def __setattr__(self, item, value):
        if item != '_dtable' and self.has_attribute(item):
            if item.endswith('_all'):
                return self.set_attribute(item[:-4], value, flagged_only = False)
            else:
                return self.set_attribute(item, value, flagged_only = True)
        else:
            super(PeakList, self).__setattr__(item, value)

    # for proper pickling
    def __setstate__(self, state):
        super(PeakList, self).__setattr__('__dict__', state)

    def __getstate__(self):
        return self.__dict__

    # properties
    @property
    def ID(self):
        return self._id

    @property
    def size(self):
        return np.sum(self._flags)

    @property
    def full_size(self):
        return len(self._flags)

    @property
    def shape(self):
        return self.size, len(self.attributes)

    @property
    def full_shape(self):
        return self.full_size, len(self.attributes)

    @property
    def metadata(self):
        return self._metadata

    @property
    def tags(self):
        return self._tags

    @property
    def attributes(self):
        return self._dtable.dtype.names

    @property
    def flag_attributes(self):
        return tuple(self._flag_attrs[:])

    @property
    def peaks(self):
        return self._dtable[self._flags].copy()

    @property
    def flags(self):
        return self._flags.copy()

    @property
    def dtable(self):
        return self._dtable # ref for faster access, using carefully

    # publics
    # manual maintainance
    # these functions are aotumatically called internally, think twice before you manually call them
    def sort_peaks_order(self):
        sids = np.argsort(self._dtable['mz'])  # do not use self.mz, recursive calling
        self._dtable = self._dtable[sids]
        self._flags = self._flags[sids]

    def calculate_flags(self):
        self._flags = np.ones_like(self._flags) if len(self._flag_attrs) == 0 else \
                      np.sum(self._dtable[self._flag_attrs].tolist(), axis = 1) == len(self._flag_attrs)
        if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]' % str(self.ID))
        return self.flags

    # attribute operations
    def has_attribute(self, attr_name):
        return attr_name in self.attributes

    def add_attribute(self, attr_name, attr_value, attr_dtype = None, is_flag = False,
                      on_index = None, flagged_only = True, invalid_value = np.nan):
        if self.__dict__.has_key(attr_name):
            raise AttributeError('attribute name already been used by property')
        if attr_name in ('mz', 'intensity', 'flags'):
            raise AttributeError('cannot add reserved attribute [%s]' % attr_name)
        if self.has_attribute(attr_name):
            raise AttributeError('attribute [%s] already exists' % attr_name)
        if is_flag and self.size > 0 and not (set(attr_value) in ({0}, {1}, {0, 1})):
            raise ValueError('flag attribute can only contain True / False values')
        if on_index is not None and not (-self.shape[1]+1 < on_index < 0 or 1 < on_index < self.shape[1]):
            raise IndexError('index [%d] out of (insertable) range' % on_index)

        adt = bool if is_flag else \
              attr_dtype if attr_dtype is not None else \
              attr_value.dtype.str if hasattr(attr_value, 'dtype') else \
              ('S%d' % max(map(len, attr_value))) if type(attr_value[0]) in (unicode, str) else \
              type(attr_value[0])
        if adt in (bool, 'bool', '|b1'): adt = 'b'  # fix numpy dtype bug

        if flagged_only:
            if len(attr_value) != self.size:
                raise ValueError('input attibute value size not match')
            nattr = np.array([invalid_value] * self._dtable.shape[0]).astype(adt)
            nattr[self._flags] = attr_value
        else:
            if len(attr_value) != self.full_size:
                raise ValueError('input attibute value size not match')
            nattr = np.array(attr_value).astype(adt)

        if on_index is None: on_index = self.shape[1]
        anames, atypes = map(list, zip(*self._dtable.dtype.descr))
        prevnm, restnm, resttp = anames[:on_index], anames[on_index:], atypes[on_index:]

        prevtb = np.array(nattr, dtype = [(attr_name, adt)]) if len(prevnm) == 0 else \
            rfn.append_fields(self._fields_view(self._dtable, prevnm), attr_name, nattr, dtypes = adt, usemask = False)
        self._dtable = prevtb if len(restnm) == 0 else \
            rfn.append_fields(prevtb, restnm, zip(*self._fields_view(self._dtable, restnm)), dtypes = resttp, usemask = False)

        if is_flag:
            self._flag_attrs += [attr_name]
            self.calculate_flags()

        return self

    def drop_attribute(self, attr_name):
        if attr_name in ('mz', 'intensity', 'flags'):
            raise AttributeError('cannot drop reserved attribute [%s]' % attr_name)
        if not self.has_attribute(attr_name):
            raise AttributeError('attribute [%s] not exists' % attr_name)

        self._dtable = self._dtable[list(filter(lambda x: x != attr_name, self.attributes))]

        if attr_name in self._flag_attrs:
            logging.warning('flags recalculated, unflagged peaks may contain incorrect values')
            self._flag_attrs = filter(lambda x: x != attr_name, self._flag_attrs)
            self.calculate_flags()

        return self

    def set_attribute(self, attr_name, attr_value, flagged_only = True, unsorted_mz = False):
        if attr_name == 'flags':
            raise AttributeError('cannot assign read-only attribute [flag]')
        if not self.has_attribute(attr_name):
            raise AttributeError('attribute [%s] not exists' % attr_name)
        if attr_name == 'mz' and not (unsorted_mz or self._is_ordered(attr_value)):
            raise ValueError('attribute [mz] not in ascending order')
        if attr_name != 'mz' and unsorted_mz == True:
            logging.warning('setting unsorted_mz flag for non-mz attribute, ignore')

        # May raise FutureWarning, but that shold be a false alarm, because attr_name is a string rather than a list
        # we are not accessing multiple fields by their names here
        self._dtable[attr_name][self._flags if flagged_only else slice(None)] = attr_value
        if attr_name == 'mz' and unsorted_mz: self.sort_peaks_order()
        if attr_name in self._flag_attrs: self.calculate_flags()

        return self

    def get_attribute(self, attr_name, flagged_only = True):
        if not self.has_attribute(attr_name):
            raise AttributeError("cannot find attribute '%s'" % attr_name)
        return self._dtable[attr_name][self._flags if flagged_only else slice(None)] # slice to create data copy

    # peaks operations
    def set_peak(self, peak_index, peak_value, flagged_only = True):
        self._dtable[np.where(self._flags)[0][peak_index] if flagged_only else peak_index] = peak_value
        self.sort_peaks_order()
        self.calculate_flags()
        return self

    def get_peak(self, peak_index, flagged_only = True):
        return self._dtable[self._flags][peak_index] if flagged_only else self._dtable[peak_index]

    def insert_peak(self, peak_value):
        pid = np.where(peak_value[0] < self.mz_all)[0][0]
        self._dtable = np.insert(self._dtable, pid, peak_value)
        self._flags = np.insert(self._flags, pid, False)
        self.calculate_flags()
        return self

    def remove_peak(self, peak_index, flagged_only = True):
        if isinstance(peak_index, Iterable): peak_index = list(peak_index)
        rmid = np.where(self._flags)[0][peak_index] if flagged_only else peak_index
        self._dtable = np.delete(self._dtable, rmid)
        self._flags = np.delete(self._flags, rmid)
        if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]' % str(self.ID))
        return self

    def cleanup_unflagged_peaks(self, flag_name = None):
        if not (flag_name is None or flag_name == 'flags' or flag_name in self._flag_attrs):
            raise AttributeError('[%s] is not a flag attribution name')
        rmids = np.where((self._flags if flag_name is None else self._dtable[flag_name]) == 0)
        self._dtable = np.delete(self._dtable, rmids)
        self._flags = np.delete(self._flags, rmids)
        if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]' % str(self.ID))
        return self

    # exports
    def to_list(self):
        return zip(*self._dtable.tolist()) + [self._flags.tolist()]

    def to_dict(self, dict_type = OrderedDict):
        _conv = lambda x: (x.astype(int) if x.dtype == np.bool else x).tolist()
        return dict_type([(n, _conv(self._dtable[n])) for n in self.attributes] + [('flags', _conv(self._flags))])

    def to_str(self, delimiter = ','):
        title, data = zip(*self.to_dict().items())
        return join(map(lambda x: join(map(str, x), delimiter), [title] + zip(*data)), '\n')

    # utils
    def copy(self):
        return deepcopy(self)


