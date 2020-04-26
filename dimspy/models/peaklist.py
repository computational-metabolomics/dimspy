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


import logging
import warnings
from collections import OrderedDict
from collections.abc import Iterable
from copy import deepcopy
from typing import Callable, Sequence, Mapping, Union, Type

import numpy as np
import numpy.lib.recfunctions as rfn
import pandas as pd

from .peaklist_metadata import PeakList_Metadata
from .peaklist_tags import PeakList_Tags


class PeakList(object):
    """
    The PeakList class.

    Stores mass spectrometry peaks list data. It requires an ID, mz values, and intensities. It can store extra peak
    attributes e.g. SNRs, and peaklist tags and metadata. It utilises the automatically managed flags to "remove" or
    "retain" peaks without actually delete them. Therefore the filterings on the peaks are traceable.

    :param ID: The ID of the peaklist data, unique string or integer value is recommended
    :param mz: Mz values of all the peaks. Must in the ascending order
    :param intensity: Intensities of all the peaks. Must have the same size as mz
    :param kwargs: Key-value pairs of the peaklist metadata

    >>> mz_values = np.random.uniform(100, 1200, size = 100)
    >>> int_values = np.random.normal(60, 10, size = 100)
    >>> peaks = PeakList('dummy', mz_values, int_values, description = 'a dummy peaklist')

    Internally the peaklist data is stored by using numpy structured array namely the attribute talbe (this may change in the future):

    +-------+-----------+------+----------+-----+----------+
    | mz    | intensity | snr  | snr_flag | ... | *flags** |
    +=======+===========+======+==========+=====+==========+
    | 102.5 | 21.7      | 10.5 | True     |     | True     |
    +-------+-----------+------+----------+     +----------+
    | 111.7 | 12.3      | 5.1  | False    | ... | False    |
    +-------+-----------+------+----------+     +----------+
    | 126.3 | 98.1      | 31.7 | True     |     | True     |
    +-------+-----------+------+----------+     +----------+
    | 133.1 | 68.9      | 12.6 | True     |     | True     |
    +-------+-----------+------+----------+     +----------+
    | ...   |           |      |          |     |          |
    +-------+-----------+------+----------+-----+----------+

    Each column is called an attribute. The first two attributes are fixed as "mz" and "intensity". They cannot be added or
    removed as the others. The last "attribute" is the "flags", which is fact stored separately. The "flags" column is
    calculated automatically according to all the manually set flag attributes, e.g., the "snr_flag". It can only be changed
    by the class itself. The unflagged peaks are considered as "removed". They are kept internally mainly for visualization
    and tracing purposes.

    .. warning::
        Removing a flag attribute may change the "flags" column, and cause the unflagged peaks to be flagged again. As
        most the processes are applied only on the flagged peaks, these peaks, if the others have gone through such process,
        may have incorrect values.

        In principle, setting a flag attribute should be considered as an irreversible process.

    """

    _is_ordered = staticmethod(lambda vals: all([x[0] - x[1] >= 0 for x in zip(vals[1:], vals[:-1])]))

    def __init__(self, ID: str, mz: Sequence[float], intensity: Sequence[float], **metadata):
        if not self._is_ordered(mz):
            raise ValueError('mz values not in ascending order, check input data')
        if len(intensity) != len(mz):
            raise ValueError('mz values and intensities must have the same size')

        self._dtable = np.array(list(zip(mz, intensity)), dtype=[('mz', 'f8'), ('intensity', 'f8')])

        self._id = ID
        self._metadata = PeakList_Metadata(**metadata)

        self._tags = PeakList_Tags()
        self._flags = np.ones_like(mz, dtype=np.bool)
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

    def __repr__(self):
        return self.to_str('\t')

    def __getitem__(self, item: Union[str, int, slice, list, np.ndarray]):
        if type(item) in (int, slice, list, np.ndarray):
            return self.get_peak(item)
        else:
            return self.get_attribute(item)

    def __setitem__(self, item: Union[str, int, slice, list, np.ndarray], value):
        if type(item) in (int, slice, list, np.ndarray):
            self.set_peak(item, value)
        else:
            self.set_attribute(item, value)

    def __getattr__(self, item: str):
        if item.endswith('_all'):
            return self.get_attribute(item[:-4], flagged_only=False)
        else:
            return self.get_attribute(item, flagged_only=True)

    def __setattr__(self, item: str, value: Sequence):
        if item != '_dtable' and self.has_attribute(item):
            if item.endswith('_all'):
                return self.set_attribute(item[:-4], value, flagged_only=False)
            else:
                return self.set_attribute(item, value, flagged_only=True)
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
        """
        Property of the peaklist ID.

        :getter: Returns the peaklist ID
        :setter: Set the peaklist ID
        :type: Same as input ID

        """
        return self._id

    @ID.setter
    def ID(self, value: str):
        self._id = value

    @property
    def size(self):
        """
        Property of the peaklist size.

        :getter: Returns the flagged peaklist size
        :type: int

        """
        return np.sum(self._flags)

    @property
    def full_size(self):
        """
        Property of the peaklist full size.

        :getter: Returns the full peaklist size, i.e., including the unflagged peaks
        :type: int

        """
        return len(self._flags)

    @property
    def shape(self):
        """
        Property of the peaklist attributes table shape.

        :getter: Returns the attibutes table shape, i.e., peaks number x attributes number. The "flags" column does not count
        :type: tuple

        """
        return self.size, len(self.attributes)

    @property
    def full_shape(self):
        """
        Property of the peaklist full attributes table shape.

        :getter: Returns the full attibutes table shape, including the unflagged peaks
        :type: tuple

        """
        return self.full_size, len(self.attributes)

    @property
    def metadata(self):
        """
        Property of the peaklist metadata.

        :getter: Returns an access interface to the peaklist metadata object
        :type: PeakList_Metadata object

        """
        return self._metadata

    @property
    def tags(self):
        """
        Property of the peaklist tags.

        :getter: Returns an access interface to the peaklist tags object
        :type: PeakList_Tags object

        """
        return self._tags

    @property
    def attributes(self):
        """
        Property of the attribute names.

        :getter: Returns a tuple of the attribute names
        :type: tuple

        """
        return self._dtable.dtype.names

    @property
    def flag_attributes(self):
        """
        Property of the flag attribute names.

        :getter: Returns a tuple of the flag attribute names
        :type: tuple

        """
        return tuple(self._flag_attrs[:])

    @property
    def peaks(self):
        """
        Property of the attribute table.

        :getter: Returns a deep copy of the flagged attribute table
        :type: numpy structured array

        """
        return self._dtable[self._flags].copy()

    @property
    def flags(self):
        """
        Property of the flags.

        :getter: Returns a deep copy of the flags array
        :type: numpy array

        """
        return self._flags.copy()

    @property
    def dtable(self):
        """
        Property of the overall attribute table.

        :getter: Returns the original attribute table
        :type: numpy structured array

        .. warning::
            This property directly accesses the internal attribute table. Be careful when manipulating the data,
            particularly pay attention to the potential side-effects.

        """
        return self._dtable  # ref for faster access, using carefully

    # publics
    # manual maintainance
    # these functions are aotumatically called internally, think twice before you manually call them
    def sort_peaks_order(self):
        """
        Sorts peaklist mz values into ascending order.

        .. note::
            This method will be called automatically every time the mz values are changed.

        """
        sids = np.argsort(self._dtable['mz'])  # do not use self.mz, recursive calling
        self._dtable = self._dtable[sids]
        self._flags = self._flags[sids]

    def calculate_flags(self):
        """
        Re-calculates the flags according to the flag attributes.

        :rtype: numpy array

        .. note::
            This method will be called automatically every time a flag attribute is added, removed, or changed.

        """
        self._flags = np.ones_like(self._flags) if len(self._flag_attrs) == 0 else \
            np.sum(self._dtable[self._flag_attrs].tolist(), axis=1) == len(self._flag_attrs)
        if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]', str(self.ID))
        return self.flags

    # attribute operations
    def has_attribute(self, attr_name: str):
        """
        Checks whether there exists an attribute in the table.

        :param attr_name: The attribute name for checking
        :rtype: bool

        """
        return attr_name in self.attributes

    def add_attribute(self, attr_name: str, attr_value: Sequence, attr_dtype: Union[Type, str, None] = None,
                      is_flag: bool = False,
                      on_index: Union[int, None] = None, flagged_only: bool = True, invalid_value=np.nan):
        """
        Adds an new attribute to the PeakList attribute table.

        :param attr_name: The name of the new attribute, must be a string
        :param attr_value: The values of the new attribute. It's size must equals to PeakList.size
            (if flagged_only == True), or PeakList.full_size (if flagged_only == False)
        :param attr_dtype: The data type of the new attribute. If it is set to None, the PeakList will
            try to detect the data type based on attr_value. If the detection failed it will take the "object" type. Default = None
        :param is_flag: Whether the new attribute is a flag attribute, i.e., will be used in flags calculation. Default = False
        :param on_index: Insert the new attribute on a specific column. It can't be 0 or 1, as the first two
            attributes are fixed as mz and intensity. Setting to None means to put it to the last column. Default = None
        :param flagged_only: Whether the attr_value is set to the flagged peaks or all peaks. Default = True
        :param invalid_value: If flagged_only is set to True, this value will be assigned to the unflagged peaks.
            The actual value depends on the attribute data type. For instance, on a boolean attribute invalid_value = 0 will
            be converted to False. Default = numpy.nan
        :rtype: PeakList object (self)

        """
        if attr_name in self.__dict__:
            raise AttributeError('attribute name already been used by property')
        if attr_name in ('mz', 'intensity', 'flags'):
            raise AttributeError('cannot add reserved attribute [%s]' % attr_name)
        if self.has_attribute(attr_name):
            raise AttributeError('attribute [%s] already exists' % attr_name)
        if on_index is not None and not (-self.shape[1] + 1 < on_index < 0 or 1 < on_index < self.shape[1]):
            raise IndexError('index [%d] out of (insertable) range' % on_index)

        attr_name = str(attr_name)  # rfn.append_fields doesn't recognise unicode

        adt = bool if is_flag else \
            attr_dtype if attr_dtype is not None else \
                attr_value.dtype.str if hasattr(attr_value, 'dtype') else \
                    ('U%d' % max(map(len, attr_value))) if isinstance(attr_value[0], str) else \
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

        if is_flag and self.size > 0 and not (set(nattr) in ({0}, {1}, {0, 1})):
            raise ValueError('flag attribute can only contain True / False values')

        if on_index is None: on_index = self.shape[1]
        anames, atypes = map(list, zip(*self._dtable.dtype.descr))
        prevnm, restnm, resttp = anames[:on_index], anames[on_index:], atypes[on_index:]

        # suppress numpy's future warning regarding the structure array indexing
        # "Numpy has detected that you (may be) writing to an array returned
        #  by numpy.diagonal or by selecting multiple fields in a structured
        #  array. This code will likely break in a future numpy release ..."
        warnings.simplefilter(action='ignore', category=FutureWarning)
        prevtb = np.array(nattr, dtype=[(attr_name, adt)]) if len(prevnm) == 0 else \
            rfn.append_fields(self._fields_view(self._dtable, prevnm), attr_name, nattr, dtypes=adt, usemask=False)
        self._dtable = prevtb if len(restnm) == 0 else \
            rfn.append_fields(prevtb, restnm, list(zip(*self._fields_view(self._dtable, restnm))), dtypes=resttp,
                              usemask=False)
        warnings.resetwarnings()

        if is_flag:
            self._flag_attrs += [attr_name]
            self.calculate_flags()

        return self

    def drop_attribute(self, attr_name: str):
        """
        Drops an existing attribute.

        :param attr_name: The attribute name to drop. It cannot be mz, intensity, or flags
        :rtype: PeakList object (self)

        """
        if attr_name in ('mz', 'intensity', 'flags'):
            raise AttributeError('cannot drop reserved attribute [%s]' % attr_name)
        if not self.has_attribute(attr_name):
            raise AttributeError('attribute [%s] does not exist' % attr_name)

        self._dtable = self._dtable[[x for x in self.attributes if x != attr_name]]

        if attr_name in self._flag_attrs:
            logging.warning('flags recalculated, unflagged peaks may contain incorrect values')
            self._flag_attrs = [x for x in self._flag_attrs if x != attr_name]
            self.calculate_flags()

        return self

    def set_attribute(self, attr_name: str, attr_value: Sequence, flagged_only: bool = True, unsorted_mz: bool = False):
        """
        Sets values to an existing attribute.

        :param attr_name: The attribute to set values
        :param attr_value: The new attribute values, It's size must equals to PeakList.size
            (if flagged_only == True), or PeakList.full_size (if flagged_only == False)
        :param flagged_only: Whether the attr_value is set to the flagged peaks or all peaks. Default = True
        :param unsorted_mz: Whether the attr_value contains unsorted mz values. This parameter is valid only when
            attr_name == "mz". Default = False
        :rtype: PeakList object (self)

        """
        if attr_name == 'flags':
            raise AttributeError('cannot assign read-only attribute [flag]')
        if not self.has_attribute(attr_name):
            raise AttributeError('attribute [%s] does not exist' % attr_name)
        if attr_name == 'mz' and not (unsorted_mz or self._is_ordered(attr_value)):
            raise ValueError('attribute [mz] not in ascending order')
        if attr_name != 'mz' and unsorted_mz == True:
            logging.warning('setting unsorted_mz flag for non-mz attribute, ignore')

        # May raise FutureWarning, but that shold be a false alarm, because attr_name is a string rather than a list
        # we are not accessing multiple fields by their names here
        warnings.simplefilter(action='ignore', category=FutureWarning)
        self._dtable[attr_name][self._flags if flagged_only else slice(None)] = attr_value
        warnings.resetwarnings()

        if attr_name == 'mz' and unsorted_mz: self.sort_peaks_order()
        if attr_name in self._flag_attrs: self.calculate_flags()
        return self

    def get_attribute(self, attr_name: str, flagged_only: bool = True):
        """
        Gets values of an existing attribute.

        :param attr_name: The attribute to get values
        :param flagged_only: Whether to return the values of flagged peaks or all peaks. Default = True
        :rtype: numpy array

        """
        if not self.has_attribute(attr_name):
            raise AttributeError("cannot find attribute '%s'" % attr_name)
        return self._dtable[attr_name][self._flags if flagged_only else slice(None)]  # slice to create data copy

    # peaks operations
    def set_peak(self, peak_index: int, peak_value: Sequence, flagged_only: bool = True):
        """
        Sets values to a peak.

        :param peak_index: The index of the peak to set values
        :param peak_value: The new peak values. Must contain values for all the attributes (not including flags)
        :param flagged_only: Whether the peak_value is set to the index of flagged peaks or all peaks. Default = True
        :rtype: PeakList object (self)

        >>> print(peaks)
        mz, intensity, snr, flags
        10, 10, 10, True
        20, 20, 20, True
        30, 30, 30, False
        40, 40, 40, True
        >>> print(peaks.set_peak(2, [50, 50, 50], flagged_only = True))
        mz, intensity, snr, flags
        10, 10, 10, True
        20, 20, 20, True
        30, 30, 30, False
        50, 50, 50, True
        >>> print(peaks.set_peak(2, [40, 40, 40], flagged_only = False))
        mz, intensity, snr, flags
        10, 10, 10, True
        20, 20, 20, True
        40, 40, 40, False
        50, 50, 50, True

        """
        if isinstance(peak_index, Iterable): peak_index = list(peak_index)
        self._dtable[np.where(self._flags)[0][peak_index] if flagged_only else peak_index] = peak_value
        self.sort_peaks_order()
        self.calculate_flags()
        return self

    def get_peak(self, peak_index: Union[int, Sequence[int]], flagged_only: bool = True):
        """
        Gets values of a peak.

        :param peak_index: The index of the peak to get values
        :param flagged_only: Whether the values are taken from the index of flagged peaks or all peaks. Default = True
        :rtype: numpy array

        """
        if isinstance(peak_index, Iterable): peak_index = list(peak_index)
        return self._dtable[self._flags][peak_index] if flagged_only else self._dtable[peak_index]

    def insert_peak(self, peak_value: Sequence):
        """
        Insert a new peak.

        :param peak_value: The values of the new peak. Must contain values for all the attributes. It's position depends
            on the mz value, i.e., the 1st value of the input
        :rtype: PeakList object (self)

        """
        pid = np.where(peak_value[0] < self.mz_all)[0][0]
        self._dtable = np.insert(self._dtable, pid, peak_value)
        self._flags = np.insert(self._flags, pid, False)
        self.calculate_flags()
        return self

    def remove_peak(self, peak_index: Union[int, Sequence[int]], flagged_only: bool = True):
        """
        Remove an existing peak.

        :param peak_index: The index of the peak to remove
        :param flagged_only: Whether the index is for flagged peaks or all peaks. Default = True
        :rtype: PeakList object (self)

        """
        if isinstance(peak_index, Iterable): peak_index = list(peak_index)
        rmid = np.where(self._flags)[0][peak_index] if flagged_only else peak_index
        self._dtable = np.delete(self._dtable, rmid)
        self._flags = np.delete(self._flags, rmid)
        if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]', str(self.ID))
        return self

    def cleanup_unflagged_peaks(self, flag_name: Union[str, None] = None):
        """
        Remove unflagged peaks.

        :param flag_name: Remove peaks unflagged by this flag attribute. Setting None means to remove peaks unflagged by
            the overall flags. Default = None
        :rtype: PeakList object (self)

        >>> print(peaks)
        mz, intensity, intensity_flag, snr, snr_flag, flags
        10, 70, True, 10, False, False
        20, 60, True, 20, True, True
        30, 50, False, 30, True, False
        40, 40, False, 40, True, False
        >>> print(peaks.cleanup_unflagged_peaks('snr_flag'))
        mz, intensity, intensity_flag, snr, snr_flag, flags
        20, 60, True, 20, True, True
        30, 50, False, 30, True, False
        40, 40, False, 40, True, False
        >>> print(peaks.cleanup_unflagged_peaks())
        mz, intensity, intensity_flag, snr, snr_flag, flags
        20, 60, True, 20, True, True

        """
        if not (flag_name is None or flag_name == 'flags' or flag_name in self._flag_attrs):
            raise AttributeError('[%s] is not a flag attribution name')
        rmids = np.where((self._flags if flag_name is None else self._dtable[flag_name]) == 0)
        self._dtable = np.delete(self._dtable, rmids)
        self._flags = np.delete(self._flags, rmids)
        if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]', str(self.ID))
        return self

    # exports
    def to_list(self):
        """
        Exports peaklist attribute table to a list, including the flags.

        :rtype: list

        """
        return list(zip(*self._dtable.tolist())) + [self._flags.tolist()]

    def to_dict(self, dict_type: Callable[[Sequence], Mapping] = OrderedDict) -> Mapping:
        """
        Exports peaklist attribute table to a dictionary (mappable object), including the flags.

        :param dict_type: Result dictionary type, Default = OrderedDict
        :rtype: list

        """
        _conv = lambda x: (x.astype(int) if x.dtype == np.bool else x).tolist()
        return dict_type([(n, _conv(self._dtable[n])) for n in self.attributes] + [('flags', _conv(self._flags))])

    def to_str(self, delimiter: str = ','):
        """
        Exports peaklist attribute table to a string, including the flags. It can also be used inexplicitly.

        :rtype: str

        """
        title, data = zip(*self.to_dict().items())
        return str.join('\n', [str.join(delimiter, map(str, x)) for x in [title] + list(zip(*data))])

    def to_df(self):
        """
        Exports peaklist attribute table to Pandas DataFrame, including the flags.

        :rtype: pd.DataFrame

        """
        title, data = zip(*self.to_dict().items())
        return pd.DataFrame(list(zip(*data)), columns=title)

    # utils
    def copy(self):
        """
        Returns a deep copy of the peaklist.

        :rtype: PeakList object

        """
        return deepcopy(self)
