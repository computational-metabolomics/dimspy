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


# internal classes
# -----------------------------------------------------------
# Internal use only. DO NOT try metadata.metadata.attr. All attribute methods overrided
class _Metadata(dict):
    def __init__(self, *args, **kwargs):
        super(_Metadata, self).__init__(*args, **kwargs)

    def __getattr__(self, item):
        if self.has_key(item):
            return self[item]
        else:
            super(_Metadata, self).__getattribute__(item)  # TODO: fix no_such_method_in_root_class after txt loader

    def __setattr__(self, item, value):
        assert item != '__dict__', '"__dict__" is not an acceptable metadata key'
        assert type(value) is not _Metadata, 'metadata object is not an acceptable metadata value'
        if not self.__dict__.has_key(item):
            self[item] = value
        else:
            super(_Metadata, self).__setattr__(item, value)

    def __delattr__(self, item):
        if self.has_key(item):
            del self[item]
        else:
            super(_Metadata, self).__delattr__(item)


class _Tags(object):
    def __init__(self, *args, **kwargs):
        self._tags = {'_untyped': [], '_typed': {}}
        self.add_tags(*args, **kwargs)

    # build-ins
    def __str__(self):
        return str(self.to_list())[1:-1]

    # properties
    @property
    def tag_types(self):
        return self._tags['_typed'].keys()

    @property
    def tag_values(self):
        return self._tags['_untyped'] + self._tags['_typed'].values()

    # methods
    def tag_of(self, tag_type=None):
        assert tag_type is None or self.has_tag_type(tag_type), 'unknown tag type [%s]' % tag_type
        return self._tags['_untyped'] if tag_type is None else self._tags['_typed'][tag_type]

    def has_tags(self, *args, **kwargs):
        assert None not in args or None not in kwargs.values(), '[None] is an acceptable tag value'
        assert not any(map(lambda x: x in kwargs.values(), args)), 'searching tags have duplication'
        return all(map(lambda x: x in self.tag_values, args)) and \
               all(map(lambda x: self._tags['_typed'].get(x[0], None) == x[1], kwargs.items()))

    def has_any_tags(self, *args, **kwargs):
        assert None not in args or None not in kwargs.values(), '[None] is an acceptable tag value'
        assert not any(map(lambda x: x in kwargs.values(), args)), 'searching tags have duplication'
        return any(map(lambda x: x in self.tag_values, args)) or \
               any(map(lambda x: self._tags['_typed'].get(x[0], None) == x[1], kwargs.items()))

    def has_tag_type(self, tag_type):
        return tag_type in self.tag_types

    def add_tags(self, *args, **kwargs):
        assert 'NA' not in kwargs.keys(), '[NA] is not an acceptable tag type' # reserve for hdf5 protal
        assert None not in args or None not in kwargs.values(), '[None] is not an acceptable tag value'
        assert not any(map(lambda x: x in self.tag_values, list(args) + kwargs.values())), 'tag already exists'
        self._tags['_untyped'] += args
        self._tags['_typed'].update(kwargs)

    def drop_tags(self, *args):
        self._tags['_untyped'] = filter(lambda x: x not in args, self._tags['_untyped'])
        self._tags['_typed'] = dict(filter(lambda x: x[1] not in args, self._tags['_typed'].items()))

    def drop_tag_types(self, *args):
        self._tags['_typed'] = dict(filter(lambda x: x[0] not in args, self._tags['_typed'].items()))

    def drop_all_tags(self):
        self._tags = {'_untyped': [], '_typed': {}}

    def to_list(self):
        return self._tags['_untyped'] + self._tags['_typed'].items()


# public classes
# -----------------------------------------------------------
class PeakList(object):
    def __init__(self, ID, mz, intensity, **metadata):
        assert all(
            map(lambda x: x[0] - x[1] >= 0, zip(mz[1:], mz[:-1]))), 'mz values not in ascending order, check input data'
        assert len(intensity) == len(mz), 'mz values and intensities must have the same size'
        self.__dict__['_dtable'] = np.array(zip(mz, intensity), dtype=[('mz', 'f8'), ('intensity', 'f8')])

        self._id = str(ID)
        self._tags = _Tags()
        self._metadata = _Metadata(**metadata)

        self._flags = np.ones_like(mz, dtype=np.bool)
        self._flag_attrs = []

    # built-ins
    def __len__(self):
        return self.size

    def __str__(self):
        return self.to_str(', ')

    def __getitem__(self, item):
        if type(item) in (int, slice):
            return self.get_peak(item)
        else:
            return self.get_attribute(item)

    def __setitem__(self, item, value):
        if type(item) in (int, slice):
            self.set_peak(item, value)
        else:
            self.set_attribute(item, value)

    def __getattr__(self, item):
        if item.endswith('_all'):
            return self.get_attribute(item[:-4], flagged_only=False)
        else:
            return self.get_attribute(item, flagged_only=True)

    def __setattr__(self, item, value):
        try:
            return self.set_attribute(item[:-4] if item.endswith('_all') else item, value, flagged_only=True)
        except AttributeError:
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
    def tags(self):
        return self._tags

    @property
    def tag_types(self):
        return self._tags.tag_types

    @property
    def tag_values(self):
        return self._tags.tag_values

    @property
    def metadata(self):
        return self._metadata

    @property
    def attributes(self):
        return self._dtable.dtype.names

    @property
    def flag_attributes(self):
        return self._flag_attrs[:] 

    @property
    def peaks(self):
        return self._dtable[self._flags].copy()

    @property
    def flags(self):
        return self._flags.copy()

    @property
    def dtable(self):
        return self._dtable.copy()

    # publics
    # manual maintainance
    def sort_peaks_order(self):
        sids = np.argsort(self._dtable['mz'])  # do not use self.mz, recursive calling
        self._dtable = self._dtable[sids]
        self._flags = self._flags[sids]

    def calculate_flags(self):
        if len(self._flag_attrs) == 0:
            self._flags = np.ones_like(self._flags)
        else:
            self._flags = np.sum(self._dtable[self._flag_attrs].tolist(), axis=1) == len(self._flag_attrs)
            if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]' % str(self.ID))
        return self._flags

    # attribute operations
    def has_attribute(self, attr_name):
        return attr_name in self.attributes

    def add_attribute(self, attr_name, attr_value, attr_dtype=None, is_flag=False, flagged_only=True,
                      invalid_value=np.nan):
        assert attr_name not in self.__dict__.keys(), 'attribute name already been used by property'
        assert attr_name not in ('mz', 'intensity', 'flags'), 'cannot add reserved attribute [%s]' % attr_name
        assert not self.has_attribute(attr_name), 'attribute [%s] already exists' % attr_name
        if is_flag and self.size > 0: assert set(attr_value) in (
        {0}, {1}, {0, 1}), 'flag attribute can only contain True/False values'

        adt = bool if is_flag else \
              attr_dtype if attr_dtype is not None else \
              attr_value.dtype.str if hasattr(attr_value, 'dtype') else type(attr_value[0])
        if adt in (bool, 'bool', '|b1'): adt = 'b'  # fix numpy dtype bug

        nattr = np.array([invalid_value] * self._dtable.shape[0]).astype([(attr_name, adt)])
        nattr[attr_name][self._flags if flagged_only else slice(None)] = attr_value
        self._dtable = rfn.merge_arrays((self._dtable, nattr), flatten=True)

        if not is_flag:
            return
        self._flag_attrs += [attr_name]
        self.calculate_flags()

        return self

    def drop_attribute(self, attr_name):
        assert attr_name not in ('mz', 'intensity', 'flags'), 'cannot drop reserved attribute [%s]' % attr_name
        assert self.has_attribute(attr_name), 'attribute [%s] not exists' % attr_name

        self._dtable = self._dtable[list(filter(lambda x: x != attr_name, self.attributes))]

        if attr_name not in self._flag_attrs: return
        logging.warning('flags recalculated, unflagged peaks may contain incorrect values')
        self._flag_attrs = filter(lambda x: x != attr_name, self._flag_attrs)
        self.calculate_flags()

        return self

    def set_attribute(self, attr_name, attr_value, flagged_only=True):
        assert attr_name not in ('flags',), 'cannot assign read-only attribute [%s]' % attr_name
        if not self.has_attribute(attr_name):
            raise AttributeError("cannot find attribute '%s'" % attr_name)

        if flagged_only:
            self._dtable[attr_name][self._flags] = attr_value
        else:
            self._dtable[attr_name] = attr_value

        if attr_name == 'mz':
            self.sort_peaks_order()
        if attr_name in self._flag_attrs:
            self.calculate_flags()

    def get_attribute(self, attr_name, flagged_only=True):
        if not self.has_attribute(attr_name):
            raise AttributeError("cannot find attribute '%s'" % attr_name)
        return self._dtable[attr_name][self._flags] if flagged_only else self._dtable[attr_name]

    # peaks operations
    def insert_peak(self, peak_value):
        pid = np.where(peak_value[0] < self.mz_all)[0][0]
        self._dtable = np.insert(self._dtable, pid, peak_value)
        self._flags = np.insert(self._flags, pid, False)
        self.calculate_flags()

    def remove_peak(self, peak_index):
        rmid = np.where(self._flags)[0][list(peak_index) if isinstance(peak_index, Iterable) else peak_index]
        self._dtable = np.delete(self._dtable, rmid)
        self._flags = np.delete(self._flags, rmid)
        if self.size == 0: logging.warning('all peaks are removed for peaklist [%s]' % str(self.ID))

    def set_peak(self, peak_index, peak_value, flagged_only=True):
        if flagged_only:
            self._dtable[np.where(self._flags)[0][peak_index]] = peak_value
        else:
            self._dtable[peak_index] = peak_value

        self.sort_peaks_order()
        self.calculate_flags()

    def get_peak(self, peak_index, flagged_only=True):
        return self._dtable[self._flags][peak_index] if flagged_only else self._dtable[peak_index]

    # tags operations (wrapper of _Tags)
    def tag_of(self, tag_type=None):
        return self._tags.tag_of(tag_type)

    def has_tags(self, *args, **kwargs):
        return self._tags.has_tags(*args, **kwargs)

    def has_any_tags(self, *args, **kwargs):
        return self._tags.has_any_tags(*args, **kwargs)

    def has_tag_type(self, tag_type):
        return self._tags.has_tag_type(tag_type)

    def add_tags(self, *args, **kwargs):
        return self._tags.add_tags(*args, **kwargs)

    def drop_tags(self, *args):
        return self._tags.drop_tags(*args)

    def drop_tag_types(self, *args):
        return self._tags.drop_tag_types(*args)

    def drop_all_tags(self):
        return self._tags.drop_all_tags()

    # exports
    def to_list(self):
        return zip(*self._dtable.tolist()) + [self._flags.tolist()]

    def to_dict(self, dict_type=OrderedDict):
        _conv = lambda x: (x.astype(int) if x.dtype == np.bool else x).tolist()
        return dict_type([(n, _conv(self._dtable[n])) for n in self.attributes] + [('flags', _conv(self._flags))])

    def to_str(self, delimiter=','):
        title, data = zip(*self.to_dict().items())
        return join(map(lambda x: join(map(str, x), delimiter), [title] + zip(*data)), '\n')

    # utils
    def copy(self):
        return deepcopy(self)


# test
if __name__ == '__main__':
    # generate random peak
    size = 100
    mzs = sorted(np.random.uniform(100, 1000, size=size))
    ints = np.abs(np.random.normal(10, 3, size=size))
    snr = np.abs(np.random.normal(1000, 400, size=size))

    pl = PeakList('sample_peaklist', mzs, ints, mz_range=(100, 1000))
    pl.add_tags('sample', 'passed_qc', main_class='treatment_1')
    pl.add_attribute('snr', snr)
    pl.metadata.type = 'blank'

    import pdb; pdb.set_trace()

    import cPickle as cp

    s = cp.dumps(pl)
    pl = cp.loads(s)

    # magic functions
    print len(pl)
    print pl

    # tags
    print pl.tags
    print pl.tag_types
    print pl.tag_values
    pl.add_tags('not_QC')
    pl.add_tags(time_point='1hr')
    print pl.tag_of('time_point')
    print pl.has_tags('1hr')
    print pl.has_any_tags(main_class='1hr')
    print pl.has_tag_type('main_class')
    pl.drop_tags('not_QC')
    pl.drop_tag_types('time_point')
    print pl.tags.to_list()

    # metadata
    pl.metadata['replicate'] = '4'
    pl.metadata.sample = 'S01'
    print pl.metadata
    del pl.metadata.sample
    print pl.metadata
    print pl.metadata.keys()

    # add flag
    print pl.flags
    pl.add_attribute('odd_flags', np.arange(len(pl)) % 2 == 0, is_flag=True)
    print pl.flags

    # attrs access
    print pl.mz  # internally == getAttribute
    print pl.snr
    print pl.mz_all
    print pl.snr_all

    pl.mz = np.random.uniform(100, 1000, size=len(pl))  # internally == setAttribute
    print pl[:10]

    print pl['mz']  # internally == getPeak
    print pl['snr']

    print pl[0]  # internally == getPeak
    pl[0] = (0, 100, 200, True)  # internally == setPeak
    print pl[:10]

    pl[0][0] = 1  # not working
    print pl[:10]

    pl['mz'][0] = 2  # not working
    print pl[:10]

    # properties
    print pl.ID
    print pl.size
    print pl.shape
    print pl.attributes
    print pl.peaks
    print pl.dtable

    # publics
    print pl.has_attribute('snr')
    print pl.has_attribute('no_such_thing')
    pl.add_attribute('snr_filter_flag', pl.snr > 20, is_flag=True)  # must be all full size vector
    print pl
    pl.drop_attribute('snr_filter_flag')
    print pl

    pl.insert_peak((900, 20, 30, True))  # auto sort
    print pkl
    pkl.remove_peak(-2)
    print pkl

    print pkl.to_list()[0]
    print pkl.to_dict().keys()
