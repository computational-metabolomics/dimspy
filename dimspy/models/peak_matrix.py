#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
peaksAlignment: the PeaksAlignment class.

author(s): Albert Zhou, Ralf Weber
origin: Nov. 10, 2016

"""

import logging
import numpy as np
from string import join
from peaklist import PeakList


# np.average accepts weights, set to 0 to exclude value
def _masknan(x, axis, _func):
    bwt = x.astype(bool)
    nid = (np.sum(~np.isclose(bwt, 0), axis=axis) == 0)
    bwt[:, nid] = 1
    val = _func(bwt)
    val[nid] = np.nan  # all values = 0 -> nanzero average / std = nan
    return val


def _nzmean(x, axis):
    _func = lambda b: np.average(x, axis=axis, weights=b)
    return _masknan(x, axis, _func)


def _nzstd(x, axis):
    _func = lambda b: np.sqrt(np.average(np.power(x - np.average(x, axis=axis, weights=b), 2), axis=axis, weights=b))
    return _masknan(x, axis, _func)


class PeakMatrix(object):
    def __init__(self, peaklist_ids, peaklist_tags, peaklist_attributes, mask=None):
        assert len(peaklist_ids) == len(peaklist_tags) and \
               all(map(lambda x: len(peaklist_ids) == x.shape[0],
                       peaklist_attributes.values())), 'alignment input data shape not match'
        assert 'mz' in peaklist_attributes.keys() and 'intensity' in peaklist_attributes.keys(), 'required attribute fields not available'

        self._pids = np.array(peaklist_ids)
        self._tags = np.array(peaklist_tags)
        self._attr_dict = {k: np.array(v) for k,v in peaklist_attributes.items()}
        self.mask = mask

        if self.is_empty(): logging.warning('alignment input data is empty')

    # built-ins
    def __len__(self):
        return self.shape

    def __str__(self):
        return self.to_str(',')

    # properties
    @property
    def mask(self):  # mask is different from flags in that they will not be masked by ifself, i.e., more temporary
        return self._mask.copy()

    @mask.setter
    def mask(self, value):
        self._mask = np.zeros_like(self._pids, dtype=bool)
        self._mask[slice(None) if value is None else value] = True

    @property
    def peaklist_ids(self):
        return self._pids[self._mask]

    @property
    def peaklist_tags(self):
        return self._tags[self._mask]

    @property
    def peaklist_tag_types(self):
        return list(set(reduce(lambda x, y: x + y, [t.tag_types for t in self.peaklist_tags])))

    @property
    def peaklist_tag_values(self):
        return list(set(reduce(lambda x, y: x + y, [t.tag_values for t in self.peaklist_tags])))

    @property
    def attributes(self):
        return self._attr_dict.keys()

    @property
    def shape(self):
        return np.sum(self._mask), self._attr_dict['mz'].shape[1]

    @property
    def full_shape(self):
        return self._attr_dict['mz'].shape

    @property
    def present(self):
        return np.sum(self.mz_matrix > 0, axis=0)

    @property
    def missing(self):
        return np.sum(self.mz_matrix == 0, axis=1)

    @property
    def rsd(self):
        assert self.shape[0] > 1, 'calculating RSD on less than 2 samples'
        rsd = (lambda m: _nzstd(m, 0) / _nzmean(m, 0) * 100)(self.intensity_matrix)
        rsd[np.where(map(lambda x: len(set(x[np.nonzero(x)])) == 1, self.intensity_matrix.T))] = np.nan
        return rsd

    @property
    def mz_matrix(self):
        return self.attr_matrix('mz')

    @property
    def intensity_matrix(self):
        return self.attr_matrix('intensity')

    @property
    def mz_mean_vector(self):
        return self.attr_mean_vector('mz')

    @property
    def ints_mean_vector(self):
        return self.attr_mean_vector('intensity')

    # publics
    def tags_of(self, tag_type=None):
        assert tag_type is None or all(
            map(lambda x: x.has_tag_type(tag_type), self._tags)), 'not all samples has tag type [%s]' % tag_type
        tlst = map(lambda x:
            x.tag_of(tag_type), self.peaklist_tags)
        if tag_type is None:
            tlst = reduce(lambda x, y: x + y, tlst)
        return list(set(tlst))

    def mask_tags(self, *args, **kwargs):  # any
        if kwargs.has_key('override'):
            override = kwargs['override']
            del kwargs['override']
        else:
            override = False
        mask = map(lambda x: x.has_any_tags(*args, **kwargs), self._tags)
        self.mask = np.logical_and(True if override else self._mask, mask)
        return self

    def unmask_tags(self, *args, **kwargs):  # all
        if kwargs.has_key('override'):
            override = kwargs['override']; del kwargs['override']
        else:
            override = False
        mask = map(lambda x: not x.has_tags(*args, **kwargs), self._tags)
        self.mask = np.logical_and(True if override else self._mask, mask)
        return self

    def attr_matrix(self, attr_name, masked_only=True):
        return self._attr_dict[attr_name][self._mask if masked_only else slice(None)]

    def attr_mean_vector(self, attr_name, masked_only=True):
        return _nzmean(self.attr_matrix(attr_name, masked_only), 0)

    def to_peaklist(self, ID):
        return PeakList(ID, self.mz_mean_vector, self.ints_mean_vector,
                        rsd=self.rsd, present=self.present, missing=self.missing, aligned_ids=self.peaklist_ids)

    def remove_samples(self, indices, remove_empty_peaks=True, masked_only=True):
        ids = np.arange(self._pids.shape[0])[self.mask][indices] if masked_only else indices
        self._pids = np.delete(self._pids, ids, axis=0)
        self._attr_dict = {k: np.delete(v, ids, axis=0) for k, v in self._attr_dict.items()}
        self._tags = np.delete(self._tags, ids, axis=0)
        self._mask = np.delete(self._mask, ids, axis=0)
        if remove_empty_peaks: self.remove_peaks(np.where(np.sum(self.mz_matrix, axis=0) == 0), False)
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def remove_peaks(self, indices, remove_empty_samples=True):
        self._attr_dict = {k: np.delete(v, indices, axis=1) for k, v in self._attr_dict.items()}
        if remove_empty_samples: self.remove_samples(np.where(np.sum(self.mz_matrix, axis=1) == 0), False, False)
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def is_empty(self):
        return 0 in self.shape

    # exports
    def to_str(self, delimiter='\t', transpose=False):
        dm = [[''] + map(str, self.mz_mean_vector)] # + ['']]
        dm += zip(*([map(str, self.peaklist_ids)] +
                    [map(str, ln) for ln in self.intensity_matrix.T])) #+
                    #[map(str, self.peaklist_tags)]))
        if transpose: dm = zip(*dm)
        return join(map(lambda x: join(x, delimiter), dm), '\n')


# with statements
class mask_peakmatrix:
    def __init__(self, pm, *args, **kwargs):
        self._pm = pm
        self._utags = args
        self._ttags = kwargs
        if not self._ttags.has_key('override'): self._ttags['override'] = True  # default for with statement
        self._oldmask = pm.mask

    def __enter__(self):
        self._pm.mask_tags(*self._utags, **self._ttags)
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = self._oldmask


class unmask_peakmatrix:
    def __init__(self, pm, *args, **kwargs):
        self._pm = pm
        self._utags = args
        self._ttags = kwargs
        if not self._ttags.has_key('override'): self._ttags['override'] = True  # default for with statement
        self._oldmask = pm.mask

    def __enter__(self):
        self._pm.unmask_tags(*self._utags, **self._ttags)
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = self._oldmask

class unmask_all_peakmatrix:
        def __init__(self, pm):
            self._pm = pm
            self._oldmask = pm.mask

        def __enter__(self):
            self._pm.mask  = None
            return self._pm

        def __exit__(self, exc_type, exc_val, exc_tb):
            self._pm.mask = self._oldmask


# test
if __name__ == '__main__':
    # generate random peaklists
    _mzs = lambda: sorted(np.random.uniform(100, 1000, size=100))
    _ints = lambda: np.abs(np.random.normal(10, 3, size=100))

    pkls = [
        PeakList('sample_1_1', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('sample_1_2', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('QC_1', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('sample_2_1', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('sample_2_2', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('QC_2', _mzs(), _ints(), mz_range=(100, 1000)),
    ]

    pkls[0].add_tags('sample', treatment='compound_1', time_point='1hr', plate=1)
    pkls[1].add_tags('sample', treatment='compound_1', time_point='6hr', plate=1)
    pkls[2].add_tags('qc', plate=1)
    pkls[3].add_tags('sample', treatment='compound_2', time_point='1hr', plate=2)
    pkls[4].add_tags('sample', treatment='compound_2', time_point='6hr', plate=2)
    pkls[5].add_tags('qc', plate=2)

    # create matrix
    pm = PeakMatrix(
        [p.ID for p in pkls],
        [p.tags for p in pkls],
        {a: np.vstack([p.get_attribute(a) for p in pkls]) for a in pkls[0].attributes}
    )

    import pdb;
    pdb.set_trace()

    # properties
    print pm.mask
    print pm.peaklist_ids
    print pm.peaklist_tag_types
    print pm.peaklist_tag_values
    print pm.attributes
    print pm.shape
    print pm.present
    print pm.missing
    print pm.rsd

    # mask tags
    print pm.tags_of('plate')
    print pm.tags_of(None)  # special usage, use only when you know what you are doing
    print pm.mz_matrix
    # pm.mask_tags('sample').unmask_tags(treatment = 'compound_1')
    with mask_peakmatrix(pm, 'sample'):
        print pm.mask
        print pm.peaklist_ids
        pm.remove_peaks(np.where(pm.rsd > 30))
        print pm.mask
    print pm.mask

    pm.remove_samples(np.where(pm.peaklist_ids == 'sample_2_2')[0])
    print pm.peaklist_ids
    pm.remove_peaks([0, 1, 2])
    print pm.mz_matrix

    pm.mask = None
    print pm.peaklist_ids

    # export
    print pm.to_peaklist('merged')
    print pm.to_str(transpose=True)
