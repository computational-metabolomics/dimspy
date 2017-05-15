#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
peaksAlignment: the PeaksAlignment class.

author(s): Albert Zhou, Ralf Weber
origin: Nov. 10, 2016

"""


from __future__ import division

import logging
import numpy as np
from string import join
from peaklist import PeakList


class PeakMatrix(object):
    def __init__(self, peaklist_ids, peaklist_tags, **peaklist_attributes):
        if len(peaklist_attributes) < 2 or not (peaklist_attributes.has_key('mz') and peaklist_attributes.has_key('intensity')):
            raise ValueError('required attribute fields "mz" and "intensity" not available')

        pnum, psize = peaklist_attributes.values()[0].shape
        if not len(peaklist_ids) == len(peaklist_tags) == pnum:
            raise ValueError('input peaklists number not match')
        if not all(map(lambda x: x.shape == (pnum, psize), peaklist_attributes.values())):
            raise ValueError('input attribute matrix shape not match')

        self._pids = np.array(peaklist_ids)
        self._tags = np.array(peaklist_tags)
        self._attr_dict = {k: np.array(v) for k,v in peaklist_attributes.items()}
        self._mask = np.ones(pnum, dtype = bool)

    # privates
    # np.average accepts weights, set to 0 to exclude value
    @staticmethod
    def _masknan(x, axis, _func):
        bwt = x.astype(bool)
        nid = (np.sum(~np.isclose(bwt, 0), axis = axis) == 0)
        bwt[:, nid] = 1
        val = _func(bwt)
        val[nid] = np.nan  # all values = 0 -> nanzero average / std = nan
        return val

    def _nzmean(self, x, axis):
        _func = lambda b: np.average(x, axis = axis, weights = b)
        return self._masknan(x, axis, _func)

    def _nzstd(self, x, axis):
        _func = lambda b: np.sqrt(np.average(np.power(x - np.average(x, axis = axis, weights = b), 2), axis = axis, weights = b))
        return self._masknan(x, axis, _func)

    # built-ins
    def __len__(self):
        return self.shape[0]

    def __str__(self):
        return self.to_str(delimiter = ', ', transpose = True, extend = False)

    # properties
    @property
    def mask(self):  # mask is different from flags in that they will not be masked by themselves, i.e., more temporary
        return self._mask.copy()

    @mask.setter
    def mask(self, value):
        if value is None:
            self._mask = np.ones(self.full_shape[0], dtype = bool)
        else:
            self._mask = np.zeros(self.full_shape[0], dtype = bool)
            self._mask[value] = True

    @property
    def peaklist_ids(self):
        return tuple(self._pids[self._mask])

    @property
    def peaklist_tags(self):
        return tuple(self._tags[self._mask])

    @property
    def peaklist_tag_types(self):
        return tuple(set(reduce(lambda x, y: x + y, [t.tag_types for t in self.peaklist_tags], ())))

    @property
    def peaklist_tag_values(self):
        return tuple(set(reduce(lambda x, y: x + y, [t.tag_values for t in self.peaklist_tags], ())))

    @property
    def attributes(self):
        return tuple(self._attr_dict.keys())

    @property
    def shape(self):
        return np.sum(self._mask), self._attr_dict['mz'].shape[1]

    @property
    def full_shape(self):
        return self._attr_dict['mz'].shape

    @property
    def present(self):
        return np.sum(self.intensity_matrix > 0, axis = 0)

    @property
    def fraction(self):
        return self.present / self.shape[0]

    @property
    def missing_values(self):
        return np.sum(self.intensity_matrix == 0, axis = 1)

    @property
    def rsd(self):
        if self.shape[0] < 2:
            #logging.warning('calculating RSD on less than 2 samples')
            return np.ones(self.shape[1]) * np.nan
        rsd = (lambda m: self._nzstd(m, 0) / self._nzmean(m, 0) * 100)(self.intensity_matrix)
        rsd[np.where(map(lambda x: len(set(x[np.nonzero(x)])) == 1, self.intensity_matrix.T))] = np.nan
        return rsd

    @property
    def occurance(self):
        if not self._attr_dict.has_key('intra_count'):
            logging.warning("attribute matrix ['intra_count'] not available")
            return np.ones(self.shape[1])
        return np.sum(self.attr_matrix('intra_count'), axis = 0)

    @property
    def impure(self):
        if not self._attr_dict.has_key('intra_count'):
            logging.warning("attribute matrix ['intra_count'] not available")
            return np.zeros(self.shape[1])
        return np.sum(self.attr_matrix('intra_count') > 1, axis = 0) / self.shape[0]

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
    # tags and mask
    def tags_of(self, tag_type = None):
        if tag_type is not None and not all(map(lambda x: x.has_tag_type(tag_type), self.peaklist_tags)):
            raise ValueError('not all samples has tag type [%s]' % tag_type)
        tlst = [t.tag_of(tag_type) for t in self.peaklist_tags]
        if tag_type is None: tlst = reduce(lambda x, y: x + y, tlst)
        return tuple(set(tlst))

    def mask_tags(self, *args, **kwargs):  # match to any
        override = kwargs.pop('override') if kwargs.has_key('override') else False
        mask = map(lambda x: any(map(lambda t: x.has_tag(t), args)) or
                             any(map(lambda t: x.has_tag(**dict([t])), kwargs.items())), self._tags)
        self.mask = np.logical_and(True if override else self._mask, mask)
        return self

    def unmask_tags(self, *args, **kwargs):  # match to all
        override = kwargs.pop('override') if kwargs.has_key('override') else False
        mask = map(lambda x: not (all(map(lambda t: x.has_tag(t), args)) and
                                  all(map(lambda t: x.has_tag(**dict([t])), kwargs.items()))), self._tags)
        self.mask = np.logical_and(True if override else self._mask, mask)
        return self

    # access
    def attr_matrix(self, attr_name):
        if not self._attr_dict.has_key(attr_name):
            raise KeyError('attribute matrix [%s] not available' % attr_name)
        return self._attr_dict[attr_name][self._mask]

    def attr_mean_vector(self, attr_name):
        return self._nzmean(self.attr_matrix(attr_name), 0)

    def to_peaklist(self, ID):
        pl = PeakList(ID, self.mz_mean_vector, self.ints_mean_vector, aligned_ids = self.peaklist_ids)
        pl.add_attribute("present", self.present)
        pl.add_attribute("fraction", self.fraction)
        pl.add_attribute("rsd", self.rsd)
        pl.add_attribute("occurance", self.occurance)
        pl.add_attribute("impure", self.impure)
        return pl

    def remove_samples(self, ids, remove_empty_peaks = True, masked_only = True):
        rmids = np.arange(self.full_shape[0])[self.mask][list(ids)] if masked_only else ids
        self._pids = np.delete(self._pids, rmids, axis = 0)
        self._tags = np.delete(self._tags, rmids, axis = 0)
        self._attr_dict = {k: np.delete(v, rmids, axis = 0) for k, v in self._attr_dict.items()}
        self._mask = np.delete(self._mask, rmids, axis = 0)
        if remove_empty_peaks: self.remove_peaks(np.where(np.sum(self.intensity_matrix, axis = 0) == 0), False)
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def remove_peaks(self, ids, remove_empty_samples = True):
        self._attr_dict = {k: np.delete(v, ids, axis = 1) for k, v in self._attr_dict.items()}
        if remove_empty_samples:
            with unmask_all_peakmatrix(self) as pm: rmsids = np.where(np.sum(pm.intensity_matrix, axis = 1) == 0)
            self.remove_samples(rmsids, False, False)
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def is_empty(self):
        return 0 in self.shape

    # exports
    def to_str(self, attr_name = 'intensity', delimiter = '\t', transpose = False, extend = True):
        hd = ['m/z'] + map(str, self.mz_mean_vector)
        dm = [map(str, self.peaklist_ids)] + [map(str, ln) for ln in self.attr_matrix(attr_name).T]

        if extend:
            ttypes = set(reduce(lambda x,y: x+y, map(lambda x: x.tag_types, self.peaklist_tags)))
            tnum = len(ttypes)
            hd = [hd[0]] + ['missing values'] + map(lambda x: 'tags_' + x, ttypes) + ['tags_untyped'] + hd[1:]
            dm = [dm[0]] + \
                 [map(str, self.missing_values)] + \
                 [map(lambda x: str(x.tag_of(t)) if x.has_tag_type(t) else '', self.peaklist_tags) for t in ttypes] + \
                 [map(lambda x: join(map(str, x.tag_of(None)), ';'), self.peaklist_tags)] + \
                 dm[1:]
            prelst = ['present'] + ([''] * (tnum + 2)) + map(str, self.present)
            rsdlst = ['rsd'] + ([''] * (tnum + 2)) + map(str, self.rsd)
            ocrlst = ['occurance'] + ([''] * (tnum + 2)) + map(str, self.occurance)
            implst = ['impure'] + ([''] * (tnum + 2)) + map(str, self.impure)
            dm = zip(*([prelst, rsdlst, ocrlst, implst] + zip(*dm)))

        lm = [hd] + zip(*dm)
        return join(map(lambda x: join(x, delimiter), zip(*lm) if transpose else lm), '\n')


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
        if len(self._oldmask) != self._pm.full_shape[0]:
            raise ReferenceError('sample number changed within with... statement')
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
        if len(self._oldmask) != self._pm.full_shape[0]:
            raise ReferenceError('sample number changed within with... statement')
        self._pm.mask = self._oldmask

class unmask_all_peakmatrix:
        def __init__(self, pm):
            self._pm = pm
            self._oldmask = pm.mask

        def __enter__(self):
            self._pm.mask  = None
            return self._pm

        def __exit__(self, exc_type, exc_val, exc_tb):
            if len(self._oldmask) != self._pm.full_shape[0]:
                raise ReferenceError('sample number changed within with... statement')
            self._pm.mask = self._oldmask
