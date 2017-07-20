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
from collections import OrderedDict, Iterable
from string import join
from peaklist import PeakList


class PeakMatrix(object):
    def __init__(self, peaklist_ids, peaklist_tags, peaklist_attributes):
        if len(set(peaklist_ids)) != len(peaklist_ids):
            raise ValueError('duplicate peaklist IDs found')
        self._pids = np.array(peaklist_ids)
        self._tags = np.array(peaklist_tags)

        self._attr_dict = OrderedDict((str(k), np.array(v)) for k,v in peaklist_attributes)
        if not (self._attr_dict.keys()[0] == 'mz' and self._attr_dict.keys()[1] == 'intensity'):
            raise ValueError('required attribute fields "mz" and "intensity" not available')

        pnum, psize = self._attr_dict.values()[0].shape
        if not len(peaklist_ids) == len(peaklist_tags) == pnum:
            raise ValueError('input peaklists number not match')
        if not all(map(lambda x: x.shape == (pnum, psize), self._attr_dict.values())):
            raise ValueError('input attribute matrix shape not match')

        self._mask = np.zeros(pnum, dtype=bool)
        self._flags_dict = OrderedDict()

    # privates
    # np.average accepts weights, set to 0 to exclude value
    def _present_mask(self, x, _func, flagged_only):
        pmx = self.present_matrix if flagged_only else self.full_present_matrix
        mskx = np.ma.masked_array(x, mask = ~pmx)
        vals = _func(mskx)
        rval = np.array(vals)
        rval[vals.mask] = np.nan
        return rval

    def _present_mean(self, x, axis, flagged_only):
        return self._present_mask(x, lambda v: np.mean(v, axis = axis), flagged_only)

    def _present_std(self, x, axis, flagged_only):
        return self._present_mask(x, lambda v: np.std(v, axis = axis), flagged_only)

    # built-ins
    def __len__(self):
        return self.shape[0]

    def __str__(self):
        return self.to_str(attr_name='intensity', delimiter=', ', transpose=True, comprehensive=False)

    # properties
    @property
    def mask(self):  # mask is different from flags in that they will not be masked by themselves, i.e., more temporary
        return self._mask.copy()

    @mask.setter
    def mask(self, value):
        self._mask = np.zeros(self.full_shape[0], dtype = bool)
        if value is not None: self._mask[value] = True

    @property
    def flag_names(self):
        return tuple(self._flags_dict.keys())

    @property
    def flags(self):
        """
        Property of the flags.

        :getter: returns a deep copy of the flags array
        :type: numpy array

        """
        return np.ones(self.full_shape[1], dtype = bool) if len(self._flags_dict) == 0 else \
               np.logical_and.reduce(self._flags_dict.values())

    @property
    def attributes(self):
        return tuple(self._attr_dict.keys())

    @property
    def peaklist_ids(self):
        return tuple(self._pids[~self._mask])

    @property
    def peaklist_tags(self):
        return tuple(self._tags[~self._mask])

    @property
    def peaklist_tag_types(self):
        return tuple(set(reduce(lambda x, y: x + y, [t.tag_types for t in self.peaklist_tags], ())))

    @property
    def peaklist_tag_values(self):
        return tuple(set(reduce(lambda x, y: x + y, [t.tag_values for t in self.peaklist_tags], ())))

    @property
    def shape(self):
        return np.sum(~self._mask), np.sum(self.flags)

    @property
    def full_shape(self):
        return self._attr_dict['mz'].shape

    @property
    def present(self):
        return np.sum(self.present_matrix, axis=0)

    @property
    def present_matrix(self):
        pmx = self.attr_matrix('intra_count') if 'intra_count' in self.attributes else self.mz_matrix
        return pmx > 0

    @property
    def full_present_matrix(self):
        pmx = self.attr_matrix('intra_count', flagged_only = False) if 'intra_count' in self.attributes else \
              self.attr_matrix('mz', flagged_only = False)
        return pmx > 0

    @property
    def fraction(self):
        return self.present / self.shape[0]

    @property
    def missing_values(self):
        return np.sum(~self.present_matrix, axis=1)

    @property
    def rsd(self):
        if self.shape[0] < 2:
            # logging.warning('calculating RSD on less than 2 samples')
            return np.ones(self.shape[1]) * np.nan
        ints = self.intensity_matrix
        rsd = self._present_std(ints, 0, True) / self._present_mean(ints, 0, True) * 100
        rsd[np.where(map(lambda x: len(set(x[np.nonzero(x)])) == 1, ints.T))] = np.nan # only one valid value
        return rsd

    @property
    def occurrence(self):
        if not self._attr_dict.has_key('intra_count'):
            logging.warning("attribute matrix ['intra_count'] not available")
            return np.ones(self.shape[1])
        return np.sum(self.attr_matrix('intra_count'), axis=0)

    @property
    def purity(self):
        if not self._attr_dict.has_key('intra_count'):
            logging.warning("attribute matrix ['intra_count'] not available")
            return np.zeros(self.shape[1])
        return 1 - (np.sum(self.attr_matrix('intra_count') > 1, axis=0) / self.shape[0])

    # attribute matrix
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
    def intensity_mean_vector(self):
        return self.attr_mean_vector('intensity')

    # publics
    # tags and mask
    def tags_of(self, tag_type=None):
        if not (tag_type is None or all(map(lambda x: x.has_tag_type(tag_type), self.peaklist_tags))):
            raise ValueError('not all samples has tag type [%s]' % tag_type)
        tlst = [t.tag_of(tag_type) for t in self.peaklist_tags]
        if tag_type is None: tlst = reduce(lambda x, y: x + y, tlst)
        return tuple(set(tlst))

    def mask_tags(self, *args, **kwargs):  # match to all
        override = kwargs.pop('override') if kwargs.has_key('override') else False
        mask = map(lambda x: all(map(lambda t: x.has_tag(t), args)) and
                             all(map(lambda t: x.has_tag(**dict([t])), kwargs.items())), self._tags)
        self.mask = np.logical_or(False if override else self._mask, mask)
        return self

    def unmask_tags(self, *args, **kwargs):  # match to all
        override = kwargs.pop('override') if kwargs.has_key('override') else False
        mask = map(lambda x: not (all(map(lambda t: x.has_tag(t), args)) and
                                  all(map(lambda t: x.has_tag(**dict([t])), kwargs.items()))), self._tags)
        self.mask = np.logical_and(True if override else self._mask, mask)
        return self

    # flags
    def add_flag(self, flag_name, flag_values, flagged_only = True):
        if flag_name == 'flags':
            raise KeyError('reserved flag name [flags] cannot be added')
        if self._flags_dict.has_key(flag_name):
            raise KeyError('flag name [%s] already exists' % flag_name)
        if flagged_only:
            if not len(flag_values) == self.shape[1]: raise ValueError('flag values and peak matrix shape not match')
            newf = np.zeros(self.full_shape[1], dtype = bool)
            newf[self.flags] = flag_values
            self._flags_dict[flag_name] = newf
        else:
            if not len(flag_values) == self.full_shape[1]: raise ValueError('flag values and peak matrix shape not match')
            self._flags_dict[flag_name] = np.array(flag_values, dtype = bool)

    def drop_flag(self, flag_name):
        if flag_name == 'flags':
            raise KeyError('reserved flag name [flags] cannot be droppped')
        if not self._flags_dict.has_key(flag_name):
            raise KeyError('flag name [%s] not exists' % flag_name)
        del self._flags_dict[flag_name]

    def flag_values(self, flag_name):
        if flag_name == 'flags':
            raise KeyError('use PeakMatrix.flags to access reserved flag name [flags]')
        if not self._flags_dict.has_key(flag_name):
            raise KeyError('flag name [%s] not exists' % flag_name)
        return self._flags_dict[flag_name]

    # access
    def attr_matrix(self, attr_name, flagged_only = True):
        if not self._attr_dict.has_key(attr_name):
            raise KeyError('attribute matrix [%s] not available' % attr_name)
        aM = self._attr_dict[attr_name][~self._mask]
        return aM[:,self.flags] if flagged_only else aM

    def attr_mean_vector(self, attr_name, flagged_only = True):
        aM = self.attr_matrix(attr_name, flagged_only)
        aV = self._present_mean(aM, 0, flagged_only) if aM.dtype.kind in ('i', 'u', 'f') else \
             np.array([join([str(v) for v,p in zip(ln,self.present_matrix[:,j]) if p],';')
                       for j,ln in enumerate(zip(*aM))]) # for strings
        return aV

    def remove_samples(self, sample_ids, masked_only=True):
        if isinstance(sample_ids, Iterable): sample_ids = list(sample_ids)
        rmids = np.where(~self.mask)[0][sample_ids] if masked_only else sample_ids
        self._pids = np.delete(self._pids, rmids, axis=0)
        self._tags = np.delete(self._tags, rmids, axis=0)
        self._attr_dict = OrderedDict((k, np.delete(v, rmids, axis=0)) for k,v in self._attr_dict.items())
        self._mask = np.delete(self._mask, rmids, axis=0)
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def remove_peaks(self, peak_ids, flagged_only = True):
        if isinstance(peak_ids, Iterable): peak_ids = list(peak_ids)
        rmids = np.where(self.flags)[0][peak_ids] if flagged_only else peak_ids
        self._attr_dict = OrderedDict((k, np.delete(v, rmids, axis=1)) for k,v in self._attr_dict.items())
        self._flags_dict = OrderedDict((k, np.delete(v, rmids)) for k,v in self._flags_dict.items())
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def is_empty(self):
        return 0 in self.shape

    # exports
    def extract_peaklist(self, peaklist_id):
        if peaklist_id not in self.peaklist_ids:
            raise ValueError('peaklist id has to match those in the peak matrix')
        idx = self.peaklist_ids.index(peaklist_id)

        nzero_idx = self.present_matrix[idx]
        pl = PeakList(peaklist_id, self.mz_matrix[idx,nzero_idx], self.intensity_matrix[idx,nzero_idx])
        for attr in self.attributes[2:]: pl.add_attribute(attr, self.attr_matrix(attr)[idx,nzero_idx])

        return pl

    def extract_peaklists(self):
        return map(self.extract_peaklist, self.peaklist_ids)

    def to_peaklist(self, ID):
        presids = self.present > 0 # presented peaks only
        if False in presids:
            logging.warning('[%d] empty peaks removed when exporting PeakMatrix to PeakList' % np.sum(~presids))
        pl = PeakList(ID, self.mz_mean_vector[presids], self.intensity_mean_vector[presids], aligned_ids = self.peaklist_ids)
        pl.add_attribute('present', self.present[presids])
        pl.add_attribute('fraction', self.fraction[presids])
        pl.add_attribute('rsd', self.rsd[presids])
        pl.add_attribute('occurrence', self.occurrence[presids])
        pl.add_attribute('purity', self.purity[presids])
        return pl

    def to_str(self, attr_name='intensity', delimiter='\t', transpose=False, comprehensive=True):
        hd = ['m/z'] + map(str, self.attr_mean_vector('mz', flagged_only = False))
        dm = [map(str, self.peaklist_ids)] + \
             [map(str, ln) for ln in self.attr_matrix(attr_name, flagged_only = False).T]

        if comprehensive:
            ttypes = set(reduce(lambda x, y: x + y, map(lambda x: x.tag_types, self.peaklist_tags)))
            tnum = len(ttypes)
            hd = [hd[0]] + ['missing values'] + map(lambda x: 'tags_' + x, ttypes) + ['tags_untyped'] + hd[1:]
            dm = [dm[0]] + \
                 [map(str, self.missing_values)] + \
                 [map(lambda x: str(x.tag_of(t)) if x.has_tag_type(t) else '', self.peaklist_tags) for t in ttypes] + \
                 [map(lambda x: join(map(str, x.tag_of(None)), ';'), self.peaklist_tags)] + \
                 dm[1:]

            def _refill(vals):
                vect = np.array(['--'] * self.full_shape[1], dtype = vals.dtype)
                vect[np.where(self.flags)] = vals
                return list(vect)

            prelst = ['present']    + ([''] * (tnum + 2)) + _refill(self.present.astype(str))
            rsdlst = ['rsd']        + ([''] * (tnum + 2)) + _refill(self.rsd.astype(str))
            ocrlst = ['occurrence'] + ([''] * (tnum + 2)) + _refill(self.occurrence.astype(str))
            puplst = ['purity']     + ([''] * (tnum + 2)) + _refill(self.purity.astype(str))
            flgmtx = [[fn] + ([''] * (tnum + 2)) + map(str, self.flag_values(fn)) for fn in self.flag_names]
            flglst = ['flags']      + ([''] * (tnum + 2)) + map(str, self.flags)
            dm = zip(*([prelst, rsdlst, ocrlst, puplst] + flgmtx + [flglst] + zip(*dm)))

        lm = [hd] + zip(*dm)
        return join(map(lambda x: join(x, delimiter), zip(*lm) if transpose else lm), '\n')


# with statements
class mask_peakmatrix:
    def __init__(self, pm, *args, **kwargs):
        self._pm = pm
        self._utags = args
        self._ttags = kwargs
        if not self._ttags.has_key('override'): self._ttags['override'] = True  # default for with statement
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.mask_tags(*self._utags, **self._ttags)
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]


class unmask_peakmatrix:
    def __init__(self, pm, *args, **kwargs):
        self._pm = pm
        self._utags = args
        self._ttags = kwargs
        if not self._ttags.has_key('override'): self._ttags['override'] = True  # default for with statement
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.unmask_tags(*self._utags, **self._ttags)
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]


class mask_all_peakmatrix:
    def __init__(self, pm):
        self._pm = pm
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.mask = [True] * self._pm.full_shape[0]
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]


class unmask_all_peakmatrix:
    def __init__(self, pm):
        self._pm = pm
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.mask = None
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]

