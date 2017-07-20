#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
filters: PeakList and PeakMatrix filters

author(s): Albert Zhou, Ralf Weber
origin: Nov. 2, 2016

"""

from __future__ import division

import logging
import numpy as np
from dimspy.models.peak_matrix import mask_peakmatrix, unmask_peakmatrix


# peaklist filters
def filter_attr(pl, attr_name, max_threshold=None, min_threshold=None, flag_name=None, flag_index=None):
    if min_threshold is None and max_threshold is None:
        raise ValueError('must specify minimum or maximum threshold value')
    flt = lambda x: np.logical_and((min_threshold <= x) if min_threshold is not None else True,
                                   (x <= max_threshold) if max_threshold is not None else True)
    if flag_name is None: flag_name = attr_name + '_flag'
    return pl.add_attribute(flag_name, flt(pl[attr_name]), is_flag=True, on_index=flag_index)


def filter_ringing(pl, threshold, bin_size=1.0, flag_name='ringing_flag', flag_index=None):
    if not 0 <= threshold <= 1:
        raise ValueError('mzr_remove: Provide a value in the range [0.0, 1.0]')
    inds = np.digitize(pl.mz, np.arange(np.floor(np.min(pl.mz)), np.ceil(np.max(pl.mz)) + bin_size, bin_size) - 0.5)
    blks = [(inds == i) for i in np.unique(inds)]
    mask = np.array(reduce(lambda x, y: x + y, [[np.max(pl.intensity[c])] * np.sum(c) for c in blks]))
    return pl.add_attribute(flag_name, pl.intensity > (mask * threshold), is_flag=True, on_index=flag_index)


def filter_mz_ranges(pl, mzrs, flag_name='mzrs_remove_flag', flag_index=None):
    mzrs_removed_flags = np.ones(pl.shape[0], dtype=bool)
    for mzr in mzrs:
        if len(mzr) != 2:
            raise ValueError('mzr_remove: Provide a list of "start" and "end" values for each m/z range that needs to be removed.')
        if mzr[0] >= mzr[1]:
            raise ValueError('mzr_remove: Start value cannot be larger then end value.')
        mzrs_removed_flags[(pl.mz >= mzr[0]) & (pl.mz <= mzr[1])] = False
    pl.add_attribute(flag_name, mzrs_removed_flags, is_flag=True, on_index=flag_index)
    return pl


# PeakMatrix filters
def filter_rsd(pm, rsd_threshold, qc_label='qc', flag_name='rsd_flag'):
    if qc_label not in pm.peaklist_tag_values:
        raise AttributeError('peaklist object does not have QC label [%s]' % str(qc_label))

    with unmask_peakmatrix(pm, qc_label) as m: rsd_values = m.rsd
    if np.any(np.isnan(rsd_values)):
        logging.warning('nan found in QC rsd values, filter might not work properly')

    pm.add_flag(flag_name, [not (np.isnan(v) or v > rsd_threshold) for v in rsd_values])
    return pm


def filter_fraction(pm, fraction_threshold, within_classes=False, class_tag_type=None, flag_name='fraction_flag'):
    if not within_classes:
        pm.add_flag(flag_name, pm.fraction >= fraction_threshold)
    else:
        if not all(map(lambda t: t.has_tag_type(class_tag_type), pm.peaklist_tags)):
            raise AttributeError('not all tags have tag type [%s]' % class_tag_type)
        flg = np.ones(pm.shape[1])
        for tag in pm.tags_of(class_tag_type):
            with unmask_peakmatrix(pm, **{class_tag_type: tag}) as m:
                flg = np.logical_and(flg, (m.fraction >= fraction_threshold))
        pm.add_flag(flag_name, flg)
    return pm


def filter_blank_peaks(pm, blank_label, fraction_threshold=1, fold_threshold=1, method='mean', rm_blanks=True, flag_name='blank_flag'):
    if blank_label not in pm.peaklist_tag_values:
        raise ValueError('blank label [%s] does not exist' % blank_label)
    if method not in ('mean', 'median', 'max'):
        raise ValueError('filter method must be mean, median or max')

    with unmask_peakmatrix(pm, blank_label) as m:
        ints = m.intensity_matrix if m.shape[0] == 1 else \
               np.max(m.intensity_matrix, axis=0) if method == 'max' else \
               np.array(map(lambda x: getattr(np, method)(x), m.intensity_matrix.T))
        ints *= fold_threshold

    with mask_peakmatrix(pm, blank_label) as m:
        faild_int = np.sum(m.intensity_matrix >= ints, axis=0) < (fraction_threshold * m.shape[0])
        m.add_flag(flag_name, ~((ints > 0) & faild_int))

    if rm_blanks:
        pm = pm.remove_samples(np.where(map(lambda x: x.has_tag(blank_label), pm.peaklist_tags)))
    return pm

