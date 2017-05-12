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
def filter_attr(peaks, threshold, snr_attr_name, flag_name, compa = None, on_index = None):
    if not peaks.has_attribute(snr_attr_name):
        raise AttributeError('peaklist object does not have SNR attribute [%s]' % snr_attr_name)
    if compa is None: compa = lambda x, t: x > t
    peaks.add_attribute(flag_name, compa(peaks[snr_attr_name], threshold), is_flag = True, on_index = on_index)
    return peaks

# PeakMatrix filters
def filter_rsd(pm, rsd_threshold = None, qc_label = None):
    if rsd_threshold is None and qc_label is None:
        raise ValueError('must provide rsd threshold or qc label')

    if qc_label is None:
        pm.remove_peaks(np.where([np.isnan(v) or v > rsd_threshold for v in pm.rsd]))
    else:
        if not pm.has_attribute(qc_label):
            raise AttributeError('peaklist object does not have QC label [%s]' % qc_label)
        with mask_peakmatrix(pm, qc_label):
            rsd_values = pm.rsd
        with unmask_peakmatrix(pm, qc_label):
            rmids = np.where([np.isnan(v) or v > rsd_values for v in pm.rsd]) # cannot remove samples inside with... statement
        pm.remove_peaks(rmids)
    return pm

def filter_fraction(pm, fraction_threshold, within_classes = False, class_tag_type = None):
    if not within_classes:
        rmids = np.where(pm.fraction < fraction_threshold)
        pm.remove_peaks(rmids)
    else:
        if not all(map(lambda t: t.has_tag_type(class_tag_type), pm.peaklist_tags)):
            raise AttributeError('not all tags have tag type [%s]' % class_tag_type)
        for tag in pm.tags_of(class_tag_type):
            with mask_peakmatrix(pm, **{class_tag_type: tag}):
                rmids = np.where(pm.fraction < fraction_threshold)
            pm.remove_peaks(rmids)
    return pm

def filter_blank_peaks(pm, blank_label, fraction_threshold = 1.0, fold_threshold = 1.0, method = 'mean', rm_samples = True):
    if blank_label not in pm.peaklist_tag_values:
        raise ValueError('blank label [%s] does not exist' % blank_label)
    if method not in ('mean', 'median', 'max'):
        raise ValueError('filter method must be mean, median or max')

    with mask_peakmatrix(pm, blank_label):
        ints = pm.intensity_matrix if pm.shape[0] == 1 else \
               np.max(pm.intensity_matrix, axis = 0) if method == 'max' else \
               np.array(map(lambda x: getattr(np, method)(x), pm.intensity_matrix.T))
               # note: quick fix of unexpected dtype conversion in apply_along_axis (float64 -> int64)
               # np.apply_along_axis(lambda x: _skipempty(getattr(np, function), x[np.nonzero(x)]), 0, pm.intensity_matrix)
        ints *= fold_threshold

    with unmask_peakmatrix(pm, blank_label):
        rmids = np.where(np.logical_and(ints > 0, pm.fraction < fraction_threshold))
    pm.remove_peaks(rmids)

    if rm_samples:
        pm = pm.remove_samples(np.where(map(lambda x: x.has_tag(blank_label), pm.peaklist_tags)))
    return pm

def filter_sparsity(pm, ppm_threshold):
    mmzs = pm.mzs_mean_vector
    rmids = np.where(np.abs((mmzs[1:] - mmzs[:-1]) / mmzs[1:]) * 1e+6 < ppm_threshold)
    return pm.remove_peaks(rmids)

