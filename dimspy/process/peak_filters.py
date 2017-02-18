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
from ..models.peak_matrix import mask_peakmatrix, unmask_peakmatrix


# peaklist filters
def filter_abs_intensity(peaks, threshold, flag_name='cutoff_ints_filter'):
    peaks.add_attribute(flag_name, peaks.ints > threshold, is_flag=True)
    return peaks


def filter_rel_intensity(peaks, threshold, flag_name='relative_ints_filter'):
    return filter_abs_intensity(peaks, np.max(peaks.ints) * threshold, flag_name)


def filter_snr(peaks, threshold, snr_attr_name='snr', flag_name='snr_flag'):
    assert peaks.has_attribute(snr_attr_name), 'peak list does not have SNR attribute [%s]' % snr_attr_name
    peaks.add_attribute(flag_name, peaks[snr_attr_name] > threshold, is_flag=True)
    return peaks


# PeakMatrix filters
def filter_rsd(pm, threshold, qc_label):
    if qc_label is not None:
        assert qc_label in pm.peaklist_tag_values, "QC label does not exist"
        with mask_peakmatrix(pm, qc_label):
            rsd_values = pm.rsd
        with unmask_peakmatrix(pm, qc_label):
            pm.remove_peaks(np.where(pm.rsd > rsd_values))  # TODO: compare when rsd values are nan?
    else:
        pm.remove_peaks(np.where(pm.rsd > threshold))
    return pm


def filter_across_classes(pm, min_fraction):
    rmids = np.where(np.sum(pm.intensity_matrix > 0, axis=0) / pm.shape[0] < min_fraction)
    return pm.remove_peaks(rmids)


def filter_within_classes(pm, tag_type, min_fraction):
    for tag in pm.tags_of(tag_type):
        with mask_peakmatrix(pm, **{tag_type:tag}):
            pm.remove_peaks(np.where(np.sum(pm.intensity_matrix > 0, axis=0) / pm.shape[0] < min_fraction))
    return pm


def filter_blank_peaks(pm, blank_label, min_fraction=1.0, min_fold=1.0, function="mean", rm_samples=True):
    print blank_label, type(blank_label)
    print pm.peaklist_tag_values
    for t in pm.peaklist_tags: print t

    assert blank_label in pm.peaklist_tag_values, "Blank label does not exist"
    assert 0 < min_fraction <= 1, "Provide a value between 0. and 1."
    assert min_fold >= 0, "Provide a value larger than zero."
    assert function in ("mean", "median", "max"), "Mean, median or max intensity"

    with mask_peakmatrix(pm, blank_label):
        ints = pm.intensity_matrix if pm.shape[0] == 1 else \
               np.max(pm.intensity_matrix, axis=0) if function == "max" else \
               np.array(map(lambda x: getattr(np, function)(x), pm.intensity_matrix.T))
               # quick fix of unexpected dtype conversion in apply_along_axis (float64 -> int64)
               # np.apply_along_axis(lambda x: _skipempty(getattr(np, function), x[np.nonzero(x)]), 0, pm.intensity_matrix)
        ints *= min_fold

    with unmask_peakmatrix(pm, blank_label):
        has_blank = ints > 0
        faild_int = np.sum(pm.intensity_matrix >= ints, axis = 0) < (min_fraction * pm.shape[0])
        pm.remove_peaks(np.where(np.logical_and(has_blank, faild_int)))

    # min_fraction of the non-zero values
    # pm.intensity_matrix.shape[0] should non-zero value only instead
    # column is zero

    if rm_samples:
        pm = pm.remove_samples(np.where(map(lambda x: x.has_tags(blank_label), pm.peaklist_tags)))
    return pm


def filter_sparsity(pm, ppm):
    mmzs = pm.mzs_mean_vector
    rmids = np.where(np.abs((mmzs[1:] - mmzs[:-1]) / mmzs[1:]) * 1e+6 < ppm)
    return pm.remove_peaks(rmids)

