#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017-2019 Ralf Weber, Albert Zhou.
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
from functools import reduce
from typing import Union, Sequence, Tuple, Any

import numpy as np
from ..models.peak_matrix import PeakMatrix, mask_peakmatrix, unmask_peakmatrix
from ..models.peaklist import PeakList


# peaklist filters
def filter_attr(pl: PeakList, attr_name: str, max_threshold: Union[int, float, None] = None,
                min_threshold: [int, float, None] = None, flag_name: Union[str, None] = None,
                flag_index: Union[int, None] = None):
    """
    Peaklist attribute values filter.

    :param pl: the target peaklist
    :param attr_name: name of the target attribute
    :param max_threshold: maximum threshold. A peak will be unflagged if the value of it's attr_name is larger than the
        threshold. Default = None, indicating no threshold
    :param min_threshold: Minimum threshold. A peak will be unflagged if the value of it's attr_name is smaller than the
        threshold. Default = None, indicating no threshold
    :param flag_name: name of the new flag attribute. Default = None, indicating using attr_name + '_flag'
    :param flag_index: index of the new flag to be inserted into the peaklist. Default = None
    :rtype: PeakList object

    This filter accepts real value attributes only.

    """
    if min_threshold is None and max_threshold is None:
        raise ValueError('must specify minimum or maximum threshold value')
    flt = lambda x: np.logical_and((min_threshold <= x) if min_threshold is not None else True,
                                   (x <= max_threshold) if max_threshold is not None else True)
    if flag_name is None: flag_name = attr_name + '_flag'
    return pl.add_attribute(flag_name, flt(pl[attr_name]), is_flag=True, on_index=flag_index)


def filter_ringing(pl: PeakList, threshold: float, bin_size: Union[int, float] = 1.0, flag_name: str = 'ringing_flag',
                   flag_index: Union[int, None] = None):
    """
    Peaklist ringing filter.

    :param pl: the target peaklist
    :param threshold: intensity threshold ratio
    :param bin_size: size of the mz chunk for intensity filtering. Default = 1.0 ppm
    :param flag_name: name of the new flag attribute. Default = 'ringing_flag'
    :param flag_index: index of the new flag to be inserted into the peaklist. Default = None
    :rtype: PeakList object

    This filter will split the mz values into bin_size chunks, and search the highest intensity value for each chunk.
    All other peaks, if it's intensity is smaller than threshold x the highest intensity in that chunk, will be unflagged.

    """
    if not 0 <= threshold <= 1:
        raise ValueError('mzr_remove: Provide a value in the range [0.0, 1.0]')
    inds = np.digitize(pl.mz, np.arange(np.floor(np.min(pl.mz)), np.ceil(np.max(pl.mz)) + bin_size, bin_size) - 0.5)
    blks = [(inds == i) for i in np.unique(inds)]
    mask = np.array(reduce(lambda x, y: x + y, [[np.max(pl.intensity[c])] * np.sum(c) for c in blks]))
    return pl.add_attribute(flag_name, pl.intensity > (mask * threshold), is_flag=True, on_index=flag_index)


def filter_mz_ranges(pl: PeakList, mz_ranges: Sequence[Tuple[float, float]], flag_name: str = 'mz_ranges_flag',
                     flagged_only: bool = False, flag_index: Union[int, None] = None):
    """
    Peaklist mz range filter.

    :param pl: the target peaklist
    :param mz_ranges: the mz ranges to remove. Must be in the format of [(mz_min1, mz_max2), (mz_min2, mz_max2), ...]
    :param flag_name: name of the new flag attribute. Default = 'mz_range_remove_flag'
    :param flag_index: index of the new flag to be inserted into the peaklist. Default = None
    :rtype: PeakList

    This filter will remove all the peaks whose mz values are within any of the ranges in the mz_remove_rngs.

    """

    if flagged_only:
        flags = np.ones(pl.shape[0], dtype=bool)
    else:
        flags = np.ones(pl.full_size, dtype=bool)

    for mzr in mz_ranges:
        if len(mzr) != 2:
            raise ValueError(
                'mzr_remove: Provide a list of "start" and "end" values for each m/z range that needs to be removed.')
        if mzr[0] >= mzr[1]:
            raise ValueError('mzr_remove: Start value cannot be larger then end value.')
        flags[
            (pl.get_attribute("mz", flagged_only) >= mzr[0]) & (pl.get_attribute("mz", flagged_only) <= mzr[1])] = False
    pl.add_attribute(flag_name, flags, flagged_only=flagged_only, is_flag=True, on_index=flag_index)
    return pl


# PeakMatrix filters
def filter_rsd(pm: PeakMatrix, rsd_threshold: Union[int, float], qc_tag: Any, on_attr: str = 'intensity',
               flag_name: str = 'rsd_flag'):
    """
    PeakMatrix RSD filter.

    :param pm: the target peak matrix
    :param rsd_threshold: threshold of the RSD of the QC samples
    :param qc_tag: tag (label) to unmask qc samples
    :param on_attr: calculate RSD on given attribute. Default = "intensity"
    :param flag_name: name of the new flag. Default = 'rsd_flag'
    :rtype: PeakMatrix

    This filter will calculate the RSD values of the QC samples. A peak with a QC RSD value larger than the
    threshold will be unflagged.

    """
    rsd_values = pm.rsd(qc_tag, on_attr=on_attr)
    if np.any(np.isnan(rsd_values)):
        logging.warning('nan found in QC rsd values, filter might not work properly')

    pm.add_flag(flag_name, [not (np.isnan(v) or v > rsd_threshold) for v in rsd_values])
    return pm


def filter_fraction(pm: PeakMatrix, fraction_threshold: float, within_classes: bool = False, class_tag_type: Any = None,
                    flag_name: str = 'fraction_flag'):
    """
    PeakMatrix fraction filter.

    :param pm: the target peak matrix
    :param fraction_threshold: threshold of the sample fractions
    :param within_classes: whether to calculate the fraction array within each class. Default = False
    :param class_tag_type: tag type to unmask samples within the same class (e.g. "classLabel"). Default = None
    :param flag_name: name of the new flag. Default = 'fraction_flag'
    :rtype: PeakMatrix object

    This filter will calculate the fraction array over all samples or within each class (based on class_tag_type).
    The peaks with a fraction value smaller than the threshold will be unflagged.

    """
    if not within_classes:
        pm.add_flag(flag_name, pm.fraction >= fraction_threshold)
    else:
        if class_tag_type is None:
            raise KeyError('must provide class tag type for within classes filtering')
        if not all([t.has_tag_type(class_tag_type) for t in pm.peaklist_tags]):
            raise AttributeError('not all tags have tag type [%s]' % class_tag_type)
        flg = np.zeros(pm.shape[1])
        for tag in pm.tags_of(class_tag_type):
            with unmask_peakmatrix(pm, tag) as m:
                flg = np.logical_or(flg, (m.fraction >= fraction_threshold))
        pm.add_flag(flag_name, flg)
    return pm


def filter_blank_peaks(pm: PeakMatrix, blank_tag: Any, fraction_threshold: Union[int, float] = 1,
                       fold_threshold: Union[int, float] = 1,
                       method: str = 'mean', rm_blanks: bool = True, flag_name: str = 'blank_flag'):
    """
    PeakMatrix blank filter.

    :param pm: the target peak matrix
    :param blank_tag: tag (label) to mask blank samples. e.g Tag("blank", "classLabel")
    :param fraction_threshold: threshold of the sample fractions. Default = 1
    :param fold_threshold: threshold of the blank sample intensity folds. Default = 1
    :param method: method to calculate blank sample intensity array. Valid values include 'mean', 'median', and 'max'.
        Default = 'mean'
    :param rm_blanks: whether to remove (not mask) blank samples after filtering
    :param flag_name: name of the new flag. Default = 'blank_flag'
    :rtype: PeakMatrix object

    This filter will calculate the intensity array of the blanks using the "method", and compare with the
    intensities of the other samples. If fraction_threshold% of the intensity values of a peak are smaller than the
    blank intensities x fold_threshold, this peak will be unflagged.

    """
    if not any([blank_tag in x for x in pm.peaklist_tags]):
        raise ValueError('blank tag [%s] does not exist' % blank_tag)
    if method not in ('mean', 'median', 'max'):
        raise ValueError('filter method must be mean, median or max')

    with unmask_peakmatrix(pm, blank_tag) as m:
        mm = np.ma.masked_array(m.intensity_matrix, mask = ~(m.intensity_matrix > 0))
        ints = mm[0] if mm.shape[0] == 1 else getattr(np, method)(mm, axis = 0)
        imsk = ints.mask
        ints = np.array(ints) * fold_threshold

    with mask_peakmatrix(pm, blank_tag) as m:
        faild_int = np.sum(m.intensity_matrix >= ints, axis=0) < (fraction_threshold * m.shape[0])
        m.add_flag(flag_name, ~(~imsk & faild_int))

    if rm_blanks:
        pm = pm.remove_samples(np.where([x.has_tag(blank_tag) for x in pm.peaklist_tags])[0])
    return pm

