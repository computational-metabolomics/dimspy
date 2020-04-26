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


import collections
import csv
import os
import re
import warnings
from typing import Sequence, Dict

import numpy as np

from .models.peak_matrix import PeakMatrix
from .models.peaklist import PeakList


def mz_range_from_header(h: str) -> Sequence[float]:
    """
    Extract m/z range from header or filter string

    :param h: Header or filter string
    :return: m/z range
    """

    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


def ms_type_from_header(h: str) -> str:
    """
    Extract the ms type from header or filter string

    :param h: header or filter string
    :return: ms type (e.g. FTMS and ITMS)
    """

    return h.split(" ")[0]


def scan_type_from_header(h: str) -> str:
    """
    Extract the scan type from the header of filter string

    :param h: header or filter string
    :return: Scan type (e.g. full or sim)
    """

    if " full " in h.lower():
        return "Full"
    elif " sim " in h.lower():
        return "SIM"
    else:
        return None


def mode_type_from_header(h: str) -> str:
    """
    Extract scan mode from the header of filter string

    :param h: header or filter string
    :return: Scan type (e.g. p = profile, c = centroid)
    """

    if " p " in h.lower():
        return "p"
    elif " c " in h.lower():
        return "c"
    else:
        return None


def count_scan_types(hs: list) -> int:
    """
    Count the number of unique scan types

    :param hs: List of headers or filter strings
    :return: Count
    """

    return len(set([scan_type_from_header(h) for h in hs]))


def count_ms_types(hs: list) -> int:
    """
    Count the number of unique ms types

    :param hs: List of headers or filter strings
    :return: Count
    """

    return len(set([ms_type_from_header(h) for h in hs]))


def _partially_overlapping_windows(mzrs: list) -> list:
    """
    Select all adjacent m/z windows that partially overlap
    For example: [100-200] and [185-285] (Valid for SIM-stitch)

    :param mzrs: Nested list of mz ranges / windows
    :return: Nested list of m/z ranges / windows
    """

    assert type(mzrs) == list, "List required"
    temp = []
    for i in range(0, len(mzrs) - 1):
        if mzrs[i][0] < mzrs[i + 1][0] and mzrs[i][1] > mzrs[i + 1][0] and mzrs[i][1] < mzrs[i + 1][1]:
            if mzrs[i] not in temp:
                temp.append(mzrs[i])
            if mzrs[i + 1] not in temp:
                temp.append(mzrs[i + 1])
    return temp


def _first_fully_overlapping_windows(mzrs: list) -> list:
    """
    Select m/z windows that fall within another window and have a different mass ranges
    For example: [100-200] and [125-175] (Invalid)

    :param mzrs: Nested list of m/z ranges / windows
    :return: Nested list of m/z ranges / windows
    """

    assert type(mzrs) == list, "List required"

    for i in range(0, len(mzrs) - 1):
        if mzrs[i][0] <= mzrs[i + 1][0] and mzrs[i][1] >= mzrs[i + 1][1]:
            return mzrs[i], mzrs[i + 1]  # Temporary print
    return []


def _non_overlapping_windows(mzrs: list) -> list:
    """
    Select windows that do not overlap with other windows.
    For example: [100-200] and [200-400] (Valid for merging)

    :param mzrs: Nested list of m/z ranges / windows
    :return: Nested list of m/z ranges / windows
    """

    assert type(mzrs) == list, "List required"
    temp = []

    for i in range(0, len(mzrs)):
        c = 0
        for j in range(0, len(mzrs)):
            if mzrs[i][0] <= mzrs[j][0] and mzrs[i][1] <= mzrs[j][0]:
                c += 1
            elif mzrs[i][0] >= mzrs[j][1] and mzrs[i][1] >= mzrs[j][1]:
                c += 1
        if c == len(mzrs) - 1:
            temp.append(mzrs[i])
    return temp


def interpret_method(mzrs: list):
    """
    Interpret and define type of method

    :param mzrs: Nested list of m/z ranges / windows
    :return: Type of MS method
    """

    mzrs.sort(key=lambda x: x[1])

    now = _non_overlapping_windows(mzrs)
    pow = _partially_overlapping_windows(mzrs)

    if len(mzrs) == 1:
        print("Single m/z window.....")
        method = "single"
    elif len(now) == len(mzrs):
        print("Adjacent m/z windows.....")
        method = "adjacent"
    elif len(pow) == len(mzrs):
        print("SIM-Stitch method - Overlapping m/z windows.....")
        method = "overlapping"
    else:
        raise IOError("SIM-Stitch cannot be applied; 'filter_scan_events' required or set 'skip_stitching' to False")

    return method

def to_int(x):
    """
    :param x: Value to convert to int
    :return: Value as int (or False if conversion not possible)
    """
    try:
        i = int(x)
        return i
    except ValueError as e:
        return False


def validate_metadata(fn_tsv: str) -> collections.OrderedDict:
    """
    Check and validate metadata within a tab-separated file

    :param fn_tsv: Path to tab-separated file
    :return: Dictionary
    """

    assert os.path.isfile(fn_tsv.encode('unicode_escape')), "{} does not exist".format(fn_tsv)
    with open(fn_tsv.encode('unicode_escape')) as tsv:
        fm_dict = collections.OrderedDict()
        for row in csv.DictReader(tsv, delimiter="\t"):
            for k, v in row.items():
                fm_dict.setdefault(k, []).append(v)

    if "filename" not in fm_dict:
        raise IOError("Column 'filename' missing.")

    unique, counts = np.unique(fm_dict["filename"], return_counts=True)
    if len(unique) != sum(counts):
        raise ValueError("Duplicate filename in list")

    # convert relevant columns to int
    for h in ['replicate', 'batch', 'injectionOrder', 'multilist']:
        if h in fm_dict:
            int_l = []
            for c, x in enumerate(fm_dict[h]):
                i = to_int(x)
                assert to_int(i), "Column '{}' values should be integers, see row {}".format(h, c+1)
                int_l.append(i)
            fm_dict[h] = int_l

    if "replicate" in fm_dict.keys():

        if 0 in fm_dict["replicate"]:
            raise IOError("Incorrect replicate number in list. Row {}".format(list(fm_dict["replicate"]).index(0)))

        idxs_replicates = idxs_reps_from_filelist(fm_dict["replicate"])
        counts = {}
        for idxs in idxs_replicates:
            if len(idxs) not in counts:
                counts[len(idxs)] = 1
            else:
                counts[len(idxs)] += 1
        for k, v in list(counts.items()):
            print("{} sample(s) with {} replicate(s)".format(v, k))
    else:
        print("Column for replicate numbers missing. Only required for replicate filter.")

    if "batch" in fm_dict.keys():
        unique_batches, counts = np.unique(fm_dict["batch"], return_counts=True)
        print("Batch numbers:", unique_batches)
        print("Number of samples in each Batch:", dict(list(zip(unique_batches, counts))))
    else:
        print("Column for batch number missing. Not required.")

    if "injectionOrder" in fm_dict:
        assert np.array_equal(fm_dict["injectionOrder"], sorted(
            fm_dict["injectionOrder"])), "Check the injectionOrder column - samples not in order"
    else:
        print("Column for sample injection order missing. Not required.")

    if "classLabel" in fm_dict:
        if "replicate" in fm_dict:
            for i in range(len(idxs_replicates)):
                assert len(np.unique(fm_dict["classLabel"][min(idxs_replicates[i]):max(
                    idxs_replicates[i]) + 1])) == 1, "class names do not match with number of replicates"
        unique, counts = np.unique(fm_dict["classLabel"], return_counts=True)
        cls = dict(list(zip(unique, counts)))
        print("Classes:", cls)
    else:
        warnings.warn("Column 'classLabel' for class labels missing. Not required.")

    if "multilist" not in fm_dict:
        print("Column 'multilist' for spliting peaklists is missing. Not required.")

    return fm_dict


def update_metadata_and_labels(peaklists: Sequence[PeakList], fl: Dict):
    """
    Update metadata

    :param peaklists: List of peaklist Objects
    :param fl: Dictionary with meta data
    :return: List of peaklist objects
    """

    if not isinstance(peaklists[0], PeakList):
        raise IOError("PeakList object required")

    for k in list(fl.keys()):
        for pl in peaklists:
            if pl.ID not in fl[list(fl.keys())[0]]:
                raise IOError("filelist and peaklist do not match {}".format(pl.ID))

            index = fl[list(fl.keys())[0]].index(pl.ID)
            pl.metadata[k] = fl[k][index]
            # pl.metadata["filelist"] = {k:fl[k][index] for k in fl.keys()}

            for tag_name in ["replicate", "replicates", "batch", "injectionOrder", "classLabel"]:
                if tag_name in list(fl.keys()):
                    if pl.tags.has_tag_type(tag_name):
                        pl.tags.drop_tag_type(tag_name)
                    pl.tags.add_tag(fl[tag_name][index], tag_name)

    return peaklists


def update_labels(pm: PeakMatrix, fn_tsv: str) -> PeakMatrix:
    """
    Update Sample labels PeakMatrix object
    :param pm: peakMatrix Object
    :param fn_tsv: Path to tab-separated file
    :return: peakMatrix Object
    """

    assert os.path.isfile(fn_tsv.encode('unicode_escape')), "{} does not exist".format(fn_tsv)

    with open(fn_tsv.encode('unicode_escape')) as tsv:
        fm_dict = collections.OrderedDict()
        for row in csv.DictReader(tsv, delimiter="\t"):
            for k, v in row.items():
                fm_dict.setdefault(k, []).append(v)

    assert "sample_id" == list(fm_dict.keys())[0] or "filename" == list(fm_dict.keys())[
        0], "Column for class labels not available"
    assert "classLabel" in fm_dict.keys(), "Column for class label (classLabel) not available"
    assert (fm_dict[list(fm_dict.keys())[0]] == pm.peaklist_ids).all(), "Sample ids do not match {}".format(
        np.setdiff1d(fm_dict[list(fm_dict.keys())[0]], pm.peaklist_ids))

    for tag_name in ["replicate", "replicates", "batch", "injectionOrder", "classLabel"]:
        if tag_name in fm_dict:
            for i in range(len(fm_dict[tag_name])):
                if pm.peaklist_tags[i].has_tag_type(tag_name):
                    pm.peaklist_tags[i].drop_tag_type(tag_name)
                pm.peaklist_tags[i].add_tag(fm_dict[tag_name][i], tag_name)
    return pm


def idxs_reps_from_filelist(replicates: list):
    """

    :param replicates:
    :return:
    """

    idxs, temp = [], [0]
    replicates = [int(r) for r in replicates]
    for i in range(1, len(replicates)):
        if (replicates[i - 1] == replicates[i] or replicates[i - 1] > replicates[i]) and replicates[i] == 1:
            idxs.append(temp)
            temp = [i]
        elif replicates[i - 1] < replicates[i] and replicates[i - 1] - replicates[i] == -1:
            temp.append(i)
        else:
            raise ValueError("Incorrect numbering for replicates. Row {}".format(i))
    idxs.append(temp)
    return idxs
