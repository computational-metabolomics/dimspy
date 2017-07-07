#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import warnings
import collections
import re
import numpy as np
from models.peaklist import PeakList


def mz_range_from_header(h):
    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


def ms_type_from_header(h):
    return h.split(" ")[0]


def scan_type_from_header(h):
    if " full " in h.lower():
        return "Full"
    elif " sim " in h.lower():
        return "SIM"
    else:
        # Need to check if there are any other scan types (e.g. Thermo and Waters)
        return None


def mode_type_from_header(h):
    if " p " in h.lower():
        return "p"
    elif " c " in h.lower():
        return "c"
    else:
        return None


def count_scan_types(hs):
    return len(set([scan_type_from_header(h) for h in hs]))


def count_ms_types(hs):
    return len(set([ms_type_from_header(h) for h in hs]))


def _partially_overlapping_windows(mzrs):
    """
    Select adjecent windows that partially overlap
    For example: [100-200] and [185-285] (Valid for SIM-stitch)
    """
    assert type(mzrs) == list, "List required"
    temp = []
    for i in range(0, len(mzrs) - 1):
        if mzrs[i][0] < mzrs[i + 1][0] and mzrs[i][1] > mzrs[i + 1][0] and mzrs[i][1] < mzrs[i + 1][1]:
            if mzrs[i] not in temp:
                temp.append(mzrs[i])
            if mzrs[i + 1] not in temp:
                temp.append(mzrs[i+1])
    return temp


def _first_fully_overlapping_windows(mzrs):
    """
    Select windows that fall within another window but do not have identical mass ranges
    For example: [100-200] and [125-175] (Invalid)
    """
    assert type(mzrs) == list, "List required"

    for i in range(0, len(mzrs) - 1):
        if mzrs[i][0] <= mzrs[i + 1][0] and mzrs[i][1] >= mzrs[i + 1][1]:
            return mzrs[i], mzrs[i + 1] # Temporary print
    return []


def _non_overlapping_windows(mzrs):
    """
    Select windows that do not overlap with other windows.
    For example: [100-200] and [200-400] (Valid for merging)
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
        if c == len(mzrs)-1:
            temp.append(mzrs[i])
    return temp


def interpret_experiment(mzrs):

    mzrs.sort(key=lambda x: x[1])

    now = _non_overlapping_windows(mzrs)
    pow = _partially_overlapping_windows(mzrs)

    if len(mzrs) == 1:
        print "Single m/z window....."
        experiment = "single"
    elif len(now) == len(mzrs):
        print "Adjacent m/z windows....."
        experiment = "adjacent"
    elif len(pow) == len(mzrs):
        print "SIM-Stitch experiment - Overlapping m/z windows....."
        experiment = "overlapping"
    else:
        raise IOError("SIM-Stitch cannot be applied; 'filter_scan_events' required or set 'skip_stitching' to False")

    return experiment


def check_metadata(fn_tsv):

    assert os.path.isfile(fn_tsv.encode('string-escape')), "{} does not exist".format(fn_tsv)

    fm = np.genfromtxt(fn_tsv.encode('string-escape'), dtype=None, delimiter="\t", names=True)
    if len(fm.shape) == 0:
        fm = np.array([fm])

    fm_dict = collections.OrderedDict()
    for k in fm.dtype.names:
        fm_dict[k] = list(fm[k])

    unique, counts = np.unique(fm_dict["filename"], return_counts=True)
    if len(unique) != sum(counts):
        raise ValueError("Duplicate filenames in filelist")

    if "replicate" in fm.dtype.names:

        if 0 in fm["replicate"]:
            raise IOError("Incorrect replicate number in filelist. Row {}".format(list(fm["replicate"]).index(0)))

        idxs_replicates = idxs_reps_from_filelist(fm["replicate"])
        counts = {}
        for idxs in idxs_replicates:
            if len(idxs) not in counts:
                counts[len(idxs)] = 1
            else:
                counts[len(idxs)] += 1
        for k, v in counts.items():
            print "{} sample(s) with {} replicate(s)".format(v, k)
    else:
        print "Column for replicate numbers missing. Only required for replicate filter."

    if "batch" in fm.dtype.names:
        unique_batches, counts = np.unique(fm["batch"], return_counts=True)
        print "Batch numbers:", unique_batches
        print "Number of samples in each Batch:", dict(zip(unique_batches, counts))
    else:
        print "Column for batch number missing. Not required."

    if "order" in fm.dtype.names:
        assert np.array_equal(fm["order"], sorted(fm["order"])), "Check the order column - samples not in order"
    else:
        print "Column for sample order missing. Not required."

    if "class" in fm.dtype.names:
        if "replicate" in fm.dtype.names:
            for i in range(len(idxs_replicates)):
                assert len(np.unique(fm["class"][min(idxs_replicates[i]):max(idxs_replicates[i])+1])) == 1, "class names do not match with number of replicates"
        unique, counts = np.unique(fm["class"], return_counts=True)
        cls = dict(zip(unique, counts))
        print "Classes:", cls
    else:
        warnings.warn("Column 'class' for class labels missing.")

    return fm_dict


def update_metadata(peaklists, fl):
    assert isinstance(peaklists[0], PeakList), "PeakList object required"
    for k in fl.keys():  # Update metadata
        for pl in peaklists:
            assert pl.ID in fl[fl.keys()[0]], "filelist and peaklist do not match {}".format(pl.ID)
            index = fl[fl.keys()[0]].index(pl.ID)
            pl.metadata[k] = fl[k][index]
            if "class" in fl.keys():
                if pl.tags.has_tag_type("class_label"):
                    pl.tags.drop_tag_types("class_label")
                pl.tags.add_tags(class_label=fl["class"][index])
    return peaklists


def idxs_reps_from_filelist(replicates):
    idxs, temp = [], [0]
    for i in range(1, len(replicates)):
        if (replicates[i-1] == replicates[i] or replicates[i-1] > replicates[i]) and replicates[i] == 1:
            idxs.append(temp)
            temp = [i]
        elif replicates[i-1] < replicates[i] and replicates[i-1] - replicates[i] == -1:
            temp.append(i)
        else:
            raise ValueError("Incorrect numbering for replicates. Row {}".format(i))
    idxs.append(temp)
    return idxs


def update_class_labels(pm, fn_tsv):

    assert os.path.isfile(fn_tsv.encode('string-escape')), "{} does not exist".format(fn_tsv)

    fm = np.genfromtxt(fn_tsv.encode('string-escape'), dtype=None, delimiter="\t", names=True)
    if len(fm.shape) == 0:
        fm = np.array([fm])

    assert "sample_id" == fm.dtype.names[0] or "filename" == fm.dtype.names[0], "Column for class labels not available"
    assert "class" in fm.dtype.names, "Column for class label (class) not available"
    assert (fm[fm.dtype.names[0]] == pm.peaklist_ids).all(), "Sample ids do not match {}".format(np.setdiff1d(fm[fm.dtype.names[0]], pm.peaklist_ids))
    # TODO: class_labels
    for i in range(len(fm["class"])):
        if pm.peaklist_tags[i].has_tag_type("class_label"):
            pm.peaklist_tags[i].drop_tag_types("class_label")
            pm.peaklist_tags[i].add_tags(class_label=fm["class"][i])
    return pm

