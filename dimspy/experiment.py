#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Define DIMS experiment

author(s): Ralf Weber
origin: Oct. 2016
"""
import os
import warnings
import collections
import re
import copy
import numpy as np
import zipfile


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


def sort_mz_ranges(mz_ranges):
    return collections.OrderedDict(sorted(mz_ranges.iteritems(), key=lambda y: y[1][0]))


def _partially_overlapping_windows(mzrs):
    """
    Select adjecent windows that partially overlap
    For example: [100-200] and [185-285] (Valid for SIM-stitch)
    """
    assert type(mzrs) == collections.OrderedDict, "OrderedDict requried"
    temp = []
    mzrs_val = mzrs.values()
    for i in range(0, len(mzrs_val) - 1):
        if mzrs_val[i][0] < mzrs_val[i + 1][0] and mzrs_val[i][1] > mzrs_val[i + 1][0] and mzrs_val[i][1] < mzrs_val[i + 1][1]:
            if mzrs.keys()[i] not in temp:
                temp.append(mzrs.keys()[i])
            if mzrs.keys()[i + 1] not in temp:
                temp.append(mzrs.keys()[i+1])
    return temp


def _first_fully_overlapping_windows(mzrs):
    """
    Select windows that fall within another window but do not have identical mass ranges
    For example: [100-200] and [125-175] (Invalid)
    """
    assert type(mzrs) == collections.OrderedDict, "OrderedDict requried"
    mzrs_val= mzrs.values()
    for i in range(0, len(mzrs_val) - 1):
        if mzrs_val[i][0] <= mzrs_val[i + 1][0] and mzrs_val[i][1] >= mzrs_val[i + 1][1]:
            return mzrs.keys()[i], mzrs.keys()[i + 1] # Temporary print
    return []


def _non_overlapping_windows(mzrs):
    """
    Select windows that do not overlap with other windows.
    For example: [100-200] and [200-400] (Valid for merging)
    """
    assert type(mzrs) == collections.OrderedDict, "OrderedDict requried"
    temp = []
    mzrs_val= mzrs.values()
    for i in range(0, len(mzrs_val)):
        c = 0
        for j in range(0, len(mzrs_val)):
            if mzrs_val[i][0] <= mzrs_val[j][0] and mzrs_val[i][1] <= mzrs_val[j][0]:
                c += 1
            elif mzrs_val[i][0] >= mzrs_val[j][1] and mzrs_val[i][1] >= mzrs_val[j][1]:
                c += 1
        if c == len(mzrs_val)-1:
            temp.append(mzrs.keys()[i])
    return temp


def _match_header_description(h, v, e):

    if v == [e["start"], e["end"]]:
        if "ms_type" in e:
            if e["ms_type"].lower() != ms_type_from_header(h).lower():
                return False
        if "scan_type" in e:
            if e["scan_type"].lower() != scan_type_from_header(h).lower():
                return False
        if "mode" in e:
            if e["mode"].lower() != mode_type_from_header(h).lower():
                return False
        return True
    return False


def remove_headers(exp, mzrs):

    for w in exp:

        if "start" in w:
            assert type(w["start"]) == float, "Wrong format for start"
        if "end" in w:
            assert type(w["end"]) == float, "Wrong format for end"
        if "mode" in w:
            assert w["mode"].lower() in ["p", "c"], "Wrong mode (p or c)"
        if "scan_type" in w:
            assert w["scan_type"].lower() in ["sim", "full"], "Scan type does not exist"
        if "ms_type" in w:
            assert w["ms_type"].lower() in ["ftms", "itms"], "MS type does not exist"

    for h, r in copy.copy(mzrs).items():
        incl = False
        for e in exp:
            if _match_header_description(h, r, e):
                incl = True
                break
        if not incl:
            del mzrs[h]

    assert len(mzrs) == len(exp), "No matching window descriptions"
    return mzrs


def read_ms_exp_from_file(fn_tsv):

    fn_tsv = fn_tsv.encode('string-escape')
    assert os.path.isfile(fn_tsv), "{} does not exist".format(fn_tsv)

    with open(fn_tsv) as fn_exp:

        exp = []

        lines =  fn_exp.readlines()
        header = [name.lower() for name in lines[0].rstrip().split("\t")]

        assert len(lines) >= 2, "too few rows"
        assert len(header) <= 6 , "too many columns"
        assert len(header) >= 2 , "too few columns"
        assert "start" in header, "provide start value for each window"
        assert "end" in header, "provide end value for each window"
        assert False not in [False for h in header if h not in ["start", "end", "scan_type", "ms_type", "mode"]], "unknown column name"

        for line in lines[1:]:

            d = collections.OrderedDict(zip(header, line.rstrip().lower().split("\t")))

            if "start" in d:
                d["start"] = float(d["start"])
            if "end" in d:
                d["end"] = float(d["end"])
            if "mode" in d:
                assert d["mode"].lower() in ["p", "c"], "wrong mode (p or c)"
            if "scan_type" in d:
                assert d["scan_type"].lower() in ["sim", "full"], "scan type does not exist"
            if "ms_type" in d:
                assert d["ms_type"].lower() in ["ftms", "itms"], "ms type does not exist"
            exp.append(d)
        return exp
    return []


def read_ms_exp_from_headers(mz_ranges):

    mzrs = sort_mz_ranges(mz_ranges)

    now = _non_overlapping_windows(mzrs)
    pow = _partially_overlapping_windows(mzrs)
    ffow = _first_fully_overlapping_windows(mzrs)

    #print len(mzrs), len(now), len(pow), len(ffow)

    if len(now) == len(mzrs):
        print "Experiment: Join non-overlapping windows"
    elif len(pow) == len(mzrs):
        print "Experiment: Stitch overlapping windows"
    elif len(ffow) > 0:
        del mzrs[ffow[0]]
        pow2 = _partially_overlapping_windows(mzrs)
        if len(pow2) == len(mzrs) - 1:
            print "Experiment: Stitch overlapping windows (Fully overlapping window removed)"
        else:
            print "Please, describe DIMS experiment - experiment too complex"
    else:
        print "Please, describe DIMS experiment - experiment too complex"

    return mzrs.keys()


def read_filelist(fn_tsv, source):

    source = source.encode('string-escape')
    fn_tsv = fn_tsv.encode('string-escape')

    assert os.path.isfile(fn_tsv), "{} does not exist".format(fn_tsv)
    assert zipfile.is_zipfile(source) or os.path.isdir(source), "path or .zip archive does not exist"

    fm = np.genfromtxt(fn_tsv, dtype=None, delimiter="\t", names=True)
    if len(fm.shape) == 0: # TODO: Added to check if filelist has a single row
        fm = np.array([fm])

    col_names = fm.dtype.fields.keys()

    fm_dict = collections.OrderedDict()
    for k in col_names:
        fm_dict[k] = list(fm[k])

    if zipfile.is_zipfile(source):
        l = zipfile.ZipFile(source)
        for fn in fm["filename"]:
            assert fn in l.namelist(), "{} does not exist in .zip file".format(fn)

    elif os.path.isdir(source):
        l = os.listdir(source)
        for i in range(len(fm["filename"])):
            assert os.path.basename(fm["filename"][i]) in l, "{} does not exist in directory".format(os.path.basename(fm["filename"][i]))
            fn = str(fm["filename"][i])
            if os.path.basename(fn) == fm["filename"][i]:
                fm_dict["filename"][i] = os.path.join(source, fn).replace('\\', r'\\')

    #if "blank" not in col_names and "Blank" not in col_names:
    #    warnings.warn("No samples marked as blank. Column missing.")
    #else:
    #    unique, counts = np.unique(fm["blank"], return_counts=True)
    #    print "Blank samples:", counts

    #if "qc" not in col_names and "QC" not in col_names:
    #    warnings.warn("No samples marked as QC. Column missing.")
    #else:
    #    unique, counts = np.unique(fm["qc"], return_counts=True)
    #    print "QC samples:", counts

    if "replicate" in col_names:
        unique_reps, counts = np.unique(fm["replicate"], return_counts=True)
        assert len(np.unique(counts)) == 1, "Incorrect numbering of replicates"
        assert len(unique_reps) == max(unique_reps), "Incorrect numbering of replicates"
        print "Replicates:", dict(zip(unique_reps, counts))
    else:
        print "Column for replicate numbers missing"

    if "batch" in col_names:
        unique_batches, counts = np.unique(fm["batch"], return_counts=True)
        #assert np.array_equal(fm["batch"], sorted(fm["batch"])), "mixed order of batches"
        print "Batch numbers:", unique_batches
        print "Number of samples in each Batch:", dict(zip(unique_batches, counts))
    else:
        print "Column for batch number missing. Not required."

    if "order" in col_names:
        assert np.array_equal(fm["order"], sorted(fm["order"])), "Check the order column - samples not in order"
    else:
        print "Column for sample order missing. Not required."

    if "class" in col_names:

        ra = range(0, len(fm["class"]), max(unique_reps))
        for i in range(0, len(ra)-1):
            assert len(np.unique(fm["class"][ra[i]:ra[i+1]])) == 1, "class names do not match with number of replicates"

        unique, counts = np.unique(fm["class"], return_counts=True)

        cls = dict(zip(unique, counts))
        print "Classes:", cls
    else:
        warnings.warn("Column for class labels missing.")

    return fm_dict


if __name__ == '__main__':

    print

    """
    fn_mzml = "../data/Pol_pos_Lung_FeOx_QC1_Rep_3.mzML"
    msrun = pymzml.run.Reader(fn_mzml)

    scan_ids, mz_ranges = collections.OrderedDict(), collections.OrderedDict()
    for scan in msrun:
        if "MS:1000512" in scan: # Need to check how the header for a waters files looks like
            scan_ids.setdefault(scan['MS:1000512'], []).append(scan['id'])
            mzr = mz_range_from_header(scan['MS:1000512'])
            mz_ranges.setdefault(scan['MS:1000512'], mzr)

    print
    print "TEST - TEMPLATE PROVIDED"

    # EXPERIMENT TEXT FILE PROVIDED
    mzrs = sort_mz_ranges(mz_ranges)
    experiment = read_ms_experiment("exp_design_stitch_SIM.txt")
    print len(mzrs), mzrs
    mzrs = remove_headers(experiment, mzrs)
    print len(mzrs), mzrs

    # EXPERIMENT TEXT FILE PROVIDED
    mzrs = sort_mz_ranges(mz_ranges)
    experiment = read_ms_experiment("exp_design_full_simple.txt")
    print len(mzrs), mzrs
    mzrs = remove_headers(experiment, mzrs)
    print len(mzrs), mzrs

    # EXPERIMENT TEXT FILE PROVIDED
    #mzrs = sort_mz_ranges(mz_ranges)
    #experiment = read_experiment("exp_design_full_simple_wrong_scan_type.txt")
    #print len(mzrs), mzrs
    #mzrs = remove_headers(experiment, mzrs)
    #print len(mzrs), mzrs




    print
    print "TEST - WITHOUT TEMPLATE (1)"
    fn_mzml = "../data/Pol_pos_Lung_FeOx_QC1_Rep_3.mzML"
    msrun = pymzml.run.Reader(fn_mzml)

    scan_ids, mz_ranges = collections.OrderedDict(), collections.OrderedDict()
    for scan in msrun:
        if "MS:1000512" in scan:  # Need to check how the header for a waters files looks like
            scan_ids.setdefault(scan['MS:1000512'], []).append(scan['id'])
            mzr = mz_range_from_header(scan['MS:1000512'])
            mz_ranges.setdefault(scan['MS:1000512'], mzr)

    mzrs = sort_mz_ranges(mz_ranges)

    # ------
    #     -------
    #          -------

    now = non_overlapping_windows(mzrs)
    pow = partially_overlapping_windows(mzrs)
    ffow = first_fully_overlapping_windows(mzrs)
    print len(mzrs), len(now), len(pow), len(ffow)

    if len(now) == len(mzrs):
        print "READY TO PROCESS - MERGE"
        print
    elif len(pow) == len(mzrs):
        print "READY TO PROCESS - SIM STITCH"
        print
    elif len(ffow) > 0:
        print "REMOVE FULLY OVERLAPPING WINDOW"
        del mzrs[ffow[0]]
        pow2 = partially_overlapping_windows(mzrs)
        if len(pow2) == len(mzrs) - 1:
            print "READY TO PROCESS - SIM STITCH (AFTER OVERLAPPING WINDOW REMOVED (LIKELY FULL MS)"
            print
        else:
            print "Please, describe DIMS experiment - experiment too complex"
            print
    else:
        print "Please, describe DIMS experiment - experiment too complex"
        print




    print
    print "TEST - WITHOUT TEMPLATE (2)"
    mzrs = sort_mz_ranges(mz_ranges)
    del mzrs['FTMS + p NSI SIM ms [215.00-290.00]'] # TO CREAT
    # FOR EXAMPLE (TWO SECTIONS WITH OVERLAPPING WINDOWS)
    # ------
    #     -------
    #          -------
    #                       -------
    #                            -------
    now = non_overlapping_windows(mzrs)
    pow = partially_overlapping_windows(mzrs)
    ffow = first_fully_overlapping_windows(mzrs)
    print len(mzrs), len(now), len(pow), len(ffow)

    if len(now) == len(mzrs):
        print "READY TO PROCESS - MERGE"
        print
    elif len(pow) == len(mzrs):
        print "READY TO PROCESS - SIM STITCH"
        print
    elif len(ffow) > 0:
        print "REMOVE FULLY OVERLAPPING WINDOW"
        del mzrs[ffow[0]]
        now2 = non_overlapping_windows(mzrs)
        pow2 = partially_overlapping_windows(mzrs)
        ffow2 = first_fully_overlapping_windows(mzrs)
        if len(pow2) == len(mzrs) - 1:
            print "READY TO PROCESS - SIM STITCH (AFTER OVERLAPPING WINDOW REMOVED (LIKELY FULL MS)"
            print
        elif len(pow2) < len(mzrs) - 1 and len(now2) == 0 and len(ffow2) == 0:
            print "READY TO PROCESS - SIM STITCH TWO SECTIONS (AFTER OVERLAPPING WINDOW REMOVED (LIKELY FULL MS)"
            print
        else:
            print "Please, describe DIMS experiment - experiment too complex"
            print
    else:
        print "Please, describe DIMS experiment - experiment too complex"
        print

    print "TEST - WITHOUT TEMPLATE (3)"
    # TO COMPLEX TO PROCESS
    fn_mzml = "../data/SimStitchOptimisation_megamix_1e6_1.mzML"
    msrun = pymzml.run.Reader(fn_mzml)

    scan_ids, mz_ranges = collections.OrderedDict(), collections.OrderedDict()
    for scan in msrun:
        if "MS:1000512" in scan:  # Need to check how the header for a waters files looks like
            scan_ids.setdefault(scan['MS:1000512'], []).append(scan['id'])
            mzr = mz_range_from_header(scan['MS:1000512'])
            mz_ranges.setdefault(scan['MS:1000512'], mzr)

    mzrs = sort_mz_ranges(mz_ranges)
    now = non_overlapping_windows(mzrs)
    pow = partially_overlapping_windows(mzrs)
    ffow = first_fully_overlapping_windows(mzrs)
    print len(mzrs), len(now), len(pow), len(ffow)

    if len(now) == len(mzrs):
        print "READY TO PROCESS - MERGE"
        print
    elif len(pow) == len(mzrs):
        print "READY TO PROCESS - SIM STITCH"
        print
    elif len(ffow) > 0:
        print "REMOVE FULLY OVERLAPPING WINDOW"
        del mzrs[ffow[0]]
        pow2 = partially_overlapping_windows(mzrs)
        if len(pow2) == len(mzrs) - 1:
            print "READY TO PROCESS - SIM STITCH (AFTER OVERLAPPING WINDOW REMOVED (LIKELY FULL MS)"
            print
        else:
            print "Please, describe DIMS experiment - experiment too complex"
            print
    else:
        print "Please, describe DIMS experiment - experiment too complex"
        print
"""