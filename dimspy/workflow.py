#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
author(s): Ralf Weber
origin: Nov. 2016
"""
import time
import os
import sys
import zipfile
import cPickle as pickle
import numpy as np

from models.peak_matrix import PeakMatrix
from models.peaklist import PeakList

from process.peak_alignment import align_peaks
from process.scan_processing import read_scan_data
from process.scan_processing import process_replicate_scans
from process.scan_processing import remove_edges
from process.scan_processing import join_peaklists

from process.peak_filters import filter_blank_peaks
from process.peak_filters import filter_across_classes
from process.peak_filters import filter_within_classes
from process.peak_filters import filter_rsd

from experiment import read_filelist


def process_scans(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, min_fraction=None, rsd_thres=None, block_size=2000, ncpus=None):

    if filelist is None and os.path.isdir(source):
        fl = {"filename":[os.path.join(source,fn) for fn in os.listdir(source) if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]}
    elif filelist is None and zipfile.is_zipfile(source):
        with zipfile.ZipFile(source) as zf:
            assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported"
            fl = {"filename": [fn for fn in zf.namelist() if fn.lower().endswith(".mzml")]}
    elif os.path.isfile(filelist) and zipfile.is_zipfile(source):
        with zipfile.ZipFile(source) as zf:
            assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported"
        fl = read_filelist(filelist, source)
    elif os.path.isfile(filelist) and os.path.isdir(source):
        fl = read_filelist(filelist, source)
    else:
        print "Can not read and parse {} and {}".format(source, filelist)
        sys.exit()

    pkls = []
    for i in range(len(fl["filename"])):

        print
        print fl["filename"][i]

        scans = read_scan_data(fl["filename"][i], source, function_noise, nscans, fn_exp)
        scans_re = remove_edges(scans)

        prs = process_replicate_scans(scans_re, snr_thres, ppm, min_fraction, rsd_thres, block_size, ncpus)

        cl = ()
        if "class" in fl: cl = (fl["class"][i],)

        pkl = join_peaklists(os.path.basename(fl["filename"][i]), prs, cl)
        for k in fl.keys(): pkl.metadata[k] = fl[k][i]
        pkls.append(pkl)
    return pkls


# placeholder (synonym)
def stitch(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres=None, rsd_thres=None, block_size=2000, ncpus=None):
    return process_scans(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres, rsd_thres, block_size, ncpus)


def replicate_filter(peaklists, ppm, reps, minpeaks, rsd_thres=None, block_size=2000, ncpus=None):

    assert type(peaklists) == list, "Provide list of peaklists"

    unique, counts = np.unique([pl.metadata.replicate for pl in peaklists], return_counts=True)
    assert max(unique) == reps, "replicates missing (1)"
    assert len(counts) == reps, "replicates missing (2)"
    assert sum(counts) == len(peaklists), "replicates missing (3)"
    assert list(unique) == range(1, reps+1), "replicates missing (4)"
    assert len(counts) > 1, "No technical replicates available (single) - No need to run replicate filter "

    idx_peaklists = range(0, len(peaklists) + reps, reps)
    pkls_rep_filt = []
    for i in range(len(idx_peaklists)-1):

        pls = peaklists[idx_peaklists[i]:idx_peaklists[i+1]]
        pm = align_peaks(pls, ppm, block_size, byunique=False, ncpus=ncpus)
        #############################################################
        # Write some sort of merging function for multiple replicate file names
        #############################################################
        prefix = os.path.commonprefix([p.ID for p in pls])
        merged_id = "{}{}".format(prefix, "_".join(map(str, [p.ID.replace(prefix, "").split(".")[0] for p in pls])))
        #############################################################

        pkl = pm.to_peaklist(ID = merged_id)
        pkl.add_tags(class_name = pls[0].metadata.tags[0]) # TODO: fix raw parser to remove tags from metadata
        pkl.add_attribute("present", pm.present)
        pkl.add_attribute("rsd", pm.rsd)
        pkl.add_attribute("present_flag", pm.present >= minpeaks, is_flag=True)
        if rsd_thres is not None: pkl.add_attribute("rsd_flag", pm.rsd <= rsd_thres, flagged_only=False, is_flag=True)
        for k in pls[0].metadata:
            if k != "filename": pkl.metadata[k] = pls[0].metadata[k]
        pkls_rep_filt.append(pkl)
    return pkls_rep_filt


def align_samples(peaklists, ppm, block_size=2000, ncpus=None):
    return align_peaks(peaklists, ppm=ppm, block_size=block_size, byunique=False, ncpus=ncpus)


def blank_filter(peak_matrix, blank_label, min_fraction=1.0, min_fold_change=1.0, function="mean", rm_samples=True):
    return filter_blank_peaks(peak_matrix, blank_label, min_fraction, min_fold_change, function, rm_samples)


def sample_filter(peak_matrix, min_fraction, within=False, rsd=None, qc_label=None):

    if not within:
        peak_matrix = filter_across_classes(peak_matrix, min_fraction)
    elif within:
        peak_matrix = filter_within_classes(peak_matrix, None, min_fraction) # TODO: use tag_type instead of None, currently it's temp workaround
    if rsd is not None:
        peak_matrix = filter_rsd(peak_matrix, rsd, qc_label)
    return peak_matrix


def to_readable(pickle_file, path_out, separator, transpose=False):
    assert os.path.isfile(pickle_file), "Pickle file does not exist"
    assert separator in ["tab", "comma"], "Incorrect separator [tab, comma]"
    seps = {"comma":",", "tab":"\t"}
    with open(pickle_file, "rb") as fn_pkl_in:
        pkl = pickle.load(fn_pkl_in)
        if type(pkl) == list:
            assert isinstance(pkl[0], PeakList), "Not compatible with {}".format(type(pkl[0]))
            for pl in pkl:
                with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                    pk_out.write(pl.to_str(seps[separator]))
                time.sleep(1)

        elif isinstance(pkl, PeakMatrix):
            assert os.path.isfile(path_out), "Provide filename for peak matrix"
            with open(os.path.join(path_out), "w") as pk_out:
                pk_out.write(pkl.to_str(seps[separator], transpose))
    return


