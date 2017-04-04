#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
author(s): Ralf Weber
origin: Nov. 2016
"""
import time
import os
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

from experiment import validate_source
from experiment import validate_metadata


def process_scans(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, min_fraction=None, rsd_thres=None, block_size=2000, ncpus=None):

    files = validate_source(filelist, source)
    fl = validate_metadata(filelist)

    pkls = []
    for i in range(len(files)):

        print
        print files[i]

        scans = read_scan_data(files[i], source, function_noise, nscans, fn_exp)
        scans_re = remove_edges(scans)

        prs = process_replicate_scans(scans_re, snr_thres, ppm, min_fraction, rsd_thres, block_size, ncpus)

        cl = ()
        if "class" in fl: cl = (fl["class"][i],) # TODO: Tags from metadata

        pkl = join_peaklists(os.path.basename(files[i]), prs, cl)
        for k in fl.keys(): pkl.metadata[k] = fl[k][i]
        pkls.append(pkl)
    return pkls


# placeholder (synonym)
def stitch(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres=None, rsd_thres=None, block_size=2000, ncpus=None):
    return process_scans(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres, rsd_thres, block_size, ncpus)


def replicate_filter(peaklists, ppm, reps, minpeaks, rsd_thres=None, filelist=None, block_size=2000, ncpus=None):

    if filelist is not None:
        files = validate_source(filelist, peaklists)
        fl = validate_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in files]
        for k in fl.keys():  # Update metadata
            for pl in peaklists:
                pl.metadata[k] = fl[k][files.index(pl.ID)]

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

        pkl = pm.to_peaklist(ID=merged_id)
        pkl.add_tags(class_name=pls[0].metadata.tags[0]) # TODO: fix raw parser to remove tags from metadata
        pkl.add_attribute("present", pm.present)
        pkl.add_attribute("rsd", pm.rsd)
        pkl.add_attribute("present_flag", pm.present >= minpeaks, is_flag=True)
        if rsd_thres is not None: pkl.add_attribute("rsd_flag", pm.rsd <= rsd_thres, flagged_only=False, is_flag=True)
        for k in pls[0].metadata:
            if k != "filename": pkl.metadata[k] = pls[0].metadata[k]
        pkls_rep_filt.append(pkl)
    return pkls_rep_filt


def align_samples(peaklists, ppm, filelist=None, block_size=2000, ncpus=None):

    assert filelist is None or os.path.isfile(filelist), "Provide list of peaklists"
    assert type(peaklists) == list, "Provide list of peaklists"

    if filelist is not None:
        files = validate_source(filelist, peaklists)
        fl = validate_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in files]
        for k in fl.keys():  # Update metadata
            for pl in peaklists:
                pl.metadata[k] = fl[k][files.index(pl.ID)]
    return align_peaks(peaklists, ppm=ppm, block_size=block_size, byunique=False, ncpus=ncpus)


def blank_filter(peak_matrix, blank_label, min_fraction=1.0, min_fold_change=1.0, function="mean", rm_samples=True):
    return filter_blank_peaks(peak_matrix, blank_label, min_fraction, min_fold_change, function, rm_samples)


def sample_filter(peak_matrix, min_fraction, within=False, rsd=None, qc_label=None, fn_classes=None):
    if not within:
        peak_matrix = filter_across_classes(peak_matrix, min_fraction)
    elif within:
        peak_matrix = filter_within_classes(peak_matrix, None, min_fraction) # TODO: use tag_type instead of None, currently it's temp workaround
    if rsd is not None:
        peak_matrix = filter_rsd(peak_matrix, rsd, qc_label)
    return peak_matrix


