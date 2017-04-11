#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import numpy as np
import collections
from process.peak_alignment import align_peaks
from process.scan_processing import read_scans
from process.scan_processing import average_replicate_scans
from process.scan_processing import remove_edges
from process.scan_processing import join_peaklists
from process.peak_filters import filter_blank_peaks
from process.peak_filters import filter_across_classes
from process.peak_filters import filter_within_classes
from process.peak_filters import filter_rsd
from portals import check_paths
from experiment import check_metadata
from experiment import update_metadata
from experiment import update_class_labels
from portals import load_peaklists


def process_scans(source, function_noise, snr_thres, nscans, ppm, min_fraction=None, rsd_thres=None, filelist=None, subset_mzrs=None, block_size=2000, ncpus=None):

    filenames = check_paths(filelist, source)
    assert len([fn for fn in filenames if not fn.lower().endswith(".mzml") or not fn.lower().endswith(".raw")]) > 0, "Incorrect file format. Provide .mzML and .raw files"

    if filelist is not None:
        fl = check_metadata(filelist)
    else:
        fl = collections.OrderedDict()

    pls = []
    for i in range(len(filenames)):

        print
        print os.path.basename(filenames[i])

        scans = read_scans(filenames[i], source, function_noise, nscans, subset_mzrs)
        scans_er = remove_edges(scans)

        prs = average_replicate_scans(scans_er, snr_thres, ppm, min_fraction, rsd_thres, block_size, ncpus)
        pl = join_peaklists(os.path.basename(filenames[i]), prs)

        if "class" in fl:
            pl.add_tags(class_label=fl["class"][i])  # TODO: Tags from metadata
            # pl.add_tags(class_label2=fl["class"][i]) # TODO: Tags from metadata
            # TODO: assert not any(map(lambda x: x in self.tag_values, list(args) + kwargs.values())), 'tag already exists'

        for k in fl.keys():
            pl.metadata[k] = fl[k][i]

        pls.append(pl)
    return pls


# placeholder (synonym)
def stitch(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres=None, rsd_thres=None, block_size=2000, ncpus=None):
    return process_scans(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres, rsd_thres, block_size, ncpus)


def replicate_filter(source, ppm, reps, min_peaks, rsd_thres=None, filelist=None, block_size=2000, ncpus=None):

    filenames = check_paths(filelist, source)
    assert len(filenames) > 0, "Provide a filelist that list all the text files (columnname:filename) and assign replicate numbers to each filename/sample (columnname:replicate)"
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata(peaklists, fl)

    assert hasattr(peaklists[0].metadata, "replicate"), "Provide a filelist and assign replicate numbers (columnname:replicate) to each filename/sample"

    unique, counts = np.unique([pl.metadata.replicate for pl in peaklists], return_counts=True)
    assert max(unique) == reps, "replicates missing (1)"
    assert len(counts) == reps, "replicates missing (2)"
    assert sum(counts) == len(peaklists), "replicates missing (3)"
    assert list(unique) == range(1, reps+1), "replicates missing (4)"
    assert len(counts) > 1, "No technical replicates available (single) - Skip replicate filter "

    idx_peaklists = range(0, len(peaklists) + reps, reps)
    pls_rep_filt = []
    for i in range(len(idx_peaklists)-1):

        pls = peaklists[idx_peaklists[i]:idx_peaklists[i+1]]
        pm = align_peaks(pls, ppm, block_size, byunique=False, ncpus=ncpus)
        #############################################################
        # TODO: Write some sort of merging function for multiple replicate file names
        #############################################################
        prefix = os.path.commonprefix([p.ID for p in pls])
        merged_id = "{}{}".format(prefix, "_".join(map(str, [p.ID.replace(prefix, "").split(".")[0] for p in pls])))
        #############################################################

        pl = pm.to_peaklist(ID=merged_id)
        for j, t in enumerate(pls[0].tag_types):  # TODO:
            pl.add_tags(t, pls[0].tag_values[j])
        pl.add_attribute("present", pm.present)
        pl.add_attribute("rsd", pm.rsd)
        pl.add_attribute("present_flag", pm.present >= min_peaks, is_flag=True)
        if rsd_thres is not None:
            pl.add_attribute("rsd_flag", pm.rsd <= rsd_thres, flagged_only=False, is_flag=True)
        for k in pls[0].metadata:
            if k != "filename":
                pl.metadata[k] = pls[0].metadata[k]
        pls_rep_filt.append(pl)
    return pls_rep_filt


def align_samples(source, ppm, filelist=None, block_size=2000, ncpus=None):

    filenames = check_paths(filelist, source)
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata(peaklists, fl)

    return align_peaks(peaklists, ppm=ppm, block_size=block_size, byunique=False, ncpus=ncpus)


def blank_filter(peak_matrix, blank_label, min_fraction=1.0, min_fold_change=1.0, function="mean", rm_samples=True, tsv_labels=None):

    assert 0 < min_fraction <= 1, "Provide a value between 0. and 1."
    assert min_fold_change >= 0, "Provide a value larger than zero."
    assert function in ("mean", "median", "max"), "Mean, median or max intensity"

    if tsv_labels is not None:
        peak_matrix = update_class_labels(peak_matrix, tsv_labels)

    assert blank_label in peak_matrix.peaklist_tag_values, "Blank label ({}) does not exist".format(blank_label)

    return filter_blank_peaks(peak_matrix, blank_label, min_fraction, min_fold_change, function, rm_samples)


def sample_filter(peak_matrix, min_fraction, within=False, rsd=None, qc_label=None, tsv_labels=None):

    if tsv_labels is not None:
        assert os.path.isfile(tsv_labels), "File with class labels not available"
        peak_matrix = update_class_labels(peak_matrix, tsv_labels)

    if qc_label is not None:
        assert qc_label in peak_matrix.peaklist_tag_values

    if not within:
        peak_matrix = filter_across_classes(peak_matrix, min_fraction)
    elif within:
        peak_matrix = filter_within_classes(peak_matrix, None, min_fraction)  # TODO: use tag_type instead of None, currently it's temp workaround
    if rsd is not None:
        peak_matrix = filter_rsd(peak_matrix, rsd, qc_label)
    return peak_matrix
