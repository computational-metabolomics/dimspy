#!/usr/bin/python
# -*- coding: utf-8 -*-


import collections
import os
import h5py
import numpy as np
from dimspy.models.peaklist import PeakList
from dimspy.models.peak_matrix import PeakMatrix
from portals import hdf5_portal
from portals.paths import check_paths
from portals.txt_portal import load_peak_matrix_from_txt
from experiment import check_metadata
from experiment import update_class_labels
from experiment import update_metadata
from process.peak_alignment import align_peaks
from process.peak_filters import filter_fraction
from process.peak_filters import filter_blank_peaks
from process.peak_filters import filter_rsd
from process.scan_processing import average_replicate_scans
from process.scan_processing import join_peaklists
from process.scan_processing import read_scans
from process.scan_processing import remove_edges


def process_scans(source, function_noise, snr_thres, nscans, ppm, min_fraction=None, rsd_thres=None, filelist=None, subset_scan_events=None, block_size=2000, ncpus=None):

    filenames = check_paths(filelist, source)
    if len([fn for fn in filenames if not fn.lower().endswith(".mzml") or not fn.lower().endswith(".raw")]) == 0:
        raise IOError("Incorrect file format. Provide .mzML and .raw files")

    if filelist is not None:
        fl = check_metadata(filelist)
    else:
        fl = collections.OrderedDict()

    pls = []
    for i in range(len(filenames)):

        print
        print os.path.basename(filenames[i])

        scans = read_scans(filenames[i], source, function_noise, nscans, subset_scan_events)
        scans_er = remove_edges(scans)

        prs = average_replicate_scans(scans_er, snr_thres, ppm, min_fraction, rsd_thres, block_size, ncpus)
        pl = join_peaklists(os.path.basename(filenames[i]), prs)

        if "class" in fl:
            pl.add_tags(class_label=fl["class"][i])
            # pl.add_tags(class_label2=fl["class"][i])

        for k in fl.keys():
            pl.metadata[k] = fl[k][i]

        pls.append(pl)
    return pls


# placeholder (synonym)
def stitch(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres=None, rsd_thres=None, block_size=2000, ncpus=None):
    return process_scans(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres, rsd_thres, block_size, ncpus)


def replicate_filter(source, ppm, reps, min_peaks, rsd_thres=None, filelist=None, block_size=2000, ncpus=None):

    filenames = check_paths(filelist, source)
    if len(filenames) == 0:
        raise IOError("Provide a filelist that list all the text files (columnname:filename) and assign replicate numbers to each filename/sample (columnname:replicate)")
    peaklists = hdf5_portal.load_peaklists_from_hdf5(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata(peaklists, fl)

    if not hasattr(peaklists[0].metadata, "replicate"):
        raise IOError("Provide a filelist and assign replicate numbers (columnname:replicate) to each filename/sample")

    unique, counts = np.unique([pl.metadata.replicate for pl in peaklists], return_counts=True)
    if max(unique) != reps:
        raise ValueError("replicates missing (1)")
    if len(counts) != reps:
        raise ValueError("replicates missing (2)")
    if sum(counts) != len(peaklists):
        raise ValueError("replicates missing (3)")
    if list(unique) != range(1, reps+1):
        raise ValueError("replicates missing (4)")
    if len(counts) <= 1:
        raise ValueError("No technical replicates available (single) - Skip replicate filter")

    idx_peaklists = range(0, len(peaklists) + reps, reps)
    pls_rep_filt = []
    for i in range(len(idx_peaklists)-1):

        pls = peaklists[idx_peaklists[i]:idx_peaklists[i+1]]
        pm = align_peaks(pls, ppm, block_size, ncpus=ncpus)
        #############################################################
        # TODO: Write some sort of merging function for multiple replicate file names
        #############################################################
        prefix = os.path.commonprefix([p.ID for p in pls])
        merged_id = "{}{}".format(prefix, "_".join(map(str, [p.ID.replace(prefix, "").split(".")[0] for p in pls])))
        #############################################################

        pl = pm.to_peaklist(ID=merged_id)
        pl.add_tags(*pls[0].tag_of(None), **{t: pls[0].tag_of(t) for t in pls[0].tag_types})
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
    peaklists = hdf5_portal.load_peaklists_from_hdf5(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata(peaklists, fl)

    return align_peaks(peaklists, ppm=ppm, block_size=block_size, ncpus=ncpus)


def blank_filter(peak_matrix, blank_label, min_fraction=1.0, min_fold_change=1.0, function="mean", rm_samples=True, class_labels=None):

    if min_fraction < 0.0 or min_fraction > 1.0:
        raise ValueError("Provide a value between 0. and 1.")
    if min_fold_change < 0:
        raise ValueError("Provide a value larger than zero.")
    if function not in ("mean", "median", "max"):
        raise ValueError("Mean, median or max intensity")
    if not os.path.isfile(peak_matrix):
        raise IOError("{} does not exist".format(peak_matrix))

    if h5py.is_hdf5(peak_matrix):
        peak_matrix = hdf5_portal.load_peak_matrix_from_hdf5(peak_matrix)
    else:
        peak_matrix = load_peak_matrix_from_txt(peak_matrix)

    if class_labels is not None:
        peak_matrix = update_class_labels(peak_matrix, class_labels)

    if blank_label not in peak_matrix.peaklist_tag_values:
        raise IOError("Blank label ({}) does not exist".format(blank_label))

    return filter_blank_peaks(peak_matrix, blank_label, min_fraction, min_fold_change, function, rm_samples)


def sample_filter(peak_matrix, min_fraction, within=False, rsd=None, qc_label=None, class_labels=None):

    if not os.path.isfile(peak_matrix):
        raise IOError("{} does not exist".format(peak_matrix))

    if h5py.is_hdf5(peak_matrix):
        peak_matrix = hdf5_portal.load_peak_matrix_from_hdf5(peak_matrix)
    else:
        peak_matrix = load_peak_matrix_from_txt(peak_matrix)

    if class_labels is not None:
        if not os.path.isfile(class_labels):
            raise IOError("{} does not exist".format(class_labels))
        peak_matrix = update_class_labels(peak_matrix, class_labels)

    if qc_label is not None:
        if qc_label not in peak_matrix.peaklist_tag_values:
            raise IOError("QC label ({}) does not exist".format(qc_label))

    peak_matrix = filter_fraction(peak_matrix, min_fraction, within_classes=within, class_tag_type=None)

    if rsd is not None:
        peak_matrix = filter_rsd(peak_matrix, rsd, qc_label)
    return peak_matrix


def hdf5_to_text(fname, path_out, separator="\t", transpose=False):
    assert os.path.isfile(fname), 'HDF5 database [%s] not exists' % fname
    assert h5py.is_hdf5(fname), 'input file [%s] is not a valid HDF5 database' % fname
    seps = {"comma": ",", "tab": "\t"}
    if separator in seps: separator = seps[separator]
    assert separator in [",", "\t"], "Incorrect separator ('tab', 'comma', ',', '\t')"
    f = h5py.File(fname, 'r')
    if "mz" in f:
        obj = hdf5_portal.load_peak_matrix_from_hdf5(fname)
        assert isinstance(obj, PeakMatrix)
        obj = hdf5_portal.load_peak_matrix_from_hdf5(fname)
        with open(os.path.join(path_out), "w") as pk_out:
            pk_out.write(obj.to_str(separator, transpose))
    else:
        assert os.path.isdir(path_out), "File or Directory does not exist:".format(path_out)
        obj = hdf5_portal.load_peaklists_from_hdf5(fname)
        assert isinstance(obj[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
        for pl in obj:
            with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                pk_out.write(pl.to_str(separator))
    return