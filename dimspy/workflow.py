#!/usr/bin/python
# -*- coding: utf-8 -*-


import os
import logging
import operator
import collections
import h5py
import numpy as np
from itertools import combinations
import zipfile
from models.peaklist import PeakList
from models.peak_matrix import PeakMatrix
from portals import hdf5_portal
from portals import txt_portal
from portals.paths import check_paths
from experiment import check_metadata
from experiment import update_class_labels
from experiment import update_metadata
from experiment import idxs_reps_from_filelist
from process.peak_alignment import align_peaks
from process.peak_filters import filter_fraction
from process.peak_filters import filter_blank_peaks
from process.peak_filters import filter_rsd
from process.peak_filters import filter_mz_ranges
from process.peak_filters import filter_attr
from process.scan_processing import average_replicate_scans
from process.scan_processing import join_peaklists
from process.scan_processing import read_scans
from process.scan_processing import remove_edges



def process_scans(source, function_noise, snr_thres, nscans, ppm, min_fraction=None, rsd_thres=None, filelist=None, mzrs_to_remove=[], scan_events=[], block_size=2000, ncpus=None):

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

        scans = read_scans(filenames[i], source, function_noise, nscans, scan_events)

        if scan_events == "all":
            pls_filt = average_replicate_scans(scans, snr_thres, ppm, min_fraction, rsd_thres, block_size, ncpus)
            for h in pls_filt:
                pl_filt_copy = join_peaklists(os.path.basename(filenames[i]), {h: pls_filt[h]})
                for k in fl.keys():
                    pl_filt_copy.metadata[k] = fl[k][i]
                pl_filt_copy.metadata["header"] = h
                pl_filt_copy.metadata["scan_ids"] = [int(pl_scan.ID) for pl_scan in scans[h]]
                pls.append(pl_filt_copy)

        elif type(scan_events) is list or scan_events is None:

            scans_er = remove_edges(scans)
            prs = average_replicate_scans(scans_er, snr_thres, ppm, min_fraction, rsd_thres, block_size, ncpus)
            pl_filt = join_peaklists(os.path.basename(filenames[i]), prs)

            if type(mzrs_to_remove) == list:
                if len(mzrs_to_remove) > 0:
                    pl_filt = filter_mz_ranges(pl_filt, mzrs_to_remove)
                else:
                    pass

            elif mzrs_to_remove is not None:
                raise ValueError("mzr_remove: Provide a list of 'start' and 'end' values for each m/z range that needs to be removed.")
            else:
                pass

            if "class" in fl:
                pl_filt.tags.add_tags(class_label=fl["class"][i])

            for k in fl.keys():
                pl_filt.metadata[k] = fl[k][i]
            pls.append(pl_filt)

        else:
            raise ValueError("scan_events: Provide a list e.g. [[50.0, 1000.0, full]] or the string 'all'")

    return pls


# placeholder (synonym)
def stitch(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres=None, rsd_thres=None, block_size=2000, ncpus=None):
    return process_scans(source, filelist, fn_exp, nscans, function_noise, snr_thres, ppm, presence_thres, rsd_thres, block_size, ncpus)


def replicate_filter(source, ppm, replicates, min_peaks, rsd_thres=None, filelist=None, block_size=2000, ncpus=None):

    if replicates < min_peaks:
        raise IOError("Provide realistic values for the number of replicates and minimum number of peaks present (min_peaks)")

    filenames = check_paths(filelist, source)
    if len(filenames) == 0:
        raise IOError("Provide a filelist that list all the text files (columnname:filename) and assign replicate numbers to each filename/sample (columnname:replicate)")
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata(peaklists, fl)

    if not hasattr(peaklists[0].metadata, "replicate"):
        raise IOError("Provide a filelist and assign replicate numbers (columnname:replicate) to each filename/sample")

    idxs_peaklists = idxs_reps_from_filelist([pl.metadata.replicate for pl in peaklists])
    unique, counts = np.unique([pl.metadata.replicate for pl in peaklists], return_counts=True)

    if len(counts) <= 1:
        raise ValueError("No technical replicates available (single) - Skip 'replicate filter'")
    if max(unique) < replicates:
        raise ValueError("Replicates incorrectly labeled")
    if sum(counts) != len(peaklists):
        raise ValueError("Replicates incorrectly labeled")

    pls_rep_filt = []

    if len([True for idxs_pls in idxs_peaklists if len(idxs_pls) > replicates]) > 0:
        print
        print "All combinations (n={}) for each each set of replicates will be processed to calculate the most reproducible set".format(replicates)
        print
        print "rank\tID\tpeaks\tmedian_RSD({}/{})".format(replicates, replicates)

    for idxs_pls in idxs_peaklists:

        temp = []

        for j, pls_comb in enumerate(combinations(peaklists[idxs_pls[0]:idxs_pls[-1] + 1], replicates)):

            pm = align_peaks(pls_comb, ppm, block_size, ncpus=ncpus)

            #############################################################
            # TODO: Write some sort of merging function for multiple replicate file names
            #############################################################
            prefix = os.path.commonprefix([p.ID for p in pls_comb])
            merged_id = "{}{}".format(prefix, "_".join(map(str, [p.ID.replace(prefix, "").split(".")[0] for p in pls_comb])))
            #############################################################

            pl = pm.to_peaklist(ID=merged_id)
            if "snr" in pm.attributes:
                pl.add_attribute("snr", pm.attr_mean_vector("snr"), on_index=2)

            pl.tags.add_tags(*pls_comb[0].tags.tag_of(None), **{t: pls_comb[0].tags.tag_of(t) for t in pls_comb[0].tags.tag_types})
            pl.add_attribute("present_flag", pm.present >= min_peaks, is_flag=True)

            if rsd_thres is not None:
                rsd_flag = map(lambda x: not np.isnan(x) and x < rsd_thres, pm.rsd)
                pl.add_attribute("rsd_flag", rsd_flag, flagged_only=False, is_flag=True)

            for k in pls_comb[0].metadata:
                if k != "filename":
                    pl.metadata[k] = pls_comb[0].metadata[k]

            pl_filt = filter_attr(pl.copy(), attr_name="present", min_threshold=replicates, flag_name="pres_rsd")
            temp.append([pl, pl_filt.shape[0], np.median(pl_filt.rsd)])

        temp.sort(key=operator.itemgetter(2, 1))
        pls_rep_filt.append(temp[0][0]) # Most reproducible set of replicates

        for p in range(0, len(temp)):
            print "{}\t{}\t{}\t{}\t".format(p+1, temp[p][0].ID, temp[p][1], temp[p][2])
        print

    return pls_rep_filt


def align_samples(source, ppm, filelist=None, block_size=2000, ncpus=None):

    filenames = check_paths(filelist, source)
    peaklists = load_peaklists(source)

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

    if not isinstance(peak_matrix, PeakMatrix):
        if h5py.is_hdf5(peak_matrix):
            peak_matrix = hdf5_portal.load_peak_matrix_from_hdf5(peak_matrix)
        else:
            peak_matrix = txt_portal.load_peak_matrix_from_txt(peak_matrix)

    if class_labels is not None:
        peak_matrix = update_class_labels(peak_matrix, class_labels)

    if blank_label not in peak_matrix.peaklist_tag_values:
        raise IOError("Blank label ({}) does not exist".format(blank_label))

    return filter_blank_peaks(peak_matrix, blank_label, min_fraction, min_fold_change, function, rm_samples)


def sample_filter(peak_matrix, min_fraction, within=False, rsd=None, qc_label=None, class_labels=None):

    if not isinstance(peak_matrix, PeakMatrix):
        if h5py.is_hdf5(peak_matrix):
            peak_matrix = hdf5_portal.load_peak_matrix_from_hdf5(peak_matrix)
        else:
            peak_matrix = txt_portal.load_peak_matrix_from_txt(peak_matrix)

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


def hdf5_to_txt(fname, path_out, attr_name="intensity", separator="\t", transpose=False, comprehensive=False):
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
            pk_out.write(obj.to_str(attr_name=attr_name, delimiter=separator, transpose=transpose, comprehensive=comprehensive))
    else:
        assert os.path.isdir(path_out), "File or Directory does not exist:".format(path_out)
        obj = hdf5_portal.load_peaklists_from_hdf5(fname)
        assert isinstance(obj[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
        for pl in obj:
            with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                pk_out.write(pl.to_str(delimiter=separator))
    return


def load_peaklists(source):

    if type(source) == str:
        source = source.encode('string-escape')
        if h5py.is_hdf5(source):
            peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
        elif zipfile.is_zipfile(source):
            zf = zipfile.ZipFile(source)
            filenames = zf.namelist()
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [txt_portal.load_peaklist_from_txt(zf.open(fn), ID=os.path.basename(fn), has_flag_col=True) for fn in filenames]
        elif os.path.isdir(source):
            filenames = os.listdir(source)
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [txt_portal.load_peaklist_from_txt(os.path.join(source, fn), ID=os.path.basename(fn), delimiter = "\t", has_flag_col=False) for fn in filenames]
        else:
            raise TypeError("Incorrect format. Process .mzML and .raw files first using the 'process scans' function")
    elif type(source) == list or type(source) == tuple:
        if not isinstance(source[0], PeakList):
            raise TypeError("List has incorrect format. PeakList objects required.")
        else:
            peaklists = source
    else:
        raise IOError("Inccorrect input: list with peaklist objects or path")

    return peaklists
