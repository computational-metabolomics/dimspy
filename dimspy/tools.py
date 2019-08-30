#!/usr/bin/python
# -*- coding: utf-8 -*-

import collections
import logging
import operator
import os
import zipfile
from itertools import combinations
from typing import Sequence, Dict

import h5py
import numpy as np

from .experiment import check_metadata
from .experiment import idxs_reps_from_filelist
from .experiment import interpret_experiment
from .experiment import mz_range_from_header
from .experiment import update_labels
from .experiment import update_metadata_and_labels
from .models.peak_matrix import PeakMatrix
from .models.peaklist import PeakList
from .models.peaklist_tags import Tag
from .portals import hdf5_portal
from .portals import txt_portal
from .portals.paths import check_paths
from .process.peak_alignment import align_peaks
from .process.peak_filters import filter_attr
from .process.peak_filters import filter_blank_peaks
from .process.peak_filters import filter_fraction
from .process.peak_filters import filter_mz_ranges
from .process.peak_filters import filter_ringing
from .process.peak_filters import filter_rsd
from .process.replicate_processing import average_replicate_peaklists
from .process.replicate_processing import average_replicate_scans
from .process.replicate_processing import join_peaklists
from .process.replicate_processing import read_scans
from .process.replicate_processing import remove_edges


def process_scans(source: str, function_noise: str, snr_thres: float, ppm: float, min_fraction: float or None = None,
                  rsd_thres: float or None = None, min_scans: int = 1, filelist: str or None = None,
                  skip_stitching: bool = False, remove_mz_range: list or None = None,
                  ringing_thres: float or None = None, filter_scan_events: Dict or None = None,
                  report: str or None = None, block_size: int = 5000, ncpus: int or None = None):
    """

    :param source:
    :param function_noise:
    :param snr_thres:
    :param ppm:
    :param min_fraction:
    :param rsd_thres:
    :param min_scans:
    :param filelist:
    :param skip_stitching:
    :param remove_mz_range:
    :param ringing_thres:
    :param filter_scan_events:
    :param report:
    :param block_size:
    :param ncpus:
    :return:
    """

    if filter_scan_events is None:
        filter_scan_events = {}
    if remove_mz_range is None:
        remove_mz_range = []

    filenames = check_paths(filelist, source)

    if len([fn for fn in filenames if not fn.lower().endswith(".mzml") or not fn.lower().endswith(".raw")]) == 0:
        raise IOError("Incorrect file format. Provide .mzML and .raw files")

    if filelist is not None:
        fl = check_metadata(filelist)
    else:
        fl = collections.OrderedDict()

    if report is not None:
        out = open(report, "w")
        out.write("filename\tevent\tscans\tpeaks\tmedian_rsd\n")

    pls = []
    for i in range(len(filenames)):

        print()
        print((os.path.basename(filenames[i])))

        if type(source) is not str:
            source = ""

        print("Reading scans....")
        pls_scans = read_scans(filenames[i], source, function_noise, min_scans, filter_scan_events)

        if type(remove_mz_range) == list and len(remove_mz_range) > 0:
            print("Removing m/z ranges.....")
            for h in pls_scans:
                pls_scans[h] = [filter_mz_ranges(pl, remove_mz_range) if len(pl.mz) > 0 else pl
                                for pl in pls_scans[h]]

        if not skip_stitching:
            mz_ranges = [mz_range_from_header(h) for h in pls_scans]
            exp = interpret_experiment(mz_ranges)
            if exp == "overlapping":
                print("Removing 'edges' from SIM windows.....")
                pls_scans = remove_edges(pls_scans)

        if ringing_thres is not None and float(ringing_thres) > 0.0:
            print("Removing ringing artifacts.....")
            for h in pls_scans:
                pls_scans[h] = [filter_ringing(pl, threshold=ringing_thres, bin_size=1.0) if len(pl.mz) > 0 else pl
                                for pl in pls_scans[h]]

        print("Removing noise.....")
        for h in pls_scans:
            pls_scans[h] = [filter_attr(pl, "snr", min_threshold=snr_thres) if len(pl.mz) > 0 else pl
                            for pl in pls_scans[h]]

        print("Aligning, averaging and filtering peaks.....")
        pls_avg = []

        for h in pls_scans:

            nscans, n_peaks, median_rsd = len(pls_scans[h]), 0, "NA"
            # pls_scans[h] = [pl for pl in pls_scans[h] if len(pl.mz) > 0]

            if len(pls_scans[h]) >= 1:
                if sum(pl.shape[0] for pl in pls_scans[h]) == 0:
                    logging.warning("No scan data available for {}".format(h))
                else:
                    pl_avg = average_replicate_scans(h, pls_scans[h], ppm, min_fraction, rsd_thres, "intensity",
                                                     block_size, ncpus)
                    pls_avg.append(pl_avg)
                    n_peaks, median_rsd = pl_avg.shape[0], np.nanmedian(pl_avg.rsd)
            else:
                logging.warning("No scan data available for {}".format(h))

            if report is not None:
                out.write("{}\t{}\t{}\t{}\t{}\n".format(os.path.basename(filenames[i]), h, nscans, n_peaks, median_rsd))

        if len(pls_avg) == 0:
            raise IOError("No peaks remaining after filtering. Remove file from Study (filelist).")

        if not skip_stitching or len(list(pls_scans.keys())) == 1:
            pl = join_peaklists(os.path.basename(filenames[i]), pls_avg)
            pl = update_metadata_and_labels([pl], fl)
            pls.extend(pl)
            if len(list(pls_scans.keys())) > 1 and report is not None:
                out.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(os.path.basename(filenames[i]), "SIM-Stitch", "NA", pl[0].shape[0],
                                                  np.nanmedian(pl[0].rsd)))
        else:
            for pl in pls_avg:
                pl = update_metadata_and_labels([pl], fl)
                pl = join_peaklists("{}#{}".format(os.path.basename(filenames[i]), pl[0].metadata["header"][0]), pl)
                pls.append(pl)

    if report is not None:
        out.close()

    return pls


# placeholder (synonym)
def sim_stitch(source: str, function_noise: str, snr_thres: float, ppm: float, min_fraction: float or None = None,
               rsd_thres: float or None = None, min_scans: int = 1, filelist: str or None = None,
               skip_stitching: bool = False, remove_mz_range: list or None = None, ringing_thres: float or None = None,
               filter_scan_events: Dict or None = None, report: str or None = None, block_size: int = 5000,
               ncpus: int or None = None):
    """

    :param source:
    :param function_noise:
    :param snr_thres:
    :param ppm:
    :param min_fraction:
    :param rsd_thres:
    :param min_scans:
    :param filelist:
    :param skip_stitching:
    :param remove_mz_range:
    :param ringing_thres:
    :param filter_scan_events:
    :param report:
    :param block_size:
    :param ncpus:
    :return:
    """
    if filter_scan_events is None:
        filter_scan_events = {}
    if remove_mz_range is None:
        remove_mz_range = []

    return process_scans(source, function_noise, snr_thres, ppm, min_fraction, rsd_thres, min_scans, filelist,
                         skip_stitching, remove_mz_range, ringing_thres, filter_scan_events, report, block_size, ncpus)


def replicate_filter(source: str or Sequence[PeakList], ppm: float, replicates: int, min_peaks: int,
                     rsd_thres: float or None = None, filelist: str or None = None, report: str or None = None,
                     block_size: int = 5000, ncpus: int or None = None):
    """

    :param source:
    :param ppm:
    :param replicates:
    :param min_peaks:
    :param rsd_thres:
    :param filelist:
    :param report:
    :param block_size:
    :param ncpus:
    :return:
    """

    if replicates < min_peaks:
        raise IOError(
            "Provide realistic values for the number of replicates and minimum number of peaks present (min_peaks)")

    filenames = check_paths(filelist, source)
    if len(filenames) == 0:
        raise IOError(
            "Provide a filelist that list all the text files (columnname:filename) and assign replicate numbers to "
            "each filename/sample (columnname:replicate)")
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata_and_labels(peaklists, fl)

    if not hasattr(peaklists[0].metadata, "replicate"):
        raise IOError("Provide a filelist and assign replicate numbers (columnname:replicate) to each filename/sample")

    if report is not None:
        out = open(report, "w")

    idxs_peaklists = idxs_reps_from_filelist([pl.metadata.replicate for pl in peaklists])
    unique, counts = np.unique([int(pl.metadata.replicate) for pl in peaklists], return_counts=True)

    if len(counts) <= 1:
        raise ValueError("No technical replicates available (single) - Skip 'replicate filter'")
    if max(unique) < replicates:
        raise ValueError("Replicates incorrectly labeled")
    if sum(counts) != len(peaklists):
        raise ValueError("Replicates incorrectly labeled")

    reps_each_sample = [len(idxs_pls) for idxs_pls in idxs_peaklists]
    if min(reps_each_sample) < replicates:
        raise ValueError("Not enough (technical) replicates available for each sample.")

    if max(reps_each_sample) > replicates:
        print(("NOTE: All combinations (n={}) for each each set of replicates are "
               "processed to calculate the most reproducible set.".format(replicates)))
        if report is not None:
            out.write("set\trank\tname\tpeaks\tpeaks_{}oo{}\tmedian_rsd_{}oo{}\tscore\n".format(replicates, replicates,
                                                                                                replicates, replicates))
    else:
        if report is not None:
            out.write(
                "name\tpeaks\tpeaks_{}oo{}\tmedian_rsd_{}oo{}\n".format(replicates, replicates, replicates, replicates))

    pls_rep_filt = []
    for idxs_pls in range(len(idxs_peaklists)):

        temp = []
        max_peak_count, max_peak_count_present = 0, 0

        for pls_comb in combinations(peaklists[idxs_peaklists[idxs_pls][0]:idxs_peaklists[idxs_pls][-1] + 1],
                                     replicates):

            pl = average_replicate_peaklists(pls_comb, ppm, min_peaks, rsd_thres, block_size=block_size, ncpus=ncpus)

            if hasattr(pls_comb[0].metadata, "injectionOrder"):
                pl.metadata["injectionOrder"] = int(pls_comb[0].metadata["injectionOrder"])
                pl.tags.add_tag(int(pls_comb[0].metadata["injectionOrder"]), "injectionOrder")

            reps = [_pl.metadata["replicate"] for _pl in pls_comb]
            pl.metadata["replicates"] = reps
            pl.tags.add_tag("-".join(map(str, reps)), "replicates")

            for t in pls_comb[0].tags.tags:
                if t.ttype != "replicate":
                    if not pl.tags.has_tag_type(t.ttype):
                        pl.tags.add_tag(t)

            pl_filt = filter_attr(pl.copy(), attr_name="present", min_threshold=replicates, flag_name="pres_rsd")
            median_rsd = np.median(pl_filt.get_attribute("rsd", flagged_only=True))

            temp.append([pl, pl.shape[0], pl_filt.shape[0], median_rsd])

            if pl_filt.shape[0] > max_peak_count_present:
                max_peak_count_present = pl_filt.shape[0]

        # find the RSD category for the median RDS
        bins = np.array(list(range(0, 55, 5)))
        rsd_scores = [0.1 * b for b in reversed(list(range(len(bins))))]
        for i, comb in enumerate(temp):

            if np.isnan(comb[3]):
                rsd_score = 0.0
            else:
                inds = np.digitize([comb[3]], bins)
                rsd_score = rsd_scores[inds[0] - 1]

            # score 1: peak count / peak count present in n-out-n (e.g. 3-out-of-3)
            # score 2: peak count present in n-out-n (e.g. 3-out-of-3) / MAX peak count present in n-out-n across replicates
            # score 3: RSD categories (0-5 (score=1.0), 5-10 (score=0.9), 10-15 (score=0.8), etc)
            scores = [comb[2] / float(comb[1]), comb[2] / float(max_peak_count_present), rsd_score]
            if np.isnan(sum(scores)):
                scores.append(0)
            else:
                scores.append(sum(scores) / 3.0)
            temp[i].extend(scores)

        if sum([comb[-1] for comb in temp]) == 0.0:
            logging.warning(
                "insufficient data available to calculate scores for {}".format(str([comb[0].ID for comb in temp])))

        # sort the scores from high to low
        temp.sort(key=operator.itemgetter(-1), reverse=True)
        # select the replicate filtered peaklist that is ranked first
        pls_rep_filt.append(temp[0][0])

        if report is not None:
            for p in range(0, len(temp)):
                if max(reps_each_sample) > replicates:
                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(idxs_pls + 1, p + 1, temp[p][0].ID, temp[p][1],
                                                                    temp[p][2], temp[p][3], temp[p][-1]))
                else:
                    out.write("{}\t{}\t{}\t{}\n".format(temp[p][0].ID, temp[p][1], temp[p][2], temp[p][3]))

    if report is not None:
        out.close()

    return pls_rep_filt


def align_samples(source: str or Sequence[PeakList], ppm: float, filelist: str or None = None, block_size: int = 5000,
                  ncpus: int or None = None):
    """

    :param source:
    :param ppm:
    :param filelist:
    :param block_size:
    :param ncpus:
    :return:
    """

    filenames = check_paths(filelist, source)
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = check_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata_and_labels(peaklists, fl)

    return align_peaks(peaklists, ppm=ppm, block_size=block_size, ncpus=ncpus)


def blank_filter(peak_matrix: str or PeakMatrix, blank_label: str, min_fraction: float = 1.0,
                 min_fold_change: float = 1.0, function: str = "mean", rm_samples: bool = True,
                 labels: str or None = None):
    """

    :param peak_matrix:
    :param blank_label:
    :param min_fraction:
    :param min_fold_change:
    :param function:
    :param rm_samples:
    :param labels:
    :return:
    """

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

    if labels is not None:
        peak_matrix = update_labels(peak_matrix, labels)

    if not any([Tag(blank_label, 'classLabel') in x for x in peak_matrix.peaklist_tags]):
        raise IOError("Blank label ({}) does not exist".format(blank_label))

    return filter_blank_peaks(peak_matrix, Tag(blank_label, 'classLabel'), min_fraction, min_fold_change, function,
                              rm_samples)


def sample_filter(peak_matrix: str or PeakMatrix, min_fraction: float, within: bool = False, rsd: float or None = None,
                  qc_label: str or None = None, labels: str or None = None):
    """

    :param peak_matrix:
    :param min_fraction:
    :param within:
    :param rsd:
    :param qc_label:
    :param labels:
    :return:
    """

    if not isinstance(peak_matrix, PeakMatrix):
        if h5py.is_hdf5(peak_matrix):
            peak_matrix = hdf5_portal.load_peak_matrix_from_hdf5(peak_matrix)
        else:
            peak_matrix = txt_portal.load_peak_matrix_from_txt(peak_matrix)

    if labels is not None:
        if not os.path.isfile(labels):
            raise IOError("{} does not exist".format(labels))
        peak_matrix = update_labels(peak_matrix, labels)

    if qc_label is not None:
        if Tag(qc_label, 'classLabel') not in peak_matrix.peaklist_tags:
            raise IOError("QC label ({}) does not exist".format(qc_label))

    peak_matrix = filter_fraction(peak_matrix, min_fraction, within_classes=within, class_tag_type="classLabel")

    if rsd is not None:
        peak_matrix = filter_rsd(peak_matrix, rsd, Tag(qc_label, "classLabel"))
    return peak_matrix


def missing_values_sample_filter(peak_matrix: PeakMatrix, max_fraction: float):
    """

    :param peak_matrix:
    :param max_fraction:
    :return:
    """

    return peak_matrix.remove_samples(
        np.where([(x / float(peak_matrix.shape[1]) >= max_fraction) for x in peak_matrix.missing_values]))


def remove_samples(obj: PeakList or PeakMatrix, sample_names: list):
    """

    :param obj:
    :param sample_names:
    :return:
    """

    if isinstance(obj, PeakMatrix):
        return obj.remove_samples(np.where([(x in sample_names) for x in obj.peaklist_ids]))
    elif isinstance(obj[0], PeakList):
        return [pl for pl in obj if pl.ID not in sample_names]
    else:
        raise IOError("Incorrect format - PeakMatrix object or list of PeakList objects")


def hdf5_peak_matrix_to_txt(filename: str, path_out: str, attr_name: str = "intensity", rsd_tags: tuple = (),
                            delimiter: str = "\t", samples_in_rows: bool = True, comprehensive: bool = False,
                            compatibility_mode: bool = False):
    """

    :param filename:
    :param path_out:
    :param attr_name:
    :param rsd_tags:
    :param delimiter:
    :param samples_in_rows:
    :param comprehensive:
    :param compatibility_mode:
    :return:
    """

    if not os.path.isfile(filename):
        raise IOError('HDF5 database [%s] does not exist' % filename)
    if not h5py.is_hdf5(filename):
        raise IOError('input file [%s] is not a valid HDF5 database' % filename)

    obj = hdf5_portal.load_peak_matrix_from_hdf5(filename, compatibility_mode=compatibility_mode)
    with open(os.path.join(path_out), "w") as pk_out:
        pk_out.write(obj.to_str(attr_name=attr_name, delimiter=delimiter,
                                samples_in_rows=samples_in_rows, rsd_tags=rsd_tags,
                                comprehensive=comprehensive))
    return


def hdf5_peaklists_to_txt(filename: str, path_out: str, delimiter: str = "\t", compatibility_mode: bool = False):
    """

    :param filename:
    :param path_out:
    :param delimiter:
    :param compatibility_mode:
    :return:
    """

    if not os.path.isfile(filename):
        raise IOError('HDF5 database [%s] does not exist' % filename)
    if not h5py.is_hdf5(filename):
        raise IOError('input file [%s] is not a valid HDF5 database' % filename)

    if not os.path.isdir(path_out):
        raise IOError("File or Directory does not exist:".format(path_out))

    obj = hdf5_portal.load_peaklists_from_hdf5(filename, compatibility_mode=compatibility_mode)
    if "#" in obj[0].ID:
        fns = set([pl.ID.split("#")[0] for pl in obj])
        sub_ids = [pl.ID.split("#")[1] for pl in obj]
        for fn in fns:
            with open(os.path.join(path_out, os.path.splitext(fn)[0] + ".txt"), "w") as pk_out:
                for i, pl in enumerate(obj):
                    if fn in pl.ID:
                        pl.add_attribute("window", pl.full_shape[0] * [sub_ids[i]], flagged_only=False, on_index=3)
                        pk_out.write(pl.to_str(delimiter=delimiter))
                        pl.drop_attribute("window")
    else:
        for pl in obj:
            with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                pk_out.write(pl.to_str(delimiter=delimiter))
    return


def merge_peaklists(source: Sequence[PeakList], filelist: str or None = None):
    """

    :param source:
    :param filelist:
    :return:
    """

    if not isinstance(source, list):
        raise IOError(
            "Incorrect input: list of lists of peaklists, list of peak matrix objects or list of HDF5 files expected.")

    pls_merged = []
    for s in source:
        if isinstance(s, list) or isinstance(s, tuple):
            if isinstance(s[0], PeakList):
                pls_merged.extend(s)
            else:
                raise IOError("Incorrect Object in list. Peaklist Object expected.")
        elif isinstance(s, PeakMatrix):
            pls = s.extract_peaklists()
            pls_merged.extend(pls)
        elif h5py.is_hdf5(s):
            f = h5py.File(s, 'r')
            if "mz" in f:
                pm = txt_portal.load_peak_matrix_from_txt(s)
                pls = pm.extract_peaklists()
            else:
                pls = hdf5_portal.load_peaklists_from_hdf5(s)
            f.close()
            pls_merged.extend(pls)
        else:
            raise IOError(
                "Incorrect input: list of lists of peaklists, list of peak matrix objects or list of HDF5 files expected.")

    if filelist is not None:
        fl = check_metadata(filelist)
        pls_merged = update_metadata_and_labels(pls_merged, fl)

        if 'multilist' in list(fl.keys()):
            # make sure the peaklists are in the correct order (need to be sorted ascending)
            order_indx = np.argsort([i.metadata['multilist'] for i in pls_merged])
            nlists = [fl['multilist'][i] for i in order_indx]
            pls_merged = [pls_merged[i] for i in order_indx]

            # get the break points of the different lists to join together
            bp = list(np.cumsum(np.unique(nlists, return_counts=True)[1]))
            bp = bp[:-1]

            # break up the list into a list of lists
            pls_merged = partition(pls_merged, bp)

    return pls_merged


def partition(alist: list, indices: list):
    """

    :param alist:
    :param indices:
    :return:
    """

    return [alist[i:j] for i, j in zip([0] + indices, indices + [None])]


def load_peaklists(source: str or Sequence[PeakList]):
    """

    :param source:
    :return:
    """

    if type(source) == str:
        source = source.encode('string-escape')
        if h5py.is_hdf5(source):
            peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
        elif zipfile.is_zipfile(source):
            zf = zipfile.ZipFile(source)
            filenames = zf.namelist()
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0, \
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [txt_portal.load_peaklist_from_txt(zf.open(fn), ID=os.path.basename(fn), has_flag_col=True) for
                         fn in filenames]
        elif os.path.isdir(source):
            filenames = os.listdir(source)
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0, \
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [
                txt_portal.load_peaklist_from_txt(os.path.join(source, fn), ID=os.path.basename(fn), delimiter="\t",
                                                  has_flag_col=False) for fn in filenames]
        else:
            raise IOError("Incorrect format. Process .mzML and .raw files first using the 'process scans' function")
    elif type(source) == list or type(source) == tuple:
        if not isinstance(source[0], PeakList):
            raise IOError("List has incorrect format. PeakList objects required.")
        else:
            peaklists = source
    else:
        raise IOError("Inccorrect input: list with peaklist objects or path")

    return peaklists


def create_sample_list(source, path_out, delimiter="\t", qc_label="QC"):
    """

    :param source:
    :param path_out:
    :param delimiter:
    :param qc_label:
    :return:
    """

    if isinstance(source, list) or isinstance(source, tuple):

        if not isinstance(source[0], PeakList):
            raise IOError("Inccorrect input: list with peaklist objects")

        pls_tags = [pl.tags for pl in source]
        pls_ids = [pl.ID for pl in source]

    elif isinstance(source, PeakMatrix):

        pls_tags = source.peaklist_tags
        pls_ids = source.peaklist_ids

    else:
        raise IOError("Incorrect format")

    header = ["filename"]
    header.extend([tn.ttype for tn in pls_tags[0].typed_tags])
    with open(path_out, "w") as out:
        out.write("{}".format(delimiter).join(map(str, header)) + "\n")
        for i in range(len(pls_tags)):
            row = ["NA"] * len(header)
            row[0] = pls_ids[i]
            for t in pls_tags[i].typed_tags:
                row[header.index(t.ttype)] = pls_tags[i].tag_of(t.ttype).value
            out.write("{}".format(delimiter).join(map(str, row)) + "\n")

    return
