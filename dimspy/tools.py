#!/usr/bin/python
# -*- coding: utf-8 -*-

import collections
import logging
import operator
import os
from itertools import combinations
from typing import Sequence, Dict, Union

import h5py
import numpy as np

from .metadata import validate_metadata
from .metadata import idxs_reps_from_filelist
from .metadata import interpret_method
from .metadata import mz_range_from_header
from .metadata import update_labels
from .metadata import update_metadata_and_labels
from .models.peak_matrix import PeakMatrix
from .models.peaklist import PeakList
from .models.peaklist_tags import Tag
from .portals import hdf5_portal
from .portals import txt_portal
from .portals.paths import validate_and_sort_paths
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


def process_scans(source: str, function_noise: str, snr_thres: float, ppm: float, min_fraction: Union[float, None] = None,
                  rsd_thres: Union[float, None] = None, min_scans: int = 1, filelist: Union[str, None] = None,
                  skip_stitching: bool = False, remove_mz_range: list or None = None,
                  ringing_thres: Union[float, None] = None, filter_scan_events: Dict or None = None,
                  report: Union[str, None] = None, block_size: int = 5000, ncpus: int or None = None):
    """
    Extract, filter and average spectral data from input .RAW or .mzML files and generate a single mass
    spectral peaklist (object) for each of the data files defined in the ‘filelist’ (see below).

    .. warning::
        When using .mzML files generated using the Proteowizard tool, SIM-type scans will only be treated
        as spectra if the ‘simAsSpectra’ filter was set to true during the conversion process:
        *msconvert.exe example.raw* **--simAsSpectra** *--64 --zlib --filter "peakPicking true 1-”*

    :param source: Path to a set of input .RAW or .mzML files
    :param function_noise: Function to calculate the noise from each scan. The following options are available:

        * **median** - the median of all peak intensities within a given scan is used as the noise value.

        * **mean** - the unweighted mean average of all peak intensities within a given scan is used as the noise value.

        * **mad (Mean Absolute Deviation)** - the noise value is set as the mean of the absolute differences between peak
          intensities and the mean peak intensity (calculated across all peak intensities within a given scan).

        * **noise_packets** - the noise value is calculated using the proprietary algorithms contained in Thermo Fisher
          Scientific’s msFileReader library. This option should only be applied when you are processing .RAW files.

    :param snr_thres: Peaks with a signal-to-noise ratio (SNR) less-than or equal-to this value will be removed
        from the output peaklist.

    :param ppm: A positive numerical value equal-to or greater-than zero. This option impacts the clustering of peaks
        extracted from an input file. If the mass-to-charge ratios of two peaks, when divided by the average of
        their mass-to-charge ratios and then multiplied by 1 × 10\ :sup:`6`, is equal-to or less-than this user-defined value,
        then these peaks are clustered together as a single peak. Clustering is applied across all replicates of a given
        scan event type i.e. with a given input file, all peaks detected in the three replicates of a 50-400 *m/z* scan event
        would undergo assessment for the need for clustering.

    :param min_fraction: A numerical value from 0 to 1 that specifies the minimum proportion of scans a given mass
        spectral peak must be detected in, in order for it to be kept in the output peaklist. Here, scans refers to
        replicates of the same scan event type, i.e. if set to 0.33, then a peak would need to be detected in at least
        1 of the 3 replicates of a given scan event type.

    :param rsd_thres: Relative standard deviation threshold - A numerical value equal-to or greater-than 0.
        If greater than 0, then peaks whose intensity values have a percent relative standard deviation (otherwise termed
        the percent coefficient of variation) greater-than this value are excluded from the output peaklist.

    :param min_scans: Minimum number of scans required for each *m/z* window or event within a raw/mzML data file.

    :param filelist: A tab-delimited text file containing **filename** and **classLabel** information for each
        experimental sample. These column headers MUST be included in the first row of the table. For a standard DIMS
        experiment, users are advised to also include the following additional columns:

        * **injectionOrder** - integer values ranging from 1 to i, where i is the total number of independent injections
          performed as part of a DIMS experiment. e.g. if a study included 20 samples, each of which was injected as four
          independent replicates, there would be at least 20 * 4 injections, so i = 80 and the range for injection
          order would be from 1 to 80 in steps of 1.

        * **replicate** - integer value from 1 to r, indicating the order in which technical replicates of each study
          sample were injected in to the mass spectrometer, e.g. if study samples were analysed in quadruplicate,
          r = 4 and integer values are accordingly 1, 2, 3, 4.

        * **batch** - integer value from 1 to b, where b corresponds to the total number of batches analysed under
          define analysis conditions, for any given experiment. e.g. : if 4 independent plates of polar extracts were
          analysed in the positive ionisation mode, then valid values for batch are 1, 2, 3 and 4.

        This filelist may include additional columns, e.g. additional metadata relating to study samples. Ensure that columns
        names do not conflict with existing column names.

    :param skip_stitching: Selected Ion Monitoring (SIM) scans with overlapping scan ranges can be "stitched" together
        in to a pseudo-spectrum. This is achieved by setting this parameter to False (default).
    :param remove_mz_range: this option allows for specific m/z regions of the output peaklist to be deleted, this
        option may be useful for removing sections of a spectrum known to correspond to system noise peaks.

    :param ringing_thres: Fourier transform-based mass spectra often contain peaks (ringing artefacts) around
        spectral features that require removal. This threshold is a positive float indicating the required relative
        intensity a peak must exceed (with reference to the largest peak in a cluster of peaks) in order to be retained.

    :param filter_scan_events: Include or exclude specific scan events, by default all ALL scan events will be
        included. To include or exclude specific scan events use the following format of a dictionary.

        >>> {"include":[[100, 300, "sim"]]} or {"include":[[100, 1000, "full"]]}

    :param report: A tab-delimited text file to write measures of quality (e.g. RSD, number of peaks, etc) for each scan event processed in each .RAW or .mzML files.
    :param block_size: Number peaks in each centre clustering block.
    :param ncpus: Number of CPUs for parallel clustering. Default = None, indicating using all CPUs that are available

    :return: List of peaklist objects
    """

    if filter_scan_events is None:
        filter_scan_events = {}
    if remove_mz_range is None:
        remove_mz_range = []

    filenames = validate_and_sort_paths(source, filelist)

    if len([fn for fn in filenames if not fn.lower().endswith(".mzml") or not fn.lower().endswith(".raw")]) == 0:
        raise IOError("Incorrect file format. Provide .mzML and .raw files")

    if filelist is not None:
        fl = validate_metadata(filelist)
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
            exp = interpret_method(mz_ranges)
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


def replicate_filter(source: Union[Sequence[PeakList], str], ppm: float, replicates: int, min_peaks: int,
                     rsd_thres: Union[float, None] = None, filelist: Union[str, None] = None, report: Union[str, None] = None,
                     block_size: int = 5000, ncpus: int or None = None):
    """
    Peaks from each technical replicate (for a given study sample) are aligned using a one-dimensional hierarchical
    clustering procedure (applied on the mass-to-charge level).
    Peaks are aligned only if the difference in their mass-to-charge ratios, when divided by the average of their
    mass-to-charge ratios and multiplied by 1 × 10\ :sup:`6` \ (i.e. when measured in units of parts-per-million, ppm),
    is less-than or equal-to the user-defined ‘ppm error tolerance’. After alignment, a set of user-defined filters are
    applied to retain only those peaks that:

        * occur in equal-to or more-than the user-defined 'Number of technical replicates a peak has to be present
          in', i.e. if set to 2, then a peak must be detected in at least two of the replicate analyses, **and/or**

        * have relative standard deviation (measured in %; may otherwise be referred to as the percent coefficient
          of variation) of intensity values, across technical replicates, that is equal-to or less-than the user-defined
          ‘relative standard deviation threshold’ (if defined, otherwise ignored).

    .. warning::
        When the parameter “number of technical replicates for each sample” is set to a value less-than the total
        number of technical replicates actually acquired for each study sample, this tool will automatically determine
        which combination of technical replicates to combine. See the parameter description (below) for further details.

    :param source: A list of processed peaklist objects generated by 'process_scans' or path to hdf5 file
    :param ppm:
    :param replicates: Number of technical replicates for each sample - the total number of technical replicates
        acquired for each study sample. This value must be set to the lowest number of technical replicates acquired
        for ANY of the study samples, or alternatively, may be set to the minimum number of replicates the user would
        like to select from the total number of technical replicates for a biological sample.

    :param min_peaks: Minimum number of technical replicates a peak has to be present in.
        For a given biological sample, the number of replicates that will be used to generate the replicate-filtered
        peaklist. If this parameter is set to a value less-than the total number of technical replicates acquired for
        each biological sample, it will automatically determines which combination of technical replicates yields
        the best overall rank. Otherwise, all technical replicates are used. Ranking of the combinations of
        technical replicates is based on the average of the following three scores:

        * score 1: peak count / peak count present in n-out-n (e.g. 3-out-of-3)

        * score 2: peak count present in x-out-of-n (e.g. 3-out-of-3) / MAX peak count present in x-out-of-n across
          sets of replicates

        * score 3: RSD categories (0-5 (score=1.0), 5-10 (score=0.9), 10-15 (score=0.8), etc)

    :param rsd_thres: Relative standard deviation threshold - a numerical value from 0 upwards that defines the
        acceptable percentage relative standard deviation (otherwise termed the percent coefficient of variation)
        of a peak’s intensity across technical replicates. Peaks are removed from the output ‘replicate-filtered’
        peaklist if this condition is not met. Set to None to skipe this filter.

    :param filelist:
    :param report: A tab-delimited text file to write measures of quality (e.g. RSD, number of peaks, etc) for each
        processed 'replicate-filtered' peaklist.
    :param block_size: Number peaks in each centre clustering block.
    :param ncpus: Number of CPUs for parallel clustering. Default = None, indicating using all CPUs that are available

    :return: List of peaklist objects
    """

    if replicates < min_peaks:
        raise IOError(
            "Provide realistic values for the number of replicates and minimum number of peaks present (min_peaks)")

    filenames = validate_and_sort_paths(source, filelist)
    if len(filenames) == 0:
        raise IOError(
            "Provide a filelist that list all the text files (columnname:filename) and assign replicate numbers to "
            "each filename/sample (columnname:replicate)")
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = validate_metadata(filelist)
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


def align_samples(source: Union[Sequence[PeakList], str], ppm: float, filelist: Union[str, None] = None, block_size: int = 5000,
                  ncpus: int or None = None):
    """

    :param source: A list of processed peaklist objects generated by 'process_scans' and/or 'replicate_filter', or path to hdf5 file.
    :param ppm:
    :param filelist:
    :param block_size: Number peaks in each centre clustering block.
    :param ncpus: Number of CPUs for parallel clustering. Default = None, indicating using all CPUs that are available

    :return:
    """

    filenames = validate_and_sort_paths(source, filelist)
    peaklists = load_peaklists(source)

    if filelist is not None:
        fl = validate_metadata(filelist)
        peaklists = [pl for pl in peaklists if pl.ID in [os.path.basename(fn) for fn in filenames]]
        peaklists = update_metadata_and_labels(peaklists, fl)

    return align_peaks(peaklists, ppm=ppm, block_size=block_size, ncpus=ncpus)


def blank_filter(peak_matrix: Union[PeakMatrix, str], blank_label: str, min_fraction: float = 1.0,
                 min_fold_change: float = 1.0, function: str = "mean", rm_samples: bool = True,
                 labels: Union[str, None] = None):
    """

    :param peak_matrix: PeakMatrix object
    :param blank_label: Label for the blank samples - a string indicating the name of the class to be used for
        filtering (e.g. blank), i.e. the “reference” class. This string must have been included in the “classLabel”
        column of the metadata file associated with the process_sans or replicate_filter function(s).
    :param min_fraction: A numeric value ranging from 0 to 1. Setting this value to None or 0 will skip this
        filtering step. A value greater than 0 requires that for each peak in the peak intensity matrix,
        at least this proportion of non-reference samples have to have an intensity value that exceeds the product
        of: (A) the average intensity of “reference” class intensities and (B) the user-defined “min_fold_change”.
        If this condition is not met, the peak is removed from the peak intensity matrix.

    :param min_fold_change: A numeric value from 0 upwards. When minimum fraction filtering is enabled, this value
        defines the minimum required ratio between the intensity of a peak in a “non-reference” sample and the average
        intensity of the “reference” sample(s). Peaks with ratios exceeding this threshold are considered to have been
        reliably detected in a “non-reference” sample.

    :param function: Function to calculate the 'reference' intensity

        * **mean** - corresponds to using the non-weighted average of “reference” sample peak intensities
          (NA values are ignored) in calculating the “reference” to “non-reference” peak intensity ratio.

        * **median** - corresponds to using the median of “reference” sample peak intensities (NA values are ignored)
          in calculating the “reference” to “non-reference” peak intensity ratio.

        * **max** corresponds to the use of the maximum intensity among “reference” sample peak intensities
          (NA values are ignored) in calculating the “reference” to “non-reference” peak intensity ratio.

    :param rm_samples: 	Remove blank samples from the output peak matrix:
        * **True** - samples belonging to the user-defined “reference” class are removed from the output peak matrix
        * **False** - samples belonging to the user-defined “reference” class are retained in the output peak matrix.

    :param labels: Path to the metadata file
    :return: PeakMatrix object
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


def sample_filter(peak_matrix: Union[PeakMatrix, str], min_fraction: float, within: bool = False,
                  rsd_thres: Union[float, None] = None, qc_label: Union[str, None] = None, labels: Union[str, None] = None):
    """

    :param peak_matrix: PeakMatrix object or path to hdf5 file
    :param min_fraction: Minimum fraction - a numeric value between 0 and 1 indicating the proportion of study
        samples in which a peak must have a recorded intensity value in order for it to be retained in the output peak
        intensity matrix; e.g. 0.5 means that at least 50% of samples (whether assessed across all classes, or within
        each class individually) must have a recorded intensity value for a specific peak in order for it to be retained
        in the output peak matrix.
    :param within: Apply sample filter within each sample class

        * **False** - check across ALL classes simultaneously whether greater-than the user-defined “Minimum fraction”
          of samples contained an intensity value for a specific mass spectral peak.
        * **True** - check within EACH class separately whether greater-than the user-defined “Minimum fraction” of
          samples contained an intensity value for a specific mass spectral peak.

        .. warning::
            if in ANY class a peak is detected in greater-than the user-defined minimum fraction of samples, then
            the peak is retained in the output peak matrix. For classes in which this condition is not met, the
            peak intensity recorded for that peak (if any) will still be presented in the output peak matrix.
            If no peak intensity was recorded in a sample, then a ‘0’ is inserted in to the peak matrix.

    :param rsd_thres: Relative standard deviation threshold - A numerical value equal-to or greater-than 0.
        If greater than 0, then peaks whose intensity values have a percent relative standard deviation (otherwise termed
        the percent coefficient of variation) greater-than this value are excluded from the output PeakMatrix object.
    :param qc_label: Label for the QC samples - a string indicating the name of the class to be used for
        filtering, i.e. the “reference” class. This string must have been included in the “classLabel”
        column of the metadata file associated with the process_sans or replicate_filter function(s).
    :param labels: Path to a metadata file
    :return: PeakMatrix object
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

    if rsd_thres is not None:
        peak_matrix = filter_rsd(peak_matrix, rsd_thres, Tag(qc_label, "classLabel"))
    return peak_matrix


def missing_values_sample_filter(peak_matrix: PeakMatrix, max_fraction: float):
    """

    :param peak_matrix: PeakMatrix object
    :param max_fraction:
    :return: PeakMatrix object
    """

    return peak_matrix.remove_samples(
        np.where([(x / float(peak_matrix.shape[1]) >= max_fraction) for x in peak_matrix.missing_values]))


def remove_samples(obj: Union[PeakMatrix, Sequence[PeakList]], sample_names: list):
    """

    :param obj: Peaklist Object
    :param sample_names: List of sample names
    :return: List of Peaklist Objects or PeakMatrix object
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


def merge_peaklists(source: Sequence[PeakList], filelist: Union[str, None] = None):
    """

    :param source: List or typle of Peaklist objects, or hdf5 file
    :param filelist:
    :return: Nested lists of Peaklist objects (e.g. [[pl_01, pl_02], [pl_03, pl_04, pl05]]
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
                pm = hdf5_portal.load_peak_matrix_from_hdf5(s)
                pls = pm.extract_peaklists()
            else:
                pls = hdf5_portal.load_peaklists_from_hdf5(s)
            f.close()
            pls_merged.extend(pls)
        else:
            raise IOError(
                "Incorrect input: list of lists of peaklists, list of peak matrix objects or list of HDF5 files expected.")

    if filelist is not None:
        fl = validate_metadata(filelist)
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

    :param alist: list
    :param indices: indices
    :return: Nested list
    """

    return [alist[i:j] for i, j in zip([0] + indices, indices + [None])]


def load_peaklists(source: Union[Sequence[PeakList] or str]):
    """

    :param source: hdf5 file or list of Peaklist objects
    :return:
    :rtype: Sequence[PeakList]
    """

    if type(source) == str:
        source = source.encode('unicode_escape')
        if h5py.is_hdf5(source):
            peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
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


def create_sample_list(source: Union[Sequence[PeakList], PeakMatrix], path_out: str, delimiter: str = "\t"):
    """

    :param source: List of PeakList objects or PeakMatrix object
    :param path_out: Path to a text file
    :param delimiter:
    :param qc_label:
    :return:
    :rtype: None
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
