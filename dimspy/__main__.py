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


import argparse
import os
import zipfile

import h5py
from dimspy import __version__

from . import tools
from .portals import hdf5_portal


def map_delimiter(delimiter: str):  # pragma: no cover
    seps = {"comma": ",", "tab": "\t"}
    if delimiter in seps:
        return seps[delimiter]
    else:
        return delimiter


def main():  # pragma: no cover

    print(("Executing dimspy version %s." % __version__))

    parser = argparse.ArgumentParser(description='Python package for processing direct-infusion mass spectrometry-based metabolomics and lipidomics data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')

    parser_ps = subparsers.add_parser('process-scans', help='Process scans and/or stitch SIM windows.')
    parser_rf = subparsers.add_parser('replicate-filter',
                                      help='Filter irreproducible peaks from technical replicate peaklists.')
    parser_as = subparsers.add_parser('align-samples', help='Align peaklists across samples.')
    parser_bf = subparsers.add_parser('blank-filter',
                                      help='Filter peaks across samples that are present in the blank samples.')
    parser_sf = subparsers.add_parser('sample-filter',
                                      help='Filter peaks based on certain reproducibility and sample class criteria.')
    parser_rs = subparsers.add_parser('remove-samples',
                                      help='Remove sample(s) from a peak matrix object or list of peaklist objects.')
    parser_mvsf = subparsers.add_parser('mv-sample-filter',
                                        help='Filter samples based on the percentage of missing values.')
    parser_mp = subparsers.add_parser('merge-peaklists',
                                      help='Merge peaklists from multiple lists of peaklist or peak matrix objects.')
    parser_gp = subparsers.add_parser('get-peaklists', help='Get peaklists from a peak matrix object.')
    parser_gap = subparsers.add_parser('get-average-peaklist',
                                       help='Get an average peaklist from a peak matrix object.')
    parser_hpmt = subparsers.add_parser('hdf5-pm-to-txt', help='Write HDF5 output (peak matrix) to text format.')
    parser_hplt = subparsers.add_parser('hdf5-pls-to-txt', help='Write HDF5 output (peak lists) to text format.')
    parser_csl = subparsers.add_parser('create-sample-list',
                                       help='Create a sample list from a peak matrix object or list of peaklist objects.')
    parser_un = subparsers.add_parser('unzip', help='Extract files from zip file')

    parser_lic = subparsers.add_parser('licenses', help='Show licenses DIMSpy and RawFileReader')


    #################################
    # PROCESS SCANS
    #################################

    parser_ps.add_argument('-i', '--input',
                           type=str, action='append', required=True, metavar='source',
                           help="Directory (*.raw, *.mzml or tab-delimited peaklist files), single *.mzml/*.raw file or zip archive (*.mzml only)")

    parser_ps.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peaklist objects to.")

    parser_ps.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-delimited file that include the name of the data files (*.raw or *.mzml) and meta data. "
                                "Column names: filename, replicate, batch, injectionOrder, classLabel.")

    parser_ps.add_argument('-m', '--function-noise',
                           choices=["median", "mean", "mad", "noise_packets"], required=True,
                           help="Select function to calculate noise.")

    parser_ps.add_argument('-s', '--snr-threshold',
                           default=3.0, type=float, required=True,
                           help="Signal-to-noise threshold")

    parser_ps.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in Parts per million to group peaks across scans / mass spectra.")

    parser_ps.add_argument('-n', '--min_scans',
                           default=1, type=int, required=False,
                           help="Minimum number of scans required for each m/z range or event.")

    parser_ps.add_argument('-a', '--min-fraction',
                           default=0.5, type=float, required=False,
                           help="Minimum fraction a peak has to be present. Use 0.0 to not apply this filter.")

    parser_ps.add_argument('-d', '--rsd-threshold',
                           default=None, type=float, required=False,
                           help="Maximum threshold - relative standard deviation (Calculated for peaks that have been measured across a minimum of two scans).")

    parser_ps.add_argument('-k', '--skip-stitching',
                           action='store_true', required=False,
                           help="Skip the step where (SIM) windows are 'stitched' or 'joined' together. Individual peaklists are generated for each window.")

    parser_ps.add_argument('-r', '--ringing-threshold',
                           default=None, type=float, required=False,
                           help="Ringing")

    parser_ps.add_argument('-e', '--include-scan-events',
                           action='append', nargs=3, required=False,
                           metavar=('start', 'end', 'scan_type'), default=[],
                           help="Scan events to select. E.g. 100.0 200.0 sim or 50.0 1000.0 full")

    parser_ps.add_argument('-x', '--exclude-scan-events',
                           action='append', nargs=3, required=False,
                           metavar=('start', 'end', 'scan_type'), default=[],
                           help="Scan events to select. E.g. 100.0 200.0 sim  or  50.0 1000.0 full")

    parser_ps.add_argument('-z', '--remove-mz-range',
                           action='append', nargs=2, required=False, metavar=('start', 'end'), default=[],
                           help="M/z range(s) to remove. E.g. 100.0 102.0  or  140.0 145.0.")

    parser_ps.add_argument('-u', '--report',
                           type=str, required=False, default=None,
                           help="Summary/Report of processed mass spectra")

    parser_ps.add_argument('-b', '--block-size',
                           default=5000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_ps.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units (CPUs).")

    #################################
    # REPLICATE FILTER
    #################################

    parser_rf.add_argument('-i', '--input',
                           type=str, required=True,
                           help="HDF5 file (Peaklist objects) from step 'process-scans' or directory path that contains tab-delimited peaklists.")

    parser_rf.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peaklist objects to.")

    parser_rf.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in Parts per million to group peaks across scans / mass spectra.")

    parser_rf.add_argument('-r', '--replicates',
                           default=2,
                           type=int,
                           required=True,
                           help="Number of technical replicates.")

    parser_rf.add_argument('-m', '--min-peak-present',
                           default=2,
                           type=int,
                           required=True,
                           help="Minimum number of times a peak has to be present (number) across technical replicates.")

    parser_rf.add_argument('-d', '--rsd-threshold',
                           default=None, type=float, required=False,
                           help="Maximum threshold - Relative Standard Deviation.")

    parser_rf.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-delimited file that list all the data files (*.raw or *.mzml) and meta data (filename, technical replicate, class, batch).")

    parser_rf.add_argument('-u', '--report',
                           type=str, required=False, default=None,
                           help="Summary/Report of processed mass spectra")

    parser_rf.add_argument('-b', '--block-size',
                           default=5000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_rf.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units (CPUs).")

    #################################
    # Align Samples
    #################################

    parser_as.add_argument('-i', '--input',
                           type=str, required=True,
                           help="HDF5 file (Peaklist objects) from step 'process-scans / replicate-filter' "
                                "or directory path that contains tab-delimited peaklists.")

    parser_as.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peak matrix object to.")

    parser_as.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in parts per million to group peaks across scans / mass spectra.")

    parser_as.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-delimited file that include the name of the samples and meta data."
                                "Column names: filename, replicate, batch, injectionOrder, classLabel.")

    parser_as.add_argument('-b', '--block-size',
                           default=5000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_as.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units (CPUs).")

    #################################
    # Blank Filter
    #################################

    parser_bf.add_argument('-i', '--input',
                           type=str, required=True,
                           help="HDF5 file or tab-delimited file that contains a peak matrix (object).")

    parser_bf.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peak matrix object to.")

    parser_bf.add_argument('-l', '--blank-label',
                           default="blank", type=str, required=True,
                           help="Class label for blanks.")

    parser_bf.add_argument('-m', '--min-fraction',
                           default=1.0, type=float, required=False,
                           help="Minium fold change blank versus sample.")

    parser_bf.add_argument('-f', '--function',
                           default="mean", choices=["mean", "median", "max"], required=False,
                           help="Select function to calculate blank intenstiy.")

    parser_bf.add_argument('-c', '--min-fold-change',
                           default=1.0, type=float, required=False,
                           help="Minium fold change blank versus sample.")

    parser_bf.add_argument('-r', '--remove-blank-samples',
                           action='store_true', required=False,
                           help="Remove blank samples from peak matrix.")

    parser_bf.add_argument('-a', '--labels',
                           type=str, required=False,
                           help="Tab delimited file with at least two columns named 'filename' and 'classLabel'.")

    #################################
    # Sample filter
    #################################

    parser_sf.add_argument('-i', '--input',
                           type=str, required=True,
                           help="HDF5 file or tab-delimited file that contains a peak matrix.")

    parser_sf.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peak matrix object to.")

    parser_sf.add_argument('-p', '--min-fraction',
                           type=float, required=False, default=0.5,
                           help="Minimum percentage of samples a peak has to be present.")

    parser_sf.add_argument('-w', '--within',
                           action='store_true', required=False,
                           help="Apply sample filter within each sample class.")

    parser_sf.add_argument('-d', '--rsd-threshold',
                           default=None, type=float, required=False,
                           help="Peaks where the associated QC peaks are above this threshold will be removed.")

    parser_sf.add_argument('-q', '--qc-label',
                           default=None, type=str, required=False,
                           help="Class label for QCs")

    parser_sf.add_argument('-a', '--labels',
                           type=str, required=False,
                           help="Tab delimited file with at least two columns named 'filename' and 'classLabel'.")

    #################################
    # Missing Values Sample filter
    #################################

    parser_mvsf.add_argument('-i', '--input',
                             type=str, required=True,
                             help="HDF5 file file that contains a peak matrix object.")

    parser_mvsf.add_argument('-o', '--output',
                             type=str, required=True,
                             help="HDF5 file to save the peak matrix object to.")

    parser_mvsf.add_argument('-m', '--max-fraction',
                             type=float, required=False,
                             help="Maximum percentage of missing values allowed across a sample.")

    #################################
    # Remove Samples
    #################################

    parser_rs.add_argument('-i', '--input',
                           type=str, action='append', required=True, metavar='source',
                           help="HDF5 file that contains a peak matrix object or list of peaklist objects from one of the processing steps.")

    parser_rs.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peak matrix object or peaklist objects to.")

    parser_rs.add_argument('-s', '--sample-names',
                           type=str, action='append', required=True,
                           help="Sample name(s)")

    #################################
    # Merge peaklists
    #################################

    parser_mp.add_argument('-i', '--input',
                           required=True, default=[], action='append',
                           help="Multiple HDF5 files that contain peaklists or peak matrix from one of the processing steps.")
    parser_mp.add_argument('-o', '--output',
                           type=str, required=True,
                           help="Directory (if using multilist column in filelist) or HDF5 file to write to.")
    parser_mp.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-delimited file that list all the data files (*.raw or *.mzml) and meta data "
                                "(filename, technical replicate, class, batch, multiList).")

    #################################
    # Get peaklists
    #################################

    parser_gp.add_argument('-i', '--input',
                           required=True, default=[], action='append',
                           help="Single or Multiple HDF5 files that contain a peak matrix object from one of the processing steps.")

    parser_gp.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peaklist objects to.")

    #########################################
    # Get average peaklist from peak matrix
    #########################################

    parser_gap.add_argument('-i', '--input',
                            type=str, required=True,
                            help="Single or Multiple HDF5 files that contain a peak matrix object from one of the processing steps.")

    parser_gap.add_argument('-o', '--output',
                            type=str, required=True,
                            help="HDF5 file to save the peaklist object to.")

    parser_gap.add_argument('-n', '--name-peaklist',
                            type=str, required=True,
                            help="Name of the peaklist.")

    #################################
    # HDF5 peak matrix to text
    #################################

    parser_hpmt.add_argument('-i', '--input',
                             type=str, required=True,
                             help="HDF5 file that contains a peak matrix object from one of the processing steps.")

    parser_hpmt.add_argument('-o', '--output',
                             type=str, required=True,
                             help="Directory (peaklists) or text file (peak matrix) to write to.")

    parser_hpmt.add_argument('-a', '--attribute_name',
                             default="intensity", choices=["intensity", "mz", "snr"], required=False,
                             help="Type of matrix to print.")

    parser_hpmt.add_argument('-l', '--class-label-rsd',
                             action='append', required=False,
                             default=(),
                             help="Class label to select samples for RSD calculatons (e.g. QC).")

    parser_hpmt.add_argument('-d', '--delimiter',
                             default="tab", choices=["tab", "comma"],
                             help="Values on each line of the file are separated by this character.")

    parser_hpmt.add_argument('-s', '--representation-samples',
                             default="rows", choices=["rows", "columns"],
                             help="Should the rows or columns respresent the samples?")

    parser_hpmt.add_argument('-c', '--comprehensive',
                             action='store_true', required=False,
                             help="Comprehensive version of the peak matrix")

    #################################
    # HDF5 peaklists to text
    #################################

    parser_hplt.add_argument('-i', '--input',
                             type=str, required=True,
                             help="HDF5 file that contains a list of peaklist objects from one of the processing steps.")

    parser_hplt.add_argument('-o', '--output',
                             type=str, required=True,
                             help="Directory to write to.")

    parser_hplt.add_argument('-d', '--delimiter',
                             default="tab", choices=["tab", "comma"],
                             help="Values on each line of the file are separated by this character.")

    #################################
    # Create Sample List
    #################################

    parser_csl.add_argument('-i', '--input',
                            type=str, required=True,
                            help="HDF5 file that contains a peak matrix object from one of the processing steps.")

    parser_csl.add_argument('-o', '--output',
                            type=str, required=True,
                            help="Text file to write to.")

    parser_csl.add_argument('-d', '--delimiter',
                            default="tab", choices=["tab", "comma"],
                            help="Values on each line of the file are separated by this character.")

    #################################
    # Unzip archive
    #################################

    parser_un.add_argument('-i', '--input',
                            type=str, required=True,
                            help="file[.zip]")

    parser_un.add_argument('-o', '--output',
                            type=str, required=True,
                            help="Directory to write to.")

    args = parser.parse_args()

    print(args)

    if args.step == "process-scans":

        filter_scan_events = {}
        if args.exclude_scan_events != [] and args.include_scan_events != []:
            raise argparse.ArgumentTypeError(
                "-e/--include-scan-events and -x/--exclude-scan-events can not be used together.")
        elif args.exclude_scan_events != []:
            for se in args.exclude_scan_events:
                if "exclude" not in filter_scan_events:
                    filter_scan_events["exclude"] = []
                filter_scan_events["exclude"].append([se[0], se[1], se[2]])
        elif args.include_scan_events != []:
            for se in args.include_scan_events:
                if "include" not in filter_scan_events:
                    filter_scan_events["include"] = []
                filter_scan_events["include"].append([se[0], se[1], se[2]])

        remove_mz_range = [[float(mzr[0]), float(mzr[1])] for mzr in args.remove_mz_range]

        if len(args.input) == 1:  # Directory / zipfile / single filename
            args.input = args.input[0]

        peaklists = tools.process_scans(source=args.input,
                                        function_noise=args.function_noise,
                                        snr_thres=args.snr_threshold,
                                        ppm=args.ppm,
                                        min_fraction=args.min_fraction,
                                        rsd_thres=args.rsd_threshold,
                                        min_scans=args.min_scans,
                                        filelist=args.filelist,
                                        skip_stitching=args.skip_stitching,
                                        ringing_thres=args.ringing_threshold,
                                        filter_scan_events=filter_scan_events,
                                        remove_mz_range=remove_mz_range,
                                        report=args.report,
                                        block_size=args.block_size,
                                        ncpus=args.ncpus)
        hdf5_portal.save_peaklists_as_hdf5(peaklists, args.output)

    elif args.step == "replicate-filter":
        peaklists_rf = tools.replicate_filter(source=args.input,
                                              ppm=args.ppm,
                                              replicates=args.replicates,
                                              min_peaks=args.min_peak_present,
                                              rsd_thres=args.rsd_threshold,
                                              filelist=args.filelist,
                                              report=args.report,
                                              block_size=args.block_size,
                                              ncpus=args.ncpus)
        hdf5_portal.save_peaklists_as_hdf5(peaklists_rf, args.output)

    elif args.step == "align-samples":
        pm = tools.align_samples(source=args.input,
                                 ppm=args.ppm,
                                 filelist=args.filelist,
                                 block_size=args.block_size,
                                 ncpus=args.ncpus)
        hdf5_portal.save_peak_matrix_as_hdf5(pm, args.output)

    elif args.step == "blank-filter":
        pm_bf = tools.blank_filter(peak_matrix=args.input,
                                   blank_label=args.blank_label,
                                   min_fraction=args.min_fraction,
                                   min_fold_change=args.min_fold_change,
                                   function=args.function,
                                   rm_samples=args.remove_blank_samples,
                                   labels=args.labels)
        hdf5_portal.save_peak_matrix_as_hdf5(pm_bf, args.output)

    elif args.step == "sample-filter":
        pm_sf = tools.sample_filter(peak_matrix=args.input,
                                    min_fraction=args.min_fraction,
                                    within=args.within,
                                    rsd_thres=args.rsd_threshold,
                                    qc_label=args.qc_label,
                                    labels=args.labels)
        hdf5_portal.save_peak_matrix_as_hdf5(pm_sf, args.output)

    elif args.step == "mv-sample-filter":
        pm = hdf5_portal.load_peak_matrix_from_hdf5(args.input)
        pm_mvf = tools.missing_values_sample_filter(pm, args.max_fraction)
        hdf5_portal.save_peak_matrix_as_hdf5(pm_mvf, args.output)

    elif args.step == "remove-samples":
        f = h5py.File(args.input, 'r')
        if "mz" in f:
            pm = hdf5_portal.load_peak_matrix_from_hdf5(args.input)
            pm_rs = tools.remove_samples(pm, args.sample_names)
            hdf5_portal.save_peak_matrix_as_hdf5(pm_rs, args.output)
        else:
            pls = hdf5_portal.load_peaklists_from_hdf5(args.input)
            pls_rs = tools.remove_samples(pls, args.sample_names)
            hdf5_portal.save_peak_matrix_as_hdf5(pls_rs, args.output)

    elif args.step == "merge-peaklists":

        pls_merged = tools.merge_peaklists(source=args.input, filelist=args.filelist)

        # if a list of lists, save each as separate list of peaklists
        if any(isinstance(l, list) for l in pls_merged):
            for i in range(len(pls_merged)):
                if 'multilistLabel' in list(pls_merged[i][0].metadata.keys()):
                    m_label = pls_merged[i][0].metadata['multilistLabel']
                    hdf5_portal.save_peaklists_as_hdf5(pls_merged[i],
                                                       os.path.join(args.output, '{}.hdf5'.format(m_label)))
                else:
                    m_nm = pls_merged[i][0].metadata['multilist']
                    hdf5_portal.save_peaklists_as_hdf5(pls_merged[i],
                                                       os.path.join(args.output,
                                                                    'merged_peaklist_{:03d}.hdf5'.format(m_nm)))
        else:
            hdf5_portal.save_peaklists_as_hdf5(pls_merged, args.output)

    elif args.step == "get-peaklists":
        pls = []
        for s in args.input:
            pm = hdf5_portal.load_peak_matrix_from_hdf5(s)
            pls.extend(pm.extract_peaklists())
        hdf5_portal.save_peaklists_as_hdf5(pls, args.output)

    elif args.step == "get-average-peaklist":
        pls = [hdf5_portal.load_peak_matrix_from_hdf5(args.input).to_peaklist(ID=args.name_peaklist)]
        hdf5_portal.save_peaklists_as_hdf5(pls, args.output)

    elif args.step == "hdf5-pm-to-txt":
        if args.representation_samples == "rows":
            samples_in_rows = True
        else:
            samples_in_rows = False
        tools.hdf5_peak_matrix_to_txt(args.input,
                                      path_out=args.output,
                                      attr_name=args.attribute_name,
                                      delimiter=map_delimiter(args.delimiter),
                                      rsd_tags=args.class_label_rsd,
                                      samples_in_rows=samples_in_rows,
                                      comprehensive=args.comprehensive)

    elif args.step == "hdf5-pls-to-txt":
        tools.hdf5_peaklists_to_txt(args.input, path_out=args.output, delimiter=map_delimiter(args.delimiter))

    elif args.step == "create-sample-list":
        try:
            inp = hdf5_portal.load_peaklists_from_hdf5(args.input)
        except:
            inp = hdf5_portal.load_peak_matrix_from_hdf5(args.input)
        tools.create_sample_list(inp, args.output, delimiter=map_delimiter(args.delimiter))

    elif args.step == "unzip":
        with zipfile.ZipFile(args.input, 'r') as zip_ref:
            zip_ref.extractall(args.output)

    elif args.step == "licenses":
        print("""
DIMSpy is licensed under the GNU General Public License v3.0. Copyright © 2017 - 2020 Ralf Weber, Albert Zhou

RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved. 

Using DIMSpy software for processing Thermo Fisher Scientific *.raw files implies the acceptance of the RawFileReader license terms.

Anyone receiving RawFileReader as part of a larger software distribution (in the current context, as part of DIMSpy) is considered an "end user" under
section 3.3 of the RawFileReader License, and is not granted rights to redistribute RawFileReader.

This license (see 'SOFTWARE LICENSE AGREEMENT' below) covers the following files which are distributed with the DIMSpy software package:
 - ThermoFisher.CommonCore.BackgroundSubtraction.dll
 - ThermoFisher.CommonCore.BackgroundSubtraction.XML
 - ThermoFisher.CommonCore.Data.dll
 - ThermoFisher.CommonCore.Data.XML
 - ThermoFisher.CommonCore.MassPrecisionEstimator.dll
 - ThermoFisher.CommonCore.MassPrecisionEstimator.XML
 - ThermoFisher.CommonCore.RawFileReader.dll
 - ThermoFisher.CommonCore.RawFileReader.XML

**SOFTWARE LICENSE AGREEMENT ("License") FOR RawFileReader**
----------------------------------------------------------------------
These License terms are an agreement between you and Thermo Finnigan LLC ("Licensor"). They apply to Licensor's MSFileReader software program ("Software"), which includes documentation and any media on which you received it. These terms also apply to any updates or supplements for this Software, unless other terms accompany those items, in which case those terms apply. **If you use this Software, you accept this License. If you do not accept this License, you are prohibited from using this software.  If you comply with these License terms, you have the rights set forth below.**

1. Rights Granted:

1.1. You may install and use this Software on any of your computing devices.

1.2. You may distribute this Software to others, but only in combination with other software components and/or programs that you provide and subject to the distribution requirements and restrictions below.

2.  Use Restrictions:

2.1. You may not decompile, disassemble, reverse engineer, use reflection or modify this Software.

3. Distribution Requirements:

If you distribute this Software to others, you agree to:

3.1. Indemnify, defend and hold harmless the Licensor from any claims, including attorneys&#39; fees, related to the distribution or use of this Software;

3.2. Display the following text in your software's "About"; box: "**RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved**.";

3.3. Require your end users to agree to a license agreement that prohibits them from redistributing this Software to others.

4.  Distribution Restrictions:

4.1. You may not use the Licensor&#39;s trademarks in a way that suggests your software components and/or programs are provided by or are endorsed by the Licensor; and

4.2. You may not commercially exploit this Software or products that incorporate this Software without the prior written consent of Licensor. Commercial exploitation includes, but is not limited to, charging a purchase price, license fee, maintenance fee, or subscription fee; or licensing, transferring or redistributing the Software in exchange for consideration of any kind.

4.3. Your rights to this Software do not include any license, right, power or authority to subject this Software in whole or in part to any of the terms of an Excluded License. &quot;Excluded License&quot; means any license that requires as a condition of use, modification and/or distribution of software subject to the Excluded License, that such software or other software combined and/or distributed with such software be (a) disclosed or distributed in source code form; or (b) licensed for the purpose of making derivative works.  Without limiting the foregoing obligation, you are specifically prohibited from distributing this Software with any software that is subject to the General Public License (GPL) or similar license in a manner that would create a combined work.

5.  Additional Terms Applicable to Software:

5.1. This Software is licensed, not sold. This License only gives you some rights to use this Software; the Licensor reserves all other rights. Unless applicable law gives you more rights despite this limitation, you may use this Software only as expressly permitted in this License.

5.2. Licensor has no obligation to fix, update, supplement or support this Software.

5.3. This Software is not designed, manufactured or intended for any use requiring fail-safe performance in which the failure of this Software could lead to death, serious personal injury or severe physical and environmental damage (&quot;High Risk Activities&quot;), such as the operation of aircraft, medical or nuclear facilities. You agree not to use, or license the use of, this Software in connection with any High Risk Activities.

5.4. Your rights under this License terminate automatically if you breach this License in any way. Termination of this License will not affect any of your obligations or liabilities arising prior to termination. The following sections of this License shall survive termination: 2.1, 3.1, 3.2, 3.3, 4.1, 4.2, 4.3, 5.1, 5.2, 5.3, 5.5, 5.6, 5.7, 5.8, and 5.9.

5.5. This Software is subject to United States export laws and regulations. You agree to comply with all domestic and international export laws and regulations that apply to this Software. These laws include restrictions on destinations, end users and end use.

5.6. This License shall be construed and controlled by the laws of the State of California, U.S.A., without regard to conflicts of law. You consent to the jurisdiction of the state and federal courts situated in the State of California in any action arising under this License. The application of the U.N. Convention on Contracts for the International Sale of Goods to this License is hereby expressly excluded. If any provision of this License shall be deemed unenforceable or contrary to law, the rest of this License shall remain in full effect and interpreted in an enforceable manner that most nearly captures the intent of the original language.

5.7. THIS SOFTWARE IS LICENSED &quot;AS IS&quot;. YOU BEAR ALL RISKS OF USING IT. LICENSOR GIVES NO AND DISCLAIMS ALL EXPRESS AND IMPLIED WARRANTIES, REPRESENTATIONS OR GUARANTEES.  YOU MAY HAVE ADDITIONAL CONSUMER RIGHTS UNDER YOUR LOCAL LAWS WHICH THIS LICENSE CANNOT CHANGE. TO THE EXTENT PERMITTED UNDER YOUR LOCAL LAWS, LICENSOR EXCLUDES THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT.

5.8. LICENSOR&#39;S TOTAL LIABILITY TO YOU FOR DIRECT DAMAGES ARISING UNDER THIS LICENSE IS LIMITED TO U.S. $1.00. YOU CANNOT RECOVER ANY OTHER DAMAGES, INCLUDING CONSEQUENTIAL, LOST PROFITS, SPECIAL, INDIRECT OR INCIDENTAL DAMAGES, EVEN IF LICENSOR IS EXPRESSLY MADE AWARE OF THE POSSIBILITY THEREOF OR IS NEGLIGENT. THIS LIMITATION APPLIES TO ANYTHING RELATED TO THIS SOFTWARE, SERVICES, CONTENT (INCLUDING CODE) ON THIRD PARTY INTERNET SITES, OR THIRD PARTY PROGRAMS, AND CLAIMS FOR BREACH OF CONTRACT, BREACH OF WARRANTY, GUARANTEE  OR CONDITION, STRICT LIABILITY, NEGLIGENCE, OR OTHER TORT TO THE EXTENT PERMITTED BY APPLICABLE LAW.

5.9. Use, duplication or disclosure of this Software by the U.S. Government is subject to the restricted rights applicable to commercial computer software (under FAR 52.227019 and DFARS 252.227-7013 or parallel regulations). The manufacturer for this purpose is Thermo Finnigan LLC, 355 River Oaks Parkway, San Jose, California 95134, U.S.A.
""")


if __name__ == "__main__":
    main()
