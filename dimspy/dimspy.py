#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import workflow
from portals import hdf5_portal
from . import __version__
from experiment import check_metadata
from experiment import update_metadata


def main():

    print("Executing dimspy version %s." % __version__)

    parser = argparse.ArgumentParser(description='Python package to process DIMS data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')

    parser_cf = subparsers.add_parser('check-filelist', help='Validate samplelist (filename, replicates, order, batch, class, qc, blank).')
    parser_ps = subparsers.add_parser('process-scans', help='Process scans and/or stitch windows within MS data file.')
    parser_rf = subparsers.add_parser('replicate-filter', help='Filter irreproducible peaks from technical replicate peaklists.')
    parser_as = subparsers.add_parser('align-samples', help='Align mass spectra across samples.')
    parser_bf = subparsers.add_parser('blank-filter', help='Filter peaks present in the blank samples.')
    parser_sf = subparsers.add_parser('sample-filter', help='Filter peaks based on certain reproducibility and sample class criteria.')
    parser_mp = subparsers.add_parser('merge-peaklists', help='Merge multiple HDF5 files that contain peaklist objects.')
    parser_ht = subparsers.add_parser('hdf5-to-txt', help='Write HDF5 output to text format.')

    parser_cf.add_argument('-l', '--filelist',
                           type=str, required=True,
                           help="Tab-delimited file that list all the MS data files (*.raw, *.mzml or tab-delimited files) and associated meta data (e.g. replicate, order, batch, class, qc, blank).")

    parser_cf.add_argument('-i', '--input',
                           type=str, required=True,
                           help="Directory (*.raw, *.mzml or tab-delimited peaklist files) or zip archive (*.mzml or tab-delimited peaklist files)")

    parser_cf.add_argument('-r', '--replicates',
                           type=int, required=False,
                           help="Number of technical replicates for each sample.")

    parser_cf.add_argument('-a', '--batches',
                           type=int, required=False,
                           help="Number of batches")

    parser_cf.add_argument('-q', '--name-QC',
                           type=str,  required=False,
                           help="Only required when QCs are not marked in the filelist (i.e. QC column).")

    parser_cf.add_argument('-b', '--name-blank',
                           type=str, required=False,
                           help="Only required when 'blanks' are not marked in the filelist (i.e. blank column).")

    #################################
    # PROCESS SCANS
    #################################

    parser_ps.add_argument('-i', '--input',
                           type=str, action='append', required=True, metavar=('source'),
                           help="Directory (*.raw, *.mzml or tab-delimited peaklist files), single *.mzml/*.raw file or zip archive (*.mzml only)")

    parser_ps.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peaklist objects.")

    parser_ps.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-separated file that list all the data files (*.raw or *.mzml) and meta data (filename, technical replicate, class, batch). "
                                "HIGHLY RECOMMENDED when directory or zip file is provided.")

    parser_ps.add_argument('-m', '--function-noise',
                           choices=["median", "mean", "mad", "msfilereader"], required=True,
                           help="Select function to calculate noise.")

    parser_ps.add_argument('-s', '--snr-threshold',
                           default=3.0,  type=float,  required=True,
                           help="Signal-to-noise threshold")

    parser_ps.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in Parts per million to group peaks across scans / mass spectra.")

    parser_ps.add_argument('-n', '--min_scans',
                           default=1, type=int, required=False,
                           help="Minimum number of scans required for each m/z range or event (header).")

    parser_ps.add_argument('-a', '--min-fraction',
                           default=1, type=float, required=False,
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

    parser_ps.add_argument('-u', '--include-scan-events',
                           action='append', nargs=3, required=False, 
                           metavar=('start', 'end', 'scan_type'), default=[],
                           help="Scan events to select. E.g. 100.0 200.0 sim  or  50.0 1000.0 full")

    parser_ps.add_argument('-x', '--exclude-scan-events',
                           action='append', nargs=3, required=False, 
                           metavar=('start', 'end', 'scan_type'), default=[],
                           help="Scan events to select. E.g. 100.0 200.0 sim  or  50.0 1000.0 full")

    parser_ps.add_argument('-z', '--remove-mz-range',
                           action='append', nargs=2, required=False, metavar=('start','end'), default=[],
                           help="M/z range(s) to remove. E.g. 100.0 102.0  or  140.0 145.0.")

    parser_ps.add_argument('-b', '--block-size',
                           default=2000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_ps.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units (CPUs).")

    #################################
    # REPLICATE FILTER
    #################################

    parser_rf.add_argument('-i', '--input',
                           type=str, required=False,
                           help="HDF5 file (Peaklist objects) from step 'process-scans' or directory path that contains tab-delimited peaklists.")

    parser_rf.add_argument('-o', '--output',
                           type=str, required=False,
                           help="HDF5 file to save the peaklist objects.")

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
                           default=None, type=float,  required=False,
                           help="Maximum threshold - Relative Standard Deviation.")

    parser_rf.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-separated file that list all the data files (*.raw or *.mzml) and meta data (filename, technical replicate, class, batch).")

    parser_rf.add_argument('-b', '--block-size',
                           default=2000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_rf.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units (CPUs).")

    #################################
    # Align Samples
    #################################

    parser_as.add_argument('-i', '--input',
                           type=str, required=False,
                           help="HDF5 file (Peaklist objects) from step 'process-scans / replicate-filter' or directory path that contains tab-delimited peaklists.")

    parser_as.add_argument('-o', '--output',
                           type=str, required=False,
                           help="HDF5 file to save the peak matrix object.")

    parser_as.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in parts per million to group peaks across scans / mass spectra.")

    parser_as.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-separated file that list all the data files (*.raw or *.mzml) and meta data (filename, technical replicate, class, batch).")

    parser_as.add_argument('-b', '--block-size',
                           default=2000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_as.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units (CPUs).")

    #################################
    # Blank Filter
    #################################

    parser_bf.add_argument('-i', '--input',
                           type=str, required=False,
                           help="HDF5 file or tab-delimited file that contains a peak matrix (object).")

    parser_bf.add_argument('-o', '--output',
                           type=str, required=False,
                           help="HDF5 file to save the peak matrix object.")

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

    parser_bf.add_argument('-a', '--class-labels',
                           type=str, required=False,
                           help="Tab delimited file (two columns) with the filenames / sample ids in the first columns and class label in the second column.")

    #################################
    # Sample filter
    #################################

    parser_sf.add_argument('-i', '--input',
                           type=str, required=False,
                           help="HDF5 file or tab-delimited file that contains a peak matrix.")

    parser_sf.add_argument('-o', '--output',
                           type=str, required=True,
                           help="HDF5 file to save the peak matrix object.")

    parser_sf.add_argument('-p', '--min-fraction',
                           type=float, required=False,
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

    parser_sf.add_argument('-a', '--class-labels',
                           type=str, required=False,
                           help="Tab delimited file (two columns) with the filenames / sample ids in the first columns and class label in the second column.")

    #################################
    # Merge peaklists
    #################################

    parser_mp.add_argument('-i', '--input',
                           required=True, default=[], action='append',
                           help="Multiple HDF5 files that contain peaklists from one of the processing steps.")
    parser_mp.add_argument('-o', '--output',
                           required=True, type=str,
                           help="HDF5 file to save the merged peaklist objects.")
    parser_mp.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="Tab-separated file that list all the data files (*.raw or *.mzml) and meta data (filename, technical replicate, class, batch).")

    #################################
    # HDF5 to text
    #################################

    parser_ht.add_argument('-i', '--input',
                           type=str, required=True,
                           help="HDF5 file that contains peaklists or a peak matrix from one of the processing steps.")

    parser_ht.add_argument('-o', '--output',
                           type=str, required=True,
                           help="Directory (peaklists) or text file (peak matrix) to write to.")

    parser_ht.add_argument('-a', '--attribute_name',
                           default="intensity", choices=["intensity", "mz"], required=False,
                           help="Type of matrix to print.")

    parser_ht.add_argument('-s', '--separator',
                           default="tab", choices=["tab", "comma"],
                           help="Format of the file for further data processing or data analysis.")

    parser_ht.add_argument('-t', '--transpose',
                           action='store_true', required=False,
                           help="Transpose peak matrix")

    parser_ht.add_argument('-c', '--comprehensive',
                           action='store_true', required=False,
                           help="Comprehensive output of the peak matrix")


    args = parser.parse_args()
    print args

    if args.step == "process-scans":
        filter_scan_events = {}
        if args.exclude_scan_events != []:
            for se in args.exclude_scan_events: 
                if "exclude" not in filter_scan_events:
                    filter_scan_events["exclude"] = []
                filter_scan_events["exclude"].append([se[0], se[1], se[2]])
        elif args.include_scan_events != []:
            for se in args.include_scan_events: 
                if "exclude" not in filter_scan_events:
                    filter_scan_events["include"] = []
                filter_scan_events["include"].append([se[0], se[1], se[2]])
        elif args.exclude_scan_events != [] and args.include_scan_events != []:
            raise argparse.ArgumentTypeError("-u/--include-scan-events and -x/--exclude-scan-events can not be used together.")

        remove_mz_range = [[float(mzr[0]), float(mzr[1])] for mzr in args.remove_mz_range]

        if len(args.input) == 1: # Directory / zipfile / single filename
            args.input = args.input[0]

        peaklists = workflow.process_scans(source=args.input,
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
                                           block_size=args.block_size,
                                           ncpus=args.ncpus)
        hdf5_portal.save_peaklists_as_hdf5(peaklists, args.output)

    elif args.step == "mass-calibrate":
        # TODO implement mass calibration
        # TODO methods available in python
        print args.method

    elif args.step == "replicate-filter":
        peaklists_rf = workflow.replicate_filter(source=args.input,
                                                 ppm=args.ppm,
                                                 replicates=args.replicates,
                                                 min_peaks=args.min_peak_present,
                                                 rsd_thres=args.rsd_threshold,
                                                 filelist=args.filelist,
                                                 block_size=args.block_size,
                                                 ncpus=args.ncpus)
        hdf5_portal.save_peaklists_as_hdf5(peaklists_rf, args.output)

    elif args.step == "align-samples":
        pm = workflow.align_samples(source=args.input,
                                    ppm=args.ppm,
                                    filelist=args.filelist,
                                    block_size=args.block_size,
                                    ncpus=args.ncpus)
        hdf5_portal.save_peak_matrix_as_hdf5(pm, args.output)
    elif args.step == "blank-filter":
        pm_bf = workflow.blank_filter(peak_matrix=args.input,
                                      blank_label=args.blank_label,
                                      min_fraction=args.min_fraction,
                                      min_fold_change=args.min_fold_change,
                                      function=args.function,
                                      rm_samples=args.remove_blank_samples,
                                      class_labels=args.class_labels)
        hdf5_portal.save_peak_matrix_as_hdf5(pm_bf, args.output)

    elif args.step == "sample-filter":
        pm_sf = workflow.sample_filter(peak_matrix=args.input,
                                       min_fraction=args.min_fraction,
                                       within=args.within,
                                       rsd=args.rsd_threshold,
                                       qc_label=args.qc_label,
                                       class_labels=args.class_labels)
        hdf5_portal.save_peak_matrix_as_hdf5(pm_sf, args.output)

    elif args.step == "merge-peaklists":
        pls_merged = []
        for fn_hdf5 in args.input:
            pls = hdf5_portal.load_peaklists_from_hdf5(fn_hdf5)
            if args.filelist is not None:
                fl = check_metadata(args.filelist)
                pls = update_metadata(pls, fl)
            pls_merged.extend(pls)
        hdf5_portal.save_peaklists_as_hdf5(pls_merged, args.output)

    elif args.step == "hdf5-to-txt":
        workflow.hdf5_to_txt(args.input,
                             path_out=args.output,
                             attr_name=args.attribute_name,
                             separator=args.separator,
                             transpose=args.transpose,
                             comprehensive=args.comprehensive)
