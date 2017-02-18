#!/usr/bin/python

"""dimspy:dimspy provides entry point main()."""

import argparse
import sys
import workflow
import cPickle as pickle
from . import __version__

from models.peak_matrix import PeakMatrix
from models.peaklist import PeakList


def main():

    print("Executing dimspy version %s." % __version__)

    parser = argparse.ArgumentParser(description='Python package to process DIMS data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')

    parser_cf = subparsers.add_parser('check-filelist', help='Validate samplelist (filename, replicates, order, batch, class, qc, blank)')
    parser_ps = subparsers.add_parser('process-scans', help='Process scans and/or stitch windows within MS data file')
    parser_rf = subparsers.add_parser('replicate-filter', help='Filter irreproducible peaks from technical replicate peaklists')
    parser_as = subparsers.add_parser('align-samples', help='Align mass spectra across samples')
    parser_bf = subparsers.add_parser('blank-filter', help='Filter peaks present in the blank samples')
    parser_sf = subparsers.add_parser('sample-filter', help='Filter peaks based on certain reproducibility and sample class criteria')
    parser_pr = subparsers.add_parser('pickle-to-readable', help='Write output to human-readable format')

    parser_cf.add_argument('-l', '--filelist',
                           type=str, required=True,
                           help="csv or tab-separated file that list all the MS data files (*.raw or *.mzml) and associated meta data (replicate, order, batch, class, qc, blank)")

    parser_cf.add_argument('-i', '--source',
                           type=str, required=True,
                           help="Directory (*.raw or *.mzml), archive (*.zip) or MS data file (*.raw or *.mzml)")

    parser_cf.add_argument('-r', '--replicates',
                           type=int, required=False,
                           help="Number of technical replicates for each sample")

    parser_cf.add_argument('-a', '--batches',
                           type=int, required=False,
                           help="Number of batches")

    parser_cf.add_argument('-q', '--name-QC',
                           type=str,  required=False,
                           help="Only required when QCs are not marked in the filelist (i.e. QC column)")

    parser_cf.add_argument('-b', '--name-blank',
                           type=str, required=False,
                           help="Only required when 'blanks' are not marked in the filelist (i.e. blank column)")

    #################################
    # PROCESS SCANS
    #################################

    parser_ps.add_argument('-i', '--source',
                           type=str, required=True,
                           help="Directory (*.raw or *.mzml), archive (*.zip) or MS data file (*.raw or *.mzml)")

    parser_ps.add_argument('-o', '--pickle-file-out',
                           type=str, required=True,
                           help="Pickle file to store the peaklist objects")

    parser_ps.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="File (csv or tab-separated ) that list all the data files (*.raw or *.mzml) and meta data (filename, technical replicate, group, batch). "
                                "HIGHLY RECOMMENDED when directory or zip file is provided")

    parser_ps.add_argument('-m', '--function-noise',
                           choices=["median", "mean", "mad", "msfilereader"], required=True,
                           help="Select function to calculate noise")

    parser_ps.add_argument('-s', '--snr-threshold',
                           default=3.0,  type=float,  required=True,
                           help="Signal-to-noise threshold")

    parser_ps.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in Parts per million to group peaks across scans / mass spectra")

    parser_ps.add_argument('-n', '--nscans',
                           default=0, type=int, required=False,
                           help="Maximum number of scans to select for each scan event. Use zero for all scans.")

    parser_ps.add_argument('-a', '--min-fraction',
                           default=1, type=float, required=False,
                           help="Minimum number of scans a peak has to be present.")

    parser_ps.add_argument('-d', '--rsd-threshold',
                           default=None, type=float, required=False,
                           help="Maximum threshold - relative standard deviation (Only applied to peaks that have been measured across a minimum of three scans).")

    parser_ps.add_argument('-e', '--filename-experiment',
                           type=str, required=False,
                           help="Filename that contains a description of the DIMS experiment.")

    parser_ps.add_argument('-b', '--block-size',
                           default=2000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_ps.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units.")

    #################################
    # REPLICATE FILTER
    #################################

    parser_rf.add_argument('-i', '--pickle-file-in',
                           type=str, required=False,
                           help="Pickle file that contains peaklist objects")

    parser_rf.add_argument('-o', '--pickle-file-out',
                           type=str, required=False,
                           help="Pickle file to store the peaklist objects")

    parser_rf.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in Parts per million to group peaks across scans / mass spectra.")

    parser_rf.add_argument('-r', '--replicates',
                           default=2,
                           type=int,
                           required=True,
                           help="Number of technical replicates.")

    parser_rf.add_argument('-m', '--min-peaks',
                           default=2,
                           type=int,
                           required=True,
                           help="Minimum number of times a peak has to be present (number) across technical replicates.")

    parser_rf.add_argument('-d', '--rsd-threshold',
                           default=None, type=float,  required=False,
                           help="Maximum threshold - Relative Standard Deviation.")

    parser_rf.add_argument('-b', '--block-size',
                           default=2000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_rf.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units.")

    #################################
    # Align Samples
    #################################

    parser_as.add_argument('-i', '--pickle-file-in',
                           type=str, required=False,
                           help="Pickle file that contains peaklists from step 'process-scans' or 'replicate filter'")

    parser_as.add_argument('-o', '--pickle-file-out',
                           type=str, required=False,
                           help="Pickle file to store the peaklist objects")

    parser_as.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in parts per million to group peaks across scans / mass spectra")

    parser_as.add_argument('-b', '--block-size',
                           default=2000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_as.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units.")

    #################################
    # Blank Filter
    #################################

    parser_bf.add_argument('-i', '--pickle-file-in',
                           type=str, required=False,
                           help="Pickle file that contains a peak matrix from step 'align samples'")

    parser_bf.add_argument('-o', '--pickle-file-out',
                           type=str, required=False,
                           help="Pickle file to store the peak matrix object")

    parser_bf.add_argument('-l', '--blank-label',
                           default="blank", type=str, required=True,
                           help="Class label for blanks")

    parser_bf.add_argument('-m', '--min-fraction',
                           default=1.0, type=float, required=False,
                           help="Minium fold change blank verus sample.")

    parser_bf.add_argument('-f', '--function',
                           choices=["mean", "median", "max"], required=False,
                           help="Select function to calculate blank intenstiy")

    parser_bf.add_argument('-c', '--min-fold-change',
                           default=1.0, type=float, required=False,
                           help="Minium fold change blank verus sample.")
                           
    parser_bf.add_argument('-r', '--remove-blank-samples',
                           action='store_true', required=False,
                           help="Remove blank samples from peak matrix")

    #################################
    # Sample Filter
    #################################

    parser_sf.add_argument('-i', '--pickle-file-in',
                           type=str, required=False,
                           help="Pickle file that contains a peak matrix from step 'align samples'.")

    parser_sf.add_argument('-o', '--pickle-file-out',
                           type=str, required=True,
                           help="Pickle file that contains a peak matrix from step 'align samples'.")

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

    #################################
    # Pickle to readable
    #################################

    parser_pr.add_argument('-i', '--pickle-file-in',
                           type=str, required=True,
                           help="Pickle file that contains a list of peaklists or peak matrix from one of the processing steps.")

    parser_pr.add_argument('-o', '--path-out',
                           type=str, required=True,
                           help="Directory (peaklists) or output file (peak matrix).")

    parser_pr.add_argument('-s', '--separator',
                           default="tab", choices=["tab", "comma"],
                           help="Format of the file for further data processing or data analysis.")

    parser_pr.add_argument('-t', '--transpose',
                           action='store_true', required=False,
                           help="Transpose peak matrix")

    args = parser.parse_args()
    print args

    if args.step == "check-filelist":
        workflow.read_filelist(args.filelist, args.source)

    elif args.step == "process-scans":
        with open(args.pickle_file_out, "wb") as fn_pkl:
            peaklists = workflow.process_scans(source=args.source,
                filelist=args.filelist,
                fn_exp=args.filename_experiment,
                nscans=args.nscans,
                function_noise=args.function_noise,
                snr_thres=args.snr_threshold,
                ppm=args.ppm,
                min_fraction=args.min_fraction,
                rsd_thres=args.rsd_threshold,
                block_size=args.block_size,
                ncpus=args.ncpus)
            pickle.dump(peaklists, fn_pkl, -1)

    elif args.step == "mass-calibrate":
        # TODO implement mass calibration
        # TODO methods available in python
        print args.method

    elif args.step == "replicate-filter":
        with open(args.pickle_file_in, "rb") as fn_pkl_in:
            peaklists = pickle.load(fn_pkl_in)
            print peaklists
            with open(args.pickle_file_out, "wb") as fn_pkl_out:
                peaklists_rf = workflow.replicate_filter(peaklists, args.ppm, args.replicates, args.min_peaks, rsd_thres=args.rsd_threshold, block_size=args.block_size, ncpus=args.ncpus)
                pickle.dump(peaklists_rf, fn_pkl_out, -1)

    elif args.step == "align-samples":
        with open(args.pickle_file_in, "rb") as fn_pkl_in:
            with open(args.pickle_file_out, "wb") as fn_pkl_out:
                peaklists = pickle.load(fn_pkl_in)
                pm = workflow.align_samples(peaklists, args.ppm, block_size=args.block_size, ncpus=args.ncpus)
                pickle.dump(pm, fn_pkl_out, -1)

    elif args.step == "blank-filter":
        with open(args.pickle_file_in, "rb") as fn_pkl_in:
            with open(args.pickle_file_out, "wb") as fn_pkl_out:
                peak_matrix = pickle.load(fn_pkl_in)
                pm_bf = workflow.blank_filter(peak_matrix, blank_label=args.blank_label, min_fold_change=args.min_fold_change, function=args.function, min_fraction=args.min_fraction, rm_samples=args.remove_blank_samples)
                pickle.dump(pm_bf, fn_pkl_out, -1)

    elif args.step == "sample-filter":
        with open(args.pickle_file_in, "rb") as fn_pkl_in:
            with open(args.pickle_file_out, "wb") as fn_pkl_out:
                peak_matrix = pickle.load(fn_pkl_in)
                pm_sf = workflow.sample_filter(peak_matrix, min_fraction=args.min_fraction, within=args.within, rsd=args.rsd_threshold, qc_label=args.qc_label)
                pickle.dump(pm_sf, fn_pkl_out, -1)

    elif args.step == "pickle-to-readable":
        workflow.to_readable(args.pickle_file_in, path_out=args.path_out, separator=args.separator, transpose=args.transpose)
