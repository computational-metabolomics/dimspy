#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import workflow
import portals
from . import __version__


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
    parser_pr = subparsers.add_parser('pickle-to-readable', help='Write output to human-readable format.')

    parser_cf.add_argument('-l', '--filelist',
                           type=str, required=True,
                           help="Tab-delimited file that list all the MS data files (*.raw, *.mzml or tab-delimited files) and associated meta data (replicate, order, batch, class, qc, blank).")

    parser_cf.add_argument('-i', '--source',
                           type=str, required=True,
                           help="Directory (*.raw, *.mzml or tab-delimited files) or zip archive (*.mzml or tab-delimited files)")

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

    parser_ps.add_argument('-i', '--source',
                           type=str, required=True,
                           help="Directory (*.raw, *.mzml or tab-delimited files) or zip archive (*.mzml or tab-delimited files).")

    parser_ps.add_argument('-o', '--hdf5-file-out',
                           type=str, required=True,
                           help="HDF5 file to save the peaklist objects.")

    parser_ps.add_argument('-l', '--filelist',
                           type=str, required=False,
                           help="File (csv or tab-separated ) that list all the data files (*.raw or *.mzml) and meta data (filename, technical replicate, group, batch). "
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

    parser_ps.add_argument('-n', '--nscans',
                           default=0, type=int, required=False,
                           help="Maximum number of scans to select for each scan event. Use zero for all scans.")

    parser_ps.add_argument('-a', '--min-fraction',
                           default=1, type=float, required=False,
                           help="Minimum number of scans a peak has to be present.")

    parser_ps.add_argument('-d', '--rsd-threshold',
                           default=None, type=float, required=False,
                           help="Maximum threshold - relative standard deviation (Only applied to peaks that have been measured across a minimum of three scans).")

    parser_ps.add_argument('-e', '--subset-mz-ranges',
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

    parser_rf.add_argument('-i', '--hdf5-file-in',
                           type=str, required=False,
                           help="HDF5 file (Peaklist objects) or directory (tab-delimited files) that contains peaklists from step 'process-scans' or 'replicate filter'.")

    parser_rf.add_argument('-o', '--hdf5-file-out',
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

    parser_as.add_argument('-i', '--hdf5-file-in',
                           type=str, required=False,
                           help="HDF5 file (Peaklist objects) or directory (tab-delimited files) that contains peaklists from step 'process-scans' or 'replicate filter'.")

    parser_as.add_argument('-o', '--hdf5-file-out',
                           type=str, required=False,
                           help="HDF5 file to save the peak matrix object.")

    parser_as.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in parts per million to group peaks across scans / mass spectra.")

    parser_as.add_argument('-b', '--block-size',
                           default=2000, type=int, required=False,
                           help="The size of each block of peaks to perform clustering on.")

    parser_as.add_argument('-c', '--ncpus',
                           default=None, type=int, required=False,
                           help="Number of central processing units.")

    #################################
    # Blank Filter
    #################################

    parser_bf.add_argument('-i', '--file-in',
                           type=str, required=False,
                           help="HDF5 file or tab-delimited file that contains a peak matrix (object).")

    parser_bf.add_argument('-o', '--hdf5-file-out',
                           type=str, required=False,
                           help="HDF5 file to save the peak matrix object.")

    parser_bf.add_argument('-l', '--blank-label',
                           default="blank", type=str, required=True,
                           help="Class label for blanks.")

    parser_bf.add_argument('-m', '--min-fraction',
                           default=1.0, type=float, required=False,
                           help="Minium fold change blank versus sample.")

    parser_bf.add_argument('-f', '--function',
                           choices=["mean", "median", "max"], required=False,
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
    # Sample Filter
    #################################

    parser_sf.add_argument('-i', '--file-in',
                           type=str, required=False,
                           help="HDF5 file or tab-delimited file that contains a peak matrix (object).")

    parser_sf.add_argument('-o', '--hdf5-file-out',
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
    # HDF5 to text
    #################################

    parser_pr.add_argument('-i', '--hdf5-file-in',
                           type=str, required=True,
                           help="HDF5 file that contains peaklists or a peak matrix from one of the processing steps.")

    parser_pr.add_argument('-o', '--path-out',
                           type=str, required=True,
                           help="Directory (peaklists) or text file (peak matrix) to write to.")

    parser_pr.add_argument('-s', '--separator',
                           default="tab", choices=["tab", "comma"],
                           help="Format of the file for further data processing or data analysis.")

    parser_pr.add_argument('-t', '--transpose',
                           action='store_true', required=False,
                           help="Transpose peak matrix")

    args = parser.parse_args()
    print args

    if args.step == "process-scans":
        peaklists = workflow.process_scans(source=args.source,
            function_noise=args.function_noise,
            snr_thres=args.snr_threshold,
            nscans=args.nscans,
            ppm=args.ppm,
            min_fraction=args.min_fraction,
            rsd_thres=args.rsd_threshold,
            filelist=args.filelist,
            subset_mzrs=args.subset_mz_ranges,
            block_size=args.block_size,
            ncpus=args.ncpus)
        portals.save_peaklists_as_hdf5(peaklists, args.hdf5_file_out)

    elif args.step == "mass-calibrate":
        # TODO implement mass calibration
        # TODO methods available in python
        print args.method

    elif args.step == "replicate-filter":
        peaklists_rf = workflow.replicate_filter(args.source,
            ppm=args.ppm,
            reps=args.replicates,
            min_peaks=args.min_peaks,
            rsd_thres=args.rsd_threshold,
            filelist=args.filelist,
            block_size=args.block_size,
            ncpus=args.ncpus)
        portals.save_peaklists_as_hdf5(peaklists_rf, args.hdf5_file_out)

    elif args.step == "align-samples":
        pm = workflow.align_samples(args.source,
            ppm=args.ppm,
            filelist=args.filelist,
            block_size=args.block_size,
            ncpus=args.ncpus)
        portals.save_peak_matrix_as_hdf5(pm, args.hdf5_file_out)
    elif args.step == "blank-filter":
        pm_bf = workflow.blank_filter(args.source,
            blank_label=args.blank_label,
            min_fraction=args.min_fraction,
            min_fold_change=args.min_fold_change,
            function=args.function,
            rm_samples=args.remove_blank_samples,
            tsv_labels=args.class_labels)
        portals.save_peak_matrix_as_hdf5(pm_bf, args.hdf5_file_out)

    elif args.step == "sample-filter":
        pm_sf = workflow.sample_filter(args.source,
            min_fraction=args.min_fraction,
            within=args.within,
            rsd=args.rsd_threshold,
            qc_label=args.qc_label,
            tsv_labels=args.class_labels)
        portals.save_peak_matrix_as_hdf5(pm_sf, args.hdf5_file_out)

    elif args.step == "hdf5-to-text":
        portals.hdf5_to_text(args.hdf5_file_in,
            path_out=args.path_out,
            separator=args.separator,
            transpose=args.transpose)
