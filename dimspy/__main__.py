#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import argparse
import h5py
import tools
from portals import hdf5_portal
from . import __version__


def map_delimiter(delimiter):
    seps = {"comma": ",", "tab": "\t"}
    if delimiter in seps:
        return seps[delimiter]
    else:
        return delimiter


def main():

    print("Executing dimspy version %s." % __version__)

    parser = argparse.ArgumentParser(description='Python package to process DIMS data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')

    parser_ps = subparsers.add_parser('process-scans', help='Process scans and/or stitch SIM windows.')
    parser_rf = subparsers.add_parser('replicate-filter', help='Filter irreproducible peaks from technical replicate peaklists.')
    parser_as = subparsers.add_parser('align-samples', help='Align peaklists across samples.')
    parser_bf = subparsers.add_parser('blank-filter', help='Filter peaks across samples that are present in the blank samples.')
    parser_sf = subparsers.add_parser('sample-filter', help='Filter peaks based on certain reproducibility and sample class criteria.')
    parser_rs = subparsers.add_parser('remove-samples', help='Remove sample(s) from a peak matrix object or list of peaklist objects.')
    parser_mvsf = subparsers.add_parser('mv-sample-filter', help='Filter samples based on the percentage of missing values.')
    parser_mp = subparsers.add_parser('merge-peaklists', help='Merge peaklists from multiple lists of peaklist or peak matrix objects.')
    parser_gp = subparsers.add_parser('get-peaklists', help='Get peaklists from a peak matrix object.')
    parser_gap = subparsers.add_parser('get-average-peaklist', help='Get an average peaklist from a peak matrix object.')
    parser_hpmt = subparsers.add_parser('hdf5-pm-to-txt', help='Write HDF5 output (peak matrix) to text format.')
    parser_hplt = subparsers.add_parser('hdf5-pls-to-txt', help='Write HDF5 output (peak lists) to text format.')
    parser_csl = subparsers.add_parser('create-sample-list', help='Create a sample list from a peak matrix object or list of peaklist objects.')

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
                           default=3.0,  type=float,  required=True,
                           help="Signal-to-noise threshold")

    parser_ps.add_argument('-p', '--ppm',
                           default=2.0, type=float, required=False,
                           help="Mass tolerance in Parts per million to group peaks across scans / mass spectra.")

    parser_ps.add_argument('-n', '--min_scans',
                           default=1, type=int, required=False,
                           help="Minimum number of scans required for each m/z range or event (header).")

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
                           help="Scan events to select. E.g. 100.0 200.0 sim  or  50.0 1000.0 full")

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
                           default=None, type=float,  required=False,
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

    args = parser.parse_args()
    
    print args

    if args.step == "process-scans":

        filter_scan_events = {}
        if args.exclude_scan_events != [] and args.include_scan_events != []:
            raise argparse.ArgumentTypeError("-u/--include-scan-events and -x/--exclude-scan-events can not be used together.")
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
                                    rsd=args.rsd_threshold,
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
                if 'multilistLabel' in pls_merged[i][0].metadata.keys():
                    m_label = pls_merged[i][0].metadata['multilistLabel']
                    hdf5_portal.save_peaklists_as_hdf5(pls_merged[i],
                                                       os.path.join(args.output, '{}.hdf5'.format(m_label)))
                else:
                    m_nm = pls_merged[i][0].metadata['multilist']
                    hdf5_portal.save_peaklists_as_hdf5(pls_merged[i],
                                                       os.path.join(args.output,
                                                       'merged_peaklist_{:03}.hdf5'.format(m_nm)))
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

if __name__ == "__main__":
    main()

