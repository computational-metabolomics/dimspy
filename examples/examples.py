#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from dimspy.tools import process_scans
from dimspy.tools import replicate_filter
from dimspy.tools import create_sample_list
from dimspy.tools import align_samples
from dimspy.tools import blank_filter
from dimspy.tools import sample_filter
from dimspy.portals.hdf5_portal import save_peaklists_as_hdf5, save_peak_matrix_as_hdf5
from dimspy.portals.hdf5_portal import load_peaklists_from_hdf5, load_peak_matrix_from_hdf5


def main():

    # Example 1 - mzML files (zip file)
    source = os.path.join("Y:\users\zhangcy\polar_positive_data")
    fn_filelist = os.path.join("Y:\users\zhangcy\polar_positive_data\pos_filelist.txt")
    output = os.path.join("Y:\\users\\zhangcy\\polar_positive_data_results\\")

    source = os.path.join("Y:\users\zhangcy\polar_positive_data")
    fn_filelist = os.path.join("E:\\raw_109\pos_filelist.txt")
    output = os.path.join("E:\\raw_109")

    print "Process Scans....."
    pls = process_scans(source, min_scans=1, function_noise="noise_packets",
                        snr_thres=3.0, ppm=2.0, min_fraction=None, rsd_thres=None,
                        filelist=fn_filelist, remove_mz_range=[], filter_scan_events={"exclude":[[50.0, 620.0, "full"]]}, block_size=2000, ncpus=None)
    print "Finished"

    sample_list = os.path.join("E:\\raw_109\\sample_list.txt")
    create_sample_list(pls, sample_list, delimiter="\t")

    print
    print pls
    for pl in pls:
        print pl.ID, pl.shape
        with open(os.path.join(output, pl.ID + ".txt"), "w") as out: out.write(pl.to_str("\t"))

    """
    save_peaklists_as_hdf5(pls, os.path.join(output, "pls.h5"))
    pls = load_peaklists_from_hdf5(os.path.join(output, "pls.h5"))
    """

    print
    print "Replicate Filter....."
    logfile = os.path.join("E:\\raw_109\\log_replicate_filter.txt")
    pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None, quality_logfile=logfile)
    print "Finished"
    print

    sample_list = os.path.join("E:\\raw_109\\sample_list.txt")
    create_sample_list(pls_rf, sample_list, delimiter="\t")
    """
    save_peaklists_as_hdf5(pls, os.path.join(output, "pls_rf.h5"))

    print "Align Samples...."
    pm = align_samples(pls_rf, ppm=3.0)
    print "Finished", pm.shape
    print

    save_peak_matrix_as_hdf5(pm, os.path.join(output, "pm.h5"))

    print "Blank Filter"
    pm_bf = blank_filter(pm, "blank", min_fraction=1.0, min_fold_change=10.0, function="mean", rm_samples=True)
    print "Finished", pm_bf.shape
    print

    save_peak_matrix_as_hdf5(pm_bf, os.path.join(output, "pm_bf.h5"))

    print "Sample Filter"
    pm_bf_sf = sample_filter(pm, 0.8, within=False)
    print "Finished", pm_bf_sf.shape
    print

    save_peak_matrix_as_hdf5(pm_bf_sf, os.path.join(output, "pm_bf_sf.h5"))
    """

if __name__ == '__main__':
    main()
