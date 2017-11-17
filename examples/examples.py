#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from dimspy.tools import process_scans
from dimspy.tools import replicate_filter
from dimspy.tools import create_sample_list
from dimspy.tools import align_samples
from dimspy.tools import hdf5_peak_matrix_to_txt
from dimspy.tools import blank_filter
from dimspy.tools import sample_filter
from dimspy.portals.hdf5_portal import save_peaklists_as_hdf5, save_peak_matrix_as_hdf5
from dimspy.portals.hdf5_portal import load_peaklists_from_hdf5, load_peak_matrix_from_hdf5


def main():

    source = os.path.join("..", "tests", "data", "MTBLS79_subset", "MTBLS79_mzml_triplicates.zip")
    fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_mzml_triplicates.txt")
    output = os.path.join("..", "tests", "test_results")

    print "Process Scans....."
    pls = process_scans(source, min_scans=1, function_noise="median",
                        snr_thres=3.0, ppm=2.0, min_fraction=None, rsd_thres=None,
                        filelist=fn_filelist, remove_mz_range=[], block_size=5000, ncpus=None)
    print "Finished"
    print

    print "Replicate Filter....."
    logfile = os.path.join(output, "log_replicate_filter.txt")
    pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None, report=logfile, block_size=5000)
    print "Finished"
    print

    print "Create a new sample list"
    sample_list = os.path.join(output, "sample_list.txt")
    create_sample_list(pls_rf, sample_list, delimiter="\t")
    print "Finished"
    print

    print "Align Samples...."
    pm = align_samples(pls_rf, ppm=3.0, ncpus=1, block_size=5000)
    print "Finished", pm.shape
    print

    save_peak_matrix_as_hdf5(pm, os.path.join(output, "pm.h5"))
    hdf5_peak_matrix_to_txt(os.path.join(output, "pm.h5"), path_out=os.path.join(output, "pm.txt"), attr_name="intensity", comprehensive=True)

    pm = load_peak_matrix_from_hdf5(os.path.join(output, "pm.h5"))

    print "Blank Filter"
    pm_bf = blank_filter(pm, "blank", min_fraction=1.0, min_fold_change=10.0, function="mean", rm_samples=True)
    print "Finished", pm_bf.shape
    print

    print "Sample Filter"
    pm_bf_sf = sample_filter(pm, 0.8, within=False)
    print "Finished", pm_bf_sf.shape
    print

if __name__ == '__main__':
    main()
