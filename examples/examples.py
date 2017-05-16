#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from dimspy.workflow import process_scans
from dimspy.workflow import replicate_filter
from dimspy.workflow import align_samples
from dimspy.workflow import blank_filter
from dimspy.workflow import sample_filter


def main():

    # Example 1 - mzML files (zip file)
    source = os.path.join("..", "tests", "data", "MTBLS79_subset", "MTBLS79_subset.zip")
    fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_mzML_subset.txt")
    output = os.path.join("..", "tests", "data", "MTBLS79_subset", "output")

    print
    print "Process Scans....."
    pls = process_scans(source, nscans=3, function_noise="median",
        snr_thres=10.0, ppm=2.0, min_fraction=None, rsd_thres=None,
        filelist=fn_filelist, subset_scan_events=None, block_size=2000, ncpus=None)
    print "Finished"

    for pl in pls[0:6]:
        print pl.ID, pl.shape
        with open(os.path.join(output, pl.ID + ".txt"), "w") as out: out.write(pl.to_str("\t"))

    print
    print "Replicate Filter....."
    pls_rf = replicate_filter(pls, ppm=2.0, reps=3, min_peaks=2, rsd_thres=20.0)
    print "Finished"
    print
    print "Align Samples...."
    pm = align_samples(pls_rf, ppm=2.0)
    print "Finished", pm.shape
    print
    print "Blank Filter"
    pm_bf = blank_filter(pm, "blank", min_fraction=1.0, min_fold_change=10.0, function="mean", rm_samples=True)
    print "Finished", pm_bf.shape
    print
    print "Sample Filter"
    pm_bf_sf = sample_filter(pm_bf, 0.8, within=False)
    print "Finished", pm_bf_sf.shape


    # Example 2 - RAW files (Directory)
    source = os.path.join("..", "tests", "data", "MTBLS79_subset", "raw")
    fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_raw.txt")
    print
    print "Process Scans....."
    pls = process_scans(source, nscans=1, function_noise="noise_packets",
        snr_thres=3.0, ppm=2.0, filelist=fn_filelist)

    print
    print "Finished"
    print
    print "Replicate Filter....."
    pls_rf = replicate_filter(pls, ppm=2.0, reps=3, min_peaks=2, rsd_thres=None)
    print "Finished"
    print
    print "Align Samples...."
    pm = align_samples(pls_rf, ppm=2.0)
    print "Finished", pm.shape
    print
    print "Sample Filter"
    pm_sf = sample_filter(pm, 0.8, within=False)
    print "Finished", pm_sf.shape

    # Example 3 - Subset m/z ranges
    source = os.path.join("..", "tests", "data", "MTBLS79_subset", "raw")
    fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_raw.txt")
    fn_mz_ranges = os.path.join("..", "tests", "data", "MTBLS79_subset", "mz_ranges.txt")

    print
    print "Process Scans....."
    pls = process_scans(source, nscans=3, function_noise="noise_packets",
        snr_thres=3.0, ppm=2.0, filelist=fn_filelist, subset_scan_events=[[70.0, 170.0, "sim"]])
    print
    print
    print "Process Scans......"
    pls = process_scans(source, nscans=3, function_noise="noise_packets",
        snr_thres=3.0, ppm=2.0, filelist=fn_filelist, subset_scan_events=fn_mz_ranges)

    # Example 3 - Replicate filter using text files
    source = os.path.join('..', 'tests', 'data', "txt")  # , "peaklists_txt.zip")# "peaklists_txt")
    fn_filelist_txt = os.path.join("..", 'tests', 'data', 'filelist_txt_subset.txt')

    print "Replicate Filter...."
    pls_rf_txt = replicate_filter(source, ppm=2.0, reps=3, min_peaks=2, rsd_thres=None, filelist=fn_filelist_txt)
    print "Finished"
    print
    print "Align Samples...."
    pm = align_samples(pls_rf_txt, ppm=2.0)
    print "Finished", pm.shape

    print
    print "Align Samples...."
    pm = align_samples(source, ppm=2.0, filelist=fn_filelist_txt)
    print "Finished", pm.shape

if __name__ == '__main__':
    main()

