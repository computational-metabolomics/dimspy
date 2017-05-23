#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Ralf Weber, Albert Zhou
origin: 05-23-2017

"""


import unittest
import numpy as np
from dimspy.models.peaklist import PeakList
from dimspy.process.peak_alignment import align_peaks
from dimspy.workflow import *
from dimspy.portals.hdf5_portal import save_peaklists_as_hdf5

class WorkflowTestCase(unittest.TestCase):

    def test_process_scans_mzml_zip(self):
        pls = process_scans(os.path.join("data", "MTBLS79_subset", "MTBLS79_mzml_subset.zip"), "median", 3.0, 1, 2.0, min_fraction=None, rsd_thres=None, filelist=None,
                      subset_scan_events=None, block_size=2000, ncpus=None)
        #pls_comp = load_peaklists(os.path.join("data", "MTBLS79_subset", "MTBLS79_mzml_subset.hdf5"))
        save_peaklists_as_hdf5(pls, os.path.join("data", "MTBLS79_subset", "MTBLS79_mzml_subset.hdf5"))

    def test_process_scans_raw_path(self):
        pls = process_scans(os.path.join("data", "MTBLS79_subset", "raw"), "median", 3.0, 1, 2.0, min_fraction=None, rsd_thres=None, filelist=None,
                      subset_scan_events=None, block_size=2000, ncpus=None)
        #pls_comp = load_peaklists(os.path.join("data", "MTBLS79_raw_subset.hdf5"))
        save_peaklists_as_hdf5(pls, os.path.join("data", "MTBLS79_subset", "MTBLS79_raw_subset.hdf5"))

"""
    def test_replicate_filter(self):
        replicate_filter(source, ppm, reps, min_peaks, rsd_thres=None, filelist=None, block_size=2000, ncpus=None)

    def test_align_samples(self):
        align_samples(source, ppm, filelist=None, block_size=2000, ncpus=None)

    def test_blank_filter(self):
        blank_filter(peak_matrix, blank_label, min_fraction=1.0, min_fold_change=1.0, function="mean", rm_samples=True,
                     class_labels=None)
    def test_sample_filter(self):
        sample_filter(peak_matrix, min_fraction, within=False, rsd=None, qc_label=None, class_labels=None)

    def test_load_peaklists(self):
        load_peaklists(source)

    def test_hdf5_to_txt(self):
        hdf5_to_txt(fname, path_out, separator="\t", transpose=False, extend=False)
"""