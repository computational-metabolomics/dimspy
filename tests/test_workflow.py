#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Ralf Weber, Albert Zhou
origin: 05-23-2017

"""


import unittest
import numpy as np
from dimspy.workflow import *
from dimspy.portals.hdf5_portal import save_peaklists_as_hdf5
from dimspy.portals.hdf5_portal import load_peaklists_from_hdf5


class WorkflowTestCase(unittest.TestCase):
    
    path_test_data = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset")

    def test_process_scans_mzml_zip(self):
        pls = process_scans(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate.zip"), "median", 3.0, 1, 2.0, min_fraction=None, rsd_thres=None, filelist=None, subset_scan_events=None, block_size=2000, ncpus=None)
        save_peaklists_as_hdf5(pls, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_single_scan.hdf5"))
        pls_comp = load_peaklists_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_single_scan.hdf5"))
        print pls[0].mz
        print pls_comp[0].mz
        self.assertTrue(np.allclose([pl.mz for pl in pls][0], [pl.mz for pl in pls_comp][0]))
        #self.assertEqual([pl.mz for pl in pls][0], [pl.mz for pl in pls_comp][0])

"""
        pls = process_scans(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate.zip"), "median", 3.0, 3, 2.0, min_fraction=None, rsd_thres=None, filelist=None, subset_scan_events=None, block_size=2000, ncpus=None)
        save_peaklists_as_hdf5(pls, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_multiple_scans.hdf5"))
        pls_comp = load_peaklists(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_multiple_scans.hdf5"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

    def test_process_scans_raw_path(self):
        pls = process_scans(os.path.join(self.path_test_data, "raw"), "median", 3.0, 1, 2.0, min_fraction=None, rsd_thres=None, filelist=None, subset_scan_events=None, block_size=2000, ncpus=None)
        save_peaklists_as_hdf5(pls, os.path.join(self.path_test_data, "MTBLS79_raw_triplicate_single_scan.hdf5"))
        pls_comp = load_peaklists(os.path.join(self.path_test_data, "MTBLS79_raw_triplicate_single_scan.hdf5"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

    def test_replicate_filter(self):
        pls = load_peaklists(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_single_scan.hdf5"))
        pls_rf = replicate_filter(pls, 2.0, 3, 2, rsd_thres=None, filelist=os.path.join(self.path_test_data, "filelist_mzml.txt"), block_size=2000, ncpus=None)
        save_peaklists_as_hdf5(pls_rf, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_single_scan_rf.hdf5"))
        pls_comp = load_peaklists(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_single_scan_rf.hdf5"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])
        
        pls = load_peaklists(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_multiple_scans.hdf5"))
        pls_rf = replicate_filter(pls, 2.0, 3, 2, rsd_thres=None, filelist=os.path.join(self.path_test_data, "filelist_mzml.txt"), block_size=2000, ncpus=None)
        save_peaklists_as_hdf5(pls_rf, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_multiple_scans_rf.hdf5"))
        pls_comp = load_peaklists(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicate_multiple_scans_rf.hdf5"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])


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
