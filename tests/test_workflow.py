#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Ralf Weber, Albert Zhou
origin: 05-23-2017

"""


import unittest
from dimspy.workflow import *
from dimspy.portals.hdf5_portal import *


class WorkflowTestCase(unittest.TestCase):

    path_test_data = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset")

    def test_process_scans_mzml_zip(self):
        pls = process_scans(os.path.join(self.path_test_data, "MTBLS79_mzml_single.zip"), function_noise="median",
                            snr_thres=0.5, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,
                            filelist=os.path.join(self.path_test_data, "filelist_mzml_single.txt"),
                            filter_scan_events=None, block_size=2000, ncpus=None)

        # save_peaklists_as_hdf5(pls, os.path.join(self.path_test_data, "MTBLS79_mzml_single.hdf5"))
        pls_comp = load_peaklists_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_single.hdf5"))
        self.assertEqual([pl.to_str()[0:1000] for pl in pls], [pl.to_str()[0:1000] for pl in pls_comp])
        #with open(os.path.join("test_pm_comp.txt"), "w") as out: out.write(pls[0].to_str("\t"))
        #with open(os.path.join("test_pm.txt.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

    def test_process_scans_raw_path(self):
        pls = process_scans(os.path.join(self.path_test_data, "raw", "batch04_QC17_rep01_262.RAW"), function_noise="noise_packets",
                            snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None, filelist=None,
                            filter_scan_events=None, block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls, os.path.join(self.path_test_data, "MTBLS79_raw_batch04_QC17_rep01_262.hdf5"))
        pls_comp = load_peaklists_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_raw_batch04_QC17_rep01_262.hdf5"))
        self.assertEqual([pl.to_str()[0:1000] for pl in pls], [pl.to_str()[0:1000] for pl in pls_comp])
        #with open(os.path.join("test_pm_comp.txt"), "w") as out: out.write(pls[0].to_str("\t"))
        #with open(os.path.join("test_pm.txt.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

    def test_replicate_filter(self):
        # pls = process_scans(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates.zip"), function_noise="median", snr_thres=3.0, min_scans=1,
        #                     ppm=2.0, min_fraction=None, rsd_thres=None,
        #                     filelist=os.path.join(self.path_test_data, "filelist_mzml_triplicates.txt"),
        #                     filter_scan_events=None, block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates.hdf5"))
        pls = load_peaklists_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates.hdf5"))
        pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None,
                                  filelist=None,
                                  block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls_rf, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates_rf.hdf5"))
        pls_rf_comp = load_peaklists_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates_rf.hdf5"))
        self.assertEqual([pl.to_str()[0:1000] for pl in pls_rf], [pl.to_str()[0:1000] for pl in pls_rf_comp])
        self.assertEqual([pl.to_str() for pl in pls_rf], [pl.to_str() for pl in pls_rf_comp])

        pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None,
                                  filelist=os.path.join(self.path_test_data, "filelist_mzml_triplicates.txt"),
                                  block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls_rf, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates_rf_2.hdf5"))
        pls_rf_comp = load_peaklists_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates_rf_2.hdf5"))
        self.assertEqual([pl.to_str()[0:1000] for pl in pls_rf], [pl.to_str()[0:1000] for pl in pls_rf_comp])
        self.assertEqual([pl.to_str() for pl in pls_rf], [pl.to_str() for pl in pls_rf_comp])

    def test_align_samples(self):
        # pls = process_scans(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates.zip"), "median",
        #                     3.0, 1, 2.0, min_fraction=None, rsd_thres=None,
        #                     filelist=os.path.join(self.path_test_data, "filelist_mzml_triplicates.txt"),
        #                     filter_scan_events=None, block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls, os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates.hdf5"))
        pls = load_peaklists_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_triplicates.hdf5"))
        pm = align_samples(pls, ppm=2.0, filelist=None, block_size=2000, ncpus=None)
        # save_peak_matrix_as_hdf5(pm, os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix.hdf5"))
        pm_comp = load_peak_matrix_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix.hdf5"))
        self.assertEqual(pm.to_str(), pm_comp.to_str())
        self.assertEqual(pm.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

    def test_blank_filter(self):
        pm = load_peak_matrix_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix.hdf5"))
        pm_bf = blank_filter(pm, blank_label="blank", min_fraction=1.0, min_fold_change=1.0, function="mean",
                             rm_samples=True, class_labels=None)
        # save_peak_matrix_as_hdf5(pm_bf, os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix_bf.hdf5"))
        pm_comp = load_peak_matrix_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix_bf.hdf5"))
        self.assertEqual(pm_bf.to_str(), pm_comp.to_str())
        self.assertEqual(pm_bf.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

    def test_sample_filter(self):
        pm = load_peak_matrix_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix_bf.hdf5"))
        pm_sf = sample_filter(pm, min_fraction=0.8, within=False, rsd=None, qc_label=None, class_labels=None)
        # save_peak_matrix_as_hdf5(pm_sf, os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix_sf.hdf5"))
        pm_comp = load_peak_matrix_from_hdf5(os.path.join(self.path_test_data, "MTBLS79_mzml_peak_matrix_sf.hdf5"))
        self.assertEqual(pm_sf.to_str(), pm_comp.to_str())
        self.assertEqual(pm_sf.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())
