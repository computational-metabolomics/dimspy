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


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset", *args)


def to_test_result(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)


class WorkflowTestCase(unittest.TestCase):

    def test_process_scans_mzml_zip(self):
        pls = process_scans(to_test_data("MTBLS79_mzml_single.zip"), function_noise="median",
                            snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,
                            filelist=to_test_data("filelist_mzml_single.txt"),
                            filter_scan_events=None, block_size=2000, ncpus=None)

        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_single.hdf5"))
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_single.hdf5"))
        # with open(to_test_result("test_pm_mzml.txt"), "w") as out: out.write(pls[0].to_str("\t"))
        # with open(to_test_result("test_pm_mzml_comp.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

        pls = process_scans(to_test_data("MTBLS79_mzml_single.zip"), function_noise="median",
                            snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=0.5, rsd_thres=30.0,
                            filelist=to_test_data("filelist_mzml_single.txt"),
                            filter_scan_events=None, block_size=2000, ncpus=None)

        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_frac_rsd.hdf5"))
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_frac_rsd.hdf5"))
        # with open(to_test_result("test_pl_frac_rsd.txt"), "w") as out: out.write(pls[0].to_str("\t"))
        # with open(to_test_result("test_pl_frac_rsd_comp.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])


    def test_process_scans_raw_path(self):
        pls = process_scans(to_test_data("raw", "batch04_QC17_rep01_262.RAW"), function_noise="noise_packets",
                            snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None, filelist=None,
                            filter_scan_events=None, block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_raw_batch04_QC17_rep01_262.hdf5"))
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_raw_batch04_QC17_rep01_262.hdf5"))
        # with open(to_test_result("test_pl_raw.txt"), "w") as out: out.write(pls[0].to_str("\t"))
        # with open(to_test_result("test_pl_comp_raw.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

    def test_replicate_filter(self):
        # pls = process_scans(to_test_data("MTBLS79_mzml_triplicates.zip"), function_noise="median",
        #                     snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,
        #                     filelist=to_test_data("filelist_mzml_triplicates.txt"),
        #                     filter_scan_events=None, block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None,
                                  filelist=None,
                                  block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls_rf, to_test_data("MTBLS79_mzml_triplicates_rf.hdf5"))
        pls_rf_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates_rf.hdf5"))
        # with open(to_test_result("test_pl_rf.txt"), "w") as out: out.write(pls_rf[0].to_str("\t"))
        # with open(to_test_result("test_pl_rf_comp.txt"), "w") as out: out.write(pls_rf_comp[0].to_str("\t"))
        self.assertEqual([pl.to_str() for pl in pls_rf], [pl.to_str() for pl in pls_rf_comp])

        pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None,
                                  filelist=to_test_data("filelist_mzml_triplicates.txt"),
                                  block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls_rf, to_test_data("MTBLS79_mzml_triplicates_rf_2.hdf5"))
        pls_rf_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates_rf_2.hdf5"))
        # with open(to_test_result("test_pl_rf_2.txt"), "w") as out: out.write(pls_rf[0].to_str("\t"))
        # with open(to_test_result("test_pl_rf_2_comp.txt"), "w") as out: out.write(pls_rf_comp[0].to_str("\t"))
        self.assertEqual([pl.to_str() for pl in pls_rf], [pl.to_str() for pl in pls_rf_comp])

    def test_align_samples(self):
        # pls = process_scans(to_test_data("MTBLS79_mzml_triplicates.zip"), function_noise="median",
        #                     snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,
        #                     filelist=to_test_data("filelist_mzml_triplicates.txt"),
        #                     filter_scan_events=None, block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        pm = align_samples(pls, ppm=2.0, filelist=None, block_size=2000, ncpus=None)
        # save_peak_matrix_as_hdf5(pm, to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        pm_comp = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        # with open(to_test_result("test_pm_as.txt"), "w") as out: out.write(pm.to_str())
        # with open(to_test_result("test_pm_as_comp.txt"), "w") as out: out.write(pm_comp.to_str())
        self.assertEqual(pm.to_str(), pm_comp.to_str())
        self.assertEqual(pm.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

    def test_blank_filter(self):
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        pm_bf = blank_filter(pm, blank_label="blank", min_fraction=1.0, min_fold_change=1.0, function="mean",
                             rm_samples=True, class_labels=None)
        # save_peak_matrix_as_hdf5(pm_bf, to_test_data("MTBLS79_mzml_peak_matrix_bf.hdf5"))
        pm_comp = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix_bf.hdf5"))
        # with open(to_test_result("test_pm_bf.txt"), "w") as out: out.write(pm.to_str())
        # with open(to_test_result("test_pm_bf_comp.txt"), "w") as out: out.write(pm_comp.to_str())
        self.assertEqual(pm_bf.to_str(), pm_comp.to_str())
        self.assertEqual(pm_bf.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

    def test_sample_filter(self):
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix_bf.hdf5"))
        pm_sf = sample_filter(pm, min_fraction=0.8, within=False, rsd=None, qc_label=None, class_labels=None)
        # save_peak_matrix_as_hdf5(pm_sf, to_test_data("MTBLS79_mzml_peak_matrix_sf.hdf5"))
        pm_comp = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix_sf.hdf5"))
        # with open(to_test_result("test_pm_sf.txt"), "w") as out: out.write(pm.to_str())
        # with open(to_test_result("test_pm_sf_comp.txt"), "w") as out: out.write(pm_comp.to_str())
        self.assertEqual(pm_sf.to_str(), pm_comp.to_str())
        self.assertEqual(pm_sf.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

    def test_merge_peaklists(self):
        # pls_rep_01 = process_scans(to_test_data("raw", "batch04_QC17_rep01_262.RAW"),
        #                            function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,
        #                            min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,
        #                            block_size=2000, ncpus=None)
        #
        # pls_rep_02 = process_scans(to_test_data("raw", "batch04_QC17_rep02_263.RAW"),
        #                            function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,
        #                            min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,
        #                            block_size=2000, ncpus=None)
        # save_peaklists_as_hdf5(pls_rep_01, to_test_data("batch04_QC17_rep01_262.hdf5"))
        # save_peaklists_as_hdf5(pls_rep_02, to_test_data("batch04_QC17_rep02_263.hdf5"))
        pls_rep_01 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep01_262.hdf5"))
        pls_rep_02 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep02_263.hdf5"))
        pls_merged = merge_peaklists([pls_rep_01, pls_rep_02])
        # save_peaklists_as_hdf5(pls_merged, to_test_data("MTBLS79__01_262_02_263.hdf5"))
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79__01_262_02_263.hdf5"))
        self.assertEqual([pl.to_str() for pl in pls_merged], [pl.to_str() for pl in pls_comp])

    def test_hdf5_peaklists_to_txt(self):
        hdf5_peaklists_to_txt(to_test_data("MTBLS79_mzml_triplicates.hdf5"), to_test_result(""), delimiter="\t")
        for fn in ["batch04_QC17_rep01_262.txt", "batch04_QC17_rep02_263.txt", "batch04_QC17_rep03_264.txt"]:
            with open(to_test_result(fn), "r") as test_result:
                with open(to_test_data(fn), "r") as comp:
                    self.assertEqual(test_result.read(), comp.read())

    def test_hdf5_peak_matrix_to_txt(self):
        hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"), to_test_result("pm_mzml_triplicates_comprehensive.txt"),
                                attr_name="intensity", rsd_tags=(), delimiter="\t", transpose=False, comprehensive=False)
        hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"), to_test_result("pm_mzml_triplicates.txt"),
                                attr_name="intensity", rsd_tags=("QC",), delimiter="\t", transpose=False, comprehensive=True)
        hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"), to_test_result("pm_mzml_triplicates_T.txt"),
                                attr_name="intensity", rsd_tags=("QC",), delimiter="\t", transpose=True, comprehensive=True)
        for fn in ["pm_mzml_triplicates.txt", "pm_mzml_triplicates_comprehensive.txt", "pm_mzml_triplicates_T.txt"]:
            with open(to_test_result(fn), "r") as test_result:
                with open(to_test_data(fn), "r") as comp:
                    self.assertEqual(test_result.read(), comp.read())

    def test_create_sample_list(self):
        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        create_sample_list(pls, to_test_result("filelist_mzml_triplicates_test_result.txt"))
        with open(to_test_result("filelist_mzml_triplicates_test_result.txt"), "r") as test_result:
            with open(to_test_data("filelist_mzml_triplicates_comp.txt"), "r") as comp:
                self.assertEqual(test_result.read(), comp.read())
