#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Ralf Weber, Albert Zhou
origin: 05-23-2017

"""


import unittest
import numpy as np
from dimspy.tools import *
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
                            filter_scan_events=None, report=to_test_result("MTBLS79_mzml_single_report.txt"), block_size=5000, ncpus=None)

        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_single.hdf5"))  # creating test set
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_single.hdf5"))
        # with open(to_test_result("test_pm_mzml.txt"), "w") as out: out.write(pls[0].to_str("\t"))  # creating test set
        # with open(to_test_result("test_pm_mzml_comp.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))  # creating test set
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])
        with open(to_test_result("MTBLS79_mzml_single_report.txt"), "rU") as test_result:
            with open(to_test_data("MTBLS79_mzml_single_report.txt"), "rU") as comp:
                self.assertEqual(test_result.read(), comp.read())

        pls = process_scans(to_test_data("MTBLS79_mzml_single.zip"), function_noise="median",
                            snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=0.5, rsd_thres=30.0,
                            filelist=to_test_data("filelist_mzml_single.txt"),
                            filter_scan_events=None, report=None, block_size=5000, ncpus=None)

        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_frac_rsd.hdf5"))  # creating test set
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_frac_rsd.hdf5"))
        # with open(to_test_result("test_pl_frac_rsd.txt"), "w") as out: out.write(pls[0].to_str("\t"))  # creating test set
        # with open(to_test_result("test_pl_frac_rsd_comp.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))  # creating test set
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

    def test_process_scans_raw_path(self):
        pls = process_scans(to_test_data("raw", "batch04_QC17_rep01_262.RAW"), function_noise="noise_packets",
                            snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None, filelist=None,
                            filter_scan_events=None, report=None, block_size=5000, ncpus=None)
        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_raw_batch04_QC17_rep01_262.hdf5"))  # creating test set
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_raw_batch04_QC17_rep01_262.hdf5"))
        # with open(to_test_result("test_pl_raw.txt"), "w") as out: out.write(pls[0].to_str("\t"))  # creating test set
        # with open(to_test_result("test_pl_comp_raw.txt"), "w") as out: out.write(pls_comp[0].to_str("\t"))  # creating test set
        self.assertEqual([pl.to_str() for pl in pls], [pl.to_str() for pl in pls_comp])

    def test_replicate_filter(self):
        # pls = process_scans(to_test_data("MTBLS79_mzml_triplicates.zip"), function_noise="median",  # creating test set
        #                     snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,  # creating test set
        #                     filelist=to_test_data("filelist_mzml_triplicates.txt"),  # creating test set
        #                     filter_scan_events=None, report=None, block_size=5000, ncpus=None)  # creating test set
        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_triplicates.hdf5"))  # creating test set
        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None,
                                  filelist=None, report=to_test_result("MTBLS79_mzml_triplicates_report.txt"),
                                  block_size=5000, ncpus=None)

        # save_peaklists_as_hdf5(pls_rf, to_test_data("MTBLS79_mzml_triplicates_rf.hdf5"))  # creating test set
        pls_rf_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates_rf.hdf5"))
        # with open(to_test_result("test_pl_rf.txt"), "w") as out: out.write(pls_rf[0].to_str("\t"))  # creating test set
        # with open(to_test_result("test_pl_rf_comp.txt"), "w") as out: out.write(pls_rf_comp[0].to_str("\t"))  # creating test set
        self.assertEqual([pl.to_str() for pl in pls_rf], [pl.to_str() for pl in pls_rf_comp])
        with open(to_test_result("MTBLS79_mzml_triplicates_report.txt"), "rU") as test_result:
            with open(to_test_data("MTBLS79_mzml_triplicates_report.txt"), "rU") as comp:
                self.assertEqual(test_result.read(), comp.read())

        pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None,
                                  filelist=to_test_data("filelist_mzml_triplicates.txt"),
                                  report=None, block_size=5000, ncpus=None)
        # save_peaklists_as_hdf5(pls_rf, to_test_data("MTBLS79_mzml_triplicates_rf_2.hdf5"))
        pls_rf_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates_rf_2.hdf5"))
        # with open(to_test_result("test_pl_rf_2.txt"), "w") as out: out.write(pls_rf[0].to_str("\t"))  # creating test set
        # with open(to_test_result("test_pl_rf_2_comp.txt"), "w") as out: out.write(pls_rf_comp[0].to_str("\t"))  # creating test set
        self.assertEqual([pl.to_str() for pl in pls_rf], [pl.to_str() for pl in pls_rf_comp])

    def test_align_samples(self):
        # pls = process_scans(to_test_data("MTBLS79_mzml_triplicates.zip"), function_noise="median",  # creating test set
        #                     snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,  # creating test set
        #                     filelist=to_test_data("filelist_mzml_triplicates.txt"),  # creating test set
        #                     filter_scan_events=None, block_size=5000, ncpus=None)  # creating test set
        # save_peaklists_as_hdf5(pls, to_test_data("MTBLS79_mzml_triplicates.hdf5"))  # creating test set
        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        pm = align_samples(pls, ppm=2.0, filelist=None, block_size=5000, ncpus=None)
        # save_peak_matrix_as_hdf5(pm, to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))  # creating test set
        pm_comp = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        # with open(to_test_result("test_pm_as.txt"), "w") as out: out.write(pm.to_str())  # creating test set
        # with open(to_test_result("test_pm_as_comp.txt"), "w") as out: out.write(pm_comp.to_str())  # creating test set
        self.assertEqual(pm.to_str(), pm_comp.to_str())
        self.assertEqual(pm.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates_rf.hdf5"))
        pm = align_samples(pls, ppm=2.0, filelist=None, block_size=5000, ncpus=None)
        # save_peak_matrix_as_hdf5(pm, to_test_data("MTBLS79_mzml_peak_matrix_rf.hdf5"))  # creating test set
        pm_comp = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix_rf.hdf5"))
        # with open(to_test_result("test_pm_as_rf.txt"), "w") as out: out.write(pm.to_str())  # creating test set
        # with open(to_test_result("test_pm_as_rf_comp.txt"), "w") as out: out.write(pm_comp.to_str())  # creating test set
        self.assertEqual(pm.to_str(), pm_comp.to_str())
        self.assertEqual(pm.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())


    def test_blank_filter(self):
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        pm_bf = blank_filter(pm, blank_label="blank", min_fraction=1.0, min_fold_change=1.0, function="mean",
                             rm_samples=True, labels=None)
        # save_peak_matrix_as_hdf5(pm_bf, to_test_data("MTBLS79_mzml_peak_matrix_bf.hdf5"))  # creating test set
        pm_comp = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix_bf.hdf5"))
        # with open(to_test_result("test_pm_bf.txt"), "w") as out: out.write(pm.to_str())  # creating test set
        # with open(to_test_result("test_pm_bf_comp.txt"), "w") as out: out.write(pm_comp.to_str())  # creating test set
        self.assertEqual(pm_bf.to_str(), pm_comp.to_str())
        self.assertEqual(pm_bf.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

    def test_sample_filter(self):
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix_bf.hdf5"))
        pm_sf = sample_filter(pm, min_fraction=0.8, within=False, rsd=None, qc_label=None, labels=None)
        # save_peak_matrix_as_hdf5(pm_sf, to_test_data("MTBLS79_mzml_peak_matrix_sf.hdf5"))  # creating test set
        pm_comp = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix_sf.hdf5"))
        # with open(to_test_result("test_pm_sf.txt"), "w") as out: out.write(pm.to_str())  # creating test set
        # with open(to_test_result("test_pm_sf_comp.txt"), "w") as out: out.write(pm_comp.to_str())  # creating test set
        self.assertEqual(pm_sf.to_str(), pm_comp.to_str())
        self.assertEqual(pm_sf.to_peaklist("pl").to_str(), pm_comp.to_peaklist("pl").to_str())

    def test_remove_samples(self):
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        pm_rs = remove_samples(pm, ["batch04_B02_rep01_301.mzML", "batch04_B02_rep02_302.mzML", "batch04_B02_rep03_303.mzML"])
        self.assertEqual(pm_rs.shape, (6, 2692))

        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates_rf.hdf5"))
        pls_rs = remove_samples(pls, ["batch04_B02_rep01_301_2_302_3_303", "batch04_QC17_rep01_262_2_263_3_264"])
        self.assertEqual(len(pls_rs), 1)

    def test_missing_values_sample_filter(self):
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        pm_mv_sf = missing_values_sample_filter(pm, 1.0)
        self.assertEqual(pm_mv_sf.shape, (9, 4617))
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        pm_mv_sf = missing_values_sample_filter(pm, 0.6)
        self.assertEqual(pm_mv_sf.shape, (3, 2595))
        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        pm_mv_sf = missing_values_sample_filter(pm, 0.3)
        self.assertEqual(pm_mv_sf.shape, (0, 0))

    def test_merge_peaklists(self):
        # pls_rep_01 = process_scans(to_test_data("raw", "batch04_QC17_rep01_262.RAW"),  # creating test set
        #                            function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,  # creating test set
        #                            min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,  # creating test set
        #                            block_size=5000, ncpus=None)  # creating test set
        #
        # pls_rep_02 = process_scans(to_test_data("raw", "batch04_QC17_rep02_263.RAW"),  # creating test set
        #                            function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,  # creating test set
        #                            min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,  # creating test set
        #                            block_size=5000, ncpus=None)  # creating test set
        # pls_rep_03 = process_scans(to_test_data("raw", "batch04_QC17_rep03_264.RAW"),  # creating test set
        #                            function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,  # creating test set
        #                            min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,  # creating test set
        #                            block_size=5000, ncpus=None)  # creating test set
        # save_peaklists_as_hdf5(pls_rep_01, to_test_data("batch04_QC17_rep01_262.hdf5"))  # creating test set
        # save_peaklists_as_hdf5(pls_rep_02, to_test_data("batch04_QC17_rep02_263.hdf5"))  # creating test set
        # save_peaklists_as_hdf5(pls_rep_03, to_test_data("batch04_QC17_rep03_264.hdf5"))  # creating test set

        pls_rep_01 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep01_262.hdf5"))
        pls_rep_02 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep02_263.hdf5"))
        pls_rep_03 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep03_264.hdf5"))

        pls_merged = merge_peaklists([pls_rep_01, pls_rep_02, pls_rep_03])
        # save_peaklists_as_hdf5(pls_merged, to_test_data("MTBLS79__01_262_02_263.hdf5"))  # creating test set
        pls_comp = load_peaklists_from_hdf5(to_test_data("MTBLS79__01_262_02_263.hdf5"))
        self.assertEqual([pl.to_str() for pl in pls_merged], [pl.to_str() for pl in pls_comp])

    def test_merge_peaklist_multilist(self):

        pls_01 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep01_262.hdf5"))
        pls_02 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep02_263.hdf5"))
        pls_03 = load_peaklists_from_hdf5(to_test_data("batch04_QC17_rep03_264.hdf5"))

        merged_peaklists = merge_peaklists([pls_01, pls_02, pls_03], to_test_data("filelist_multi.txt"))

        # for i in range(len(merged_peaklists)):  # creating test set
        #    hdf5_portal.save_peaklists_as_hdf5(merged_peaklists[i], to_test_data('merged_{}.hdf5'.format(i)))  # creating test set

        merged_comp_01 = load_peaklists_from_hdf5(to_test_data("merged_0.hdf5"))
        merged_comp_02 = load_peaklists_from_hdf5(to_test_data("merged_1.hdf5"))

        self.assertEqual(len([pl.to_str() for pl in merged_peaklists[0]]), 2)
        self.assertEqual(len([pl.to_str() for pl in merged_peaklists[1]]), 1)
        self.assertEqual([pl.to_str() for pl in merged_peaklists[0]], [pl.to_str() for pl in merged_comp_01])
        self.assertEqual([pl.to_str() for pl in merged_peaklists[1]], [pl.to_str() for pl in merged_comp_02])

    def test_hdf5_peaklists_to_txt(self):
        hdf5_peaklists_to_txt(to_test_data("MTBLS79_mzml_triplicates.hdf5"), to_test_result(""), delimiter="\t")
        for fn in ["batch04_QC17_rep01_262.txt", "batch04_QC17_rep02_263.txt", "batch04_QC17_rep03_264.txt"]:
            with open(to_test_result(fn), "rU") as test_result:
                with open(to_test_data(fn), "rU") as comp:
                    self.assertEqual(test_result.read(), comp.read())


    def test_hdf5_peak_matrix_to_txt(self):

        hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"), to_test_result("pm_mzml_triplicates.txt"),
                                attr_name="intensity", rsd_tags=(), delimiter="\t", samples_in_rows=True, comprehensive=False)
        with open(to_test_result("pm_mzml_triplicates.txt"), "rU") as test_result:
            ln = test_result.readline().split("\t")[:5]
            self.assertEqual(ln[0], "mz")
            self.assertTrue(np.allclose(map(float,ln[1:]), [74.0166655257, 74.0198337519, 74.0200238089, 74.0202012645], atol = 1e-10))

        hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"), to_test_result("pm_mzml_triplicates_comprehensive.txt"),
                                attr_name="intensity", rsd_tags=(Tag("QC", "classLabel"),), delimiter="\t", samples_in_rows=True, comprehensive=True)
        with open(to_test_result("pm_mzml_triplicates_comprehensive.txt"), "rU") as test_result:
            ln = test_result.readline().split("\t")[:8]
            self.assertEquals(ln[:-2], ['mz', 'missing values', 'tags_batch', 'tags_replicate', 'tags_injectionOrder', 'tags_classLabel'])
            self.assertTrue(np.isclose(float(ln[-1]), 74.0166655257))

        hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"), to_test_result("pm_mzml_triplicates_snr.txt"),
                                attr_name="snr", rsd_tags=("QC",), delimiter="\t", samples_in_rows=True, comprehensive=False)
        with open(to_test_result("pm_mzml_triplicates_snr.txt"), "rU") as test_result:
            ln = test_result.readlines()[1].split("\t")[:10]
            self.assertEqual(ln[0], "batch04_B02_rep01_301.mzML")
            self.assertTrue(np.allclose(map(float,ln[1:]), [0., 0., 0., 0., 0., 0., 3.60960180872, 0., 4.35180213987], atol = 1e-10))

        hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"), to_test_result("pm_mzml_triplicates_comprehensive_T.txt"),
                                attr_name="intensity", rsd_tags=(Tag("QC", "classLabel"),), delimiter="\t", samples_in_rows=False, comprehensive=True)
        with open(to_test_result("pm_mzml_triplicates_comprehensive_T.txt"), "rU") as test_result:
            self.assertEquals(test_result.readline().split("\t")[0:5],
                              ['mz', 'present', 'occurrence', 'purity', 'rsd_QC'])


    def test_create_sample_list(self):
        pls = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates.hdf5"))
        create_sample_list(pls, to_test_result("filelist_csl_MTBLS79_mzml_triplicates.txt"))
        with open(to_test_result("filelist_csl_MTBLS79_mzml_triplicates.txt"), "rU") as test_result:
            with open(to_test_data("filelist_csl_MTBLS79_mzml_triplicates.txt"), "rU") as comp:
                self.assertEqual(test_result.read(), comp.read())

        pm = load_peak_matrix_from_hdf5(to_test_data("MTBLS79_mzml_peak_matrix.hdf5"))
        create_sample_list(pm, to_test_result("filelist_csl_MTBLS79_mzml_peak_matrix.txt"))
        with open(to_test_result("filelist_csl_MTBLS79_mzml_peak_matrix.txt"), "rU") as test_result:
            with open(to_test_data("filelist_csl_MTBLS79_mzml_peak_matrix.txt"), "rU") as comp:
                self.assertEqual(test_result.read(), comp.read())

        pm = load_peaklists_from_hdf5(to_test_data("MTBLS79_mzml_triplicates_rf.hdf5"))
        create_sample_list(pm, to_test_result("filelist_csl_MTBLS79_mzml_triplicates_rf.txt"))
        with open(to_test_result("filelist_csl_MTBLS79_mzml_triplicates_rf.txt"), "rU") as test_result:
            with open(to_test_data("filelist_csl_MTBLS79_mzml_triplicates_rf.txt"), "rU") as comp:
                self.assertEqual(test_result.read(), comp.read())