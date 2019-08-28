#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
test_Peaklist

author(s): Ralf Weber, Albert Zhou
origin: 05-23-2017

"""


import unittest
import copy
import numpy as np
from dimspy.tools import *
from dimspy.portals.hdf5_portal import *


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "MTBLS79_subset", *args)


def to_test_result(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)


class WorkflowTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        zip_ref = zipfile.ZipFile(to_test_data("MTBLS79_mzml_triplicates.zip"), 'r')
        zip_ref.extractall(to_test_result("zip_data"))
        zip_ref.close()

        zip_ref = zipfile.ZipFile(to_test_data("MTBLS79_mzml_single.zip"), 'r')
        zip_ref.extractall(to_test_result("zip_data"))
        zip_ref.close()

        cls.pls_master = process_scans(to_test_result("zip_data"), function_noise="median",
                                       snr_thres=10.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,
                                       filelist=to_test_data("filelist_mzml_triplicates.txt"),
                                       filter_scan_events=None, report=None, block_size=5000, ncpus=None)

        cls.pm_master = align_samples(cls.pls_master, ppm=2.0, filelist=None, block_size=5000, ncpus=None)

    def test_process_scans(self):

        pls = process_scans(to_test_result("zip_data"), function_noise="median",
                            snr_thres=3.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None,
                            filelist=to_test_data("filelist_mzml_single.txt"),
                            filter_scan_events=None, report=to_test_result("MTBLS79_mzml_single_report.txt"),
                            block_size=5000, ncpus=None)

        self.assertEqual(pls[0].mz[0], 74.0208368189612)
        self.assertEqual(pls[0].mz[-1], 571.3620355932388)
        self.assertEqual(list(pls[0].__dict__.keys()), ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls[0].intensity[0], 7441.7158203125)
        self.assertEqual(pls[0].intensity[-1], 3059.92041015625)
        self.assertEqual(len(pls[0]), 1800)
        self.assertEqual(pls[0].to_str()[0:37], "mz,intensity,snr,present,fraction,rsd")
        self.assertEqual(len(pls[0].to_str()), 172278)

        with open(to_test_data("MTBLS79_mzml_single_report.txt"), "r") as test_result:
            with open(to_test_data("MTBLS79_mzml_single_report.txt"), "r") as comp:
                self.assertEqual(test_result.read(), comp.read())

        pls = process_scans(to_test_result("zip_data"), function_noise="median",
                            snr_thres=10.0, min_scans=1, ppm=2.0, min_fraction=0.5, rsd_thres=30.0,
                            filelist=to_test_data("filelist_mzml_single.txt"),
                            filter_scan_events=None, report=None, block_size=5000, ncpus=None)

        self.assertEqual(pls[0].mz[0], 98.99550184665189)
        self.assertEqual(pls[0].mz[-1], 570.367357910245)
        self.assertEqual(list(pls[0].__dict__.keys()), ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls[0].intensity[0], 82339.65546875)
        self.assertEqual(pls[0].intensity[-1], 12712.233173076924)
        self.assertEqual(len(pls[0]), 369)
        self.assertEqual(pls[0].to_str()[0:37], "mz,intensity,snr,present,fraction,rsd")
        self.assertEqual(len(pls[0].to_str()), 48888)

        pls = process_scans(to_test_data("raw", "batch04_QC17_rep01_262.RAW"), function_noise="noise_packets",
                            snr_thres=10.0, min_scans=1, ppm=2.0, min_fraction=None, rsd_thres=None, filelist=None,
                            filter_scan_events=None, report=None, block_size=5000, ncpus=None)

        self.assertEqual(pls[0].mz[0], 97.00499917179926)
        self.assertEqual(pls[0].mz[-1], 570.367357910245)
        self.assertEqual(list(pls[0].__dict__.keys()), ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls[0].intensity[0], 20609.436279296875)
        self.assertEqual(pls[0].intensity[-1], 12712.233173076924)
        self.assertEqual(len(pls[0]), 982)
        self.assertEqual(pls[0].to_str()[0:37], "mz,intensity,snr,present,fraction,rsd")
        self.assertEqual(len(pls[0].to_str()), 94221)

    def test_replicate_filter(self):

        pls_rf = replicate_filter(self.pls_master, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None,
                                  filelist=None, report=to_test_result("MTBLS79_mzml_triplicates_report.txt"),
                                  block_size=5000, ncpus=None)

        with open(to_test_result("MTBLS79_mzml_triplicates_report.txt"), "r") as test_result:
            with open(to_test_data("MTBLS79_mzml_triplicates_report.txt"), "r") as comp:
                self.assertEqual(test_result.read(), comp.read())

        self.assertEqual(pls_rf[0].mz[0], 95.08548671366543)
        self.assertEqual(pls_rf[0].mz[-1], 564.5550854783114)
        self.assertEqual(list(pls_rf[0].__dict__.keys()), ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls_rf[0].intensity[0], 8190.606998697917)
        self.assertEqual(pls_rf[0].intensity[-1], 38388.22486394899)
        self.assertEqual(len(pls_rf[0]), 650)
        self.assertEqual(pls_rf[0].to_str()[0:37], "mz,intensity,snr,present,fraction,rsd")
        self.assertEqual(len(pls_rf[0].to_str()), 72877)

        pls_rf = replicate_filter(self.pls_master, ppm=2.0, replicates=3, min_peaks=3, rsd_thres=None,
                                  filelist=to_test_data("filelist_mzml_triplicates.txt"),
                                  report=None, block_size=5000, ncpus=None)

        self.assertEqual(pls_rf[0].mz[0], 95.08548671366543)
        self.assertEqual(pls_rf[0].mz[-1], 564.5550854783114)
        self.assertEqual(list(pls_rf[0].__dict__.keys()), ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls_rf[0].intensity[0], 8190.606998697917)
        self.assertEqual(pls_rf[0].intensity[-1], 38388.22486394899)
        self.assertEqual(len(pls_rf[0]), 527)
        self.assertEqual(pls_rf[0].to_str()[0:37], "mz,intensity,snr,present,fraction,rsd")
        self.assertEqual(len(pls_rf[0].to_str()), 72877)

        pls_rf = replicate_filter(self.pls_master, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=10.0,
                                  filelist=to_test_data("filelist_mzml_triplicates.txt"),
                                  report=None, block_size=5000, ncpus=None)

        self.assertEqual(pls_rf[0].mz[0], 95.08548671366543)
        self.assertEqual(pls_rf[0].mz[-1], 559.5163918132778)
        self.assertEqual(list(pls_rf[0].__dict__.keys()), ['_dtable', '_id', '_metadata', '_tags', '_flags', '_flag_attrs'])
        self.assertEqual(pls_rf[0].intensity[0], 8190.606998697917)
        self.assertEqual(pls_rf[0].intensity[-1], 33800.578125)
        self.assertEqual(len(pls_rf[0]), 317)
        self.assertEqual(pls_rf[0].to_str()[0:37], "mz,intensity,snr,present,fraction,rsd")
        self.assertEqual(len(pls_rf[0].to_str()), 74438)

    def test_align_samples(self):

        self.assertEqual(self.pm_master.shape, (9, 1383))
        self.assertListEqual(list(self.pm_master.attr_mean_vector('mz')[0:4]),
                             [95.08548671366543, 96.04772304614822, 97.00500386290739, 97.06475291143175])
        self.assertListEqual(list(self.pm_master.attr_mean_vector('mz')[-4:]),
                             [571.369498269902, 571.3707401959643, 572.3643110643018, 572.3655998454249])

        self.assertListEqual(list(self.pm_master.attr_mean_vector('intensity')[0:4]),
                             [8190.606998697917, 5827.436712032904, 24371.84765625, 5479.027180989583])
        self.assertListEqual(list(self.pm_master.attr_mean_vector('intensity')[-4:]),
                             [18668.55898980035, 10166.205078125, 13053.23828125, 13837.34375])

        pm = align_samples(self.pls_master, ppm=50.0, filelist=None, block_size=5000, ncpus=None)

        self.assertEqual(pm.shape, (9, 915))
        self.assertListEqual(list(pm.attr_mean_vector('mz')[0:4]),
                             [95.08548671366543, 96.04772304614822, 97.00500386290739, 97.06475291143175])
        self.assertListEqual(list(pm.attr_mean_vector('mz')[-4:]),
                             [569.3635826997123, 570.3674042061052, 571.3662372673122, 572.3649554548633])

        self.assertListEqual(list(self.pm_master.attr_mean_vector('intensity')[0:4]),
                             [8190.606998697917, 5827.436712032904, 24371.84765625, 5479.027180989583])
        self.assertListEqual(list(self.pm_master.attr_mean_vector('intensity')[-4:]),
                             [18668.55898980035, 10166.205078125, 13053.23828125, 13837.34375])

    def test_blank_and_sample_filter(self):

        self.assertEqual(self.pm_master.shape, (9, 1383))
        pm_sf = sample_filter(copy.deepcopy(self.pm_master), min_fraction=0.8, within=False, rsd=None, qc_label=None, labels=None)
        self.assertEqual(pm_sf.shape, (9, 60))

        self.assertEqual(self.pm_master.shape, (9, 1383))
        pm_bf = blank_filter(copy.deepcopy(self.pm_master), blank_label="blank", min_fraction=1.0, min_fold_change=1.0, function="mean",
                             rm_samples=True, labels=None)
        self.assertEqual(pm_bf.shape, (6, 638))

        pm_bf_sf = sample_filter(pm_bf, min_fraction=0.8, within=False, rsd=None, qc_label=None, labels=None)
        self.assertEqual(pm_bf_sf.shape, (6, 306))

    def test_remove_samples(self):
        self.assertEqual(self.pm_master.shape, (9, 1383))
        pm_rs = remove_samples(copy.deepcopy(self.pm_master), ["batch04_B02_rep01_301.mzML",
                                                               "batch04_B02_rep02_302.mzML",
                                                               "batch04_B02_rep03_303.mzML"])
        self.assertEqual(pm_rs.shape, (6, 769))

    def test_missing_values_sample_filter(self):

        self.assertEqual(self.pm_master.shape, (9, 1383))
        pm_mv_sf = missing_values_sample_filter(copy.deepcopy(self.pm_master), 1.0)
        self.assertEqual(pm_mv_sf.shape, (9, 1383))

        self.assertEqual(self.pm_master.shape, (9, 1383))
        pm_mv_sf = missing_values_sample_filter(copy.deepcopy(self.pm_master), 0.6)
        self.assertEqual(pm_mv_sf.shape, (3, 776))

        self.assertEqual(self.pm_master.shape, (9, 1383))
        pm_mv_sf = missing_values_sample_filter(copy.deepcopy(self.pm_master), 0.3)
        self.assertEqual(pm_mv_sf.shape, (0, 0))

    def test_merge_peaklists(self):
        pls_rep_01 = process_scans(to_test_data("raw", "batch04_QC17_rep01_262.RAW"),
                                   function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,
                                   min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,
                                   block_size=5000, ncpus=None)
        pls_rep_02 = process_scans(to_test_data("raw", "batch04_QC17_rep02_263.RAW"),
                                   function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,
                                   min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,
                                   block_size=5000, ncpus=None)
        pls_rep_03 = process_scans(to_test_data("raw", "batch04_QC17_rep03_264.RAW"),
                                   function_noise="noise_packets", snr_thres=10.0, min_scans=1, ppm=2.0,
                                   min_fraction=None, rsd_thres=None, filelist=None, filter_scan_events=None,
                                   block_size=5000, ncpus=None)

        pls_merged_01 = merge_peaklists([pls_rep_01, pls_rep_02, pls_rep_03])
        pls_merged_02 = merge_peaklists([pls_rep_01, pls_rep_02, pls_rep_03], to_test_data("filelist_multi.txt"))

        self.assertEqual(len(pls_merged_01), 3)
        self.assertEqual(len(pls_merged_02[0]), 2)
        self.assertListEqual([pls_merged_02[0][0].ID, pls_merged_02[0][1].ID], ["batch04_QC17_rep01_262.RAW",
                                                                                "batch04_QC17_rep02_263.RAW"])
        self.assertEqual(len(pls_merged_02[1]), 1)

    def test_hdf5_peaklists_to_txt(self):

        for compatibility in [False, True]:

            version = "v2"
            if compatibility: version = "v1"

            hdf5_peaklists_to_txt(to_test_data("MTBLS79_mzml_triplicates_{}.hdf5".format(version)), to_test_result(""),
                                  delimiter="\t", compatibility_mode=compatibility)
            with open(to_test_result("batch04_QC17_rep01_262.txt"), "r") as test_result:
                with open(to_test_data("batch04_QC17_rep01_262_{}.txt".format(version)), "r") as comp:
                    self.assertEqual(test_result.read(), comp.read())

    def test_hdf5_peak_matrix_to_txt(self):

        for compatibility in [False, True]:

            version = "v2"
            if compatibility: version = "v1"

            hdf5_peak_matrix_to_txt(to_test_data("MTBLS79_mzml_peak_matrix_{}.hdf5".format(version)),
                                    to_test_result("pm_mzml_triplicates_{}.txt".format(version)),
                                    attr_name="intensity", rsd_tags=(),
                                    delimiter="\t", samples_in_rows=True, comprehensive=False,
                                    compatibility_mode=compatibility)

            with open(to_test_result("pm_mzml_triplicates_{}.txt".format(version)), "r") as test_result:
                with open(to_test_data("pm_mzml_triplicates_{}.txt".format(version)), "r") as comp:
                    self.assertEqual(test_result.read(), comp.read())

    def test_hdf5_before_after(self):

        save_peaklists_as_hdf5(self.pls_master, to_test_result("MTBLS79_mzml_triplicates.hdf5"))
        pls = load_peaklists_from_hdf5(to_test_result("MTBLS79_mzml_triplicates.hdf5"))
        self.assertEqual(len(pls), len(self.pls_master))
        self.assertTrue(np.all(pls[0].mz == self.pls_master[0].mz))
        self.assertTrue(np.all(pls[0].intensity == self.pls_master[0].intensity))
        self.assertTrue(np.all(pls[0].snr == self.pls_master[0].snr))

        save_peak_matrix_as_hdf5(self.pm_master, to_test_result("MTBLS79_mzml_peak_matrix.hdf5"))
        pm = load_peak_matrix_from_hdf5(to_test_result("MTBLS79_mzml_peak_matrix.hdf5"))
        self.assertEqual(pm.shape, self.pm_master.shape)
        self.assertTrue(np.all(pm.attr_mean_vector('mz') == self.pm_master.attr_mean_vector('mz')))
        self.assertTrue(np.all(pm.attr_mean_vector('intensity') == self.pm_master.attr_mean_vector('intensity')))
        self.assertTrue(np.all(pm.attr_mean_vector('snr') == self.pm_master.attr_mean_vector('snr')))

    def test_create_sample_list(self):
        create_sample_list(self.pls_master, to_test_result("filelist_csl_MTBLS79_mzml_triplicates.txt"))
        with open(to_test_result("filelist_csl_MTBLS79_mzml_triplicates.txt"), "r") as test_result:
            with open(to_test_data("filelist_csl_MTBLS79_mzml_triplicates.txt"), "r") as comp:
                self.assertEqual(test_result.read(), comp.read())

        create_sample_list(self.pm_master, to_test_result("filelist_csl_MTBLS79_mzml_peak_matrix.txt"))
        with open(to_test_result("filelist_csl_MTBLS79_mzml_peak_matrix.txt"), "r") as test_result:
            with open(to_test_data("filelist_csl_MTBLS79_mzml_peak_matrix.txt"), "r") as comp:
                self.assertEqual(test_result.read(), comp.read())

    @classmethod
    def tearDownClass(cls):
        for fn in os.listdir(to_test_result("")):
            if os.path.isfile(to_test_result(fn)):
                os.remove(to_test_result(fn))

        if os.path.isdir(to_test_result("zip_data")):
            for fn in os.listdir(to_test_result("zip_data")):
                if os.path.isfile(to_test_result(os.path.join("zip_data", fn))):
                    os.remove(to_test_result(os.path.join("zip_data", fn)))
            os.rmdir(to_test_result("zip_data"))


if __name__ == '__main__':
    unittest.main()
