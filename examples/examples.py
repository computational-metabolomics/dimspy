#!/usr/bin/python
# -*- coding: utf-8 -*-

from dimspy.tools import *
from dimspy.portals.hdf5_portal import *
import zipfile


def main():

    source = os.path.join("..", "tests", "data", "MTBLS79_subset", "MTBLS79_mzml_triplicates.zip")
    fn_filelist = os.path.join("..", "tests", "data", "MTBLS79_subset", "filelist_mzml_triplicates.txt")
    output = os.path.join("results")
    if not os.path.exists(output):
        os.mkdir(output)

    print("Unzip mzml files.....")
    zip_ref = zipfile.ZipFile(source, 'r')
    zip_ref.extractall(os.path.join("data"))
    zip_ref.close()
    print("Completed")

    print("Process Scans.....")
    pls = process_scans("data", min_scans=1, function_noise="median",
                        snr_thres=3.0, ppm=2.0, min_fraction=None, rsd_thres=None,
                        filelist=fn_filelist, remove_mz_range=[], block_size=5000, ncpus=None)
    print("Completed")

    print("Replicate Filter.....")
    logfile = os.path.join(output, "log_replicate_filter.txt")
    pls_rf = replicate_filter(pls, ppm=2.0, replicates=3, min_peaks=2, rsd_thres=None, report=logfile, block_size=5000)
    print("Completed")

    print("Write each replicate filtered peaklist to a text file")
    for pl in pls_rf:
        with open(pl.ID + ".txt", "w") as out:
            out.write(os.path.join("results", pl.to_str("\t")))
    print("Completed")

    # print("Save, write and load peaklists")
    # save_peaklists_as_hdf5(pls_rf, os.path.join(output, "pls_rf.h5"))
    # hdf5_peaklists_to_txt(os.path.join(output, "pls_rf.h5"), path_out=output)
    # pls_rf = load_peaklists_from_hdf5(os.path.join(output, "pls_rf.h5"))
    # print("Completed")

    # print("Create a new sample list.....")
    # sample_list = os.path.join(output, "sample_list.txt")
    # create_sample_list(pls_rf, sample_list, delimiter="\t")
    # print("Completed")
    # print("")

    print("Align Samples.....")
    pm = align_samples(pls_rf, ppm=3.0, ncpus=1, block_size=5000)
    print("Completed", pm.shape)

    # print("Save, write and load peak matrix")
    # save_peak_matrix_as_hdf5(pm, os.path.join(output, "pm.h5"))
    # hdf5_peak_matrix_to_txt(os.path.join(output, "pm.h5"), path_out=os.path.join(output, "pm.txt"), attr_name="intensity", comprehensive=True)
    # pm = load_peak_matrix_from_hdf5(os.path.join(output, "pm.h5"))
    # print("Completed")

    print("Blank Filter.....")
    pm_bf = blank_filter(pm, "blank", min_fraction=1.0, min_fold_change=10.0, function="mean", rm_samples=True)
    print("Completed", pm_bf.shape)

    print("Sample Filter.....")
    pm_bf_sf = sample_filter(pm, 0.8, within=False)
    print("Completed", pm_bf_sf.shape)


if __name__ == '__main__':
    main()
