#!/usr/bin/env bash

dimspy --help

dimspy unzip \
--input ../tests/data/MTBLS79_subset/MTBLS79_mzml_triplicates.zip \
--output results/mzml

dimspy process-scans \
--input results/mzml \
--output results/peaklists.hdf5 \
--filelist ../tests/data/MTBLS79_subset/filelist_mzml_triplicates.txt \
--function-noise median \
--snr-threshold 3.0 \
--ppm 2.0 \
--min_scans 1 \
--min-fraction 0.5 \
--block-size 5000 \
--ncpus 2

dimspy replicate-filter \
--input results/peaklists.hdf5 \
--output results/peaklists_rf.hdf5 \
--ppm 2.0 \
--replicates 3 \
--min-peak-present 2

dimspy align-samples \
--input results/peaklists.hdf5 \
--output results/pm_a.hdf5 \
--ppm 2.0

dimspy blank-filter \
--input results/pm_a.hdf5 \
--output results/pm_a_bf.hdf5 \
--blank-label blank \
--remove

dimspy sample-filter \
--input results/pm_a_bf.hdf5 \
--output results/pm_a_bf_sf.hdf5 \
--min-fraction 0.8

dimspy hdf5-pls-to-txt \
--input results/peaklists_rf.hdf5 \
--output results \
--delimiter tab

dimspy hdf5-pm-to-txt \
--input results/pm_a_bf_sf.hdf5 \
--output results/pm_a_bf_sf.txt \
--delimiter tab

dimspy merge-peaklists \
--input results/peaklists_rf.hdf5 \
--input results/peaklists.hdf5 \
--output results/peaklists_merged.hdf5
