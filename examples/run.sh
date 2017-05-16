#!/usr/bin/env bash

python -m dimspy process-scans \
--input tests/data/MTBLS79_subset/MTBLS79_subset.zip \
--output tests/data/temp/peaklists.hdf5 \
--filelist tests/data/MTBLS79_subset/filelist_mzML.txt \
--subset-scan-events None \
--function-noise median \
--snr-threshold 3.0 \
--ppm 2.0 \
--nscans 1 \
--min-fraction 1 \
--rsd-threshold 30.0 \
--block-size 2000 \
--ncpus 2

python -m dimspy replicate-filter \
--input tests/data/temp/peaklists.hdf5 \
--output tests/data/temp/peaklists_rf.hdf5 \
--ppm 2.0 \
--replicates 3 \
--min-peak-present 2 \
--rsd-threshold 30.0

python -m dimspy align-samples \
--input tests/data/temp/peaklists_rf.hdf5 \
--output tests/data/temp/pm_a.hdf5 \
--ppm 2.0

python -m dimspy blank-filter \
--input tests/data/temp/pm_a.hdf5 \
--output tests/data/temp/pm_a_bf.hdf5 \
--blank-label blank \
--remove

python -m dimspy sample-filter \
--input tests/data/temp/pm_a_bf.hdf5 \
--output tests/data/temp/pm_a_bf_sf.hdf5 \
--min-fraction 0.8 \
--rsd-threshold 30.00

python -m dimspy hdf5-to-txt \
--input tests/data/temp/peaklists_rf.hdf5 \
--output tests/data/temp

python -m dimspy hdf5-to-txt \
--input tests/data/temp/pm_a_bf_sf.hdf5 \
--output tests/data/temp/pm_a_bf_sf.txt
