python -m dimspy process-scans^
 --input tests/data/MTBLS79_subset/MTBLS79_mzml_triplicates.zip^
 --output tests/test_results/peaklists.hdf5^
 --filelist tests/data/MTBLS79_subset/filelist_mzml_triplicates.txt^
 --function-noise median^
 --snr-threshold 3.0^
 --ppm 2.0^
 --min_scans 1^
 --min-fraction 0.5^
 --block-size 5000^
 --ncpus 2

python -m dimspy replicate-filter^
 --input tests/test_results/peaklists.hdf5^
 --output tests/test_results/peaklists_rf.hdf5^
 --ppm 2.0^
 --replicates 3^
 --min-peak-present 2

python -m dimspy align-samples^
 --input tests/test_results/peaklists_rf.hdf5^
 --output tests/test_results/pm_a.hdf5^
 --ppm 2.0

python -m dimspy blank-filter^
 --input tests/test_results/pm_a.hdf5^
 --output tests/test_results/pm_a_bf.hdf5^
 --blank-label blank^
 --remove

python -m dimspy sample-filter^
 --input tests/test_results/pm_a_bf.hdf5^
 --output tests/test_results/pm_a_bf_sf.hdf5^
 --min-fraction 0.8

python -m dimspy hdf5-pls-to-txt^
 --input tests/test_results/peaklists_rf.hdf5^
 --output tests/test_results^
 --delimiter tab

python -m dimspy hdf5-pm-to-txt^
 --input tests/test_results/pm_a_bf_sf.hdf5^
 --output tests/test_results/pm_a_bf_sf.txt^
 --delimiter tab
