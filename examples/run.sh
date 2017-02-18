#!/usr/bin/env bash
python -m dimspy check-filelist --source /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data --filelist /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/filelist.tsv

python -m dimspy process-scans \
--source /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data \
--pickle-file-out /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/peaklists.pkl \
--filelist /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/filelist.tsv \
--filename-experiment None \
--mode-noise median \
--snr-threshold 3.0 \
--ppm 2.0 \
--nscans 1 \
--min-peaks 1 \
--rsd-threshold 30.0 \
--block-size 2000 \
--ncpus 2

python -m dimspy replicate-filter \
--pickle-file-in /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/peaklists.pkl \
--pickle-file-out /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/peaklist_rf.pkl \
--ppm 2.0 \
--replicates 3 \
--min-peaks 2

python -m dimspy align-samples \
--pickle-file-in /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/peaklist_rf.pkl \
--pickle-file-out /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/pm.pkl \
--ppm 2.0

python -m dimspy blank-filter \
--pickle-file-in /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/pm.pkl \
--pickle-file-out /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/pm_bf.pkl \
--blank-label blank \
--remove

python -m dimspy sample-filter \
--pickle-file-in /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/pm_bf.pkl \
--pickle-file-out /Users/RW/Dropbox/Projects/Bioinformatics/Development/DIMS-workflow-development/computational/dimspy/tests/data/pm_sf.pkl \
--min-fraction 0.8
--rsd-threshold 30.00
