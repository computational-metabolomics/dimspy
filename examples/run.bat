python -m dimspy check-filelist --source E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data --filelist E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data\filelist.tsv

python -m dimspy process-scans^
 --source E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data^
 --pickle-file-out E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data\peaklist.pkl^
 --filelist E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data\filelist.tsv^
 --filename-experiment None^
 --mode-noise median^
 --snr-threshold 3.0^
 --ppm 2.0^
 --nscans 1^
 --min-peaks 1^
 --rsd-threshold 30.0^
 --block-size 2000^
 --ncpus 2

python -m dimspy replicate-filter^
 --pickle-file-in E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data\peaklist.pkl^
 --pickle-file-out E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data\peaklist_rf.pkl^
 --ppm 2.0^
 --min-replicates 2

python -m dimspy sample-filter^
 --pickle-file-in E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data\pm.pkl^
 --pickle-file-out E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\tests\data\pm_sf.pkl^
 --min-fraction 0.25^
 --rsd-threshold 30.00^
 --within