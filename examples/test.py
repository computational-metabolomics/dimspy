import cPickle as cp
from dimspy.workflow import process_scans
from dimspy.workflow import replicate_filter
from dimspy.workflow import align_samples
from dimspy.workflow import blank_filter
from dimspy.workflow import sample_filter
from dimspy.models.peaklist import PeakList
from dimspy.models.peak_matrix import PeakMatrix


def main():
    import os

    #my_path = '/Users/Albert/Repository/Shoukaku/data/'
    my_path = os.path.dirname(os.path.realpath(__file__))


    # load raw files
    # ----------------------------------------------
    # fn_fls = os.path.join(my_path, 'tests', 'MS filename_polar_DIMS_ZF_May2016.txt')
    # source = os.path.join(my_path, 'tests', 'zebrafish_raw')

    # fn_fls = os.path.join(my_path, 'tests', 'data', 'filelist.tsv')
    # source = os.path.join(my_path, 'tests', 'data', "MTBLS79_subset")
    #
    # # ZIP FILE WORKS!
    fn_fls = os.path.join(my_path, 'tests', 'data', 'filelist_mzML.tsv')
    source = os.path.join(my_path, 'tests', 'data', "MTBLS79_subset", "MTBLS79_subset.zip")

    # porder = [pdct[i] for i in pids]
    # ALL SCANS AND NO RSD FILTER
    #pls = process_scans(source, fn_fls, nscans=None, fn_exp=None, function_noise="median",
    #                    snr_thres=3.0, ppm=2.0, min_fraction=0.5, rsd_thres=None, block_size=200, ncpus=None)

    # File: Batch04_S01_rep02_248.RAW
    # ValueError: negative dimensions are not allowed
    # THREE SCANS AND NO RSD FILTER
    # pls = process_scans(source, fn_fls, nscans=3, fn_exp=None, function_noise="median",
    #                    snr_thres=3.0, ppm=2.0, min_fraction=0.5, rsd_thres=None, block_size=200, ncpus=None)

    # SINGLE SCANS
    pls = process_scans(source, fn_fls, nscans=5, fn_exp=None, function_noise="median",
                        snr_thres=3.0, ppm=2.0, min_fraction=0.5, rsd_thres=30.0, block_size=2000, ncpus=None)
    # export
    for pl in pls:
        print pl.ID, pl.shape
        #cp.dump(pl, open(os.path.join(my_path, 'tests', 'data', pl.ID + ".pkl"), "w"))
        with open(os.path.join(my_path, 'tests', 'data', pl.ID + ".txt"), "w") as out: out.write(pl.to_str("\t"))

    # # import
    # src_path = os.path.join(my_path, 'tests', 'data')
    # files = filter(os.path.isfile, [os.path.join(src_path, fn) for fn in os.listdir(src_path) if fn.endswith('.pkl') and not fn.startswith('.')])
    # pls = map(lambda x: cp.load(open(x, 'r')), files)
    # print '[%d] peaklists loaded' % len(pls)


    # replicate Filter
    # ----------------------------------------------
    print
    print "Replicate Filter...."
    pls_rf = replicate_filter(pls, ppm=2.0, reps=3, minpeaks=2, rsd_thres=20.0)
    print "Finished"

    # export
    for pl in pls_rf:
        print pl.ID, pl.shape
        #cp.dump(pl, open(os.path.join(my_path, 'tests', 'data', pl.ID + ".pkl"), "w"))
        with open(os.path.join(my_path, 'tests', 'data', pl.ID + ".txt"), "w") as out: out.write(pl.to_str("\t"))

    # # import
    # src_path = os.path.join(my_path, 'tests', 'data')
    # files = filter(os.path.isfile, [os.path.join(src_path, fn) for fn in os.listdir(src_path) if fn.endswith('.pkl') and not fn.startswith('.')])
    # pls_rf = map(lambda x: cp.load(open(x, 'r')), files)
    # print '[%d] peaklists loaded' % len(pls_rf)


    # alignment
    # ----------------------------------------------
    print
    print "Align Samples...."
    pm = align_samples(pls_rf, ppm=2.0)
    print "Finished", pm.shape

    cp.dump(pm, open(os.path.join(my_path, 'tests', 'data', "pm.pkl"), "w"))
    with open(os.path.join(my_path, 'tests', 'data', "pm.txt"), "w") as out: out.write(pm.to_str('\t'))

    # pm = cp.load(open(os.path.join(my_path, 'tests', 'data', "peak_matrix.pkl"), "r"))
    # print '[%s] peak matrix loaded' % str(pm.shape)[1:-1]


    # blank Filter
    # ----------------------------------------------
    print
    print "Blank Filter"

    # E:\Dropbox\Projects\Bioinformatics\Development\DIMS-workflow-development\computational\dimspy\dimspy\models\peak_matrix.py:57: VisibleDeprecationWarning: boolean index did n
    # ot match indexed array along dimension 0; dimension is 7 but corresponding boolean dimension is 6
    # return self._tags[self._mask]
    pm_bf = blank_filter(pm, "blank", min_fraction=1.0, min_fold_change=10.0, function="mean", rm_samples=True)
    print "Finished", pm_bf.shape
    print pm_bf.peaklist_ids
    print pm_bf.peaklist_tag_values

    cp.dump(pm, open(os.path.join(my_path, 'tests', 'data', "pm_bf.pkl"), "w"))
    with open(os.path.join(my_path, 'tests', 'data', "pm_bf.txt"), "w") as out: out.write(pm.to_str('\t'))

    # pm_bf = cp.load(open(os.path.join(my_path, 'tests', 'data', "pm_bf.pkl"), "r"))
    # print '[%s] peak matrix loaded' % str(pm_bf.shape)[1:-1]


    # sample Filter
    # ----------------------------------------------
    print
    print "Sample Filter"
    import pdb; pdb.set_trace()
    pm_bf_sf = sample_filter(pm, 0.8, within=False)
    print "Finished(1)", pm_bf_sf.shape
    pm_bf_sf = sample_filter(pm, 0.8, within=False, rsd=30.0, qc_label=None)
    print "Finished(2)", pm_bf_sf.shape
    
    # ERROR: rsd_values = pm.rsd
    # "Weights sum to zero, can't be normalized")
    pm_bf_sf = sample_filter(pm, 0.8, within=False, rsd=30.0, qc_label="QC")
    print "Finished(3)", pm_bf_sf.shape
    
    pm_bf_sf = sample_filter(pm, 0.8, within=True, rsd=30.0, qc_label=None)
    print "Finished(4)", pm_bf_sf.shape

    cp.dump(pm, open(os.path.join(my_path, 'tests', 'data', "pm_sf.pkl"), "w"))
    with open(os.path.join(my_path, 'tests', 'data', "pm_sf.txt"), "w") as out: out.write(pm.to_str('\t'))

    # pm_sf = cp.load(open(os.path.join(my_path, 'tests', 'data', "pm_sf.pkl"), "r"))
    # print '[%s] peak matrix loaded' % str(pm_sf.shape)[1:-1]


# main
if __name__ == '__main__':
    main()

