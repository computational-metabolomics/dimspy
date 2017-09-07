from dimspy.workflow import *
from dimspy.portals.hdf5_portal import *
from dimspy.models import *




peaks1 = process_scans(source='/mnt/hgfs/DATA/nearline_dims_check/test', function_noise="median",
                            snr_thres=3.0, ppm=2.0, min_fraction=0.7, rsd_thres=None, filelist="/mnt/hgfs/DATA/nearline_dims_check/filelist_01.txt",
                            block_size=2000, ncpus=None, skip_stitching=False,
                       filter_scan_events={'include': [[50.0, 1000.0, 'full']]})
pm = align_samples(peaks1, ppm=2)
pm_bf = blank_filter(pm, blank_label='blank')
print pm_bf
#
#
# peaks2 = process_scans(source='/mnt/hgfs/DATA/nearline_dims_check/blank', function_noise="median",
#                             snr_thres=3.0,  ppm=2.0, min_fraction=0.7, rsd_thres=None, filelist="/mnt/hgfs/DATA/nearline_dims_check/filelist_02.txt",
#                             block_size=2000, ncpus=None, skip_stitching=False,
#                        filter_scan_events={'include': [[50.0, 1000.0, 'full']]})

# temp = []
# for p in peaks1:
#     temp.append(p)
# for p in peaks2:
#     temp.append(p)
# pm = align_samples(temp,  ppm=2)
# pm_bf = blank_filter(pm, blank_label='blank')

# merged_peaklists = merge_peaklists([peaks1, peaks2], "/mnt/hgfs/DATA/nearline_dims_check/filelist_multi_test.txt")
# print merged_peaklists

# pm = align_samples(merged_peaklists, ppm=2)
# pm_bf = blank_filter(pm, blank_label='blank')

# for i in range(len(merged_peaklists)):
#     print '####################'
#     print merged_peaklists[i][0].ID
#     print merged_peaklists[i][1].ID
#     pm = align_samples(merged_peaklists[i],  ppm=2)
#     pm_bf = blank_filter(pm, blank_label='blank')
#     print pm_bf


# for i in range(len(merged_peaklists)):
#     hdf5_portal.save_peaklists_as_hdf5(merged_peaklists[i], 'merged_{}.hdf5'.format(i))