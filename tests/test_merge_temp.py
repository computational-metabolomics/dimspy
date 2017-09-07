import unittest
from dimspy.workflow import *
from dimspy.portals.hdf5_portal import *
from dimspy.models import *



pls_01 = load_peaklists_from_hdf5("/home/tnl495/galaxy_testing/dimspy2/tests/data/MTBLS79_subset/batch04_QC17_rep01_262.hdf5")
pls_02 = load_peaklists_from_hdf5("/home/tnl495/galaxy_testing/dimspy2/tests/data/MTBLS79_subset/batch04_QC17_rep02_263.hdf5")
pls_03 = load_peaklists_from_hdf5("/home/tnl495/galaxy_testing/dimspy2/tests/data/MTBLS79_subset/batch04_QC17_rep03_264.hdf5")

print [pls_01[0], pls_02[0], pls_03[0]]


merged_peaklists = merge_peaklists([pls_01, pls_02, pls_03], "/home/tnl495/galaxy_testing/dimspy2/tests/data/MTBLS79_subset/filelist_multi.txt")
print merged_peaklists
for i in range(len(merged_peaklists)):
    hdf5_portal.save_peaklists_as_hdf5(merged_peaklists[i], 'merged_{}.hdf5'.format(i))

