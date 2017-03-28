import collections
import copy
import os
import re
import zipfile

import numpy as np

from dimspy.models.peaklist import PeakList
from dimspy.parsers import Mzml
from dimspy.parsers import ThermoRaw, ThermoRawWine
from dimspy.process.peak_alignment import align_peaks
from dimspy.process.peak_filters import filter_snr
from dimspy.experiment import read_ms_exp_from_file
from dimspy.experiment import read_ms_exp_from_headers
from dimspy.experiment import remove_headers


def _calculate_edges(mz_ranges):
    s_mz_ranges = map(sorted, mz_ranges)
    if len(s_mz_ranges) == 1: return s_mz_ranges
    
    s_min,s_max = zip(*s_mz_ranges)
    assert all(map(lambda x: x[0] < x[1], zip(s_min[:-1], s_min[1:]))), 'start values not in order'
    assert all(map(lambda x: x[0] < x[1], zip(s_max[:-1], s_max[1:]))), 'end values not in order'
    
    s_zip = zip(s_min[1:], s_max[:-1])
    e_size = map(lambda x: (x[1]-x[0]) * 0.5, s_zip)
    assert all(map(lambda x: x > 0, e_size)), 'incorrect overlap'
    
    merged = (s_min[0],) + reduce(lambda x,y: x+y, [(z[0]+e, z[1]-e) for z,e in zip(s_zip, e_size)]) + (s_max[-1],)
    return zip(merged[::2], merged[1::2])


# TODO: following is temporarily retained in case the new codes fail to cover all the situations
# mz_ranges.sort(key=lambda x: x[0])
#
# import pdb; pdb.set_trace()
#
# if len(mz_ranges) == 1:
#     return mz_ranges
# else:
#     st_mzrs = map(lambda x: np.min(x), mz_ranges)
#     assert map(lambda x: x[0] < x[1], zip(st_mzrs[:-1], st_mzrs[1:])), 'start values not in order'
#
#     ed_mzrs = map(lambda x: np.max(x), mz_ranges)
#     assert map(lambda x: x[0] < x[1], zip(ed_mzrs[:-1], ed_mzrs[1:])), 'end values not in order'
#
#     new_mzrs = copy.deepcopy(mz_ranges)
#
#     for i in range(0, len(new_mzrs) - 1):
#
#         overlap = max(0, min(new_mzrs[i][1], new_mzrs[i+1][1]) - max(new_mzrs[i][0], new_mzrs[i+1][0]))
#         assert overlap > 0, "incorrect overlap"
#
#         if i == 0:
#             # first window
#             new_mzrs[i] = (new_mzrs[i][0] + (overlap * 0.5), new_mzrs[i][1] - (overlap * 0.5))
#             new_mzrs[i + 1] = (new_mzrs[i + 1][0] + (overlap * 0.5), new_mzrs[i + 1][1])
#         elif i == len(new_mzrs):
#             # final window
#             new_mzrs[i] = (new_mzrs[i][0], new_mzrs[i][1] - (overlap * 0.5))
#             new_mzrs[i + 1] = (new_mzrs[i + 1][0] + (overlap * 0.5), new_mzrs[i + 1][1] - (overlap * 0.5))
#         else:
#             # windows in between
#             new_mzrs[i] = (new_mzrs[i][0], new_mzrs[i][1] - (overlap * 0.5))
#             new_mzrs[i + 1] = (new_mzrs[i + 1][0] + (overlap * 0.5), new_mzrs[i + 1][1])
#
# return new_mzrs


def mz_range_from_header(h):
    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


def remove_edges(pkls_sd):

    assert type(pkls_sd) == dict or type(pkls_sd) == collections.OrderedDict, "Incorrect format ()"

    mzrs = [mz_range_from_header(h) for h in pkls_sd]
    new_mzrs = _calculate_edges(mzrs)
    for h in pkls_sd.keys():
        mz_ranges = len(pkls_sd[h]) * [new_mzrs[pkls_sd.keys().index(h)]]
        for i in range(len(pkls_sd[h])):
            remove = [np.where(pkls_sd[h][i].mz == mz)[0][0] for mz in pkls_sd[h][i].mz if mz < mz_ranges[i][0] or mz >= mz_ranges[i][1]]
            pkls_sd[h][i].remove_peak(remove)
    return pkls_sd


def read_scan_data(fn, source, function_noise, nscans, fn_exp):

    fn = fn.encode('string-escape')
    source = source.encode('string-escape')
    if fn_exp:
        fn_exp = fn_exp.encode('string-escape')

    # assert os.path.isfile(fn), "File does not exist"
    assert fn.lower().endswith(".mzml") or fn.lower().endswith(".raw"), "Check format raw data (.RAW or .mzML)"

    assert type(nscans) == int and nscans >= 0, "User a integer >= 0"

    # TEMPORARY, BETTER TO CHECK PART OF THE CONTENT
    if fn.lower().endswith(".mzml"):
        if zipfile.is_zipfile(source):
            run = Mzml(fn, source)
        else:
            run = Mzml(fn)

    elif fn.lower().endswith(".raw"):
        if os.name == "nt":
            run = ThermoRaw(fn, source)
        else:
            run = ThermoRawWine(fn, source)

    h_sids = run.headers_scan_ids()
    mzrs = collections.OrderedDict(zip(h_sids.keys(), [mz_range_from_header(h) for h in h_sids]))

    if os.path.isfile(str(fn_exp)):
        exp_desc = read_ms_exp_from_file(fn_exp)
        h_rm = remove_headers(exp_desc, mzrs)
        h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)
    else:
        h_rm = read_ms_exp_from_headers(mzrs)
        h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)

    # Validate that there are enough scans for each window
    assert min([len(scans) for h, scans in h_sids.items()]) >= nscans, "not enough scans for each window, nscans = {}".format(nscans)
    #retireve scan data / create a peaklist class for each scan

    scans = collections.OrderedDict()
    for h, sids in h_sids.iteritems():
        if nscans > 0:
            sids = sids[0:nscans]
        scans[h] = run.peaklists(sids, function_noise)

    return scans


def process_replicate_scans(pkls, snr_thres=3.0, ppm=2.0, min_fraction=0.8, rsd_thres=30.0, block_size=2000, ncpus=None):
    print "Removing noise....."
    for h in pkls:
        pkls[h] = [filter_snr(pk, snr_thres) for pk in pkls[h] if len(pk.mz) > 0] # TODO: empty scans

    print "Align and filtering peaks from scans....."
    for h in pkls:
        print h
        if len(pkls[h]) > 1:
            pm = align_peaks(pkls[h], ppm=ppm, block_size=block_size, ncpus=ncpus)
            # remove clusters that have a higher number of peaks than samples
            # OR we can take the most accurate group of peaks and remove remaining peaks
            # Better to first remove clusters of higher number of peaks and log it
            pkls[h] = pm.to_peaklist(ID=h)
            pkls[h].add_attribute("snr", pm.attr_mean_vector('snr'))
            pkls[h].add_attribute("snr_flag", np.ones(pkls[h].full_size), flagged_only=False, is_flag=True)
            pkls[h].add_attribute("present", pm.present) # np.zeros(len(pkls[h].mz)))
            pkls[h].add_attribute("fraction", pm.present / float(pm.shape[0]))
            pkls[h].add_attribute("rsd", pm.rsd)
            pkls[h].add_attribute("fraction_flag", (pm.present / float(pm.shape[0])) >= min_fraction, flagged_only=False, is_flag=True)
            if rsd_thres is None:
                pkls[h].add_attribute("rsd_flag", np.ones(pkls[h].full_size), flagged_only=False, is_flag=True)
            else:
                #pkls[h].add_attribute("rsd_flag", np.logical_or(np.isnan(pm.rsd), pm.rsd < rsd_thres), flagged_only=False, is_flag=True)
                pkls[h].add_attribute("rsd_flag", pm.rsd <= rsd_thres, flagged_only=False, is_flag=True)

        elif len(pkls[h]) == 1:
            pkls[h] = pkls[h][0]
            pkls[h].add_attribute("present", np.ones(pkls[h].full_size), flagged_only = False)
            pkls[h].add_attribute("fraction", np.ones(pkls[h].full_size), flagged_only = False)
            pkls[h].add_attribute("fraction_flag", np.ones(pkls[h].full_size), flagged_only=False, is_flag=True)
        else:
            print "No scans available for {}".format(h) # TODO: check if it is valid to remove header
            del pkls[h]
    return pkls


def join_peaklists(ID, pkls, class_label):

    def _join_atrtributes(pls):
        attrs_out = collections.OrderedDict()
        for pl in pls:
            for atr in pl.attributes:
                attrs_out.setdefault(atr, []).extend(list(pl.get_attribute(atr, flagged_only=False)))
            assert list(pl.attributes) == attrs_out.keys(), "different attributes"
        return attrs_out

    attrs = _join_atrtributes(pkls.values())
    pl_j = PeakList(ID=ID, mz=attrs["mz"], intensity=attrs["intensity"], tags=class_label)
    pl_j.metadata.ion_injection_time = 1.0
    del attrs["mz"], attrs["intensity"]  # default attributes
    for a in attrs:
        #pl_j.add_attribute(a, attrs[a])
        pl_j.add_attribute(a, attrs[a], is_flag=(a in pkls.values()[0].flag_attributes), flagged_only=False)
    return pl_j
