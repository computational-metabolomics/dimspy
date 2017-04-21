#!/usr/bin/python
# -*- coding: utf-8 -*-

import logging
import collections
import os
import zipfile
import numpy as np
from dimspy.models.peaklist import PeakList
from dimspy.portals import Mzml
from dimspy.portals import ThermoRaw
from dimspy.process.peak_alignment import align_peaks
from dimspy.process.peak_filters import filter_snr
from dimspy.experiment import define_mz_ranges
from dimspy.experiment import interpret_experiment_from_headers
from dimspy.experiment import remove_headers
from dimspy.experiment import mz_range_from_header


def _calculate_edges(mz_ranges):
    s_mz_ranges = map(sorted, mz_ranges)
    if len(s_mz_ranges) == 1: return s_mz_ranges
    
    s_min, s_max = zip(*s_mz_ranges)
    assert all(map(lambda x: x[0] < x[1], zip(s_min[:-1], s_min[1:]))), 'start values not in order'
    assert all(map(lambda x: x[0] < x[1], zip(s_max[:-1], s_max[1:]))), 'end values not in order'
    
    s_zip = zip(s_min[1:], s_max[:-1])
    e_size = map(lambda x: (x[1]-x[0]) * 0.5, s_zip)
    assert all(map(lambda x: x > 0, e_size)), 'incorrect overlap'
    
    merged = (s_min[0],) + reduce(lambda x, y: x+y, [(z[0]+e, z[1]-e) for z, e in zip(s_zip, e_size)]) + (s_max[-1],)
    return zip(merged[::2], merged[1::2])


def remove_edges(pls_sd):

    assert type(pls_sd) == dict or type(pls_sd) == collections.OrderedDict, "Incorrect format ()"

    mzrs = [mz_range_from_header(h) for h in pls_sd]
    new_mzrs = _calculate_edges(mzrs)
    for h in pls_sd.keys():
        mz_ranges = len(pls_sd[h]) * [new_mzrs[pls_sd.keys().index(h)]]
        for i in range(len(pls_sd[h])):
            remove = [np.where(pls_sd[h][i].mz == mz)[0][0] for mz in pls_sd[h][i].mz if mz < mz_ranges[i][0] or mz >= mz_ranges[i][1]]
            pls_sd[h][i].remove_peak(remove)
    return pls_sd


def read_scans(fn, source, function_noise, nscans, subset_scan_events=None):

    fn = fn.encode('string-escape')
    source = source.encode('string-escape')

    # assert os.path.isfile(fn), "File does not exist"
    assert fn.lower().endswith(".mzml") or fn.lower().endswith(".raw"), "Check format raw data (.RAW or .mzML)"
    assert type(nscans) == int and nscans >= 0, "User a integer >= 0"

    if fn.lower().endswith(".mzml"):
        if zipfile.is_zipfile(source):
            run = Mzml(fn, source)
        else:
            run = Mzml(fn)
    elif fn.lower().endswith(".raw"):
        run = ThermoRaw(fn)
    else:
        pass

    h_sids = run.headers_scan_ids()
    mzrs = collections.OrderedDict(zip(h_sids.keys(), [mz_range_from_header(h) for h in h_sids]))

    if subset_scan_events is None:
        h_rm = interpret_experiment_from_headers(mzrs)
        h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)
    elif type(subset_scan_events) == list:
        subset = define_mz_ranges(subset_scan_events)
        h_rm = remove_headers(subset, mzrs)
        h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)
    elif os.path.isfile(subset_scan_events.encode('string-escape')):
        with open(subset_scan_events.encode('string-escape'), 'r') as f:
            mzrs_from_fn = [line.strip().split("\t") for line in f]
            subset = define_mz_ranges(mzrs_from_fn)
            h_rm = remove_headers(subset, mzrs)
            h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)
    elif subset_scan_events == "all":
        pass

    # Validate that there are enough scans for each window
    assert min([len(scans) for h, scans in h_sids.items()]) >= nscans, "not enough scans for each window, nscans = {}".format(nscans)
    #retireve scan data / create a peaklist class for each scan

    scans = collections.OrderedDict()
    for h, sids in h_sids.iteritems():
        if nscans > 0:
            sids = sids[0:nscans]
        scans[h] = run.peaklists(sids, function_noise)
    return scans


def average_replicate_scans(pls, snr_thres=3.0, ppm=2.0, min_fraction=0.8, rsd_thres=30.0, block_size=2000, ncpus=None):

    print "Removing noise....."
    for h in pls:
        pls[h] = [filter_snr(pl, snr_thres) for pl in pls[h] if len(pl.mz) > 0]

    print "Align, averaging and filtering peaks....."
    for h in pls:
        print h
        if len(pls[h]) > 1:
            pm = align_peaks(pls[h], ppm=ppm, block_size=block_size, ncpus=ncpus)
            # TODO: remove clusters that have a higher number of peaks than samples
            # OR we can take the most accurate group of peaks and remove remaining peaks
            # Better to first remove clusters of higher number of peaks and log it
            pls[h] = pm.to_peaklist(ID=h)

            pls[h].add_attribute("snr", pm.attr_mean_vector('snr'))
            pls[h].add_attribute("snr_flag", np.ones(pls[h].full_size), flagged_only=False, is_flag=True)
            pls[h].add_attribute("present", pm.present) # np.zeros(len(pls[h].mz)))
            pls[h].add_attribute("fraction", pm.present / float(pm.shape[0]))
            pls[h].add_attribute("rsd", pm.rsd)

            if min_fraction is not None:
                pls[h].add_attribute("fraction_flag", (pm.present / float(pm.shape[0])) >= min_fraction, flagged_only=False, is_flag=True)
            if rsd_thres is not None:
                pls[h].add_attribute("rsd_flag", pm.rsd <= rsd_thres, flagged_only=False, is_flag=True)

        elif len(pls[h]) == 1:
            pls[h] = pls[h][0]
            # snr and snr_flag attribute already available
            pls[h].remove_unflagged_peaks('snr_flag')
            pls[h].add_attribute("present", np.ones(pls[h].full_size), flagged_only=False)
            pls[h].add_attribute("fraction", np.ones(pls[h].full_size), flagged_only=False)
            pls[h].add_attribute("rsd", np.nan * np.ones(pls[h].full_size), flagged_only=False)
            if min_fraction is not None:
                pls[h].add_attribute("fraction_flag", np.ones(pls[h].full_size), flagged_only=False, is_flag=True)
            if rsd_thres is not None:
                logging.warning('applying RSD filter on single scan, all peaks removed')
                pls[h].add_attribute("rsd_flag", np.zeros(pls[h].full_size), flagged_only=False, is_flag=True)
        else:
            print "No scans available for {}".format(h)
            del pls[h]
    return pls


def join_peaklists(ID, pls):

    def _join_atrtributes(pls):
        attrs_out = collections.OrderedDict()
        for pl in pls:
            for atr in pl.attributes:
                attrs_out.setdefault(atr, []).extend(list(pl.get_attribute(atr, flagged_only=False)))
            assert list(pl.attributes) == attrs_out.keys(), "different attributes"
        return attrs_out

    attrs = _join_atrtributes(pls.values())
    pl_j = PeakList(ID=ID, mz=attrs["mz"], intensity=attrs["intensity"])
    del attrs["mz"], attrs["intensity"]  # default attributes
    for a in attrs:
        pl_j.add_attribute(a, attrs[a], is_flag=(a in pls.values()[0].flag_attributes), flagged_only=False)
    return pl_j
