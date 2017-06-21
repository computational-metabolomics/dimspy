#!/usr/bin/python
# -*- coding: utf-8 -*-

import logging
import collections
import os
import zipfile
import numpy as np
from dimspy.models.peaklist import PeakList
from dimspy.portals import mzml_portal
from dimspy.portals import thermo_raw_portal
from dimspy.process.peak_alignment import align_peaks
from dimspy.process.peak_filters import filter_attr
from dimspy.experiment import define_mz_ranges
from dimspy.experiment import interpret_experiment_from_headers
from dimspy.experiment import remove_headers
from dimspy.experiment import mz_range_from_header
from string import join

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

    if type(pls_sd) is not dict and type(pls_sd) is not collections.OrderedDict:
        raise TypeError("Incorrect format - dict or collections.OrderedDict required")

    mzrs = [mz_range_from_header(h) for h in pls_sd]
    new_mzrs = _calculate_edges(mzrs)
    for h in pls_sd.keys():
        mz_ranges = len(pls_sd[h]) * [new_mzrs[pls_sd.keys().index(h)]]
        for i in range(len(pls_sd[h])):
            remove = [np.where(pls_sd[h][i].mz == mz)[0][0] for mz in pls_sd[h][i].mz if mz < mz_ranges[i][0] or mz >= mz_ranges[i][1]]
            pls_sd[h][i].remove_peak(remove)
    return pls_sd


def read_scans(fn, source, function_noise, nscans, scan_events={}):

    fn = fn.encode('string-escape')
    source = source.encode('string-escape')
    
    # assert os.path.isfile(fn), "File does not exist"
    if not fn.lower().endswith(".mzml") and not fn.lower().endswith(".raw"):
        raise IOError("Check format raw data (.RAW or .mzML)")

    if nscans is not None and type(nscans) is not int:
        raise ValueError("Integer (>= 0) or None required for nscans")

    if zipfile.is_zipfile(source):
        if fn.lower().endswith(".mzml"):
            run = mzml_portal.Mzml(fn, source)
        elif fn.lower().endswith(".raw"):
            raise IOError("Zip file with raw files not supported")
        else:
            raise IOError("Incorrect format: {}".format(os.path.basename(fn)))
    else:
        if fn.lower().endswith(".mzml"):
            run = mzml_portal.Mzml(fn)
        elif fn.lower().endswith(".raw"):
            run = thermo_raw_portal.ThermoRaw(fn)
        else:
            raise IOError("Incorrect format: {}".format(os.path.basename(fn)))

    h_sids = run.headers_scan_ids()
    mzrs = collections.OrderedDict(zip(h_sids.keys(), [mz_range_from_header(h) for h in h_sids]))

    if scan_events is None or scan_events == {}:
        h_rm = interpret_experiment_from_headers(mzrs)
        h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)
    elif type(scan_events) is dict:
        subset = define_mz_ranges(scan_events)
        h_rm = remove_headers(subset, mzrs)
        h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)
    elif scan_events == "all":
        pass

    #elif os.path.isfile(scan_events.encode('string-escape')):
    #    print "Reading scans....."
    #    with open(scan_events.encode('string-escape'), 'r') as f:
    #        mzrs_from_fn = [line.strip().split("\t") for line in f]
    #        subset = define_mz_ranges(mzrs_from_fn)
    #        h_rm = remove_headers(subset, mzrs)
    #        h_sids = collections.OrderedDict((key, value) for key, value in h_sids.items() if key in h_rm)

    # Validate that there are enough scans for each window
    if nscans is not None:
        if min([len(scans) for h, scans in h_sids.items()]) < nscans:
            raise IOError("not enough scans for each window, nscans = {}".format(nscans))
    #retireve scan data / create a peaklist class for each scan

    scans = collections.OrderedDict()
    for h, sids in h_sids.iteritems():
        if nscans is not None:
            if nscans > 0:
                sids = sids[0:nscans]
        scans[h] = run.peaklists(sids, function_noise)
    return scans


def average_replicate_scans(pls, snr_thres=3.0, ppm=2.0, min_fraction=0.8, rsd_thres=30.0, block_size=2000, ncpus=None):

    print "Removing noise....."
    pls_c = collections.OrderedDict()
    for h in pls:
        pls_c[h] = [filter_attr(pl.copy(), "snr", min_threshold=snr_thres) for pl in pls[h] if len(pl.mz) > 0]

    print "Aligning, averaging and filtering peaks....."
    for h in pls_c:
        print h
        emlst = np.array(map(lambda x: x.size == 0, pls_c[h]))
        if np.sum(emlst) > 0:
            logging.warning('droping empty peaklist(s) [%s]' % join(map(str, [p.ID for e, p in zip(emlst,  pls_c[h]) if e]), ','))
            pls_c[h] = [p for e, p in zip(emlst,  pls_c[h]) if not e]

        if len(pls_c[h]) >= 1:
            pm = align_peaks(pls_c[h], ppm=ppm, block_size=block_size, ncpus=ncpus)
            # TODO: remove clusters that have a higher number of peaks than samples
            # OR we can take the most accurate group of peaks and remove remaining peaks
            # Better to first remove clusters of higher number of peaks and log it
            pls_c[h] = pm.to_peaklist(ID=h)

            pls_c[h].add_attribute("snr", pm.attr_mean_vector('snr'), on_index=2)
            pls_c[h].add_attribute("snr_flag", np.ones(pls_c[h].full_size), flagged_only=False, is_flag=True)

            if min_fraction is not None:
                pls_c[h].add_attribute("fraction_flag", (pm.present / float(pm.shape[0])) >= min_fraction, flagged_only=False, is_flag=True)
            if rsd_thres is not None:
                if pm.shape[0] == 1:
                    logging.warning('applying RSD filter on single scan, all peaks removed')
                rsd_flag = map(lambda x: not np.isnan(x) and x < snr_thres, pm.rsd)
                pls_c[h].add_attribute("rsd_flag", rsd_flag, flagged_only=False, is_flag=True)
        else:
            logging.warning("No scan data available for {}".format(h))
            del pls_c[h]
    return pls_c


def join_peaklists(ID, pls):
    def _join_atrtributes(pls):
        attrs_out = collections.OrderedDict()
        for pl in pls:
            for atr in pl.attributes:
                attrs_out.setdefault(atr, []).extend(list(pl.get_attribute(atr, flagged_only=False)))
            if list(pl.attributes) != attrs_out.keys():
                raise IOError("Different attributes")
        return attrs_out

    attrs = _join_atrtributes(pls.values())
    pl_j = PeakList(ID=ID, mz=attrs["mz"], intensity=attrs["intensity"])
    del attrs["mz"], attrs["intensity"]  # default attributes
    for a in attrs:
        pl_j.add_attribute(a, attrs[a], is_flag=(a in pls.values()[0].flag_attributes), flagged_only=False)
    return pl_j
