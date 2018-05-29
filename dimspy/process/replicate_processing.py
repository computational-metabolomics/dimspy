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
from dimspy.experiment import scan_type_from_header
from dimspy.experiment import mz_range_from_header
from string import join


def _calculate_edges(mz_ranges):
    s_mz_ranges = map(sorted, mz_ranges)
    if len(s_mz_ranges) == 1:
        return s_mz_ranges
    
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
            for mz in pls_sd[h][i].mz:
                if mz < mz_ranges[i][0] or mz >= mz_ranges[i][1]:
                    remove.extend(list(np.where(pls_sd[h][i].mz == mz)[0]))
            pls_sd[h][i].remove_peak(remove)
    return pls_sd


def read_scans(fn, source, function_noise, min_scans=1, filter_scan_events=None):

    if filter_scan_events is None:
        filter_scan_events = {}
    if not fn.lower().endswith(".mzml") and not fn.lower().endswith(".raw"):
        raise IOError("Check format raw data (.RAW or .mzML)")

    if min_scans is not None and type(min_scans) is not int:
        raise ValueError("Integer (>= 1) or None required for min_scans")

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

    h_sids = run.headers()

    if type(filter_scan_events) is dict and len(filter_scan_events) > 0:

        if ("include" in filter_scan_events and "exclude" in filter_scan_events) or \
                ("include" not in filter_scan_events and "exclude" not in filter_scan_events):
            raise ValueError("Use 'exclude' or 'include' for filter_scan_events not both. E.g {'include': [[70.0, 170.0, 'sim']]}")

        if len([True for fse in filter_scan_events.values()[0] if len(fse) == 3]) != len(filter_scan_events.values()[0]):
            raise ValueError("Provide a start, end and scan type (sim or full) for filter_scan_events.")

        filter_scan_events = {filter_scan_events.keys()[0]:
                              [[float(fse[0]), float(fse[1]), str(fse[2])] for fse in filter_scan_events.values()[0]]}

        h_descs = {}
        for h in h_sids.copy():
            mzr = mz_range_from_header(h)
            h_descs[h] = [mzr[0], mzr[1], scan_type_from_header(h).lower()]

        incl_excl = filter_scan_events.keys()[0]
        for hd in filter_scan_events[incl_excl]:
            if hd not in h_descs.values():
                logging.warning("Event {} doest not exist".format(str(hd)))

        for hd in h_descs:
            if filter_scan_events.keys()[0] == "include":
                if h_descs[hd] not in filter_scan_events["include"]:
                    del h_sids[hd]
            elif filter_scan_events.keys()[0] == "exclude":
                if h_descs[hd] in filter_scan_events["exclude"]:
                    del h_sids[hd]

    if len(h_sids) == 0:
        raise Exception("No scan data to process. Check filter_scan_events")

    scans = collections.OrderedDict()
    for h, sids in h_sids.iteritems():
        if len(sids) >= min_scans:
            scans[h] = run.peaklists(sids, function_noise)
        else:
            logging.warning('Not enough scans for [{}] [{} < {}]. Scan event {} has been removed.'.format(h, len(scans), min_scans, h))
    return scans


def average_replicate_scans(name, pls, ppm=2.0, min_fraction=0.8, rsd_thres=30.0, rsd_on="intensity", block_size=5000, ncpus=None):

    emlst = np.array(map(lambda x: x.size == 0, pls))
    if np.sum(emlst) > 0:
        logging.warning('No scan data available for [%s]' % join(map(str, [p.ID for e, p in zip(emlst,  pls) if e]), ','))
        pls = [p for e, p in zip(emlst,  pls) if not e]

    pm = align_peaks(pls, ppm=ppm, block_size=block_size, ncpus=ncpus)

    pl_avg = pm.to_peaklist(ID=name)
    # meta data
    for pl in pls:
        for k, v in pl.metadata.items():
            if k not in pl_avg.metadata:
                pl_avg.metadata[k] = []
            if v is not None:
                pl_avg.metadata[k].append(v)

    if rsd_on != "intensity":
        pl_avg.add_attribute(rsd_on, pm.attr_mean_vector(rsd_on), on_index=2)
        rsd_label = "rsd_{}".format(rsd_on)
        shift = 1
    else:
        rsd_label = "rsd"
        shift = 0

    pl_avg.add_attribute("snr", pm.attr_mean_vector('snr'), on_index=2+shift)
    pl_avg.add_attribute("snr_flag", np.ones(pl_avg.full_size), flagged_only=False, is_flag=True)

    pl_avg.add_attribute(rsd_label, pm.rsd(on_attr=rsd_on, flagged_only=False), on_index=5+shift)

    if min_fraction is not None:
        pl_avg.add_attribute("fraction_flag", (pm.present / float(pm.shape[0])) >= min_fraction, flagged_only=False, is_flag=True)
    if rsd_thres is not None:
        if pm.shape[0] == 1:
            logging.warning('applying RSD filter on single scan, all peaks removed')
        rsd_flag = map(lambda x: not np.isnan(x) and x < rsd_thres, pl_avg.get_attribute(rsd_label, flagged_only=False))
        pl_avg.add_attribute("{}_flag".format(rsd_label), rsd_flag, flagged_only=False, is_flag=True)
    return pl_avg


def average_replicate_peaklists(pls, ppm, min_peaks, rsd_thres=None, block_size=5000, ncpus=None):

    pm = align_peaks(pls, ppm, block_size, ncpus)

    prefix = os.path.commonprefix([p.ID for p in pls])
    merged_id = "{}{}".format(prefix, "_".join(map(str, [p.ID.replace(prefix, "").split(".")[0] for p in pls])))

    pl = pm.to_peaklist(ID=merged_id)
    if "snr" in pm.attributes:
        pl.add_attribute("snr", pm.attr_mean_vector("snr"), on_index=2)

    pl.add_attribute("rsd", pm.rsd(flagged_only=False), on_index=5)
    pl.add_attribute("present_flag", pm.present >= min_peaks, is_flag=True)

    if rsd_thres is not None:
        rsd_flag = map(lambda x: not np.isnan(x) and x < rsd_thres, pl.get_attribute("rsd", flagged_only=False))
        pl.add_attribute("rsd_flag", rsd_flag, flagged_only=False, is_flag=True)

    return pl


def join_peaklists(name, pls):

    def _join_atrtributes(pls):
        attrs_out = collections.OrderedDict()
        for pl in pls:
            for atr in pl.attributes:
                attrs_out.setdefault(atr, []).extend(list(pl.get_attribute(atr, flagged_only=False)))
            if list(pl.attributes) != attrs_out.keys():
                raise IOError("Different attributes")
        return attrs_out

    def _join_meta_data(pl, pls):
        # meta data
        for pl_ in pls:
            for k, v in pl_.metadata.items():
                if k not in pl.metadata:
                    pl.metadata[k] = []
                if v is not None:
                    pl.metadata[k].extend(v)
        return pl

    attrs = _join_atrtributes(pls)
    pl_j = PeakList(ID=name, mz=attrs["mz"], intensity=attrs["intensity"])
    del attrs["mz"], attrs["intensity"]  # default attributes
    for a in attrs:
        pl_j.add_attribute(a, attrs[a], is_flag=(a in pls[0].flag_attributes), flagged_only=False)

    return _join_meta_data(pl_j, pls)
