#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import collections
import pymzml
import numpy as np
import zipfile
from dimspy.models.peaklist import PeakList


def mz_range_from_header(h):
    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


class Mzml():
    def __init__(self, fname="", archive=None):
        self.fname = fname
        self.archive = archive

    def run(self):
        if self.archive is not None:
            zf = zipfile.ZipFile(self.archive, 'r')
            return pymzml.run.Reader('', file_object=zf.open(self.fname))
        elif os.path.isfile(self.fname) and (self.fname.lower().endswith(".mzml") or self.fname.lower().endswith(".mzml.gz")):
            return pymzml.run.Reader(self.fname)
        else:
            return None

    def headers(self):
        return list(set([scan['MS:1000512'] for scan in self.run() if 'MS:1000512' in scan]))

    def peaklist(self, scan_id, mode="median"):

        assert mode in ["mean", "median", "mad"], "select a method that is available [msfilereader, mean, median, mad]"
        for scan in self.run():
            if scan["id"] == scan_id:

                mzs, ints = zip(*scan.peaks)

                scan_time = scan["MS:1000016"]
                tic = scan["total ion current"]
                if "MS:1000927" in scan:
                    ion_injection_time = scan["MS:1000927"]
                else:
                    ion_injection_time = None
                header = scan['MS:1000512']
                mz_range = mz_range_from_header(header)

                pkl = PeakList(ID=scan["id"], mz=mzs, intensity=ints,
                               mz_range=mz_range,
                               header=header,
                               ion_injection_time=ion_injection_time,
                               scan_time=scan_time,
                               tic=tic,
                               segment=None,
                               calibraton=None,
                               instrument=None)

                snr = np.divide(ints, scan.estimatedNoiseLevel(mode=mode))
                pkl.add_attribute('snr', snr)
                return pkl
        return None

    def peaklists(self, scan_ids, mode="median"):  # generator

        assert mode in ["mean", "median", "mad"], "select a method that is available"
        # somehow i can not access the scans directly when run() uses an open archive object
        # print self.run()[2] fails ... strange
        return [self.peaklist(scan["id"], mode) for scan in self.run() if scan["id"] in scan_ids]

    def headers_scan_ids(self, n=None):
        h_sids = collections.OrderedDict()
        for scan in self.run():
            if 'MS:1000512' in scan:
                if n is None:
                    h_sids.setdefault(scan['MS:1000512'], []).append(scan['id'])
                elif len(h_sids[scan['MS:1000512']]) < n:
                    h_sids.setdefault(scan['MS:1000512'], []).append(scan['id'])
        return h_sids

    def tics(self):
        # somehow i can not access the scans directly when run() uses an open archive object
        # print self.run()[2]
        for scan in self.run():
            if scan["id"] == "TIC":
                return zip(*scan.peaks)
        return None

    def extra_info(self, scan_id):
        # Not available
        return None
