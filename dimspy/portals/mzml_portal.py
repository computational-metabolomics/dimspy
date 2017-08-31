#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 1.0.0

"""

import os
import collections
import pymzml
import numpy as np
import zipfile
from dimspy.models.peaklist import PeakList
from dimspy.experiment import mz_range_from_header


class Mzml:
    def __init__(self, filename="", archive=None):
        self.filename = filename
        self.archive = archive

    def run(self):
        assert self.filename.lower().endswith(".mzml") or self.filename.lower().endswith(".mzml.gz") or self.filename.lower().endswith(".zip"), "Incorrect format for mzml parser"
        if self.archive is not None:
            assert zipfile.is_zipfile(self.archive), 'Input file [%s] is not a valid zip archive' % self.archive
            zf = zipfile.ZipFile(self.archive, 'r')
            assert self.filename in zf.namelist(), "{} does not exist in zip file".format(self.filename)
            return pymzml.run.Reader('', file_object=zf.open(self.filename))
        elif self.filename.lower().endswith(".mzml") or self.filename.lower().endswith(".mzml.gz"):
            assert os.path.isfile(self.filename), "{} does not exist".format(self.filename)
            return pymzml.run.Reader(self.filename)
        else:
            return None

    def headers(self):
        return list(set([scan['MS:1000512'] for scan in self.run() if 'MS:1000512' in scan]))

    def peaklist(self, scan_id, function_noise="median"):

        if function_noise not in ["mean", "median", "mad"]:
            raise ValueError("select a function that is available [mean, median, mad]")

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

                pl = PeakList(ID=scan["id"], mz=mzs, intensity=ints,
                              mz_range=mz_range,
                              header=header,
                              ion_injection_time=ion_injection_time,
                              scan_time=scan_time,
                              tic=tic,
                              function_noise=function_noise)

                snr = np.divide(ints, scan.estimatedNoiseLevel(mode=function_noise))
                pl.add_attribute('snr', snr)
                return pl
        return None

    def peaklists(self, scan_ids, function_noise="median"):  # generator

        if function_noise not in ["mean", "median", "mad"]:
            raise ValueError("select a function that is available [mean, median, mad]")

        # somehow i can not access the scans directly when run() uses an open archive object
        # print self.run()[2] fails ... strange
        return [self.peaklist(scan["id"], function_noise) for scan in self.run() if scan["id"] in scan_ids]

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
