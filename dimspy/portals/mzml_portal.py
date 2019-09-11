#!/usr/bin/python
# -*- coding: utf-8 -*-

import collections
import os

import numpy as np
import pymzml

from dimspy.experiment import mz_range_from_header
from dimspy.models.peaklist import PeakList


class Mzml:
    def __init__(self, filename="", **kwargs):
        self.filename = filename

        if not os.path.isfile(self.filename):
            raise IOError("{} does not exist".format(self.filename))

        if not self.filename.lower().endswith(".mzml") and not self.filename.lower().endswith(".mzml.gz"):
            raise IOError('Incorrect file format for mzML parser')

        if "ms_precisions" in kwargs:
            self.ms_precisions = kwargs["ms_precisions"]
        else:
            self.ms_precisions = dict(zip(range(3, 11), 8 * [5e-6]))

        self._sids = self._scan_ids()

        self.run = pymzml.run.Reader(self.filename)
        self.run.ms_precisions.update(self.ms_precisions)

    def headers(self):
        """

        :param n:
        :return:
        """
        h_sids = collections.OrderedDict()
        for scan_id in self._sids:
            if 'MS:1000512' in self.run[scan_id]:
                h_sids.setdefault(self.run[scan_id]['MS:1000512'], []).append(scan_id)
        return h_sids

    def _scan_ids(self):
        """

        :return:
        """
        sids_h = collections.OrderedDict()
        run = pymzml.run.Reader(self.filename)
        run.ms_precisions.update(self.ms_precisions)
        for scan in run:
            if 'MS:1000512' in scan:
                sids_h[scan.ID] = str(scan['MS:1000512'])
            else:
                sids_h[scan.ID] = None
        run.close()
        return sids_h

    def scan_ids(self):
        """

        :return:
        """
        return self._sids

    def peaklist(self, scan_id, function_noise="median"):
        """

        :param scan_id:
        :param function_noise:
        :return:
        """
        if function_noise not in ["mean", "median", "mad"]:
            raise ValueError("select a function that is available [mean, median, mad]")

        scan = self.run[scan_id]
        peaks = scan.peaks("raw")
        if len(peaks) > 0:
            mzs, ints = list(zip(*peaks))
        else:
            mzs, ints = [], []

        scan_time = scan["MS:1000016"]
        tic = scan["total ion current"]
        if "MS:1000927" in scan:
            ion_injection_time = scan["MS:1000927"]
        else:
            ion_injection_time = None
        header = scan['MS:1000512']
        mz_range = mz_range_from_header(header)
        ms_level = scan['ms level']
        pl = PeakList(ID=scan.ID, mz=mzs, intensity=ints,
                      mz_range=mz_range,
                      header=header,
                      ms_level=ms_level,
                      ion_injection_time=ion_injection_time,
                      scan_time=scan_time,
                      tic=tic,
                      function_noise=function_noise)
        snr = np.divide(ints, scan.estimated_noise_level(mode=function_noise))
        pl.add_attribute('snr', snr)
        return pl

    def peaklists(self, scan_ids, function_noise="median"):
        """

        :param scan_ids:
        :param function_noise:
        :return:
        """
        if function_noise not in ["mean", "median", "mad"]:
            raise ValueError("select a function that is available [mean, median, mad]")

        return [self.peaklist(scan_id, function_noise) for scan_id in scan_ids if scan_id in self._sids]

    def tics(self):
        """

        :return:
        """
        tic_values = collections.OrderedDict()
        for scan_id in self._sids:
            tic_values[scan_id] = self.run[scan_id].TIC
        return tic_values

    def ion_injection_times(self):
        """

        :return:
        """
        iits = collections.OrderedDict()
        for scan_id in self._sids:
            scan = self.run[scan_id]
            if "MS:1000927" in scan:
                iits[scan_id] = scan["MS:1000927"]
            else:
                iits[scan_id] = None
        return iits

    def scan_dependents(self):
        """

        :return:
        """
        l = []
        for scan_id in self._sids:
            scan = self.run[scan_id]
            if scan.selected_precursors:
                precursor = scan.element.find("./{}precursorList/{}precursor".format(scan.ns, scan.ns))
                l.append([int(precursor.get("spectrumRef").split("scan=")[1]), scan.ID])
        return l

    def close(self):
        """

        :return:
        """
        self.run.close()
