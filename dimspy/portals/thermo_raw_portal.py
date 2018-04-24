#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

.. moduleauthor:: Ralf Weber, Albert Zhou

.. versionadded:: 1.0.0

"""

import sys
import os
import re
import collections
import numpy as np
from dimspy.models.peaklist import PeakList
import clr
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "ThermoRawFileReader_3_0_41/Libraries"))
clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
clr.AddReference('ThermoFisher.CommonCore.Data')
import ThermoFisher.CommonCore.Data.Business as Business
import ThermoFisher.CommonCore.RawFileReader as RawFileReader


def mz_range_from_header(h):
    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


class ThermoRaw:

    def __init__(self, filename):
        self.run = RawFileReader.RawFileReaderAdapter.FileFactory(filename)
        self.run.SelectInstrument(Business.Device.MS, 1)
        self.filename = filename

    def headers(self):

        sids = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            sids.setdefault(str(self.run.GetFilterForScanNumber(scan_id).Filter), []).append(scan_id)
        return sids

    def scan_ids(self):

        sids = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            sids[scan_id] = str(self.run.GetFilterForScanNumber(scan_id).Filter)
        return sids

    def peaklist(self, scan_id, function_noise="noise_packets"):

        if function_noise not in ["noise_packets", "mean", "median", "mad"]:
            raise ValueError("select a function that is available [noise_packets, mean, median, mad]")

        scan = self.run.GetCentroidStream(scan_id, False)
        if scan.Masses is not None:
            mz_ibn = zip(scan.Masses, scan.Intensities, scan.Baselines, scan.Noises)  # SignalToNoise not available
            mz_ibn.sort()
            mzs, ints, baseline, noise = zip(*mz_ibn)
        else:
            mzs, ints, baseline, noise = [], [], [], []

        if function_noise == "noise_packets":
            snr = [p.SignalToNoise for p in scan.GetCentroids()]
        elif function_noise == "median" and len(ints) > 0:
            snr = ints / np.median(ints)
        elif function_noise == "mean" and len(ints) > 0:
            snr = ints / np.mean(ints)
        elif function_noise == "mad" and len(ints) > 0:
            snr = ints / np.median(np.abs(np.subtract(ints, np.median(ints))))
        else:
            snr = []

        scan_stats = self.run.GetScanStatsForScanNumber(scan_id)

        extra_values = list(self.run.GetTrailerExtraInformation(scan_id).Values)
        extra_labels = list(self.run.GetTrailerExtraInformation(scan_id).Labels)
        for i, label in enumerate(extra_labels):
            if "Ion Injection Time (ms):" == label:
                ion_injection_time = extra_values[i]
            else:
                ion_injection_time = None
            if "Elapsed Scan Time (sec):" == label:
                scan_time = extra_values[i]
            else:
                scan_time = None
            if "Micro Scan Count:" == label:
                micro_scans = extra_values[i]
            else:
                micro_scans = None

        tic = scan_stats.TIC
        segment = scan_stats.SegmentNumber
        header = str(self.run.GetScanEventStringForScanNumber(scan_id))
        ms_level = header.count("@") + 1

        pl = PeakList(ID=scan_id, mz=mzs, intensity=ints,
                      mz_range=mz_range_from_header(header),
                      header=header,
                      ms_level=ms_level,
                      micro_scans=micro_scans,
                      segment=segment,
                      ion_injection_time=ion_injection_time,
                      scan_time=scan_time,
                      tic=tic,
                      function_noise=function_noise)

        if len(pl.mz) > 0:
            pl.add_attribute('snr', snr)
            pl.add_attribute('noise', noise)
            pl.add_attribute('baseline', baseline)

        return pl

    def peaklists(self, scan_ids, function_noise="noise_packets"):
        if function_noise not in ["noise_packets", "mean", "median", "mad"]:
            raise ValueError("select a function that is available [noise_packets, mean, median, mad]")

        return [self.peaklist(scan_id, function_noise=function_noise) for scan_id in scan_ids]

    def tics(self):
        tics = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            scan_stats = self.run.GetScanStatsForScanNumber(scan_id)
            tics[scan_id].append(scan_stats.TIC)
        return tics

    def injection_times(self):
        injection_times = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            extra_values = list(self.run.GetTrailerExtraInformation(scan_id).Values)
            extra_labels = list(self.run.GetTrailerExtraInformation(scan_id).Labels)
            for i, label in enumerate(extra_labels):
                if "Ion Injection Time (ms):" == label:
                    injection_times[scan_id] = float(extra_values[i])
            if scan_id not in injection_times:
                injection_times[scan_id] = None
        return injection_times

    def scan_dependents(self):
        l = []
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            gsd = self.run.GetScanDependents(scan_id, 5)
            if gsd is not None:
                for i, d in enumerate(gsd.ScanDependentDetailArray):
                    l.append([scan_id, d.ScanIndex])
        return l
