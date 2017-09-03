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
    """
    Extract a list of headers / .
    :rtype: list
    """
    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


class ThermoRaw:
    """
    Extract a list of headers / .
    :rtype: list
    """
    def __init__(self, filename):
        self.run = RawFileReader.RawFileReaderAdapter.FileFactory(filename)
        self.run.SelectInstrument(Business.Device.MS, 1)

    def headers(self):
        """
        Extract a particular scan from a *.raw file and return a PeakList objects
        :rtype: dict
        """
        sids = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            sids.setdefault(str(self.run.GetFilterForScanNumber(scan_id).Filter), []).append(scan_id)
        return sids

    def scan_ids(self):
        """
        Extract a particular scan from a *.raw file and return a PeakList objects
        :rtype: dict
        """
        sids = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            sids[scan_id] = str(self.run.GetFilterForScanNumber(scan_id).Filter)
        return sids

    def peaklist(self, scan_id, function_noise="noise_packets"):
        """
        Extract a particular scan from a *.raw file and return a PeakList objects

        :param scan_ids:
        :rtype: list
        """
        if function_noise not in ["noise_packets", "mean", "median", "mad"]:
            raise ValueError("select a function that is available [noise_packets, mean, median, mad]")

        scan = self.run.GetCentroidStream(scan_id, False)

        mz_ibn = zip(scan.Masses, scan.Intensities, scan.Baselines, scan.Noises)  # SignalToNoise not available
        mz_ibn.sort()
        mzs, ints, baseline, noise = zip(*mz_ibn)

        if function_noise == "noise_packets":
            snr = [p.SignalToNoise for p in scan.GetCentroids()]
        elif function_noise == "median":
            snr = ints / np.median(ints)
        elif function_noise == "mean":
            snr = ints / np.mean(ints)
        elif function_noise == "mad":
            snr = ints / np.median(np.abs(np.subtract(ints, np.median(ints))))

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

        pl.add_attribute('snr', snr)
        pl.add_attribute('noise', noise)
        pl.add_attribute('baseline', baseline)
        return pl

    def peaklists(self, scan_ids, function_noise="noise_packets"):
        """
        Extract the scans from a *.raw file and return a list of PeakList objects

        :param scan_ids:
        :rtype: list

        """
        if function_noise not in ["noise_packets", "mean", "median", "mad"]:
            raise ValueError("select a function that is available [noise_packets, mean, median, mad]")

        return [self.peaklist(scan_id, function_noise=function_noise) for scan_id in scan_ids]

    def tics(self):
        # somehow i can not access the scans directly when run() uses an open archive object
        # print self.run()[2]
        tics = []
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            scan_stats = self.run.GetScanStatsForScanNumber(scan_id)
            tics.append(scan_stats.TIC)
        return tics

    def scan_dependents(self):
        l = []
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            gsd = self.run.GetScanDependents(scan_id, 5)
            if gsd is not None:
                for i, d in enumerate(gsd.ScanDependentDetailArray):
                    print scan_id, self.run.GetFilterForScanNumber(scan_id).Filter, d.ScanIndex, d.FilterString
                    l.append([scan_id, d.ScanIndex])
        return l
