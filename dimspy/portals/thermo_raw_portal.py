#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017-2020 Ralf Weber, Albert Zhou.
#
# This file is part of DIMSpy.
#
# DIMSpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIMSpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DIMSpy.  If not, see <https://www.gnu.org/licenses/>.
#


import collections
import os
from typing import Sequence, Union
import re
import sys

import numpy as np
from ..models.peaklist import PeakList

try:
    import clr
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "ThermoRawFileReader_3_0_41/Libraries"))
    clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
    clr.AddReference('ThermoFisher.CommonCore.Data')
    import ThermoFisher.CommonCore.Data.Business as Business
    import ThermoFisher.CommonCore.RawFileReader as RawFileReader
except ImportError:
    import warnings
    warnings.warn("""
                  DIMSpy requires the Mono framework in order to read and process .raw files. 
                  Install dimspy via conda (highly recommended) to automatically install Mono 
                  (see https://dimspy.readthedocs.io/en/latest/installation.html) or 
                  install Mono from (https://www.mono-project.com). 
                  You can ignore this warning if you use DIMSpy to read and process .mzML files.
                  """)


def mz_range_from_header(h: str) -> list:
    """
    Extract the m/z range from a header or filterstring

    :param h: str
    :return: Sequence[float, float]
    """
    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


class ThermoRaw:
    "ThermoRaw portal"
    def __init__(self, filename):
        """
        Initialise a object interface to a mzML file.

        :param filename: Path to the mzML file

        """
        self.run = RawFileReader.RawFileReaderAdapter.FileFactory(filename)
        self.run.SelectInstrument(Business.Device.MS, 1)
        self.filename = filename
        self.timestamp = self.run.CreationDate

    def headers(self) -> collections.OrderedDict:
        """
        Get all unique header or filter strings and associated scan ids.
        :return: Dictionary
        """
        sids = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            sids.setdefault(str(self.run.GetFilterForScanNumber(scan_id).Filter), []).append(scan_id)
        return sids

    def scan_ids(self) -> collections.OrderedDict:
        """
        Get all scan ids and associated headers or filter strings.
        :return: Dictionary
        """
        sids = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            sids[scan_id] = str(self.run.GetFilterForScanNumber(scan_id).Filter)
        return sids

    def peaklist(self, scan_id, function_noise="noise_packets") -> PeakList:
        """
        Create a peaklist object for a specific scan id.
        :param scan_id: Scan id
        :param function_noise: Function to calculate the noise from each scan. The following options are available:

        * **median** - the median of all peak intensities within a given scan is used as the noise value.

        * **mean** - the unweighted mean average of all peak intensities within a given scan is used as the noise value.

        * **mad (Mean Absolute Deviation)** - the noise value is set as the mean of the absolute differences between peak
          intensities and the mean peak intensity (calculated across all peak intensities within a given scan).

        * **noise_packets** - the noise value is calculated using the proprietary algorithms contained in Thermo Fisher
          Scientific’s msFileReader library. This option should only be applied when you are processing .RAW files.

        :return: PeakList object
        """
        if function_noise not in ["noise_packets", "mean", "median", "mad"]:
            raise ValueError("select a function that is available [noise_packets, mean, median, mad]")

        scan = self.run.GetCentroidStream(scan_id, False)
        if scan.Masses is not None:
            mz_ibn = list(
                zip(scan.Masses, scan.Intensities, scan.Baselines, scan.Noises))  # SignalToNoise not available
            mz_ibn.sort()
            mzs, ints, baseline, noise = list(zip(*mz_ibn))
        else:
            mzs, ints, baseline, noise = [], [], [], []

        if function_noise == "noise_packets" and len(ints) > 0:
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

        ion_injection_time = None
        micro_scans = None
        elapsed_scan_time = None

        extra_values = list(self.run.GetTrailerExtraInformation(scan_id).Values)
        extra_labels = list(self.run.GetTrailerExtraInformation(scan_id).Labels)
        for i, label in enumerate(extra_labels):
            if "Ion Injection Time (ms):" == label:
                ion_injection_time = float(extra_values[i])
            if "Elapsed Scan Time (sec):" == label:
                elapsed_scan_time = float(extra_values[i])
            if "Micro Scan Count:" == label:
                micro_scans = float(extra_values[i])

        scan_time = float(scan_stats.StartTime)
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
                      elapsed_scan_time=elapsed_scan_time,
                      tic=tic,
                      function_noise=function_noise)

        if len(pl.mz) > 0:
            pl.add_attribute('snr', snr)
            pl.add_attribute('noise', noise)
            pl.add_attribute('baseline', baseline)

        return pl

    def peaklists(self, scan_ids, function_noise="noise_packets") -> Sequence[PeakList]:
        """
        Create a list of peaklist objects for each scan id in the list.
        :param scan_ids: List of scan ids

        :param function_noise: Function to calculate the noise from each scan. The following options are available:

        * **median** - the median of all peak intensities within a given scan is used as the noise value.

        * **mean** - the unweighted mean average of all peak intensities within a given scan is used as the noise value.

        * **mad (Mean Absolute Deviation)** - the noise value is set as the mean of the absolute differences between peak
          intensities and the mean peak intensity (calculated across all peak intensities within a given scan).

        * **noise_packets** - the noise value is calculated using the proprietary algorithms contained in Thermo Fisher
          Scientific’s msFileReader library. This option should only be applied when you are processing .RAW files.

        :return: List of PeakList objects
        """
        if function_noise not in ["noise_packets", "mean", "median", "mad"]:
            raise ValueError("select a function that is available [noise_packets, mean, median, mad]")

        return [self.peaklist(scan_id, function_noise=function_noise) for scan_id in scan_ids]

    def tics(self) -> collections.OrderedDict:
        """
        Get all TIC values and associated scan ids
        :return: Dictionary
        """
        tics = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            scan_stats = self.run.GetScanStatsForScanNumber(scan_id)
            tics[scan_id] = scan_stats.TIC
        return tics

    def ion_injection_times(self) -> collections.OrderedDict:
        """
        Get all TIC values and associated scan ids
        :return: Dictionary
        """
        iits = collections.OrderedDict()
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            extra_values = list(self.run.GetTrailerExtraInformation(scan_id).Values)
            extra_labels = list(self.run.GetTrailerExtraInformation(scan_id).Labels)
            for i, label in enumerate(extra_labels):
                if "Ion Injection Time (ms):" == label:
                    iits[scan_id] = float(extra_values[i])
            if scan_id not in iits:
                iits[scan_id] = None
        return iits

    def scan_dependents(self) -> list:
        """
        Get a nested list of scan id pairs. Each pair represents a fragementation event.
        :return: List
        """
        l = []
        for scan_id in range(self.run.RunHeaderEx.FirstSpectrum, self.run.RunHeaderEx.LastSpectrum + 1):
            gsd = self.run.GetScanDependents(scan_id, 5)
            if gsd is not None:
                for i, d in enumerate(gsd.ScanDependentDetailArray):
                    l.append([scan_id, d.ScanIndex])
        return l

    def close(self):
        """
        Close the reader/file object
        :return: None
        """
        self.run.Close()
