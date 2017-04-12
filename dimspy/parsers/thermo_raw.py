#!/usr/bin/python
# -*- coding: utf-8 -*-
import clr
import sys
sys.path.append("ThermoRawFileReader_3_0_41/Libraries")
clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
clr.AddReference('ThermoFisher.CommonCore.Data')
import ThermoFisher.CommonCore.Data.Business as Business
import re
import collections
import numpy as np
from dimspy.models.peaklist import PeakList


def mz_range_from_header(h):
    return [float(m) for m in re.findall(r'([\w\.-]+)-([\w\.-]+)', h)[0]]


class ThermoRaw():

    def __init__(self, fname):

        self.run = Business.RawFileReaderFactory.ReadFile(fname)
        self.run.SelectInstrument(Business.Device.MS, 1)

    def headers(self):
        return [str(self.run.GetScanEventStringForScanNumber(scan_id)) for scan_id in range(1, self.run.ScanEvents.ScanEvents + 1)]

    def headers_scan_ids(self):
        sids = collections.OrderedDict()
        for scan_id in range(1, self.run.ScanEvents.ScanEvents + 1):
            sids.setdefault(str(self.run.GetScanEventStringForScanNumber(scan_id)), []).append(scan_id)
        return sids

    def peaklist(self, scan_id, mode_noise="noise_packets"):  # generator
        assert mode_noise in ["noise_packets", "mean", "median", "mad"], "select a method that is available [noise_packets, mean, median, mad]"

        scan = self.run.GetCentroidStream(scan_id, True)

        mz_ibn = zip(scan.Masses, scan.Intensities, scan.Baselines, scan.Noises)  # SignalToNoise not available
        mz_ibn.sort()
        mzs, ints, baseline, noise = zip(*mz_ibn)

        if mode_noise == "noise_packets":
            snr = [p.SignalToNoise for p in scan.GetCentroids()]
        elif mode_noise == "median":
            snr = ints / np.median(ints)
        elif mode_noise == "mean":
            snr = ints / np.mean(ints)
        elif mode_noise == "mad":
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
        instrument = self.run.GetInstrumentData().Name

        pl = PeakList(ID=scan_id, mz=mzs, intensity=ints,
                           mz_range=mz_range_from_header(header),
                           header=header,
                           micro_scans=micro_scans,
                           ion_injection_time=ion_injection_time,
                           scan_time=scan_time,
                           tic=tic,
                           segment=segment,
                           instrument=instrument,
                           mode_noise=mode_noise)

        pl.add_attribute('snr', snr)
        pl.add_attribute('noise', noise)
        pl.add_attribute('baseline', baseline)
        return pl

    def peaklists(self, scan_ids, mode_noise="noise_packets"):
        assert mode_noise in ["noise_packets", "mean", "median", "mad"], "select a method that is available [noise_packets, mean, median, mad]"
        return [self.peaklist(scan_id, mode_noise=mode_noise) for scan_id in scan_ids]


# testing
if __name__ == '__main__':
    TR = ThermoRaw("../../tests/data/raw/centroid.raw")
    with open("../../tests/data/raw/centroid_scan_03.txt", "w") as out:
        out.write(TR.peaklist(3).to_str("\t"))
    with open("../../tests/data/raw/centroid_scan_43.txt", "w") as out:
        out.write(TR.peaklist(43).to_str("\t"))

    TR = ThermoRaw("../../tests/data/raw/profile.raw")
    with open("../../tests/data/raw/profile_scan_03.txt", "w") as out:
        out.write(TR.peaklist(3).to_str("\t"))
    with open("../../tests/data/raw/profile_scan_43.txt", "w") as out:
        out.write(TR.peaklist(43).to_str("\t"))
