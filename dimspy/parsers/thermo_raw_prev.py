#!/usr/bin/python
# -*- coding: utf-8 -*-
import tempfile
import os
import subprocess
import collections
import cPickle as pickle
import numpy as np
import argparse
from dimspy.experiment import mz_range_from_header


class ThermoRaw():

    def __init__(self, fname):

        import MSFileReader
        self.run = MSFileReader.ThermoRawfile(fname)

    def headers(self):
        return [str(self.run.GetFilterForScanNum(scan_id)) for scan_id in range(self.run.GetFirstSpectrumNumber(), self.run.GetLastSpectrumNumber() + 1)]

    def headers_scan_ids(self):
        sids = collections.OrderedDict()
        for scan_id in range(self.run.GetFirstSpectrumNumber(), self.run.GetLastSpectrumNumber() + 1):
            sids.setdefault(str(self.run.GetFilterForScanNum(scan_id)), []).append(scan_id)
        return sids

    def peaklist(self, scan_id, mode_noise="noise_packets", out="class"):  # generator
        assert mode_noise in ["noise_packets", "mean", "median", "mad"], "select a method that is available [noise_packets, mean, median, mad]"

        if self.run.IsCentroidScanForScanNum(scan_id):
            msl = self.run.GetMassListRangeFromScanNum(scan_id)
            noiseData = self.run.GetNoiseData(scan_id)
            mzs, ints = msl[0][0], msl[0][1]

            noise = np.interp(mzs, noiseData[0], noiseData[1])
            baseline = np.interp(mzs, noiseData[0], noiseData[2])
        else:
            gld = self.run.GetLabelData(scan_id)
            mzs, ints, baseline, noise = gld[0][0], gld[0][1], gld[0][3], gld[0][4]

        # TODO: sort mz values.
        # Found a raw file where mz files where not in order (20170313_FlyFaeces_DIMS_Lipids_Neg_1116_V5F2c.raw).
        # Correct place to do this?
        mz_ibn = zip(mzs, ints, baseline, noise)
        mz_ibn.sort()
        mzs, ints, baseline, noise = zip(*mz_ibn)

        if mode_noise == "noise_packets":
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2871411
            # Mol Cell Proteomics. 2010 May; 9(5): 754–763.
            # TODO: Negative baseline?
            snr = np.divide(np.subtract(ints, baseline), noise)
        elif mode_noise == "median":
            snr = ints / np.median(ints)
        elif mode_noise == "mean":
            snr = ints / np.mean(ints)
        elif mode_noise == "mad":
            snr = ints / np.median(np.abs(np.subtract(ints, np.median(ints))))

        extra = self.run.GetTrailerExtraForScanNum(scan_id)

        if "Elapsed Scan Time (sec)" in extra:
            scan_time = extra["Elapsed Scan Time (sec)"]  # TODO
        else:
            scan_time = None

        tic = self.run.GetScanHeaderInfoForScanNum(scan_id)["TIC"]
        ion_injection_time = extra["Ion Injection Time (ms)"]
        header = str(self.run.GetFilterForScanNum(scan_id))
        mz_range = mz_range_from_header(header)
        calibration = dict((k.replace("Conversion Parameter ", ""), extra[k]) for k in extra if "Conversion Parameter" in k)
        segment = extra['Scan Segment']
        instrument = str(self.run.GetInstName())

        if out == "dict":
            pl = {"ID": scan_id, "mz": mzs, "intensity": ints,
                           "mz_range": mz_range,
                           "header": header,
                           "ion_injection_time": ion_injection_time,
                           "scan_time": scan_time,
                           "tic": tic,
                           "segment": segment,
                           "calibration": calibration,
                           "instrument": instrument,
                           "mode_noise": mode_noise,
                           "snr": snr,
                           "noise": noise,
                           "baseline": baseline}

        else:
            from dimspy.models.peaklist import PeakList
            pl = PeakList(ID=scan_id, mz=mzs, intensity=ints,
                           mz_range=mz_range,
                           header=header,
                           ion_injection_time=ion_injection_time,
                           scan_time=scan_time,
                           tic=tic,
                           segment=segment,
                           calibraton=calibration,
                           instrument=instrument,
                           mode_noise=mode_noise)

            pl.add_attribute('snr', snr)
            pl.add_attribute('noise', noise)
            pl.add_attribute('baseline', baseline)
        return pl

    def peaklists(self, scan_ids, mode_noise="noise_packets", out="class"):
        assert mode_noise in ["noise_packets", "mean", "median", "mad"], "select a method that is available [noise_packets, mean, median, mad]"
        return [self.peaklist(scan_id, mode_noise=mode_noise, out=out) for scan_id in scan_ids]

    def tics(self):
        return [(self.run.RTFromScanNum(scan_id), self.run.GetScanHeaderInfoForScanNum(scan_id)["TIC"]) for scan_id in range(self.run.GetFirstSpectrumNumber(), self.run.GetLastSpectrumNumber() + 1)]

    def extra_info(self, scan_id):
        return self.run.GetTrailerExtraForScanNum(scan_id)

    def close(self):
        self.run.Close()


class ThermoRawWine():
    def __init__(self, fname):

        fname_pl = tempfile.NamedTemporaryFile()
        ms_file_reader = os.path.abspath(__file__)

        p = subprocess.Popen(
            ["wine", "C:\\python27\\python.exe", ms_file_reader, "--input", fname, "--output", fname_pl.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        print out, err
        self.run = pickle.load(open(fname_pl.name, 'rb'))

    def headers(self):
        return self.run["headers_scan_ids"].keys()

    def headers_scan_ids(self):
        return self.run["headers_scan_ids"]

    def peaklist(self, scan_id, mode_noise="noise_packets", out="class"): # generator
        from dimspy.models.peaklist import PeakList
        d = self.run["peaklists"][scan_id]
        pl = PeakList(ID=d["ID"], mz=d["mz"], intensity=d["intensity"],
            mz_range=d["mz_range"],
            header=d["header"],
            ion_injection_time=d["ion_injection_time"],
            scan_time=d["scan_time"],
            tic=d["tic"],
            segment=d["segment"],
            calibraton=d["calibration"],
            instrument=d["instrument"],
            mode_noise=d["mode_noise"])
        pl.add_attribute('snr', d["snr"])
        pl.add_attribute('noise', d["noise"])
        pl.add_attribute('baseline', d["baseline"])
        return pl

    def peaklists(self, scan_ids, mode_noise="noise_packets", out="class"):
        assert mode_noise in ["noise_packets", "mean", "median", "mad"], "select a method that is available [noise_packets, mean, median, mad]"
        return [self.peaklist(scan_id, mode_noise="noise_packets", out=out) for scan_id in scan_ids]

    def tics(self):
        return self.run["tics"]

    def extra_info(self, scan_id):
        return self.run["extra"][scan_id]

    def close(self):
        return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='This is a python wrapper for MSFileReader via Wine.')
    parser.add_argument('-i', '--input', help='Input instrument .RAW file', required=True)
    parser.add_argument('-o', '--output', help='Output (temporary) .pl file', required=True)
    parser.add_argument('-m', '--mode_noise', default="noise_packets", help='["noise_packets", "mean", "median", "mad"]')
    parser.add_argument('-s', '--scan_ids', help='scan_ids')

    args = parser.parse_args()

    run = ThermoRaw(args.input)
    h_sids = run.headers_scan_ids()
    scan_ids = [sid for h in h_sids for sid in h_sids[h]]

    extras = collections.OrderedDict(zip(scan_ids, [run.extra_info(sid) for sid in scan_ids]))
    pls_dicts = collections.OrderedDict(zip(scan_ids, run.peaklists(scan_ids, args.mode_noise, "dict")))

    d = {"peaklists":pls_dicts,
                 "extra_info":extras,
                 "tics":run.tics(),
                 "headers_scan_ids":h_sids}

    pickle.dump(d, open(args.output, "wb"))
    run.close()