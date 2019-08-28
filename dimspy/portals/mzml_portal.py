#!/usr/bin/python 
# -*- coding: utf-8 -*-

"""

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 1.0.0

"""


import os
import collections
import numpy as np
import pymzml
from dimspy.models.peaklist import PeakList
from dimspy.experiment import mz_range_from_header


class Mzml:
    def __init__(self, filename="", preload=True):
        self.filename = filename

        if not os.path.isfile(self.filename):
            raise IOError("{} does not exist".format(self.filename))

        if not self.filename.lower().endswith(".mzml") and not self.filename.lower().endswith(".mzml.gz"):
            raise IOError('Incorrect file format for mzML parser')

        # run = pymzml.run.Reader(self.filename)
        # run.info["file_object"].close()

    def headers(self, n=None):
        h_sids = collections.OrderedDict()
        run = pymzml.run.Reader(self.filename)
        for scan in run:
            if 'MS:1000512' in scan:
                h_sids.setdefault(scan['MS:1000512'], []).append(scan['id'])
        run.info["file_object"].close()
        return h_sids

    def scan_ids(self):
        h_sids = collections.OrderedDict()
        run = pymzml.run.Reader(self.filename)
        for scan in run:
            if 'MS:1000512' in scan:
                h_sids[scan['id']] = str(scan['MS:1000512'])
        run.info["file_object"].close()
        return h_sids

    def peaklist(self, scan_id, function_noise="median"):

        if function_noise not in ["mean", "median", "mad"]:
            raise ValueError("select a function that is available [mean, median, mad]")

        run = pymzml.run.Reader(self.filename)
        for scan in run:
            if scan["id"] == scan_id:
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
                pl = PeakList(ID=scan["id"], mz=mzs, intensity=ints,
                              mz_range=mz_range,
                              header=header,
                              ms_level=ms_level,
                              ion_injection_time=ion_injection_time,
                              scan_time=scan_time,
                              tic=tic,
                              function_noise=function_noise)
                snr = np.divide(ints, scan.estimated_noise_level(mode=function_noise))
                pl.add_attribute('snr', snr)
                run.info["file_object"].close()
                return pl
        return None

    def peaklists(self, scan_ids, function_noise="median"):  # generator

        if function_noise not in ["mean", "median", "mad"]:
            raise ValueError("select a function that is available [mean, median, mad]")
        run = pymzml.run.Reader(self.filename)
        pls = [self.peaklist(scan["id"], function_noise) for scan in run if scan["id"] in scan_ids]
        run.info["file_object"].close()
        return pls

    def tics(self):
        run = pymzml.run.Reader(self.filename)
        for scan in run:
            if scan["id"] == "TIC":
                hit =  zip(*scan.peaks("raw"))[1]
                run.info["file_object"].close()
                return hit
        run.info["file_object"].close()
        return

    def injection_times(self):
        injection_times = {}
        run = pymzml.run.Reader(self.filename)
        for scan in run:
            injection_times[scan['id']] = None
            for element in scan.xmlTree:
                if "MS:1000927" == element.get('accession'):
                    injection_times[scan['id']] = float(element.get("value"))
                    break
            if scan['id'] not in injection_times:
                injection_times[scan['id']] = None
        run.info["file_object"].close()
        return injection_times

    def scan_dependents(self):
        l = []
        run = pymzml.run.Reader(self.filename)
        for scan in run:
            if type(scan["id"]) == int:
                scan_id = scan["id"]
                if "precursors" in list(scan.keys()):
                    spectrum_ref = None
                    for element in scan.xmlTree:
                        for e in list(element.items()):
                            if e[0] == 'spectrumRef':
                                spectrum_ref = int(e[1].split("scan=")[1])
                    if spectrum_ref is not None:
                        l.append([spectrum_ref, scan_id])
        run.info["file_object"].close()
        return l
