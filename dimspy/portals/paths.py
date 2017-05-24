#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
hdf5Portal: PeakList and PeakMatrix HDF5 IO portals.

author(s): Albert Zhou, Ralf Weber
origin: May. 15, 2017

"""


import os
import numpy as np
import zipfile
import h5py
from dimspy.portals import hdf5_portal
from dimspy.models.peaklist import PeakList


def check_paths(tsv, source):
    if tsv is None:
        if type(source) == str:
            source = source.encode('string-escape')
            if os.path.isdir(source):
                filenames = [os.path.join(source, fn) for fn in os.listdir(source) if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]
            elif zipfile.is_zipfile(source):
                with zipfile.ZipFile(source) as zf:
                    assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
                    filenames = [fn for fn in zf.namelist() if fn.lower().endswith(".mzml")]
            elif h5py.is_hdf5(source):
                peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
                filenames = [pl.ID for pl in peaklists]
#
        elif type(source) == list or type(source) == tuple:
            assert isinstance(source[0], PeakList), "Incorrect Objects in list. PeakList class required."
            filenames = [pl.ID for pl in source]
        else:
            pass
#
    elif os.path.isfile(tsv.encode('string-escape')):
        tsv = tsv.encode('string-escape')
        fm = np.genfromtxt(tsv, dtype=None, delimiter="\t", names=True)
        if len(fm.shape) == 0:
            fm = np.array([fm])
        assert fm.dtype.names[0] == "filename" or fm.dtype.names[0] == "sample_id", \
            "Incorrect header for first column. Use filename or sample_id"
#
        filenames = []
        if type(source) == list or type(source) == tuple:
            assert isinstance(source[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
            for fn in fm[fm.dtype.names[0]]:
                assert fn in [pl.ID for pl in source], "{} does not exist in list with Objects".format(fn)
                filenames.append(fn)
        elif type(source.encode('string-escape')) == str:
            source = source.encode('string-escape')
            if os.path.isdir(source):
                l = os.listdir(source)
                for fn in fm[fm.dtype.names[0]]:
                    assert os.path.basename(fn) in l, "{} does not exist in directory provided".format(os.path.basename(fn))
                    filenames.append(os.path.join(source, fn).replace('\\', r'\\'))
#
            elif zipfile.is_zipfile(source.encode('string-escape')):
                with zipfile.ZipFile(source.encode('string-escape')) as zf:
                    assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
                    for fn in fm[fm.dtype.names[0]]:
                        assert fn in zf.namelist(), "{} does not exist in .zip file".format(fn)
                        filenames.append(fn)
#
            elif h5py.is_hdf5(source):
                peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
                filenames = [pl.ID for pl in peaklists]
#
#
            else:
                raise IOError("Can not read and parse {} or {}".format(source, tsv))
    else:
        raise IOError("File {} does not exist".format(tsv))
#
    return filenames
