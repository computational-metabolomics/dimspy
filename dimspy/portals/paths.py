#!/usr/bin/python
# -*- coding: utf-8 -*-

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
            elif os.path.isfile(source):
                if source.lower().endswith(".raw") or source.lower().endswith(".mzml"):
                    filenames = [source]
                else:
                    raise IOError("Incorrect file format, provide .mzml or .raw files: {}".format(source))
            else:
                raise IOError("[Errno 2] No such file or directory: {}".format(source))

        elif type(source) == list or type(source) == tuple:
            if isinstance(source[0], PeakList):
                filenames = [pl.ID for pl in source]
            else:
                filenames = []
                for fn in source:
                    if os.path.isfile(fn):
                        if fn.lower().endswith(".raw") or fn.lower().endswith(".mzml"):
                            filenames.append(fn)
                        else:
                            raise IOError("Incorrect file format, provide .mzml or .raw files: {}".format(source))
                    else:
                        raise IOError("[Errno 2] No such file or directory: {}".format(source))
        else:
            raise IOError("[Errno 2] No such file or directory: {}".format(source))

    elif os.path.isfile(tsv.encode('string-escape')):
        tsv = tsv.encode('string-escape')
        fm = np.genfromtxt(tsv, dtype=None, delimiter="\t", names=True)
        if len(fm.shape) == 0:
            fm = np.array([fm])
        assert fm.dtype.names[0] == "filename" or fm.dtype.names[0] == "sample_id", \
            "Incorrect header for first column. Use filename or sample_id"

        filenames = []
        if type(source) == list or type(source) == tuple:
            if isinstance(source[0], PeakList):
                for fname in fm[fm.dtype.names[0]]:
                    if fname in [pl.ID for pl in source]:
            	        filenames.append(fname)
                    else:
                        raise IOError("{} does not exist in list with Peaklist objects".format(fn))
            else:
                for fname in fm[fm.dtype.names[0]]:
                    if fname not in [os.path.basename(fn) for fn in source]:
                        raise IOError("{} (row {}) does not exist in source provided".format(fname, list(fm[fm.dtype.names[0]]).index(fname)+1))
                for fn in source:
                    if os.path.isfile(fn):
                        filenames.append(fn)
                    else:
                        raise IOError("[Errno 2] No such file or directory: {}".format(fn))
  
        elif type(source.encode('string-escape')) == str:
            source = source.encode('string-escape')
            if os.path.isdir(source):
                l = os.listdir(source)
                for fn in fm[fm.dtype.names[0]]:
                    assert os.path.basename(fn) in l, "{} does not exist in directory provided".format(os.path.basename(fn))
                    filenames.append(os.path.join(source, fn).replace('\\', r'\\'))

            elif zipfile.is_zipfile(source.encode('string-escape')):
                with zipfile.ZipFile(source.encode('string-escape')) as zf:
                    assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
                    for fn in fm[fm.dtype.names[0]]:
                        assert fn in zf.namelist(), "{} does not exist in .zip file".format(fn)
                        filenames.append(fn)

            elif h5py.is_hdf5(source):
                peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
                filenames = [pl.ID for pl in peaklists]
            else:
                raise IOError("[Errno 2] No such file or directory: {} or {}".format(source, tsv))
    else:
        raise IOError("[Errno 2] No such file or directory: {} or {}".format(source, tsv))

    return filenames
