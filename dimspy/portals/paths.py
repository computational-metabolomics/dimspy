#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 1.0.0

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
            if os.path.isdir(source):
                filenames = [os.path.join(source, fn) for fn in os.listdir(source) if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]
            elif zipfile.is_zipfile(source):
                with zipfile.ZipFile(source) as zf:
                    if len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) > 0:
                        raise IOError("Archive with *.raw files not yet supported. Convert to mzML")
                    filenames = [fn for fn in zf.namelist() if fn.lower().endswith(".mzml")]
            elif h5py.is_hdf5(source):
                peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
                filenames = [os.path.join(os.path.abspath(os.path.dirname(source)), pl.ID) for pl in peaklists]
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

    elif os.path.isfile(tsv):
        fm = np.genfromtxt(tsv, dtype=None, delimiter="\t", names=True)
        if len(fm.shape) == 0:
            fm = np.array([fm])
        if fm.dtype.names[0] != "filename" and fm.dtype.names[0] != "sample_id":
            raise IOError("Incorrect header for first column. Use filename or sample_id")

        filenames = []
        if type(source) == list or type(source) == tuple:
            if isinstance(source[0], PeakList):
                for filename in fm[fm.dtype.names[0]]:
                    if filename in [pl.ID for pl in source]:
                        filenames.append(filename)
                    else:
                        raise IOError("{} does not exist in list with Peaklist objects".format(filename))
            else:
                for filename in fm[fm.dtype.names[0]]:
                    if filename not in [os.path.basename(fn) for fn in source]:
                        raise IOError("{} (row {}) does not exist in source provided".format(filename, list(fm[fm.dtype.names[0]]).index(filename)+1))
                for fn in source:
                    if os.path.isfile(fn):
                        filenames.append(fn)
                    else:
                        raise IOError("[Errno 2] No such file or directory: {}".format(fn))
  
        elif type(source) == str:
            if os.path.isdir(source):
                l = os.listdir(source)
                for fn in fm[fm.dtype.names[0]]:
                    if os.path.basename(fn) not in l:
                        raise IOError("{} does not exist in directory provided".format(os.path.basename(fn)))
                    filenames.append(os.path.join(source, fn))

            elif zipfile.is_zipfile(source):
                with zipfile.ZipFile(source) as zf:
                    if len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) > 0:
                        raise IOError("Archive with *.raw files not yet supported. Convert to mzML")
                    for fn in fm[fm.dtype.names[0]]:
                        if fn not in zf.namelist():
                            raise IOError("{} does not exist in .zip file".format(fn))
                        filenames.append(fn)

            elif h5py.is_hdf5(source):
                peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
                filenames = [pl.ID for pl in peaklists]
            else:
                raise IOError("[Errno 2] No such file or directory: {} or {}".format(source, tsv))
    else:
        raise IOError("[Errno 2] No such file or directory: {} or {}".format(source, tsv))

    return filenames
