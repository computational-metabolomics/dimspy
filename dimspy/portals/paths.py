#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017-2019 Ralf Weber, Albert Zhou.
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


import os

import h5py
import numpy as np
from datetime import datetime
import platform

from ..models.peaklist import PeakList
from ..portals import hdf5_portal
from ..portals.mzml_portal import Mzml
from ..portals.thermo_raw_portal import ThermoRaw


def sort_ms_files_by_timestamp(ps):
    s_files = {}
    for i, fn in enumerate(ps):
        if fn.lower().endswith(".raw"):
            run = ThermoRaw(fn)
            try:
                pattern = "%d/%m/%Y %H:%M:%S"
            except:
                pattern = "%m/%d/%Y %I:%M:%S %p"

        elif fn.lower().endswith(".mzml"):
            run = Mzml(fn)
            pattern = "%Y-%m-%dT%H:%M:%SZ"
        else:
            continue
        s_files[fn] = str(run.timestamp)
        run.close()
    return sorted(s_files.items(), key=lambda x: datetime.strptime(x[1], pattern), reverse=False)


def validate_and_sort_paths(source, tsv):
    """

    :param tsv:
    :param source:
    :return:
    """
    if tsv is None:
        if type(source) == str:
            if os.path.isdir(source):
                filenames = [os.path.join(source, fn) for fn in os.listdir(source) if
                             fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]
                filenames = [fd[0] for fd in sort_ms_files_by_timestamp(filenames)]

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
        fm = np.genfromtxt(tsv, dtype=None, delimiter="\t", names=True, encoding=None)
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
                for fn in source:
                    if not os.path.isfile(fn):
                        raise IOError("[Errno 2] No such file or directory: {}".format(fn))

                for filename in fm[fm.dtype.names[0]]:
                    fns = [os.path.basename(fn) for fn in source]
                    if filename in fns:
                        filenames.append(source[fns.index(filename)])
                    else:
                        raise IOError("{} (row {}) does not exist in source provided".format(filename, list(
                            fm[fm.dtype.names[0]]).index(filename) + 1))

        elif type(source) == str:
            if os.path.isdir(source):
                l = os.listdir(source)
                for fn in fm[fm.dtype.names[0]]:
                    if os.path.basename(fn) not in l:
                        raise IOError("{} does not exist in directory provided".format(os.path.basename(fn)))
                    filenames.append(os.path.join(source, fn))

            elif h5py.is_hdf5(source):
                peaklists = hdf5_portal.load_peaklists_from_hdf5(source)
                filenames = [pl.ID for pl in peaklists]
            elif os.path.isfile(source):
                if source.lower().endswith(".raw") or source.lower().endswith(".mzml"):
                    filenames.append(source)
                else:
                    raise IOError("Incorrect file format, provide .mzml or .raw files: {}".format(source))
            else:
                raise IOError("[Errno 2] No such file or directory: {} or {}".format(source, tsv))
    else:
        raise IOError("[Errno 2] No such file or directory: {} or {}".format(source, tsv))

    return filenames
