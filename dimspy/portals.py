#!/usr/bin/python
# -*- coding: utf-8 -*-
from string import strip
import os
import cPickle as pickle
import time
import zipfile
from ast import literal_eval
import numpy as np
from models import PeakList
from models import PeakMatrix
import h5py

def check_paths(tsv, source):

    if tsv is None:
        if type(source) == str:
            source.encode('string-escape')
            if os.path.isdir(source):
                filenames = [os.path.join(source, fn) for fn in os.listdir(source) if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]
            elif zipfile.is_zipfile(source):
                with zipfile.ZipFile(source) as zf:
                    assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
                    filenames = [fn for fn in zf.namelist() if fn.lower().endswith(".mzml")]
            elif h5py.is_hdf5(source):
                with h5py.file(source) as h5:
                    # filenames = h5.
                    pass
        elif type(source) == list:
            assert type(source[0]) == PeakList, "Incorrect Objects in list. PeakList class required."
            filenames = [pl.ID for pl in source]
        else:
            pass

    elif os.path.isfile(tsv.encode('string-escape')):
        tsv = tsv.encode('string-escape')
        fm = np.genfromtxt(tsv, dtype=None, delimiter="\t", names=True)
        if len(fm.shape) == 0:  # TODO: Added to check if filelist has a single row
            fm = np.array([fm])
        assert fm.dtype.names[0] == "filename" or fm.dtype.names[0] == "sample_id", \
            "Incorrect header for first column. Use filename or sample_id"

        filenames = []
        if type(source) == list:
            assert type(source[0]) == PeakList, "Incorrect Objects in list. Peaklist Object required."
            for fn in fm[fm.dtype.names[0]]:
                assert fn in [pl.ID for pl in source], "{} does not exist in list with Objects".format(fn)
                filenames.append(fn)
        elif type(source.encode('string-escape')) == str:
            if os.path.isdir(source):
                l = os.listdir(source.encode('string-escape'))
                for fn in fm[fm.dtype.names[0]]:
                    assert os.path.basename(fn) in l, "{} does not exist in directory provided".format(os.path.basename(fn))
                    filenames.append(os.path.join(source, fn).replace('\\', r'\\'))

            elif zipfile.is_zipfile(source.encode('string-escape')):
                with zipfile.ZipFile(source) as zf:
                    assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
                    for fn in fm[fm.dtype.names[0]]:
                        assert fn in zf.namelist(), "{} does not exist in .zip file".format(fn)
                        filenames.append(fn)
            elif h5py.is_hdf5(source.encode('string-escape')):
                with h5py.file(source.encode('string-escape')) as h5:
                    # filenames = h5.
                    pass

            else:
                raise IOError("Can not read and parse {} or {}".format(source, tsv))
    else:
        raise IOError("File {} does not exist".format(tsv))

    return filenames


def load_peaklists(source):

    if type(source) == str:
        source = source.encode('string-escape')
        if h5py.is_hdf5(source):
            peaklists = hdf5_to_peaklists(source)
        elif zipfile.is_zipfile(source):
            zf = zipfile.ZipFile(source)
            filenames = zf.namelist()
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [text_to_peaklist(zf.open(fn), ID=os.path.basename(fn), has_flag_col=True) for fn in filenames]
        elif os.path.isdir(source):
            filenames = os.listdir(source)
            assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
                "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
            peaklists = [text_to_peaklist(os.path.join(source, fn), ID=os.path.basename(fn), has_flag_col=False) for fn in filenames]
        else:
            raise TypeError("Incorrect format. Process .mzML and .raw files first using the 'process scans' function")
    elif type(source) == list:
        if not isinstance(source[0], PeakList):
            raise TypeError("List has incorrect format. PeakList objects required.")
        else:
            peaklists = source
    else:
        raise IOError("Inccorrect input: list with peaklist objects or path")

    return peaklists


def text_to_peaklist(file_name, ID, attr_names_dict=None, flag_names=None, delimiter='\t', has_flag_col=False):
    if hasattr(file_name, 'read'):
        rlns = filter(lambda x: x != '', map(strip, file_name.readlines()))
    else:
        with open(file_name, 'rU') as f:
            rlns = filter(lambda x: x != '', map(strip, f.readlines()))

    dlns = map(lambda x: map(strip, x.split(delimiter)), rlns)

    if has_flag_col: dlns = zip(*zip(*dlns)[:-1]) # flag_col should be the last one
    if attr_names_dict is None: attr_names_dict = {n:n for n in set(dlns[0])}
    if flag_names is None: flag_names = filter(lambda x: x.endswith('_flag'), set(dlns[0]))

    assert all(map(lambda x: x in attr_names_dict.values(), ('mz', 'intensity'))), 'required attribute(s) not exist'
    assert all(map(lambda x: x in attr_names_dict.keys(), flag_names)), 'flag attribute(s) not found'
    assert all(map(lambda x: len(x) == len(dlns[0]), dlns[1:])), 'data matrix size not match'
    assert set(attr_names_dict.keys()) == set(dlns[0]), 'attribution names dict not match'

    def _eval(val):
        try: return literal_eval(val)
        except (ValueError, SyntaxError): return val

    ddct = dict(map(lambda x: (attr_names_dict[x[0]], (map(_eval, x[1:]), x[0] in flag_names)), zip(*dlns)))

    pl = PeakList(ID, ddct['mz'][0], ddct['intensity'][0])
    for k, (v, f) in filter(lambda x: x[0] not in ('mz', 'intensity'), ddct.items()):
        pl.add_attribute(k, v, is_flag=f, flagged_only=False)
    return pl

def to_readable(pickle_file, path_out, separator, transpose=False):
    assert os.path.isfile(pickle_file), "Pickle file does not exist"
    assert separator in ["tab", "comma"], "Incorrect separator [tab, comma]"
    seps = {"comma": ",", "tab": "\t"}
    with open(pickle_file, "rb") as pkl:
        pls = pickle.load(pkl)
        if type(pls) == list:
            assert isinstance(pls[0], PeakList), "Not compatible with {}".format(type(pls[0]))
            for pl in pls:
                with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                    pk_out.write(pl.to_str(seps[separator]))
                time.sleep(1)

        elif isinstance(pls, PeakMatrix):
            assert os.path.isfile(path_out), "Provide filename for peak matrix"
            with open(os.path.join(path_out), "w") as pk_out:
                pk_out.write(pls.to_str(seps[separator], transpose))
    return


def hdf5_to_peaklists(f):
    pass


def hdf5_to_peak_matrix(f):
    pass


def peaklists_to_hdf5(peaklists):
    pass


def peak_matrix_to_hdf5(pm):
    pass
