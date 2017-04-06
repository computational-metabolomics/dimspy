#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Workflow DIMS processing

author(s): Ralf Weber
origin: Nov. 2016
"""

from string import strip
import os
import cPickle as pickle
import time
import zipfile
from ast import literal_eval
from models import PeakList
from models import PeakMatrix


def textfile_to_peaklist(file_name, ID, attr_names_dict=None, flag_names=None, delimiter='\t', has_flag_col=False):
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


def textfiles_to_peaklist(source, filenames):
    peaklists = []
    for fn in filenames:
        if zipfile.is_zipfile(source.encode('string-escape')):
            zf = zipfile.ZipFile(source.encode('string-escape'))
            pl = textfile_to_peaklist(zf.open(fn), ID=os.path.basename(fn), attr_names_dict=None, flag_names=None, delimiter='\t', has_flag_col=True)
        else:
            pl = textfile_to_peaklist(fn, ID=os.path.basename(fn), attr_names_dict=None, flag_names=None, delimiter='\t', has_flag_col=False)
        peaklists.append(pl)
    return peaklists


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
