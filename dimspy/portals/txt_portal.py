#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
txt_portal: PeakList and PeakMatrix plain text portals.

author(s): Albert Zhou, Ralf Weber
origin: Apr. 13, 2017

"""

import os
import logging
import numpy as np
from string import strip
from ast import literal_eval
from dimspy.models.peaklist_tags import PeakList_Tags
from dimspy.models.peaklist import PeakList
from dimspy.models.peak_matrix import PeakMatrix, unmask_all_peakmatrix


def _evalv(vect):
    try:
        ctype = type(literal_eval(vect[0]))
    except (ValueError, SyntaxError):
        ctype = None
    return vect if ctype is None else map(ctype, vect)


# peaklist portals
def save_peaklist_as_txt(pkl, file_name, *args, **kwargs):
    if os.path.isfile(file_name):
        logging.warning('plain text file [%s] already exists, override' % file_name)
    with open(file_name, 'w') as f: f.write(pkl.to_str(*args, **kwargs))


def load_peaklist_from_txt(file_name, ID, delimiter=',', flag_names='auto', has_flag_col=True):
    if not os.path.isfile(file_name):
        raise IOError('plain text file [%s] not exists' % file_name)
    with open(file_name, 'rU') as f:
        rlns = filter(lambda x: x != '', map(strip, f.readlines()))

    dlns = map(lambda x: map(strip, x.split(delimiter)), rlns)
    if any(map(lambda x: len(x) != len(dlns[0]), dlns[1:])):
        raise IOError('data matrix size not match')

    hd, dm = dlns[0], zip(*dlns[1:])
    if has_flag_col:
        hd, dm = hd[:-1], dm[:-1]  # flag_col must be the last one, and discarded
    if len(set(hd)) != len(hd):
        raise IOError('duplicate headers found')

    mzs, ints = np.array(dm[0], dtype=float), np.array(dm[1], dtype=float)  # first two cols must be mz and ints
    pkl = PeakList(ID, mzs, ints)

    flag_names = filter(lambda x: x.endswith('_flag'), hd) if flag_names == 'auto' else \
        [] if flag_names is None else \
            set(flag_names)
    for n, v in zip(hd[2:], dm[2:]): pkl.add_attribute(n, _evalv(v), is_flag=n in flag_names, flagged_only=False)

    return pkl


# peak matrix portals
def save_peak_matrix_as_txt(pm, file_name, *args, **kwargs):
    if os.path.isfile(file_name):
        logging.warning('plain text file [%s] already exists, override' % file_name)
    with open(file_name, 'w') as f:
        with unmask_all_peakmatrix(pm) as m: f.write(m.to_str(*args, **kwargs))


def load_peak_matrix_from_txt(file_name, delimiter='\t', transposed=False, comprehensive=False):
    if not os.path.isfile(file_name):
        raise IOError('plain text file [%s] not exists' % file_name)
    with open(file_name, 'rU') as f:
        rlns = filter(lambda x: x != '', f.readlines())

    dlns = map(lambda x: map(strip, x.split(delimiter)), rlns)
    if any(map(lambda x: len(x) != len(dlns[0]), dlns[1:])):
        raise IOError('data matrix size not match')

    if not transposed: dlns = zip(*dlns)

    def _parseflags():
        fgs = []
        for l, ln in enumerate(zip(*dlns)[5:]):
            if ln[0] == 'flags': break
            fgs += [(ln[0], map(eval, filter(lambda x: x != '', ln[1:])))]
        return fgs
    flgs = _parseflags() if comprehensive else []

    # must refactor if PeakMatrix.to_str changed
    pcol = 6+len(flgs) if comprehensive else 1
    pids = dlns[0][pcol:]

    def _parsetags(tgs):
        for l, ln in enumerate(dlns[2:]):  # line 1 = missing
            if not ln[0].startswith('tags_'): break
            tn, tv = ln[0][5:], ln[pcol:]
            tl = filter(lambda x: x[1] != '', enumerate(_evalv(tv)))
            for i, v in tl: tgs[i].add_tags(v) if tn == 'untyped' else tgs[i].add_tags(**{tn: v})
        return l, tgs
    tnum, tags = 0, [PeakList_Tags() for _ in pids]
    if comprehensive: tnum, tags = _parsetags(tags)

    rlns = zip(*dlns[2 + tnum:])
    mz = np.array([rlns[0]] * len(pids), dtype=float)
    ints = np.array(rlns[pcol:], dtype=float)

    pm = PeakMatrix(pids, tags, [('mz', mz), ('intensity', ints)])
    for fn, fv in flgs: pm.add_flag(fn, fv, flagged_only = False)
    return pm

