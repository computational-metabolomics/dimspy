#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
The PeakList and PeakMatrix plain text portals.

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 1.0.0

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
def save_peaklist_as_txt(pkl, filename, *args, **kwargs):
    """
    Saves a peaklist object to a plain text file.

    :param pkl: the target peaklist object
    :param filename: path to a new text file
    :param args: arguments to be passed to PeakList.to_str
    :param kwargs: keyword arguments to be passed to PeakList.to_str

    """
    if os.path.isfile(filename):
        logging.warning('plain text file [%s] already exists, override' % filename)
    with open(filename, 'w') as f: f.write(pkl.to_str(*args, **kwargs))


def load_peaklist_from_txt(filename, ID, delimiter=',', flag_names='auto', has_flag_col=True):
    """
    Loads a peaklist from plain text file.

    :param filename: path to an exiting text-based peaklist file
    :param ID: ID of the peaklist
    :param delimiter: delimiter of the text lines. Default = ',', i.e., CSV format
    :param flag_names: names of the flag attributes. Default = 'auto', indicating all the attribute names ends
        with "_flag" will be treated as flag attibute. Provide None to indicate no flag attributes
    :param has_flag_col: whether the text file contains the overall "flags" column. If True, it's values will be
        discarded. The overall flags of the new peaklist will be calculated automatically. Default = True
    :rtype: PeakList object

    """
    if not os.path.isfile(filename):
        raise IOError('plain text file [%s] does not exist' % filename)
    with open(filename, 'rU') as f:
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
                 [] if flag_names is None else set(flag_names)
    for n, v in zip(hd[2:], dm[2:]): pkl.add_attribute(n, _evalv(v), is_flag=n in flag_names, flagged_only=False)

    return pkl


# peak matrix portals
def save_peak_matrix_as_txt(pm, filename, *args, **kwargs):
    """
    Saves a peak matrix in plain text file.

    :param pm: the target peak matrix object
    :param filename: path to a new text file
    :param args: arguments to be passed to PeakMatrix.to_str
    :param kwargs: keyword arguments to be passed to PeakMatrix.to_str

    """
    if os.path.isfile(filename):
        logging.warning('plain text file [%s] already exists, override' % filename)
    with open(filename, 'w') as f:
        with unmask_all_peakmatrix(pm) as m: f.write(m.to_str(*args, **kwargs))


def load_peak_matrix_from_txt(filename, delimiter='\t', samples_in_rows=True, comprehensive='auto'):
    """
    Loads a peak matrix from plain text file.

    :param filename: path to an exiting text-based peak matrix file
    :param delimiter: delimiter of the text lines. Default = '\t', i.e., TSV format
    :param samples_in_rows: whether or not the samples are stored in rows. Default = True
    :param comprehensive: whether the input is a 'comprehensive' or 'simple' version of the matrix. Default = 'auto', i.e., auto detect
    :rtype: PeakMatrix object

    """
    if not os.path.isfile(filename):
        raise IOError('plain text file [%s] does not exist' % filename)
    with open(filename, 'rU') as f:
        rlns = filter(lambda x: x != '', f.readlines())

    dlns = map(lambda x: map(strip, x.split(delimiter)), rlns)
    if any(map(lambda x: len(x) != len(dlns[0]), dlns[1:])):
        raise IOError('data matrix size not match')

    if samples_in_rows: dlns = zip(*dlns)
    if comprehensive == 'auto': comprehensive = ('flags' in dlns[0])
    rdlns = zip(*dlns)
    rsdrow = filter(lambda x: x[1][0] == 'rsd_all', enumerate(rdlns))[0][0]

    def _parseflags():
        fgs = []
        for l, ln in enumerate(rdlns[rsdrow+1:]):
            if ln[0] == 'flags': break
            fgs += [(ln[0], map(eval, filter(lambda x: x != '', ln[1:])))]
        return fgs
    flgs = _parseflags() if comprehensive else []

    # must refactor if PeakMatrix.to_str changed
    pcol = rsdrow + len(flgs) + 2 if comprehensive else 1
    pids = dlns[0][pcol:]

    def _parsetags(tgs):
        for l, ln in enumerate(dlns[2:]):  # line 1 = missing
            if not ln[0].startswith('tags_'): break
            tn, tv = ln[0][5:], ln[pcol:]
            tl = filter(lambda x: x[1] != '', enumerate(_evalv(tv)))
            for i, v in tl: tgs[i].add_tag(v) if tn == 'untyped' else tgs[i].add_tag(v, tn)
        return l, tgs
    tnum, tags = 0, [PeakList_Tags() for _ in pids]
    if comprehensive: tnum, tags = _parsetags(tags)

    rlns = zip(*dlns[2 + tnum:])
    mz = np.array([rlns[0]] * len(pids), dtype=float)
    ints = np.array(rlns[pcol:], dtype=float)

    pm = PeakMatrix(pids, tags, [('mz', mz), ('intensity', ints)])
    for fn, fv in flgs: pm.add_flag(fn, fv, flagged_only = False)
    return pm

