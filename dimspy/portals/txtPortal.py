#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
txtPortal: PeakList and PeakMatrix plain text IO portals.

author(s): Albert Zhou, Ralf Weber
origin: Sep. 27, 2016

"""


import os
import numpy as np
from string import strip
from ast import literal_eval
from dimspy.models.peaklist import PeakList, _Tags
from dimspy.models.peak_matrix import PeakMatrix


def _evalv(vect):
    try: ctype = type(literal_eval(vect[0]))
    except (ValueError, SyntaxError): ctype = None
    return vect if ctype is None else map(ctype, vect)

def save_peaklist_as_txt(pkl, file_name, *args, **kwargs):
    if os.path.isfile(file_name): print '\nWARNING: plain text file [%s] already exists and overrided\n' % file_name
    with open(file_name, 'w') as f: f.write(pkl.to_str(*args, **kwargs))

def load_peaklist_from_txt(file_name, ID, delimiter=',', flag_names='auto', has_flag_col=True):
    assert os.path.isfile(file_name), 'plain text file [%s] not exists' % file_name
    with open(file_name, 'rU') as f: rlns = filter(lambda x: x != '', map(strip, f.readlines()))

    dlns = map(lambda x: map(strip, x.split(delimiter)), rlns)
    assert all(map(lambda x: len(x) == len(dlns[0]), dlns[1:])), 'data matrix size not match'

    hd, dm = dlns[0], zip(*dlns[1:])
    if has_flag_col: hd, dm = hd[:-1], dm[:-1] # flag_col must be the last one, and discarded

    mzs, ints = np.array(dm[0], dtype = float), np.array(dm[1], dtype = float) # first two cols must be mz and ints
    pkl = PeakList(ID, mzs, ints)

    assert len(set(hd)) == len(hd), 'duplicate headers found'
    flag_names = filter(lambda x: x.endswith('_flag'), hd) if flag_names == 'auto' else \
                 [] if flag_names is None else set(flag_names)

    for n, v in zip(hd[2:], dm[2:]): pkl.add_attribute(n, _evalv(v), is_flag = n in flag_names, flagged_only = False)
    return pkl

def save_peak_matrix_as_txt(pm, file_name, *args, **kwargs):
    if os.path.isfile(file_name): print '\nWARNING: plain text file [%s] already exists and overrided\n' % file_name
    with open(file_name, 'w') as f: f.write(pm.to_str(*args, **kwargs))

def load_peak_matrix_from_txt(file_name, delimiter='\t', transposed=False, extended=False):
    assert os.path.isfile(file_name), 'plain text file [%s] not exists' % file_name
    with open(file_name, 'rU') as f: rlns = filter(lambda x: x != '', f.readlines())

    dlns = map(lambda x: map(strip, x.split(delimiter)), rlns)
    assert all(map(lambda x: len(x) == len(dlns[0]), dlns[1:])), 'data matrix size not match'

    if not transposed: dlns = zip(*dlns)
    pids = dlns[0][5:] if extended else dlns[0][1:]

    def _parsetags(tgs):
        for l, ln in enumerate(dlns[2:]): # line 1 = missing
            if not ln[0].startswith('tags_'): break
            tn, tv = ln[0][5:], ln[5:]
            tl = filter(lambda x: x[1] != '', enumerate(_evalv(tv)))
            for i, v in tl: tgs[i].add_tags(v) if tn == 'untyped' else tgs[i].add_tags(**{tn: v})
        return l, tgs
    tnum, tags = 0, [_Tags() for _ in pids]
    if extended: tnum, tags = _parsetags(tags)

    rlns = zip(*dlns[2+tnum:])
    mz, ints = np.array([rlns[0]] * len(pids), dtype = float), np.array(rlns[5:] if extended else rlns[1:], dtype = float)
    return PeakMatrix(pids, tags, {'mz': mz, 'intensity': ints})


# testing
if __name__ == '__main__':
    _mzs = lambda: sorted(np.random.uniform(0, 1000, size=1000))
    _ints = lambda: np.abs(np.random.normal(10, 3, size=1000))

    pkls = [
        PeakList('sample_1_1', _mzs(), _ints(), mz_range=(0, 1000)),
        PeakList('sample_1_2', _mzs(), _ints(), mz_range=(0, 1000)),
        PeakList('QC_1', _mzs(), _ints(), mz_range=(0, 1000)),
        PeakList('sample_2_1', _mzs(), _ints(), mz_range=(0, 1000)),
        PeakList('sample_2_2', _mzs(), _ints(), mz_range=(0, 1000)),
        PeakList('QC_2', _mzs(), _ints(), mz_range=(0, 1000)),
    ]

    pkls[0].add_tags('sample', treatment='compound_1', time_point='1hr', plate=1)
    pkls[1].add_tags('sample', treatment='compound_1', time_point='6hr', plate=1)
    pkls[2].add_tags('qc', plate=1)
    pkls[3].add_tags('sample', treatment='compound_2', time_point='1hr', plate=2)
    pkls[4].add_tags('sample', treatment='compound_2', time_point='6hr', plate=2)
    pkls[5].add_tags('qc', plate=2)

    pkls[0].metadata['file_name'] = 'S1_1.txt'
    pkls[1].metadata['file_name'] = 'S1_2.txt'
    pkls[2].metadata['file_name'] = 'QC_1.txt'
    pkls[3].metadata['file_name'] = 'S2_1.txt'
    pkls[4].metadata['file_name'] = 'S2_2.txt'
    pkls[5].metadata['file_name'] = 'QC_2.txt'

    pkls[0].add_attribute('snr_flag', [0, 1] * 500, is_flag = True, flagged_only = False)

    save_peaklist_as_txt(pkls[0], 'test_pl.txt')
    import pdb; pdb.set_trace()
    pkls[0] = load_peaklist_from_txt('test_pl.txt', pkls[0].ID)

    pkls[0].add_tags('sample', treatment = 'compound_1', time_point = '1hr', plate = 1) # load from txt will lost all tags
    pkls[0].drop_attribute('snr_flag')

    from dimspy.process.peak_alignment import align_peaks
    pm = align_peaks(pkls, 5000.)

    save_peak_matrix_as_txt(pm, 'test_pm.txt', transpose=True, extend=True)
    import pdb; pdb.set_trace()
    pm = load_peak_matrix_from_txt('test_pm.txt',transposed=True, extended=True)



