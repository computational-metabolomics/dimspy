#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Workflow DIMS processing

author(s): Ralf Weber
origin: Nov. 2016
"""

from string import strip
from ast import literal_eval


def _loadPeaklist(file_name, ID, attr_names_dict=None, flag_names=None, delimiter='\t', has_flag_col=True):
    print 'loading [%s] ...' % file_name

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

    pkl = PeakList(ID, ddct['mz'][0], ddct['intensity'][0])
    for k, (v, f) in filter(lambda x: x[0] not in ('mz', 'intensity'), ddct.items()):
        pkl.add_attribute(k, v, is_flag=f, flagged_only=False)
    return pkl