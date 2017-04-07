#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Workflow DIMS processing

author(s): Ralf Weber
origin: Nov. 2016
"""


import os, h5py
import numpy as np
from string import strip
from ast import literal_eval
from models import PeakList, PeakMatrix, unmask_all_peakmatrix

# this is the only exception, should always use peaklists to create peakmatrix in other parts of the source codes
from models.peaklist import _Tags


def loadPeaklistFromTxt(file_name, ID, attr_names_dict=None, flag_names=None, delimiter='\t', has_flag_col=True):
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


def savePeaklistsToHDF(pkls, fname):
    if os.path.isfile(fname): print '\nWARNING: HDF5 database [%s] already exists and overrided\n' % fname
    f = h5py.File(fname, 'w')

    def _savepkl(i, pkl):
        assert pkl.ID not in f.keys(), 'peaklist [%s] already exists'
        dset = f.create_dataset(pkl.ID, pkl.full_shape[::-1])
        dset.attrs['class'] = 'PeakList'
        dset.attrs['order'] = i

        dset[...] = np.array(pkl.to_list()[:-1], dtype = float) # skip flags
        dset.attrs['dtable_names'], dset.attrs['dtable_types'] = zip(*pkl.dtable.dtype.descr)

        dset.attrs['flag_attrs'] = pkl.flag_attributes
        dset.attrs['tags'] = [(t if type(t) in (tuple, list) else ('NA', t)) for t in pkl.tags.to_list()]
        for k, v in pkl.metadata.items(): dset.attrs['metadata_' + k] = v

    map(lambda x: _savepkl(*x), enumerate(pkls))

def loadPeaklistsFromHDF(fname):
    assert os.path.isfile(fname), 'HDF5 database [%s] not exists' % fname
    assert h5py.is_hdf5(fname), 'input file [%s] is not a valid HDF5 database' % fname
    f = h5py.File(fname, 'r')

    def _loadpkl(ID):
        dset = f[ID]
        assert dset.attrs.get('class', '') == 'PeakList', 'unknown object found in the database'

        dm = dset.value
        dn = dset.attrs['dtable_names']
        dt = dset.attrs['dtable_types']

        assert dn[0] == 'mz' and dn[1] == 'intensity', 'PANIC: HDF5 dataset matrix not in order'
        pkl = PeakList(ID, dm[0], dm[1], **{k[9:]: v for k,v in dset.attrs.items() if k.startswith('metadata_')})

        for n, v, t in zip(dn[2:], dm[2:], dt[2:]):
            pkl.add_attribute(n, v, t, is_flag = (n in dset.attrs['flag_attrs']), flagged_only = False)

        pkl.add_tags(*[t[1] for t in dset.attrs['tags'] if t[0] == 'NA'], **dict(t for t in dset.attrs['tags'] if t[0] != 'NA'))
        return dset.attrs['order'], pkl

    return zip(*sorted(map(_loadpkl, f.keys())))[1]


def savePeakMatrixToHDF(pm, fname):
    if os.path.isfile(fname): print '\nWARNING: HDF5 database [%s] already exists and overrided\n' % fname
    f = h5py.File(fname, 'w')

    def _saveattr(attr):
        assert attr not in f.keys(), 'attribute [%s] already exists'
        dset = f.create_dataset(attr, pm.full_shape)
        dset[...] = pm.attr_matrix(attr, masked_only = False)
    map(_saveattr, pm.attributes)

    dset = f['mz'] # must exists in pm
    dset.attrs['class'] = 'PeakMatrix'
    with unmask_all_peakmatrix(pm):
        dset.attrs['peaklist_ids'] = pm.peaklist_ids
        for i, tags in enumerate(pm.peaklist_tags):
            dset.attrs['peaklist_tags_%d' % i] = [(t if type(t) in (tuple, list) else ('NA', t)) for t in tags.to_list()]
    dset.attrs['mask'] = pm.mask

def loadPeakMatrixFromHDF(fname):
    assert os.path.isfile(fname), 'HDF5 database [%s] not exists' % fname
    assert h5py.is_hdf5(fname), 'input file [%s] is not a valid HDF5 database' % fname
    f = h5py.File(fname, 'r')

    assert 'mz' in f, 'input database missing crucial attribute [mz]'
    dset = f['mz']
    assert dset.attrs.get('class', '') == 'PeakMatrix', 'input database is not a valid PeakMatrix'
    pids = dset.attrs['peaklist_ids']
    ptgs = [_Tags(*[t[1] for t in tags if t[0] == 'NA'], **dict(t for t in tags if t[0] != 'NA'))
            for n, tags in sorted(dset.attrs.items(), key = lambda x: x[0]) if n.startswith('peaklist_tags_')]
    adct = {attr: f[attr] for attr in f}
    mask = dset.attrs['mask']

    return PeakMatrix(pids, ptgs, adct, mask)


# testing
if __name__ == '__main__':
    _mzs = lambda: sorted(np.random.uniform(100, 1000, size=100))
    _ints = lambda: np.abs(np.random.normal(10, 3, size=100))

    pkls = [
        PeakList('sample_1_1', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('sample_1_2', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('QC_1', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('sample_2_1', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('sample_2_2', _mzs(), _ints(), mz_range=(100, 1000)),
        PeakList('QC_2', _mzs(), _ints(), mz_range=(100, 1000)),
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

    pkls[0].add_attribute('snr_flag', [0, 1] * 50, is_flag = True, flagged_only = False)
    pkls[3].add_attribute('snr_flag', [1, 0] * 50, is_flag = True, flagged_only = False)


    savePeaklistsToHDF(pkls, 'test_pl.hdf5')
    pkls = loadPeaklistsFromHDF('test_pl.hdf5')


    for p in pkls:
        if p.has_attribute('snr_flag'): p.drop_attribute('snr_flag')

    pm = PeakMatrix(
        [p.ID for p in pkls],
        [p.tags for p in pkls],
        {a: np.vstack([p.get_attribute(a) for p in pkls]) for a in pkls[0].attributes}
    )
    pm.unmask_tags('qc')

    import pdb; pdb.set_trace()

    savePeakMatrixToHDF(pm, 'test_pm.hdf5')
    pm = loadPeakMatrixFromHDF('test_pm.hdf5')

