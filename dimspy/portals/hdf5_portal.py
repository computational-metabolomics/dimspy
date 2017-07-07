#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
hdf5Portal: PeakList and PeakMatrix HDF5 IO portals.

author(s): Albert Zhou, Ralf Weber
origin: Apr. 13, 2017

"""


import os, logging, h5py
import numpy as np
from ast import literal_eval
from dimspy.models.peaklist_tags import PeakList_Tags
from dimspy.models.peaklist import PeakList
from dimspy.models.peak_matrix import PeakMatrix, unmask_all_peakmatrix


def _eval(v):
    try:
        return literal_eval(v)
    except (ValueError, SyntaxError):
        return str(v)


# peaklists portals
def save_peaklists_as_hdf5(pkls, fname):
    if os.path.isfile(fname):
        logging.warning('HDF5 database [%s] already exists and override' % fname)
    f = h5py.File(fname, 'w')

    def _savepkl(i, pkl):
        if pkl.ID in f.keys():
            raise IOError('peaklist [%s] already exists' % pkl.ID)
        dset = f.create_dataset(pkl.ID, pkl.full_shape[::-1], dtype=np.float64)
        dset.attrs['class'] = 'PeakList'
        dset.attrs['order'] = i

        dset[...] = np.array(pkl.to_list()[:-1], dtype=float) # skip flags
        dset.attrs['dtable_names'], dset.attrs['dtable_types'] = zip(*pkl.dtable.dtype.descr)

        dset.attrs['flag_attrs'] = pkl.flag_attributes
        dset.attrs['tags'] = [(t if type(t) in (tuple, list) else ('None', t)) for t in pkl.tags.to_list()]
        for k, v in pkl.metadata.items(): dset.attrs['metadata_' + k] = v

    map(lambda x: _savepkl(*x), enumerate(pkls))

def load_peaklists_from_hdf5(fname):
    if not os.path.isfile(fname):
        raise IOError('HDF5 database [%s] not exists' % fname)
    if not h5py.is_hdf5(fname):
        raise IOError('input file [%s] is not a valid HDF5 database' % fname)
    f = h5py.File(fname, 'r')

    def _loadpkl(ID):
        dset = f[ID]
        if dset.attrs.get('class', '') != 'PeakList':
            raise IOError('unknown object found in the database')

        dm = dset.value
        dn = dset.attrs['dtable_names']
        dt = dset.attrs['dtable_types']

        if dn[0] != 'mz' or dn[1] != 'intensity':
            raise IOError('PANIC: HDF5 dataset matrix not in order')
        pkl = PeakList(ID, dm[0], dm[1], **{k[9:]: v for k,v in dset.attrs.items() if k.startswith('metadata_')})

        for n, v, t in zip(dn[2:], dm[2:], dt[2:]):
            pkl.add_attribute(n, v, t, is_flag=(n in dset.attrs['flag_attrs']), flagged_only=False)

        pkl.tags.add_tags(*[_eval(t[1]) for t in dset.attrs['tags'] if t[0] == 'None'],
                          **{t[0]: _eval(t[1]) for t in dset.attrs['tags'] if t[0] != 'None'})
        return dset.attrs['order'], pkl

    return zip(*sorted(map(_loadpkl, f.keys())))[1]


# peak matrix portals
def save_peak_matrix_as_hdf5(pm, fname):
    if os.path.isfile(fname):
        logging.warning('HDF5 database [%s] already exists and override' % fname)
    f = h5py.File(fname, 'w')

    def _saveattr(attr):
        if attr in f.keys():
            raise IOError('attribute [%s] already exists' % attr)
        ds = f.create_dataset(attr, pm.full_shape, dtype=np.float64)
        with unmask_all_peakmatrix(pm) as m:
            ds[...] = m.attr_matrix(attr)
    map(_saveattr, pm.attributes)

    dset = f['mz']  # must exists in pm
    dset.attrs['class'] = 'PeakMatrix'
    with unmask_all_peakmatrix(pm):
        dset.attrs['peaklist_ids'] = pm.peaklist_ids
        for i, tags in enumerate(pm.peaklist_tags):
            dset.attrs['peaklist_tags_%d' % i] = [(t if type(t) in (tuple, list) else ('None', t)) for t in tags.to_list()]
    dset.attrs['mask'] = pm.mask


def load_peak_matrix_from_hdf5(fname):
    if not os.path.isfile(fname):
        raise IOError('HDF5 database [%s] not exists' % fname)
    if not h5py.is_hdf5(fname):
        raise IOError('input file [%s] is not a valid HDF5 database' % fname)
    f = h5py.File(fname, 'r')

    if 'mz' not in f:
        raise IOError('input database missing crucial attribute [mz]')
    dset = f['mz']
    if dset.attrs.get('class', '') != 'PeakMatrix':
        raise IOError('input database is not a valid PeakMatrix')
    pids = dset.attrs['peaklist_ids']
    tatt = sorted(filter(lambda x: x.startswith('peaklist_tags_'), dset.attrs.keys()), key=lambda x: int(x[14:]))
    ptgs = [PeakList_Tags(*[_eval(t[1]) for t in tags if t[0] == 'None'],
                          **{t[0]: _eval(t[1]) for t in tags if t[0] != 'None'})
            for tags in map(lambda x: dset.attrs[x], tatt)]
    adct = {attr: f[attr] for attr in f}
    mask = dset.attrs['mask']

    pm = PeakMatrix(pids, ptgs, **adct)
    pm.mask = mask
    return pm

