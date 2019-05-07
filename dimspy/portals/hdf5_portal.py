#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
The PeakList and PeakMatrix HDF5 portals.

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 1.0.0

"""


import os, logging, zlib, h5py
import pickle as cp
import numpy as np
from typing import Sequence
from ast import literal_eval
from dimspy.models.peaklist_tags import Tag, PeakList_Tags
from dimspy.models.peaklist import PeakList
from dimspy.models.peak_matrix import PeakMatrix, unmask_all_peakmatrix


def _eval(v):
    try:
        return literal_eval(v)
    except (ValueError, SyntaxError):
        return str(v)

_packMeta = lambda x: np.string_(zlib.compress(cp.dumps(x)) + b'\xFF') # numpy truncates right-side \x00
_unpackMeta = lambda x: cp.loads(zlib.decompress(x[:-1]))


# peaklists portals
def save_peaklists_as_hdf5(pkls: Sequence[PeakList], filename: str):
    """
    Saves multiple peaklists in a HDF5 file.

    :param pkls: the target list of peaklist objects
    :param filename: path to a new HDF5 file

    To incorporate with different dtypes in the attribute matrix, this portal converts all the arribute values
    into fix-length strings for HDF5 data tables storage. The order of the peaklists will be retained.

    """
    if os.path.isfile(filename):
        logging.warning('HDF5 database [%s] already exists and override' % filename)
    f = h5py.File(filename, 'w')

    def _savepkl(i, pkl):
        if pkl.ID in f.keys():
            raise IOError('peaklist [%s] already exists' % pkl.ID)

        dm = np.string_(pkl.to_list()[:-1])
        dset = f.create_dataset(pkl.ID, pkl.full_shape[::-1], dtype=dm.dtype)
        dset.attrs['class'] = 'PeakList'
        dset.attrs['order'] = i

        dset[...] = np.array(dm) # skip flags
        dset.attrs['dtable_names'], dset.attrs['dtable_types'] = list(map(np.string_, zip(*pkl.dtable.dtype.descr)))

        dset.attrs['flag_attrs'] = np.string_(pkl.flag_attributes)
        dset.attrs['tags'] = np.string_([(t or 'None', v) for v,t in pkl.tags.to_list()])
        for k, v in pkl.metadata.items(): dset.attrs['metadata_' + k] = _packMeta(v)

    for pl in enumerate(pkls): _savepkl(*pl)


def load_peaklists_from_hdf5(filename: str):
    """
    Loads a list of peaklist objects from a HDF5 file.

    :param filename: path to a HDF5 file
    :rtype: list

    The values in HDF5 data tables are automatically converted to their original dtypes before loading in the peaklist.

    """
    if not os.path.isfile(filename):
        raise IOError('HDF5 database [%s] does not exist' % filename)
    if not h5py.is_hdf5(filename):
        raise IOError('input file [%s] is not a valid HDF5 database' % filename)
    f = h5py.File(filename, 'r')

    def _loadpkl(ID):
        dset = f[ID]
        if dset.attrs.get('class', '') != 'PeakList':
            raise IOError('unknown object found in the database')

        dm = dset.value
        dn = dset.attrs['dtable_names'].astype(str)
        dt = dset.attrs['dtable_types'].astype(str)

        if dn[0] != 'mz' or dn[1] != 'intensity':
            raise IOError('PANIC: HDF5 dataset matrix not in order')
        pkl = PeakList(ID, dm[0].astype(np.float64), dm[1].astype(np.float64),
                       **{k[9:]: _unpackMeta(v) for k,v in list(dset.attrs.items()) if k.startswith('metadata_')})

        for n, v, t in zip(dn[2:], dm[2:], dt[2:]):
            pkl.add_attribute(n, v, t, is_flag=(n in dset.attrs['flag_attrs'].astype(str)), flagged_only=False)

        for t,v in map(lambda x: x.astype(str), dset.attrs['tags']): pkl.tags.add_tag(_eval(v), None if t == 'None' else t)
        return dset.attrs['order'], pkl

    return list(zip(*sorted(map(_loadpkl, f.keys()))))[1]


# peak matrix portals
def save_peak_matrix_as_hdf5(pm: PeakMatrix, filename: str):
    """
    Saves a peak matrix object to a HDF5 file.

    :param pm: the target peak matrix object
    :param filename: path to a new HDF5 file

    The order of the attributes and flags will be retained.

    """
    if os.path.isfile(filename):
        logging.warning('HDF5 database [%s] already exists and override' % filename)
    f = h5py.File(filename, 'w')

    def _saveattr(attr):
        if attr in f.keys():
            raise IOError('attribute [%s] already exists' % attr)

        with unmask_all_peakmatrix(pm) as m:
            dm = m.attr_matrix(attr, flagged_only = False)

        dt = np.float64 if dm.dtype.kind == 'f' else \
             np.int64 if dm.dtype.kind in ('i', 'u') else \
             ('S%d' % np.max([list(map(len,l)) for l in dm]))

        ds = f.create_dataset(attr, dm.shape, dtype=dt)
        ds[...] = dm.astype(dt)
        ds.attrs['dtype'] = dm.dtype.str

    for a in pm.attributes: _saveattr(a)

    dset = f['mz']  # must exists in pm
    dset.attrs['class'] = 'PeakMatrix'
    dset.attrs['attributes'] = np.string_(pm.attributes)
    dset.attrs['mask'] = pm.mask

    with unmask_all_peakmatrix(pm):
        dset.attrs['peaklist_ids'] = np.string_(pm.peaklist_ids)
        for i, tags in enumerate(pm.peaklist_tags):
            dset.attrs['peaklist_tags_%d' % i] = np.string_([(t or 'None', v) for v,t in tags.to_list()])

        dset.attrs['flag_names'] = np.string_(pm.flag_names)
        for fn in pm.flag_names:
            dset.attrs[fn] = pm.flag_values(fn)


def load_peak_matrix_from_hdf5(filename):
    """
    Loads a peak matrix from a HDF5 file.

    :param filename: path to an existing HDF5 file
    :rtype: PeakMatrix object

    """
    if not os.path.isfile(filename):
        raise IOError('HDF5 database [%s] does not exist' % filename)
    if not h5py.is_hdf5(filename):
        raise IOError('input file [%s] is not a valid HDF5 database' % filename)
    f = h5py.File(filename, 'r')

    if 'mz' not in f:
        raise IOError('input database missing crucial attribute [mz]')

    dset = f['mz']
    if dset.attrs.get('class', '') != 'PeakMatrix':
        raise IOError('input database is not a valid PeakMatrix')
    attl = dset.attrs['attributes'].astype(str)
    pids = dset.attrs['peaklist_ids'].astype(str)
    mask = dset.attrs['mask']

    tatt = sorted([x for x in dset.attrs.keys() if x.startswith('peaklist_tags_')], key=lambda x: int(x[14:]))
    ptgs = [PeakList_Tags(*[Tag(_eval(v), None if t == 'None' else t) for t,v in map(lambda x: x.astype(str), tags)]) for tags in [dset.attrs[x] for x in tatt]]

    flgs = [(fn, dset.attrs[fn]) for fn in dset.attrs['flag_names'].astype(str)]
    alst = [(attr, np.array(f[attr]).astype(f[attr].attrs['dtype'])) for attr in attl]

    pm = PeakMatrix(pids, ptgs, alst)
    pm.mask = mask
    for fn, fv in flgs: pm.add_flag(fn, fv, flagged_only = False)
    return pm

