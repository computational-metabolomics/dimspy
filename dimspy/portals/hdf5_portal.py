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
    try: return literal_eval(v)
    except (ValueError, SyntaxError): return str(v)

# peaklists portals
def save_peaklists_as_hdf5(pkls, fname):
    if os.path.isfile(fname):
        logging.warning('HDF5 database [%s] already exists, override' % fname)
    f = h5py.File(fname, 'w')

    def _savepkl(i, pkl):
        if pkl.ID in f.keys():
            raise IOError('peaklist [%s] already exists')
        dset = f.create_dataset(pkl.ID, pkl.full_shape[::-1])
        dset.attrs['class'] = 'PeakList'
        dset.attrs['order'] = i

        dset[...] = np.array(pkl.to_list()[:-1], dtype = float) # skip flags
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
            pkl.add_attribute(n, v, t, is_flag = (n in dset.attrs['flag_attrs']), flagged_only = False)

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
            raise IOError('attribute [%s] already exists')
        ds = f.create_dataset(attr, pm.full_shape)
        with unmask_all_peakmatrix(pm) as m: ds[...] = m.attr_matrix(attr)
    map(_saveattr, pm.attributes)

    dset = f['mz'] # must exists in pm
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
    ptgs = [PeakList_Tags(*[_eval(t[1]) for t in tags if t[0] == 'None'],
                          **{t[0]: _eval(t[1]) for t in tags if t[0] != 'None'})
            for n, tags in sorted(dset.attrs.items(), key = lambda x: x[0]) if n.startswith('peaklist_tags_')]
    adct = {attr: f[attr] for attr in f}
    mask = dset.attrs['mask']

    pm = PeakMatrix(pids, ptgs, **adct)
    pm.mask = mask
    return pm



# TODO: move these to a more general module?
# ----------------------------------------------------------
# def check_paths(tsv, source):
#     if tsv is None:
#         if type(source) == str:
#             source = source.encode('string-escape')
#             if os.path.isdir(source):
#                 filenames = [os.path.join(source, fn) for fn in os.listdir(source) if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]
#             elif zipfile.is_zipfile(source):
#                 with zipfile.ZipFile(source) as zf:
#                     assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
#                     filenames = [fn for fn in zf.namelist() if fn.lower().endswith(".mzml")]
#             elif h5py.is_hdf5(source):
#                 peaklists = load_peaklists_from_hdf5(source)
#                 filenames = [pl.ID for pl in peaklists]
#
#         elif type(source) == list:
#             assert isinstance(source[0], PeakList), "Incorrect Objects in list. PeakList class required."
#             filenames = [pl.ID for pl in source]
#         else:
#             pass
#
#     elif os.path.isfile(tsv.encode('string-escape')):
#         tsv = tsv.encode('string-escape')
#         fm = np.genfromtxt(tsv, dtype=None, delimiter="\t", names=True)
#         if len(fm.shape) == 0:
#             fm = np.array([fm])
#         assert fm.dtype.names[0] == "filename" or fm.dtype.names[0] == "sample_id", \
#             "Incorrect header for first column. Use filename or sample_id"
#
#         filenames = []
#         if type(source) == list:
#             assert isinstance(source[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
#             for fn in fm[fm.dtype.names[0]]:
#                 assert fn in [pl.ID for pl in source], "{} does not exist in list with Objects".format(fn)
#                 filenames.append(fn)
#         elif type(source.encode('string-escape')) == str:
#             source = source.encode('string-escape')
#             if os.path.isdir(source):
#                 l = os.listdir(source)
#                 for fn in fm[fm.dtype.names[0]]:
#                     assert os.path.basename(fn) in l, "{} does not exist in directory provided".format(os.path.basename(fn))
#                     filenames.append(os.path.join(source, fn).replace('\\', r'\\'))
#
#             elif zipfile.is_zipfile(source.encode('string-escape')):
#                 with zipfile.ZipFile(source.encode('string-escape')) as zf:
#                     assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
#                     for fn in fm[fm.dtype.names[0]]:
#                         assert fn in zf.namelist(), "{} does not exist in .zip file".format(fn)
#                         filenames.append(fn)
#
#             elif h5py.is_hdf5(source):
#                 peaklists = load_peaklists_from_hdf5(source)
#                 filenames = [pl.ID for pl in peaklists]
#
#
#             else:
#                 raise IOError("Can not read and parse {} or {}".format(source, tsv))
#     else:
#         raise IOError("File {} does not exist".format(tsv))
#
#     return filenames
#
#
# def load_peaklists(source):
#
#     if type(source) == str:
#         source = source.encode('string-escape')
#         if h5py.is_hdf5(source):
#             peaklists = load_peaklists_from_hdf5(source)
#         elif zipfile.is_zipfile(source):
#             zf = zipfile.ZipFile(source)
#             filenames = zf.namelist()
#             assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
#                 "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
#             peaklists = [text_to_peaklist(zf.open(fn), ID=os.path.basename(fn), has_flag_col=True) for fn in filenames]
#         elif os.path.isdir(source):
#             filenames = os.listdir(source)
#             assert len([fn for fn in filenames if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]) == 0,\
#                 "Incorrect format. Process .mzML and .raw files first using the \'process scans\' function"
#             peaklists = [text_to_peaklist(os.path.join(source, fn), ID=os.path.basename(fn), has_flag_col=False) for fn in filenames]
#         else:
#             raise TypeError("Incorrect format. Process .mzML and .raw files first using the 'process scans' function")
#     elif type(source) == list:
#         if not isinstance(source[0], PeakList):
#             raise TypeError("List has incorrect format. PeakList objects required.")
#         else:
#             peaklists = source
#     else:
#         raise IOError("Inccorrect input: list with peaklist objects or path")
#
#     return peaklists
#
# def hdf5_to_text(fname, path_out, separator="\t", transpose=False):
#     assert os.path.isfile(fname), 'HDF5 database [%s] not exists' % fname
#     assert h5py.is_hdf5(fname), 'input file [%s] is not a valid HDF5 database' % fname
#     seps = {"comma": ",", "tab": "\t"}
#     if separator in seps: separator = seps[separator]
#     assert separator in [",", "\t"], "Incorrect separator ('tab', 'comma', ',', '\t')"
#     f = h5py.File(fname, 'r')
#     if "mz" in f:
#         obj = load_peak_matrix_from_hdf5(fname)
#         assert isinstance(obj, PeakMatrix)
#         obj = load_peak_matrix_from_hdf5(fname)
#         with open(os.path.join(path_out), "w") as pk_out:
#             pk_out.write(obj.to_str(separator, transpose))
#     else:
#         assert os.path.isdir(path_out), "File or Directory does not exist:".format(path_out)
#         obj = load_peaklists_from_hdf5(fname)
#         assert isinstance(obj[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
#         for pl in obj:
#             with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
#                 pk_out.write(pl.to_str(separator))
#     return
# ----------------------------------------------------------

