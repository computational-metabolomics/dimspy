#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
hdf5Portal: PeakList and PeakMatrix HDF5 IO portals.

author(s): Albert Zhou, Ralf Weber
origin: Apr. 13, 2017

"""


import os, h5py, zipfile
import numpy as np
from dimspy.models.peaklist import PeakList, _Tags
from dimspy.models.peak_matrix import PeakMatrix, unmask_all_peakmatrix


# TODO: move these to a more general module?
# ----------------------------------------------------------
def check_paths(tsv, source):
    if tsv is None:
        if type(source) == str:
            source = source.encode('string-escape')
            if os.path.isdir(source):
                filenames = [os.path.join(source, fn) for fn in os.listdir(source) if fn.lower().endswith(".mzml") or fn.lower().endswith(".raw")]
            elif zipfile.is_zipfile(source):
                with zipfile.ZipFile(source) as zf:
                    assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
                    filenames = [fn for fn in zf.namelist() if fn.lower().endswith(".mzml")]
            elif h5py.is_hdf5(source):
                peaklists = load_peaklists_from_hdf5(source)
                filenames = [pl.ID for pl in peaklists]

        elif type(source) == list:
            assert isinstance(source[0], PeakList), "Incorrect Objects in list. PeakList class required."
            filenames = [pl.ID for pl in source]
        else:
            pass

    elif os.path.isfile(tsv.encode('string-escape')):
        tsv = tsv.encode('string-escape')
        fm = np.genfromtxt(tsv, dtype=None, delimiter="\t", names=True)
        if len(fm.shape) == 0:
            fm = np.array([fm])
        assert fm.dtype.names[0] == "filename" or fm.dtype.names[0] == "sample_id", \
            "Incorrect header for first column. Use filename or sample_id"

        filenames = []
        if type(source) == list:
            assert isinstance(source[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
            for fn in fm[fm.dtype.names[0]]:
                assert fn in [pl.ID for pl in source], "{} does not exist in list with Objects".format(fn)
                filenames.append(fn)
        elif type(source.encode('string-escape')) == str:
            source = source.encode('string-escape')
            if os.path.isdir(source):
                l = os.listdir(source)
                for fn in fm[fm.dtype.names[0]]:
                    assert os.path.basename(fn) in l, "{} does not exist in directory provided".format(os.path.basename(fn))
                    filenames.append(os.path.join(source, fn).replace('\\', r'\\'))

            elif zipfile.is_zipfile(source.encode('string-escape')):
                with zipfile.ZipFile(source.encode('string-escape')) as zf:
                    assert len([fn for fn in zf.namelist() if fn.lower().endswith(".raw")]) == 0, "Archive with *.raw files not yet supported. Convert to mzML"
                    for fn in fm[fm.dtype.names[0]]:
                        assert fn in zf.namelist(), "{} does not exist in .zip file".format(fn)
                        filenames.append(fn)

            elif h5py.is_hdf5(source):
                peaklists = load_peaklists_from_hdf5(source)
                filenames = [pl.ID for pl in peaklists]


            else:
                raise IOError("Can not read and parse {} or {}".format(source, tsv))
    else:
        raise IOError("File {} does not exist".format(tsv))

    return filenames


def load_peaklists(source):

    if type(source) == str:
        source = source.encode('string-escape')
        if h5py.is_hdf5(source):
            peaklists = load_peaklists_from_hdf5(source)
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

def hdf5_to_text(fname, path_out, separator="\t", transpose=False):
    assert os.path.isfile(fname), 'HDF5 database [%s] not exists' % fname
    assert h5py.is_hdf5(fname), 'input file [%s] is not a valid HDF5 database' % fname
    seps = {"comma": ",", "tab": "\t"}
    if separator in seps: separator = seps[separator]
    assert separator in [",", "\t"], "Incorrect separator ('tab', 'comma', ',', '\t')"
    f = h5py.File(fname, 'r')
    if "mz" in f:
        obj = load_peak_matrix_from_hdf5(fname)
        assert isinstance(obj, PeakMatrix)
        obj = load_peak_matrix_from_hdf5(fname)
        with open(os.path.join(path_out), "w") as pk_out:
            pk_out.write(obj.to_str(separator, transpose))
    else:
        assert os.path.isdir(path_out), "File or Directory does not exist:".format(path_out)
        obj = load_peaklists_from_hdf5(fname)
        assert isinstance(obj[0], PeakList), "Incorrect Objects in list. Peaklist Object required."
        for pl in obj:
            with open(os.path.join(path_out, os.path.splitext(pl.ID)[0] + ".txt"), "w") as pk_out:
                pk_out.write(pl.to_str(separator))
    return
# ----------------------------------------------------------



def save_peaklists_as_hdf5(pkls, fname):
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
        dset.attrs['tags'] = [(t if type(t) in (tuple, list) else ('None', t)) for t in pkl.tags.to_list()]
        for k, v in pkl.metadata.items(): dset.attrs['metadata_' + k] = v

    map(lambda x: _savepkl(*x), enumerate(pkls))


def load_peaklists_from_hdf5(fname):
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

        pkl.add_tags(*[t[1] for t in dset.attrs['tags'] if t[0] == 'None'], **dict(t for t in dset.attrs['tags'] if t[0] != 'None'))
        return dset.attrs['order'], pkl

    return zip(*sorted(map(_loadpkl, f.keys())))[1]


def save_peak_matrix_as_hdf5(pm, fname):
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
            dset.attrs['peaklist_tags_%d' % i] = [(t if type(t) in (tuple, list) else ('None', t)) for t in tags.to_list()]
    dset.attrs['mask'] = pm.mask


def load_peak_matrix_from_hdf5(fname):
    assert os.path.isfile(fname), 'HDF5 database [%s] not exists' % fname
    assert h5py.is_hdf5(fname), 'input file [%s] is not a valid HDF5 database' % fname
    f = h5py.File(fname, 'r')

    assert 'mz' in f, 'input database missing crucial attribute [mz]'
    dset = f['mz']
    assert dset.attrs.get('class', '') == 'PeakMatrix', 'input database is not a valid PeakMatrix'
    pids = dset.attrs['peaklist_ids']
    ptgs = [_Tags(*[t[1] for t in tags if t[0] == 'None'], **dict(t for t in tags if t[0] != 'None'))
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

    save_peaklists_as_hdf5(pkls, 'test_pl.hdf5')
    pkls = load_peaklists_from_hdf5('test_pl.hdf5')

    for p in pkls:
        if p.has_attribute('snr_flag'): p.drop_attribute('snr_flag')

    pm = PeakMatrix(
        [p.ID for p in pkls],
        [p.tags for p in pkls],
        {a: np.vstack([p.get_attribute(a) for p in pkls]) for a in pkls[0].attributes}
    )
    pm.unmask_tags('qc')

    #import pdb; pdb.set_trace()

    save_peak_matrix_as_hdf5(pm, 'test_pm.hdf5')
    pm = load_peak_matrix_from_hdf5('test_pm.hdf5')

    hdf5_to_text('test_pm.hdf5', "pm.txt", "\t", transpose=False)