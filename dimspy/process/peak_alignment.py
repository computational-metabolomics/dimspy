#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Cluster and align peaklists into peak matrix.

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 1.0.0

"""

from __future__ import division

import logging
import numpy as np
import fastcluster as fc
from string import join
from operator import itemgetter
from multiprocessing import Pool, cpu_count
from collections import Counter
from scipy import cluster
from scipy.spatial.distance import squareform
from dimspy.models.peak_matrix import PeakMatrix


# single cluster
def _cluster_peaks(mzs, ppm, distype='euclidean', linkmode='centroid'):
    if len(mzs) == 0:
        return np.array([])
    if len(mzs) == 1:
        return np.zeros_like(mzs, dtype=int).reshape((-1, 1))

    m = np.column_stack([mzs])
    mdist = fc.pdist(m, metric=distype)

    outer_mzs = np.add.outer(mzs, mzs)
    np.fill_diagonal(outer_mzs, 0)
    avg_mz_pair = np.divide(outer_mzs, 2)
    mdist_mz_pair = squareform(avg_mz_pair)
    relative_errors = np.multiply(mdist_mz_pair, 1e-6)

    with np.errstate(divide='ignore', invalid='ignore'):  # using errstate context to avoid seterr side effects
        m_mass_tol = np.divide(mdist, relative_errors)
        m_mass_tol[np.isnan(m_mass_tol)] = 0.0
    z = fc.linkage(m_mass_tol, method=linkmode)

    # cut tree at ppm threshold & order matches the order of mzs
    return cluster.hierarchy.cut_tree(z, height=ppm)


# multiprocess cluster
def _cluster_peaks_mp(params):
    return _cluster_peaks(*params)


def _cluster_peaks_map(mzs, ppm, block_size, fixed_block, edge_extend, ncpus):
    if not np.all(mzs[1:] >= mzs[:-1]):
        raise ValueError('mz values not in ascending order')
    if not 1 <= block_size <= len(mzs):
        # logging.warning('block size (%d) not in range [1, #peaks (%d)]' % (block_size, len(mzs)))
        block_size = min(max(block_size, 1), len(mzs))

    # split blocks
    if fixed_block:
        sids = range(block_size, len(mzs), block_size)
    else:
        umzs, urids = np.unique(mzs, return_index=True)
        sids = [urids[i] for i in range(block_size, len(umzs), block_size)]  # appx size
    if len(sids) > 0 and sids[-1] == len(mzs) - 1:
        sids = sids[:-1]  # no need to cluster the final peak separately

    # create mapping pool
    def _mmap(f, p):
        ppool = Pool(cpu_count() - 1 if ncpus is None else ncpus)
        rets = ppool.map(f, p)
        ppool.close()  # close after parallel finished
        return rets

    def _smap(f, p):
        return map(f, p)

    def _pmap(f, p):
        largechk = filter(lambda x: len(x[0]) > 1E+5, p)
        if len(largechk) > 0:
            raise RuntimeError('Some of the clustering chunks contain too many peaks: \n%s' %
                join(map(lambda x: 'mz range [%.5f - %.5f] ... [%d] peaks' % (min(x[0]),max(x[0]),len(x[0])), largechk), '\n'))
        return (_smap if ncpus == 1 or cpu_count() <= 2 else _mmap)(f, p)

    # align edges
    eeppm = edge_extend * ppm * 1e-6
    _rng = lambda x: (lambda v: np.where(np.logical_and(x - v < mzs, mzs <= x + v))[0])(eeppm * x)
    erngs = [_rng(mzs[i]) for i in sids]

    overlap = [p[-1] >= s[0] for p, s in zip(erngs[:-1], erngs[1:])]  # in case edges have overlap
    if True in overlap:
        logging.warning('[%d] edge blocks overlapped, consider increasing the block size' % (sum(overlap) + 1))
        erngs = reduce(lambda x, y: (x[:-1] + [np.unique(np.hstack((x[-1], y[0])))]) if y[1] else x + [y[0]],
                       zip(erngs[1:], overlap), [erngs[0]])
        sids = [sids[0]] + [s for s, o in zip(sids[1:], overlap) if not o]

    _cids = _pmap(_cluster_peaks_mp, [(mzs[r], ppm) for r in erngs])
    eblks = [r[c == c[r == s]] for s, r, c in zip(sids, erngs, map(lambda x: x.flatten(), _cids))]
    ecids = map(lambda x: np.zeros_like(x).reshape((-1, 1)), eblks)

    # align blocks
    # keep () in reduce in case eblks is empty
    brngs = np.array(
        (0,) + reduce(lambda x, y: x + y, map(lambda x: (x[0], x[-1] + 1), eblks), ()) + (len(mzs),)
    ).reshape((-1, 2))

    # in case edges have reached mz bounds
    bkmzs = [mzs[slice(*r)] for r in brngs]
    slimbk = map(lambda x: len(x) == 0 or abs(x[-1] - x[0]) / x[0] < eeppm * 10, bkmzs)
    if np.sum(slimbk) > 0:
        pbrngs = [map(lambda x: min(x, len(mzs) - 1), (r[0], r[-1] - 1 if r[-1] != r[0] else r[-1])) for r in brngs]
        pblns = ['block %d' % i + ': [%f, %f]' % itemgetter(*r)(mzs)
                 for i, (s, r) in enumerate(zip(slimbk, pbrngs)) if s]
        logging.warning('[%d] empty / slim clustering block(s) found, consider increasing the block size\n%s' %
                        (np.sum(slimbk), join(pblns, '\n')))
    bcids = _pmap(_cluster_peaks_mp, [(m, ppm) for m in bkmzs])

    # combine
    cids = [None] * (len(bcids) + len(ecids))
    cids[::2], cids[1::2] = bcids, ecids
    return cids


def _cluster_peaks_reduce(clusters):
    return reduce(lambda x, y: np.vstack((x, y + np.max(x) + 1)), filter(lambda x: len(x) > 0, clusters)).flatten()


# alignment
def _align_peaks(cids, pids, *attrs):
    if not all(map(lambda x: x.shape == cids.shape == pids.shape, attrs)):
        raise ValueError('attributes shape not match')

    # encode string id list to continuous values and search uniques
    def _idsmap(ids):
        _, ri, vi = np.unique(ids, return_index=True, return_inverse=True)
        sri = np.argsort(ri)  # ensure order
        return np.argsort(sri)[vi], ids[ri[sri]]

    (mcids, mpids), (ucids, upids) = zip(*map(_idsmap, (cids, pids)))

    # count how many peaks from same sample being clustered into one peak
    cM = np.zeros(map(len, (upids, ucids)))
    for pos, count in Counter(zip(mpids, mcids)).items():
        cM[pos] = count

    # fill all the attributes into matrix
    def _avg_am(a):
        aM = np.zeros(map(len, (upids, ucids)))
        for p, v in zip(zip(mpids, mcids), a):
            aM[p] += v
        with np.errstate(divide='ignore', invalid='ignore'):
            aM /= cM
        aM[np.isnan(aM)] = 0
        return aM

    def _cat_am(a):
        aM = [[[] for _ in ucids] for _ in upids]
        for (r, c), v in zip(zip(mpids, mcids), a):
            aM[r][c] += [str(v)]
        aM = [[join(val, ',') for val in ln] for ln in aM]
        return np.array(aM)

    def _fillam(a):
        alg = _avg_am if a.dtype.kind in ('i', 'u', 'f') else \
              _cat_am if a.dtype.kind in ('?', 'b', 'a', 'S', 'U') else \
              lambda x: logging.warning('undefined alignment behaviour for [%s] dtype data') # returns None
        return alg(a)

    attrMs = map(_fillam, attrs)

    # sort mz values, ensure mzs matrix be the first
    sortids = np.argsort(np.average(attrMs[0], axis=0, weights=attrMs[0].astype(bool)))
    return upids, map(lambda x: x[:, sortids], attrMs + [cM])


# interface
def align_peaks(peaks, ppm=2.0, block_size=5000, fixed_block=True, edge_extend=10, ncpus=None):
    """
    Cluster and align peaklists into a peak matrix.

    :param peaks: list of peaklists for alignment
    :param ppm: the hierarchical clustering cutting height, i.e., ppm range for each aligned mz value. Default = 2.0
    :param block_size: number peaks in each centre clustering block. This can be a exact or approximate number depends
        on the fixed_block parameter. Default = 5000
    :param fixed_block: whether the blocks contain fixed number of peaks. Default = True
    :param edge_extend: ppm range for the edge blocks. Default = 10
    :param ncpus: number of CPUs for parallel clustering. Default = None, indicating using as many as possible
    :rtype: PeakMatrix object

    .. figure::  images/alignment.png
        :align:   center

    This function uses hierarchical clustering to align the mz values of the input peaklists. The alignment "width" is
    decided by the parameter of ppm. Due to a large number of peaks, this function splits them into blocks with fixed
    or approximate length, and clusters in a parallel manner on multiple CPUs. When running, the edge blocks are
    clustered first to prevent separating the same peak into two adjacent centre blocks. The size of the edge blocks is
    decided by edge_extend. The clustering of centre blocks is conducted afterwards.

    After merging the clustering results, all the attributes (mz, intensity, snr, etc.) are aligned into matrix
    accordingly. If multiple peaks from the same sample are clustered into one mz value, their attributes are averaged
    (for real value attributes e.g. mz and intensity) or concatenated (string, unicode, or bool attributes). The flag
    attributes are ignored. The number of these overlapping peaks is recorded in a new intra_count attribute matrix.

    """
    # remove empty peaklists
    emlst = np.array(map(lambda x: x.size == 0, peaks))
    if np.sum(emlst) > 0:
        logging.warning(
            'droping empty peaklist(s) [%s]' % join(map(str, [p.ID for e, p in zip(emlst, peaks) if e]), ','))
        peaks = [p for e, p in zip(emlst, peaks) if not e]
    if len(peaks) == 0: raise ValueError('all input peaklists for alignment are empty')

    # obtain attrs
    attrs = peaks[0].attributes
    if attrs[:2] != ('mz', 'intensity'):
        raise AttributeError('PANIC: peak attributes in wrong order')
    if not all(map(lambda x: attrs == x.attributes, peaks)):
        raise ValueError('peak attributes not the same')
    if 'intra_count' in attrs:
        raise AttributeError('preserved attribute name [intra_count] already exists')
    attrs = filter(lambda x: x not in peaks[0].flag_attributes, attrs) # flags should be excluded

    # single peaklist
    if len(peaks) == 1:
        attrlst = [(a, peaks[0][a].reshape((1, -1))) for a in attrs] + \
                  [('intra_count', np.ones((1, peaks[0].size)))]
        return PeakMatrix([peaks[0].ID], [peaks[0].tags], attrlst)

    # flatten
    f_pids = np.hstack(map(lambda p: [p.ID] * p.size, peaks))
    f_attrs = map(lambda attr: np.hstack(map(lambda p: p[attr], peaks)), attrs)

    sortids = np.argsort(f_attrs[0])  # attrs[0] -> mz values
    s_pids = f_pids[sortids]
    s_attrs = map(lambda x: x[sortids], f_attrs)

    # cluster
    clusters = _cluster_peaks_map(s_attrs[0], ppm, block_size, fixed_block, edge_extend, ncpus)
    cids = _cluster_peaks_reduce(clusters)

    # align
    a_pids, a_attrms = _align_peaks(cids, s_pids, *s_attrs)
    attrs += ('intra_count',)  # for cM

    # sort by original pid
    pids = f_pids[sorted(np.unique(f_pids, return_index=True)[1])]
    pdct = dict((i, mi) for mi, i in enumerate(a_pids))
    porder = [pdct[i] for i in pids]
    o_attrms = map(lambda x: x[porder] if x is not None else None, a_attrms)

    return PeakMatrix(pids, [p.tags for p in peaks], filter(lambda x: x[1] is not None, zip(attrs, o_attrms)))
