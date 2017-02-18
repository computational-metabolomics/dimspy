#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Cluster and Align peaklists

author(s): Ralf Weber, Albert Zhou
origin: Oct. 2016

references:
http://stats.stackexchange.com/questions/3685/where-to-cut-a-dendrogram
https://people.duke.edu/~ccc14/sta-663/Optimization_Bakeoff.html
http://stackoverflow.com/questions/28687882/cutting-scipy-hierarchical-dendrogram-into-clusters-via-a-threshold-value

"""


from __future__ import division

import logging, sys
import numpy as np
import fastcluster as fc
from string import join
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from scipy import cluster
from scipy.spatial.distance import squareform
from ..models.peak_matrix import PeakMatrix


def _counter(lst):
    dct = defaultdict(lambda: 0)
    for k in set(lst): dct[k] += 1
    return dct
if sys.version_info < (2, 7):
    counter = _counter
else:
    from collections import Counter
    counter = Counter


# single cluster
def _cluster_peaks(mzs, ppm, distype='euclidean', linkmode='centroid'):
    if len(mzs) == 1:
        return np.zeros_like(mzs, dtype=int).reshape((-1, 1))

    m = np.column_stack([mzs])
    mdist = fc.pdist(m, metric=distype)

    outer_mzs = np.add.outer(mzs, mzs)
    np.fill_diagonal(outer_mzs, 0)
    avg_mz_pair = np.divide(outer_mzs, 2)
    mdist_mz_pair = squareform(avg_mz_pair)
    relative_errors = np.multiply(mdist_mz_pair, 1e-6)

    with np.errstate(divide='ignore', invalid='ignore'): # using errstate context to avoid seterr side effects
        m_mass_tol = np.divide(mdist, relative_errors)
        m_mass_tol[np.isnan(m_mass_tol)] = 0.0
    z = fc.linkage(m_mass_tol, method=linkmode)

    # cut tree at ppm threshold & order matches the order of mzs
    return cluster.hierarchy.cut_tree(z, height=ppm)


# multiprocess cluster
def _cluster_peaks_mp(params): return _cluster_peaks(*params)

def _cluster_peaks_map(mzs, ppm, block_size, byunique=False, edge_extend=10, ncpus=None):
    assert np.all(mzs[1:] >= mzs[:-1]), 'mz values not in ascending order'
    block_size = min(len(mzs), block_size)
    assert len(mzs) >= block_size > 0, 'block size not in (0, #peaks) range'

    # create mapping pool
    def _mmap(f, p):
        ppool = Pool(cpu_count()-1 if ncpus is None else ncpus)
        rets = ppool.map(f, p)
        ppool.close()  # close after parallel finished
        return rets

    def _smap(f, p):
        return map(f, p)

    pmap = _smap if ncpus == 1 or cpu_count() <= 2 else _mmap

    # split blocks
    if byunique:
        umzs, urids = np.unique(mzs, return_index=True)
        sids = [urids[i] for i in range(block_size, len(umzs), block_size)] # appx size
    else:
        sids = range(block_size, len(mzs), block_size)
    if len(sids) > 0 and sids[-1] == len(mzs)-1: sids = sids[:-1] # no need to cluster the final peak separately


    # align edges
    _rng = lambda x: (lambda v: np.where(np.logical_and(x-v < mzs, mzs <= x+v))[0])(edge_extend * ppm * x * 1e-6)
    erngs = [_rng(mzs[i]) for i in sids]
    _cids = pmap(_cluster_peaks_mp, [(mzs[r], ppm) for r in erngs])
    eblks = [r[c == c[r == s]] for s, r, c in zip(sids, erngs, map(lambda x: x.flatten(), _cids))]
    ecids = map(lambda x: np.zeros_like(x).reshape((-1, 1)), eblks)
    
    # align blocks
    brngs = np.array((0,) + reduce(lambda x,y: x+y, map(lambda x: (x[0], x[-1]+1), eblks), ()) + (len(mzs),)).reshape((-1, 2)) # keep () in reduce in case eblks is empty
    if len(set(brngs[0])) == 1:
        brngs = brngs[1:] # in case edges have reached mz bounds
    if len(set(brngs[-1])) == 1:
        brngs = brngs[:-1]
    bkmzs = [mzs[slice(*r)] for r in brngs]
    assert all(map(lambda x: abs(x[-1]-x[0]) / x[0] > edge_extend * ppm * 10 * 1e-6, bkmzs[:-1])), \
        'mz value range in cluster block is too small, consider increase the number of peaks in each block'
    bcids = pmap(_cluster_peaks_mp, [(m, ppm) for m in bkmzs])

    # combine
    cids = [None] * (len(bcids) + len(ecids))
    cids[::2], cids[1::2] = bcids, ecids
    return cids


def _cluster_peaks_reduce(clusters):
    return reduce(lambda x, y: np.vstack((x, y+np.max(x)+1)), clusters).flatten()


# alignment
def _align_peaks(cids, pids, *attrs):
    assert all(map(lambda x: x.shape == cids.shape == pids.shape, attrs)), 'attributes shape not match'

    # encode string id list to continuous values and search uniques
    def _idsmap(ids):
        ri, vi = np.unique(ids, return_index=True, return_inverse=True)[1:]
        sri = np.argsort(ri)  # ensure order
        return np.argsort(sri)[vi], ids[ri[sri]]
    (mcids, mpids), (ucids, upids) = zip(*map(_idsmap, (cids, pids)))

    # count how many peaks from same sample being clustered into one peak
    cM = np.zeros(map(len, (upids, ucids)))
    for pos, count in counter(zip(mpids, mcids)).items():
        cM[pos] = count

    # fill all the attributes into matrix
    def _fillam(a):
        aM = np.zeros(map(len, (upids, ucids)))
        for p, v in zip(zip(mpids, mcids), a):
            aM[p] += v
        with np.errstate(divide='ignore', invalid='ignore'):
            aM /= cM
        aM[np.isnan(aM)] = 0
        return aM
    attrMs = map(_fillam, attrs)

    # sort mz values, ensure mzs matrix be the first
    sortids = np.argsort(np.average(attrMs[0], axis=0, weights=attrMs[0].astype(bool)))
    return upids, map(lambda x: x[:, sortids], attrMs + [cM])


# interface
def align_peaks(peaks, ppm=2.0, block_size=2000, byunique=False, ncpus=None):
    # remove empty peaklists
    emlst = np.array(map(lambda x: x.size == 0, peaks))
    if np.sum(emlst) > 0:
        logging.warning('droping empty peaklist(s) [%s]' % join(map(str, [p.ID for e, p in zip(emlst, peaks) if e]), ','))
    peaks = [p for e, p in zip(emlst, peaks) if not e]

    # obtain attrs
    attrs = peaks[0].attributes
    assert attrs[:2] == ('mz', 'intensity'), 'PANIC: peak attributes in wrong order'
    assert all(map(lambda x: attrs == x.attributes, peaks)), 'peak attributes not the same'

    # flatten
    f_pids = np.hstack(map(lambda p: [p.ID] * p.size, peaks))
    f_attrs = map(lambda attr: np.hstack(map(lambda p: p[attr], peaks)), attrs)

    sortids = np.argsort(f_attrs[0])  # attrs[0] -> mz values
    s_pids = f_pids[sortids]
    s_attrs = map(lambda x: x[sortids], f_attrs)

    # cluster
    clusters = _cluster_peaks_map(s_attrs[0], ppm=ppm, block_size=block_size, byunique=byunique, ncpus=ncpus)
    cids = _cluster_peaks_reduce(clusters)

    # align
    a_pids, a_attrms = _align_peaks(cids, s_pids, *s_attrs)
    assert 'alignment_counts' not in attrs, 'preserved attribute name [alignment_counts] already exists'
    attrs += ('alignment_counts', )  # for cM

    # sort by original pid
    pids = f_pids[sorted(np.unique(f_pids, return_index=True)[1])]
    pdct = dict((i, mi) for mi, i in enumerate(a_pids))
    porder = [pdct[i] for i in pids]
    o_attrms = map(lambda x: x[porder], a_attrms)

    return PeakMatrix(pids, [p.tags for p in peaks], dict(zip(attrs, o_attrms)))

