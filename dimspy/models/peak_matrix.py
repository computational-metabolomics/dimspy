#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
The PeakMatrix data object class.

.. moduleauthor:: Albert Zhou, Ralf Weber

.. versionadded:: 0.1

"""

from __future__ import division

import logging
import numpy as np
from collections import OrderedDict, Iterable
from string import join
from peaklist import PeakList


class PeakMatrix(object):
    """
    The PeakMatrix class.

    Stores aligned mass spectrometry peaks matrix data. It requires IDs, tags, and attributes from the source peak
    lists. It uses tags based mask to "hide" the unrelated samples for convenient processing. It utilises the
    automatically managed flags to "remove" peaks without actually delete them. Therefore the filterings on the peaks
    are traceable. Normally, PeakMatrix object is created by functions e.g. align_peaks() rather than manual.

    :param peaklist_ids: the IDs of the source peak lists
    :param peaklist_tags: the tags (PeakList_Tags) of the source peak lists
    :param peaklist_attributes: the attributes of the source peak lists. Must be a list or tuple in the format of
        [(attr_name, attr_matrix), ...], where attr_name is name of the attribute, and attr_matrix is the vertically
        stacked arrtibute values in the shape of samples x peaks. The order of the attributes will be kept in the
        PeakMatrix. The first two attributes must be "mz" and "intensity".

    >>> pids = [pl.ID for pl in peaklists]
    >>> tags = [pl.tags for pl in peaklists]
    >>> attrs = [(attr_name, np.vstack([pl[attr_name] for pl in peaklists])) \
                 for attr_name in peaklists[0].attributes]
    >>> pm = PeakMatrix(pids, tags, attrs)

    Internally the attribute data is stored in OrderedDict as a list of matrix. An attribute matrix can be illustrated
    as follows, in which the mask and flags are the same for all attributes. The final row "flags" is automatically
    calculated based on the manually added flags. It decides which peaks are "removed" i.e. unflagged. Particularly,
    the "--" indicates no peak in that sample can be aligned into the mz value.

    **attribute: "mz"**

    +------------+--------+--------+--------+-----+
    | mask       | peak_1 | peak_2 | peak_3 | ... |
    +============+========+========+========+=====+
    | False      | 12.7   | 14.9   | 21.0   |     |
    +------------+--------+--------+--------+     +
    | True       | --     | 15.1   | 21.1   | ... |
    +------------+--------+--------+--------+     +
    | False      | 12.1   | 14.7   | --     |     |
    +------------+--------+--------+--------+     +
    | False      | 12.9   | 14.8   | 20.9   |     |
    +------------+--------+--------+--------+-----+
    | ...        |        |        |        |     |
    +------------+--------+--------+--------+-----+
    | **flag_1** | True   | False  | True   |     |
    +------------+--------+--------+--------+     +
    | **flag_2** | True   | True   | False  | ... |
    +------------+--------+--------+--------+     +
    | **flags*** | True   | False  | False  |     |
    +------------+--------+--------+--------+-----+

    .. warning::
        Removing a flag may change the overall "flags", and cause the unflagged peaks to be flagged again. As
        most the processes are applied only on the flagged peaks, these peaks, if the others have gone through such process,
        may have incorrect values.

        In principle, setting a flag attribute should be considered as an irreversible process.

    Different from the flags, mask should be considered as a more temporary way to hide the unrelated samples. A masked
    sample (row) will not be used for processing, but its data is still in the attribute matrix. For this reason,
    the mask_peakmatrix, unmask_peakmatrix, and unmask_all_peakmatrix statements are provided as a more flexible way
    to set / unset the mask.

    """

    def __init__(self, peaklist_ids, peaklist_tags, peaklist_attributes):
        if len(set(peaklist_ids)) != len(peaklist_ids):
            raise ValueError('duplicate peaklist IDs found')
        self._pids = np.array(peaklist_ids)
        self._tags = np.array(peaklist_tags)

        self._attr_dict = OrderedDict((str(k), np.array(v)) for k,v in peaklist_attributes)
        if not (self._attr_dict.keys()[0] == 'mz' and self._attr_dict.keys()[1] == 'intensity'):
            raise ValueError('required attribute fields "mz" and "intensity" not available')

        pnum, psize = self._attr_dict.values()[0].shape
        if not len(peaklist_ids) == len(peaklist_tags) == pnum:
            raise ValueError('input peaklists number not match')
        if not all(map(lambda x: x.shape == (pnum, psize), self._attr_dict.values())):
            raise ValueError('input attribute matrix shape not match')

        self._mask = np.zeros(pnum, dtype=bool)
        self._flags_dict = OrderedDict()

    # privates
    # np.average accepts weights, set to 0 to exclude value
    def _present_mask(self, x, _func, flagged_only):
        pmx = self.property('present_matrix', flagged_only)
        mskx = np.ma.masked_array(x, mask = ~pmx)
        vals = _func(mskx)
        rval = np.array(vals)
        rval[vals.mask] = np.nan
        return rval

    def _present_mean(self, x, axis, flagged_only):
        return self._present_mask(x, lambda v: np.mean(v, axis = axis), flagged_only)

    def _present_std(self, x, axis, flagged_only):
        return self._present_mask(x, lambda v: np.std(v, axis = axis), flagged_only)

    # built-ins
    def __len__(self):
        return self.shape[0]

    def __str__(self):
        return self.to_str(attr_name='intensity', delimiter=', ', transpose=True, comprehensive=False)

    # properties
    @property
    def mask(self):  # mask is different from flags in that they will not be masked by themselves, i.e., more temporary
        """
        Property of the mask.

        :getter: returns a deep copy of the mask array
        :setter: sets the mask array. Provide None to unmask all samples
        :type: numpy array

        """
        return self._mask.copy()

    @mask.setter
    def mask(self, value):
        self._mask = np.zeros(self.full_shape[0], dtype = bool)
        if value is not None: self._mask[value] = True

    @property
    def flag_names(self):
        """
        Property of the flag names.

        :getter: returns a tuple including the names of the manually set flags
        :type: tuple

        """
        return tuple(self._flags_dict.keys())

    @property
    def flags(self):
        """
        Property of the flags.

        :getter: returns a deep copy of the flags array
        :type: numpy array

        """
        return np.ones(self.full_shape[1], dtype = bool) if len(self._flags_dict) == 0 else \
               np.logical_and.reduce(self._flags_dict.values())

    @property
    def attributes(self):
        """
        Property of the attribute names.

        :getter: returns a tuple including the names of the attribute matrix
        :type: tuple

        """
        return tuple(self._attr_dict.keys())

    @property
    def peaklist_ids(self):
        """
        Property of the source peaklist IDs.

        :getter: returns a tuple including the IDs of the source peaklists
        :type: tuple

        """
        return tuple(self._pids[~self._mask])

    @property
    def peaklist_tags(self):
        """
        Property of the source peaklist tags.

        :getter: returns a tuple including the Peaklist_Tags objects of the source peaklists
        :type: tuple

        """
        return tuple(self._tags[~self._mask])

    @property
    def peaklist_tag_types(self):
        """
        Property of the source peaklist tag types.

        :getter: returns a tuple including the types of the typed tags of the source peaklists
        :type: tuple

        """
        return tuple(set(reduce(lambda x, y: x + y, [t.tag_types for t in self.peaklist_tags], ())))

    @property
    def peaklist_tag_values(self):
        """
        Property of the source peaklist tag values.

        :getter: returns a tuple including the values of the source peaklists tags, both typed and untyped
        :type: tuple

        """
        return tuple(set(reduce(lambda x, y: x + y, [t.tag_values for t in self.peaklist_tags], ())))

    @property
    def shape(self):
        """
        Property of the peak matrix shape.

        :getter: returns the shape of the attribute matrix
        :type: tuple

        """
        return np.sum(~self._mask), np.sum(self.flags)

    @property
    def full_shape(self):
        """
        Property of the peak matrix full shape.

        :getter: returns the full shape of the attribute matrix, i.e., ignore mask and flags
        :type: tuple

        """
        return self._attr_dict['mz'].shape

    @property
    def present(self):
        """
        Property of the present array.

        :getter: returns the present array, indicating how many peaks are aligned in each mz value
        :type: numpy array

        """
        return  self.property('present')

    @property
    def present_matrix(self):
        """
        Property of the present matrix.

        :getter: returns the present matrix, indicating whether a sample has peak(s) aligned in each mz value
        :type: numpy array

        >>> print pm.present_matrix
        array([[ True,  True,  True,  True, False],
               [ True,  True, False, False,  True],
               [ True,  True,  True,  True,  True],
               [False,  True, False,  True,  True],])
        >>> print pm.present
        array([3, 4, 2, 3, 3])

        """
        return self.property('present_matrix')

    @property
    def fraction(self):
        """
        Property of the fraction array.

        :getter: returns the fraction array, indicating the ratio of present peaks on each mz value
        :type: numpy array

        >>> print pm.present
        array([3, 4, 2, 3, 3])
        >>> print pm.shape[0]
        4
        >>> print pm.fraction
        array([0.75, 1.0, 0.5, 0.75, 0.75])

        """
        return self.property('fraction')

    @property
    def missing_values(self):
        """
        Property of the missing values array.

        :getter: returns the missing values array, indicating the number of unaligned peaks on each sample
        :type: numpy array

        >>> print pm.present_matrix
        array([[ True,  True,  True,  True, False],
               [ True,  True, False, False,  True],
               [ True,  True,  True,  True,  True],
               [False,  True, False,  True,  True],])
        >>> print pm.missing_values
        array([1, 2, 0, 2])

        """
        return self.property('missing_values')

    @property
    def occurrence(self):
        """
        Property of the occurrence array.

        :getter: returns the occurrence array, indicating the total number of peaks (including peaks in the same sample)
            aliged in each mz value. This property is valid only when the *intra_count* attribute matrix is available
        :type: numpy array

        >>> print pm.attr_matrix('intra_count')
        array([[ 2,  1,  1,  1,  0],
               [ 1,  1,  0,  0,  1],
               [ 1,  3,  1,  2,  1],
               [ 0,  1,  0,  1,  1],])
        >>> print pm.occurrence
        array([ 4,  6,  2,  4,  3])

        """
        return self.property('occurrence')

    @property
    def purity(self):
        """
        Property of the purity level array.

        :getter: returns the purity array, indicating the ratio of only one peak in each sample being aligned in each mz
            value. This property is valid only when the *intra_count* attribute matrix is available
        :type: numpy array

        >>> print pm.attr_matrix('intra_count')
        array([[ 2,  1,  1,  1,  0],
               [ 1,  1,  0,  0,  1],
               [ 1,  3,  1,  2,  1],
               [ 0,  1,  0,  1,  1],])
        >>> print pm.purity
        array([ 0.667,  0.75,  1.0,  0.667,  1.0])

        """
        return self.property('purity')

    # attribute matrix
    @property
    def mz_matrix(self):
        """
        Property of the mz matrix.

        :getter: returns the mz attribute matrix, unmasked and flagged values only
        :type: numpy array
        
        """       
        return self.attr_matrix('mz')

    @property
    def intensity_matrix(self):
        """
        Property of the intensity matrix.

        :getter: returns the intensity attribute matrix, unmasked and flagged values only
        :type: numpy array
        
        """        
        return self.attr_matrix('intensity')

    @property
    def mz_mean_vector(self):
        """
        Property of the mz mean values array.

        :getter: returns the mean values array of the mz attribute matrix, unmasked and flagged values only
        :type: numpy array
        
        """               
        return self.attr_mean_vector('mz')

    @property
    def intensity_mean_vector(self):
        """
        Property of the intensity mean values array.

        :getter: returns the mean values array of the intensity attribute matrix, unmasked and flagged values only
        :type: numpy array
        
        """                   
        return self.attr_mean_vector('intensity')

    # publics
    # tags and mask
    def tags_of(self, tag_type=None):
        """
        Obtains tags of the peaklist_tags with particular tag type.

        :param tag_type: the type of the returning tags. Provide None to obtain untyped tags
        :rtype: tuple

        """        
        if not (tag_type is None or all(map(lambda x: x.has_tag_type(tag_type), self.peaklist_tags))):
            raise ValueError('not all samples has tag type [%s]' % tag_type)
        tlst = [t.tag_of(tag_type) for t in self.peaklist_tags]
        if tag_type is None: tlst = reduce(lambda x, y: x + y, tlst)
        return tuple(set(tlst))

    def mask_tags(self, *args, **kwargs):  # match to all
        """
        Masks samples with particular tags.

        :param args: target tag values, both typed and untyped
        :param kwargs: target typed tag types and values
        :param override: whether to override the current mask, default = False
        :rtype: PeakMatrix object (self)

        This function will mask samples with ALL the tags. To match ANY of the tags, use cascade form instead.
        
        >>> pm.mask_tags('qc', plate = 1)
        (will mask all QC samples on plate 1)
        >>> pm.mask_tags('qc').mask_tags(plate = 1)
        (will mask QC samples and all samples on plate 1)

        """            
        override = kwargs.pop('override') if kwargs.has_key('override') else False
        mask = map(lambda x: all(map(lambda t: x.has_tag(t), args)) and
                             all(map(lambda t: x.has_tag(**dict([t])), kwargs.items())), self._tags)
        self.mask = np.logical_or(False if override else self._mask, mask)
        return self

    def unmask_tags(self, *args, **kwargs):  # match to all
        """
        Unmasks samples with particular tags.

        :param args: target tag values, both typed and untyped
        :param kwargs: target typed tag types and values
        :param override: whether to override the current mask, default = False
        :rtype: PeakMatrix object (self)

        This function will unmask samples with ALL the tags. To unmask ANY of the tags, use cascade form instead.
        
        >>> pm.mask = [True] * pm.full_shape[0]
        >>> pm.unmask_tags('qc', plate = 1)
        (will unmask all QC samples on plate 1)
        >>> pm.unmask_tags('qc').unmask_tags(plate = 1)
        (will unmask QC samples and all samples on plate 1)

        """                
        override = kwargs.pop('override') if kwargs.has_key('override') else False
        mask = map(lambda x: not (all(map(lambda t: x.has_tag(t), args)) and
                                  all(map(lambda t: x.has_tag(**dict([t])), kwargs.items()))), self._tags)
        self.mask = np.logical_and(True if override else self._mask, mask)
        return self

    # flags
    def add_flag(self, flag_name, flag_values, flagged_only = True):
        """
        Adds a flag to the peak matrix peaks.

        :param flag_name: name of the flag, it must be unique and not equal to "flags"
        :param flag_values: values of the flag. It must have a length of pm.shape[1] if flagged_only = True, or
            pm.full_shape[1] if flagged_only = False
        :param flagged_only: whether to set the flagged peaks only. Default = True, and the values of the unflagged peaks
            are set to False

        The overall flags property will be automatically recalculated.

        """                
        if flag_name == 'flags':
            raise KeyError('reserved flag name [flags] cannot be added')
        if self._flags_dict.has_key(flag_name):
            raise KeyError('flag name [%s] already exists' % flag_name)
        if flagged_only:
            if not len(flag_values) == self.shape[1]: raise ValueError('flag values and peak matrix shape not match')
            newf = np.zeros(self.full_shape[1], dtype = bool)
            newf[self.flags] = flag_values
            self._flags_dict[flag_name] = newf
        else:
            if not len(flag_values) == self.full_shape[1]: raise ValueError('flag values and peak matrix shape not match')
            self._flags_dict[flag_name] = np.array(flag_values, dtype = bool)

    def drop_flag(self, flag_name):
        """
        Drops a existing flag from the peak matrix.

        :param flag_name: name of the flag to drop. It must exist and not equal to "flags"

        The overall flags property will be automatically recalculated.

        """                
        if flag_name == 'flags':
            raise KeyError('reserved flag name [flags] cannot be droppped')
        if not self._flags_dict.has_key(flag_name):
            raise KeyError('flag name [%s] not exists' % flag_name)
        del self._flags_dict[flag_name]

    def flag_values(self, flag_name):
        """
        Obtains values of an existing flag.

        :param flag_name: name of the target flag. It must exist and not equal to "flags"
        :rtype: numpy array

        """                
        if flag_name == 'flags':
            raise KeyError('use PeakMatrix.flags to access reserved flag name [flags]')
        if not self._flags_dict.has_key(flag_name):
            raise KeyError('flag name [%s] not exists' % flag_name)
        return self._flags_dict[flag_name]

    # access
    def property(self, prop_name, flagged_only = True):
        """
        Obtains an existing attribute matrix.

        :param prop_name: name of the target property. Valid properties include 'present', 'present_matrix', 'fraction',
            'missing_values', 'occurrence', and 'purity'
        :param flagged_only: whether to return the flagged values only. Default = True
        :rtype: numpy array

        """
        if prop_name in ('occurrence', 'purity') and 'intra_count' not in self.attributes:
            logging.warning("attribute matrix ['intra_count'] not available")

        # use lambda to postpone the actual calculation
        _prop = {
            'present_matrix':
                lambda: self.attr_matrix('intra_count' if 'intra_count' in self.attributes else 'mz', flagged_only) > 0,
            'present':
                lambda: np.sum(self.property('present_matrix', flagged_only), axis = 0),
            'fraction':
                lambda: self.property('present', flagged_only) / self.shape[0],
            'missing_values':
                lambda: np.sum(~self.property('present_matrix', flagged_only), axis = 1),
            'occurrence':
                lambda: np.ones(self.shape[1] if flagged_only else self.full_shape[1]) if 'intra_count' not in self.attributes else
                        np.sum(self.attr_matrix('intra_count', flagged_only), axis = 0),
            'purity':
                lambda: np.ones(self.shape[1] if flagged_only else self.full_shape[1]) if 'intra_count' not in self.attributes else
                        1 - np.sum(self.attr_matrix('intra_count', flagged_only) > 1, axis = 0) /
                            np.sum(self.property('present_matrix', flagged_only), axis = 0),
        }.get(prop_name, None)

        if _prop is None: raise ValueError('unknown property name [%s]' % prop_name)
        return _prop()

    def attr_matrix(self, attr_name, flagged_only = True):
        """
        Obtains an existing attribute matrix.

        :param attr_name: name of the target attribute
        :param flagged_only: whether to return the flagged values only. Default = True
        :rtype: numpy array

        """                        
        if not self._attr_dict.has_key(attr_name):
            raise KeyError('attribute matrix [%s] not available' % attr_name)
        aM = self._attr_dict[attr_name][~self._mask]
        return aM[:,self.flags] if flagged_only else aM

    def attr_mean_vector(self, attr_name, flagged_only = True):
        """
        Obtains the mean array of an existing attribute matrix.

        :param attr_name: name of the target attribute
        :param flagged_only: whether to return the mean array of the flagged values only. Default = True
        :rtype: numpy array

        Noting that only the "present" peaks will be used for mean values calculation. If the attribute matrix has a
        string / unicode data type, the values in each column will be concatenated.

        """                                
        aM = self.attr_matrix(attr_name, flagged_only)
        aV = self._present_mean(aM, 0, flagged_only) if aM.dtype.kind in ('i', 'u', 'f') else \
             np.array([join([str(v) for v,p in zip(ln,self.present_matrix[:,j]) if p],';') for j,ln in enumerate(zip(*aM))]) # for strings
        return aV

    def remove_samples(self, sample_ids, masked_only=True):
        """
        Removes samples from the peak matrix.

        :param sample_ids: the indices of the samples to remove
        :param masked_only: whether the indices are for unmasked samples or all samples. Default = True
        :rtype: PeakMatrix object (self)

        """        
        if isinstance(sample_ids, Iterable): sample_ids = list(sample_ids)
        rmids = np.where(~self.mask)[0][sample_ids] if masked_only else sample_ids
        self._pids = np.delete(self._pids, rmids, axis=0)
        self._tags = np.delete(self._tags, rmids, axis=0)
        self._attr_dict = OrderedDict((k, np.delete(v, rmids, axis=0)) for k,v in self._attr_dict.items())
        self._mask = np.delete(self._mask, rmids, axis=0)
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def remove_peaks(self, peak_ids, flagged_only = True):
        """
        Removes peaks from the peak matrix.

        :param peak_ids: the indices of the peaks to remove
        :param flagged_only: whether the indices are for flagged peaks or all peaks. Default = True
        :rtype: PeakMatrix object (self)

        """
        if isinstance(peak_ids, Iterable): peak_ids = list(peak_ids)
        rmids = np.where(self.flags)[0][peak_ids] if flagged_only else peak_ids
        self._attr_dict = OrderedDict((k, np.delete(v, rmids, axis=1)) for k,v in self._attr_dict.items())
        self._flags_dict = OrderedDict((k, np.delete(v, rmids)) for k,v in self._flags_dict.items())
        if self.is_empty(): logging.warning('matrix is empty after removal')
        return self

    def is_empty(self):
        """
        Checks whether the peak matrix is empty under the current mask and flags.

        :rtype: bool

        """
        return 0 in self.shape

    # exports
    def extract_peaklist(self, peaklist_id):
        """
        Extracts one peaklist from the peak matrix.

        :param peaklist_id: ID of the peaklist to extract
        :rtype: PeakList object

        Only the "present" peaks will be included in the result peaklist.

        """
        if peaklist_id not in self.peaklist_ids:
            raise ValueError('peaklist id has to match those in the peak matrix')
        idx = self.peaklist_ids.index(peaklist_id)

        nzero_idx = self.present_matrix[idx]
        pl = PeakList(peaklist_id, self.mz_matrix[idx,nzero_idx], self.intensity_matrix[idx,nzero_idx])
        for attr in self.attributes[2:]: pl.add_attribute(attr, self.attr_matrix(attr)[idx,nzero_idx])

        return pl

    def extract_peaklists(self):
        """
        Extracts all peaklists from the peak matrix.

        :rtype: list

        """
        return map(self.extract_peaklist, self.peaklist_ids)

    def to_peaklist(self, ID):
        """
        Averages the peak matrix into a single peaklist.

        :param ID: ID of the merged peaklist
        :rtype: PeakList object

        Only the "present" peaks will be included in the result peaklist. The new peaklist will only contain the
        following attributes: mz, intensity, present, fraction, rsd, occurence, and purity.

        """
        presids = self.present > 0 # presented peaks only
        if False in presids:
            logging.warning('[%d] empty peaks removed when exporting PeakMatrix to PeakList' % np.sum(~presids))
        pl = PeakList(ID, self.mz_mean_vector[presids], self.intensity_mean_vector[presids], aligned_ids = self.peaklist_ids)
        pl.add_attribute('present', self.present[presids])
        pl.add_attribute('fraction', self.fraction[presids])
        pl.add_attribute('occurrence', self.occurrence[presids])
        pl.add_attribute('purity', self.purity[presids])
        return pl

    def to_str(self, attr_name='intensity', delimiter='\t', transpose=False, comprehensive=True):
        """
        Exports the peak matrix to a string.

        :param attr_name: name of the attribute matrix for exporting. Default = 'intensity'
        :param delimiter: delimiter to separate the matrix. Default = '\t', i.e., TSV format
        :param transpose: whether to transpose the matrix. Default = False
        :param comprehensive: whether to include comprehensive info, e.g., mask, flags, present, rsd etc. Default = True
        :rtype: str

        """
        hd = ['m/z'] + map(str, self.attr_mean_vector('mz', flagged_only = False))
        dm = [map(str, self.peaklist_ids)] + \
             [map(str, ln) for ln in self.attr_matrix(attr_name, flagged_only = False).T]

        if comprehensive:
            ttypes = set(reduce(lambda x, y: x + y, map(lambda x: x.tag_types, self.peaklist_tags)))
            tnum = len(ttypes)
            hd = [hd[0]] + ['missing values'] + map(lambda x: 'tags_' + x, ttypes) + ['tags_untyped'] + hd[1:]
            dm = [dm[0]] + \
                 [map(str, self.missing_values)] + \
                 [map(lambda x: str(x.tag_of(t)) if x.has_tag_type(t) else '', self.peaklist_tags) for t in ttypes] + \
                 [map(lambda x: join(map(str, x.tag_of(None)), ';'), self.peaklist_tags)] + \
                 dm[1:]

            prelst = ['present']    + ([''] * (tnum + 2)) + map(str, self.property('present', flagged_only = False))
            ocrlst = ['occurrence'] + ([''] * (tnum + 2)) + map(str, self.property('occurrence', flagged_only = False))
            puplst = ['purity']     + ([''] * (tnum + 2)) + map(str, self.property('purity', flagged_only = False))
            flgmtx = [[fn] + ([''] * (tnum + 2)) + map(str, self.flag_values(fn).astype(int)) for fn in self.flag_names]
            flglst = ['flags']      + ([''] * (tnum + 2)) + map(str, self.flags.astype(int))
            dm = zip(*([prelst, ocrlst, puplst] + flgmtx + [flglst] + zip(*dm)))

        lm = [hd] + zip(*dm)
        return join(map(lambda x: join(x, delimiter), zip(*lm) if transpose else lm), '\n')


# with statements
class mask_peakmatrix:
    """
    The mask_peakmatrix statement.

    Temporary mask the peak matrix with particular tags. Within the statement the samples can be motified or removed.
    After leaving the statement the original mask will be recoverd.

    :param pm: the target peak matrix
    :param override: whether to override the current mask, default = True
    :param args: target tag values, both typed and untyped
    :param kwargs: target typed tag types and values
    :rtype: PeakMatrix object

    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')
    >>> with mask_peakmatrix(pm., 'qc') as m: print m.peaklist_ids
    ('sample_1', 'sample_2', 'sample_3', 'sample_4')
    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')

    """

    def __init__(self, pm, *args, **kwargs):
        self._pm = pm
        self._utags = args
        self._ttags = kwargs
        if not self._ttags.has_key('override'): self._ttags['override'] = True  # default for with statement
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.mask_tags(*self._utags, **self._ttags)
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]


class unmask_peakmatrix:
    """
    The unmask_peakmatrix statement.

    Temporary unmask the peak matrix with particular tags. Within the statement the samples can be motified or removed.
    After leaving the statement the original mask will be recoverd.

    :param pm: the target peak matrix
    :param override: whether to override the current mask, default = True
    :param args: target tag values, both typed and untyped
    :param kwargs: target typed tag types and values
    :rtype: PeakMatrix object

    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')
    >>> with unmask_peakmatrix(pm, 'qc') as m: print m.peaklist_ids
    ('qc_1', 'qc_2') # no need to set pm.mask to True
    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')

    """

    def __init__(self, pm, *args, **kwargs):
        self._pm = pm
        self._utags = args
        self._ttags = kwargs
        if not self._ttags.has_key('override'): self._ttags['override'] = True  # default for with statement
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.unmask_tags(*self._utags, **self._ttags)
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]


class mask_all_peakmatrix:
    """
    The mask_all_peakmatrix statement.

    Temporary mask all the peak matrix samples. Within the statement the samples can be motified or removed.
    After leaving the statement the original mask will be recoverd.

    :param pm: the target peak matrix
    :rtype: PeakMatrix object

    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')
    >>> with mask_all_peakmatrix(pm) as m: print m.peaklist_ids
    ()
    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')

    """

    def __init__(self, pm):
        self._pm = pm
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.mask = [True] * self._pm.full_shape[0]
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]


class unmask_all_peakmatrix:
    """
    The unmask_all_peakmatrix statement.

    Temporary unmask all the peak matrix samples. Within the statement the samples can be motified or removed.
    After leaving the statement the original mask will be recoverd.

    :param pm: the target peak matrix
    :rtype: PeakMatrix object

    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')
    >>> with unmask_all_peakmatrix(pm) as m: print m.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')
    >>> print pm.peaklist_ids
    ('sample_1', 'sample_2', 'qc_1', 'sample_3', 'sample_4', 'qc_2')

    """

    def __init__(self, pm):
        self._pm = pm
        self._oldmask = dict(zip(pm._pids, pm._mask))

    def __enter__(self):
        self._pm.mask = None
        return self._pm

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._pm.mask = [self._oldmask[i] for i in self._pm._pids]


# RSD calculation
def peak_matrix_rsd(pm, *args, **kwargs):
    """
    Calculate the QC samples' relative standard deviation (RSD) array for PeakMatrix.

    :param pm: target peak matrix object
    :param args: untyped QC tag label for RSD calculation
    :param kwargs: typed QC tag label for RSD calculation
    :param flagged_only: whether to calculate on flagged peaks only. Default = True
    :type: numpy array

    The RSD is calculated as:

    >>> rsd = std(pm.intensity_matrix, axis = 0) / mean(pm.intensity_matrix, axis = 0) * 100

    Noting that only the "present" peaks will be used for calculation. If a column has less than 2 peaks, the
    corresponding rsd value will be set to np.nan.

    """
    flagged_only = kwargs.pop('flagged_only') if kwargs.has_key('flagged_only') else True

    if pm.shape[0] < 2:
        # logging.warning('calculating RSD on less than 2 samples')
        return np.ones(pm.shape[1] if flagged_only else pm.full_shape[1]) * np.nan

    with (unmask_all_peakmatrix(pm) if len(args) + len(kwargs) == 0 else unmask_peakmatrix(pm, *args, **kwargs)) as m:
        if m.shape[0] == 0:
            raise AttributeError('peak matrix does not have QC labels [%s]' % join(map(lambda x: str(x)[1:-1], (args,kwargs)), ', '))
        ints = m.attr_matrix('intensity', flagged_only)
        rsd = m._present_std(ints, 0, flagged_only) / m._present_mean(ints, 0, flagged_only) * 100

    rsd[np.where(map(lambda x: len(set(x[np.nonzero(x)])) == 1, ints.T))] = np.nan # only one valid value
    return rsd
