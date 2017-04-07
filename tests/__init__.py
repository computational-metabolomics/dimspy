#!/usr/bin/python
# -*- coding: utf-8 -*-
import unittest
import pprint
import numpy as np
from dimspy.models import PeakMatrix
from dimspy.models import PeakList
from dimspy.models import mask_peakmatrix

pp = pprint.PrettyPrinter(indent=4)
np.random.seed(1234)

class TestDimspy(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_models_peaklist(self):
        # generate random peak
        size = 100
        mzs = sorted(np.random.uniform(100, 1000, size=size))
        ints = np.abs(np.random.normal(10, 3, size=size))
        snr = np.abs(np.random.normal(1000, 400, size=size))

        pl = PeakList('sample_peaklist', mzs, ints, mz_range=(100, 1000))
        pl.add_tags('sample', 'passed_qc', main_class='treatment_1')
        pl.add_attribute('snr', snr)
        pl.metadata.type = 'blank'

        """
        import cPickle as cp

        s = cp.dumps(pl)
        pl = cp.loads(s)

        # magic functions
        print len(pl)
        print pl

        # tags
        print pl.tags
        print pl.tag_types
        print pl.tag_values
        pl.add_tags('not_QC')
        pl.add_tags(time_point='1hr')
        print pl.tag_of('time_point')
        print pl.has_tags('1hr')
        print pl.has_any_tags(main_class='1hr')
        print pl.has_tag_type('main_class')
        pl.drop_tags('not_QC')
        pl.drop_tag_types('time_point')
        print pl.tags.to_list()

        # metadata
        pl.metadata['replicate'] = '4'
        pl.metadata.sample = 'S01'
        print pl.metadata
        del pl.metadata.sample
        print pl.metadata
        print pl.metadata.keys()

        # add flag
        print pl.flags
        pl.add_attribute('odd_flags', np.arange(len(pl)) % 2 == 0, is_flag=True)
        print pl.flags

        # attrs access
        print pl.mz  # internally == getAttribute
        print pl.snr
        print pl.mz_all
        print pl.snr_all

        pl.mz = np.random.uniform(100, 1000, size=len(pl))  # internally == setAttribute
        print pl[:10]

        print pl['mz']  # internally == getPeak
        print pl['snr']

        print pl[0]  # internally == getPeak
        pl[0] = (0, 100, 200, True)  # internally == setPeak
        print pl[:10]

        pl[0][0] = 1  # not working
        print pl[:10]

        pl['mz'][0] = 2  # not working
        print pl[:10]

        # properties
        print pl.ID
        print pl.size
        print pl.shape
        print pl.attributes
        print pl.peaks
        print pl.dtable

        # publics
        print pl.has_attribute('snr')
        print pl.has_attribute('no_such_thing')
        pl.add_attribute('snr_filter_flag', pl.snr > 20, is_flag=True)  # must be all full size vector
        print pl
        pl.drop_attribute('snr_filter_flag')
        print pl

        pl.insert_peak((900, 20, 30, True))  # auto sort
        print pl
        pl.remove_peak(-2)
        print pl

        print pl.to_list()[0]
        print pl.to_dict().keys()
        """

    def test_models_peak_matrix(self):
        # generate random peaklists

        _mzs = lambda: sorted(np.random.uniform(100, 1000, size=100))
        _ints = lambda: np.abs(np.random.normal(10, 3, size=100))

        pls = [
            PeakList('sample_1_1', _mzs(), _ints(), mz_range=(100, 1000)),
            PeakList('sample_1_2', _mzs(), _ints(), mz_range=(100, 1000)),
            PeakList('QC_1', _mzs(), _ints(), mz_range=(100, 1000)),
            PeakList('sample_2_1', _mzs(), _ints(), mz_range=(100, 1000)),
            PeakList('sample_2_2', _mzs(), _ints(), mz_range=(100, 1000)),
            PeakList('QC_2', _mzs(), _ints(), mz_range=(100, 1000)),
        ]

        pls[0].add_tags('sample', treatment='compound_1', time_point='1hr', plate=1)
        pls[1].add_tags('sample', treatment='compound_1', time_point='6hr', plate=1)
        pls[2].add_tags('qc', plate=1)
        pls[3].add_tags('sample', treatment='compound_2', time_point='1hr', plate=2)
        pls[4].add_tags('sample', treatment='compound_2', time_point='6hr', plate=2)
        pls[5].add_tags('qc', plate=2)

        # create matrix
        pm = PeakMatrix(
            [p.ID for p in pls],
            [p.tags for p in pls],
            {a: np.vstack([p.get_attribute(a) for p in pls]) for a in pls[0].attributes}
        )

        # properties
        np.testing.assert_array_equal(pm.mask, [True,  True,  True,  True,  True,  True])
        np.testing.assert_array_equal(pm.peaklist_ids, ['sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'sample_2_2', 'QC_2'])
        np.testing.assert_array_equal(pm.peaklist_tag_types, ['plate', 'treatment', 'time_point'])
        np.testing.assert_array_equal(pm.peaklist_tag_values, ['1hr', 1, 2, '6hr', 'sample', 'qc', 'compound_1', 'compound_2'])
        np.testing.assert_array_equal(pm.attributes, ['intensity', 'mz'])
        np.testing.assert_array_equal(pm.shape, (6, 100L))

        np.testing.assert_array_equal(pm.present, [6] * 100)
        np.testing.assert_array_equal(pm.missing, [0] * 6)
        # self.assertEqual(np.isclose(pm.rsd[0:5], np.array([20.2078199, 37.63703196,
        #                                                                27.49539753, 33.07843295,
        #                                                                17.95949877])), True)  # TODO: replace isclose

        np.testing.assert_array_equal(pm.tags_of('plate'), [1, 2])
        # print pm.mz_matrix
        # pm.mask_tags('sample').unmask_tags(treatment = 'compound_1')
        with mask_peakmatrix(pm, 'sample'):
            np.testing.assert_array_equal(pm.mask, [True,  True, False,  True,  True, False])
            np.testing.assert_array_equal(pm.peaklist_ids, ['sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2'])

            pm.remove_peaks(np.where(pm.rsd > 30))

            np.testing.assert_array_equal(pm.mask, [True,  True, False,  True,  True, False])

        pm.remove_samples(np.where(pm.peaklist_ids == 'sample_2_2')[0])
        np.testing.assert_array_equal(pm.peaklist_ids, ['sample_1_1', 'sample_1_2', 'QC_1', 'sample_2_1', 'QC_2'])
        np.testing.assert_array_equal(pm.shape, (5, 70L))
        pm.remove_peaks([0, 1, 2])
        np.testing.assert_array_equal(pm.shape, (5, 67L))

        # print pm.to_peaklist('merged')  # TODO:
        # print pm.to_str(transpose=True)  # TODO:

    def test_parsers_mzml(self):
        pass

    def test_parsers_thermo_raw(self):
        pass

    def test_process(self):
        pass

    def test_peak_alignment(self):
        pass

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDimspy)
    unittest.TextTestRunner(verbosity=2).run(suite)

