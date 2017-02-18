import os, sys
from string import join

sys.path.append(os.path.abspath(os.path.dirname(__file__)).rsplit(os.sep,1)[0])
from dimspy.models.peaklist import _loadPeaklist
from dimspy.process.peak_alignment import align_peaks


if __name__ == '__main__':
    attrdct = {'M/Z': 'mzs', 'INTENSITY': 'ints', 'SNR': 'snr', 'NON-NOISE_FLAG': 'non_noise_flag'}
    pkls = map(lambda x: _loadPeaklist(x, x.rsplit('.', 1)[0], attrdct, ('NON-NOISE_FLAG',)), sys.argv[1:])

    pm = align_peaks(pkls, ppm = 2.0, block_size = 2000, byunique = True, procs = 2)

    with open('test_alignment.txt', 'w') as f:
        f.write(join(map(lambda ln: join(map(str, ln), '\t'), pm.intensity_matrix), '\n'))
