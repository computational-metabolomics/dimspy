#!/usr/bin/env python
# 2016 copyright Ralf Weber, Albert Zhou
# released under the GNU General Public License version 3.0 (GPLv3)

from setuptools import setup
from setuptools import find_packages
import sys
import os
import dimspy


def main():
    if sys.version_info[0] != 2 and sys.version_info[1] <= 7:
        sys.exit("Python-2.7.8 is required ")

    setup(name="dimspy",
          version=dimspy.__version__,
          description="Python package to process DIMS data",
          long_description=open('README.rst').read(),
          author="Ralf Weber, Albert Zhou",
          author_email="r.j.weber@bham.ac.uk, j.zhou.3@bham.ac.uk ",
          url="https://github.com/computational-metabolomics/dimspy",
          license="GPLv3",
          platforms=['Windows, UNIX'],
          keywords=['Metabolomics', 'Mass spectrometry', 'Data Processing', 'Direct-Infusion Mass Spectrometry'],
          packages=find_packages(),

          test_suite='tests',

          install_requires=open('requirements.txt').read().splitlines()
                            if os.name != "nt" \
                            else open('requirements-win.txt').read().splitlines(),
          
          include_package_data=True,

          classifiers=[
              "Programming Language :: Python :: 2",
              "Programming Language :: Python :: 2.7",
              "Topic :: Scientific/Engineering :: Bio-Informatics",
              "Topic :: Scientific/Engineering :: Chemistry",
              "Topic :: Utilities",
              "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
              "Operating System :: OS Independent",
          ])


if __name__ == "__main__":
    main()
