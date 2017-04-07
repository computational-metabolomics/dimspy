#!/usr/bin/python
# -*- coding: utf-8 -*-
import setuptools
import sys
import dimspy
import os


def main():

    if sys.version_info[0] != 2 and sys.version_info[1] <= 7:
        sys.exit("Python-2.7.8 is required ")

    required_unix = open('requirements.txt').read().splitlines()
    required_win = open('requirements_win.txt').read().splitlines()

    setuptools.setup(name="dimspy",
        version=dimspy.__version__,
        description="Python package to process DIMS data",
        long_description=open('README.rst').read(),
        author="Ralf Weber, Albert Zhou",
        author_email="r.j.weber@bham.ac.uk, j.zhou.3@bham.ac.uk ",
        url="https://github.com/computational-metabolomics/dimspy",
        license="GPLv3",
        platforms=['Windows, UNIX'],
        keywords=['Metabolomics', 'Mass spectrometry', 'Data Processing', 'Direct-Infusion Mass Spectrometry'],
        packages=setuptools.find_packages(),
        test_suite='tests',
        install_requires=required_unix
            if os.name != "nt" \
            else required_unix.extend(required_win),

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
