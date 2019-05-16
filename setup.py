#!/usr/bin/python
# -*- coding: utf-8 -*-
import setuptools
import sys
import dimspy


def main():

    setuptools.setup(name="dimspy",
        version=dimspy.__version__,
        description="Python package for data processing of direct-infusion mass spectrometry-based metabolomics and lipidomics data",
        long_description=open('README.rst').read(),
        author="Ralf Weber, Albert Zhou",
        author_email="r.j.weber@bham.ac.uk, j.zhou.3@bham.ac.uk ",
        url="https://github.com/computational-metabolomics/dimspy",
        license="GPLv3",
        platforms=['Windows, UNIX'],
        keywords=['Metabolomics', 'Lipidomics', 'Mass spectrometry', 'Data Processing', 'Direct-Infusion Mass Spectrometry'],
        packages=setuptools.find_packages(),
        test_suite='tests.suite',
        python_requires='>=3.7',
        install_requires=open('requirements.txt').read().splitlines(),
        include_package_data=True,
        classifiers=[
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.7",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
        ],
        entry_points={
         'console_scripts': [
             'dimspy = dimspy.__main__:main'
         ]
        }
    )


if __name__ == "__main__":
    main()
