#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017-2020 Ralf Weber, Albert Zhou.
#
# This file is part of DIMSpy.
#
# DIMSpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIMSpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DIMSpy.  If not, see <https://www.gnu.org/licenses/>.
#


import setuptools
import sys
import dimspy


def main():

    setuptools.setup(name="dimspy",
        version=dimspy.__version__,
        description="Python package for processing of direct-infusion mass spectrometry-based metabolomics and lipidomics data",
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
        project_urls={
            "Documentation": "https://dimspy.readthedocs.io/en/latest/",
            "Changelog": "https://dimspy.readthedocs.io/en/latest/changelog.html",
            "Bug Tracker": "https://github.com/computational-metabolomics/dimspy/issues",
        },
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
