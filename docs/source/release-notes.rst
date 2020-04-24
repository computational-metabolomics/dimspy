Release Notes
=============

All notable changes to this project will be documented.
For more details please refer to `github <https://github.com/computational-metabolomics/dimspy>`_ commit history

DIMSpy v2.0.0
-------------

**Release date: 24 April 2020**

- First stable Python 3 only release
- Refactor and improve HDF5 portal to save peaklists and/or peak matrices
- Add compatibility for previous HDF5 files (python 2 version of DIMSpy)
- Improve filelist handling
- mzML or raw files are ordered by timestamp if no filelist is provided (i.e. process_scans)
- Fix warnings (NaturalNameWarning, ResourceWarning, DeprecationWarning)
- Fix 'blank filter' bug (missing and/or zero values are excluded)
- Improve sub setting / filtering of scan events
- Optimise imports
- Increase coverage of tests
- Improve documentation (`Read the Docs <https://dimspy.readthedocs.io/en/latest/>`_), including docstring

DIMSpy v2.0.0
-------------

**Release date: 2 October 2019**

- Final Python 2 release