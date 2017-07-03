dimspy
======
Python Package to process Direct-Infusion Mass Spectrometry (DIMS) data

|Version| |Py versions| |Git| |Bioconda| |Build Status (Travis)| |Build Status (AppVeyor)| |License| |RTD doc|

Overview
--------

Install
--------


Conda_
~~~~~~~

1. Install Conda_ (For example: `Miniconda Python distribution <http://conda.pydata.org/miniconda.html>`__).
2. Run the following commands to install DIMSpy.

Linux-64

::

    $ conda create -n dimspy python=2.7 numpy=1.13.0 scipy=0.19.1 pymzml=0.7.8 pythonnet=2.3.0 h5py=2.7.0 fastcluster=1.1.23 -c conda-forge -c bioconda
    $ source activate dimspy  
    $ pip install git+https://github.com/computational-metabolomics/dimspy.git


OSx (TO WRITE & TEST)

::

    $ 

Windows-64

::

    $ conda create -n dimspy python=2.7 numpy=1.13.0 scipy=0.19.1 fastcluster=1.1.23 h5py==2.7.0 -c conda-forge -c bioconda
    $ activate dimspy
    $ pip install pythonnet==2.3.0 pymzml==0.7.8
    $ pip install git+https://github.com/computational-metabolomics/dimspy.git


Usage
------

Command line
~~~~~~~~~~~~~

::

    $ python -m dimspy --help


Workflow
---------
TODO


.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/dimspy.svg?style=flat&maxAge=3600&label=Travis-CI
   :target: https://travis-ci.org/computational-metabolomics/dimspy

.. |Build Status (AppVeyor)| image:: https://img.shields.io/appveyor/ci/computational-metabolomics/mzml2isa.svg?style=flat&maxAge=3600&label=AppVeyor
   :target: https://ci.appveyor.com/project/computational-metabolomics/dimspy

.. |Py versions| image:: https://img.shields.io/pypi/pyversions/dimspy.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/dimspy/

.. |Version| image:: https://img.shields.io/pypi/v/dimspy.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/dimspy/

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/dimspy

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: http://bioconda.github.io/recipes/dimspy/README.html

.. |License| image:: https://img.shields.io/pypi/l/dimspy.svg?style=flat&maxAge=3600
   :target: https://www.gnu.org/licenses/gpl-3.0.html

.. |RTD doc| image:: https://img.shields.io/badge/documentation-RTD-71B360.svg?style=flat&maxAge=3600
   :target: http://dimspy.readthedocs.io/en/latest/dimspy/index.html

.. _pip: https://pip.pypa.io/
.. _Conda: http://conda.pydata.org/docs/
