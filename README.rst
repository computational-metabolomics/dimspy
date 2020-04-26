DIMSpy
======
|Py versions| |Version| |Bioconda| |Galaxy-eu| |Git| |Build Status (Travis)| |Build Status (AppVeyor)| |codecov| |License| |binder| |RTD doc| |gitter|

Python package for processing direct-infusion mass spectrometry-based metabolomics and lipidomics data

- **Documentation:** https://dimspy.readthedocs.io/en/latest
- **Source:** https://github.com/computational-metabolomics/dimspy
- **Bug reports:** https://github.com/computational-metabolomics/dimspy/issues


Installation (Conda, PyPi and Galaxy)
-------------------------------------
See the `Installation page <https://dimspy.readthedocs.io/en/latest/installation.html>`__ of
the `online documentation <https://computational-metabolomics.github.io/dimspy/>`__.


Bug reports
-----------
Please report any bugs that you find `here <https://github.com/computational-metabolomics/dimspy/issues>`_.
Or fork the repository on `GitHub <https://github.com/computational-metabolomics/dimspy/>`_
and create a pull request (PR). We welcome all contributions, and we
will help you to make the PR if you are new to `git`.


Credits
-------
DIMSpy was originally written by Ralf Weber and Albert Zhou and has been developed with the help of many others.
Thanks to everyone who has improved DIMSpy by contributing code, adding features, bug reports and fixes, and documentation.

**Developers and contributers**
 - Ralf J. M. Weber (r.j.weber@bham.ac.uk) - `University of Birmingham (UK) <https://www.birmingham.ac.uk/staff/profiles/biosciences/weber-ralf.aspx>`__
 - Jiarui (Albert) Zhou (j.zhou.3@bham.ac.uk) - `University of Birmingham (UK) <http://www.birmingham.ac.uk/index.aspx>`__, `HIT Shenzhen (China) <http://www.hitsz.edu.cn>`_
 - Thomas N. Lawson (t.n.lawson@bham.ac.uk) - `University of Birmingham (UK) <http://www.birmingham.ac.uk/index.aspx>`__
 - Martin R. Jones (martin.jones@eawag.ch) - `Eawag  (Switzerland) <https://www.eawag.ch/en/aboutus/portrait/organisation/staff/profile/martin-jones/show/>`_

**DIMSpy acknowledges support from the following funders:**
 - BBSRC, grant number BB/M019985/1
 - European Commission's H2020 programme, grant agreement number 654241
 - Wellcome Trust, grant number 202952/Z/16/Z

**Citation**

To cite DIMSpy please use one of the Zenodo references listed `here <https://dimspy.readthedocs.io/en/latest/citation.html>`__.


License
--------
DIMSpy is licensed under the GNU General Public License v3.0 (see `LICENSE file <https://github.com/computational-metabolomics/dimspy/blob/master/LICENSE>`_ for licensing information). Copyright © 2017 - 2020 Ralf Weber, Albert Zhou

**Third-party licenses and copyright**

RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved. See `RawFileReaderLicense <https://github.com/computational-metabolomics/dimspy/blob/master/RawFileReaderLicense.rst>`_ for licensing information.
Using DIMSpy software for processing Thermo Fisher Scientific \*.raw files implies the acceptance of the RawFileReader license terms.
Anyone receiving RawFileReader as part of a larger software distribution (in the current context, as part of DIMSpy) is considered an "end user" under
section 3.3 of the RawFileReader License, and is not granted rights to redistribute RawFileReader.


.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/dimspy.svg?logo=travis&maxAge=600&style=flat-square
   :target: https://travis-ci.com/computational-metabolomics/dimspy

.. |Build Status (AppVeyor)| image:: https://img.shields.io/appveyor/ci/RJMW/dimspy.svg?logo=appveyor&style=flat-square&maxAge=600
   :target: https://ci.appveyor.com/project/RJMW/dimspy/branch/master

.. |Py versions| image:: https://img.shields.io/pypi/pyversions/dimspy.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/dimspy/

.. |Version| image:: https://img.shields.io/pypi/v/dimspy.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/dimspy/

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/dimspy

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/dimspy?style=flat-square&maxAge=3600
   :target: http://bioconda.github.io/recipes/dimspy/README.html
   
.. |galaxy-eu| image:: https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==
   :target: http://usegalaxy.eu

.. |License| image:: https://img.shields.io/pypi/l/dimspy.svg?style=flat&maxAge=3600
   :target: https://www.gnu.org/licenses/gpl-3.0.html

.. |RTD doc| image:: https://img.shields.io/badge/documentation-RTD-71B360.svg?style=flat&maxAge=3600
   :target: https://dimspy.readthedocs.io/en/latest/
   
.. |codecov| image:: https://codecov.io/gh/computational-metabolomics/dimspy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/computational-metabolomics/dimspy

.. |binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/computational-metabolomics/dimspy/master?filepath=notebooks%2Fworkflow.ipynb

.. |gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/computational-metabolomics/dimspy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
