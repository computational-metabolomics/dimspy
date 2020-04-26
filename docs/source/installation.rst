Installation
============

Conda (recommended)
-------------------

Install Miniconda, follow the steps described `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install>`_

Start the ``conda prompt``

* Windows: Open the ``Anaconda Prompt`` via the Start menu
* macOS or Linux: Open a ``Terminal``

Create a dimspy specific ``conda`` environment.
This will install a the dependencies required to run ``dimspy``::

    $ conda create --yes --name dimspy -c conda-forge -c bioconda -c computational-metabolomics

.. note::

    * The installation process will take a few minutes.
    * Feel free to use a different name for the Conda environment

    You can use the following command to remove a conda environment::

        $ conda env remove -y --name dimspy

    This is only required if something has gone wrong in the previous step.

Activate the ``dimspy`` environment::

    $ conda activate dimspy

To test your ``dimspy`` installation, in your Conda Prompt, run the command::

    $ dimspy --help

or::

    $ python
    import dimspy

Close and deactivate the ``dimspy`` environment when you’re done::

    $ conda deactivate


PyPi
----

Install the current release of ``dimspy`` with ``pip``::

    $ pip install dimspy

.. note::

    * The installation process will take a few minutes.

To upgrade to a newer release use the ``--upgrade`` flag::

    $ pip install --upgrade dimspy

If you do not have permission to install software systemwide, you can
install into your user directory using the ``--user`` flag::

    $ pip install --user dimspy

Alternatively, you can manually download ``dimspy`` from
`GitHub <https://github.com/computational-metabolomics/dimspy/releases>`_  or
`PyPI <https://pypi.python.org/pypi/dimspy>`_.
To install one of these versions, unpack it and run the following from the
top-level source directory using the Terminal::

    $ pip install .


Testing
-------
DIMSpy uses the Python ``pytest`` testing package.  You can learn more
about pytest on their `homepage <https://pytest.org>`_.