.. _step_by_ste_install:

********************
Package Installation
********************

Ther are different ways to install ``SDSS_Catl_Utils``. The preferred way
to install it is through `pip <https://pypi.org/>`_. Another possibility
is to install the package through the source code by cloning
the repository from `Github <https://github.com/vcalderon2009/sdss_catl_utils>`_ 
and build from the source code.

``pip`` Installation
====================

The simplest way to install ``SDSS_Catl_Utils`` is with ``pip``.
To install it with `pip`

.. code-block:: bash

    pip install sdss_catl_utils

This will install the latest official release of the code.

Upgrading via ``pip``
=====================

Whenever there is a new release, you can upgrade your current
version by running ::

    pip install --upgrade sdss_catl_utils

This will ensure that you have the most up-to-date version of
``SDSS_Catl_Utils``.

.. note::

    Consider installing ``sdss_catl_utils`` into a virtual environment.
    Setting this up is completely straightforward and takes less than
    a minute, even if this is your first time using a virtual environment.
    Using a virtual environment simplifies not just the current installation
    but also package upgrades and your subsequent workflow.
    If you use `conda <https://www.continuum.io/downloads>`_
    to manage your python distribution, you can find explicit instructions
    in the :ref:`installing_sdss_catl_utils_with_virtualenv`
    section of the documentation.

Building from Source
====================

If you don't install the latest release using ``pip``, you can instead
`clone` the source code and call the setup file.
This is the most common way to install ``SDSS_Catl_Utils`` if you want
versions of the code that have been updated since the last latest
official release. In this case, after installation, it is particularly
important that you follow the instructions in :ref:`verifying_installation`
section bellow.

The first step is to clone the ``SDSS_Catl_Utils`` repository ::

    git clone https://github.com/vcalderon2009/sdss_catl_utils
    cd sdss_catl_utils

Installing one of the official releases
----------------------------------------

All oficial releases of the code are tagged with their version name,
e.g. `v0.1.0` pertains to the ``0.1.0`` release.
To install a particular release ::

    git checkout v0.1.0
    python setup.py install

This will install the ``v0.1.0`` release of the ``SDSS_Catl_Utils``.
Other official release version can be installed similarly.

Installing the most recent master branch
----------------------------------------

If you prefer to use the most recent version of the code ::

    git checkout master
    python setup.py install

This will install the ``master`` branch of the code, which is 
the branch currently under development. While the features in the
officiail releases have a stable API, new features being developed
in the ``master`` branch may not. However, the master branch may have
*new* features and/or perfomance enhancements that you may wish to use
for your science applications. A concerted effort is made to ensure
that only thoroughly tested and documented code appears in the public
``master`` branch, though ``SDSS_Catl_Utils`` users should be awayre
of the distinction between the bleeding edge version in master
and the official release version available through ``pip``.

.. _sdss_catl_utils_dependencies:

Dependencies
============

If your install ``sdss_catl_utils`` using pip, then most of your dependencies
will be handled for you automatically. The only additional dependency
you may need is:

* `h5py <http://h5py.org/>`__ : 2.5 or later

The h5py package is used for fast I/O of galaxy and group catalogues.

If you did not use pip, then you should be aware of the following strict
requirements:

* `Python <http://www.python.org/>`_: 3x or higher 
* `Numpy <http://www.numpy.org>`_ 
* `Astropy <http://www.astropy.org/>`__ 
* `Pandas <https://pandas.pydata.org/>`_ 
* `h5py <https://www.h5py.org/>`__ 
* `GitPython <https://gitpython.readthedocs.io/en/stable/>`_ 
* `Cython <https://cython.org/>`_ 
* `requests <http://docs.python-requests.org/en/master/>`_ 
* `numexpr <https://github.com/pydata/numexpr>`_ 
* `Scipy <https://www.scipy.org/>`_ 
* `Scikit-Learn <https://scikit-learn.org>`_ 
* `BeautifulSoup <https://www.crummy.com/software/BeautifulSoup/>`_ 
* `wget <https://bitbucket.org/techtonik/python-wget/src>`_ 
* `tqdm <https://tqdm.github.io/>`_ 
* `cosmo-utils <https://github.com/vcalderon2009/cosmo_utils>`_ 

Any of the above can be installed with either `pip` or `conda`.

.. _verifying_installation:

Verifying your installation
===========================

After installing the code and its dependencies, fire up a Python interpreter
and check that the version number matches what you expect:

.. code-block:: python

    >>> import sdss_catl_utils
    >>> print(sdss_catl_utils.__version__) # doctest: +SKIP

If the version number is not what it should be, this likely means you have a 
previous installation that is superseding the version you tried to install.
This *should* be accomplished by doing ``pip uninstall sdss_catl_utils``
before your new installation, but you may need to uninstall the previous 
build "manually". Like all python packages, you can find the installation 
location as follows:

.. code-block:: python

    >>> import sdss_catl_utils
    >>> print(sdss_catl_utils.__file__) # doctest: +SKIP

This wil show where your active version is located on your machine. You 
can manually delete this copy of ``SDSS_Catl_Utils`` prior to your new
installation to avoid version conflicts. (There may be multiple copies
of ``SDSS_Catl_Utils`` in this location, depending on how many times
you have previously installed the code - all such copies my be deleted
prior to reinstallation).

