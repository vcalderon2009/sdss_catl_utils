#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018_12-18
# Last Modified: 2018_12-18
# Vanderbilt University
from __future__ import absolute_import, division, print_function 
__author__     = ['Victor Calderon']
__copyright__  = ["Copyright 2018 Victor Calderon, "]
__email__      = ['victor.calderon@vanderbilt.edu']
__maintainer__ = ['Victor Calderon']
__all__        = [  'CatlUtils']


import os
from cosmo_utils.utils import file_utils      as cfutils
from cosmo_utils.utils import work_paths      as cwpaths
from cosmo_utils.utils import web_utils       as cweb

from sdss_catl_utils.mocks_manager import mocks_defaults as md

## Functions and classes

class CatlUtils(object):
    """
    Class used to handle the galaxy/group galaxy/group catalogues, *after*
    the catalogues have been downloaded. This class has functions to read,
    modify, and analyze the galaxy/group catalogues.

    For a list of available pre-processed galaxy- and group-galaxy
    catalogues provided by ``sdss_catl_utils``, see ...
    """
    def __init__(self, **kwargs):
        """
        Parameters
        ------------
        catl_kind : {``data``, ``mocks``} `str`
            Kind of catalogues to download. This variable is set to
            ``mocks`` by default.

            Options:
                - ``data``: Downloads the SDSS DR7 real catalogues.
                - ``mocks``: Downloads the synthetic catalogues of SDSS DR7.

        hod_n : `int`, optional
            Number of the HOD model to use. This value is set to `0` by
            default.

        halotype : {'so', 'fof'}, `str`, optional
            Type of dark matter definition to use. This value is set to
            ``so`` by default.

            Options:
                - ``so``: Spherical Overdensity halo definition.
                - ``fof``: Friends-of-Friends halo definition.

        clf_method : {1, 2, 3}, `int`, optional
            Method for assigning galaxy properties to mock galaxies.
            This variable dictates how galaxies are assigned
            luminosities or stellar masses based on their galaxy type
            and host halo's mass. This variable is set to ``1`` by
            default.

            Options:
                - ``1``: Independent assignment of (g-r) colour, sersic, and specific star formation rate (`logssfr`)
                - ``2``: (g-r) colour dictates active/passive designation and draws values independently.
                - ``3``: (g-r) colour dictates active/passive designation, and assigns other galaxy properties for that given galaxy.

        clf_seed : `int`, optional
            Value of the random seed used for the conditional luminosity function.
            This variable is set to ``1235`` default.

        dv : `float`, optional
            Value for the ``velocity bias`` parameter. It is the difference
            between the galaxy and matter velocity profiles.

            .. math::
                dv = \\frac{v_{g} - v_{c}}{v_{m} - v_{c}}

            where :math:`v_g` is the galaxy's velocity; :math:`v_m`, the
            matter velocity.

        sample : {'19', '20', '21'}, `str`, optional
            Luminosity of the SDSS volume-limited sample to analyze.
            This variable is set to ``'19'`` by default.

            Options:
                - ``'19'``: :math:`M_r = 19` volume-limited sample
                - ``'20'``: :math:`M_r = 20` volume-limited sample
                - ``'21'``: :math:`M_r = 21` volume-limited sample

        catl_type : {'mr', 'mstar'}, `str`, optional
            Type of Abundance matching used in the catalogue. This
            variable is set to ``'mr'`` by default.

            Options:
                - ``'mr'``: Luminosity-based abundance matching used
                - ``'mstar'``: Stellar-mass-based abundance matching used.

        cosmo_choice : { ``'LasDamas'``, ``'Planck'``} `str`, optional
            Choice of cosmology to use. This variable is set to ``LasDamas``
            by default.

            Options:
                - ``LasDamas`` : Uses the cosmological parameters from the 
                    `LasDamas <http://lss.phy.vanderbilt.edu/lasdamas/simulations.html>`_ simulations.
                - ``Planck`` : Uses the Planck 2015 cosmology.

        perf_opt : `bool`, optional
            If True, it chooses to analyze the ``perfect`` version of
            the synthetic galaxy/group galaxy catalogues. Otherwise,
            it downloads the catalogues with group-finding errors
            included. This variable is set to ``False`` by default.

        environ_name : `str`
            Name of the environment variable to assign to ``outdir``.
            This variable is set to the default ``environ_name`` from
            `~sdss_catl_utils.mocks_manager.mocks_default`

        Examples
        ----------

        >>> from sdss_catl_utils.mocks_manager.download_manager import DownloadManager
        """
        super().__init__()
        # Assigning variables
        self.catl_kind    = kwargs.get('catl_kind', md.catl_kind)
        self.hod_n        = kwargs.get('hod_n', md.hod_n)
        self.halotype     = kwargs.get('halotype', md.halotype)
        self.clf_method   = kwargs.get('clf_method', md.clf_method)
        self.clf_seed     = kwargs.get('clf_seed', md.clf_seed)
        self.dv           = kwargs.get('dv', md.dv)
        self.sample       = kwargs.get('sample', md.sample)
        self.type_am      = kwargs.get('type_am', md.type_am)
        self.cosmo_choice = kwargs.get('cosmo_choice', md.cosmo_choice)
        self.perf_opt     = kwargs.get('perf_opt', md.perf_opt)
        self.remove_files = kwargs.get('remove_files', md.remove_files)
        self.environ_name = kwargs.get('environ_name', md.environ_name)
        # Other variables
        self.sample_Mr    = 'Mr{0}'.format(self.sample)
        self.sample_s     = str(self.sample)

    # Main directory path - Path to which all catalogues are saved
    def main_dir(self):
        """
        Path to the main directory, in which all catalogues are stored.

        Returns
        ---------
        main_dirpath : `str`
            Path to the main directory, in which all catalogues are stored.
            It uses the environment variable ``environ_name`` from
            the class.
        """
        if os.environ.get(self.environ_name):
            main_dirpath = os.path.join(os.environ[self.environ_name],
                                        'data',
                                        'external',
                                        'SDSS/')
            # Checking directory exists
            if not os.path.exists(main_dirpath):
                msg = '`main_dirpath` ({0}) does not exist!'.format(
                    main_dirpath)
                raise ValueError(msg)
        else:
            msg = '`environ_name` ({0}) is not a system variable! '
            msg += 'Please download the catalogues and add the path '
            msg += 'to the catalogues as `{0}` environment variable!'
            msg = msg.format(self.environ_name)
            raise ValueError(msg)

        return main_dirpath

    # Location prefix
    def _catl_prefix(self, catl_type='memb', catl_kind='mocks',
        perf_opt=False):
        """
        Prefix of the paths based on the type of catalogues and input
        parameters chosen.

        Parameters
        ------------
        catl_type : {'gal', 'memb', 'group'}, `str`
            Option for which kind of catalogue is being analyzed. This
            variable is set to ``memb`` by default.

            Options:
                - ``'gal'`` : Galaxy catalogue
                - ``'memb'`` : Group Member galaxy catalogue
                - ``'group'`` : Group galaxy catalogue

        catl_kind : {``data``, ``mocks``} `str`
            Kind of catalogues to download. This variable is set to
            ``mocks`` by default.

            Options:
                - ``data``: Downloads the SDSS DR7 real catalogues.
                - ``mocks``: Downloads the synthetic catalogues of SDSS DR7.

        perf_opt : `bool`, optional
            If True, it chooses to analyze the ``perfect`` version of
            the synthetic galaxy/group galaxy catalogues. Otherwise,
            it downloads the catalogues with group-finding errors
            included. This variable is set to ``False`` by default.

        Returns
        ----------
        catl_prefix : `str`
            Prefix of the paths based on the type of catalogues and
            input parameters.
        """
        ## Checking input parameters
        # `catl_type` - Input variable
        catl_type_arr = ['gal', 'memb', 'group']
        if not (catl_type in catl_type_arr):
            msg = '>>> `catl_type` ({0}) is not a valid input!'
            msg = msg.format(catl_type)
            raise ValueError(msg)
        # `catl_type` - Type
        if not (isinstance(catl_type, str)):
            msg = '>>> `catl_type` ({0}) is not a valid type!'
            msg = msg.format(type(catl_type))
            raise TypeError(msg)
        # `catl_kind` - Input variable
        catl_kind_arr = ['data', 'mocks']
        if not (catl_kind in catl_kind_arr):
            msg = '>>> `catl_kind` ({0}) is not a valid input!'
            msg = msg.format(catl_kind)
            raise ValueError(msg)
        # `catl_kind` - Type
        if not (isinstance(catl_kind, str)):
            msg = '>>> `catl_kind` ({0}) is not a valid type!'
            msg = msg.format(type(catl_kind))
            raise TypeError(msg)
        # `perf_opt` - Type
        if not (isinstance(perf_opt, bool)):
            msg = '`perf_opt` ({0}) is not a valid type!'
            msg = msg.format(type(perf_opt))
            raise TypeError(msg)
        ##
        ## Option for which type of catalogue to download
        catl_type_dict = {  'gal'  : 'galaxy_catalogues',
                            'memb' : 'member_galaxy_catalogues',
                            'group': 'group_galaxy_catalogues',
                            'perf_group': 'perfect_group_galaxy_catalogues',
                            'perf_memb' : 'perfect_member_galaxy_catalogues'}
        # Deciding which folder to use
        if (catl_type == 'gal'):
            catl_type_str = 'gal'
        elif (catl_type in ['memb', 'group']):
            if perf_opt:
                # Member galaxies
                if (catl_type == 'memb'):
                    catl_type_str = 'perf_memb'
                # Groups
                elif (catl_type == 'group'):
                    catl_type_str = 'perf_group'
            else:
                # Member galaxies
                if (catl_type == 'memb'):
                    catl_type_str = 'memb'
                # Groups
                elif (catl_type == 'group'):
                    catl_type_str = 'group'
        ##
        ## Parsing prefix path
        # `Data`
        if (catl_kind == 'data'):
            catl_prefix = os.path.join( 'data',
                                        self.type_am,
                                        self.sample_Mr,
                                        catl_type_dict[catl_type_str])
        # `Mocks`
        if (catl_kind == 'mocks'):
            catl_prefix = os.path.join(
                                    'mocks',
                                    'halos_{0}'.format(self.halotype),
                                    'dv_{0}'.format(self.dv),
                                    'hod_model_{0}'.format(self.hod_n),
                                    'clf_seed_{0}'.format(self.clf_seed),
                                    'clf_method_{0}'.format(self.clf_method),
                                    self.type_am,
                                    self.sample_Mr,
                                    catl_type_dict[catl_type_str])

        return catl_prefix

    # Location of galaxy/group catalogues with specified parameters
    def catls_dir(self, catl_type='memb', catl_kind='mocks',
        print_filedir=False):
        """
        Location of the group and group galaxy catalogues with the specified
        parameters.

        Parameters
        ------------
        catl_type : {'gal', 'memb', 'group'}, `str`
            Option for which kind of catalogue is being analyzed. This
            variable is set to ``memb`` by default.

            Options:
                - ``'gal'`` : Galaxy catalogue
                - ``'memb'`` : Group Member galaxy catalogue
                - ``'group'`` : Group galaxy catalogue

        catl_kind : {``data``, ``mocks``} `str`
            Kind of catalogues to download. This variable is set to
            ``mocks`` by default.

            Options:
                - ``data``: Downloads the SDSS DR7 real catalogues.
                - ``mocks``: Downloads the synthetic catalogues of SDSS DR7.

        print_filedir : `bool`, optional
            If True, the path of the catalogue directory is printed
            onto the screen. This variable is set to ``False`` by default.

        Returns
        ---------
        catls_dirpath : `str`
            Path to the location of the group and galaxy catalogues with
            the specified parameters.
        """
        ## Checking input parameters
        # `catl_type` - Value
        catl_type_arr = ['gal', 'memb', 'group']
        if not (catl_type in catl_type_arr):
            msg = '>>> `catl_type` ({0}) is not a valid input!'
            msg = msg.format(catl_type)
            raise ValueError(msg)
        # `catl_type` - Type
        if not (isinstance(catl_type, str)):
            msg = '>>> `catl_type` ({0}) is not a valid type!'
            msg = msg.format(type(catl_type))
            raise TypeError(msg)
        # ``catl_kind`` - Input Value
        catl_kind_arr = ['data', 'mocks']
        if not (catl_kind in catl_kind):
            msg = '>>> `catl_kind` ({0}) is not a valid input!'
            msg = msg.format(catl_kind)
            raise ValueError(msg)
        # ``catl_kind`` - Type
        if not (isinstance(catl_kind, str)):
            msg = '>>> `catl_kind` ({0}) is not a valid type!'
            msg = msg.format(type(catl_kind))
            raise TypeError(msg)
        # `print_filedir` - Type
        if not (isinstance(print_filedir, bool)):
            msg = '`print_filedir` ({0}) is not a valid type!'
            msg = msg.format(type(print_filedir))
            raise TypeError(msg)
        ##
        ## Reading in environment variable
        if os.environ.get(self.environ_name):
            base_dir = self.main_dir()
            assert(os.path.exists(base_dir))
            # Parsing path to catalogues
            catls_dirpath = os.path.join(   base_dir,
                                            self._catl_prefix(
                                                catl_type=catl_type,
                                                catl_kind=catl_kind,
                                                perf_opt=self.perf_opt))
            # Check directory exists
            if not (os.path.exists(catls_dirpath)):
                msg = '>> `catls_dirpath` ({0}) does not exist!'
                msg = msg.format(catls_dirpath)
                raise ValueError(msg)

        else:
            msg = '`environ_name` ({0}) is not a system variable! '
            msg += 'Please download the catalogues and add the path '
            msg += 'to the catalogues as `{0}` environment variable!'
            msg = msg.format(self.environ_name)
            raise ValueError(msg)

        if print_filedir:
            print(catls_dirpath)

        return catls_dirpath

    # Extract the list of catalogues
    def catl_arr_calc(self, catl_type='memb', catl_kind='mocks',
        ext='hdf5', print_filedir=False, return_len=False):
        """
        Extracts the list of galaxy/group catalogues.

        Parameters
        --------------
        catl_type : {'gal', 'memb', 'group'}, `str`
            Option for which kind of catalogue is being analyzed. This
            variable is set to ``memb`` by default.

            Options:
                - ``'gal'`` : Galaxy catalogue
                - ``'memb'`` : Group Member galaxy catalogue
                - ``'group'`` : Group galaxy catalogue

        catl_kind : {``data``, ``mocks``} `str`
            Kind of catalogues to download. This variable is set to
            ``mocks`` by default.

            Options:
                - ``data``: Downloads the SDSS DR7 real catalogues.
                - ``mocks``: Downloads the synthetic catalogues of SDSS DR7.

        ext : {'hdf5'} `str`
            File extension used for the catalogues. This variable is set
            to ``hdf5`` by default.

        print_filedir : `bool`, optional
            If True, the path of the catalogue directory is printed
            onto the screen. This variable is set to ``False`` by default.


        Returns
        ---------
        catls_arr : `numpy.ndarray`, shape (N,)
            Array of the paths of the galaxy/group catalogues with specified
            parameters. The shape of the variable is (`N`,), where `N`,
            is the number of catalogues in the directory.
        """
        ## Checking input parameters
        # `catl_type` - Value
        catl_type_arr = ['gal', 'memb', 'group']
        if not (catl_type in catl_type_arr):
            msg = '>>> `catl_type` ({0}) is not a valid input!'
            msg = msg.format(catl_type)
            raise ValueError(msg)
        # `catl_type` - Type
        if not (isinstance(catl_type, str)):
            msg = '>>> `catl_type` ({0}) is not a valid type!'
            msg = msg.format(type(catl_type))
            raise TypeError(msg)
        # ``catl_kind`` - Input Value
        catl_kind_arr = ['data', 'mocks']
        if not (catl_kind in catl_kind):
            msg = '>>> `catl_kind` ({0}) is not a valid input!'
            msg = msg.format(catl_kind)
            raise ValueError(msg)
        # ``catl_kind`` - Type
        if not (isinstance(catl_kind, str)):
            msg = '>>> `catl_kind` ({0}) is not a valid type!'
            msg = msg.format(type(catl_kind))
            raise TypeError(msg)
        # `ext` - Value
        ext_arr = ['hdf5']
        if not (ext in ext_arr):
            msg = '`ext` ({0}) is not a valid input parameter!'
            msg = msg.format(ext)
            raise ValueError(msg)
        # `print_filedir` - Type
        if not (isinstance(print_filedir, bool)):
            msg = '`print_filedir` ({0}) is not a valid type!'
            msg = msg.format(type(print_filedir))
            raise TypeError(msg)
        # `return_len` - Type
        if not (isinstance(return_len, bool)):
            msg = '`return_len` ({0}) is not a valid type!'
            msg = msg.format(type(return_len))
            raise TypeError(msg)
        ##
        ## Path to catalogues
        catl_dirpath = self.catls_dir(  catl_type=catl_type,
                                        catl_kind=catl_kind,
                                        print_filedir=print_filedir)
        # List of files
        catl_arr = cfutils.Index(catl_dirpath, '.' + ext)
        # Returning number of files if necessary
        if return_len:
            return catl_arr, len(catl_arr)
        else:
            return catl_arr



















