#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018_12-12
# Last Modified: 2018_12-18
# Vanderbilt University
from __future__ import absolute_import, division, print_function 
__author__     = ['Victor Calderon']
__copyright__  = ["Copyright 2018 Victor Calderon, "]
__email__      = ['victor.calderon@vanderbilt.edu']
__maintainer__ = ['Victor Calderon']
__all__        = [  'DownloadManager']


import os
from cosmo_utils.utils import file_utils      as cfutils
from cosmo_utils.utils import work_paths      as cwpaths
from cosmo_utils.utils import web_utils       as cweb

from sdss_catl_utils.mocks_manager import mocks_defaults as md

class DownloadManager(object):
    """
    Class used to scrape the web for galaxy and group galaxy
    catalogue data and cache the downloaded catalogues.

    For list of available pre-processed galaxy- and group-galaxy
    catalogues proved by ``sdss_catl_utils``, see

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

    # Writes the location of catalogues as environment variable
    def _environ_variable_write(self, outdir):
        """
        Writes the location of catalogues as environment variable

        Parameters
        ------------
        outdir : `str`
            Path to the output directory, to which catalogues can be saved
            and stored.
        """
        # Bashrc and bash_profile location
        bashrc_path = os.path.expanduser('~/.bashrc')
        bash_profile_path = os.path.expanduser('~/.bash_profile')
        # Checking if files exist
        if os.path.exists(bashrc_path):
            conf_file = bashrc_path
        elif os.path.exists(bash_profile):
            conf_file = bash_profile
        else:
            conf_file = bashrc_file
            open(conf_file, 'a').close()
            print('`bashrc` file created: `{0}`'.format(bashrc_path))
            cfutils.File_Exists(conf_file)
        ##
        ## Checking that ``outdir`` exists
        assert(os.path.exists(outdir))
        # Appending to ``conf_file`` file if necessary
        if not os.environ.get(self.environ_name):
            with open(conf_file, 'a') as outfile:
                # String for exporting environment variable
                export_str = 'export {0}={1}'.format(   self.environ_name,
                                                        outdir)
                # Appending to file
                outfile.write('\n\n')
                outfile.write('## Output directory for SDSS Catalogues\n')
                outfile.write(export_str)
                outfile.write('\n\n')
                # Printing out message
                msg = '>> The system variable `{0}` has been '
                msg += 'added to your `{1}` file!'
                msg = msg.format(self.environ_name, conf_file)
                print(msg)

    # Checks file and folder structures
    def directory_tree_check(self, loc='./'):
        """
        Checks the file structure of the current directory, and determines
        if it is a ``git`` repository or not. If not, it will save the
        catalogues to a location in the user's home directory, under
        ``~/user/.sdss_catls``. It also adds the location of the
        output directory as an environment variable.

        Parameters
        ------------
        loc : `str`
            Path to the directory being analyzed for possible git
            initialization directory. This variable is set to ``./``
            by default.

        Returns
        ------------
        outdir : `str`
            Path to the output directory, to which catalogues can be saved
            and stored.
        """
        # Defining output directory
        try:
            # Environment variable
            catl_maindir = os.environ[self.environ_name]
            # Output directory - Check if directory exists
            outdir = os.path.join(  catl_maindir,
                                    'data',
                                    'external',
                                    'SDSS/')
            cfutils.Path_Folder(outdir)
            assert(os.path.exists(outdir))
        except KeyError:
            # Git Repository
            try:
                base_dir = cwpaths.git_root_dir(loc)
                outdir   = os.path.join(base_dir, 'data', 'external', 'SDSS/')
                mod_dir  = cwpaths.git_root_dir(path=os.path.realpath(__file__))
                # Creating directory if necessary
                if not (os.path.realpath(mod_dir) != os.path.realpath(base_dir)):
                    # Creating new directory
                    cfutils.Path_Folder(outdir)
                    assert(os.path.exists(outdir))
                else:
                    msg = '>> Current directory ({0}) is not supported! '
                    msg += 'You must move directories!'
                    msg = msg.format(outdir)
            except:
                # Home directory location
                base_dir = os.path.join(os.path.expanduser('~'),
                                        '.sdss_catls')
                # Output directory
                outdir   = os.path.join(base_dir,
                                        'data',
                                        'external',
                                        'SDSS/')
                # Creating directory
                cfutils.Path_Folder(outdir)
                assert(os.path.exists(outdir))
            ##
            ## Adding directory as environment variable
            self._environ_variable_write(outdir)

        return outdir

    # Output directory for local copies of catalogues
    def local_outdir_path(self, loc='./', catl_type='memb', catl_kind='data',
        check_exist=False, create_dir=False, perf_opt=False):
        """
        Output directory for ``local`` copies of the catalogues.

        Parameters
        ------------
        loc : `str`
            Path to the directory being analyzed for possible git
            initialization directory.

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

        check_exist : `bool`, optional
            If `True`, it checks for whether or not the file exists.
            This variable is set to `False` by default.

        create_dir : `bool`, optional
            If `True`, it creates the directory if it does not exist.
            This variable is set to `False` by default.

        perf_opt : `bool`, optional
            If True, it chooses to analyze the ``perfect`` version of
            the synthetic galaxy/group galaxy catalogues. Otherwise,
            it downloads the catalogues with group-finding errors
            included. This variable is set to ``False`` by default.

        Returns
        ----------
        outdir_local : `str`
            Path to the local output directory, based on the type of
            catalogues and input parameters.
        """
        ## Checking input parameters
        # `loc` - Type
        if not (isinstance(loc, str)):
            msg = '`loc` ({0}) is not a valid input type!'
            msg = msg.format(loc)
            raise TypeError(loc)
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
        # `perf_opt` - Type
        if not (isinstance(perf_opt, bool)):
            msg = '`perf_opt` ({0}) is not a valid type!'
            msg = msg.format(type(perf_opt))
            raise TypeError(msg)
        # `check_exist` - Type
        if not (isinstance(check_exist, bool)):
            msg = '`check_exist` ({0}) is not a valid input type!'
            msg = msg.format(type(check_exist))
            raise TypeError(msg)
        # `create_dir` - Type
        if not (isinstance(create_dir, bool)):
            msg = '`create_dir` ({0}) is not a valid input type!'
            msg = msg.format(type(create_dir))
            raise TypeError(msg)
        ##
        ## Local base path
        local_base   = self.directory_tree_check(loc=loc)
        ## Local Output directory
        local_outdir = os.path.join(local_base,
                                    self._catl_prefix(
                                        catl_kind=catl_kind,
                                        catl_type=catl_type,
                                        perf_opt=perf_opt))
        # Creating directory
        if create_dir:
            cfutils.Path_Folder(catl_output_dir)
        # Check for its existence
        if check_exist:
            if not (os.path.exists(catl_output_dir)):
                msg = '`catl_output_dir` ({0}) was not found!'.format(
                    catl_output_dir)
                raise FileNotFoundError(msg)

        return local_outdir

    # Prefix of the catalogues being downloaded
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

    # Checks that the combination of parameters is valid
    def catl_web_params_check(self, catl_kind='mocks', catl_type='memb',
        ext='hdf5', return_files_url=False, perf_opt=False):
        """
        Checks that the combination of parameters is valid and that there
        are catalogues available to download.

        Parameters
        ------------
        catl_kind : {``data``, ``mocks``} `str`
            Kind of catalogues to download. This variable is set to
            ``mocks`` by default.

            Options:
                - ``data``: Downloads the SDSS DR7 real catalogues.
                - ``mocks``: Downloads the synthetic catalogues of SDSS DR7.

        catl_type : {'gal', 'memb', 'group'}, `str`
            Option for which kind of catalogue is being analyzed. This
            variable is set to ``memb`` by default.

            Options:
                - ``'gal'`` : Galaxy catalogue
                - ``'memb'`` : Group Member galaxy catalogue
                - ``'group'`` : Group galaxy catalogue

        ext : {'hdf5'} `str`
            File extension used for the catalogues. This variable is set
            to ``hdf5`` by default.

        return_files_url : `bool`, optional
            If True, it returns an array of files from the desired combination
            of parameters, and with a file extension ``ext``. This variable
            is set to ``False`` by default.

        perf_opt : `bool`, optional
            If True, it chooses to analyze the ``perfect`` version of
            the synthetic galaxy/group galaxy catalogues. Otherwise,
            it downloads the catalogues with group-finding errors
            included. This variable is set to ``False`` by default.

        Returns
        ------------
        catl_url : `str`
            URL of the files being downloaded.

        catl_files_url : `numpy.ndarray`, optional
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
        # ``return_files_url``
        if not (isinstance(return_files_url, bool)):
            msg = '>>> `return_files_url` ({0}) is not a valid input type!'
            msg = msg.format(type(return_files_url))
            raise TypeError(msg)
        # `perf_opt` - Type
        if not (isinstance(perf_opt, bool)):
            msg = '`perf_opt` ({0}) is not a valid type!'
            msg = msg.format(type(perf_opt))
            raise TypeError(msg)
        ## Parsing URL
        catl_url = os.path.join(    md.sdss_catl_url,
                                    self._catl_prefix(
                                        catl_kind=catl_kind,
                                        catl_type=catl_type,
                                        perf_opt=perf_opt))
        ## Check that URL exists
        cweb.url_checker(catl_url)
        # Returning files
        if not return_files_url:
            return catl_url
        else:
            # Checking that there are files available
            catl_files_url = cweb.url_file_list(catl_url, ext)
            # Returning elements
            return catl_url, catl_files_url

    # Downloads the catalogues to the corresponding Output directory
    def catl_download(self, outdir='./', download_type='all', ext='hdf5'):
        """
        Downloads the corresponding catalogues to the designated folder.

        Parameters
        ------------
        outdir : `str`, optional
            Output directory, to which catalogues will be saved. This
            variable is set to the default location of
            `~directory_tree_check`. This variable is set to ``./``
            by default.

        download_type : {``all``, ``data``, ``mocks``} `str`, optional
            Type of catalogues to download. This variable is set to ``all``
            by default.

            Options:
                - ``'data'``: Downloads the SDSS `data` catalogues.
                - ``'mocks'``: Downloads the SDSS synthetic catalogues.
                - ``'all'``: Downloads both ``data`` and ``mocks`` catalogues.

        ext : {'hdf5'} `str`
            File extension used for the catalogues. This variable is set
            to ``hdf5`` by default.
        
        Examples
        ----------

        >>> # To download the synthetic catalogues # doctest: +SKIP
        >>> from sc.mocks_manager import DownloadManager
        >>>
        >>> # Downloading catalogues
        >>> A = DownloadManager()
        >>> A.catl_download(outdir='./', download_type='all') # doctest: +SKIP
        >>>
        """
        ## Checking input parameters
        # `download_type` - Value
        download_type_arr = ['all', 'data', 'mocks']
        if not (download_type in download_type_arr):
            msg = '`download_type` ({0}) is not a valid input parameter!'
            msg = msg.format(download_type)
            raise ValueError(msg)
        # `download_type` - Type
        if not (isinstance(download_type, str)):
            msg = '`download_type` ({0}) is not a valid input type!'
            msg = msg.format(type(download_type))
            raise TypeError(msg)
        # `outdir` - Type
        if not (isinstance(outdir, str)):
            msg = '`outdir` ({0}) is not a valid input type!'
            msg = msg.format(type(outdir))
            raise TypeError(msg)
        # `ext` - Type
        if not (isinstance(ext, str)):
            msg = '`ext` ({0}) is not a valid input type!'
            msg = msg.format(type(ext))
            raise TypeError(msg)
        # `ext` - Value
        ext_arr = ['hdf5']
        if not (ext in ext_arr):
            msg = '`ext` ({0}) is not a valid input parameter!'
            msg = msg.format(ext)
            raise ValueError(msg)
        # `check_exist` - Type
        if not (isinstance(check_exist, bool)):
            msg = '`check_exist` ({0}) is not a valid input type!'
            msg = msg.format(type(check_exist))
            raise TypeError(msg)
        # `create_dir` - Type
        if not (isinstance(create_dir, bool)):
            msg = '`create_dir` ({0}) is not a valid input type!'
            msg = msg.format(type(create_dir))
            raise TypeError(msg)
        ##
        ## Downloading catalogues
        # `download_type` Array
        if (download_type == 'mocks'):
            download_arr = ['mocks']
        elif (download_type == 'data'):
            download_arr = ['data']
        elif (download_type == 'all'):
            download_arr = ['data', 'mocks']
        # `catl_type` array
        catl_type_arr = ['memb', 'group']
        # Looping over different options
        for download_opt_ii in download_arr:
            for catl_type_ii in catl_type_arr:
                # Local Path
                local_ii = self.local_outdir_path(  loc=outdir,
                                                    catl_type=catl_type_ii,
                                                    catl_kind=download_opt_ii)
                # Web URL
                url_ii   = self.catl_web_params_check(
                                                    catl_kind=download_opt_ii,
                                                    catl_type=catl_type_ii,
                                                    ext=ext,
                                                    return_files_url=False)
                ##
                ## Downloading catalogues to corresponding catalogues
                cweb.url_files_download(url_ii_perf,
                                        ext,
                                        local_ii_perf,
                                        create_dir=True,
                                        check_exist=True,
                                        remove_files=self.remove_files)
                ##
                ## Perfect option
                if (self.perf_opt) and (download_opt_ii == 'mocks'):
                # Local Path
                    local_ii_perf = self.local_outdir_path(
                                                    loc=outdir,
                                                    catl_type=catl_type_ii,
                                                    catl_kind=download_opt_ii,
                                                    perf_opt=self.perf_opt)
                    # Web URL
                    url_ii_perf   = self.catl_web_params_check(
                                                    catl_kind=download_opt_ii,
                                                    catl_type=catl_type_ii,
                                                    ext=ext,
                                                    return_files_url=False,
                                                    perf_opt=self.perf_opt)
                    ##
                    ## Downloading catalogues to corresponding catalogues
                    cweb.url_files_download(url_ii_perf,
                                            ext,
                                            local_ii_perf,
                                            create_dir=True,
                                            check_exist=True,
                                            remove_files=self.remove_files)











    



























