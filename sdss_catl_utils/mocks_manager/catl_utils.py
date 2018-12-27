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
__all__        = [  'catl_keys',
                    'catl_keys_prop',
                    'catl_clean',
                    'catl_clean_nmin',
                    'catl_prefix_main',
                    'check_input_params',
                    'CatlUtils']

import os
import numpy as np
import pandas as pd
from collections import Counter
from cosmo_utils.utils import file_utils      as cfutils
from cosmo_utils.utils import work_paths      as cwpaths
from cosmo_utils.utils import web_utils       as cweb
from cosmo_utils.utils import file_readers    as cfreaders

from sdss_catl_utils.mocks_manager import mocks_defaults as md
from sdss_catl_utils.custom_exceptions import SDSSCatlUtils_Error

## -- Functions and classes -- ##

# Catalogue keys - Main
def catl_keys(catl_kind='data', perf_opt=False, return_type='list'):
    """
    Provides a dictionary/list with the corresponding keys for 1) halo mass,
    2) Haloid/GroupID, 3) Group/Halo Galaxy-Type.

    Parameters
    ------------
    catl_kind : {``data``, ``mocks``} `str`, optional
        Type of the catalogue being analyzed. This variable corresponds
        to whether a ``real`` or ``synthetic/mock`` catalogue is being
        read/analyzed. This variable is set to ``data`` by default.

        Options:
            - ``data``: Catalogue(s) from the SDSS `real` catalogues
            - ``mocks``: Catalogue(s) from the `mock` catalogues.

    perf_opt : `bool`, optional
        If `True`, it returns the corresponding keys for a ``perfect``
        SDSS catalogue. This option only applies when ``catl_kind == 'mocks'``.

    return_type : {``list``, ``dict``} `str`, optional
        Type of output to be returned. This variable is set to `list` by
        default.

        Options:
            - ``list``: Returns the output as part of a list. The order of
                        the elements are: ``'group mass'``,
                        ``'groupid/haloid'``, and ``Group/Halo ID``.

    Returns
    ------------
    catl_objs : `dict` or `list`
        Dictionary or list of keys for ``group/halo mass``,
        ``group ID`` and ``galaxy type``
        columns in catalogues. This variable depends on the choice of
        `return_type`.

    Raises
    ------------
    SDSSCatlUtils_Error : Exception from `~sdss_catl_utils.SDSSCatlUtils_Error`
        Program exception if input parameters are `not` accepted.

    Examples
    ------------
    This function can be used for several different combinations of parameters.
    For example, if we wanted to analyze a ``data`` catalogue, and have this
    function return a ``list`` as the output, one could write::

        >>> catl_keys(catl_kind='data', return_type='list')
        ['M_h', 'groupid', 'galtype']

    This list corresponds to the 1) Group estimated mass, 2) Galaxy's group ID,
    and 3) Galaxy's group galaxy type.

    If instead, we wanted to analyze ``perfect mock catalogues``, we
    could write::

        >>> catl_keys(catl_kind='mocks', perf_opt=True, return_type='list')
        ['M_h', 'haloid', 'galtype']

    For more information and examples, please refer to
    :ref:`quickstart_getting_started`.
    """
    file_msg = cfutils.Program_Msg(__file__)
    ## Checking input parameters
    # `catl_kind` - Value
    catl_kind_arr = ['data', 'mocks']
    if not (catl_kind in catl_kind_arr):
        msg = '{0} `catl_kind` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, catl_kind)
        raise SDSSCatlUtils_Error(msg)
    # `catl_kind` - Type
    if not (isinstance(catl_kind, str)):
        msg = '{0} `catl_kind` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(catl_kind))
        raise TypeError(msg)
    # `perf_opt`  - Type
    if not (isinstance(perf_opt, bool)):
        msg = '{0} `perf_opt` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(perf_opt))
        raise TypeError(msg)
    # `return_type` - Value
    return_type_arr = ['list', 'dict']
    if not (return_type in return_type_arr):
        msg = '{0} `return_type` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, return_type)
        raise SDSSCatlUtils_Error(msg)
    # `return_type` - Type
    if not (isinstance(return_type, str)):
        msg = '{0} `return_type` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(return_type))
        raise TypeError(msg)
    ##
    ## If analyzing real catalogues, setting ``perf_opt`` to False.
    if (catl_kind == 'data'):
        perf_opt = False
    ##
    ## Keys for the different sets of catalogues
    if (catl_kind == 'data'):
        (   gm_key,
            id_key,
            galtype) = ['M_h', 'groupid', 'galtype']
    elif (catl_kind == 'mocks'):
        if perf_opt:
            (   gm_key,
                id_key,
                galtype) = ['M_h', 'haloid', 'galtype']
        else:
            (   gm_key,
                id_key,
                galtype) = ['M_group', 'groupid', 'g_galtype']
    ##
    ## Determining which type of output to return
    if (return_type == 'list'):
        catl_objs = [gm_key, id_key, galtype]
    elif (return_type == 'dict'):
        catl_objs = {'gm_key': gm_key, 'id_key': id_key,
                        'galtype_key':galtype}

    return catl_objs

# Catalogue Keys - Galaxy properties
def catl_keys_prop(catl_kind='data', catl_info='memb', return_type='list'):
    """
    Provides a dictionary/list with the corresponding keys for 1) specific
    star formation rate (sSFR) and 2) stellar mass.

    Parameters
    ------------
    catl_kind : {``data``, ``mocks``} `str`, optional
        Type of the catalogue being analyzed. This variable corresponds
        to whether a ``real`` or ``synthetic/mock`` catalogue is being
        read/analyzed. This variable is set to ``data`` by default.

        Options:
            - ``data``: Catalogue(s) from the SDSS `real` catalogues
            - ``mocks``: Catalogue(s) from the `mock` catalogues.

    catl_info : {``memb``, ``groups``} `bool`, optional
        Option for which type of catalogue is being analyzed. This variable
        correspondos to whether a ``galaxy``-catalogue or a ``group``-catalogue
        is being analyzed. This variable is set to ``memb`` by default.

        Options:
            - ``memb``: Galaxy catalogue with the `member` galaxies of groups.
            - ``groups``: Catalogues with `group` information.

    return_type : {``list``, ``dict``} `str`, optional
        Type of output to be returned. This variable is set to `list` by
        default.

        Options:
            - ``list``: Returns the output as part of a list. The order of
                        the elements are: ``'group mass'``,
                        ``'groupid/haloid'``, and ``Group/Halo ID``.
    Returns
    ------------
    catl_objs : `dict` or `list`
        Dictionary or list of keys for ``logssfr`` and ``logmstar``
        columns in catalogues. This variable depends on the choice of
        `return_type`.

    Raises
    ------------
    SDSSCatlUtils_Error : Exception from `~sdss_catl_utils.SDSSCatlUtils_Error`
        Program exception if input parameters are `not` accepted.

    Examples
    ------------
    This function can be used for several different combinations of parameters.
    For example, if we wanted to analyze a ``data`` catalogue, and have this
    function return a ``list`` as the output, one could write::

        >>> catl_keys_prop(catl_kind='data', return_type='list')
        ['logssfr', 'logMstar_JHU']

    This list corresponds to the 1) `specific star formation rate` key and
    2) `stellar mass` key.

    If instead, we wanted to analyze a ``mock`` catalogue and return the
    appropriate keys for a ``group catalogue``, we could write::

        >>> catl_keys_prop(catl_kind='mocks', catl_info='groups')
        ['logssfr', 'logMstar']

    For more information and examples, please refer to
    :ref:`quickstart_getting_started`.
    """
    file_msg = cfutils.Program_Msg(__file__)
    ## Checking input parameters
    # `catl_kind` - Value
    catl_kind_arr = ['data', 'mocks']
    if not (catl_kind in catl_kind_arr):
        msg = '{0} `catl_kind` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, catl_kind)
        raise SDSSCatlUtils_Error(msg)
    # `catl_kind` - Type
    if not (isinstance(catl_kind, str)):
        msg = '{0} `catl_kind` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(catl_kind))
        raise TypeError(msg)
    # `catl_info` - Value
    catl_info_arr = ['memb', 'groups']
    if not (catl_info in catl_info_arr):
        msg = '{0} `catl_info` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, catl_info)
        raise SDSSCatlUtils_Error(msg)
    # `catl_info` - Type
    if not (isinstance(catl_info, str)):
        msg = '{0} `catl_info` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(catl_info))
        raise TypeError(msg)
    # `return_type` - Value
    return_type_arr = ['list', 'dict']
    if not (return_type in return_type_arr):
        msg = '{0} `return_type` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, return_type)
        raise SDSSCatlUtils_Error(msg)
    # `return_type` - Type
    if not (isinstance(return_type, str)):
        msg = '{0} `return_type` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(return_type))
        raise TypeError(msg)
    ##
    ## Dictionary with entries for ssfr and mstar for different combinations.
    keys_dict = {   'data_memb'   : ['logssfr', 'logMstar_JHU'],
                    'data_groups' : ['logssfr_tot', 'logMstar_tot'],
                    'mocks_memb'  : ['logssfr', 'logMstar'],
                    'mocks_groups': ['logssfr', 'logMstar']}
    # Deciding which key to use
    catl_key_str = '{0}_{1}'.format(catl_kind, catl_info)
    #
    # Determining which type of output to return
    if (return_type == 'list'):
        catl_objs = keys_dict[catl_key_str]
    elif (return_type == 'dict'):
        catl_objs = {   'logssfr_key' : keys_dict[catl_key_str][0],
                        'logmstar_key': keys_dict[catl_key_str][1]}

    return catl_objs

# Cleaning the catalogues from `bad` inputs
def catl_clean(catl_pd, catl_kind, catl_info='memb', reindex=True):
    """
    Cleans and removes the ``bad`` rows, i.e. those that contain `failed`
    entries for sSFR and Mstar.

    Parameters
    ------------
    catl_pd : `pandas.DataFrame`
        DataFrame containing the information about galaxies or galaxy groups.

    catl_kind : {``data``, ``mocks``} `str`, optional
        Type of the catalogue being analyzed. This variable corresponds
        to whether a ``real`` or ``synthetic/mock`` catalogue is being
        read/analyzed.

        Options:
            - ``data``: Catalogue(s) from the SDSS `real` catalogues
            - ``mocks``: Catalogue(s) from the `mock` catalogues.

    catl_info : {``memb``, ``groups``} `bool`, optional
        Option for which type of catalogue is being analyzed. This variable
        correspondos to whether a ``galaxy``-catalogue or a ``group``-catalogue
        is being analyzed. This variable is set to ``memb`` by default.

        Options:
            - ``memb``: Galaxy catalogue with the `member` galaxies of groups.
            - ``groups``: Catalogues with `group` information.
    
    reindex : `bool`, optional
        If `True`, the output catalogue is reindexed from the original dataframe
        `catl_pd`. This variable is set to `True` by default.

    Returns
    ----------
    catl_pd_mod : `pandas.DataFrame`
        Modified clean version of `catl_pd`. It removes the `failed`
        values.

    Raises
    ------------
    SDSSCatlUtils_Error : Exception from `~sdss_catl_utils.SDSSCatlUtils_Error`
        Program exception if input parameters are `not` accepted.
    """
    file_msg = cfutils.Program_Msg(__file__)
    ## Checking input parameters
    # `catl_pd` - Type
    if not (isinstance(catl_pd, pd.DataFrame)):
        msg = '{0} `catl_pd` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(catl_pd))
        raise TypeError(msg)
    # `catl_kind` - Value
    catl_kind_arr = ['data', 'mocks']
    if not (catl_kind in catl_kind_arr):
        msg = '{0} `catl_kind` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, catl_kind)
        raise SDSSCatlUtils_Error(msg)
    # `catl_kind` - Type
    if not (isinstance(catl_kind, str)):
        msg = '{0} `catl_kind` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(catl_kind))
        raise TypeError(msg)
    # `catl_info` - Value
    catl_info_arr = ['memb', 'groups']
    if not (catl_info in catl_info_arr):
        msg = '{0} `catl_info` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, catl_info)
        raise SDSSCatlUtils_Error(msg)
    # `catl_info` - Type
    if not (isinstance(catl_info, str)):
        msg = '{0} `catl_info` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(catl_info))
        raise TypeError(msg)
    # `reindex  ` - Type
    if not (isinstance(reindex, bool)):
        msg = '{0} `reindex` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(reindex))
        raise TypeError(msg)
    ##
    ## List of values that are considered `failed`. These values
    ## were assigned by the MPA-JHU group that developed the VAGC catalogue.
    ssfr_fail_arr = [0, -99, -999, np.nan]
    mstar_fail_arr= [-1, 0, np.nan]
    # Label for `ssfr` and `mstar`
    logssfr_key, logmstar_key = catl_keys_prop( catl_kind=catl_kind,
                                                catl_info=catl_info,
                                                return_type='list')
    ##
    ## Cleaning catalogue entries
    if (logmstar_key in catl_pd.columns.values):
        catl_pd_mod = catl_pd.loc[~catl_pd[logssfr_key].isin(ssfr_fail_arr) &
                                  ~catl_pd[logmstar_key].isin(mstar_fail_arr)]
    else:
        catl_pd_mod = catl_pd.loc[~catl_pd[logssfr_key].isin(ssfr_fail_arr)]
    #
    # Option if to reindex the new `clean` catalogue
    if reindex:
        catl_pd_mod.reset_index(drop=True, inplace=True)

    return catl_pd_mod

# Cleaning the catalogue from `fail` values and only include galaxies
# from groups larger than the number of galaxy threshold, i.e. ``nmin``.
def catl_clean_nmin(catl_pd, catl_kind, catl_info='memb', reindex=True,
    nmin=1, perf_opt=False):
    """
    Cleans and removed the ``bad`` rows with `failed` values, i.e.
    those that contain `failed` entries. This method also includes
    galaxies from groups above the ``nmin`` galaxy number threshold.

    Parameters
    --------------
    catl_pd : `pandas.DataFrame`
        DataFrame containing the information about galaxies or galaxy groups.

    catl_kind : {``data``, ``mocks``} `str`, optional
        Type of the catalogue being analyzed. This variable corresponds
        to whether a ``real`` or ``synthetic/mock`` catalogue is being
        read/analyzed.

        Options:
            - ``data``: Catalogue(s) from the SDSS `real` catalogues
            - ``mocks``: Catalogue(s) from the `mock` catalogues.

    catl_info : {``memb``, ``groups``} `bool`, optional
        Option for which type of catalogue is being analyzed. This variable
        correspondos to whether a ``galaxy``-catalogue or a ``group``-catalogue
        is being analyzed. This variable is set to ``memb`` by default.

        Options:
            - ``memb``: Galaxy catalogue with the `member` galaxies of groups.
            - ``groups``: Catalogues with `group` information.
    
    reindex : `bool`, optional
        If `True`, the output catalogue is reindexed from the original dataframe
        `catl_pd`. This variable is set to `True` by default.

    nmin : `int`, optional
        Minimum group richness to have in the (galaxy) group catalogue.
        This variable is set to ``1`` by default, and must be larger than
        `1`.

    perf_opt : `bool`, optional
        Option for using a `perfect` mock catalogue. This variable is set
        to `False` by default.

    Returns
    --------------
    catl_pd_mod : `pandas.DataFrame`
        Version of `catl_pd` after having removed the `failed` values
        of `sSFR` and `Mstar`, and also after having chosen only galaxies
        and groups with group richnesses larger than ``nmin``.

    Raises
    ------------
    SDSSCatlUtils_Error : Exception from `~sdss_catl_utils.SDSSCatlUtils_Error`
        Program exception if input parameters are `not` accepted.

    Examples
    ------------
    Before using this function, one needs to have read one of the
    (galaxy) group catalogues. If for example, one wants to create a
    new object from the ``data`` `real` SDSS catalogue with galaxies 
    from groups with ``n > 10``, one can do:

    >>> from cosmo_utils.utils import file_readers as cfr
    >>> from sdss_catl_utils.mocks_manager.catl_utils import catl_clean_nmin
    >>> nmin = 10 # Minimum number of galaxies in file
    >>> catl_pd = cfr.read_hdf5_file_to_pandas_DF('/path/to/file') # doctest: +SKIP
    >>> catl_mod = catl_clean_nmin(catl_pd, 'data', nmin=nmin) # doctest: +SKIP

    Now, the resulting catalogue will only include galaxies from groups
    with ``n > 10``.

    """
    file_msg = cfutils.Program_Msg(__file__)
    ## Checking input parameters
    # `catl_pd` - Type
    if not (isinstance(catl_pd, pd.DataFrame)):
        msg = '{0} `catl_pd` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(catl_pd))
        raise TypeError(msg)
    # `catl_kind` - Value
    catl_kind_arr = ['data', 'mocks']
    if not (catl_kind in catl_kind_arr):
        msg = '{0} `catl_kind` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, catl_kind)
        raise SDSSCatlUtils_Error(msg)
    # `catl_kind` - Type
    if not (isinstance(catl_kind, str)):
        msg = '{0} `catl_kind` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(catl_kind))
        raise TypeError(msg)
    # `catl_info` - Value
    catl_info_arr = ['memb', 'groups']
    if not (catl_info in catl_info_arr):
        msg = '{0} `catl_info` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, catl_info)
        raise SDSSCatlUtils_Error(msg)
    # `catl_info` - Type
    if not (isinstance(catl_info, str)):
        msg = '{0} `catl_info` ({1}) is not a valid input type!!'
        msg = msg.format(file_msg, type(catl_info))
        raise TypeError(msg)
    # `reindex  ` - Type
    if not (isinstance(reindex, bool)):
        msg = '{0} `reindex` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(reindex))
        raise TypeError(msg)
    # `nmin` - Type
    nmin_type_arr = (int, float)
    if (not isinstance(nmin, nmin_type_arr)):
        msg = '{0} `nmin` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(nmin))
        raise TypeError(msg)
    # `nmin` - Value
    if not (nmin >= 1):
        msg = '{0} `nmin` ({1}) must be larger than `1`!'
        msg = msg.format(file_m, nmin)
        raise SDSSCatlUtils_Error(msg)
    # `perf_opt` - Type
    if not (isinstance(perf_opt, bool)):
        msg = '{0} `perf_opt` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(perf_opt))
        raise TypeError(msg)
    ##
    ## Types of galaxies
    cens = int(1)
    nmin = int(nmin)
    # Keys for the catalogue
    gm_key, id_key, galtype_key = catl_keys(catl_kind, perf_opt=perf_opt)
    # Cleaning catalogue entries
    catl_pd_mod = catl_clean(   catl_pd,
                                catl_kind=catl_kind,
                                catl_info=catl_info,
                                reindex=reindex)
    # Choosing only galaxies in groups of richness >= `nmin`
    if (catl_info == 'groups'):
        if ('ngals' in catl_pd_mod.columns):
            catl_pd_mod_nmin = catl_pd_mod.loc[catl_pd_mod['ngals'] >= nmin]
        else:
            msg = '{0} Key `ngals` not found in DataFrame!'.format(file_msg)
            raise SDSSCatlUtils_Error(msg)
    elif (catl_info == 'memb'):
        # List of central galaxies
        cens_pd = catl_pd_mod.loc[(catl_pd_mod[galtype_key] == cens), id_key]
        catl_pd_cen = catl_pd_mod.loc[catl_pd_mod[id_key].isin(cens_pd)]
        # Counting the number of galaxies in each galaxy group and choosing
        # those above the `nmin` threshold.
        g_counts = Counter(catl_pd_cen[id_key].values)
        g_ngals  = [xx for xx in g_counts.keys() if g_counts[xx]>=nmin]
        # Selecting only galaxies with ``group_id`` in ``g_ngals``.
        catl_pd_mod_nmin = catl_pd_cen.loc[catl_pd_cen[id_key].isin(g_ngals)]
    #
    # Resetting index if necessary
    if reindex:
        catl_pd_mod_nmin.reset_index(drop=True, inplace=True)

    return catl_pd_mod_nmin

# Catalogue prefix of the catalogues
def catl_prefix_main(catl_type='memb', catl_kind=md.catl_kind, hod_n=md.hod_n,
    halotype=md.halotype, clf_method=md.clf_method, clf_seed=md.clf_seed,
    dv=md.dv, sample=md.sample, type_am=md.type_am, perf_opt=md.perf_opt):
    """
    Prefix of the paths based on the type of catalogues and input parameters
    chosen.

    Parameters
    -----------
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

    type_am : {'mr', 'mstar'}, `str`, optional
        Type of Abundance matching used in the catalogue. This
        variable is set to ``'mr'`` by default.

        Options:
            - ``'mr'``: Luminosity-based abundance matching used
            - ``'mstar'``: Stellar-mass-based abundance matching used.

    perf_opt : `bool`, optional
        If `True`, it chooses to analyze the ``perfect`` version of
        the synthetic galaxy/group galaxy catalogues. Otherwise,
        it downloads the catalogues with group-finding errors
        included. This variable is set to ``False`` by default.

    Returns
    ---------
    catl_prefix : `str`
        Prefix of the paths based on the type of catalogues and input
        parameters.
    """
    file_msg = cfutils.Program_Msg(__file__)
    ## Checking input parameters
    # `catl_type`
    check_input_params(catl_type, 'catl_type', check_type='type')
    check_input_params(catl_type, 'catl_type', check_type='vals')
    # `catl_kind`
    check_input_params(catl_kind, 'catl_kind', check_type='type')
    check_input_params(catl_kind, 'catl_kind', check_type='vals')
    # `hod_n`
    check_input_params(hod_n, 'hod_n', check_type='type')
    check_input_params(hod_n, 'hod_n', check_type='vals')
    # `halotype`
    check_input_params(halotype, 'halotype', check_type='type')
    check_input_params(halotype, 'halotype', check_type='vals')
    # `clf_method`
    check_input_params(clf_method, 'clf_method', check_type='type')
    check_input_params(clf_method, 'clf_method', check_type='vals')
    # `clf_seed`
    check_input_params(clf_seed, 'clf_seed', check_type='type')
    # `dv`
    check_input_params(dv, 'dv', check_type='type')
    # `sample`
    check_input_params(sample, 'sample', check_type='type')
    check_input_params(sample, 'sample', check_type='vals')
    # `type_am`
    check_input_params(type_am, 'type_am', check_type='type')
    check_input_params(type_am, 'type_am', check_type='vals')
    # `perf_opt`
    check_input_params(perf_opt, 'perf_opt', check_type='type')
    ##
    ## Option for which type of catalogue to download
    catl_type_dict = {  'gal'  : 'galaxy_catalogues',
                        'memb' : 'member_galaxy_catalogues',
                        'group': 'group_galaxy_catalogues',
                        'perf_group': 'perfect_group_galaxy_catalogues',
                        'perf_memb' : 'perfect_member_galaxy_catalogues'}
    # Setting `perf_opt` to `False` if necessary
    if (catl_kind == 'data'):
        perf_opt = False
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
    # Extra parameters
    sample_Mr = 'Mr{0}'.format(sample)
    ##
    ## Parsing prefix path
    # `Data`
    if (catl_kind == 'data'):
        catl_prefix = os.path.join( 'data',
                                    type_am,
                                    sample_Mr,
                                    catl_type_dict[catl_type_str])
    # `Mocks`
    if (catl_kind == 'mocks'):
        catl_prefix = os.path.join(
                                'mocks',
                                'halos_{0}'.format(halotype),
                                'dv_{0}'.format(dv),
                                'hod_model_{0}'.format(hod_n),
                                'clf_seed_{0}'.format(clf_seed),
                                'clf_method_{0}'.format(clf_method),
                                type_am,
                                sample_Mr,
                                catl_type_dict[catl_type_str])

    return catl_prefix

# Directory with accepted input parameters
def _get_input_params_dict():
    """
    Checks the input parameters to the class object.

    Returns
    ---------
    input_dict : `dict`
        Dictionary with the accepted ``type`` and ``value`` input
        parameters.
    """
    ## Dictionary with input variables
    # Variable types
    input_dict_type = { 'catl_kind'    : (str),
                        'hod_n'        : (int),
                        'halotype'     : (str),
                        'clf_method'   : (int),
                        'clf_seed'     : (int),
                        'dv'           : (int, float),
                        'sample'       : (str),
                        'type_am'      : (str),
                        'cosmo_choice' : (str),
                        'perf_opt'     : (bool),
                        'remove_files' : (bool),
                        'environ_name' : (str),
                        'catl_type'    : (str)}
    # Variable inputs
    input_dict_vals = { 'catl_kind'    : ['data', 'mocks'],
                        'hod_n'        : list(range(10)),
                        'halotype'     : ['fof', 'so'],
                        'clf_method'   : [1, 2, 3],
                        'sample'       : ['19', '20', '21'],
                        'type_am'      : ['mr', 'mstar'],
                        'cosmo_choice' : ['LasDamas', 'Planck'],
                        'catl_type'    : ['gal', 'memb', 'group']}
    # Merging dictionaries
    input_dict = {'type': input_dict_type, 'vals': input_dict_vals}

    return input_dict

# Checking input parameters to make sure they are `expected`
def check_input_params(input_var, var_name, check_type='type'):
    """
    Checks the type and/or values for different variables.

    Parameters
    ------------
    input_var : `int`, `float`, `bool`, `str`
        Input variable to be evaluated.

    var_name : `str`
        Name of the input parameter being evaluated. This variable name
        must correspond to one of the keys in the `type` or `vals`
        dictionaries.

    check_type : {``type``, ``vals``} `str`
        Type of check to perform. This variable is set to ``type`` by default.

        Options:
            - ``type``: It checks for the `type` of `input_var`.
            - ``vals``: It checks for the `value` of `input_var`.
    
    Raises
    ------------
    SDSSCatlUtils_Error : Exception from `~sdss_catl_utils.SDSSCatlUtils_Error`
        Program exception if input parameters are `not` accepted.
    """
    file_msg = cfutils.Program_Msg(__file__)
    ##
    ## Checking input parameters
    # `input_var` - Type
    input_var_arr = (int, float, bool, str)
    if not (isinstance(input_var, input_var_arr)):
        msg = '{0} `input_var` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(input_var))
        raise TypeError(msg)
    # `check_type` - Type
    if not (isinstance(check_type, str)):
        msg = '{0} `check_type` ({1}) is not a valid input type!'
        msg = msg.format(file_msg, type(check_type))
        raise TypeError(msg)
    # `check_type` - Value
    check_type_arr = ['type', 'vals']
    if not (check_type in check_type_arr):
        msg = '{0} `check_type` ({1}) is not a valid input value!'
        msg = msg.format(file_msg, check_type)
        raise ValueError(msg)
    ##
    ## Extracting input dictinaries
    input_dict = _get_input_params_dict()
    # Obtaining list of acceptable inputs
    if (check_type == 'type'):
        # Checking for input parameters - Type
        check_dict = input_dict[check_type]
        # Checking if key exists
        if (var_name in list(check_dict.keys())):
            if not (isinstance(input_var, check_dict[var_name])):
                msg = '{0} `{1}` ({2}) is not a valid input! Valid: `{3}` type'
                msg = msg.format(file_msg, var_name, type(input_var),
                        check_dict[var_name])
                raise TypeError(msg)
        else:
            # If `var_name` is not in dictionary
            msg = '{0} `{1}` is not in dictionary! Keys: `{2}`'
            msg = msg.format(file_msg, var_name, list(check_dict.keys()))
            raise KeyError(msg)
    elif (check_type == 'vals'):
        # Checking for input parameters - Values
        check_dict = input_dict[check_type]
        # Checking if key exists
        if (var_name in list(check_dict.keys())):
            if not (input_var in check_dict[var_name]):
                msg = '{0} `{1}` ({2}) is not a valid input! Valid: `{3}` type'
                msg = msg.format(file_msg, var_name, input_var,
                        check_dict[var_name])
                raise ValueError(msg)
        else:
            # If `var_name` is not in dictionary
            msg = '{0} `{1}` is not in dictionary! Keys: `{2}`'
            msg = msg.format(file_msg, var_name, list(check_dict.keys()))
            raise KeyError(msg)

# Main class to handle catalogues
class CatlUtils(object):
    """
    Class used to handle the galaxy/group galaxy/group catalogues, *after*
    the catalogues have been downloaded. This class has functions to read,
    modify, and analyze the galaxy/group catalogues.

    For a list of available pre-processed galaxy- and group-galaxy
    catalogues provided by ``sdss_catl_utils``, see ...
    """
    def __init__(self, **kwargs):
        r"""
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

        type_am : {'mr', 'mstar'}, `str`, optional
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
            If `True`, it chooses to analyze the ``perfect`` version of
            the synthetic galaxy/group galaxy catalogues. Otherwise,
            it downloads the catalogues with group-finding errors
            included. This variable is set to ``False`` by default.

        environ_name : `str`
            Name of the environment variable to assign to ``outdir``.
            This variable is set to the default ``environ_name`` from
            `~sdss_catl_utils.mocks_manager.mocks_default`

        Examples
        ----------

        >>> from sdss_catl_utils.mocks_manager.catl_utils import CatlUtils

        If one wants to initialize an ``catalogue`` object with the default
        input parameters, one could write:

        >>> catl_obj = CatlUtils() # Initialzing catalogue object

        However, if one wants to have a `modified` version of ``catl_obj``,
        one can pass a dictionary with the proper set of input parameters:

        >>> params_d = {'catl_kind': 'mocks', 'clf_method': 3, 'clf_seed': 123}
        >>> catl_obj = CatlUtils(**params_d)

        This will create an `catalogue object` that corresponds to ``mocks``
        catalogues with ``clf_method = 3`` and ``clf_seed = 123``.
        
        Notes
        ---------

        There are many combinations of parameters that could be performed.
        `However, note that not all variations exist and would not be
        available for download/analysis`!
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
        # Checking input parameters
        self._check_input_parameters()
        # Dictionary of input parameters
        self.param_dict = self.get_params_dict()

    # Checking input parameters to make sure they are `expected`
    def _check_input_parameters(self):
        r"""
        Checks whether or not the input parameters are what is expected or not.
        """
        ## Checking for input TYPES
        # `catl_kind`
        check_input_params(self.catl_kind, 'catl_kind', check_type='type')
        check_input_params(self.catl_kind, 'catl_kind', check_type='vals')
        # `hod_n`
        check_input_params(self.hod_n, 'hod_n', check_type='type')
        check_input_params(self.hod_n, 'hod_n', check_type='vals')
        # `halotype`
        check_input_params(self.halotype, 'halotype', check_type='type')
        check_input_params(self.halotype, 'halotype', check_type='vals')
        # `clf_method`
        check_input_params(self.clf_method, 'clf_method', check_type='type')
        check_input_params(self.clf_method, 'clf_method', check_type='vals')
        # `clf_seed`
        check_input_params(self.clf_seed, 'clf_seed', check_type='type')
        # `dv`
        check_input_params(self.dv, 'dv', check_type='type')
        # `sample`
        check_input_params(self.sample, 'sample', check_type='type')
        check_input_params(self.sample, 'sample', check_type='vals')
        # `type_am`
        check_input_params(self.type_am, 'type_am', check_type='type')
        check_input_params(self.type_am, 'type_am', check_type='vals')
        # `cosmo_choice`
        check_input_params(self.cosmo_choice, 'cosmo_choice', check_type='type')
        check_input_params(self.cosmo_choice, 'cosmo_choice', check_type='vals')
        # `perf_opt`
        check_input_params(self.perf_opt, 'perf_opt', check_type='type')
        # `remove_files`
        check_input_params(self.remove_files, 'remove_files', check_type='type')
        # `environ_name`
        check_input_params(self.environ_name, 'environ_name', check_type='type')

    # Main directory path - Path to which all catalogues are saved
    def main_dir(self):
        """
        Path to the main directory, in which all catalogues are stored.
        This directory depends on your installation of `SDSS_Catl_Utils`
        and on your current working space.

        Returns
        ---------
        main_dirpath : `str`
            Path to the main directory, in which all catalogues are stored.
            It uses the environment variable ``environ_name`` from
            the class.

        Examples
        ----------
        The main directory of the catalogues can be easily accessed after
        having created an object for the combination of input parameters
        for the catalogues.
        
        >>> from sdss_catl_utils.mocks_manager.catl_utils import CatlUtils
        >>> catl_params_dict = {'catl_kind': 'mocks', 'clf_seed': 3} # Catalogue parameters
        >>> catl_obj = CatlUtils(**catl_params_dict) # Initialing object

        One can access the directory, under which the galaxy-
        and group-catalogues are saved, by typing:

        >>> catl_obj.main_dir() # doctest: +SKIP

        This will return the path of the directory.

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

    # Get dictionary of input parameters
    def get_params_dict(self):
        r"""
        Gets the dictionary of input parameters for the given model.

        Returns
        ---------
        param_dict : `dict`
            Dictionary of input parameters for the chosen combination of
            parameters.

        Examples
        ----------
        After having initializing an `CatlUtils` object, one can easily
        recover the set of input parameters used.

        >>> from sdss_catl_utils.mocks_manager.catl_utils import CatlUtils
        >>> catl_params_dict = {'catl_kind': 'mocks', 'clf_seed': 3} # Catalogue parameters
        >>> catl_obj = CatlUtils(**catl_params_dict) # Initialing object

        The dictionary of parameters is saved as ``self.param_dict`` for
        the class object. It can be easily accessed by:

        >>> catl_obj.param_dict # doctest: +SKIP

        Or, it can also be accessed as:

        >>> catl_obj.get_params_dict() # doctest: +SKIP
        
        For example, in order to print out the list of input parameters,
        one can do the following:

        >>> param_dict = catl_obj.param_dict # doctest: +SKIP
        >>> param_dict # doctest: +SKIP
        {'catl_kind': 'mocks',
         'hod_n': 0,
         'halotype': 'fof',
         'clf_method': 1,
         'clf_seed': 3,
         'dv': 1.0,
         'sample': '19',
         'type_am': 'mr',
         'cosmo_choice': 'LasDamas',
         'perf_opt': False,
         'remove_files': False,
         'environ_name': 'sdss_catl_path',
         'sample_Mr': 'Mr19',
         'sample_s': '19'}

        This dictionary can now be used in other parts of one's analysis,
        or as a reference of the parameters used.
        """
        # Initializing dictionary
        param_dict = {}
        # Populating dictionary
        param_dict['catl_kind'   ] = self.catl_kind
        param_dict['hod_n'       ] = self.hod_n
        param_dict['halotype'    ] = self.halotype
        param_dict['clf_method'  ] = self.clf_method
        param_dict['clf_seed'    ] = self.clf_seed
        param_dict['dv'          ] = self.dv
        param_dict['sample'      ] = self.sample
        param_dict['type_am'     ] = self.type_am
        param_dict['cosmo_choice'] = self.cosmo_choice
        param_dict['perf_opt'    ] = self.perf_opt
        param_dict['remove_files'] = self.remove_files
        param_dict['environ_name'] = self.environ_name
        param_dict['sample_Mr'   ] = self.sample_Mr
        param_dict['sample_s'    ] = self.sample_s

        return param_dict

    # Location prefix
    def _catl_prefix(self, catl_type='memb', catl_kind='mocks',
        perf_opt=False):
        r"""
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
            If `True`, it chooses to analyze the ``perfect`` version of
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
        # `catl_type` - Type and Value
        check_input_params(catl_type, 'catl_type', check_type='type')
        check_input_params(catl_type, 'catl_type', check_type='vals')
        # `catl_kind` - Type and Value
        check_input_params(catl_kind, 'catl_kind', check_type='type')
        check_input_params(catl_kind, 'catl_kind', check_type='vals')
        # `perf_opt`  - Type
        check_input_params(perf_opt, 'perf_opt', check_type='type')
        ##
        ## Catalogue prefix
        catl_prefix = catl_prefix_main( catl_type=catl_type,
                                        catl_kind=catl_kind,
                                        perf_opt=perf_opt,
                                        hod_n=self.hod_n,
                                        halotype=self.halotype,
                                        clf_method=self.clf_method,
                                        clf_seed=self.clf_seed,
                                        dv=self.dv,
                                        sample=self.sample,
                                        type_am=self.type_am)

        return catl_prefix

    # Location of galaxy/group catalogues with specified parameters
    def catls_dir(self, catl_type='memb', catl_kind='mocks',
        print_filedir=False):
        r"""
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
            If `True`, the path of the catalogue directory is printed
            onto the screen. This variable is set to ``False`` by default.

        Returns
        ---------
        catls_dirpath : `str`
            Path to the location of the group and galaxy catalogues with
            the specified parameters.

        Notes
        ---------
        This method can be used to access the different `kinds` of
        catalogues, i.e. galaxy-, member-, and group-catalogues for different
        combinations of input parameters.

        Examples
        ----------
        `catls_dir` can be used to explore the directories that host different
        `kinds` of catalogues, i.e. galaxy-, member-, and group-catalogues
        for various combinations of input parameters.

        For example, if one wants to read the ``member`` galaxy catalogue,
        i.e. a catalogue with both `galaxy` and `group` information, for
        ``mocks`` catalogues, one could easily write:

        >>> memb_dir = catls_dir(catl_type='memb', catl_kind='mocks') # doctest: +SKIP

        This will return the path to the `member` galaxy catalogues directory
        with the synthetic catalogues in it. This method only returns
        the path of the directory, and not the path to the actual
        files in the directory. To accomplish this, one would need to use
        the function `catl_arr_extract`.
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
    def catl_arr_extract(self, catl_type='memb', catl_kind='mocks',
        ext='hdf5', print_filedir=False, return_len=False):
        r"""
        Extracts the list of galaxy/group catalogues.

        Parameters
        --------------
        catl_type : {``gal``, ``memb``, ``group``}, `str`
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
            If `True`, the path of the catalogue directory is printed
            onto the screen. This variable is set to ``False`` by default.


        Returns
        ---------
        catls_arr : `numpy.ndarray`, shape (N,)
            Array of the paths of the galaxy/group catalogues with specified
            parameters. The shape of the variable is (`N`,), where `N`,
            is the number of catalogues in the directory.

        Examples
        -----------
        This function is able to return the path to individual kinds of
        galaxy and group catalogues that meet some criteria based on the
        different combinations of input parameters.

        For example, if one wants to get the paths to the ``mocks`` ``member``
        galaxy catalogues with a ``halotype = 'fof'`` prescription, and
        following a ``clf_method = 2`` methodology, one could retrieve them
        by:
        
        >>> params_dict = {'catl_kind': 'mocks', 'halotype': 'fof', 'clf_method': 2}
        >>> catl_obj = CatlUtils(**params_dict) # doctest: +SKIP
        >>> catl_paths = catl_arr_extract(catl_type='memb', catl_kind='mocks', halotype='fof') # doctest: +SKIP

        This will return the paths to the `member` galaxy catalogues that
        meet those prescribed parameters.

        Notes
        ---------
        If no files are present or exist in the given directory, it will
        return an ``empty`` array of files.
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

    # Merges the `member` and `group` catalogues into a single DataFrame
    def catl_merge(self, catl_idx=0, return_memb_group=False,
        print_filedir=False):
        """
        Merges the `member` and `group` catalogues for a given set of
        input parameters, and returns a `modified` version of the galaxy
        group catalogues with added info about the galaxy groups.

        Parameters
        ------------
        catl_idx : `int`, optional
            Index of the catalogue to match. The index must be smaller
            than the number of files returned by function `~catl_arr_extract`
            for the given combination of parameters. This variable
            is set to ``0`` by default.

        return_memb_group : `bool`, optional
            If `True`, the function returns the `member` and `group`
            catalogues along with the `merged` catalogue. This variable is
            set to `False` by default.

        print_filedir : `bool`, optional
            If `True`, the output directory is printed onto the screen.
            This variable is set to `False` by default.

        Returns
        ------------
        merged_pd : `pandas.DataFrame`
            Combined DataFrame containing both ``galaxy`` and ``group``
            information for the given combination of parameters and the
            ``catl_idx``-th catalogue(s).

        memb_pd : `pandas.DataFrame`
            Member galaxy catalogue of the ``catl_idx``-th catalogue.
            This catalogue contains the information about the ``member``
            galaxies. This object is only returned if
            ``return_memb_group == True``, 

        group_pd : `pandas.DataFrame`
            Group galaxy catalogue of the ``catl_idx``-th catalogue.
            This catalogue contains the information about the ``galaxy groups``
            This object is only returned if ``return_memb_group == True``, 

        Raises
        ------------
        SDSSCatlUtils_Error : Exception from `~sdss_catl_utils.SDSSCatlUtils_Error`
            Program exception if input parameters are `not` accepted.

        Examples
        -----------
        This function will merge the ``member`` and ``group`` catalogues into
        a single `pandas.DataFrame` object, and it can be used to merge
        two existing catalogues. For example, to merge the `member` and
        `group` catalogues using the *default* parameters, one can easily
        write:
        
        >>> merged_pd = catl_merge() # doctest: +SKIP

        The resulting object will consist of a merge between the member and
        group catalogues, with the columns pertaining to **group** properties
        having a ``GG_`` attached to their names.

        However, if one wants to merge two **mock** catalogues with
        a ``clf_method = 2`` and ``halotype = 'fof'`` prescription, one could
        easily write this:

        >>> params_dict = {'catl_kind': 'mocks', 'halotype': 'fof', 'clf_method': 2}
        >>> catl_obj = CatlUtils(**params_dict) # doctest: +SKIP
        >>> merged_pd = catl_obj.catl_merge()

        Additionally, one could recover the ``member`` and ``group``
        catalogues as well:

        >>> merged_pd, memb_pd, group_pd = catl_obj.catl_merge(return_memb_group=True) # doctest: +SKIP

        For more information and examples, please refer to
        :ref:`quickstart_getting_started`.
        """
        file_msg = cfutils.Program_Msg(__file__)
        ## Checking input parameters
        # `catl_idx` - Type
        catl_idx_arr = (float, int, np.int64, np.int32, np.float32, np.float64)
        if not (isinstance(catl_idx, int)):
            msg = '{0} `catl_idx` ({1}) is not a valid input type!'
            msg = msg.format(file_msg, type(catl_idx))
            raise TypeError(msg)
        else:
            catl_idx = int(catl_idx)
        # `catl_idx` - Value
        idx_range = [0, 100]
        if not ((catl_idx >= idx_range[0]) and (catl_idx < idx_range[1])):
            msg = '{0} `catl_idx` ({1}) is not within the acceptable range! '
            msg = 'Range: `{2}`'
            msg = msg.format(file_msg, catl_idx, idx_range)
            raise ValueError(msg)
        # `return_memb_group` - Type
        if not (isinstance(return_memb_group, bool)):
            msg = '{0} `return_memb_group` ({1}) is not a valid input type!'
            msg = msg.format(file_msg, type(return_memb_group))
            raise TypeError(msg)
        # `print_filedir` - Type
        if not (isinstance(print_filedir, bool)):
            msg = '{0} `print_filedir` ({1}) is not a valid input type!'
            msg = msg.format(file_msg, type(print_filedir))
            raise TypeError(msg)
        ##
        ## Extracting catalogues given input parameters
        (   memb_arr,
            memb_len) = self.catl_arr_extract(  catl_type='memb',
                                                catl_kind=self.catl_kind,
                                                return_len=True,
                                                print_filedir=print_filedir)
        #
        # Checking if `catl_idx` is less than `memb_len`
        if (catl_idx > (memb_len - 1)):
            msg = '{0} `catl_idx` ({1}) is OUT of range ({2})!'
            msg = msg.format(file_msg, catl_idx, memb_len)
            raise ValueError(msg)
        #
        # Extracting group galaxy catalogue
        memb_path = memb_arr[catl_idx]
        # ith galaxy group catalogue
        group_path_base = self.catls_dir(   catl_type='group',
                                            catl_kind=self.catl_kind,
                                            print_filedir=print_filedir)
        # Modifying path to group catalogue
        if (self.catl_kind == 'mocks'):
            group_path = os.path.join(  group_path_base,
                            os.path.basename(memb_path).replace(
                                'memb','group'))
        elif (self.catl_kind == 'data'):
            group_path = os.path.join(  group_path_base,
                            os.path.basename(memb_path).replace(
                                'Gals', 'Group'))
        # Checking if files exist
        cfutils.File_Exists(group_path)
        # Reading in member and group galaxy catalogues into DataFrames
        memb_pd  = cfreaders.read_hdf5_file_to_pandas_DF(memb_path)
        group_pd = cfreaders.read_hdf5_file_to_pandas_DF(group_path)
        # Column keys for catalogues
        (   gm_key,
            id_key,
            galtype_key) = catl_keys(   self.catl_kind,
                                        perf_opt=self.perf_opt,
                                        return_type='list')
        # Matching keys from the group catalogue
        memb_id_unq  = np.unique(memb_pd[id_key].values)
        group_id_unq = np.unique(group_pd[id_key].values)
        if np.array_equal(memb_id_unq, group_id_unq):
            # Group column names
            group_cols = np.sort(group_pd.columns.values)
            # Sorting `memb_pd` by `id_key`
            # Member galaxy catalogue
            memb_pd.sort_values(by=id_key, inplace=True)
            memb_pd.reset_index(drop=True, inplace=True)
            # Group galaxy catalogue
            group_pd.sort_values(by=id_key, inplace=True)
            group_pd.reset_index(drop=True, inplace=True)
            # Reanming columns
            g_cols_dict = {ii: 'GG_' + ii for ii in group_cols}
            group_pd.rename(columns=g_cols_dict, inplace=True)
            group_pd.rename(columns={'GG_' + id_key: id_key}, inplace=True)
            # Merging the 2 DataFrames
            merged_pd = pd.merge(   left=memb_pd,
                                    right=group_pd,
                                    how='left',
                                    left_on=id_key,
                                    right_on=id_key)
        else:
            msg  = '{0} Lengths of the 2 DataFrames (`memb_pd`, `group_pd`) '
            msg += 'do not match! `memb_pd`: {1}, `group_pd`: {2}!'
            msg  = msg.format(file_msg, len(memb_id_unq), len(group_id_unq))
            raise SDSSCatlUtils_Error(msg)
        #
        # Returning DataFrames if necessary
        if return_memb_group:
            return_obj = (merged_pd, memb_pd, group_pd)
        else:
            return_obj = merged_pd

        return return_obj

