#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018-12-24
# Last Modified: 2018-12-24
# Vanderbilt University
from __future__ import absolute_import, division, print_function 
__author__     = ['Victor Calderon']
__copyright__  = ["Copyright 2018 Victor Calderon, 2018"]
__email__      = ['victor.calderon@vanderbilt.edu']
__maintainer__ = ['Victor Calderon']
"""
Set of test functions for the `catl_utils` functions
"""

import numpy as np
import pytest
from sdss_catl_utils.mocks_manager import catl_utils
from sdss_catl_utils.custom_exceptions import SDSSCatlUtils_Error

## Functions

#### ----------------- Test `catl_keys` function - Types --------------------##

catl_keys_types_arr = [     ('data' , 'list', 3, list),
                            ('data' , 'dict', 3, dict),
                            ('mocks', 'list', 3, list),
                            ('mocks', 'dict', 3, dict) ]
@pytest.mark.parametrize('catl_kind, return_type, nelem, expected',
    catl_keys_types_arr)
def test_catl_keys_types_nelem(catl_kind, return_type, nelem, expected):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys` for input and 
    output variables.

    It verifies the `type` of the output returned by the function.

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} `str`
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues

    return_type : {'list', 'dict'} `str`
        Type of output to the be returned. This variable is set to `list`
        by default.

        Options:
            - 'list' : Returns the values as part of a list
            - 'dict' : Returns the values as part of a python dictionary

    nelem : `int`
        Expected number of elements inside the object returned by the function.

    expected : `str`
        Expected type of element from the `catl_keys` function
    """
    ## Constants
    perf_opt = False
    ## Running element
    output = catl_utils.catl_keys(catl_kind, return_type=return_type,
        perf_opt=perf_opt)
    ## Comparing against `expected` value - Type
    assert(isinstance(output, expected))
    ## Checking number of elements returned
    if isinstance(output, list):
        assert(len(output) == nelem)
    elif isinstance(output, dict):
        assert(len(output.keys()) == nelem)

#### ----------------- Test `catl_keys` function - Outputs ------------------##

catl_keys_return_arr = [ 'list' , 'dict']
catl_keys_output_arr = [('data' , False, ['M_h', 'groupid', 'galtype']),
                        ('data' , False, ['M_h', 'groupid', 'galtype']),
                        ('mocks', False, ['M_group', 'groupid', 'g_galtype']),
                        ('mocks', True, ['M_h', 'haloid', 'galtype'])]
@pytest.mark.parametrize('return_type', catl_keys_return_arr)
@pytest.mark.parametrize('catl_kind, perf_opt, expected', catl_keys_output_arr)
def test_catl_keys_outputs(catl_kind, perf_opt, return_type, expected):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys` for input and 
    output variables.

    It verifies the output returned by the function.

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues

    perf_opt : `bool`, optional
        Option for using a `perfect` mock catalogue.

    return_type : {'list', 'dict'} str
        Type of output to the be returned. This variable is set to `list`
        by default.

        Options:
            - 'list' : Returns the values as part of a list
            - 'dict' : Returns the values as part of a python dictionary

    expected : str
        Expected type of element from the `catl_keys` function
    """
    ## Running element
    output = catl_utils.catl_keys(catl_kind, perf_opt=perf_opt, 
        return_type=return_type)
    ## Comparing against `expected` value - Output
    if isinstance(output, list):
        np.testing.assert_equal(output, expected)
    elif isinstance(output, dict):
        out_keys = ['gm_key', 'id_key', 'galtype_key']
        out_vals = [output[xx] for xx in out_keys]
        np.testing.assert_equal(out_vals, expected)

#### ----------- Test `catl_keys` function - Errors - `catl_kind` -----------##

catl_keys_catl_kind_arr = [ 'data1', 'mocks1', 'NoMethod']
catl_keys_catl_perf_arr = [ True, False]
catl_keys_return_arr    = [ 'list' , 'dict']
@pytest.mark.parametrize('catl_kind', catl_keys_catl_kind_arr)
@pytest.mark.parametrize('return_type', catl_keys_return_arr)
@pytest.mark.parametrize('perf_opt', catl_keys_catl_perf_arr)
def test_catl_keys_catl_kind_errors_1(catl_kind, perf_opt, return_type):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys` for input and 
    output variables.

    It verifies if errors are raised when `catl_kind` is incorrect

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues
    """
    ## Running function
    with pytest.raises(SDSSCatlUtils_Error):
        output = catl_utils.catl_keys(catl_kind, perf_opt=perf_opt,
            return_type=return_type)

#### --------- Test `catl_keys` function - Errors - `return_type` -----------##

catl_keys_catl_kind_arr = ['data', 'mocks']
catl_keys_catl_perf_arr = [True, False]
catl_keys_return_arr    = [ 'list_no' , 'dict1', 'NoMethod']
@pytest.mark.parametrize('catl_kind', catl_keys_catl_kind_arr)
@pytest.mark.parametrize('return_type', catl_keys_return_arr)
@pytest.mark.parametrize('perf_opt', catl_keys_catl_perf_arr)
def test_catl_keys_catl_kind_errors_2(catl_kind, perf_opt, return_type):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys` for input and 
    output variables.

    It verifies if errors are raised when `catl_kind` is incorrect

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues
    """
    ## Running function
    with pytest.raises(SDSSCatlUtils_Error):
        output = catl_utils.catl_keys(catl_kind, perf_opt=perf_opt,
            return_type=return_type)

#### --------- Test `catl_keys` function - Errors - `return_type` -----------##

catl_keys_catl_kind_arr = ['data', 'mocks']
catl_keys_catl_perf_arr = [ 'NotBoolean', 1, 'mark', 1.2]
catl_keys_return_arr    = [ 'list' , 'dict']
@pytest.mark.parametrize('catl_kind', catl_keys_catl_kind_arr)
@pytest.mark.parametrize('return_type', catl_keys_return_arr)
@pytest.mark.parametrize('perf_opt', catl_keys_catl_perf_arr)
def test_catl_keys_catl_kind_errors_3(catl_kind, perf_opt, return_type):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys` for input and 
    output variables.

    It verifies if errors are raised when `catl_kind` is incorrect

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues
    """
    ## Running function
    with pytest.raises(TypeError):
        output = catl_utils.catl_keys(catl_kind, perf_opt=perf_opt,
            return_type=return_type)

#########-------------------------------------------------------------#########
#########-------------------------------------------------------------#########

#### ----------------- Test `catl_keys_prop` function - Types ---------------##

catl_keys_prop_info_arr  = ['memb', 'groups']
catl_keys_prop_types_arr = [('data' , 'list', 2, list),
                            ('data' , 'dict', 2, dict),
                            ('mocks', 'list', 2, list),
                            ('mocks', 'dict', 2, dict) ]
@pytest.mark.parametrize('catl_info', catl_keys_prop_info_arr)
@pytest.mark.parametrize('catl_kind, return_type, nelem, expected',
    catl_keys_prop_types_arr)
def test_catl_keys_prop_types_nelem(catl_kind, catl_info, return_type, nelem,
    expected):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys_prop` for input and 
    output variables.

    It verifies the `type` of the output returned by the function.

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues

    catl_info : {'memb', 'groups'} str, optional
        Option for which kind of catalogues to use.

    return_type : {'list', 'dict'} str
        Type of output to the be returned. This variable is set to `list`
        by default.

        Options:
            - 'list' : Returns the values as part of a list
            - 'dict' : Returns the values as part of a python dictionary

    nelem : int
        Expected number of elements inside the object returned by the function.

    expected : str
        Expected type of element from the `catl_keys_prop` function
    """
    ## Running element
    output = catl_utils.catl_keys_prop(catl_kind, catl_info=catl_info,
        return_type=return_type)
    ## Comparing against `expected` value - Type
    assert(isinstance(output, expected))
    ## Checking number of elements returned
    if isinstance(output, list):
        assert(len(output) == nelem)
    elif isinstance(output, dict):
        assert(len(output.keys()) == nelem)

#### ----------------- Test `catl_keys_prop` function - Output --------------##

catl_keys_prop_return_arr = [ 'list' , 'dict']
catl_keys_prop_output_arr = [('data' , 'memb', ['logssfr'    , 'logMstar_JHU']),
                        ('data' , 'groups' , ['logssfr_tot', 'logMstar_tot']),
                        ('mocks', 'memb', ['logssfr'    , 'logMstar']),
                        ('mocks', 'groups' , ['logssfr'    , 'logMstar'])]
@pytest.mark.parametrize('return_type', catl_keys_prop_return_arr)
@pytest.mark.parametrize('catl_kind, catl_info, expected', catl_keys_prop_output_arr)
def test_catl_keys_prop_outputs(catl_kind, catl_info, return_type, expected):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys_prop` for input and 
    output variables.

    It verifies the output returned by the function.

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues

    catl_info : {'memb', 'groups'} str, optional
        Option for which kind of catalogues to use.

    return_type : {'list', 'dict'} str
        Type of output to the be returned. This variable is set to `list`
        by default.

        Options:
            - 'list' : Returns the values as part of a list
            - 'dict' : Returns the values as part of a python dictionary

    expected : str
        Expected type of element from the `catl_keys_prop` function
    """
    ## Running element
    output = catl_utils.catl_keys_prop(catl_kind, catl_info=catl_info, 
        return_type=return_type)
    ## Comparing against `expected` value - Output
    if isinstance(output, list):
        np.testing.assert_equal(output, expected)
    elif isinstance(output, dict):
        out_keys = ['logssfr_key', 'logmstar_key']
        out_vals = [output[xx] for xx in out_keys]
        np.testing.assert_equal(out_vals, expected)

#### -------- Test `catl_keys_prop` function - Errors - `catl_kind` ---------##

catl_keys_prop_catl_kind_arr = [ 'data1', 'mocks1', 'NoMethod']
catl_keys_prop_return_arr    = [ 'list' , 'dict']
catl_keys_prop_catl_info_arr = [ 'memb', 'groups']
@pytest.mark.parametrize('catl_kind', catl_keys_prop_catl_kind_arr)
@pytest.mark.parametrize('return_type', catl_keys_prop_return_arr)
@pytest.mark.parametrize('catl_info', catl_keys_prop_catl_info_arr)
def test_catl_keys_prop_catl_kind_errors(catl_kind, catl_info, return_type):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys_prop` for input and 
    output variables.

    It verifies if errors are raised when `catl_kind` is incorrect

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues
    """
    ## Running function
    with pytest.raises(SDSSCatlUtils_Error):
        output = catl_utils.catl_keys_prop(catl_kind, catl_info=catl_info,
            return_type=return_type)

#### -------- Test `catl_keys_prop` function - Errors - `catl_info` ---------##

## Test `catl_keys_prop` function - Errors - `catl_info`
catl_keys_prop_catl_kind_arr = [ 'data', 'mocks']
catl_keys_prop_return_arr    = [ 'list' , 'dict']
catl_keys_prop_catl_info_arr = [ 'members_no', 'groups_Invalid', 1, 1.2]
@pytest.mark.parametrize('catl_kind', catl_keys_prop_catl_kind_arr)
@pytest.mark.parametrize('return_type', catl_keys_prop_return_arr)
@pytest.mark.parametrize('catl_info', catl_keys_prop_catl_info_arr)
def test_catl_keys_prop_catl_info_errors(catl_kind, catl_info, return_type):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys_prop` for input and 
    output variables.

    It verifies if errors are raised when `catl_info` is incorrect

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues
    """
    ## Running function
    with pytest.raises(SDSSCatlUtils_Error):
        output = catl_utils.catl_keys_prop(catl_kind, catl_info=catl_info,
            return_type=return_type)

#### ------- Test `catl_keys_prop` function - Errors - `return_type` --------##

catl_keys_prop_catl_kind_arr = [ 'data', 'mocks']
catl_keys_prop_return_arr    = [ 'list_no' , 'dict1', 'NoMethod']
catl_keys_prop_catl_info_arr = [ 'memb', 'groups']
@pytest.mark.parametrize('catl_kind', catl_keys_prop_catl_kind_arr)
@pytest.mark.parametrize('return_type', catl_keys_prop_return_arr)
@pytest.mark.parametrize('catl_info', catl_keys_prop_catl_info_arr)
def test_catl_keys_prop_return_type_errors(catl_kind, catl_info, return_type):
    """
    Tests the function:
        cosmo_utils.mock_catalogues.catl_utils.catl_keys_prop` for input and 
    output variables.

    It verifies if errors are raised when `return_type` is incorrect

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues
    """
    ## Running function
    with pytest.raises(SDSSCatlUtils_Error):
        output = catl_utils.catl_keys_prop(catl_kind, catl_info=catl_info,
            return_type=return_type)

#########-------------------------------------------------------------#########
#########-------------------------------------------------------------#########

#### --------------- Test `check_input_params` function - Types -------------##

input_arr = [   ('catl_kind', 'data'),
                ('hod_n', 1),
                ('halotype', 'fof'),
                ('clf_method', 1),
                ('clf_seed', 1234),
                ('dv', 1.),
                ('sample', '19'),
                ('type_am', 'mstar'),
                ('cosmo_choice', 'LasDamas'),
                ('perf_opt', True),
                ('remove_files', True),
                ('environ_name', 'test_name')]
@pytest.mark.parametrize('var_name, input_var', input_arr)
def test_check_input_params_types(input_var, var_name):
    """
    Checks the function `~sdss_catl_utils.mocks_manager.catl_utils.check_input_params`
    for input parameters.

    Parameters
    ------------
    input_var : `int`, `float`, `bool`, `str`
        Input variable to be evaluated.

    var_name : `str`
        Name of the input parameter being evaluated. This variable name
        must correspond to one of the keys in the `type` or `vals`
        dictionaries.
    """
    check_type = 'type'
    # Running function
    catl_utils.check_input_params(input_var, var_name, check_type=check_type)

#### --------------- Test `check_input_params` function - Values ------------##

input_arr = [   ('catl_kind', 'data'),
                ('catl_kind', 'mocks'),
                ('hod_n', 1),
                ('hod_n', 6),
                ('hod_n', 9),
                ('halotype', 'fof'),
                ('halotype', 'so'),
                ('clf_method', 1),
                ('clf_method', 2),
                ('clf_method', 3),
                ('sample', '19'),
                ('sample', '20'),
                ('sample', '21'),
                ('type_am', 'mstar'),
                ('type_am', 'mr'),
                ('cosmo_choice', 'LasDamas'),
                ('cosmo_choice', 'Planck')]
@pytest.mark.parametrize('var_name, input_var', input_arr)
def test_check_input_params_vals(input_var, var_name):
    """
    Checks the function `~sdss_catl_utils.mocks_manager.catl_utils.check_input_params`
    for input parameters.

    Parameters
    ------------
    input_var : `int`, `float`, `bool`, `str`
        Input variable to be evaluated.

    var_name : `str`
        Name of the input parameter being evaluated. This variable name
        must correspond to one of the keys in the `type` or `vals`
        dictionaries.
    """
    check_type = 'vals'
    # Running function
    catl_utils.check_input_params(input_var, var_name, check_type=check_type)

#### ---------- Test `check_input_params` function - Error - Type -----------##

input_arr = [   ('catl_kind', 1),
                ('hod_n', 'test'),
                ('halotype', None),
                ('clf_method', 'test'),
                ('clf_seed', '10'),
                ('dv', '1000'),
                ('sample', 19),
                ('type_am', 10),
                ('cosmo_choice', 123),
                ('perf_opt', 'None'),
                ('remove_files', 'True'),
                ('environ_name', 1)]
@pytest.mark.parametrize('var_name, input_var', input_arr)
def test_check_input_params_err_type(input_var, var_name):
    """
    Checks the function `~sdss_catl_utils.mocks_manager.catl_utils.check_input_params`
    for input parameters.

    Parameters
    ------------
    input_var : `int`, `float`, `bool`, `str`
        Input variable to be evaluated.

    var_name : `str`
        Name of the input parameter being evaluated. This variable name
        must correspond to one of the keys in the `type` or `vals`
        dictionaries.
    """
    check_type = 'type'
    # Running function
    with pytest.raises(TypeError):
        catl_utils.check_input_params(input_var, var_name,
            check_type=check_type)

#### ---------- Test `check_input_params` function - Errors - Values --------##

input_arr = [   ('catl_kind', 'data_no'),
                ('catl_kind', 'mocks_test'),
                ('hod_n', 11),
                ('hod_n', 63),
                ('hod_n', 103),
                ('halotype', 'fof_alt'),
                ('halotype', 'sos'),
                ('clf_method', 12),
                ('clf_method', 23),
                ('clf_method', 43),
                ('sample', '22'),
                ('sample', '34'),
                ('sample', '10'),
                ('type_am', '1_mstar'),
                ('type_am', '2_mr'),
                ('cosmo_choice', 'LasDamas_old'),
                ('cosmo_choice', 'Planck_new')]
@pytest.mark.parametrize('var_name, input_var', input_arr)
def test_check_input_params_err_vals(input_var, var_name):
    """
    Checks the function `~sdss_catl_utils.mocks_manager.catl_utils.check_input_params`
    for input parameters.

    Parameters
    ------------
    input_var : `int`, `float`, `bool`, `str`
        Input variable to be evaluated.

    var_name : `str`
        Name of the input parameter being evaluated. This variable name
        must correspond to one of the keys in the `type` or `vals`
        dictionaries.
    """
    check_type = 'vals'
    # Running function
    with pytest.raises(ValueError):
        catl_utils.check_input_params(input_var, var_name,
            check_type=check_type)

#### ---------- Test `check_input_params` function - Errors - KeyError --------##

input_arr = [   ('catl_kind_1', 'data_no'),
                ('hod_n_test', 103),
                ('_test_halotype', 'sos'),
                ('1123_clf_method', 43),
                ('_test_sample', '34'),
                ('type_type_am', '2_mr'),
                ('cosmo_choice_other_test', 'Planck_new')]
@pytest.mark.parametrize('var_name, input_var', input_arr)
def test_check_input_params_err_key(input_var, var_name):
    """
    Checks the function `~sdss_catl_utils.mocks_manager.catl_utils.check_input_params`
    for input parameters.

    Parameters
    ------------
    input_var : `int`, `float`, `bool`, `str`
        Input variable to be evaluated.

    var_name : `str`
        Name of the input parameter being evaluated. This variable name
        must correspond to one of the keys in the `type` or `vals`
        dictionaries.
    """
    check_type = 'vals'
    # Running function
    with pytest.raises(KeyError):
        catl_utils.check_input_params(input_var, var_name,
            check_type=check_type)

#########-------------------------------------------------------------#########
#########-------------------------------------------------------------#########

#### --------------- Test `catl_prefix_path` function - Types -------------##

input_arr = [
    ('data', 0, 'fof', 1, 1235, '19', 'mr', False, 'data/mr/Mr19'),
    ('mocks', 0, 'so', 1, 0, '20', 'mr', True, 'mocks/halos_so/dv_1.0/hod_model_0/clf_seed_0/clf_method_1/sigma_c_0.1417/mr/Mr20'),
    ('mocks', 6, 'so', 3, 10, '19', 'mr', False, 'mocks/halos_so/dv_1.0/hod_model_6/clf_seed_10/clf_method_3/sigma_c_0.1417/mr/Mr19')]
pytest_str  = 'catl_kind, hod_n, halotype, clf_method, clf_seed, sample, '
pytest_str += 'type_am, perf_opt, expected'
@pytest.mark.parametrize(pytest_str, input_arr)
def test_catl_prefix_path_inputs(catl_kind, hod_n, halotype, clf_method,
    clf_seed, sample, type_am, perf_opt, expected):
    """
    Checks the function
    `~sdss_catl_utils.mocks_manager.catl_utils.catl_prefix_path` for input
    parameters.

    Parameters
    -------------
    catl_kind : {``data``, ``mocks``} `str`
        Kind of catalogues to download. This variable is set to
        ``mocks`` by default.

        Options:
            - ``data``: Downloads the SDSS DR7 real catalogues.
            - ``mocks``: Downloads the synthetic catalogues of SDSS DR7.

    hod_n : `int`
        Number of the HOD model to use.

    halotype : {'so', 'fof'}, `str`
        Type of dark matter definition to use.

        Options:
            - ``so``: Spherical Overdensity halo definition.
            - ``fof``: Friends-of-Friends halo definition.

    clf_method : {1, 2, 3}, `int`
        Method for assigning galaxy properties to mock galaxies.
        This variable dictates how galaxies are assigned
        luminosities or stellar masses based on their galaxy type
        and host halo's mass.

        Options:
            - ``1``: Independent assignment of (g-r) colour, sersic, and specific star formation rate (`logssfr`)
            - ``2``: (g-r) colour dictates active/passive designation and draws values independently.
            - ``3``: (g-r) colour dictates active/passive designation, and assigns other galaxy properties for that given galaxy.

    clf_seed : `int`
        Value of the random seed used for the conditional luminosity function.

    sample : {'19', '20', '21'}, `str`
        Luminosity of the SDSS volume-limited sample to analyze.

        Options:
            - ``'19'``: :math:`M_r = 19` volume-limited sample
            - ``'20'``: :math:`M_r = 20` volume-limited sample
            - ``'21'``: :math:`M_r = 21` volume-limited sample

    type_am : {'mr', 'mstar'}, `str`
        Type of Abundance matching used in the catalogue. This

        Options:
            - ``'mr'``: Luminosity-based abundance matching used
            - ``'mstar'``: Stellar-mass-based abundance matching used.

    perf_opt : `bool`
        If `True`, it chooses to analyze the ``perfect`` version of
        the synthetic galaxy/group galaxy catalogues. Otherwise,
        it downloads the catalogues with group-finding errors
        included.

    expected : `str`
        Expected `path` to the set of catalogues
    """
    # Output path from function `catl_prefix_path`
    output_path = catl_utils.catl_prefix_path(  catl_kind=catl_kind,
                                                hod_n=hod_n,
                                                halotype=halotype,
                                                clf_method=clf_method,
                                                clf_seed=clf_seed,
                                                sample=sample,
                                                type_am=type_am,
                                                perf_opt=perf_opt)
    # Comparing expected with output
    assert(output_path == expected)

#### --------------- Test `catl_prefix_str` function - Types -------------##

input_arr = [
    ('data', 0, 'fof', 1, 1235, '19', 'mr', False, 'data_Mr19_am_mr'),
    ('mocks', 0, 'so', 1, 0, '20', 'mr', True, 'Mr20_halo_so_dv_1.0_hn_0_clfs_0_clfm_1_sigclf_0.1417_am_mr_pf_True'),
    ('mocks', 6, 'so', 3, 10, '19', 'mr', False, 'Mr19_halo_so_dv_1.0_hn_6_clfs_10_clfm_3_sigclf_0.1417_am_mr_pf_False')]
pytest_str  = 'catl_kind, hod_n, halotype, clf_method, clf_seed, sample, '
pytest_str += 'type_am, perf_opt, expected'
@pytest.mark.parametrize(pytest_str, input_arr)
def test_catl_prefix_path_inputs(catl_kind, hod_n, halotype, clf_method,
    clf_seed, sample, type_am, perf_opt, expected):
    """
    Checks the function
    `~sdss_catl_utils.mocks_manager.catl_utils.catl_prefix_str` for input
    parameters.

    Parameters
    -------------
    catl_kind : {``data``, ``mocks``} `str`
        Kind of catalogues to download. This variable is set to
        ``mocks`` by default.

        Options:
            - ``data``: Downloads the SDSS DR7 real catalogues.
            - ``mocks``: Downloads the synthetic catalogues of SDSS DR7.

    hod_n : `int`
        Number of the HOD model to use.

    halotype : {'so', 'fof'}, `str`
        Type of dark matter definition to use.

        Options:
            - ``so``: Spherical Overdensity halo definition.
            - ``fof``: Friends-of-Friends halo definition.

    clf_method : {1, 2, 3}, `int`
        Method for assigning galaxy properties to mock galaxies.
        This variable dictates how galaxies are assigned
        luminosities or stellar masses based on their galaxy type
        and host halo's mass.

        Options:
            - ``1``: Independent assignment of (g-r) colour, sersic, and specific star formation rate (`logssfr`)
            - ``2``: (g-r) colour dictates active/passive designation and draws values independently.
            - ``3``: (g-r) colour dictates active/passive designation, and assigns other galaxy properties for that given galaxy.

    clf_seed : `int`
        Value of the random seed used for the conditional luminosity function.

    sample : {'19', '20', '21'}, `str`
        Luminosity of the SDSS volume-limited sample to analyze.

        Options:
            - ``'19'``: :math:`M_r = 19` volume-limited sample
            - ``'20'``: :math:`M_r = 20` volume-limited sample
            - ``'21'``: :math:`M_r = 21` volume-limited sample

    type_am : {'mr', 'mstar'}, `str`
        Type of Abundance matching used in the catalogue. This

        Options:
            - ``'mr'``: Luminosity-based abundance matching used
            - ``'mstar'``: Stellar-mass-based abundance matching used.

    perf_opt : `bool`
        If `True`, it chooses to analyze the ``perfect`` version of
        the synthetic galaxy/group galaxy catalogues. Otherwise,
        it downloads the catalogues with group-finding errors
        included.

    expected : `str`
        Expected `path` to the set of catalogues
    """
    # Output path from function `catl_prefix_path`
    output_path = catl_utils.catl_prefix_str(  catl_kind=catl_kind,
                                                hod_n=hod_n,
                                                halotype=halotype,
                                                clf_method=clf_method,
                                                clf_seed=clf_seed,
                                                sample=sample,
                                                type_am=type_am,
                                                perf_opt=perf_opt)
    # Comparing expected with output
    assert(output_path == expected)























