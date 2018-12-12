#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018_12-12
# Last Modified: 2018_12-12
# Vanderbilt University
from __future__ import absolute_import, division, print_function 
__author__     = ['Victor Calderon']
__copyright__  = ["Copyright 2018 Victor Calderon, "]
__email__      = ['victor.calderon@vanderbilt.edu']
__maintainer__ = ['Victor Calderon']
__all__        = [  'DownloadManager']

from cosmo_utils       import mock_catalogues as cm
from cosmo_utils       import utils           as cu
from cosmo_utils.utils import file_utils      as cfutils
from cosmo_utils.utils import file_readers    as cfreaders
from cosmo_utils.utils import work_paths      as cwpaths
from cosmo_utils.utils import stats_funcs     as cstats
from cosmo_utils.utils import geometry        as cgeom
from cosmo_utils.mock_catalogues import catls_utils as cmcu

class DownloadManager(object):
    """
    Class used to scrape the web for galaxy and group galaxy
    catalogue data and cache the downloaded catalogues.

    For list of available pre-processed galaxy- and group-galaxy
    catalogues proved by ``sdss_catl_utils``, see

    """
    def __init__(self):
        ""
        ""
        self.A = 'A'

    def create_new(self):
        """
        """
        print("This is a test")