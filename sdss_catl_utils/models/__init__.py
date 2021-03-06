# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The `~sdss_catl_utils.models` sub-package is responsible
for downloading galaxy and group galaxy catalogue data, storing hdf5
binaries and keeping a persistent memory of their location on disk
and associated metadata.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

from sdss_catl_utils.models import catl_models
from sdss_catl_utils.models import catl_models_template