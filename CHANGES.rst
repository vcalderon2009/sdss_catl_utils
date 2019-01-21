0.0.3 (2019-01-21)
-----------------------

- Added new parameter ``sigma_clf_c``, which corresponds to the scatter in
  ``log(L)`` for central galaxies in the `conditional luminosity function`
  (CLF). The affected functions were:
  - `~sdss_catl_utils.models.catl_models.DownloadManager`
  - `~sdss_catl_utils.mocks_manager.catl_utils`
  - `~sdss_catl_utils.mocks_manager.catl_utils.check_input_params`
  - and more.

0.0.2 (2018-12-29)
-----------------------

- Added new functions to download catalogues to
  `~sdss_catl_utils.models.catl_models.DownloadManager`. It allows for
  the download of catalogues of specific characteristics.
- Added new function `~sdss_catl_utils.models.catl_models.CatlUtils`
  that acts as a tool for handling the synthetic catalogues.
- Added functions `~sdss_catl_utils.mocks_manager.catl_utils.catl_keys` and
  `~sdss_catl_utils.mocks_manager.catl_utils.catl_keys_prop` to facilitate
  the handling of the catalogues.
- Added unit-testing for `~sdss_catl_utils.mocks_manager.catl_utils` module.
- Added more functions for `~sdss_catl_utils.mocks_manager.catl_utils` module
  (`~sdss_catl_utils.mocks_manager.catl_utils.catl_prefix_main`,
  `~sdss_catl_utils.mocks_manager.catl_utils.check_input_params`)
- Added classes for the
  Conformity (`~sdss_catl_utils.models.catl_models.SDSSConformity`), 
  ML (`~sdss_catl_utils.models.catl_models.SDSSMLAnalysis`), and
  Catl (`~sdss_catl_utils.models.catl_models.SDSSCatlAnalysis`) analyses.


0.0.1 (2018-12-11)
-----------------------

- Initial release