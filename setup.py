#!/usr/bin/env python
# -*- coding: utf-8 -*-

############################
# Get version from git tag or RELEASE-VERSION
#
# If this package is being installed as a dependency, chances are a module
# named 'version' has already been imported. Load the modules directly.
# This will also help if dependency (e.g. pyresample) tries to import a module
# named 'version'.
import imp
git_version = imp.load_source('atrain_match_git_version', 'version.py')
dist_version = git_version.get_git_version()

from setuptools import setup

lib_modules = ['common', 'find_crosses', 'runutils', 'track_correlation',
               'calipso',
               'cloudsat_calipso_avhrr_match', 'cloudsat_calipso_avhrr_plot',
               'cloudsat_calipso_avhrr_prepare', 'cloudsat_calipso_avhrr_statistics',
               'cloudsat_calipso_process_master', 'cloudsat_cwc', 'cloudsat',
               'cloudsat5km_cwc',
               'common', 'config', 'runutils', 'track_correlation',
               'filtfunc', 'radiance_tb_tables_kgtest', 'trajectory_plot']
script_modules = ['merge_tles', 'amsr_avhrr_match', 'amsr_avhrr_validate',
                  'validate_cph', 'validate_cph_all',
                  'clean_sno_results', 'collect_last_day', 'fetch_nwp',
                  'nwp_profile', 'process_master', 'run_aapp_on_ears',
                  'compile_stats']

setup(name='atrain_match',
      description="Library modules used in matching satellite swaths",
      author=("Jakob Malm <jakob.malm@smhi.se>, "
              "Erik Johansson <erik.johansson@smhi.se>, "
              "Karl-Göran Karlsson <kgkarl@smhi.se>"),
      author_email="FoUa@smhi.se",
      url='http://nwcsaf.org',
      long_description='long description',
      license='EUMETSAT NWCSAF license (?)',
      version=dist_version,
      provides=['atrain_match'],
      py_modules=lib_modules + script_modules,
      packages=['amsr_avhrr', 'statistics'],
      install_requires=['pyephem', 'numpy', 'scipy'],
      zip_safe=False
      )