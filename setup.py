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

lib_modules = ['common',  
               'runutils', 
               'calipso',
               'config',
               'cloudsat_calipso_imager_match', 
               'cloudsat_calipso_imager_prepare', 
               'cloudsat_calipso_imager_statistics',
               'cloudsat','amsr', 'iss',
               'extract_imager_along_track',
               'get_flag_info',
               'matchobject_io',
               'validate_cph_util',
               'stat_util',
               'pps_prototyping_util',
               'read_modis_products',
               'read_cloudproducts_maia',
               'read_cloudproducts_cci',
               'read_cloudproducts_and_nwp_pps',
               ]
script_modules = ['process_master',
                  'compile_stats']

setup(name='atrain_match',
      description="Library modules used in matching satellite swaths",
      author=("Jakob Malm <jakob.malm@smhi.se>, "
              "Adam Dybbroe <adam.dybbroe@smhi.se>, "
              "Erik Johansson <erik.johansson@smhi.se>, "
              "Nina Hakansson <nina.hakansson@smhi.se>, "
              "Karl-Göran Karlsson <kgkarl@smhi.se>"),
      author_email="FoUa@smhi.se",
      long_description='Software for matching imager data with lidar/radar and other truths',
      license='(?)',
      version=dist_version,
      provides=['atrain_match'],
      py_modules=lib_modules + script_modules,
      data_files=[('cfg', ['etc/atrain_match.cfg']), ],
      packages=['amsr_imager', 'statistics', 'plotting',
                'reshaped_files_plotting',
                'reshaped_files_analyze_scripts'],
      install_requires=['numpy', 'h5py'],
      zip_safe=False
      )
