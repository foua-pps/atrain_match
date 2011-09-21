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
      py_modules=['find_crosses', 'common', 'merge_tles', 'track_correlation',
                  'runutils', 'amsr_avhrr_match'],
      packages=['amsr_avhrr'],
      install_requires=['pyephem', 'numpy', 'scipy'],
      zip_safe=False
      )