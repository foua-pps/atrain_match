#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Scripts to match level2 products to CALIOP, CloudSat (CPR), Synop, etc."""

from setuptools import find_packages, setup

try:
    # HACK: https://github.com/pypa/setuptools_scm/issues/190#issuecomment-351181286
    # Stop setuptools_scm from including all repository files
    import setuptools_scm.integration
    setuptools_scm.integration.find_files = lambda _: []
except ImportError:
    pass

NAME = "atrain_match"
README = open('README.txt', 'r').read()

setup(name=NAME,
      description="Library modules used for matching satellite swaths",
      author=("Nina Hakansson <nina.hakansson@smhi.se>, "
              "Jakob Malm <jakob.malm@smhi.se>, "
              "Adam Dybbroe <adam.dybbroe@smhi.se>, "
              "Erik Johansson <erik.johansson@smhi.se>, "
              "Karl-Göran Karlsson <kgkarl@smhi.se>"),
      author_email="FoUa@smhi.se",
      long_description=README,
      classifiers=["Development Status :: 3 - Alpha",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License v3 " +
                   "or later (GPLv3+)",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering"],
      url="https://github.com:foua-pps/atrain_match",
      packages=find_packages(),
      scripts=['atrain_match/process_master.py', #
               'atrain_match/compile_stats.py', #
               'atrain_match/process_atrain_match.py', ],
      data_files=[('cfg', ['atrain_match/etc/atrain_match.cfg']),],
      zip_safe=False,
      use_scm_version=True,
      python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*',
      install_requires=['numpy', 'h5py', 'netCDF4', 'pandas' ],
)
