.. atrain_match documentation master file, created by sphinx-quickstart on Wed Nov  3 12:14:49 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to atrain_match's documentation!
========================================

``atrain_match`` contains functionality for matching AVHRR data with Cloudsat
and Calipso data, and producing statistics summaries for PPS products Cloud
Mask, Cloud Type, and Cloud Top Temperature and Height, compared with
corresponding Cloudsat and/or Calipso products.

Additionally, the :mod:`collect_last_day` module can be set up to collect
data files needed for matching and validation from a temporary directory (e.g.
``/data/24/...`` at SMHI) to an archive directory.

Contents:

.. toctree::
   :maxdepth: 2

   utilities
   
   find_crosses
   
   cloudsat_calipso_avhrr_match

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

