

This program is used to process and output statistics for the inter-comparison
of for example PPS results and CloudSat/CALIPSO observations. 

 * The main running program is: process_atrain_match.py which will do matchups or
    process_master.py will do matchups and retrieve result foreach case.
   
 * The program truth_imager_match.py is the main library for matching
   and truth_imager_make_statistics(_lib).py is the main libaries for statistics.

 * The compile_stat.py accumulates statistics uses the module statistics to
   accumulate statistics (only run if process_master have been used and there exist
   statistics).  

 * Program is updated wo be able to use CALIOP-CALIPSO, CPR (CloudSat), AMSR_E, Synop  
   or CATS (ISS) as truth. Modules used to handle the truths (read, reshape, etc):
      cloudsat.py
      calipso.py
      amsr.py
      iss.py
      synop.py
      mora.py

 * Program can read satellite data from: PPS, CCI and MAIA etc. When satllite data 
   comes from PPS-MODIS also modis lvl-2 data can be matched.
   Files to read imager satellite data:
      read_pps.py  
      read_maia.py
      read_cci.py  
      read_oca.py
      read_patmosx.py        
      read_modis_products.py

 * Format of matchup files, and reading and writing ot these can be found in matchobject_io.py.

 * Imager cloud top height datasets have been re-calculated to heights above mean
   sea level using CloudSat and CALIPSO elevation data

 * The MODIS cloud flag has been added to the extracted CALIPSO dataset. This
   enables direct comparisons to the MODIS cloud mask! Consequently,
   corresponding MODIS Cloud Mask statistics are calculated and printed.

 * The Vertical Feature Mask parameter in the CALIPSO dataset has been used to
   subdivide results into three cloud groups: Low, Medium and High. This has
   enabled an evaluation of PPS Cloud Type results and a further sub-division
   of Cloud Top Height results

 * The National Snow and Ice Data Center (NSIDC) ice and snow mapping results
   have been added to the extracted Calipso parameters. Together with the IGBP
   land use classification it is then possible to isolate the study to focus on
   one of the several surface categories.

 * Modules: reshaped_files_scr example script to plot things from reshaped files.

 * Configuration in the etc/atrain_match.cfg. There are some enironment variables that need to be set:
 
        ATRAINMATCH_CONFIG_DIR # path to atrain_match.cfg
        AREA_CONFIG_FILE       # only for plotting
        VALIDATION_RESULTS_DIR # path to validation main directory
        ATRAIN_RESOLUTION      # 1 or 5
