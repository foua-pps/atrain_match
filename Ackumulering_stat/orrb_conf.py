'''
Created on Oct 12, 2010

@author: a001696
'''


SATELLITE = ["metop02"]
RESOLUTION = ["1km"]
STUDIED_YEAR = ["2008"]
STUDIED_MONTHS = ["08","09"]
MAP = ["arctic_super_5010"]
MAIN_DATADIR = "/data/proj/saf/ejohansson/atrain_match" # Should contain the Results directory
CASES = [{'satname': 'noaa18', 'year': 2009, 'month': 1},
         {'satname': 'noaa18', 'year': 2009, 'month': 7},
         {'satname': 'noaa19', 'year': 2009, 'month': 7}]
OUTPUT_DIR = "%s/Ackumulering_stat/Results/%s" % (MAIN_DATADIR, SATELLITE[0])
SURFACES = ["ICE_COVER_SEA","ICE_FREE_SEA","SNOW_COVER_LAND","SNOW_FREE_LAND","COASTAL_ZONE"]
