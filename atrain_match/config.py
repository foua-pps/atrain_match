# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
"""
Configuration file for ``atrain_match``. 
Configuration options not meant ot be changed.
The user configurable options are set in atrain_match.cfg
"""
def str2bool(v):
  return str(v).lower() in ("yes", "true", "t", "1")
import os

#========== Basic settings set as ENVIRON VARIABLES==========#
#: Resolution, in km, to use for data files. This setting is used 
#: to specify file names, sub-directories mm. Supported 1 or 5
RESOLUTION = int(os.environ.get('ATRAIN_RESOLUTION', 5))
#: Region configuaration file with area definitons, needed for plotting
AREA_CONFIG_FILE = os.environ.get('AREA_CONFIG_FILE', './areas.def')
#: Base directory for validation results   
_validation_results_dir = os.environ.get(
  'VALIDATION_RESULTS_DIR', 
  "/nobackup/smhid12/atrain_match_test_CALIPSOv4")
ATRAIN_MATCH_CONFIG_PATH = os.environ.get('ATRAINMATCH_CONFIG_DIR', './etc')

#========== Process modes ==========#
#: Processing modes which can be handled for all resolutions
ALLOWED_MODES = ['BASIC', 'STANDARD']
#if AMSR_MATCHING:
ALLOWED_MODES.append('SATZ_HIGH')
ALLOWED_MODES.append('SATZ_LOW')
if RESOLUTION == 1:
  #if CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA:
  ALLOWED_MODES.append('OPTICAL_DEPTH_THIN_IS_CLEAR')    
    # Filter thin CALIPSO pixels. Define OPTICAL_DETECTION_LIMIT
elif RESOLUTION == 5:
  ALLOWED_MODES.append('OPTICAL_DEPTH')      
  # Filter out cases with the thinnest topmost CALIPSO layers. 
  # Define MIN_OPTICAL_DEPTH above
else:
  raise ValueError("RESOLUTION == %s not supported" % str(RESOLUTION))  
PROCESS_SURFACES = [
  "ICE_COVER_SEA", "ICE_FREE_SEA", "SNOW_COVER_LAND", "SNOW_FREE_LAND", 
  "COASTAL_ZONE", 
  "TROPIC_ZONE", "TROPIC_ZONE_SNOW_FREE_LAND", "TROPIC_ZONE_ICE_FREE_SEA", 
  "SUB_TROPIC_ZONE", 
  "SUB_TROPIC_ZONE_SNOW_FREE_LAND", "SUB_TROPIC_ZONE_ICE_FREE_SEA", 
  "HIGH-LATITUDES", 
  "HIGH-LATITUDES_SNOW_FREE_LAND", "HIGH-LATITUDES_ICE_FREE_SEA", 
  "HIGH-LATITUDES_SNOW_COVER_LAND", "HIGH-LATITUDES_ICE_COVER_SEA", 
  "POLAR", "POLAR_SNOW_FREE_LAND", "POLAR_ICE_FREE_SEA", 
  "POLAR_SNOW_COVER_LAND", "POLAR_ICE_COVER_SEA"]
for surface in PROCESS_SURFACES:
  ALLOWED_MODES.append(surface)
# add the _NIGHT, _DAY and _TWILIGTH versions of all modes
for mode in list(ALLOWED_MODES):
  for DNT in ['_DAY', '_NIGHT', '_TWILIGHT']:
    ALLOWED_MODES.append(mode + DNT) 
    

COMPRESS_LVL = 6 #: Compresssion level for generated matched files (h5)
NODATA = -9
#: Recommended cloud threshold for the CloudSat cloud mask. In 5km data this
#: threshold has already been applied, so there is no reason to change it for
#: this data set.
CLOUDSAT_CLOUDY_THR = 30.0 
#: File lenghts, normally no need to update
CALIPSO_FILE_LENGTH = 60*60 #s calipso files are shorter 60 minutes
CLOUDSAT_FILE_LENGTH = 120*60 #s cloudsat files are shorter 120 minutes
ISS_FILE_LENGTH = 60*60 #s iss files are shorter 60 minutes 
AMSR_FILE_LENGTH = 60*60 #AMSR-Es  files are shorter 60 minutes 
SYNOP_FILE_LENGTH = 24*60 #s Our synop data comes in 1 day files
MORA_FILE_LENGTH = 24*60 #s Our MORA data comes in 1 day files

#: Maybe depricated?:
PPS_FORMAT_2012_OR_EARLIER = False
#: Surfaces for which statistics should be summarized 
SURFACES = PROCESS_SURFACES
DNT_FLAG = ['ALL', 'DAY', 'NIGHT', 'TWILIGHT']




