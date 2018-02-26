"""
Configuration file for ``atrain_match``. Most configuration options and
constants used in ``atrain_match`` are set in this file. However, there may
still be some modules which have internal constants defined. 
"""
def str2bool(v):
  return str(v).lower() in ("yes", "true", "t", "1")
import os

#========== Basic settings ==========#
# supported: AVHRR/VIIRS/MODIS/SEVIRI data:
IMAGER_INSTRUMENT = os.environ.get('IMAGER_INSTRUMENT', 'avhrr')
#: Resolution, in km, to use for data files. This setting is used 
#: to specify file names, sub-directories mm. Supported 1 or 5
RESOLUTION = int(os.environ.get('ATRAIN_RESOLUTION', 5))
#: Notice to get any matches in --sno_file processing
#: USE_ORBITS_THAT_STARTS_EXACTLY_AT_CROSS need to be False.
#: Save globbing time by using True for other processing
USE_ORBITS_THAT_STARTS_EXACTLY_AT_CROSS = True 
#: If set to True, no creation of reshaped files will be done.
#: Program will fail it there are no reshaped files already created.
#: This is useful for processing on SURFACES when original data are not
#: available any more. For normal processing let it be False.
USE_EXISTING_RESHAPED_FILES =  str2bool(
  os.environ.get('USE_EXISTING_RESHAPED_FILES', False))
WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE = True #to reduce disk usage

#: Use cloudmask to build cloud cover statistics. Traditionally CT was used.
USE_CMA_FOR_CFC_STATISTICS = True 

#========== Important time settings ==========#
#: Constant: Approximate duration of a satellite orbit in seconds 
SAT_ORBIT_DURATION = 50*60 #Not to large 
#: Allowed time deviation in seconds between AVHRR and CALIPSO/CloudSat matchup
sec_timeThr = 60*3

#========== Some paths and stuff ==========#
#: Best set as ENVIRON VARIABLES
#: Base directory for validation results   
_validation_results_dir = os.environ.get(
  'VALIDATION_RESULTS_DIR', 
  "/nobackup/smhid12/sm_kgkar/atrain_match_test_CALIPSOv4")
CONFIG_PATH = os.environ.get('ATRAINMATCH_CONFIG_DIR', './etc')
H4H5_EXECUTABLE = os.environ.get('H4H5_EXECUTABLE','h4toh5')

#========== Select imager and truth to validate ==========#
#: Choose one to validate
PPS_VALIDATION = str2bool(os.environ.get('PPS_VALIDATION', True))
CCI_CLOUD_VALIDATION = str2bool(os.environ.get('CCI_CLOUD_VALIDATION', False))
MAIA_CLOUD_VALIDATION = str2bool(os.environ.get('MAIA_CLOUD_VALIDATION', False))
CMA_PROB_VALIDATION = str2bool(os.environ.get('CMA_PROB_VALIDATION', False))
#: Turn off ISS and CLOUDSAT matching if never used
CALIPSO_MATCHING = True   #Notice can not be False if CALIPSO_REQUIRED = True
CLOUDSAT_MATCHING = True #Notice can not be False if CLOUDSAT_REQUIRED = True
ISS_MATCHING = False      #Notice can not be False if ISS_REQUIRED = True
AMSR_MATCHING = False     #Notice can not be False if AMSR_REQUIRED = True

#: Require matching. It is OK to have all False. Matching is still done
#: but program will not crach if it finds only CALIPSO data but not CloudSat.
CALIPSO_REQUIRED = False # Make progam fail if there is no CALIPSO match 
CLOUDSAT_REQUIRED = False # Make progam fail if there is no CloudSat match
ISS_REQUIRED = False # Make progam fail if there is no ISS match
AMSR_REQUIRED = False # Make progam fail if there is no ISS match

#========== Extra settings, to get things to reshaped file  ==========#
#: Save imager data also for warmest and coldest pixels:
SAVE_NEIGHBOUR_INFO = False
#: Search also for MODIS lvl2 data 
MATCH_MODIS_LVL2 = False #Only possible for MOIDS matching
#: To be able to match several PPS CTTH products in one file.
CTTH_TYPES = ["CTTH"] #["CTTHnn","CTTHold"]
#: Search also for calipso 5km aerosol data
MATCH_AEROSOL_CALIPSO = False

#========== Select cloudy/clear limits ==========#
#: Decide what is cloudy and what is not (CALIPSO)
#: It is used for CloudSat and ISS, their cloud_fraction is always 0 or 1    
#: For the combined 1km + 5km dataset cloud_fraction can only have values 
#: (0.0, 0.2, 0.4, 0.6, 0.8, 1.0). So the threshold should
#: really be set to 0.4, i.e., at least two 1 km columns should be cloudy!.
#: For calipso-v4  0, 1/15, 2/15 .. 1.0 is possible     
#: CLOUDS: pixels with cloud_fraction higher or equal to CALIPSO_CLOUDY_MIN_CFC
#: CLEAR:  pixels with cloud_fraction lower than CALIPSO_CLEAR_MAX_CFC
CALIPSO_CLOUDY_MIN_CFC = 0.5 #Tradition 0.64, KG used 0.5 for v2014 validation
CALIPSO_CLEAR_MAX_CFC = 0.5  #Tradition 0.33, KG used 0.5, 
                             #PPS development use low value 
CMA_PROB_CLOUDY_LIMIT = 50   #50% cloudy => cloud

#========== Select calipso version ==========#
#: Choose CALIPSO-CALIOP version
CALIPSO_version4 = False
CALIPSO_version3 = True

#: Only relevent for RESOLUTION-5. 
#: Both should not be True, but both can be False
ALSO_USE_1KM_FILES = False # Can be used with both CALIPSO-v3 and CALIPSO-v4
ALSO_USE_SINGLE_SHOT_CLOUD_CLEARED = True # Relevant for CALIPSO-v4

#========== Optical depth filtering ==========#
#: In the modes STANDARD (Ninas) and OPTICAL-DEPTH (KGs):
#: We consider cloud top to be OPTICAL_LIMIT_CLOUD_TOP down in the cloud layer. 
#: For clouds thinner than OPTICAL_LIMIT_CLOUD_TOP:
#: STANDARD: Use cloud base as cloud top
#: OPTICAL-DEPTH: Remove cloud and treat as clear
#: Set KG_OLD_METHOD_CLOUD_CENTER_CTTH_VALIDATION_HEIGHT to use cloud center as 
#: validation height insted.
#: In the 1km OPTICAL_DEPTH_THIN_IS_CLEAR pixels with total optical depth 
#: below OPTICAL_DETECTION_LIMIT are considered clear.
#: Settings that affect CTTH statistics should not affect CFC/CPY/CTY statistics
KG_OLD_METHOD_CLOUD_CENTER_AS_HEIGHT = False #MODE:OPTICAL_DEPTH, STATS:CTTH
           #If true use cloud layer center as validation height.
OPTICAL_LIMIT_CLOUD_TOP = 0.1 #STATS:CTTH, MODES:STANDARD and 
           # OPTICAL_DEPTH if KG_OLD_METHOD_CLOUD_CENTER_AS_HEIGHT = False
           # 1.0 used by others, we do it for each layer take something smaller
           # 1.0 gives quite bad results, for NN-CTTH 

COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC = True #MODE: ALL
                                                        #STATS: CTTH extra categories
OPTICAL_DETECTION_LIMIT = 0.2 # MODE: OPTICAL_DEPTH_THIN_IS_CLEAR 1km. STATS: CFC
                              # And to split CTTH in extra categories
                              # This should be the optical detection limit
#: Threshold for optical thickness. If optical thickness is below this 
#: value it will be filtered out in mode OPTICAL_DEPTH 
#MIN_OPTICAL_DEPTH = [0.20] # New formulation - allowing a set of values
MIN_OPTICAL_DEPTH = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                     0.50, 0.60, 0.70, 0.80, 0.90, 1.00] #MODE:OPTICAL_DEPTH, STATS:CFC 

#========== RESOLUTION-1 settings ==========#
#: For 5km processing these settings should have no effect
CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA = True #MODE: STANDARD-1km
ALSO_USE_5KM_FILES = True # Read variables from 5km, save in reshaped file
if (COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC or 
    CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA):
    ALSO_USE_5KM_FILES = True #5km data is required also for 1km processing

#========== Only for the separate CPP validation ==========#
#: important for cph_validate.py
VALIDATE_FOR_CPP_PIXELS = True # means validating for cloudtype only for 
                               # pixels where we also got cpp.cph values
CPP_REDUCE_PIXELS = int(os.environ.get('CPP_REDUCE_PIXELS', 0))
#: 1 means validate on a sub-set of the pixels, according to settings in
#:   util.py/reduce_cpp_data

#========== Process modes ==========#
#: Processing modes which can be handled for all resolutions
ALLOWED_MODES = ['BASIC', 'STANDARD']
if RESOLUTION == 1:
  if CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA:
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
    
#========== Often not changed ==========#
AREA = "no_area" #matching are no longer done for an area
COMPRESS_LVL = 6 #: Compresssion level for generated matched files (h5)
NODATA=-9
#: Recommended cloud threshold for the CloudSat cloud mask. In 5km data this
#: threshold has already been applied, so there is no reason to change it for
#: this data set.
CLOUDSAT_CLOUDY_THR = 30.0 
CALIPSO_FILE_LENGTH = 60*60 #calipso files are shorter 60 minutes
CLOUDSAT_FILE_LENGTH = 120*60 #cloudsat files are shorter 120 minutes
ISS_FILE_LENGTH = 60*60 #iss files are shorter 60 minutes 
AMSR_FILE_LENGTH = 60*60 #iss files are shorter 60 minutes 
#: Cloudsat data type (currently 'GEOPROF' are supported)
#: Traditionally also  'CWC-RVOD' where supported
CLOUDSAT_TYPE = 'GEOPROF'
#: Region configuaration file with area definitons
AREA_CONFIG_FILE = os.environ.get('AREA_CONFIG_FILE', './areas.def')
PPS_FORMAT_2012_OR_EARLIER = False
#: Set to true if you always want an avhrr orbit that starts before the cross
ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS = False

#========== Plotting ==========#
#: Processing modes for which plotting should also be performed
PLOT_MODES = ['No Plot']
#PLOT_MODES = ['BASIC']
PLOT_ONLY_PNG = True
#: MAXHEIGHT is used in the plotting. 
#: If None the maxheight is calculated from the highest cloud
#MAXHEIGHT = None
MAXHEIGHT = 18000

#========== Statistics setup ==========#
COMPILE_STATISTICS_TRUTH = ['iss','calipso','cloudsat', 'amsr']
#: DNT which statistics should be summarized 
#DNT_FLAG = ['ALL']
DNT_FLAG = ['ALL', 'DAY', 'NIGHT', 'TWILIGHT']
#: Surfaces for which statistics should be summarized 
SURFACES = PROCESS_SURFACES
#: In CASES list of dictionaries containing *satname*, *year*, 
#: and *month*, for which statistics should be summarized
#: Lats keep these last in teh config file
CASES_noaaa = [{'satname': 'noaa18', 'year': 2013, 'month': 3},
               {'satname': 'noaa18', 'year': 2013, 'month': 4},
               {'satname': 'noaa18', 'year': 2013, 'month': 5},
               {'satname': 'noaa18', 'year': 2013, 'month': 6},
               {'satname': 'noaa18', 'year': 2013, 'month': 7},
               {'satname': 'noaa18', 'year': 2013, 'month': 8},
               {'satname': 'noaa18', 'year': 2013, 'month': 9},
               {'satname': 'noaa18', 'year': 2013, 'month': 10},
               {'satname': 'noaa19', 'year': 2013, 'month': 3},
               {'satname': 'noaa19', 'year': 2013, 'month': 4},
               {'satname': 'noaa19', 'year': 2013, 'month': 5},
               {'satname': 'noaa19', 'year': 2013, 'month': 6},
               {'satname': 'noaa19', 'year': 2013, 'month': 7},
               {'satname': 'noaa19', 'year': 2013, 'month': 8},
               {'satname': 'noaa19', 'year': 2013, 'month': 9},
               {'satname': 'noaa19', 'year': 2013, 'month': 10}]

CASES = [ {'satname': 'noaa18', 'year': 2013, 'month': 6},
         {'satname': 'noaa18', 'year': 2013, 'month': 7},
         {'satname': 'noaa18', 'year': 2013, 'month': 8},
         {'satname': 'noaa18', 'year': 2013, 'month': 9},
         {'satname': 'noaa18', 'year': 2013, 'month': 10},
         {'satname': 'noaa18', 'year': 2013, 'month': 11},
         {'satname': 'noaa18', 'year': 2013, 'month': 12},
         {'satname': 'noaa18', 'year': 2008, 'month': 1},
         {'satname': 'noaa18', 'year': 2008, 'month': 2},
         {'satname': 'noaa18', 'year': 2008, 'month': 3},
         {'satname': 'noaa18', 'year': 2008, 'month': 4},
         {'satname': 'noaa18', 'year': 2008, 'month': 5},
         {'satname': 'noaa18', 'year': 2008, 'month': 6},
         {'satname': 'noaa18', 'year': 2008, 'month': 7},
         {'satname': 'noaa18', 'year': 2008, 'month': 8},
         {'satname': 'noaa18', 'year': 2008, 'month': 9},
         {'satname': 'noaa18', 'year': 2008, 'month': 10},
         {'satname': 'noaa18', 'year': 2008, 'month': 11},
         {'satname': 'noaa18', 'year': 2008, 'month': 12},
         {'satname': 'noaa18', 'year': 2009, 'month': 1},
         {'satname': 'noaa18', 'year': 2009, 'month': 2},
         {'satname': 'noaa18', 'year': 2009, 'month': 3},
         {'satname': 'noaa18', 'year': 2009, 'month': 4},
         {'satname': 'noaa18', 'year': 2009, 'month': 5},
         {'satname': 'noaa18', 'year': 2009, 'month': 6},
         {'satname': 'noaa18', 'year': 2009, 'month': 7},
         {'satname': 'noaa18', 'year': 2009, 'month': 8},
         {'satname': 'noaa18', 'year': 2009, 'month': 9},
         {'satname': 'noaa18', 'year': 2009, 'month': 10},
         {'satname': 'noaa18', 'year': 2009, 'month': 11},
         {'satname': 'noaa18', 'year': 2009, 'month': 12}]

CASES_npp =  [{'satname': 'npp', 'year': 2012, 'month': 6},
                       {'satname': 'npp', 'year': 2012, 'month': 7},
                       {'satname': 'npp', 'year': 2012, 'month': 8},
                       {'satname': 'npp', 'year': 2012, 'month': 9},
                       {'satname': 'npp', 'year': 2012, 'month': 10},
                       {'satname': 'npp', 'year': 2012, 'month': 11},
                       {'satname': 'npp', 'year': 2012, 'month': 12},
                       {'satname': 'npp', 'year': 2013, 'month': 01},
                       {'satname': 'npp', 'year': 2013, 'month': 02},
                       {'satname': 'npp', 'year': 2013, 'month': 03},
                       {'satname': 'npp', 'year': 2013, 'month': 04},
                       {'satname': 'npp', 'year': 2013, 'month': 05}]
CASES_V4test_NOAA18 =[{'satname': 'noaa18', 'year': 2006, 'month': 10},
         {'satname': 'noaa18', 'year': 2006, 'month': 11},
         {'satname': 'noaa18', 'year': 2006, 'month': 12}]
CASES_V4test_NOAA17 =[{'satname': 'noaa17', 'year': 2006, 'month': 10},
         {'satname': 'noaa17', 'year': 2006, 'month': 11},
         {'satname': 'noaa17', 'year': 2006, 'month': 12}]
CASES_modis =[{'satname': 'eos2', 'year': 2010, 'month': 1},
              {'satname': 'eos2', 'year': 2010, 'month': 2},
              {'satname': 'eos2', 'year': 2010, 'month': 3},
              {'satname': 'eos2', 'year': 2010, 'month': 4},
              {'satname': 'eos2', 'year': 2010, 'month': 5},
              {'satname': 'eos2', 'year': 2010, 'month': 6},
              {'satname': 'eos2', 'year': 2010, 'month': 7},
              {'satname': 'eos2', 'year': 2010, 'month': 8},
              {'satname': 'eos2', 'year': 2010, 'month': 9},
              {'satname': 'eos2', 'year': 2010, 'month': 10},
              {'satname': 'eos2', 'year': 2010, 'month': 11},
              {'satname': 'eos2', 'year': 2010, 'month': 12}]

CASES =[{'satname': 'noaa18', 'year': 2009, 'month': 1},
         {'satname': 'noaa18', 'year': 2009, 'month': 2},
         {'satname': 'noaa18', 'year': 2009, 'month': 3},
         {'satname': 'noaa18', 'year': 2009, 'month': 4},
         {'satname': 'noaa18', 'year': 2009, 'month': 5},
         {'satname': 'noaa18', 'year': 2009, 'month': 6},
         {'satname': 'noaa18', 'year': 2009, 'month': 7},
         {'satname': 'noaa18', 'year': 2009, 'month': 8},
         {'satname': 'noaa18', 'year': 2009, 'month': 9},
         {'satname': 'noaa18', 'year': 2009, 'month': 10},
         {'satname': 'noaa18', 'year': 2009, 'month': 11},
         {'satname': 'noaa18', 'year': 2009, 'month': 12}]

CASES = CASES_modis


