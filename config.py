"""
Configuration file for ``atrain_match``. Most configuration options and
constants used in ``atrain_match`` are set in this file. However, there may
still be some modules which have internal constants defined.
 
"""

import os
PPS_FORMAT_2012_OR_EARLIER = False

#Set to true if you always want an avhrr orbit that starts before the cross
ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS = False

#Choose one to validate
PPS_VALIDATION = True
print "PPS_VALIDATION", PPS_VALIDATION
CCI_CLOUD_VALIDATION = False


ALSO_USE_5KM_FILES = True
COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC = False
COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE = False
OPTICAL_DETECTION_LIMIT = 0.2
if COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC:
    ALSO_USE_5KM_FILES = True #5km data is needed to split result on optical depth of top layer.

#This is only for 1km RESOULTION:
#In the modes STANDARD: the 1km calipso data will be filtered 
#using info from 5km data.
#Data devided on surfaces POLAR and so on will also be filtered.
#Results in BASIC are always unfiltered.
#If 5km data is used, results for surfaces will be based on the 
#STANDARD filtered results.
#If 5km data is not to filter data, results for surfaces will be 
#based on the BASIC unfiltered results.
#If USE_5KM_FILES_TO_FILTER_CALIPSO_DATA is True use 5km data to 
#filter fluffy clouds
# We consider the cloud top to be OPTICAL_LIMIT_CLOUD_TOP down 
#in the cloud. For clouds thinner than
# OPTICAL_LIMIT_CLOUD_TOP we use the cloud base as cloud top.
USE_5KM_FILES_TO_FILTER_CALIPSO_DATA = False # to get filtered cloudheight results in mode STANDARD and all modes that is not BASIC 
OPTICAL_LIMIT_CLOUD_TOP = 1.0 #also used by xxx in EUMETSAT
if USE_5KM_FILES_TO_FILTER_CALIPSO_DATA:
    ALSO_USE_5KM_FILES = True


H4H5_EXECUTABLE = os.environ.get('H4H5_EXECUTABLE','/local_disk/opt/h4h5tools-2.2.1-linux-x86_64-static/bin/h4toh5')
#H4H5_EXECUTABLE = '/software/apps/h4h5tools/2.2.1/i1214-hdf4-4.2.8-i1214-hdf5-1.8.9-i1214/bin/h4toh5'
PLOT_ONLY_PNG = True
#Nina 20150831 I have never needed the files printed here:
#They need lot of space, lets make it optional to have them!
DO_WRITE_COVERAGE = False 
DO_WRITE_DATA = False
#important for cph_validate.py
VAL_CPP = os.environ.get('VAL_CPP', False)
VALIDATE_FOR_CPP_PIXELS = True #means validating for cloudtype only for pixels where we also got cpp.cph values
CPP_REDUCE_PIXELS = int(os.environ.get('CPP_REDUCE_PIXELS', 0))
# 1 means validate on a sub-set of the pixels, according to settings in
#   util.py/reduce_cpp_data

# Imager Instrument on which PPS has been run (currently you can only run the
# atrain match on either AVHRR data or VIIRS data, not both):
IMAGER_INSTRUMENT = os.environ.get('IMAGER_INSTRUMENT', 'avhrr')

#: Resolution, in km, to use for data files. This setting is used throughout
#: ``atrain_match`` to specify file names, sub-directories, and data handling.
#: Currently, 1 or 5 is supported
RESOLUTION = int(os.environ.get('ATRAIN_RESOLUTION', 1))
if RESOLUTION == 1:
    AVHRR_SAT = 'NPP' #'pps'
    ALSO_USE_1KM_FILES = False
elif RESOLUTION == 5:
    AVHRR_SAT = 'NOAA18'
    ALSO_USE_1KM_FILES = True

#: Base directory for validation results
#_validation_results_dir = os.environ['VALIDATION_RESULTS_DIR']    
_validation_results_dir = os.environ.get('VALIDATION_RESULTS_DIR', "/nobackup/smhid9/sm_kgkar/atrain_match_2013")
CONFIG_PATH = os.environ.get('ATRAINMATCH_CONFIG_DIR', './etc')

#: Don't know how this directory is used...
MAIN_RUNDIR = os.getcwd()

# Region configuaration file with area definitons
AREA_CONFIG_FILE = os.environ.get('AREA_CONFIG_FILE', './areas.def')

#: Cloudsat data type (currently 'GEOPROF' and 'CWC-RVOD' are supported)
#CLOUDSAT_TYPE = 'CWC-RVOD'
CLOUDSAT_TYPE = 'GEOPROF'

#: Constant: Approximate duration of a satellite orbit in seconds
SAT_ORBIT_DURATION = 110*60 #Not to large
# If to large, cloudsat_calipso_avhrr_match.py takes wrong swath
# sometimes when swaths are close in time

#: Allowed time deviation in seconds between AVHRR and CALIPSO/CloudSat matchup
sec_timeThr = 60*20

#: Recommended cloud threshold for the CloudSat cloud mask. In 5km data this
#: threshold has already been applied, so there is no reason to change it for
#: this data set.
CLOUDSAT_CLOUDY_THR = 30.0 

#: MAXHEIGHT is used in the plotting. 
#: If None the maxheight is calculated from the highest cloud
#MAXHEIGHT = None
MAXHEIGHT = 18000

#: Range of allowed (AVHRR) satellite azimuth angles, in degrees
AZIMUTH_RANGE = (0., 360.)

#: Processing modes which can be handled
# TODO: Split into latitude dependent area, snow-ice-land-sea area and day-night-twilight
ALLOWED_MODES = ['BASIC',
                 'BASIC_DAY',
                 'BASIC_NIGHT',
                 'BASIC_TWILIGHT',    
                 'ICE_COVER_SEA',   # Restrict to ice cover over sea using NSIDC and IGBP data
                 'ICE_COVER_SEA_DAY',    
                 'ICE_COVER_SEA_NIGHT',    
                 'ICE_COVER_SEA_TWILIGHT',    
                 'ICE_FREE_SEA',    # Restrict to ice-free sea using NSIDC and IGBP data
                 'ICE_FREE_SEA_DAY',     
                 'ICE_FREE_SEA_NIGHT',     
                 'ICE_FREE_SEA_TWILIGHT',  
                 'SNOW_COVER_LAND', # Restrict to snow over land using NSIDC and IGBP data
                 'SNOW_COVER_LAND_DAY',
                 'SNOW_COVER_LAND_NIGHT', 
                 'SNOW_COVER_LAND_TWILIGHT',
                 'SNOW_FREE_LAND',  # Restrict to snow-free land using NSIDC and IGBP data
                 'SNOW_FREE_LAND_DAY',
                 'SNOW_FREE_LAND_NIGHT',
                 'SNOW_FREE_LAND_TWILIGHT',
                 'COASTAL_ZONE',    # Restrict to coastal regions using NSIDC data (mixed microwave region)
                 'COASTAL_ZONE_DAY',
                 'COASTAL_ZONE_NIGHT',
                 'COASTAL_ZONE_TWILIGHT',
                 'TROPIC_ZONE',     # Restrict to tropical regions lat <= +- 10
                 'TROPIC_ZONE_DAY',
                 'TROPIC_ZONE_NIGHT',
                 'TROPIC_ZONE_TWILIGHT',
                 'TROPIC_ZONE_SNOW_FREE_LAND',      # Restrict to tropical regions  +-10 < lat <= +-45
                 'TROPIC_ZONE_SNOW_FREE_LAND_DAY',
                 'TROPIC_ZONE_SNOW_FREE_LAND_NIGHT',
                 'TROPIC_ZONE_SNOW_FREE_LAND_TWILIGHT',
                 'TROPIC_ZONE_ICE_FREE_SEA',        # Restrict to tropical regions  +-10 < lat <= +-45
                 'TROPIC_ZONE_ICE_FREE_SEA_DAY',
                 'TROPIC_ZONE_ICE_FREE_SEA_NIGHT',
                 'TROPIC_ZONE_ICE_FREE_SEA_TWILIGHT',
                 'SUB_TROPIC_ZONE',     # Restrict to sub tropical regions +-10 < lat <= +-45
                 'SUB_TROPIC_ZONE_DAY',
                 'SUB_TROPIC_ZONE_NIGHT',
                 'SUB_TROPIC_ZONE_TWILIGHT',
                 'SUB_TROPIC_ZONE_SNOW_FREE_LAND',      # Restrict to tropical regions  +-10 < lat <= +-45
                 'SUB_TROPIC_ZONE_SNOW_FREE_LAND_DAY',
                 'SUB_TROPIC_ZONE_SNOW_FREE_LAND_NIGHT',
                 'SUB_TROPIC_ZONE_SNOW_FREE_LAND_TWILIGHT',
                 'SUB_TROPIC_ZONE_ICE_FREE_SEA',        # Restrict to tropical regions  +-10 < lat <= +-45
                 'SUB_TROPIC_ZONE_ICE_FREE_SEA_DAY',
                 'SUB_TROPIC_ZONE_ICE_FREE_SEA_NIGHT',
                 'SUB_TROPIC_ZONE_ICE_FREE_SEA_TWILIGHT',
                 'HIGH-LATITUDES',     # Restrict to tropical regions  +-45 < lat <= +-75
                 'HIGH-LATITUDES_DAY',
                 'HIGH-LATITUDES_NIGHT',
                 'HIGH-LATITUDES_TWILIGHT',
                 'HIGH-LATITUDES_SNOW_FREE_LAND',      # Restrict to tropical regions  +-45 < lat <= +-75
                 'HIGH-LATITUDES_SNOW_FREE_LAND_DAY',
                 'HIGH-LATITUDES_SNOW_FREE_LAND_NIGHT',
                 'HIGH-LATITUDES_SNOW_FREE_LAND_TWILIGHT',
                 'HIGH-LATITUDES_ICE_FREE_SEA',        # Restrict to tropical regions  +-45 < lat <= +-75
                 'HIGH-LATITUDES_ICE_FREE_SEA_DAY', 
                 'HIGH-LATITUDES_ICE_FREE_SEA_NIGHT',
                 'HIGH-LATITUDES_ICE_FREE_SEA_TWILIGHT',
                 'HIGH-LATITUDES_SNOW_COVER_LAND',      # Restrict to tropical regions  +-45 < lat <= +-75
                 'HIGH-LATITUDES_SNOW_COVER_LAND_DAY',
                 'HIGH-LATITUDES_SNOW_COVER_LAND_NIGHT',
                 'HIGH-LATITUDES_SNOW_COVER_LAND_TWILIGHT',
                 'HIGH-LATITUDES_ICE_COVER_SEA',        # Restrict to tropical regions  +-45 < lat <= +-75
                 'HIGH-LATITUDES_ICE_COVER_SEA_DAY',
                 'HIGH-LATITUDES_ICE_COVER_SEA_NIGHT',
                 'HIGH-LATITUDES_ICE_COVER_SEA_TWILIGHT',
                 'POLAR',     # Restrict to tropical regions  +-75
                 'POLAR_DAY',
                 'POLAR_NIGHT',
                 'POLAR_TWILIGHT', 
                 'POLAR_SNOW_FREE_LAND',      # Restrict to tropical regions  +-75
                 'POLAR_SNOW_FREE_LAND_DAY',
                 'POLAR_SNOW_FREE_LAND_NIGHT',
                 'POLAR_SNOW_FREE_LAND_TWILIGHT',
                 'POLAR_ICE_FREE_SEA',        # Restrict to tropical regions  +-75
                 'POLAR_ICE_FREE_SEA_DAY', 
                 'POLAR_ICE_FREE_SEA_NIGHT',
                 'POLAR_ICE_FREE_SEA_TWILIGHT',
                 'POLAR_SNOW_COVER_LAND',      # Restrict to tropical regions  +-75
                 'POLAR_SNOW_COVER_LAND_DAY',
                 'POLAR_SNOW_COVER_LAND_NIGHT',
                 'POLAR_SNOW_COVER_LAND_TWILIGHT',
                 'POLAR_ICE_COVER_SEA',        # Restrict to tropical regions  +-75
                 'POLAR_ICE_COVER_SEA_DAY',
                 'POLAR_ICE_COVER_SEA_NIGHT',
                 'POLAR_ICE_COVER_SEA_TWILIGHT']

#: Threshold for optical thickness. If optical thickness is below this value it will be filtered out.
#MIN_OPTICAL_DEPTH = 0.35 # Original formulation - only allowing one value
#MIN_OPTICAL_DEPTH = [0.35] # New formulation - allowing a set of values
MIN_OPTICAL_DEPTH = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]


print "RESOLUTION=", RESOLUTION
if RESOLUTION == 1:
    if IMAGER_INSTRUMENT == 'viirs':
        # VIIRS scan period is 1.7864 - see
        # D34862-07-01_C_CDFCB-X_Volume_VII-Part_1_ECR1017D_PR.pdf
        # Adam, 2012-10-21
        DSEC_PER_AVHRR_SCALINE = 1.7864 / 16.
        #DSEC_PER_AVHRR_SCALINE = 60 / 40. # ??? AD, 2012-Oct
        SWATHWD=3200
    else:
        DSEC_PER_AVHRR_SCALINE = 1.0/6. # Full scan period, i.e. the time
                                        # interval between two consecutive
                                        # lines (sec)
        SWATHWD=2048
    AREA = "no_area"
    #AREA = "npole"
    #AREA = "europa"
    #AREA = "cea1km_test"
    #: CloudSat sampling frequency in km (or rather the most likely
    #: resolution difference between CALIPSO 1 km datasets and
    #: the CloudSat 2B-GEOPROF dataset). Nominally CloudSat sampling
    #: rate is said to be 1.1 km but it seems as the CALIPSO sampling
    #: rate is not exactly 1 km but slightly above - as based on the
    #: optimised matching of plots of the two datasets.
    # TODO: Try too remove this one. Should be much better if it could be \
    # calculated instead inside atrain_match sugestion: \
    # cloudsat_track_resolution = len(calipso.longitude) / len(cloudsat.longitude)
    CLOUDSAT_TRACK_RESOLUTION = 1.076
    if USE_5KM_FILES_TO_FILTER_CALIPSO_DATA:
         ALLOWED_MODES.append('OPTICAL_DEPTH_THIN_IS_CLEAR')      # Filter out cases with the thinnest topmost CALIPSO layers. Define MIN_OPTICAL_DEPTH above
         ALLOWED_MODES.append('OPTICAL_DEPTH_THIN_IS_CLEAR_DAY')
         ALLOWED_MODES.append('OPTICAL_DEPTH_THIN_IS_CLEAR_NIGHT')
         ALLOWED_MODES.append('OPTICAL_DEPTH_THIN_IS_CLEAR_TWILIGHT')
         ALLOWED_MODES.append('STANDARD')      # Filter out cases with the thinnest topmost CALIPSO layers. Define MIN_OPTICAL_DEPTH above
         ALLOWED_MODES.append('STANDARD_DAY')
         ALLOWED_MODES.append('STANDARD_NIGHT')
         ALLOWED_MODES.append('STANDARD_TWILIGHT')

elif RESOLUTION == 5:
    DSEC_PER_AVHRR_SCALINE = 1.0/6. * 4 # A "work for the time being" solution.
    SWATHWD=409
    #AREA = "no_area"
    AREA = "cea5km_test"#"arctic_super_1002_5km"
    #: See :data:`CLOUDSAT_TRACK_RESOLUTION` in RESOLUTION = 1km
    CLOUDSAT_TRACK_RESOLUTION = 1.076#*5.0
    CLOUDSAT5KM_TRACK_RESOLUTION = CLOUDSAT_TRACK_RESOLUTION#*5.0


    ALLOWED_MODES.append('OPTICAL_DEPTH')      # Filter out cases with the thinnest topmost CALIPSO layers. Define MIN_OPTICAL_DEPTH above
    ALLOWED_MODES.append('OPTICAL_DEPTH_DAY')
    ALLOWED_MODES.append('OPTICAL_DEPTH_NIGHT')
    ALLOWED_MODES.append('OPTICAL_DEPTH_TWILIGHT')   #Let's run these cases separately but comment them all out if you always want to do filtering/KG
else:
    raise ValueError("RESOLUTION == %s not supported" % str(RESOLUTION))


#: TODO: No description
COMPRESS_LVL = 6

#: TODO: No description
NLINES=6000

#: TODO: No description
NODATA=-9

#: Processing modes for which plotting should also be performed
#PLOT_MODES = ['BASIC']
PLOT_MODES = ['No Plot']

#========== Statistics setup ==========#
#: List of dictionaries containing *satname*, *year*, and *month*, for which
#: statistics should be summarized
#CASES = [{'satname': 'npp', 'year': 2012, 'month': 06}]


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

#CASES =  CASES_npp
#CASES = CASES_noaaa + CASES_npp
#CASES = CASES_noaaa    
#: Surfaces for which statistics should be summarized
SURFACES = ["ICE_COVER_SEA", "ICE_FREE_SEA", "SNOW_COVER_LAND", "SNOW_FREE_LAND", \
            "COASTAL_ZONE", "TROPIC_ZONE", "TROPIC_ZONE_SNOW_FREE_LAND", "TROPIC_ZONE_ICE_FREE_SEA", \
            "SUB_TROPIC_ZONE", "SUB_TROPIC_ZONE_SNOW_FREE_LAND", "SUB_TROPIC_ZONE_ICE_FREE_SEA", \
            "HIGH-LATITUDES", "HIGH-LATITUDES_SNOW_FREE_LAND", "HIGH-LATITUDES_ICE_FREE_SEA", \
            "HIGH-LATITUDES_SNOW_COVER_LAND", "HIGH-LATITUDES_ICE_COVER_SEA", \
            "POLAR", "POLAR_SNOW_FREE_LAND", "POLAR_ICE_FREE_SEA", \
            "POLAR_SNOW_COVER_LAND", "POLAR_ICE_COVER_SEA"]


#: DAY NIGHT TWILIGHT FLAG that will be used
#DNT_FLAG = ['ALL']
DNT_FLAG = ['ALL', 'DAY', 'NIGHT', 'TWILIGHT']
    


