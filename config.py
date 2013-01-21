"""
Configuration file for ``atrain_match``. Most configuration options and
constants used in ``atrain_match`` are set in this file. However, there may
still be some modules which have internal constants defined.

"""

import os
#When using 1km data, use the 5km data to filterout clouds to thin for VIIRS/AVHRR to see.
ALSO_USE_5KM_FILES = False
OPTICAL_DETECTION_LIMIT = 0.3
H4H5_EXECUTABLE = os.environ.get('H4H5_EXECUTABLE','/local_disk/opt/h4h5tools-2.2.1-linux-x86_64-static/bin/h4toh5')

VAL_CPP = os.environ.get('VAL_CPP', False)

# Imager Instrument on which PPS has been run (currently you can only run the
# atrain match on either AVHRR data or VIIRS data, not both):
IMAGER_INSTRUMENT = os.environ.get('IMAGER_INSTRUMENT', 'viirs')

#: Resolution, in km, to use for data files. This setting is used throughout
#: ``atrain_match`` to specify file names, sub-directories, and data handling.
#: Currently, 1 or 5 is supported
RESOLUTION = int(os.environ.get('ATRAIN_RESOLUTION', 1))

#: Base directory for validation results
_validation_results_dir = os.environ['VALIDATION_RESULTS_DIR']

#: Base directory where matchup files are stored. (TODO: Are these still used?)
SUB_DIR = "%s/Matchups" %_validation_results_dir

#: Base directory for files containing matched data from PPS and Calipso/Cloudsat
RESHAPE_DIR = "%s/Reshaped_Files" %_validation_results_dir

#: TODO: How is this directory used?
DATA_DIR = "%s/Data" %_validation_results_dir

#: Base directory for plots
PLOT_DIR = "%s/Plot" %_validation_results_dir

#: Base directory for statistics results
RESULT_DIR = "%s/Results" %_validation_results_dir

# Region configuaration file with area definitons
AREA_CONFIG_FILE = os.environ.get('AREA_CONFIG_FILE', './etc/areas.def')

#_satellite_data_dir = '/data/arkiv/proj/safworks/data'
_satellite_data_dir = '/local_disk/data/atrain_validation/data'
#: Base dir for PPS data
PPS_DATA_DIR = os.environ.get('PPS_DATA_DIR', _satellite_data_dir + '/pps')

#: Base dir for Cloudsat data
CLOUDSAT_DIR = os.environ.get('CLOUDSAT_DIR', _satellite_data_dir + '/cloudsat')

#: Cloudsat data type (currently 'GEOPROF' and 'CWC-RVOD' are supported)
#CLOUDSAT_TYPE = 'CWC-RVOD'
CLOUDSAT_TYPE = 'GEOPROF'

#: Base dir for Calipso data
CALIPSO_DIR = os.environ.get('CALIPSO_DIR', _satellite_data_dir + '/calipso')

#: Constant: Approximate duration of a satellite orbit in seconds
SAT_ORBIT_DURATION = 60*60 #Not to large
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

#: TODO: No description...
EMISS_MIN_HEIGHT = 2000.0

#: A value of 0.2-0.3 in cloud emissivity seems reasonable
EMISS_LIMIT = 0.2

CALIPSO_CLOUD_FRACTION = False
#: Processing modes which can be handled
# TODO: Split into latitude dependent area, snow-ice-land-sea area and day-night-twilight
ALLOWED_MODES = ['BASIC',
                 'BASIC_DAY',
                 'BASIC_NIGHT',
                 'BASIC_TWILIGHT', 
                 'EMISSFILT',       # Filter out cases with the thinnest topmost CALIPSO layers
                 'EMISSFILT_DAY',
                 'EMISSFILT_NIGHT',
                 'EMISSFILT_TWILIGHT',    
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
#latitude_filter = ['TROPICAL', 'SUB-TROPIC-ZONE', 'HIGH-LATITUDES', 'POLAR']
#land_filter = ['SNOW-COVER-LAND', 'SNOW-FREE-LAND', 'ICE-FREE-SEA', 'ICE-COVER-SEA']
#day_night_filter = ['DAY', 'NIGHT', 'TWILIGHT']
#modes = set()
#modes.add('BASIC') # no filtering
#for lat in latitude_filter:
#    modes.add(lat)
#    for land in land_filter:
#        modes.add('_'.join([lat, land]))
#        modes.add(land)
#        for day in day_night_filter:
#            modes.add('_'.join([lat, land, day]))
#            modes.add('_'.join([lat, day]))
#            modes.add(day)
    

#: Threshold for optical thickness. If optical thickness is below this value it
# will be filtered out.
MIN_OPTICAL_DEPTH = 0.5

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
    AREA = "arctic_europe_1km"
#    AREA = "arctic_super_5010"
#    AREA = "cea1km_test"
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


elif RESOLUTION == 5:
    DSEC_PER_AVHRR_SCALINE = 1.0/6. * 4 # A "work for the time being" solution.
    SWATHWD=409
    AREA = "cea5km_test"#"arctic_super_1002_5km"
    #: See :data:`CLOUDSAT_TRACK_RESOLUTION` in RESOLUTION = 1km
    CLOUDSAT_TRACK_RESOLUTION = 1.076#*5.0
    CLOUDSAT5KM_TRACK_RESOLUTION = CLOUDSAT_TRACK_RESOLUTION#*5.0
    ALLOWED_MODES.append('OPTICAL_DEPTH')      # Filter out cases with the thinnest topmost CALIPSO layers. Define MIN_OPTICAL_DEPTH above
    ALLOWED_MODES.append('OPTICAL_DEPTH_DAY')
    ALLOWED_MODES.append('OPTICAL_DEPTH_NIGHT')
    ALLOWED_MODES.append('OPTICAL_DEPTH_TWILIGHT')
else:
    raise ValueError("RESOLUTION == %s not supported" % str(RESOLUTION))


#: TODO: No description
COMPRESS_LVL = 6

#: TODO: No description
NLINES=6000

#: TODO: No description
NODATA=-9

#: Processing modes for which plotting should also be performed
PLOT_MODES = ['BASIC']
#PLOT_MODES = ['No Plot']

def subdir(self, date, *args, **kwargs):
    """
    This method is used in all ``file_finders`` instances created in 
    :mod:`cloudsat_calipso_avhrr_match`, for finding the correct sub-directory.
    
    It is required to handle the listed arguments *self*, and *date*, as well as
    any other arguments (:emphasis:`*args`) or keyword arguments 
    (:emphasis:`**kwargs`) passed to it. It should return a string containing
    the sub-directory where the file can be found.
    
    Note that :func:`subdir` is used in several different ``file_finders``
    classes (sub-classes of :class:`file_finders.SatelliteDataFileFinder`).
    ``self.__class__`` can be used to check for a particular ``file_finders``
    class, e.g.::
    
        def subdir(self, date, *args, **kwargs):
            from file_finders import PpsFileFinder
            if self.__class__ is PpsFileFinder:
                return "some/subdir/specific/to/your/pps/files"
            else:
                # Just use default subdir for the instance's class
                return self.__class__.subdir(self, date, *args, **kwargs)
    
    """
    
    # If default subdirs should be used, the following line should be uncommented
    return self.__class__.subdir(self, date, *args, **kwargs)
    
    # Example of how to set non-default subdirs for PpsFileFinder instances:
    from file_finders import PpsFileFinder #@UnresolvedImport
    if self.__class__ is PpsFileFinder:
        ending = kwargs.get('ending', None)
        if ending is None:
            ending = self.ending
        _dir = "%dkm/%d/%02d" % (RESOLUTION, date.year, date.month)
        if 'avhrr' in ending or 'viirs' in ending:
            return os.path.join(_dir,"import/PPS_data")                
        if 'sunsatangles' in ending:
            return os.path.join(_dir, "import/ANC_data")
        if 'nwp' in ending:
            return os.path.join(_dir,"import/NWP_data")
        for export_ending in ['cloudmask', 'cloudtype', 'ctth', 'precip']:
            if export_ending in ending:
                return os.path.join(_dir, "export")
    else:
        return self.__class__.subdir(self, date, *args, **kwargs)

#========== Statistics setup ==========#
#: List of dictionaries containing *satname*, *year*, and *month*, for which
#: statistics should be summarized
CASES = [{'satname': 'noaa18', 'year': 2006, 'month': 10},
         {'satname': 'noaa18', 'year': 2006, 'month': 11},
         {'satname': 'noaa18', 'year': 2006, 'month': 12},
         {'satname': 'noaa18', 'year': 2007, 'month': 1},
         {'satname': 'noaa18', 'year': 2007, 'month': 2},
         {'satname': 'noaa18', 'year': 2007, 'month': 3},
         {'satname': 'noaa18', 'year': 2007, 'month': 4},
         {'satname': 'noaa18', 'year': 2007, 'month': 5},
         {'satname': 'noaa18', 'year': 2007, 'month': 6},
         {'satname': 'noaa18', 'year': 2007, 'month': 7},
         {'satname': 'noaa18', 'year': 2007, 'month': 8},
         {'satname': 'noaa18', 'year': 2007, 'month': 9},
         {'satname': 'noaa18', 'year': 2007, 'month': 10},
         {'satname': 'noaa18', 'year': 2007, 'month': 11},
         {'satname': 'noaa18', 'year': 2007, 'month': 12},
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

CASES = [{'satname': 'npp', 'year': 2012, 'month': 6},
         {'satname': 'npp', 'year': 2012, 'month': 7},
         {'satname': 'npp', 'year': 2012, 'month': 8},
         {'satname': 'npp', 'year': 2012, 'month': 9},
         {'satname': 'npp', 'year': 2012, 'month': 10}]

#: PPS area definition (from ``acpg/cfg/region_config.cfg``) for which
#: statistics should be summarized
MAP = [AREA]

#: Base directory for ``atrain_match`` output to use when summarizing statistics.
#: Should contain the ``Results`` directory
MAIN_DATADIR = _validation_results_dir

#: TODO: No description yet...
# Seems like this can't be set to True if we run without Calipso 5km data!
# ??? Adam 2012-10-14
if CALIPSO_CLOUD_FRACTION == True:
    COMPILED_STATS_FILENAME = '%s/Results/compiled_stats_CCF' %MAIN_DATADIR
else:
    COMPILED_STATS_FILENAME = '%s/Results/compiled_stats' %MAIN_DATADIR
    
#: Surfaces for which statistics should be summarized
SURFACES = ["ICE_COVER_SEA", "ICE_FREE_SEA", "SNOW_COVER_LAND", "SNOW_FREE_LAND", \
            "COASTAL_ZONE", "TROPIC_ZONE", "TROPIC_ZONE_SNOW_FREE_LAND", "TROPIC_ZONE_ICE_FREE_SEA", \
            "SUB_TROPIC_ZONE", "SUB_TROPIC_ZONE_SNOW_FREE_LAND", "SUB_TROPIC_ZONE_ICE_FREE_SEA", \
            "HIGH-LATITUDES", "HIGH-LATITUDES_SNOW_FREE_LAND", "HIGH-LATITUDES_ICE_FREE_SEA", \
            "HIGH-LATITUDES_SNOW_COVER_LAND", "HIGH-LATITUDES_ICE_COVER_SEA", \
            "POLAR", "POLAR_SNOW_FREE_LAND", "POLAR_ICE_FREE_SEA", \
            "POLAR_SNOW_COVER_LAND", "POLAR_ICE_COVER_SEA"]

#: Filter types for which statistics should be summerized
FILTERTYPE = ['EMISSFILT']
if RESOLUTION == 5:
    FILTERTYPE.append('OPTICAL_DEPTH')
#: DAY NIGHT TWILIGHT FLAG that will be used
DNT_FLAG = ['ALL', 'DAY', 'NIGHT', 'TWILIGHT']
    
# The following are used in the old-style script interface
SATELLITE = ['noaa18', 'noaa19']
STUDIED_YEAR = ["2009"]
STUDIED_MONTHS = ['01', '07']
OUTPUT_DIR = "%s/Ackumulering_stat/Results/%s" % (MAIN_DATADIR, SATELLITE[0])


