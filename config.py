"""
Configuration file for ``atrain_match``. Most configuration options and constants
used in ``atrain_match`` are set in this file. However, there may
still be some modules which have internal constants defined.

"""

import os

def get_environ(name, default=None):
    """Get the environment variable *name*. If it is not defined, return 
    *default*."""
    try:
        return os.environ[name]
    except KeyError:
        return default


#: Resolution, in km, to use for data files. This setting is used throughout
#: ``atrain_match`` to specify file names, sub-directories, and data handling.
#: Currently, 1 or 5 is supported
RESOLUTION = 1

#: Base directory for ``atrain_match`` data
SAT_DIR = get_environ('SAT_DIR', default="/data/proj/safworks/jakob/data/atrain_validation")

#: Don't know how this directory is used...
MAIN_RUNDIR = os.getcwd()

#: Base directory for validation results
MAIN_DIR = get_environ('VALIDATION_RESULTS_DIR', "%s/results" % SAT_DIR)

#: Base directory where matchup files are stored. (TODO: Are these still used?)
SUB_DIR = "%s/Matchups" %MAIN_DIR

#: Base directory for files containing matched data from PPS and Calipso/Cloudsat
RESHAPE_DIR = "%s/Reshaped_Files" %MAIN_DIR

#: TODO: How is this directory used?
DATA_DIR = "%s/Data" %MAIN_DIR

#: Base directory for plots
PLOT_DIR = "%s/Plot" %MAIN_DIR

#: Base directory for statistics results
RESULT_DIR = "%s/Results" %MAIN_DIR


_satellite_data_dir = '/data/arkiv/proj/safworks/data'
#: Base dir for PPS data
PPS_DATA_DIR = get_environ('PPS_DATA_DIR', _satellite_data_dir + '/pps')

#: Base dir for Cloudsat data
CLOUDSAT_DIR = get_environ('CLOUDSAT_DIR', _satellite_data_dir + '/cloudsat')

#: Cloudsat data type (currently 'GEOPROF' and 'CWC-RVOD' are supported)
CLOUDSAT_TYPE = 'GEOPROF'

#: Base dir for Calipso data
CALIPSO_DIR = get_environ('CALIPSO_DIR', _satellite_data_dir + '/calipso')

#: Constant: Duration of a satellite orbit in seconds
SAT_ORBIT_DURATION = 90*60

#: CTTH file type to use in processing (One of 'ctth', 'ctth_opaque', and 'ctth_semitransparent')
CTTH_FILE = get_environ('CTTH_FILE', 'ctth')

#: Allowed time deviation in seconds between AVHRR and CALIPSO/CloudSat matchup
sec_timeThr = 60*20

#: Recommended cloud threshold for the CloudSat cloud mask. In 5km data this
#: threshold has already been applied, so there is no reason to change it for
#: this data set.
CLOUDSAT_CLOUDY_THR = 30.0 

#: TODO: This should be taken care of in plot function. See plot function for 1 km
MAXHEIGHT = 25000.0

#: Range of allowed (AVHRR) satellite azimuth angles, in degrees
AZIMUTH_RANGE = (0., 360.)

#: CloudSat sampling frequency in km (or rather the most likely
#: resolution difference between CALIPSO 1 km datasets and
#: the CloudSat 2B-GEOPROF dataset). Nominally CloudSat sampling
#: rate is said to be 1.1 km but it seems as the CALIPSO sampling
#: rate is not exactly 1 km but slightly above - as based on the
#: optimised matching of plots of the two datasets.
CLOUDSAT_TRACK_RESOLUTION = 1.076

#: See :data:`CLOUDSAT_TRACK_RESOLUTION`
CLOUDSAT5KM_TRACK_RESOLUTION = 1.076#*5.0

#: TODO: No description...
EMISS_MIN_HEIGHT = 2000.0

#: A value of 0.2-0.3 in cloud emissivity seems reasonable
EMISS_LIMIT = 0.2

#: Processing modes which can be handled
ALLOWED_MODES = ['BASIC',
             'EMISSFILT',        # Filter out cases with the thinnest topmost CALIPSO layers
             'ICE_COVER_SEA',    # Restrict to ice cover over sea using NSIDC and IGBP data
             'ICE_FREE_SEA',     # Restrict to ice-free sea using NSIDC and IGBP data
             'SNOW_COVER_LAND',  # Restrict to snow over land using NSIDC and IGBP data
             'SNOW_FREE_LAND',   # Restrict to snow-free land using NSIDC and IGBP data
             'COASTAL_ZONE']      # Restrict to coastal regions using NSIDC data (mixed microwave region)
             
if RESOLUTION == 1:
    DSEC_PER_AVHRR_SCALINE = 0.1667 # Full scan period, i.e. the time interval between two consecutive lines (sec)
    SWATHWD=2048
    AREA = "arctic_super_5010"

elif RESOLUTION == 5:
    DSEC_PER_AVHRR_SCALINE = 1.0/6*4 # A "work for the time being" solution.
    SWATHWD=409
    AREA = "arctic_super_1002_5km"
    ALLOWED_MODES.append('OPTICAL_DEPTH')      # Filter out cases with the thinnest topmost CALIPSO layers. Define MIN_OPTICAL_DEPTH below
else:
    raise ValueError("RESOLUTION == %s not supported" % str(RESOLUTION))

#: Threshold for optical thickness. If optical thickness is below this value it will be filtered out.
MIN_OPTICAL_DEPTH = 1

#: TODO: No description
COMPRESS_LVL = 6

#: TODO: No description
NLINES=6000

#: TODO: No description
NODATA=-9

#: Processing modes for which plotting should also be performed
PLOT_MODES = ['BASIC']

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
    from file_finders import PpsFileFinder
    if self.__class__ is PpsFileFinder:
        ending = kwargs.get('ending', None)
        if ending is None:
            ending = self.ending
        dir = "%dkm/%d/%02d" % (RESOLUTION, date.year, date.month)
        if 'avhrr' in ending:
            return os.path.join(dir,"import/PPS_data")                
        if 'sunsatangles' in ending:
            return os.path.join(dir, "import/ANC_data")
        if 'nwp' in ending:
            return os.path.join(dir,"import/NWP_data")
        for export_ending in ['cloudmask', 'cloudtype', 'ctth', 'precip']:
            if export_ending in ending:
                return os.path.join(dir, "export")
    else:
        return self.__class__.subdir(self, date, *args, **kwargs)


#========== Statistics setup ==========#
#: List of dictionaries containing *satname*, *year*, and *month*, for which
#: statistics should be summarized
CASES = [{'satname': 'noaa18', 'year': 2009, 'month': 1},
         {'satname': 'noaa18', 'year': 2009, 'month': 7},
         {'satname': 'noaa19', 'year': 2009, 'month': 7}]

#: PPS area definition (from ``acpg/cfg/region_config.cfg``) for which
#: statistics should be summarized
MAP = [AREA]

#: Base directory for ``atrain_match`` output to use when summarizing statistics.
#: Should contain the ``Results`` directory
MAIN_DATADIR = MAIN_DIR

#: TODO: No description yet...
COMPILED_STATS_FILENAME = '%s/Results/compiled_stats' %MAIN_DATADIR

#: Surfaces for which statistics should be summarized
SURFACES = ["ICE_COVER_SEA","ICE_FREE_SEA","SNOW_COVER_LAND","SNOW_FREE_LAND","COASTAL_ZONE"]

# The following are used in the old-style script interface
SATELLITE = ['noaa18', 'noaa19']
STUDIED_YEAR = ["2009"]
STUDIED_MONTHS = ['01', '07']
OUTPUT_DIR = "%s/Ackumulering_stat/Results/%s" % (MAIN_DATADIR, SATELLITE[0])
