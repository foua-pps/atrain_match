'''
Created on Oct 19, 2010

@author: a001696
'''
import os

def get_environ(name, default=None):
    """Get the environment variable *name*. If it is not defined, return 
    *default*."""
    try:
        return os.environ[name]
    except KeyError:
        return default


RESOLUTION = 1 #1 or 5
clsat_type = 1 #1=GEOPROG 2=CWC_RVOD
SAT_DIR = get_environ('SAT_DIR', default="/data/proj/saf/ejohansson/Satellite_Data")

match_file = "/data/proj/saf/ejohansson/SNO_tools/Snotimes/08/matchups_augsep_2008_mod.dat"

MAIN_RUNDIR = os.getcwd()
MAIN_DIR = get_environ('VALIDATION_RESULTS_DIR', MAIN_RUNDIR)
SUB_DIR = "%s/Matchups" %MAIN_DIR
RESHAPE_DIR = "%s/Reshaped_Files" %MAIN_DIR
DATA_DIR = "%s/Data" %MAIN_DIR
PLOT_DIR = "%s/Plot" %MAIN_DIR
RESULT_DIR = "%s/Results" %MAIN_DIR
PPS_DATA_DIR = get_environ('PPS_DATA_DIR', SAT_DIR) # The naming of the PPS env variable DATA_DIR is a bit unfortunate...
CLOUDSAT_DIR = "%s/CloudSat" % SAT_DIR
CLOUDSAT_TYPE = 'GEOPROF'
CALIPSO_DIR = "%s/Calipso" % SAT_DIR

SAT_ORBIT_DURATION = 90*60 # Duration of a satellite orbit in seconds

CTTH_FILE = get_environ('CTTH_FILE', 'ctth') # One of 'ctth', 'ctth_opaque', and 'ctth_semitransparent'

AREA5KM = "arctic_super_1002_5km"
AREA1KM = "baltrad1km"
sec_timeThr = 60*20 # Allowed time deviation in seconds between AVHRR and CALIPSO/CloudSat matchup

CLOUDSAT_CLOUDY_THR = 30.0  # Recommended cloud threshold for the CloudSat cloud mask. In 5km data this threshold have already been aplied therfore no reason to change it for thise data set. 
MAXHEIGHT = 25000.0
MAXHEIGHT = 12000.0
    
AZIMUTH_THR = 360.0

# CloudSat sampling frequency in km (or rather the most likely
# resolution difference between CALIPSO 1 km datasets and
# the CloudSat 2B-GEOPROF dataset). Nominally CloudSat sampling
# rate is said to be 1.1 km but it seems as the CALIPSO sampling
# rate is not exactly 1 km but slightly above - as based on the
# optimised matching of plots of the two datasets. /KG
CLOUDSAT_TRACK_RESOLUTION = 1.076
CLOUDSAT5KM_TRACK_RESOLUTION = 1.076#*5.0

EMISS_MIN_HEIGHT = 2000.0
EMISS_LIMIT = 0.2 # A value of 0.2-0.3 in cloud emissivity seems reasonable

DSEC_PER_AVHRR_SCALINE = 0.1667 # Full scan period, i.e. the time interval between two consecutive lines (sec)

ALLOWED_MODES = ['BASIC',
                 'EMISSFILT',        # Filter out cases with the thinnest topmost CALIPSO layers
                 'ICE_COVER_SEA',    # Restrict to ice cover over sea using NSIDC and IGBP data
                 'ICE_FREE_SEA',     # Restrict to ice-free sea using NSIDC and IGBP data
                 'SNOW_COVER_LAND',  # Restrict to snow over land using NSIDC and IGBP data
                 'SNOW_FREE_LAND',   # Restrict to snow-free land using NSIDC and IGBP data
                 'COASTAL_ZONE']      # Restrict to coastal regions using NSIDC data (mixed microwave region)
PLOT_MODES = ['BASIC']

if False: # Do we use this def?
    def subdir(self, satname, date):
        """This method is used by FileFinders for finding the correct subdir."""
        dir = "%s/%dkm/%d/%02d" % (satname, RESOLUTION, date.year, date.month)
        try:
            for ending in ['avhrr', 'sunsatangles', 'nwp_tsur']:
                if ending in self.ending:
                    dir = os.path.join(dir, 'import')
            for ending in ['cloudtype', 'ctth']:
                if ending in self.ending:
                    dir = os.path.join(dir, 'export')
        except (AttributeError, TypeError):
            # We're dealing with some other satellite data, e.g. Calipso or Cloudsat
            pass
        return dir