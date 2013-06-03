"""
Use this script to validate the CPP cloud phase (cph) product

"""

import os
from amsr_avhrr.util import get_avhrr_lonlat, get_cpp_product, get_avhrr_lonlat_and_time_cci,  get_cpp_product_cci
from amsr_avhrr.match import match_lonlat
from runutils import process_scenes
from cloudsat_calipso_avhrr_match import (find_calipso_files,
                                          find_files_from_avhrr,
                                          get_satid_datetime_orbit_from_fname)

import numpy as np
from amsr_avhrr.util import get_avhrr_lonlat
from amsr_avhrr.util import get_avhrr_time
from read_cloudproducts_cci import cci_read_prod
import time
import h5py
import numpy as np
import logging
from config import _validation_results_dir
from config import RESOLUTION, PPS_VALIDATION, CCI_CLOUD_VALIDATION
from amsr_avhrr.util import write_data, write_data_to_open_file
logger = logging.getLogger(__name__)

from datetime import datetime
TAI93 = datetime(1993, 1, 1)

#: Directory for mapper files
MATCH_DIR = os.environ.get('MATCH_DIR', _validation_results_dir+'/CPP_MATCH_DIR_CT')
if not os.path.exists(MATCH_DIR):
    logger.info("Creating cpp match dir: %s"%(MATCH_DIR ))
    os.makedirs(MATCH_DIR)

#Time threshold, i.e. max time diff to be considered as a match
TIME_THR=60*60.0
#NB! I'm not sure there really is a sorting due to TIME_THR. Do check in
#    the time_diff dataset after running! /Sara Hornquist 2012-06-18
#NB! This should now be working, and time threshold should be used.
# Data are still present in .h5 files but nog in .h5--values.h5 files.
# The problem was in amsr_avhrr/match.py from_file method 
# The saved time_threshold were never read again /Nina Hakansson 2013-02-06
#: h5py compression settings (True, or an integer in range(10))
_COMPRESSION = True

#: Calipso cloud phase bits
CALIPSO_PHASE_BITS = range(5, 7)

#: Calipso cloud phase values
CALIPSO_PHASE_VALUES = dict(unknown=0,
                            ice=1,
                            water=2,
                            horizontal_oriented_ice=3)

#: Water (no mixed) value
CALIPSO_WATER_VALUE = 2

#: Calipso quality bits
CALIPSO_QUAL_BITS = range(7, 9)

#: Calipso quality values
CALIPSO_QUAL_VALUES = dict(none=0,
                           low=1,
                           medium=2,
                           high=3)


#: CPP cph value meanings
CPP_PHASE_VALUES = dict(no_cloud=0,
                        liquid=1,
                        ice=2,
                        mixed=3,
                        non_opice=4,
                        non_opwater=5,
                        uncertain=6,
                        no_observation=-1)


#: Cloud type phase value equivalents
CTYPE_PHASE_BITS = {'Not processed or undefined': 1,
                    'Water': 2,
                    'Ice': 4,
                    'Tb11 below 260K': 8}



def get_calipso_lonlat(calipso_filename):
    with h5py.File(calipso_filename, 'r') as f:
        if RESOLUTION==1:    
            lon = f['Longitude'][:].ravel()
            lat = f['Latitude'][:].ravel()
        if RESOLUTION==5:
            lon = f['Longitude'][:,1].ravel()
            lat = f['Latitude'][:,1].ravel()
    if f:
        logger.info("Closing file %s"%calipso_filename)
        f.close()        
    return lon, lat


def get_calipso_time(filename):
    """
    Get time of each scanline in Calipso file *filename*.
    
    Note::
    
        Time is converted from seconds since 1993-01-01 to seconds since
        1970-01-01.
    
    """
    with h5py.File(filename, 'r') as f:
        if RESOLUTION==1: 
            sec1993 = f['Profile_Time'][:]
        if RESOLUTION==5:   
            sec1993 = f['Profile_Time'][:,1]
    if f:
        logger.info("Closing file %s"%filename)
        f.close()
    from calendar import timegm
    epoch_diff = timegm(TAI93.utctimetuple())
    
    sec1970 = sec1993.round().astype(np.int64) + epoch_diff
    
    return sec1970

def get_bits(value, bits, shift=False):
    """
    Returns value for bits *bits* in *value*.
    
    Examples
    
    >>> get_bits(6, [0, 1])
    2
    >>> get_bits(6, [1, 2])
    6
    
    If *shift* is True, shift the obtained value by min(bits) bits:
    
    >>> get_bits(6, [1, 2], shift=True)
    3
    
    """
    selected = value & sum([2**i for i in bits])
    if shift:
        return selected >> min(bits)
    return selected


def get_calipso_phase(calipso_filename, qual_min=CALIPSO_QUAL_VALUES['medium'],
                      max_layers=1):
    """
    Returns Calipso cloud phase.
    
    Pixels with quality lower than *qual_min* are masked out.
    
    Screen out pixels with more than *max_layers* layers.
    
    """
    with h5py.File(calipso_filename,'r') as f:
        features = f['Feature_Classification_Flags'][:]
    if f:
        f.close()
    # Reduce to single layer, masking any multilayer pixels
    features = np.ma.array(features[:, 0],
                           mask=(features[:, max_layers:] > 1).any(axis=-1))
    
    phase = get_bits(features, CALIPSO_PHASE_BITS, shift=True)
    qual = get_bits(features, CALIPSO_QUAL_BITS, shift=True)
    # Don't care about pixels with lower than *qual_min* quality
    return np.ma.array(phase, mask=qual < qual_min)

def get_calipso_igbp(calipso_filename):
    """
    Returns Calipso igbp.
    
    """
    with h5py.File(calipso_filename,'r') as f:
        igbp = f['IGBP_Surface_Type'][:]
    if f:
        f.close()
    return igbp.ravel()

def find_calipso_files_from_avhrr_filename(avhrr_filename, options):
    """
    Find Calipso files matching *avhrr_filename*. Returns a list of file paths.
    
    """
    import config
    sl_ = os.path.basename(avhrr_filename).split('_')
    satname = sl_[0]
    date_time = datetime.strptime(sl_[1] + sl_[2], '%Y%m%d%H%M')
    return find_calipso_files(date_time, options, values={})

    # Limit matching to AMSR-E files starting 45 min (duration of one half
    # orbit) before up to 20 min (duration of one EARS AVHRR swath) after the
    # start of the AVHRR swath
    #amsr_finder = CalipsoFileFinder(time_window=(-45 * 60, 20 * 60))
    #return amsr_finder.find(parsed['datetime'])


def match_with_calipso(calipso_filename, avhrr_filename, sunsat_filename, radius_of_influence=1e3,
          time_threshold=None, n_neighbours=1):
    """
    Find matching indices in AVHRR array for each element in Calipso swath.
    
    Arguments:
    
        cal_filename: string
            full path of Calipso HDF5 file
        avhrr_filename: string
            full path of AVHRR PPS HDF5 file
        radius_of_influence: float
            radius of influence in meters in pixel-pixel matching (default:
            1000 m)
        time_threshold: float
            largest absolute time difference to include in match
        n_neighbours: int
            number of nearest AVHRR neighbours to use
    
    Returns:
    
        mapper: `MatchMapper` instance.
    
    """

    #get_calipso_lonlat and get_calipso_time in this file!!
    logger.debug("Getting AVHRR lon/lat")
    if PPS_VALIDATION:
        avhrr_lonlat = get_avhrr_lonlat(sunsat_filename)
        avhrr_time = get_avhrr_time(sunsat_filename)
    if CCI_CLOUD_VALIDATION:
        avhrr_lonlat, avhrr_time  = get_avhrr_lonlat_and_time_cci(avhrr_filename)    
    logger.debug("Getting Calipso lon/lat")
    logger.info("Get caliop lonlat:") 
    calipso_lonlat = get_calipso_lonlat(calipso_filename)         
    logger.debug("Matching AVHRR to Calipso lon/lat")
    mapper = match_lonlat(avhrr_lonlat, calipso_lonlat, n_neighbours=1)

    #Also calculate the time diff 

    logger.info("Get caliop time:") 
    calipso_time = get_calipso_time(calipso_filename)

    tim1 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(np.min(calipso_time)))
    tim2 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(np.max(calipso_time)))
    logger.info("Starttime caliop: %s, end time: %s"%(tim1, tim2))

    time_diff = np.abs(avhrr_time[mapper.rows] -
                       calipso_time.reshape((calipso_time.size,1))).astype(np.float32)
    time_diff =np.where(mapper._pixel_mask, np.inf, time_diff) #Needed to avoid nonses values outside swath

    mapper.time_diff = time_diff
    mapper.time_threshold = time_threshold

    logger.debug("Time diff (min, max): %r" % ((time_diff.min(),
                                                time_diff.max()),))
    
    return mapper



def process_noaa_scene(avhrr_file, options, cloudtype=False, **kwargs):
    """
    Match this noaa scene with cloudsat scenes and process.
    
    """
    from pps_basic_configure import AVHRR_DIR, OUTPUT_DIR
    import pps_arguments 
    #ppsarg = pps_arguments.create_pps_argument(pps_arguments.USE_SUNSATANGLES_GLOB,
    #                                           pps_arguments.SATELLITE_PROJ,
    #                                           satname,orbit)
    #avhrr_filename = os.path.join(AVHRR_DIR, ppsarg.files.avhrr)
    logger.info("Processing file: %s"%avhrr_file)    
    if (PPS_VALIDATION):
            pps_files = find_files_from_avhrr(avhrr_file, options)
            sunsat_filename  = pps_files.sunsatangles
            calipso_filenames = find_calipso_files_from_avhrr_filename(avhrr_file, options)

    if (CCI_CLOUD_VALIDATION):
        #avhrr_file = "20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
        values = get_satid_datetime_orbit_from_fname(avhrr_file)
        date_time = values["date_time"]
        calipso_filenames = find_calipso_files(date_time, options, values)
        cpp_filename = None#avhrr_filename
        sunsat_filename = None

    logger.debug("Found Calipso files: %r" % calipso_filenames)
    if cloudtype and PPS_VALIDATION:
        cpp_filename = os.path.join(OUTPUT_DIR, pps_files.cloudtype)
        
        # Replace get_cpp_product(), imported from amsr_avhrr.util with
        # a function which reads cloud phase from the cloud type product
        def get_ctype_phase(filename, *args, **kwargs):
            from ppshdf_cloudproducts import read_cloudtype
            # read_cloudtype has to have wtype=1, or it will segfault,
            # at least with ACPG <= r2_46
            cloudtype = read_cloudtype(filename, wtype=1, wqual=0, wphase=1)
            # Initialize phase array to 'no_observation'
            phase = cloudtype.phaseflag * 0 + CPP_PHASE_VALUES['no_observation']
            phase[(cloudtype.phaseflag & CTYPE_PHASE_BITS['Water']) != 0] = \
                CPP_PHASE_VALUES['liquid']
            phase[(cloudtype.phaseflag & CTYPE_PHASE_BITS['Ice']) != 0] = \
                CPP_PHASE_VALUES['ice']
            phase[((cloudtype.phaseflag & CTYPE_PHASE_BITS['Water']) != 0) &
                  ((cloudtype.phaseflag & CTYPE_PHASE_BITS['Ice']) != 0)] = \
                CPP_PHASE_VALUES['mixed']
            
            return phase
        
        global get_cpp_product
        _get_cpp_product_orig = get_cpp_product
        get_cpp_product = get_ctype_phase
    elif PPS_VALIDATION:
        cpp_filename = pps_files.cpp
            
    
    if len(calipso_filenames)==0:
        raise ValueError("Found no matching calipso files to: %r"% avhrr_file)
    for calipso_filename in calipso_filenames:       
        process_case(calipso_filename, avhrr_file, sunsat_filename, cpp_filename, **kwargs)
    
    if cloudtype and PPS_VALIDATION:
        # Restore get_cpp_product()
        get_cpp_product = _get_cpp_product_orig


def _filename_base(avhrr_filename, calipso_filename):
    return "match--%s--%s" % (os.path.basename(avhrr_filename),
                              os.path.basename(calipso_filename))

def get_mapper(avhrr_filename, calipso_filename, sunsat_filename):
    """
    Restore mapper object from file or create a new mapper and write it to file.
    
    Returns :class:`MatchMapper` instance mapping AVHRR to Calipso.
    
    """
    match_file = _filename_base(avhrr_filename, calipso_filename) + '.h5'
    match_path = os.path.join(MATCH_DIR, match_file)
    if os.path.exists(match_path):
        logger.info("Reading match from %r" % match_path)
        from amsr_avhrr.match import MatchMapper
        mapper = MatchMapper.from_file(match_path)
    else:
        mapper = match_with_calipso(calipso_filename, avhrr_filename,sunsat_filename,
                                    radius_of_influence=1e3,
                                    time_threshold=TIME_THR,n_neighbours=1)
        logger.info("Got mapper")  
        mapper.write(match_path, compression=_COMPRESSION)
        logger.info("Match written to %r" % match_path)
    
    return mapper


#def write_matched_values(filename, cpp_phase, cal_phase, cal_igbp):
#    import h5py
#    selected = (~cpp_phase.mask & ~cal_phase.mask)
#    with h5py.File(filename, 'w') as f:
#        f.create_dataset('cpp_phase', data=cpp_phase[selected],
#                         compression=_COMPRESSION)
#        f.create_dataset('calipso_phase', data=cal_phase[selected],
#                         compression=_COMPRESSION)
#        f.create_dataset('calipso_igbp', data=cal_igbp[selected],
#                         compression=_COMPRESSION)


def validate(cpp_phase, cal_phase, verbose=False):
    """
    Validate CPP (*cpp_phase*) and Calipso (*cal_phase*) cloud phase.
    Prints results to stdout.
    
    If *verbose* is True, print extended matrix.
    
    """
    import sys
    if hasattr(cpp_phase, 'mask'):
        n_pixels = (~(cpp_phase.mask | cal_phase.mask)).sum()
    else:
        n_pixels = cpp_phase.size
        print("Number of matching pixels: %d" % n_pixels)
    cpp_keys = ['liquid', 'ice', 'mixed', 'uncertain', 'non_opice', 'non_opwater']
    
    cal_name_len = 40
    cpp_name_len = 15
    print((' ' * (cal_name_len) +
           ''.join([k.rjust(cpp_name_len) for k in cpp_keys])).rjust(cal_name_len))
    
    def print_values(cal_ice_or_water, cpp_keys=cpp_keys):
        N = {}
        for k in cpp_keys:
            n = (cal_ice_or_water & (cpp_phase == CPP_PHASE_VALUES[k])).sum()
            N[k] = n
            sys.stdout.write("% 15d" % n)
        sys.stdout.write('\n')
        return N
    
    sys.stdout.write("Calipso water: ".rjust(cal_name_len))
    water_N = print_values(cal_phase == CALIPSO_PHASE_VALUES['water'])
    sys.stdout.write("Calipso ice: ".rjust(cal_name_len))
    ice_N = print_values((cal_phase == CALIPSO_PHASE_VALUES['ice']) |
                         (cal_phase == CALIPSO_PHASE_VALUES['horizontal_oriented_ice']))
    
    pod_water = 1. * water_N['liquid'] / (water_N['liquid'] + water_N['ice'])
    pod_ice = 1. * ice_N['ice'] / (ice_N['liquid'] + ice_N['ice'])
    far_water = 1. * ice_N['liquid'] / (water_N['liquid'] + ice_N['liquid'])
    far_ice = 1. * water_N['ice'] / (water_N['ice'] + ice_N['ice'])
    
    print('')
    print("POD: ".rjust(cal_name_len) + "% 15.2f% 15.2f" % (pod_water, pod_ice))
    print("FAR: ".rjust(cal_name_len) + "% 15.2f% 15.2f" % (far_water, far_ice))
    
    if verbose:
        print('')
        print((' ' * (cal_name_len) +
               ''.join([k.rjust(cpp_name_len) for k in CPP_PHASE_VALUES.keys()])).rjust(cal_name_len))
        for k in CALIPSO_PHASE_VALUES.keys():
            sys.stdout.write(("Calipso %s: " % k).rjust(cal_name_len))
            print_values(cal_phase == CALIPSO_PHASE_VALUES[k], CPP_PHASE_VALUES.keys())


def process_case(calipso_filename, avhrr_filename,sunsat_filename, cpp_filename=None,
                 ctype_filename=None, verbose=False,
                 max_layers=9, qual_min=CALIPSO_QUAL_VALUES['none']):
    """
    This is the work horse.
    
    """
    #mapper = get_mapper(avhrr_filename, calipso_filename, sunsat_filename)
    #if mapper goes here, program crashes because of some h5py bug. Nina 2013-06-01
    mapper = get_mapper(avhrr_filename, calipso_filename, sunsat_filename) 
    logger.info("Get caliop time:")
    calipso_time = get_calipso_time(calipso_filename)

    tim1 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(np.min(calipso_time)))
    tim2 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(np.max(calipso_time)))
    logger.info("Starttime caliop: %s, end time: %s"%(tim1, tim2))
    logger.info("Get caliop lonlat: ") 
    lon, lat = get_calipso_lonlat(calipso_filename)  

    #mapper = get_mapper(avhrr_filename, calipso_filename, sunsat_filename)    
    #mapper = get_mapper(avhrr_filename, calipso_filename, sunsat_filename) 
    logger.debug("Getting CPP water")
    if PPS_VALIDATION:
        cpp_phase = mapper(get_cpp_product(cpp_filename, 'cph'))
    if CCI_CLOUD_VALIDATION:
        cpp_cci_phase = get_cpp_product_cci(avhrr_filename, 'cph') 
        cpp_phase = mapper(cpp_cci_phase)
    cpp_phase = cpp_phase[..., 0] # mapper returns extra neighbours dimension
    #print cpp_phase
    logger.debug("Getting Calipso water")
    cal_phase = get_calipso_phase(calipso_filename, max_layers=max_layers,
                                  qual_min=qual_min)
    cal_igbp = get_calipso_igbp(calipso_filename)

    selection = ~cpp_phase.mask & ~cal_phase.mask

    if np.ma.isMaskedArray(selection):
        selection = selection.filled(False)
    restrictions = {'max number of CALIOP layers': str(max_layers),
                    'minimum quality': str(qual_min)}
 
    if PPS_VALIDATION:
        avhrr_time =  get_avhrr_time(sunsat_filename)
        cph_tmp = get_cpp_product(cpp_filename, 'cph')
    if CCI_CLOUD_VALIDATION:
        avhrr_lonlat,avhrr_time  = get_avhrr_lonlat_and_time_cci(avhrr_filename) 
        cph_tmp = get_cpp_product_cci(avhrr_filename, 'cph')
    #expand avhrr_time to the same shape as input avhrr_cph
    avhrr_time_expanded = np.ones(cph_tmp.shape)
    for i in range(cph_tmp.shape[0]):
        avhrr_time_expanded[i] = np.ones(cph_tmp.shape[1])*avhrr_time[i]


    avhrr_time_remapped = mapper(avhrr_time_expanded)
    #time_diff = avhrr_time_remapped - calipso_time
    #Why the above, why not the same:
    time_diff = mapper.time_diff
    time_diff = np.array(time_diff.ravel())


    if selection.any():
 
        filename = os.path.join(MATCH_DIR,_filename_base(avhrr_filename, calipso_filename) + '--values.h5')
        logger.info("Writing values to file %s"%filename)
        with h5py.File(filename, 'w') as f:
            write_data_to_open_file(np.array(cpp_phase[selection].ravel()), 'cpp_phase', f,
                                    attributes=restrictions)
            write_data_to_open_file(np.array(cal_phase[selection].ravel()), 'calipso_phase', f,
                       attributes=restrictions)
            write_data_to_open_file(np.array(cal_igbp[selection].ravel()), 'calipso_igbp', f,
                       attributes=restrictions)
            write_data_to_open_file(np.array(calipso_time[selection].ravel()), 'calipso_time', f)
            write_data_to_open_file(np.array(avhrr_time_remapped[selection].ravel()), 'avhrr_time', f,
                   attributes=restrictions)
            write_data_to_open_file(np.array(time_diff[selection].ravel()), 'time_diff', f)
            write_data_to_open_file(np.array(lon[selection].ravel()), 'longitudes', f)#, attributes=restrictions)
            write_data_to_open_file(np.array(lat[selection].ravel()), 'latitudes', f)#, attributes=restrictions)
            write_data_to_open_file(np.array(selection.ravel()), 'selection', f)#), attributes=restrictions)
        if f:
            f.close()
        #write_data(np.array(cpp_phase[selection].ravel()), 'cpp_phase', filename, mode='w',
        #           attributes=restrictions)
        #write_data(np.array(cal_phase[selection].ravel()), 'calipso_phase', filename,
        #           attributes=restrictions)
        #write_data(np.array(cal_igbp[selection].ravel()), 'calipso_igbp', filename,
        #           attributes=restrictions)
        #write_data(np.array(calipso_time[selection].ravel()), 'calipso_time', filename)
        #write_data(np.array(avhrr_time_remapped[selection].ravel()), 'avhrr_time', filename,
        #           attributes=restrictions)
        #write_data(np.array(time_diff[selection].ravel()), 'time_diff', filename)
        #write_data(np.array(lon[selection].ravel()), 'longitudes', filename)#, attributes=restrictions)
        #write_data(np.array(lat[selection].ravel()), 'latitudes', filename)#, attributes=restrictions)
        #write_data(np.array(selection.ravel()), 'selection', filename)#), attributes=restrictions)
    
        logger.info("Validating")
        validate(cpp_phase[selection], cal_phase[selection], verbose)
    else:
        logger.info("No matching points within time between %s and %s" %(
                calipso_filename, avhrr_filename))


#def get_frac_of_land(match_file):
#    """
#    Given a *match_file*, return mapped fraction of land in range [0, 1].
#    
#    """
#    from amsr_avhrr.match import MatchMapper
#    mapper = MatchMapper.from_file(match_file)
#    
#    from runutils import parse_scene
#    scene = os.path.basename(match_file).replace('match--', '')
#    satname, _datetime, orbit = parse_scene(scene)
#    
#    from pps_basic_configure import AUX_DIR
#    import pps_arguments 
#    phppsarg = pps_arguments.create_pps_argument(pps_arguments.USE_SUNSATANGLES_GLOB,
#                                               pps_arguments.SATELLITE_PROJ,
#                                               satname,orbit)
#    physiography_filename = os.path.join(AUX_DIR, ppsarg.files.physiography)
#    from epshdf import read_physiography
#    phys = read_physiography(physiography_filename, 0, 1, 0)
#    fraction_of_land = phys.fraction_of_land / 255.
#    
#    return mapper(fraction_of_land).ravel()
#
#
#def get_land_sea_selection(landsea=None, match_file=None, lonlat=None):
#    if not landsea:
#        return None
#    
#    if landsea not in ['land', 'sea']:
#        raise ValueError("*landsea* must be either 'land' or 'sea'")
#    
#    if lonlat is not None:
#        lon, lat = lonlat
#        from ppsFracOfLandOnSatellite import doFracOfLand
#        frac_of_land = doFracOfLand(lon, lat)
#        fol = frac_of_land.data * float(frac_of_land.info.info['gain'])
#    elif match_file:
#        fol = get_frac_of_land(match_file)
#    else:
#        raise ValueError("Either *match_file* or *lonlat* must be provided")
#    
#    if landsea == 'land':
#        return fol > .9
#    else:
#        return fol < .1
#
#
def get_land_sea_selection_cal(landsea, land_sea_cal):
        return None

def validate_all(matched_values_files, verbose=False, landsea=None):
    """
    Read all *matched_values_files* and perform validation on the concatenated
    arrays.
    
    If *landsea* is 'land' ('sea'), only land (sea) pixels will be used.
    
    """
    cpp_phases = []
    cal_phases = []
    lons = []
    lats = []
    for filename in matched_values_files:
        match_file = filename.replace('--values', '')
        #landsea_select = get_land_sea_selection(landsea, match_file=match_file)
        logger.info("File: %s"%filename)
        with h5py.File(filename, 'r') as f:
            land_sea_cal = f['calipso_igbp'][:]
            landsea_select = get_land_sea_selection_cal(landsea, land_sea_cal)
            if landsea is not None:
                selection = f['selection'][:]
                _slice = landsea_select#??
                #_slice = landsea_select[selection]
            else:
                _slice = slice(None)
            cpp_phases.append(f['cpp_phase'][_slice])
            cal_phases.append(f['calipso_phase'][_slice])
            lons.append(f['longitudes'][_slice])
            lats.append(f['latitudes'][_slice])
        if f:
            f.close()
    cpp_phase = np.concatenate(cpp_phases)
    cal_phase = np.concatenate(cal_phases)
    lon = np.concatenate(lons)
    lat = np.concatenate(lats)
    
    validate(cpp_phase, cal_phase, verbose)
    
    from amsr_avhrr.plotting import distribution_map
    fig = distribution_map(lon, lat)
    fig.suptitle("Distribution of valid pixels in cloud phase validation\n" +
                  "Number of Pixels: %d" % lon.size)
    fig.savefig('cph_distribution_all.pdf')


if __name__ == '__main__':
    import ConfigParser
    CONFIG_PATH = os.environ.get('ATRAINMATCH_CONFIG_DIR', './etc')
    CONF = ConfigParser.ConfigParser()
    CONF.read(os.path.join(CONFIG_PATH, "atrain_match.cfg"))
    OPTIONS = {}
    for option, value in CONF.items('general', raw = True):
        OPTIONS[option] = value

    from optparse import OptionParser
    parser = OptionParser(usage="Usage: %prog [options] satproj "
                                "<satname> <orbit> \n"
                                "or:    %prog [options] CASE [...]\n"
                                "          where CASE ~ 'possible/path/"
                                "noaa19_20110916_0959_03456*'")
    parser.add_option('-v', '--verbose', action='store_true')
    parser.add_option('-d', '--debug', action='store_true',
                      help="Don't ignore errors")
    parser.add_option('-l', '--max-layers', type='int', metavar='N',
                      help="Screen out pixels with more than N layers (default 1)")
    parser.add_option('-q', '--quality', type='int',
                      help="Screen out pixels with lower quality than N. "
                      "N must be in 0..3 (default 2)")
    parser.add_option('-t', '--cloudtype', action='store_true',
                      help="Get cloud phase from the cloud type product rather "
                      "than from CPP cloud phase")
    opts, args = parser.parse_args()
    
    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
        logger.debug("Verbose")
    else:
        logging.basicConfig(level=logging.WARNING)
    
    processing_kwargs = {}
    if opts.verbose is not None:
        processing_kwargs['verbose'] = opts.verbose
    if opts.max_layers is not None:
        processing_kwargs['max_layers'] = opts.max_layers
    if opts.quality is not None:
        processing_kwargs['qual_min'] = opts.quality
    if opts.cloudtype is not None:
        processing_kwargs['cloudtype'] = opts.cloudtype
    
    # Command line handling
    if args[0] == 'satproj':
        satname, orbit = args[1:]
        process_noaa_scene(filename, OPTIONS, **processing_kwargs)
    else:
        process_scenes(args, process_noaa_scene, OPTIONS, ignore_errors=not opts.debug,
                       **processing_kwargs)

