"""
TODO: The following description needs updating. atrain_match is now run via
process_master.py.

Program cloudsat_calipso_avhrr_match.py

This program is used to process and output statistics for the inter-comparison
of AVHRR PPS results and CloudSat/CALIPSO observations. It may be run
repeatedly and supervised by program process_master.py.

This particular version of Adam's original CloudSat/CALIPSO matchup and
analysis program has been complemented with the following:

 * Adjusted scales between CloudSat and CALIPSO datasets. The previous
   assumption that both datasets had 1 km resolution resulted in that datasets
   went out of phase for distances longer than about 1000 km. An empirical
   scale factor (CLOUDSAT_TRACK_RESOLUTION) of 1.076 is used to get the most
   optimal match.

 * AVHRR cloud top height datasets have been recalculated to heights above mean
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
   one of the following categories:

       ICE_COVER_SEA, ICE_FREE_SEA, SNOW_COVER_LAND, SNOW_FREE_LAND or COASTAL_ZONE

RUNNING INSTRUCTIONS
--------------------

The program is capable of running in a wide range of modes. 
 
Every exection of the program prints out statistics for PPS Cloud Mask, Cloud
Type and Cloud Top Height directly to file

Don't forget to set all parameters needed for standard ACPG/AHAMAP execution
since the matchup software uses parts of the ACPG/AHAMAP software.

Dependencies: For a successful run of the program the following supporting
              python modules must be available in the default run directory:

              cloudsat.py
              calipso.py
              cloudsat_calipso_avhrr_matchup.py
Updated 20150930 
Nina och KG


"""

#change log i found in git
import os
import sys

import numpy as np
from numpy import NaN
from config import (VAL_CPP,
                    NODATA,
                    PLOT_ONLY_PNG,
                    CCI_CLOUD_VALIDATION,
                    MAIA_CLOUD_VALIDATION,
                    PPS_VALIDATION,
                    CALIPSO_version4,
                    CALIPSO_version3,
                    ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS,
                    USE_5KM_FILES_TO_FILTER_CALIPSO_DATA,
                    ALSO_USE_1KM_FILES,
                    ALSO_USE_SINGLE_SHOT_CLOUD_CLEARED,
                    PPS_FORMAT_2012_OR_EARLIER,
                    MATCH_MODIS_LVL2, 
                    RESOLUTION)
import logging
logger = logging.getLogger(__name__)

import config
from common import attach_subdir_from_config, MatchupError

from cloudsat_calipso_avhrr_statistics import (CalculateStatistics)
from trajectory_plotting import plotSatelliteTrajectory
from cloudsat_calipso_avhrr_prepare import *

from read_cloudproducts_and_nwp_pps import NWPObj
from cloudsat import reshapeCloudsat, match_cloudsat_avhrr
from cloudsat import writeCloudsatAvhrrMatchObj, readCloudsatAvhrrMatchObj
from calipso import (reshapeCalipso, 
                     discardCalipsoFilesOutsideTimeRange,
                     match_calipso_avhrr, 
                     find_break_points, 
                     time_reshape_calipso)
from matchobject_io import (writeCaliopAvhrrMatchObj, 
                            readCaliopAvhrrMatchObj,
                            DataObject)
from calipso import  (detection_height_from_5km_data,
                      add1kmTo5km,
                      addSingleShotTo5km,
                      add5kmVariablesTo1kmresolution,
                      adjust5kmTo1kmresolution)
import inspect
import numpy as np
from cloudsat_calipso_avhrr_plot import (drawCalClsatAvhrrPlotTimeDiff,
                                         drawCalClsatGEOPROFAvhrrPlot, 
                                         drawCalClsatAvhrrPlotSATZ,
                                         drawCalClsatCWCAvhrrPlot)
#All non-avhrr satellites need to be here. Avhrr is default.
INSTRUMENT = {'npp': 'viirs',
              'noaa18': 'avhrr',
              'eos1': 'modis',
              'eos2': 'modis'} 
from datetime import datetime, timedelta
from glob import glob

class ppsFiles(object):
    def __init__(self, file_name_dict):
        self.cloudtype = None
        self.ctth = None
        self.cpp = None
        self.cma = None
        self.sunsatangles = None
        self.nwp_tsur = None
        self.nwp_t500 = None
        self.nwp_t700 = None
        self.nwp_t850 = None
        self.nwp_t950 = None
        self.nwp_ttro = None
        self.nwp_ciwv = None
        self.text_r06 = None
        self.text_t11 = None
        self.text_t37t12 = None
        self.text_t37 = None
        self.thr_t11ts_inv = None
        self.thr_t11t37_inv = None
        self.thr_t37t12_inv = None
        self.thr_t11t12_inv = None
        self.thr_t85t11_inv = None
        self.thr_t11ts = None
        self.thr_t11t37 = None
        self.thr_t37t12 = None
        self.thr_t11t12 = None
        self.thr_t85t11 = None
        self.thr_r06 = None
        self.thr_r09 = None
        self.nwp_segments = None
        self.__dict__.update(file_name_dict)    

test = 0




def get_time_list(cross_time, time_window, delta_t_in_seconds):
    tlist = []
    delta_t = timedelta(seconds=delta_t_in_seconds) #search per minute!delta_t_in_seconds=60
    tobj1 = cross_time
    tobj2 = cross_time - delta_t
    while (tobj1 <= cross_time + time_window[1] or 
           tobj2 >= cross_time - time_window[0]):
        if tobj1 <= cross_time + time_window[1]:
            tlist.append(tobj1)
            tobj1 = tobj1 + delta_t
        if   tobj2 >= cross_time - time_window[0]:   
            tlist.append(tobj2)
            tobj2 = tobj2 - delta_t  
    return tlist 

def find_calipso_files_inner(date_time, time_window, options, values):
    """Find the matching Calipso file"""
    tlist = get_time_list(date_time, time_window, 600)
    flist = []
    for tobj in tlist:    
        calipso_dir = insert_info_in_filename_or_path(
            options['calipso_dir'],
            values, datetime_obj=tobj)
        calipso_file_pattern = insert_info_in_filename_or_path(
            options['calipso_file'],
            values, 
            datetime_obj=tobj)
        tmplist = glob(os.path.join(calipso_dir, calipso_file_pattern))
        print "globbing", os.path.join(calipso_dir, calipso_file_pattern)
        flist.extend([ s for s in tmplist if s not in flist ])      
    return flist

def find_cloudsat_files_inner(date_time, time_window, options, values):
    """Find the matching Cloudsat file"""
    tlist = get_time_list(date_time, time_window, 600)
    flist = []
    for tobj in tlist:    
        cloudsat_dir = insert_info_in_filename_or_path(
            options['cloudsat_dir'],
            values, datetime_obj=tobj)
        cloudsat_file_pattern = insert_info_in_filename_or_path(
            options['cloudsat_file'],
            values, 
            datetime_obj=tobj)
        tmplist = glob(os.path.join(cloudsat_dir, cloudsat_file_pattern))
        flist.extend([ s for s in tmplist if s not in flist ])      
    return flist

def get_satid_datetime_orbit_from_fname(filename,as_oldstyle=False):
    from read_cloudproducts_and_nwp_pps import get_satid_datetime_orbit_from_fname_pps
    from read_cloudproducts_cci import get_satid_datetime_orbit_from_fname_cci
    from read_cloudproducts_maia import get_satid_datetime_orbit_from_fname_maia
    #Get satellite name, time, and orbit number from avhrr_file
    if PPS_VALIDATION:
        values = get_satid_datetime_orbit_from_fname_pps(filename, as_oldstyle=as_oldstyle)  
    if CCI_CLOUD_VALIDATION:
        values = get_satid_datetime_orbit_from_fname_cci(filename)
    if MAIA_CLOUD_VALIDATION:
        values = get_satid_datetime_orbit_from_fname_maia(filename)
    return values


def insert_info_in_filename_or_path(file_or_name_path, values, datetime_obj=None):
    #file_or_name_path = file_or_name_path.format(**values)
    satellite=values.get("satellite","*")
    file_or_name_path = file_or_name_path.format(
        satellite=satellite,
        orbit=values.get("orbit","*"),
        instrument = INSTRUMENT.get(satellite,"avhrr"),
        resolution=config.RESOLUTION,
        area=config.AREA,
        lines_lines=values.get("lines_lines", "*"),
        val_dir=config._validation_results_dir,
        year=values.get('year',"unknown"),
        month=values.get('month',"unknown"),
        mode=values.get('mode',"unknown"),
        #min_opt_depth=values.get('min_opt_depth',""),
        atrain_datatype=values.get("atrain_datatype","atrain_datatype"))

    if datetime_obj is None:
        return file_or_name_path
    name = datetime_obj.strftime(file_or_name_path)
    return name

def find_calipso_files(date_time, options, values):
    if config.CLOUDSAT_TYPE == 'CWC-RVOD':
        logger.info("Do not use CALIPSO in CWC-RWOD mode, Continue")
        calipso_files = []
    else:
        #might need to geth this in before looking for matchups
        tdelta_before = timedelta(seconds = (config.CALIPSO_FILE_LENGTH + 
                                             config.sec_timeThr))
        tdelta = timedelta(seconds = (config.SAT_ORBIT_DURATION + config.sec_timeThr))
        time_window = (tdelta_before, tdelta)
        calipso_files = find_calipso_files_inner(date_time, time_window, options, values)
        if len(calipso_files) > 1:
            logger.info("More than one Calipso file found within time window!")
        elif len(calipso_files) == 0:
            raise MatchupError("Couldn't find calipso matchup!")
        calipso_files = sorted(require_h5(calipso_files))
        calipso_basenames = [ os.path.basename(s) for s in calipso_files ]
        logger.info("Calipso files: " + str(calipso_basenames))
        return calipso_files

def find_cloudsat_files(date_time, options, values):
    #might need to geth this in before looking for matchups
    tdelta_before = timedelta(seconds = (config.CLOUDSAT_FILE_LENGTH + 
                                         config.sec_timeThr))
    tdelta = timedelta(seconds = (config.SAT_ORBIT_DURATION + config.sec_timeThr))
    time_window = (tdelta_before, tdelta)
    cloudsat_files = find_cloudsat_files_inner(date_time, time_window, options, values)
    if len(cloudsat_files) > 1:
        logger.info("More than one Cloudsat file found within time window!")
    elif len(cloudsat_files) == 0:
        logger.info("No Cloudsat file found within time window!")
        #raise MatchupError("Couldn't find cloudsat matchup!")
    cloudsat_files = sorted(require_h5(cloudsat_files))
    cloudsat_basenames = [ os.path.basename(s) for s in cloudsat_files ]
    logger.info("Cloudsat files: " + str(cloudsat_basenames))
    return cloudsat_files
   
# -----------------------------------------------------
def get_time_list_and_cross_time(cross):
    """
    Find the *satellite* avhrr file closest to *datetime*.
    """
    if cross.satellite1.lower() in ['cloudsat', 'calipso']:
        cross_satellite = cross.satellite2.lower()
        cross_time = cross.time2
    else:
        cross_satellite = cross.satellite1.lower()
        cross_time = cross.time1   
    ddt = timedelta(seconds=config.SAT_ORBIT_DURATION + config.sec_timeThr)
    ddt2 = ddt
    if ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS:
        ddt2=timedelta(seconds=0)
    if config.USE_ORBITS_THAT_STARTS_EXACTLY_AT_CROSS:
        ddt=timedelta(seconds=0)
        ddt2=timedelta(seconds=0)
    time_low = cross_time - ddt
    time_high = cross_time + ddt2
    logger.info("Searching for avhrr/viirs file with start time  between" 
              ": %s and %s  "%(time_low.strftime("%d %H:%M"),time_high.strftime("%d %H:%M"))) 
    time_window = (ddt, ddt2)
    # Make list of file times to search from:
    tlist = get_time_list(cross_time, time_window, delta_t_in_seconds=60)
    return tlist, cross_time, cross_satellite

def find_radiance_file(cross, options):
    found_file, tobj= find_avhrr_file(cross, 
                                      options['radiance_dir'], 
                                      options['radiance_file'])
    if not found_file:
        raise MatchupError("No dir or file found with radiance data!\n" + 
                           "Searching for %s %s" % (options['radiance_dir'],options['radiance_file']))
    return found_file, tobj
def find_cci_cloud_file(cross, options):
    found_file, tobj= find_avhrr_file(cross, 
                                      options['cci_dir'], 
                                      options['cci_file'])
    if not found_file:
        raise MatchupError("No dir or file found with cci cloud data!\n" + 
                           "Searching under %s" % options['cci_dir'])
    return found_file, tobj
def find_maia_cloud_file(cross, options):
    found_file, tobj= find_avhrr_file(cross, 
                                      options['maia_dir'], 
                                      options['maia_file'])
    if not found_file:
        raise MatchupError("No dir or file found with maia cloud data!\n" + 
                           "Searching under %s" % options['maia_dir'])
    return found_file, tobj


def find_avhrr_file(cross, filedir_pattern, filename_pattern, values={}):
    (tlist, cross_time, cross_satellite) = get_time_list_and_cross_time(cross)
    time_window=cross.time_window
    logger.info("Time window: %s", time_window)
    logger.info("Cross time: {cross_time}".format(cross_time=cross_time))
    values["satellite"] = cross_satellite
    found_dir = None
    no_files_found = True
    checked_dir = {}
    for tobj in tlist:
        #print values
        avhrr_dir = insert_info_in_filename_or_path(filedir_pattern,
                                                  values, datetime_obj=tobj)
        #print avhrr_dir
        if os.path.exists(avhrr_dir):
            found_dir = avhrr_dir
            if avhrr_dir not in checked_dir.keys():
                checked_dir[avhrr_dir] = 1
                logger.info("Found directory "
                          "{dirctory} ".format(dirctory=found_dir))
        else:
            if not found_dir and avhrr_dir not in checked_dir.keys():
                checked_dir[avhrr_dir] = 1
                logger.info("This directory does not exist, pattern:"
                          " {directory}".format(directory=avhrr_dir))
            continue

    
        file_pattern = insert_info_in_filename_or_path(filename_pattern,
                                                       values, datetime_obj=tobj)
        files = glob(os.path.join(found_dir, file_pattern))
        if len(files) > 0:
            no_files_found = False
            logger.info("Found files: " + os.path.basename(str(files[0])))
            return files[0], tobj
    if not found_dir:       
        return None, None
    if no_files_found:       
         logger.info( "Found no files for patterns of type:"
                   " {pattern}".format(pattern=file_pattern)  )
    return None, None


def require_h5(files):
    """
    Convert any '.hdf' files in *files* to '.h5'. Returns a list of '.h5' files.
    
    """
    from config import H4H5_EXECUTABLE #@UnresolvedImport
    h5_files = []
    for f in files:
        if f.endswith('.h5'):
            h5_files.append(f)
        elif f.endswith('.hdf'):
            h5_file = f.replace('.hdf', '.h5')
            if not h5_file in files:
                print(f)
                #sys.exit()
                command = '%s %s' % (H4H5_EXECUTABLE, f)
                logger.info("Converting %r to HDF5" % f)
                if os.system(command):
                    raise RuntimeError("Couldn't convert %r to HDF5" % f)
                h5_files.append(h5_file)
        else:
            raise ValueError("File format of %r not recognized" % f)
    
    return h5_files

def get_pps_file(avhrr_file, options, values, type_of_file, file_dir):
    if not type_of_file in options:              
        logger.info("No %s file in cfg-file."%(type_of_file))
        return None
    date_time = values["date_time"]
    cloudtype_name = insert_info_in_filename_or_path(options[type_of_file], 
                                                     values, datetime_obj=date_time)                        
    path = insert_info_in_filename_or_path(options[file_dir], 
                                           values, datetime_obj=date_time)  
    try:
        file_name = glob(os.path.join(path, cloudtype_name))[0]
        #logger.debug( type_of_file +": ", file_name)
        return file_name
    except IndexError:
        logger.info("No %s file found corresponding to %s." %( type_of_file, avhrr_file))
        return None

def find_files_from_avhrr(avhrr_file, options, as_oldstyle=False):
    """
    Find all files needed to process matchup from source data files.
    """
    # Let's get the satellite and the date-time of the pps radiance
    # (avhrr/viirs) file:
    logger.info("IMAGER: %s" % avhrr_file)
    values = get_satid_datetime_orbit_from_fname(avhrr_file,
                                                 as_oldstyle=as_oldstyle)
    date_time = values["date_time"]

    cloudtype_name = insert_info_in_filename_or_path(options['cloudtype_file'], 
                                                     values, datetime_obj=date_time)
    path = insert_info_in_filename_or_path(options['cloudtype_dir'], 
                                           values, datetime_obj=date_time)
    try:
        #print os.path.join(path, cloudtype_name)
        cloudtype_file = glob(os.path.join(path, cloudtype_name))[0]
    except IndexError:
        raise MatchupError("No cloudtype file found corresponding to %s." % avhrr_file)
    logger.info("CLOUDTYPE: " + cloudtype_file)
    cma_name = insert_info_in_filename_or_path(options['cma_file'], 
                                                     values, datetime_obj=date_time)
    path = insert_info_in_filename_or_path(options['cma_dir'], 
                                           values, datetime_obj=date_time)
    try:
        #print os.path.join(path, cma_name)
        cma_file = glob(os.path.join(path, cma_name))[0]
    except IndexError:
        raise MatchupError("No cma file found corresponding to %s." % avhrr_file)
    logger.info("CMA: " + cma_file)
    ctth_name = insert_info_in_filename_or_path(options['ctth_file'], 
                                                values, datetime_obj=date_time)
    path =  insert_info_in_filename_or_path(options['ctth_dir'], 
                                                values, datetime_obj=date_time)  
    try:
        ctth_file = glob(os.path.join(path, ctth_name))[0]
    except IndexError:
        raise MatchupError("No ctth file found corresponding to %s." % avhrr_file)
    logger.info("CTTH: " + ctth_file)
    if VAL_CPP: 
        cpp_name = insert_info_in_filename_or_path(options['cpp_file'],
                                                   values, datetime_obj=date_time)
        path = insert_info_in_filename_or_path(options['cpp_dir'],
                                               values, datetime_obj=date_time)
        try:
            cpp_file = glob(os.path.join(path, cpp_name))[0]
        except IndexError:
            raise MatchupError("No cpp file found corresponding to %s." % avhrr_file)
        logger.info("CPP: " + cpp_file)
    else:     
       logger.info("Not validation of CPP ")
       cpp_file = None
    sunsatangles_name = insert_info_in_filename_or_path(options['sunsatangles_file'],
                                                        values, datetime_obj=date_time)
    path =insert_info_in_filename_or_path(options['sunsatangles_dir'], 
                                          values, datetime_obj=date_time)
    try:
        sunsatangles_file = glob(os.path.join(path, sunsatangles_name))[0]
    except IndexError:
        sunsatangles_file = None
        raise MatchupError("No sunsatangles file found corresponding to %s." % avhrr_file)
    logger.info("SUNSATANGLES: " + sunsatangles_file)

    if not 'physiography_file' in options:
        logger.warning("No physiography file searched for!")
        physiography_file = None
    else:
        physiography_name = insert_info_in_filename_or_path(options['physiography_file'],
                                                            values, datetime_obj=date_time)
        path =insert_info_in_filename_or_path(options['physiography_dir'], 
                                              values, datetime_obj=date_time)
        pattern = os.path.join(path, physiography_name)
        try:
            physiography_file = glob(pattern)[0]
        except IndexError:
            physiography_file = None
            raise MatchupError("No physiography file found corresponding to %s." % pattern)
        logger.info("PHYSIOGRAPHY: " + physiography_file)

    if not 'nwp_tsur_file' in options:
        logger.warning("No surface temperature file searched for!")
        nwp_tsur_file=None
    else:
        nwp_tsur_name = insert_info_in_filename_or_path(options['nwp_tsur_file'],
                                                        values, datetime_obj=date_time)
        path = insert_info_in_filename_or_path(options['nwp_nwp_dir'],
                                               values, datetime_obj=date_time)
        pattern = os.path.join(path, nwp_tsur_name)
        try:
            nwp_tsur_file = glob(pattern)[0]
        except IndexError:
            raise MatchupError("No nwp_tsur file found corresponding to %s." % pattern)
        logger.info("NWP_TSUR: " + nwp_tsur_file)

    file_name_dict={}
    for nwp_file in ['nwp_tsur','nwp_t500','nwp_t700',
                     'nwp_t850','nwp_t950', 'nwp_ciwv', 'nwp_ttro']:   
        file_name_dict[nwp_file] = get_pps_file(avhrr_file, options, values, 
                                                 nwp_file+'_file', 'nwp_nwp_dir')

    emis_file = get_pps_file(avhrr_file, options, values, 
                                           'emis_file', 'emis_dir')
    file_name_dict['emis'] = emis_file
    seaice_file = get_pps_file(avhrr_file, options, values, 
                               'seaice_file', 'seaice_dir')
    file_name_dict['seaice'] = seaice_file

    for text_file in ['text_r06', 'text_t11', 'text_t37t12', 'text_t37']:
        file_name_dict[text_file] = get_pps_file(avhrr_file, options, values, 
                                                 text_file+'_file', 'text_dir')

    for thr_file in ['thr_t11ts_inv', 'thr_t11t37_inv', 
                     'thr_t37t12_inv', 'thr_t11t12_inv', 
                     'thr_t11ts', 'thr_t11t37', 'thr_t37t12', 'thr_t11t12',
                     'thr_r09', 'thr_r06']:
        file_name_dict[thr_file] = get_pps_file(avhrr_file, options, values, 
                                                thr_file+'_file', 'thr_dir')
    if (values['satellite'] in ('npp', 'eos1', 'eos2')):
        for thr_file in ['thr_t85t11_inv', 'thr_t85t11']:
            file_name_dict[thr_file] = get_pps_file(avhrr_file, options,
                                                    values, 
                                                    thr_file+'_file', 'thr_dir')

    file_name_dict['nwp_segments'] = get_pps_file(avhrr_file, options,
                                                  values, 
                                                  'segment_file', 'segment_dir') 
    file_name_dict.update({'cloudtype': cloudtype_file,
                           'cma': cma_file,
                           'ctth': ctth_file,
                           'cpp': cpp_file,
                           'nwp_tsur': nwp_tsur_file,
                           'sunsatangles': sunsatangles_file,
                           'physiography': physiography_file})

    ppsfiles = ppsFiles(file_name_dict)
    return  ppsfiles

def get_cloudsat_matchups(cloudsat_files, cloudtype_file, avhrrGeoObj, avhrrObj,
                          ctype, cma, ctth, nwp_obj, avhrrAngObj, cpp, config_options):
    """
    Read Cloudsat data and match with the given PPS data.
    """
    if config.CLOUDSAT_TYPE == 'CWC-RVOD':
        if config.RESOLUTION == 1:      
            if test == 1:
                from cloudsat_cwc import match_cloudsatCwc_avhrr
                from cloudsat_cwc import reshapeCloudsat1kmCwc
                reshape_fun = reshapeCloudsat1kmCwc
                match_fun = match_cloudsatCwc_avhrr
            else:
                reshape_fun = reshapeCloudsat
                match_fun = match_cloudsat_avhrr 
        elif config.RESOLUTION == 5:
            from cloudsat5km_cwc import match_cloudsat5kmCwc_avhrr5km
            from cloudsat5km_cwc import reshapeCloudsat5kmCwc
            reshape_fun = reshapeCloudsat5kmCwc
            match_fun = match_cloudsat5kmCwc_avhrr5km
    else:
        reshape_fun = reshapeCloudsat
        match_fun = match_cloudsat_avhrr
    cloudsat = reshape_fun(cloudsat_files, avhrrGeoObj, cloudtype_file)
    cl_matchup = match_fun(cloudtype_file, cloudsat,
                           avhrrGeoObj, avhrrObj, ctype, cma,
                           ctth, nwp_obj, avhrrAngObj, cpp)
    return cl_matchup

def total_and_top_layer_optical_depth_5km(calipso, resolution=5):
    logger.info("Find total optical depth from 5km data")
    optical_depth_in = calipso.feature_optical_depth_532
    o_depth_top_layer = -9.0 + 0*calipso.number_layers_found.ravel()
    total_o_depth = -9.0 + 0*calipso.number_layers_found.ravel()
    if resolution==5:
        pixels = np.logical_and(
            calipso.number_layers_found.ravel()>0,
            optical_depth_in[:,0].ravel() >= 0)   
        o_depth_top_layer[pixels] = optical_depth_in[pixels, 0]
        total_o_depth[pixels] =  optical_depth_in[pixels, 0]       
        for lay in range(1, np.max(calipso.number_layers_found[pixels]), 1):  
            pixels = np.logical_and(
                pixels, 
                optical_depth_in[:, lay]>=0)
            total_o_depth[pixels] +=  optical_depth_in[pixels, lay]
    else:
        print "ERROR this fuction is only for 5km data!"
        print "These features can then added to 1km data set"
    calipso.feature_optical_depth_532_top_layer_5km = o_depth_top_layer
    calipso.total_optical_depth_5km = total_o_depth       
  
    return calipso 

 

def get_calipso_matchups(calipso_files, values, 
                         avhrrGeoObj, avhrrObj, ctype, cma,  ctth, 
                         nwp_obj, avhrrAngObj, options, cpp=None, 
                         nwp_segments=None, cafiles1km=None, cafiles5km=None, 
                         cafiles5km_aerosol=None):
    """
    Read Calipso data and match with the given PPS data.
    """
    surftemp = nwp_obj.surftemp
    #First remove files clearly outside time limit from the lists
    #When combinating 5km and 1km data some expensive calculations are done
    #before cutting the data that fits the time condition.
    #It is unnessecary to do this for files where no-pixel will match!
    calipso_files = discardCalipsoFilesOutsideTimeRange(
        calipso_files, avhrrGeoObj, values)
    if cafiles1km is not None:
        cafiles1km = discardCalipsoFilesOutsideTimeRange(
            cafiles1km, avhrrGeoObj, values, res=1)
    if cafiles5km is not None:
        cafiles5km = discardCalipsoFilesOutsideTimeRange(
            cafiles5km, avhrrGeoObj, values, res=5)
    if cafiles5km_aerosol is not None:
        cafiles5km_aerosol = discardCalipsoFilesOutsideTimeRange(
            cafiles5km_aerosol, avhrrGeoObj, values, res=5, ALAY=True)
        
    calipso  = reshapeCalipso(calipso_files)
    #find time breakpoints, but don't cut the data yet ..
    startBreak, endBreak = find_break_points(calipso,  avhrrGeoObj, values)
    if cafiles1km is not None and CALIPSO_version3:
        #RESOLUTION 5km also have 1km data
        logger.info("Calipso version 3 data used and old 1 km restore method!")
        calipso1km = reshapeCalipso(cafiles1km, res=1)
        calipso5km = calipso
        calipso = add1kmTo5km(calipso1km, calipso5km)
        calipso = time_reshape_calipso(calipso, 
                                       start_break=startBreak, end_break=endBreak) 
        calipso = total_and_top_layer_optical_depth_5km(calipso, resolution=5)

    elif cafiles5km is not None and CALIPSO_version4:
        #RESOLUTION 1km also have 5km data calipso version 4
        logger.info("Calipso version 4, single shot fraction and "
                    "old 5km restored optical depth method used!")
        calipso1km = calipso
        calipso5km = reshapeCalipso(cafiles5km, res=5)
        calipso5km = addSingleShotTo5km(calipso5km) 
        calipso5km = total_and_top_layer_optical_depth_5km(calipso5km, resolution=5)
        calipso1km = add5kmVariablesTo1kmresolution(calipso1km, calipso5km)
        if USE_5KM_FILES_TO_FILTER_CALIPSO_DATA:
            logger.info("Find detection height using 5km data")
            calipso1km = detection_height_from_5km_data(calipso1km, calipso5km)
        calipso = time_reshape_calipso(calipso1km,  
                                       start_break=startBreak, 
                                       end_break=endBreak) 

    elif cafiles5km is not None and CALIPSO_version3:
        #RESOLUTION 1km also have 5km data calipso version 3
        logger.info("Calipso version 3 data used and old 5 km restored optical depth method!")
        calipso1km = calipso
        calipso5km = reshapeCalipso(cafiles5km, res=5)
        calipso5km = total_and_top_layer_optical_depth_5km(calipso5km, resolution=5)
        calipso1km = add5kmVariablesTo1kmresolution(calipso1km, calipso5km)
        if USE_5KM_FILES_TO_FILTER_CALIPSO_DATA:
            logger.info("Find detection height using 5km data")
            calipso1km = detection_height_from_5km_data(calipso1km, calipso5km)
        calipso = time_reshape_calipso(calipso1km,  
                                       start_break=startBreak, 
                                       end_break=endBreak) 
    elif CALIPSO_version4 and RESOLUTION == 5 and ALSO_USE_SINGLE_SHOT_CLOUD_CLEARED:

        #RESOLUTION exclusively 5km data but additional clouds taken from 330 m single shot resolution
        logger.info("Calipso version 4 data used and new single shot restore method!")
        #calipso5km = reshapeCalipso(cafiles5km, res=5)
        calipso5km  = reshapeCalipso(calipso_files)
        calipso = addSingleShotTo5km(calipso5km) 
        calipso = time_reshape_calipso(calipso,  
                                       start_break=startBreak, 
                                       end_break=endBreak) 
        calipso = total_and_top_layer_optical_depth_5km(calipso, resolution=5)
    elif CALIPSO_version4 and RESOLUTION == 5 and ALSO_USE_1KM_FILES:

        #RESOLUTION exclusively 5km data but additional clouds taken from 1 km data
        logger.info("Calipso version 4 data used but old method combining 1 km and 5 km data!")
        #calipso5km = reshapeCalipso(cafiles5km, res=5)
        calipso5km  = reshapeCalipso(calipso_files)
        calipso1km = reshapeCalipso(cafiles1km, res=1)
        calipso = add1kmTo5km(calipso1km, calipso5km)
        calipso = time_reshape_calipso(calipso,  
                                       start_break=startBreak, 
                                       end_break=endBreak) 
        calipso = total_and_top_layer_optical_depth_5km(calipso, resolution=5)
    else:
        logger.info("Standard method used (no additional resolutions or layers)!")
        calipso = time_reshape_calipso(calipso,  
                                       start_break=startBreak, 
                                       end_break=endBreak)
        if RESOLUTION == 5:
            calipso = total_and_top_layer_optical_depth_5km(calipso, resolution=5)

    #aerosol-data
    calipso_aerosol = None
    if cafiles5km_aerosol is not None:
        calipso5km_aerosol = reshapeCalipso(cafiles5km_aerosol, res=5, ALAY=True)
        if RESOLUTION == 1:
            calipso_aerosol = adjust5kmTo1kmresolution(calipso5km_aerosol)
        elif RESOLUTION == 5:
            calipso_aerosol = calipso5km_aerosol
        calipso_aerosol = time_reshape_calipso(calipso_aerosol,  
                                               start_break=startBreak, 
                                               end_break=endBreak)
    # free some memory    
    calipso1km = None
    calipso5km = None
        
    logger.info("Matching with avhrr")
    ca_matchup = match_calipso_avhrr(values, calipso, calipso_aerosol,
                              avhrrGeoObj, avhrrObj, ctype, cma,
                              ctth, cpp, nwp_obj, avhrrAngObj, 
                              nwp_segments, options)
    return ca_matchup
def read_cloud_cci(avhrr_file):
    from read_cloudproducts_cci import cci_read_all
    return cci_read_all(avhrr_file)

def read_cloud_maia(avhrr_file):
    from read_cloudproducts_maia import maia_read_all
    return maia_read_all(avhrr_file)

def read_pps_data(pps_files, avhrr_file, cross):
    from read_cloudproducts_and_nwp_pps import pps_read_all
    return pps_read_all(pps_files, avhrr_file, cross)

def get_additional_calipso_files_if_requested(calipso_files):
    import glob
    calipso5km = None
    calipso1km = None
    calipso5km_aerosol=None

    if config.RESOLUTION == 5:
        if config.ALSO_USE_1KM_FILES == True:
            logger.info("Search for CALIPSO 1km data too")
            calipso1km = []
            calipso5km = []
            for file5km in calipso_files:
                file1km = file5km.replace('/5km/', '/1km/').\
                          replace('05kmCLay', '01kmCLay').\
                          replace('-Prov-V3-01.', '*').\
                          replace('-Prov-V3-02.', '*').\
                          replace('-Prov-V3-30.', '*').\
                          replace('.hdf', '.h5')
                files_found = glob.glob(file1km)
                if len(files_found)==0:
                    #didn't find h5 file, might be hdf file instead
                    file1km = file1km.replace('.h5','.hdf')
                    files_found = glob.glob(file1km)
                if len(files_found)>0:    
                    calipso1km.append(files_found[0])
                    calipso5km.append(file5km)
            calipso1km = sorted(require_h5(calipso1km))
            calipso_files = sorted(calipso5km)
            calipso5km = None

            if len(calipso_files) == 0:
                raise MatchupError("Did not find any matching 5km and 1km calipso files")

            if len(calipso_files) != len(calipso1km):
                raise MatchupError("Inconsistent number of calipso files...\n" + 
                                   "\tlen(calipso_files) = %d\n" % len(calipso_files) + 
                                   "\tlen(calipso1km) = %d" % len(calipso1km))
    if config.RESOLUTION == 1:
        if config.ALSO_USE_5KM_FILES == True:
            logger.info("Search for CALIPSO 5km data too")
            calipso5km = []
            calipso1km = []
            for file1km in calipso_files:
                file5km = file1km.replace('/1km/', '/5km/').\
                          replace('01kmCLay', '05kmCLay').\
                          replace('-ValStage1-V3-30.', '*').\
                          replace('-ValStage1-V3-01.', '*').\
                          replace('-ValStage1-V3-02.', '*')
                files_found = glob.glob(file5km)
                if len(files_found)==0:
                    #didn't find h5 file, might be hdf file instead
                    file5km = file5km.replace('.h5','.hdf')
                    files_found = glob.glob(file5km)
                if len(files_found)>0: 
                    calipso5km.append(files_found[0])
                    calipso1km.append(file1km)
            calipso5km = sorted(require_h5(calipso5km))
            calipso_files = sorted(calipso1km)
            calipso1km = None
            if len(calipso_files) == 0:
                raise MatchupError("Did not find any matching 5km and 1km calipso files")
            if len(calipso_files) != len(calipso5km):
                raise MatchupError("Inconsistent number of calipso files...\n" + 
                                   "\tlen(calipso_files) = %d\n" % len(calipso_files) + 
                                   "\tlen(calipso1km) = %d" % len(calipso5km))
    if config.MATCH_AEROSOL_CALIPSO:
        calipso5km_aerosol=[]
        for cfile in calipso_files:
            file5km_aerosol = cfile.replace('/CLAY/', '/ALAY/').\
                              replace('CLay', 'ALay').\
                              replace('/1km/', '/5km/').\
                              replace('01km', '05km').\
                              replace('-ValStage1-V3-30.', '*').\
                              replace('-ValStage1-V3-01.', '*').\
                              replace('-ValStage1-V3-02.', '*').\
                              replace('-Prov-V3-01.', '*').\
                              replace('-Prov-V3-02.', '*').\
                              replace('-Prov-V3-30.', '*')
            files_found_aerosol = glob.glob(file5km_aerosol)
            if len(files_found_aerosol)==0:
                #didn't find h5 file, might be hdf file instead
                file5km_aerosol = file5km_aerosol.replace('.h5','.hdf')
                files_found_aerosol = glob.glob(file5km_aerosol)
            if len(files_found_aerosol)>0: 
                calipso5km_aerosol.append(files_found_aerosol[0]) 
        print "found these aerosol files", calipso5km_aerosol             
        calipso5km_aerosol = sorted(require_h5(calipso5km_aerosol))
        print "found these aerosol files", calipso5km_aerosol
    return calipso5km, calipso1km, calipso5km_aerosol              


def add_additional_clousat_calipso_index_vars(clsatObj, caObj):
    #add cloudsat modisflag to calipso obj
    caObj.calipso.cal_MODIS_cflag = None
    if clsatObj is not None:
        #map cloudsat to calipso and the other way around!
        from amsr_avhrr.match import match_lonlat
        source = (clsatObj.cloudsat.longitude.astype(np.float64).reshape(-1,1), clsatObj.cloudsat.latitude.astype(np.float64).reshape(-1,1))
        target = (caObj.calipso.longitude.astype(np.float64).reshape(-1,1), caObj.calipso.latitude.astype(np.float64).reshape(-1,1))
        mapper = match_lonlat(source, target, radius_of_influence=1000, n_neighbours=1)
        caObj.calipso.cloudsat_index = mapper.rows.filled(NODATA).ravel()
        target = (clsatObj.cloudsat.longitude.astype(np.float64).reshape(-1,1), clsatObj.cloudsat.latitude.astype(np.float64).reshape(-1,1))
        source = (caObj.calipso.longitude.astype(np.float64).reshape(-1,1), caObj.calipso.latitude.astype(np.float64).reshape(-1,1))
        mapper = match_lonlat(source, target, radius_of_influence=1000, n_neighbours=1)
        clsatObj.cloudsat.calipso_index = mapper.rows.filled(NODATA).ravel()

        # Transfer CloudSat MODIS cloud flag to CALIPSO representation
        caObj.calipso.cal_MODIS_cflag = np.zeros(len(caObj.calipso.elevation),'b')
        caObj.calipso.cal_MODIS_cflag[caObj.calipso.cloudsat_index] = clsatObj.cloudsat.MODIS_cloud_flag
    return clsatObj, caObj

def add_elevation_corrected_imager_ctth(clsatObj, caObj):
    ## Cloudsat ##
    clsatObj.avhrr.avhrr_ctth_csat_ok = None
    if clsatObj is not None:
        # First make sure that PPS cloud top heights are converted to height
        # above sea level just as CloudSat height are defined. Use
        # corresponding DEM data.
        elevation = np.where(np.less_equal(clsatObj.cloudsat.elevation,0),
                             0,clsatObj.cloudsat.elevation)		
        num_csat_data_ok = len(clsatObj.cloudsat.elevation)
        logger.info("Length of CLOUDSAT array: %d", num_csat_data_ok )
        avhrr_ctth_csat_ok = np.array(clsatObj.avhrr.ctth_height).copy().ravel()
        if CCI_CLOUD_VALIDATION: 
            #ctth already relative mean sea level
            pass
        else: #ctth relative topography
            got_height = avhrr_ctth_csat_ok>=0                    
            avhrr_ctth_csat_ok[got_height] += elevation[got_height]*1.0
        clsatObj.avhrr.avhrr_ctth_csat_ok = avhrr_ctth_csat_ok
        if num_csat_data_ok == 0:
            logger.info("Processing stopped: Zero lenght of matching arrays!")
            print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit()
    ## Calipso ##        
    # First make sure that PPS cloud top heights are converted to height above sea level
    # just as CALIPSO heights are defined. Use corresponding DEM data.
    cal_elevation = np.where(np.less_equal(caObj.calipso.elevation,0),
                             0,caObj.calipso.elevation)
    num_cal_data_ok = len(caObj.calipso.elevation)
    logger.info("Length of CALIOP array: %d", num_cal_data_ok)
    avhrr_ctth_cal_ok = np.array(caObj.avhrr.ctth_height).copy().ravel()
    if CCI_CLOUD_VALIDATION: 
        #ctth relative mean sea level
        pass
    else: #ctth relative topography
        got_height = avhrr_ctth_cal_ok>=0                    
        avhrr_ctth_cal_ok[got_height] += cal_elevation[got_height]*1.0
    caObj.avhrr.avhrr_ctth_cal_ok = avhrr_ctth_cal_ok
    return clsatObj, caObj


def get_matchups_from_data(cross, config_options):
    """
    Retrieve Cloudsat- and Calipso-AVHRR matchups from Cloudsat, Calipso, and
    PPS files.
    """
    import os #@Reimport
    if (PPS_VALIDATION ):
        avhrr_file, tobj = find_radiance_file(cross, config_options)
        values = get_satid_datetime_orbit_from_fname(avhrr_file)
        if not avhrr_file:
            raise MatchupError("No avhrr file found!\ncross = " + str(cross))
        pps_files = find_files_from_avhrr(avhrr_file, config_options)   
        retv =read_pps_data(pps_files, avhrr_file, cross)
        (avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, 
         nwp_obj, cpp, nwp_segment, cma )= retv
        date_time = values["date_time"]
    if (CCI_CLOUD_VALIDATION):
        avhrr_file, tobj = find_cci_cloud_file(cross, config_options)
        #avhrr_file = "20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
        values = get_satid_datetime_orbit_from_fname(avhrr_file)
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surftemp, cpp =read_cloud_cci(avhrr_file)
        nwp_segment = None
        nwp_obj = NWPObj({'surftemp':surftemp})        
        avhrrGeoObj.satellite = values["satellite"];
        date_time = values["date_time"]
    if (MAIA_CLOUD_VALIDATION):
        avhrr_file, tobj = find_maia_cloud_file(cross, config_options)
        #viiCT_npp_DB_20120817_S035411_E035535_DES_N_La052_Lo-027_00001.h5
        values = get_satid_datetime_orbit_from_fname(avhrr_file)
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surftemp, cpp, cma =read_cloud_maia(avhrr_file)
        nwp_segment = None
        nwp_obj = NWPObj({'surftemp':surftemp})        
        avhrrGeoObj.satellite = values["satellite"];
        date_time = values["date_time"]

    #CLOUDSAT:  
    cl_matchup = None
    if (PPS_VALIDATION):
        cloudsat_files = find_cloudsat_files(date_time, config_options, values)
        print cloudsat_files
        if (isinstance(cloudsat_files, str) == True or 
            (isinstance(cloudsat_files, list) and len(cloudsat_files) != 0)):
            logger.info("Read CLOUDSAT %s data" % config.CLOUDSAT_TYPE)
            cl_matchup = get_cloudsat_matchups(cloudsat_files, 
                                               pps_files.cloudtype,
                                               avhrrGeoObj, avhrrObj, ctype, cma,
                                               ctth, nwp_obj, avhrrAngObj, cpp, config_options)
        else:
            logger.info("NO CLOUDSAT File, Continue")
    else:
        logger.info("NO CLOUDSAT File,"
                  "CCI-cloud validation only for calipso, Continue")
    #CALIPSO:
    ca_matchup = None
    calipso_files = find_calipso_files(date_time, config_options, values)
    if (isinstance(calipso_files, str) == True or 
        (isinstance(calipso_files, list) and len(calipso_files) != 0)):
        extra_files = get_additional_calipso_files_if_requested(calipso_files)
        calipso5km, calipso1km, calipso5km_aerosol  = extra_files
      
        logger.info("Read CALIPSO data")        
        ca_matchup= get_calipso_matchups(calipso_files, 
                                         values,
                                         avhrrGeoObj, avhrrObj, 
                                         ctype, cma, ctth, 
                                         nwp_obj, avhrrAngObj, 
                                         config_options, cpp,
                                         nwp_segment,
                                         calipso1km, calipso5km, calipso5km_aerosol)
    else:
        logger.info("NO CALIPSO File, Continue")


    # Get satellite name, time, and orbit number from avhrr_file
    date_time = values["date_time"]
    #basename = '_'.join(os.path.basename(avhrr_file).split('_')[:4])
    values = get_satid_datetime_orbit_from_fname(avhrr_file)
    basename = values["basename"]
    rematched_path = date_time.strftime(config_options['reshape_dir'].format(
            val_dir=config._validation_results_dir,
            satellite=values["satellite"],
            resolution=str(config.RESOLUTION),
            area=config.AREA,
            ))
    rematched_file = date_time.strftime(config_options['reshape_file'].format(
            satellite=values["satellite"],
            orbit=values["orbit"],
            resolution=str(config.RESOLUTION),
            instrument=INSTRUMENT.get(values["satellite"],'avhrr'),
            atrain_datatype="atrain_datatype"
            ))
    rematched_file_base = rematched_path + rematched_file
    
    # Create directories if they don't exist yet
    if not os.path.exists(os.path.dirname(rematched_path)):
        logger.info("Creating dir %s:"%(rematched_path))
        os.makedirs(os.path.dirname(rematched_path))

    #add modis lvl2    
    if config.MATCH_MODIS_LVL2 and config.IMAGER_INSTRUMENT.lower() in ['modis']:
        from read_modis_products import add_modis_06    
        ca_matchup = add_modis_06(ca_matchup, avhrr_file, config_options) 
    #add additional vars to cloudsat and calipso objects:
    cl_matchup, ca_matchup = add_additional_clousat_calipso_index_vars(cl_matchup, ca_matchup)

    # Write cloudsat matchup    
    if cl_matchup is not None:
        cl_match_file = rematched_file_base.replace(
            'atrain_datatype', 'cloudsat-%s' % config.CLOUDSAT_TYPE)
        writeCloudsatAvhrrMatchObj(cl_match_file, cl_matchup)
    else:
        logger.info('CloudSat is not defined. No CloudSat Match File created')
     
    # Write calipso matchup
    if ca_matchup is not None:
        ca_match_file = rematched_file_base.replace('atrain_datatype', 'caliop')
        avhrr_obj_name = 'pps'
        if config.CCI_CLOUD_VALIDATION:
            avhrr_obj_name = 'cci'
        if config.MAIA_CLOUD_VALIDATION:
            avhrr_obj_name = 'maia'
        writeCaliopAvhrrMatchObj(ca_match_file, ca_matchup, 
                                 avhrr_obj_name = avhrr_obj_name) 
    
    return {'cloudsat': cl_matchup, 'calipso': ca_matchup, 
            'basename': basename, 'values':values}


def get_matchups(cross, options, reprocess=False):
    """
    Retrieve Cloudsat- and Calipso-AVHRR matchups. If *reprocess* is False, and if
    matchup files exist, get matchups directly from the processed files.
    Otherwise process Cloudsat, Calipso, and PPS files first.
    """
    caObj = None
    clObj = None
    values = {}
    try:
        values["satellite"] = cross.satellite1.lower()
        date_time_cross = cross.time1
    except AttributeError:
        raise ValueError('cross is not a valid SNO cross. (cross: %s)' % cross)
    
    if values["satellite"] in ['calipso', 'cloudsat']:
        values["satellite"] = cross.satellite2.lower()
        date_time_cross = cross.time2
        
    if reprocess is False:
        import os #@Reimport
        diff_avhrr_seconds=None
        avhrr_file=None
        if (PPS_VALIDATION ):
            avhrr_file, tobj = find_radiance_file(cross, options)
            values_avhrr = get_satid_datetime_orbit_from_fname(avhrr_file)
        if (CCI_CLOUD_VALIDATION):
            avhrr_file, tobj = find_cci_cloud_file(cross, options)
        if avhrr_file is not None:
            values_avhrr = get_satid_datetime_orbit_from_fname(avhrr_file)
            date_time_avhrr = values_avhrr["date_time"]
            td = date_time_avhrr-date_time_cross
            diff_avhrr_seconds=abs(td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6


        #need to pUt in the info res, atrain data type before go inte find avhrr??
        # or change read of files
        values["atrain_sat"] = "cloudsat-%s" % config.CLOUDSAT_TYPE
        values["atrain_datatype"] = "cloudsat-%s" % config.CLOUDSAT_TYPE
        cl_match_file, tobj = find_avhrr_file(cross, options['reshape_dir'], options['reshape_file'], values=values)
        if not cl_match_file:
            logger.info("No processed CloudSat match files found." + 
                      " Generating from source data.")
            clObj = None
            date_time=date_time_cross
        else:
            date_time=tobj
            td = tobj- date_time_cross
            matchup_diff_seconds=abs(td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6
            clObj = readCloudsatAvhrrMatchObj(cl_match_file) 
            basename = '_'.join(os.path.basename(cl_match_file).split('_')[1:5])
            if (diff_avhrr_seconds is None or 
                matchup_diff_seconds<=diff_avhrr_seconds or 
                abs(matchup_diff_seconds-diff_avhrr_seconds)<300 or
                abs(matchup_diff_seconds-diff_avhrr_seconds) <300):
                logger.info("CloudSat Matchups read from previously " + 
                          "processed data.")
                date_time=tobj
            else:
                logger.info("CloudSat Matchups will be processed for better match" + 
                          " %s."%values_avhrr["basename"]) 
                logger.info("CloudSat Matchups not read from previously " + 
                          "processed data %s."%basename)  
                clObj = None
                
        values["atrain_sat"] = "caliop"
        values["atrain_datatype"] = "caliop"
        ca_match_file, tobj = find_avhrr_file(cross, options['reshape_dir'], options['reshape_file'], values=values)
        if not  ca_match_file:
            logger.info( 
                      ("No processed CALIPSO match files found. "+
                       "Generating from source data."))
            caObj = None
            date_time=date_time_cross
        else:
            #print ca_match_file
            date_time=tobj
            #print date_time
            #print tobj
            td = tobj- date_time_cross
            matchup_diff_seconds = abs(td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6
            caObj = readCaliopAvhrrMatchObj(ca_match_file)
            basename = '_'.join(os.path.basename(ca_match_file).split('_')[1:5])
            print matchup_diff_seconds, diff_avhrr_seconds
            if (diff_avhrr_seconds is None or 
                matchup_diff_seconds<=diff_avhrr_seconds or 
                abs(matchup_diff_seconds-diff_avhrr_seconds)<300 or
                matchup_diff_seconds<300):
                logger.info( 
                          "CALIPSO Matchups read from previously processed data.")
                logger.info( 'Filename: ' + ca_match_file)
            else:
                logger.info("Calipso Matchups will be processed for better match" + 
                           " %s."%values_avhrr["basename"]) 
                logger.info("Calipso Matchups not read from previously " + 
                           "processed data %s."%basename)  
                caObj = None
        if None in [caObj] and None in [clObj]:
            pass
        else:
            values['date_time'] = date_time
            values['year'] = tobj.year
            values['basename'] = basename
            values['month']="%02d"%(tobj.month)

    #TODO: Fix a better solution for below so it can handle missing cloudsat better.
    if None in [caObj] and None in [clObj]:
        return get_matchups_from_data(cross, options)
    elif None in [clObj] and config.CLOUDSAT_TYPE == 'CWC-RVOD':
        return get_matchups_from_data(cross, options)
    elif caObj is None and config.CLOUDSAT_TYPE == 'GEOPROF':
        return get_matchups_from_data(cross, options)


    return {'calipso': caObj, 'cloudsat': clObj, 
            'basename': basename,'values':values}

def plot_some_figures(clsatObj, caObj, sensor, values, basename, process_mode, 
                      config_options):

    logger.info("Plotting")
    file_type = ['eps', 'png']
    if PLOT_ONLY_PNG==True:
        file_type = ['png']
        
    plotpath = insert_info_in_filename_or_path(config_options['plot_dir'], values,
                                               datetime_obj=values['date_time'])  
    ##TRAJECTORY
    if caObj is not None:
        """
        if clsatObj is not None and config.CLOUDSAT_TYPE == 'GEOPROF':
            cllon = clsatObj.cloudsat.longitude.copy()
            cllat = clsatObj.cloudsat.latitude.copy()
        elif clsatObj is not None and config.CLOUDSAT_TYPE == 'CWC-RVOD':
            cllon = clsatObj.cloudsat.longitude.copy()
            cllat = clsatObj.cloudsat.latitude.copy()
            return
        """
        calon = caObj.calipso.longitude.copy()
        calat = caObj.calipso.latitude.copy()
        trajectorypath = os.path.join(plotpath, "trajectory_plot")
        if not os.path.exists(trajectorypath):
            os.makedirs(trajectorypath)
        trajectoryname = os.path.join(trajectorypath, 
                                      "%skm_%s_trajectory" % (int(config.RESOLUTION),
                                                              values['basename']))
        plotSatelliteTrajectory(calon, 
                                calat,
                                trajectoryname, 
                                config.AREA_CONFIG_FILE, 
                                file_type,
                                **config_options)

    if clsatObj is None or config.CLOUDSAT_TYPE=='GEOPROF':
        #HEIGHT
        drawCalClsatGEOPROFAvhrrPlot(clsatObj, 
                                     caObj, 
                                     caObj.avhrr.avhrr_ctth_cal_ok, 
                                     plotpath,
                                     basename, 
                                     process_mode, 
                                     file_type,
                                     instrument=sensor)
        #TIME DIFF SATZ 
        drawCalClsatAvhrrPlotTimeDiff(clsatObj, 
                                      caObj,
                                      plotpath, basename, 
                                      config.RESOLUTION,
                                      instrument=sensor)
        drawCalClsatAvhrrPlotSATZ(clsatObj, 
                                  caObj,
                                  plotpath, basename, 
                                  config.RESOLUTION, file_type,
                                  instrument=sensor)
    if clsatObj is not None and config.CLOUDSAT_TYPE=='CWC-RVOD':       
        if config.RESOLUTION == 1:
            elevationcwc = np.where(np.less_equal(clsatObj.cloudsatcwc.elevation,0),
                                    -9, clsatObj.cloudsatcwc.elevation)
            data_okcwc = np.ones(clsatObj.cloudsatcwc.elevation.shape,'b')                
            ### 5 KM DATA CWC-RVOD ###                       
        elif config.RESOLUTION == 5: 
            elevationcwc = np.where(np.less_equal(clsatObj.cloudsat5kmcwc.elevation,0),
                                    -9, clsatObj.cloudsat5kmcwc.elevation,-9)
            data_okcwc = np.ones(clsatObj.cloudsat5kmcwc.elevation.shape,'b')
        phase='LW'  
        drawCalClsatCWCAvhrrPlot(clsatObj, 
                                 elevationcwc, 
                                 data_okcwc, 
                                 plotpath, basename, 
                                 phase,
                                 instrument=sensor)
        phase='IW'  
        drawCalClsatCWCAvhrrPlot(clsatObj, 
                                 elevationcwc, 
                                 data_okcwc, 
                                 plotpath, basename, phase,
                                 instrument=sensor)
                
                
def run(cross, process_mode_dnt, config_options, min_optical_depth, reprocess=False):

    """
    The main work horse.
    
    """    
    logger.info("Case: %s" % str(cross))
    logger.info("Process mode: %s" % process_mode_dnt)
    # split process_mode_dnt into two parts. One with process_mode and one dnt_flag
    mode_dnt = process_mode_dnt.split('_')
    if len(mode_dnt) == 1:
        process_mode = process_mode_dnt
        dnt_flag = None
    elif mode_dnt[-1] in ['DAY', 'NIGHT', 'TWILIGHT']:
        process_mode = '_'.join(mode_dnt[0:-1])
        dnt_flag = mode_dnt[-1]
    else:
        process_mode = process_mode_dnt
        dnt_flag = None

    sno_satname = cross.satellite1.lower()
    if sno_satname in ['calipso', 'cloudsat']:
        sno_satname = cross.satellite2.lower()


    sensor = INSTRUMENT.get(sno_satname, 'avhrr')
    logger.info("Sensor = " + sensor)

    #pdb.set_trace()


    # Now fetch all the datasets for the section of the AREA where all
    # three datasets match. Also get maximum and minimum time differences to AVHRR (in seconds)
    matchup_results = get_matchups(cross, config_options, reprocess)
    caObj = matchup_results['calipso']
    clsatObj = matchup_results['cloudsat']
    values = matchup_results['values']
    #import pdb;pdb.set_trace()

    basename = matchup_results['basename']
    base_sat = basename.split('_')[0]
    base_year = basename.split('_')[1][:4]
    base_month = basename.split('_')[1][4:6]
    
    num_cal_data_ok = len(caObj.calipso.elevation)
    if (num_cal_data_ok == 0):
        logger.info("Processing stopped: Zero lenght of matching arrays!")
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit()
    if caObj.calipso.cloudsat_index is None:
        logger.info("Adding some stuff that might not be in older reshaped files")
        clsatObj, caObj = add_additional_clousat_calipso_index_vars(clsatObj, caObj)
    #Calculate hight from sea surface    
    clsatObj, caObj = add_elevation_corrected_imager_ctth( clsatObj, caObj)
    # If mode = OPTICAL_DEPTH -> Change cloud -top and -base profile
    if process_mode == 'OPTICAL_DEPTH':#Remove this if-statement if you always want to do filtering!/KG
        (new_cloud_top, new_cloud_base, new_cloud_fraction, new_fcf) = \
                        CloudsatCloudOpticalDepth(caObj.calipso.layer_top_altitude, 
                                                  caObj.calipso.layer_base_altitude, \
                                                  caObj.calipso.feature_optical_depth_532, 
                                                  caObj.calipso.cloud_fraction, 
                                                  caObj.calipso.feature_classification_flags, 
                                                  min_optical_depth)
        caObj.calipso.layer_top_altitude = new_cloud_top
        caObj.calipso.layer_base_altitude = new_cloud_base
        caObj.calipso.cloud_fraction = new_cloud_fraction
        caObj.calipso.feature_classification_flags = new_fcf

    if (config.COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC and
        (config.ALSO_USE_5KM_FILES or config.RESOLUTION==5) and 
        caObj.calipso.total_optical_depth_5km is None):
        logger.warning("rematched_file is missing total_optical_depth_5km field")
        logger.info("Consider reprocessing with: ")
        logger.info("COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC=True")
        logger.info("ALSO_USE_5KM_FILES=True or RESOLUTION==5")

    check_total_optical_depth_and_warn(caObj)

    if process_mode != 'BASIC' and RESOLUTION==1:
        caObj = CalipsoOpticalDepthHeightFiltering1km(caObj)

    if process_mode == 'OPTICAL_DEPTH_THIN_IS_CLEAR' and RESOLUTION==1:
        logger.info("Setting thin clouds to clear"
                  ", using 5km data in mode OPTICAL_DEPTH_THIN_IS_CLEAR")
        caObj = CalipsoOpticalDepthSetThinToClearFiltering1km(caObj) 

    #=============================================================
    #Get result filename    
    min_depth_to_file_name=""
    if process_mode == 'OPTICAL_DEPTH':
        min_depth_to_file_name="-%.2f"%(min_optical_depth)
    values['mode']= process_mode_dnt+min_depth_to_file_name
    result_path = insert_info_in_filename_or_path(config_options['result_dir'], values, datetime_obj=values['date_time'])
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    result_file = config_options['result_file'].format(
            resolution=str(config.RESOLUTION),
            basename=values['basename'] )
    statfilename = os.path.join(result_path, result_file)                           
    #=============================================================
    # Draw plot
    if process_mode_dnt in config.PLOT_MODES:
        plot_some_figures(clsatObj, caObj, sensor, values, basename, process_mode, 
                          config_options)
    #==============================================================
    #Calculate Statistics
    logger.info("Calculating statistics")
    CalculateStatistics(process_mode, clsatObj, statfilename, caObj, 
                        dnt_flag)
    #=============================================================
