"""
Program cloudsat_calipso_avhrr_match.py is run via process_master.py.

This program is used to process and output statistics for the inter-comparison
of AVHRR PPS results and CloudSat/CALIPSO observations. It may be run
repeatedly and supervised by program process_master.py.

This particular version of Adam's original CloudSat/CALIPSO matchup and
analysis program has been complemented with the following:

 * Program is updated wo be able to use CALIOP-CALIPSO, CPR (CloudSat), AMSR_E 
   or CATS (ISS) as truth. Modules used to handle the truths:
      cloudsat.py
      calipso.py
      amsr.py
      iss.py

 * Program can read satellite data from: PPS, CCI and MAIA. When satllite data 
   comes from PPS-MODIS also modis lvl-2 data can be matched.
   Files to read imager satellite data:
      read_cloudproducts_and_nwp_pps.py  
      read_cloudproducts_maia.py
      read_cloudproducts_cci.py          
      read_modis_products.py

 * Format of matchup files, and reading and writing can be found in matchobject_io.py.

 * The main running program is: process_master.py and compile_stat.py will 
   accumulate statistics.

 * Iamger cloud top height datasets have been recalculated to heights above mean
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


RUNNING INSTRUCTIONS
--------------------

The program is capable of running in a wide range of modes. 
 
Every exection of the program prints out statistics for PPS Cloud Mask, Cloud
Type and Cloud Top Height and CPP directly to file

Software now fully independent of ACPG/AHAMAP

Dependencies: For a successful run of the program the following supporting
              python modules must be available in the default run directory:
              cloudsat.py
              calipso.py
              cloudsat_calipso_avhrr_match.py

Updated 20181001 
Nina


"""

#change log i found in git
import os
import sys

import numpy as np
from config import (IMAGER_INSTRUMENT,
                    NODATA,
                    PLOT_TYPES, 
                    CCI_CLOUD_VALIDATION,
                    MAIA_CLOUD_VALIDATION,
                    PPS_VALIDATION,
                    CALIPSO_version4,
                    CALIPSO_version3,
                    ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS,
                    CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA,
                    ALSO_USE_1KM_FILES,
                    ALSO_USE_SINGLE_SHOT_CLOUD_CLEARED,
                    PPS_FORMAT_2012_OR_EARLIER,
                    MATCH_MODIS_LVL2, 
                    RESOLUTION,
                    USE_EXISTING_RESHAPED_FILES)
import logging
logger = logging.getLogger(__name__)

import config
print config.__file__
from common import MatchupError, ProcessingError

from cloudsat_calipso_avhrr_statistics import (CalculateStatistics)
from plotting.trajectory_plotting import plotSatelliteTrajectory
from plotting.along_track_plotting import (drawCalClsatAvhrrPlotTimeDiff,
                                           drawCalClsatGEOPROFAvhrrPlot, 
                                           drawCalClsatAvhrrPlotSATZ,
                                           drawCalClsatCWCAvhrrPlot)
from cloudsat_calipso_avhrr_prepare import (CalipsoCloudOpticalDepth_new,
                                            check_total_optical_depth_and_warn,
                                            CalipsoOpticalDepthHeightFiltering1km,
                                            detection_height_from_5km_data,
                                            CalipsoOpticalDepthSetThinToClearFiltering1km)

from read_cloudproducts_and_nwp_pps import NWPObj
from cloudsat import (reshapeCloudsat, 
                      match_cloudsat_avhrr ,
                      add_validation_ctth_cloudsat,
                      add_cloudsat_cloud_fraction,
                      mergeCloudsat)
from amsr import (reshapeAmsr, match_amsr_avhrr)
from synop import (reshapeSynop, match_synop_avhrr)
from iss import reshapeIss, match_iss_avhrr
from calipso import (reshapeCalipso, 
                     discardCalipsoFilesOutsideTimeRange,
                     match_calipso_avhrr, 
                     find_break_points, 
                     time_reshape_calipso)
from matchobject_io import (CalipsoObject,
                            writeCaliopAvhrrMatchObj, 
                            readCaliopAvhrrMatchObj,
                            writeCloudsatAvhrrMatchObj, 
                            readCloudsatAvhrrMatchObj,
                            writeIssAvhrrMatchObj, 
                            readIssAvhrrMatchObj,
                            writeAmsrAvhrrMatchObj, 
                            readAmsrAvhrrMatchObj,
                            writeSynopAvhrrMatchObj, 
                            readSynopAvhrrMatchObj
                        )
from calipso import  (add1kmTo5km,
                      addSingleShotTo5km,
                      add5kmVariablesTo1kmresolution,
                      adjust5kmTo1kmresolution,
                      add_validation_ctth_calipso)


#All non-avhrr satellites need to be here. Avhrr is default.
INSTRUMENT = {'npp': 'viirs',
              'noaa18': 'avhrr',
              'meteosat9': 'seviri',
              'noaa20': 'viirs',
              'eos1': 'modis',
              'eos2': 'modis'} 

from datetime import timedelta
from glob import glob

class ppsFiles(object):
    def __init__(self, file_name_dict):
        self.cloudtype = None
        self.ctth = None
        self.cpp = None
        self.cma = None
        self.cmaprob = None
        self.sunsatangles = None
        self.physiography = None
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

def find_truth_files_inner(date_time, time_window, options, values, truth='calipso'):
    """Find the matching Calipso file"""
    tlist = get_time_list(date_time, time_window, 600)
    flist = []
    for tobj in tlist:    
        calipso_dir = insert_info_in_filename_or_path(
            options[truth + '_dir'],
            values, datetime_obj=tobj)
        calipso_file_pattern = insert_info_in_filename_or_path(
            options[truth + '_file'],
            values, 
            datetime_obj=tobj)
        tmplist = glob(os.path.join(calipso_dir, calipso_file_pattern))
        #print "globbing", os.path.join(calipso_dir, calipso_file_pattern)
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
        ctth_type=values.get("ctth_type",""),
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

def find_truth_files(date_time, options, values, truth='calipso'):
    my_sec_THR = config.sec_timeThr
    TRUTH_FILE_LENGTH = config.CALIPSO_FILE_LENGTH 
    if truth in ['cloudsat']:
        TRUTH_FILE_LENGTH = config.CLOUDSAT_FILE_LENGTH
    if truth in ['iss']:
        TRUTH_FILE_LENGTH = config.ISS_FILE_LENGTH
    if truth in ['amsr']:
        TRUTH_FILE_LENGTH = config.AMSR_FILE_LENGTH
    if truth in ['synop']:
        TRUTH_FILE_LENGTH = config.SYNOP_FILE_LENGTH
        my_sec_THR = config.sec_timeThr_synop
    #might need to geth this in before looking for matchups
    tdelta_before = timedelta(seconds = (TRUTH_FILE_LENGTH +
                                         my_sec_THR))
    tdelta = timedelta(
        seconds = (config.SAT_ORBIT_DURATION +  my_sec_THR))

    time_window = (tdelta_before, tdelta)
    t_files = find_truth_files_inner(date_time, time_window, options, 
                                     values, truth=truth)
    if len(t_files) > 1:
        logger.debug("More than one %s file found within time window!", 
                     truth.upper())
    elif len(t_files) == 0:
        logger.info("No %s file found within time window!", 
                    truth.upper())
        return None
    if truth in ['calipso']:    
        t_files = require_h5(t_files)
    t_files = sorted(t_files)
    if truth in ['iss']:
        t_files.sort(key=lambda x: x.rsplit('.')[-2])
    truth_basenames = [ "\n          " + os.path.basename(s) for s in t_files ]
    logger.info("%s files:  %s",truth.upper(), "\n ".join(truth_basenames))
    return t_files


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


def find_avhrr_file(cross, filedir_pattern, filename_pattern, values=None):
    if values is None:
        values = {}
        
    #tlist = get_time_list(cross.time, 
    #                      (timedelta(seconds=0), timedelta(seconds=0)), 
    #                      delta_t_in_seconds=60)
    tlist = [cross.time]
    values["satellite"] = cross.satellite1
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
                logger.info("Found directory %s", found_dir)
        else:
            if not found_dir and avhrr_dir not in checked_dir.keys():
                checked_dir[avhrr_dir] = 1
                logger.info("This directory does not exist, pattern: %s", avhrr_dir)
            continue    
        file_pattern = insert_info_in_filename_or_path(filename_pattern,
                                                       values, datetime_obj=tobj)
        files = glob(os.path.join(found_dir, file_pattern))
        if len(files) > 0:
            no_files_found = False
            logger.info("Found file: %s", os.path.basename(str(files[0])))
            return files[0], tobj
    if not found_dir:       
        return None, None
    if no_files_found:       
        logger.info( "Found no files for patterns of type: %s", file_pattern)
    return None, None


def require_h5(files):
    """
    Convert any '.hdf' files in *files* to '.h5'. Returns a list of '.h5' files.
    
    """
    from config import H4H5_EXECUTABLE
    if not H4H5_EXECUTABLE:
        return files
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
                logger.info("Converting %s to HDF5", f)
                if os.system(command):
                    raise RuntimeError("Couldn't convert %r to HDF5" % f)
                h5_files.append(h5_file)
        else:
            raise ValueError("File format of %r not recognized" % f)    
    return h5_files

def get_pps_file(avhrr_file, options, values, type_of_file, file_dir, 
                 OnlyPrintInDebugMode=False, 
                 FailIfRequestedAndMissing=False):
    if not type_of_file in options and OnlyPrintInDebugMode:
        logger.debug("No %s file in cfg-file.", type_of_file)
        return None  
    elif not type_of_file in options:              
        logger.info("No %s file in cfg-file.", type_of_file)
        return None
    date_time = values["date_time"]
    pps_file_name = insert_info_in_filename_or_path(options[type_of_file], 
                                                     values, datetime_obj=date_time)                        
    path = insert_info_in_filename_or_path(options[file_dir], 
                                           values, datetime_obj=date_time)  
    try:
        file_name = glob(os.path.join(path, pps_file_name))[0]
        if FailIfRequestedAndMissing:
            #This is and important file, tell the user where it was read:
            logger.info("%s-dir: %s", type_of_file.upper().replace('_FILE',''), 
                        os.path.dirname(file_name))
            logger.info("%s-file: %s", type_of_file.upper().replace('_FILE',''), 
                        os.path.basename(file_name))
        return file_name
    except IndexError:
        if FailIfRequestedAndMissing:
            raise MatchupError("No {:s} file found.".format(type_of_file))
        else:
            logger.info("No %s file found corresponding to %s.", type_of_file, avhrr_file)
            return None

def check_cfc_configuration(file_name_dict):
    if  (file_name_dict['cma'] is None and 
         file_name_dict['cloudtype'] is None and 
         file_name_dict['cmaprob'] is None):
        raise MatchupError("No cma, cloudtype or cmaprob file "
                           "found atrain_match.cfg")
    if (config.USE_CT_FOR_CFC_STATISTICS and 
        file_name_dict['cloudtype'] is None):   
        logger.error(
            "\n\tError: USE_CT_FOR_CFC_STATISTICS=True, but ..."
            "\n\t... no cloudtype file in atrain_match.cfg")
        raise MatchupError("Configure problems, see messages above.")
    if (config.USE_CMA_FOR_CFC_STATISTICS and 
        file_name_dict['cma'] is None):   
        logger.error(
            "\n\tError: USE_CMA_FOR_CFC_STATISTICS=True, but..."
            "\n\t... no cma file in atrain_match.cfg")
        raise MatchupError("Configure problems, see messages above.")
    if (config.USE_CMAPROB_FOR_CFC_STATISTICS and 
        file_name_dict['cmaprob'] is None):   
        logger.error(
            "\n\tError: USE_CMAPROB_FOR_CFC_STATISTICS=True, but..."
            "\n\t... no cmaprob file in atrain_match.cfg")
        raise MatchupError("Configure problems, see messages above.")

def find_files_from_avhrr(avhrr_file, options, as_oldstyle=False):
    """
    Find all files needed to process matchup from source data files.
    """
    # Let's get the satellite and the date-time of the pps radiance
    # (avhrr/viirs) file:
    logger.debug("IMAGER: %s", avhrr_file)
    values = get_satid_datetime_orbit_from_fname(avhrr_file,
                                                 as_oldstyle=as_oldstyle)
    date_time = values["date_time"]
    file_name_dict={}
    file_name_dict['cma'] =  get_pps_file(avhrr_file, options, values, 
                             'cma_file', 'cma_dir', 
                             FailIfRequestedAndMissing=True)
    file_name_dict['cloudtype'] = get_pps_file(avhrr_file, options, values, 
                                  'cloudtype_file', 'cloudtype_dir', 
                                  FailIfRequestedAndMissing=True)
    file_name_dict['cmaprob'] = get_pps_file(avhrr_file, options, values, 'cmaprob_file', 'cmaprob_dir', 
                                         FailIfRequestedAndMissing=True) 
    #Check cfc configuration
    check_cfc_configuration(file_name_dict)
    # For CTTH can have several files:    
    ctth_files = {}
    if 'ctth_file' in options.keys():   
        for ctth_type in config.CTTH_TYPES:
            values['ctth_type'] = ctth_type
            ctth_files[ctth_type] = get_pps_file(avhrr_file, options, values, 
                                                 'ctth_file', 'ctth_dir', 
                                                 FailIfRequestedAndMissing=True)
    file_name_dict.update({'ctth': ctth_files})

    file_name_dict['cpp'] = get_pps_file(avhrr_file, options, values, 'cpp_file', 'cpp_dir', 
                                         FailIfRequestedAndMissing=True)                         
    file_name_dict['sunsatangles'] = get_pps_file(avhrr_file, options, values, 'sunsatangles_file', 
                                                  'sunsatangles_dir', 
                                                  FailIfRequestedAndMissing=True)
    file_name_dict['physiography'] = get_pps_file(avhrr_file, options, values, 'physiography_file', 
                                                  'physiography_dir', 
                                                  FailIfRequestedAndMissing=True)
    file_name_dict['r37'] = get_pps_file(avhrr_file, options, values, 'r37_file', 'r37_dir', 
                                         FailIfRequestedAndMissing=True)

    file_name_dict['nwp_tsur'] = get_pps_file(avhrr_file, options, values, 
                                              'nwp_tsur_file', 'nwp_nwp_dir', 
                                              FailIfRequestedAndMissing=True)
    file_name_dict['emis'] = get_pps_file(avhrr_file, options, values, 
                                          'emis_file', 'emis_dir', 
                                          FailIfRequestedAndMissing=True)

    file_name_dict['seaice'] = get_pps_file(avhrr_file, options, values, 
                                            'seaice_file', 'seaice_dir', 
                                            FailIfRequestedAndMissing=True)
    #Textures and Thresholds:
    file_name_dict['text_t11'] = get_pps_file(avhrr_file, options, values, 
                                              'text_t11_file', 'text_dir', 
                                              FailIfRequestedAndMissing=True)
    file_name_dict['thr_t11ts'] = get_pps_file(avhrr_file, options, values, 
                                               'thr_t11ts_file', 'thr_dir', 
                                               FailIfRequestedAndMissing=True)
    ####More NWP, textures and thresholds for v2014:
    for nwp_file in ['nwp_t500','nwp_t700',
                     'nwp_t850','nwp_t950', 'nwp_ciwv', 'nwp_ttro']:   
        file_name_dict[nwp_file] = get_pps_file(avhrr_file, options, values, 
                                                 nwp_file+'_file', 'nwp_nwp_dir', 
                                                 OnlyPrintInDebugMode=True)
    for text_file in ['text_r06',  'text_t37t12', 'text_t37']:
        file_name_dict[text_file] = get_pps_file(avhrr_file, options, values, 
                                                 text_file+'_file', 'text_dir', 
                                                 OnlyPrintInDebugMode=True)

    for thr_file in ['thr_t11ts_inv', 'thr_t11t37_inv', 
                     'thr_t37t12_inv', 'thr_t11t12_inv', 
                     'thr_t11ts', 'thr_t11t37', 'thr_t37t12', 'thr_t11t12',
                     'thr_r09', 'thr_r06', 'thr_t85t11_inv', 'thr_t85t11']:
        file_name_dict[thr_file] = get_pps_file(avhrr_file, options, values, 
                                                thr_file+'_file', 'thr_dir', 
                                                OnlyPrintInDebugMode=True)   
    file_name_dict['nwp_segments'] = get_pps_file(avhrr_file, options,
                                                  values, 
                                                  'segment_file', 'segment_dir') 
    ######################################## 
    ppsfiles = ppsFiles(file_name_dict)
    return  ppsfiles

def get_cloudsat_matchups(cloudsat_files, cloudsat_files_lwp, avhrrGeoObj, avhrrObj,
                          ctype, cma, ctth, nwp_obj, avhrrAngObj, 
                          cpp, nwp_segments,  config_options):
    """
    Read Cloudsat data and match with the given PPS data.
    """
    cloudsat_lwp = None
    cloudsat = None
    if cloudsat_files is not None:   
        logger.debug("Reading cloudsat for type GEOPROF.")
        cloudsat = reshapeCloudsat(cloudsat_files, avhrrGeoObj)
    if cloudsat_files_lwp is not None: 
        logger.debug("Reading cloudsat for type CWC-RVOD.")
        cloudsat_lwp = reshapeCloudsat(cloudsat_files_lwp, avhrrGeoObj)    
    if cloudsat is not None and cloudsat_lwp is not None:
        logger.info("Merging CloudSat GEOPROF and CWC-RVOD data to one object")
        cloudsat = mergeCloudsat(cloudsat, cloudsat_lwp)
    elif cloudsat is None:
        cloudsat = cloudsat_lwp        
    logger.debug("Matching CloudSat with avhrr")
    cl_matchup = match_cloudsat_avhrr(cloudsat,
                                      avhrrGeoObj, avhrrObj, ctype, cma,
                                      ctth, nwp_obj, avhrrAngObj, cpp, nwp_segments)
    return cl_matchup

def get_iss_matchups(iss_files, avhrrGeoObj, avhrrObj,
                     ctype, cma, ctth, nwp_obj, avhrrAngObj, 
                     cpp, nwp_segments,  config_options):
    """
    Read Iss data and match with the given PPS data.
    """
    iss = reshapeIss(iss_files, avhrrGeoObj)
    cl_matchup = match_iss_avhrr(iss,
                                 avhrrGeoObj, avhrrObj, ctype, cma,
                                 ctth, nwp_obj, avhrrAngObj, cpp, nwp_segments)
    return cl_matchup

def get_amsr_matchups(amsr_files, avhrrGeoObj, avhrrObj,
                     ctype, cma, ctth, nwp_obj, avhrrAngObj, 
                     cpp, nwp_segments,  config_options):
    """
    Read Amsr data and match with the given PPS data.
    """
    amsr = reshapeAmsr(amsr_files, avhrrGeoObj)
    am_matchup = match_amsr_avhrr(amsr,
                                  avhrrGeoObj, avhrrObj, ctype, cma,
                                  ctth, nwp_obj, avhrrAngObj, cpp, nwp_segments)
    return am_matchup

def get_synop_matchups(synop_files, avhrrGeoObj, avhrrObj,
                       ctype, cma, ctth, nwp_obj, avhrrAngObj, 
                       cpp, nwp_segments,  config_options):
    """
    Read Synop data and match with the given PPS data.
    """
    synop = reshapeSynop(synop_files, avhrrGeoObj)
    synop_matchup = match_synop_avhrr(synop,
                                  avhrrGeoObj, avhrrObj, ctype, cma,
                                  ctth, nwp_obj, avhrrAngObj, cpp, nwp_segments)
    return synop_matchup

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
    startBreak, endBreak = find_break_points(calipso,  avhrrGeoObj)
    if cafiles1km is not None and CALIPSO_version3 and RESOLUTION == 5:
        #RESOLUTION 5km also have 1km data
        logger.info("Calipso version 3 data used and old 1 km restore method!")
        calipso1km = reshapeCalipso(cafiles1km, res=1)
        calipso5km = calipso
        calipso = add1kmTo5km(calipso1km, calipso5km)
        calipso = time_reshape_calipso(calipso, 
                                       start_break=startBreak, end_break=endBreak) 
        calipso = total_and_top_layer_optical_depth_5km(calipso, resolution=5)

    elif cafiles5km is not None and CALIPSO_version4  and RESOLUTION == 1:
        #RESOLUTION 1km also have 5km data calipso version 4
        logger.info("Calipso version 4, single shot fraction and "
                    "old 5km restored optical depth method used!")
        calipso1km = calipso
        calipso5km = reshapeCalipso(cafiles5km, res=5)
        calipso5km = addSingleShotTo5km(calipso5km) 
        calipso5km = total_and_top_layer_optical_depth_5km(calipso5km, resolution=5)
        calipso1km = add5kmVariablesTo1kmresolution(calipso1km, calipso5km)
        if CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA:
            logger.info("Find detection height using 5km data")
            calipso1km = detection_height_from_5km_data(calipso1km, calipso5km)
        calipso = time_reshape_calipso(calipso1km,  
                                       start_break=startBreak, 
                                       end_break=endBreak) 

    elif cafiles5km is not None and CALIPSO_version3 and RESOLUTION == 1:
        #RESOLUTION 1km also have 5km data calipso version 3
        logger.info("Calipso version 3 data used and old 5 km restored optical depth method!")
        calipso1km = calipso
        calipso5km = reshapeCalipso(cafiles5km, res=5)
        calipso5km = total_and_top_layer_optical_depth_5km(calipso5km, resolution=5)
        calipso1km = add5kmVariablesTo1kmresolution(calipso1km, calipso5km)
        if CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA:
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
        logger.warning("Old metod, only one resolution used, expect bad results!")
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
        
    logger.debug("Matching CALIPSO with avhrr")
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

def read_pps_data(pps_files, avhrr_file):
    from read_cloudproducts_and_nwp_pps import pps_read_all
    return pps_read_all(pps_files, avhrr_file)

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

            if len(calipso1km) == 0:
                raise MatchupError("Did not find any matching 1km calipso files")

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
            if len(calipso5km) == 0:
                raise MatchupError("Did not find any matching 5km calipso files")
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
    
    if clsatObj is not None and caObj is not None:
        caObj.calipso.cal_MODIS_cflag = None
        #map cloudsat to calipso and the other way around!
        from amsr_avhrr.match import match_lonlat
        source = (clsatObj.cloudsat.longitude.astype(np.float64).reshape(-1,1), 
                  clsatObj.cloudsat.latitude.astype(np.float64).reshape(-1,1))
        target = (caObj.calipso.longitude.astype(np.float64).reshape(-1,1), 
                  caObj.calipso.latitude.astype(np.float64).reshape(-1,1))
        mapper, dummy = match_lonlat(source, target, radius_of_influence=1000, n_neighbours=1)
        caObj.calipso.cloudsat_index = mapper.rows.filled(NODATA).ravel()
        target = (clsatObj.cloudsat.longitude.astype(np.float64).reshape(-1,1), 
                  clsatObj.cloudsat.latitude.astype(np.float64).reshape(-1,1))
        source = (caObj.calipso.longitude.astype(np.float64).reshape(-1,1), 
                  caObj.calipso.latitude.astype(np.float64).reshape(-1,1))
        mapper, dummy = match_lonlat(source, target, radius_of_influence=1000, n_neighbours=1)
        clsatObj.cloudsat.calipso_index = mapper.rows.filled(NODATA).ravel()

        # Transfer CloudSat MODIS cloud flag to CALIPSO representation
        index = caObj.calipso.cloudsat_index.copy()
        index[index<0] = 0
        caObj.calipso.cal_MODIS_cflag = np.where(
            caObj.calipso.cloudsat_index>=0, 
            clsatObj.cloudsat.MODIS_cloud_flag[index],
            -9)
        
    if clsatObj is not None and caObj is not None:
        index = clsatObj.cloudsat.calipso_index.copy()
        index[index<0] = 0
        clsatObj.cloudsat.calipso_feature_classification_flags= np.where(
            clsatObj.cloudsat.calipso_index>=0,
            caObj.calipso.feature_classification_flags[index,0],
            -9)
        # first bas layer use height not pressure as cloudsat uses height
        clsatObj.cloudsat.calipso_layer_base_altitude = np.where(
            clsatObj.cloudsat.calipso_index>=0,
            caObj.calipso.layer_base_altitude[index,0],
            -9)
        for layer in xrange(1,10):
            clsatObj.cloudsat.calipso_layer_base_altitude = np.where(
                np.logical_and(np.logical_and(clsatObj.cloudsat.calipso_index>=0,
                                              caObj.calipso.layer_base_altitude[index,layer]>0),
                               caObj.calipso.layer_base_altitude[index,layer]< 
                               caObj.calipso.layer_base_altitude[index,layer-1]),
                
                caObj.calipso.layer_base_altitude[index,layer],
            clsatObj.cloudsat.calipso_layer_base_altitude)
        clsatObj.cloudsat.calipso_layer_base_altitude[clsatObj.cloudsat.calipso_layer_base_altitude<-999] = -9.0
        # first bas layer use height not pressure as cloudsat uses height
        clsatObj.cloudsat.calipso_layer_top_altitude = np.where(
            clsatObj.cloudsat.calipso_index>=0,
            caObj.calipso.layer_top_altitude[index,0],
            -9)
        for layer in xrange(1,10):
            clsatObj.cloudsat.calipso_layer_top_altitude = np.where(
                np.logical_and(np.logical_and(clsatObj.cloudsat.calipso_index>=0,
                                              caObj.calipso.layer_top_altitude[index,layer]>0),
                               caObj.calipso.layer_top_altitude[index,layer]> 
                               caObj.calipso.layer_top_altitude[index,layer-1]),
                
                caObj.calipso.layer_top_altitude[index,layer],
            clsatObj.cloudsat.calipso_layer_top_altitude)
        clsatObj.cloudsat.calipso_layer_top_altitude[clsatObj.cloudsat.calipso_layer_top_altitude<-999] = -9.0
    return clsatObj, caObj

def add_modis_lvl2_clousat_(clsatObj, caObj):
    #add modis lvl2 to cloudsat
    if clsatObj is not None and caObj is not None:
        clsatObj.modis.all_arrays["height"] = caObj.modis.all_arrays["height"][clsatObj.cloudsat.calipso_index]
        clsatObj.modis.all_arrays["temperature"] = caObj.modis.all_arrays["temperature"][clsatObj.cloudsat.calipso_index]
        clsatObj.modis.all_arrays["pressure"] = caObj.modis.all_arrays["pressure"][clsatObj.cloudsat.calipso_index]
        clsatObj.modis.all_arrays["cloud_emissivity"] = caObj.modis.all_arrays["cloud_emissivity"][clsatObj.cloudsat.calipso_index]
        clsatObj.modis.all_arrays["latitude_5km"] = caObj.modis.all_arrays["latitude_5km"][clsatObj.cloudsat.calipso_index]
        clsatObj.modis.all_arrays["longitude_5km"] = caObj.modis.all_arrays["longitude_5km"][clsatObj.cloudsat.calipso_index]
    return clsatObj


def add_elevation_corrected_imager_ctth(clsatObj, caObj, issObj):
    ## Cloudsat ##
    if clsatObj is None or clsatObj.avhrr.ctth_height is None:
        pass
    elif  clsatObj.avhrr.imager_ctth_m_above_seasurface is None:
        # First make sure that PPS cloud top heights are converted to height
        # above sea level just as CloudSat height are defined. Use
        # corresponding DEM data.
        elevation = np.where(np.less_equal(clsatObj.cloudsat.elevation,0),
                             0,clsatObj.cloudsat.elevation)		
        num_csat_data_ok = len(clsatObj.cloudsat.elevation)
        logger.debug("Length of CLOUDSAT array: %d", num_csat_data_ok )
        imager_ctth_m_above_seasurface = np.array(clsatObj.avhrr.ctth_height).copy().ravel()
        if CCI_CLOUD_VALIDATION: 
            #ctth already relative mean sea level
            imager_ctth_m_above_seasurface = caObj.avhrr.ctth_height
        else: #ctth relative topography
            got_height = imager_ctth_m_above_seasurface>=0                    
            imager_ctth_m_above_seasurface[got_height] += elevation[got_height]*1.0
        clsatObj.avhrr.imager_ctth_m_above_seasurface = imager_ctth_m_above_seasurface
        if num_csat_data_ok == 0:
            logger.info("Processing stopped: Zero lenght of matching arrays!")
            raise ProcessingError("Zero lenght of matching arrays!")
    ## Calipso ##        
    # First make sure that PPS cloud top heights are converted to height above sea level
    # just as CALIPSO heights are defined. Use corresponding DEM data.
    if caObj is None or caObj.avhrr.ctth_height is None:
        pass
    elif  caObj.avhrr.imager_ctth_m_above_seasurface is None:
        cal_elevation = np.where(np.less_equal(caObj.calipso.elevation,0),
                                 0,caObj.calipso.elevation)
        num_cal_data_ok = len(caObj.calipso.elevation)
        logger.debug("Length of CALIOP array: %d", num_cal_data_ok)
        imager_ctth_m_above_seasurface = np.array(caObj.avhrr.ctth_height).copy().ravel()
        logger.debug("CCI_CLOUD_VALIDATION %s", str(CCI_CLOUD_VALIDATION))
        if CCI_CLOUD_VALIDATION: 
            #ctth relative mean sea level
            imager_ctth_m_above_seasurface = caObj.avhrr.ctth_height
        else: #ctth relative topography
            got_height = imager_ctth_m_above_seasurface>=0                    
            imager_ctth_m_above_seasurface[got_height] += cal_elevation[got_height]*1.0
        caObj.avhrr.imager_ctth_m_above_seasurface = imager_ctth_m_above_seasurface
    if issObj is None or issObj.avhrr.ctth_height is None:
        pass
    elif  issObj.avhrr.imager_ctth_m_above_seasurface is None:
        iss_elevation = np.where(np.less_equal(issObj.iss.elevation,0),
                                 0,issObj.iss.elevation)
        num_iss_data_ok = len(issObj.iss.elevation)
        logger.info("Length of ISS array: %d", num_iss_data_ok)
        imager_ctth_m_above_seasurface = np.array(issObj.avhrr.ctth_height).copy().ravel()
        if CCI_CLOUD_VALIDATION: 
            #ctth relative mean sea level
            pass
        else: #ctth relative topography
            got_height = imager_ctth_m_above_seasurface>=0                    
            imager_ctth_m_above_seasurface[got_height] += iss_elevation[got_height]*1.0
        issObj.avhrr.imager_ctth_m_above_seasurface = imager_ctth_m_above_seasurface
    return clsatObj, caObj, issObj

def add_validation_ctth(clsatObj, caObj):
    if clsatObj is not None:
        if clsatObj.cloudsat.validation_height is None:
            clsatObj.cloudsat = add_validation_ctth_cloudsat(clsatObj.cloudsat)
        if clsatObj.cloudsat.cloud_fraction is None:   
            clsatObj.cloudsat = add_cloudsat_cloud_fraction(clsatObj.cloudsat) 
    if caObj is not None:
        if caObj.calipso.validation_height is None:
            caObj.calipso = add_validation_ctth_calipso(caObj.calipso)
    return clsatObj, caObj
    
def get_matchups_from_data(cross, config_options):
    """
    Retrieve Cloudsat- and Calipso-AVHRR matchups from Cloudsat, Calipso, and
    PPS files.
    """
    #STEP 1 get imager files
    if (PPS_VALIDATION):
        avhrr_file, tobj = find_radiance_file(cross, config_options)
        pps_files = find_files_from_avhrr(avhrr_file, config_options) 
    if (CCI_CLOUD_VALIDATION):
        avhrr_file, tobj = find_cci_cloud_file(cross, config_options)
    if (MAIA_CLOUD_VALIDATION):
        avhrr_file, tobj = find_maia_cloud_file(cross, config_options)
    if not avhrr_file:
        raise MatchupError("No avhrr file found!\ncross = " + str(cross))
    values = get_satid_datetime_orbit_from_fname(avhrr_file)
    date_time = values["date_time"]

    #Step 2 get truth satellite files
    #CLOUDSAT:  
    cloudsat_files = None
    if (PPS_VALIDATION and config.CLOUDSAT_MATCHING):
        cloudsat_files = find_truth_files(date_time, config_options, values, truth='cloudsat')
    elif CCI_CLOUD_VALIDATION:
        logger.info("\nCCI-cloud validation only for calipso, Continue"
                    "\nIt might be working for CloudSat though ...")
    elif not config.CLOUDSAT_MATCHING:
        logger.info("NO CLOUDSAT File, CLOUDSAT matching not requested "
                    "config.CLOUDSAT_MATCHING=False")  
    #CLOUDSAT LWP from CWC-RVOD:  
    cloudsat_files_lwp = None
    if (PPS_VALIDATION and config.CLOUDSAT_MATCHING and 'cloudsat_lwp_file' in config_options.keys()):
        cloudsat_files_lwp = find_truth_files(date_time, config_options, values, truth='cloudsat_lwp')
    elif CCI_CLOUD_VALIDATION:
        logger.info("\nCCI-cloud validation only for calipso, Continue"
                    "\nIt might be working for CloudSat though ...")
    elif not config.CLOUDSAT_MATCHING:
        logger.info("NO CLOUDSAT File, CLOUDSAT matching not requested "
                    "config.CLOUDSAT_MATCHING=False")  
    #ISS:  
    iss_files = None
    if (PPS_VALIDATION and config.ISS_MATCHING):
        iss_files = find_truth_files(date_time, config_options, values, truth='iss')
    elif CCI_CLOUD_VALIDATION:
        logger.info("\nCCI-cloud validation only for calipso, Continue" 
                    "\nIt might be working for ISS though ...")
    elif not config.ISS_MATCHING:
        logger.info("ISS matching not requested config.ISS_MATCHING=False")     
    #AMSR:  
    amsr_files = None
    if (PPS_VALIDATION and config.AMSR_MATCHING):
        amsr_files = find_truth_files(date_time, config_options, values, truth='amsr')
    elif CCI_CLOUD_VALIDATION:
        logger.info("\nCCI-cloud validation only for calipso, Continue" 
                    "\nIt might be working for AMSR though ...")
    elif not config.AMSR_MATCHING:
        logger.info("AMSR matching not requested config.AMSR_MATCHING=False")   
    #AMSR:  
    synop_files = None
    if (PPS_VALIDATION and config.SYNOP_MATCHING):
        synop_files = find_truth_files(date_time, config_options, values, truth='synop')
    elif CCI_CLOUD_VALIDATION:
        logger.info("\nCCI-cloud validation only for calipso, Continue" 
                    "\nIt might be working for SYNOP though ...")
    elif not config.SYNOP_MATCHING:
        logger.info("SYNOP matching not requested config.SYNOP_MATCHING=False")   
    #CALIPSO:
    calipso_files = None
    if config.CALIPSO_MATCHING:
        calipso_files = find_truth_files(date_time, config_options, values)
        if calipso_files is not None:
            extra_files = get_additional_calipso_files_if_requested(calipso_files)
            calipso5km, calipso1km, calipso5km_aerosol  = extra_files
    else:
        logger.info("CALIPSO matching not requested config.CALIPSO_MATCHING=False")     

    if (calipso_files is None and cloudsat_files is None and 
        cloudsat_files_lwp is None and synop_files is None and
        iss_files is None and amsr_files is None):      
        raise MatchupError(
                "Couldn't find any matching CALIPSO/CLoudSat/ISS data")

    #STEP 3 Read imager data:    
    if (PPS_VALIDATION ):
        retv =read_pps_data(pps_files, avhrr_file)
        (avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, 
         nwp_obj, cpp, nwp_segments, cma )= retv
    if (CCI_CLOUD_VALIDATION):
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surftemp, cpp, cma = read_cloud_cci(avhrr_file)
        nwp_segments = None
        nwp_obj = NWPObj({'surftemp':surftemp}) 
        avhrrGeoObj.satellite = values["satellite"]
    if (MAIA_CLOUD_VALIDATION):
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surftemp, cpp, cma = read_cloud_maia(avhrr_file)
        nwp_segments = None
        nwp_obj = NWPObj({'surftemp':surftemp})   
        avhrrGeoObj.satellite = values["satellite"]

    #STEP 4 get matchups 
    #ClloudSat
    cl_matchup = None
    if (PPS_VALIDATION and config.CLOUDSAT_MATCHING and cloudsat_files is not None):
        logger.info("Read CLOUDSAT data")
        cl_matchup = get_cloudsat_matchups(cloudsat_files, cloudsat_files_lwp, 
                                           avhrrGeoObj, avhrrObj, ctype, cma,
                                           ctth, nwp_obj, avhrrAngObj, cpp, 
                                           nwp_segments, config_options) 
    #ISS:  
    iss_matchup = None
    if (PPS_VALIDATION and config.ISS_MATCHING and iss_files is not None):
        logger.info("Read ISS data")
        iss_matchup = get_iss_matchups(iss_files, 
                                       avhrrGeoObj, avhrrObj, ctype, cma,
                                       ctth, nwp_obj, avhrrAngObj, cpp, 
                                       nwp_segments, config_options)
    #AMSR
    amsr_matchup = None
    if (PPS_VALIDATION and config.AMSR_MATCHING and amsr_files is not None):
        logger.info("Read AMSR data")
        amsr_matchup = get_amsr_matchups(amsr_files, 
                                         avhrrGeoObj, avhrrObj, ctype, cma,
                                         ctth, nwp_obj, avhrrAngObj, cpp, 
                                         nwp_segments, config_options)
    #SYNOP
    synop_matchup = None
    if (PPS_VALIDATION and config.SYNOP_MATCHING and synop_files is not None):
        logger.info("Read SYNOP data")
        synop_matchup = get_synop_matchups(synop_files, 
                                           avhrrGeoObj, avhrrObj, ctype, cma,
                                           ctth, nwp_obj, avhrrAngObj, cpp, 
                                           nwp_segments, config_options)
    #CALIPSO:
    ca_matchup = None
    if config.CALIPSO_MATCHING and calipso_files is not None:
        logger.info("Read CALIPSO data")        
        ca_matchup= get_calipso_matchups(calipso_files, 
                                         values,
                                         avhrrGeoObj, avhrrObj, 
                                         ctype, cma, ctth, 
                                         nwp_obj, avhrrAngObj, 
                                         config_options, cpp,
                                         nwp_segments,
                                         calipso1km, calipso5km, calipso5km_aerosol)
 
    if ca_matchup is None and config.CALIPSO_REQUIRED:
        raise MatchupError("No matches with CALIPSO.")
    elif cl_matchup is None and config.CLOUDSAT_REQUIRED:
        raise MatchupError("No matches with CLOUDAT.")
    elif iss_matchup is None and config.ISS_REQUIRED:
        raise MatchupError("No matches with ISS.")
    elif amsr_matchup is None and config.AMSR_REQUIRED:
        raise MatchupError("No matches with AMSR.")
    elif synop_matchup is None and config.SYNOP_REQUIRED:
        raise MatchupError("No matches with SYNOP.")
    elif (ca_matchup is None and cl_matchup is None and iss_matchup is None 
          and amsr_matchup is None and synop_matchup is None):
        raise MatchupError("No matches with any truth.")


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
        logger.info("Creating dir %s:", rematched_path)
        os.makedirs(os.path.dirname(rematched_path))

    #add modis lvl2    
    if config.MATCH_MODIS_LVL2 and config.IMAGER_INSTRUMENT.lower() in ['modis']:
        from read_modis_products import add_modis_06  
        if ca_matchup is not None:
            ca_matchup = add_modis_06(ca_matchup, avhrr_file, config_options) 
        if cl_matchup is not None:
            cl_matchup = add_modis_06(cl_matchup, avhrr_file, config_options) 
        if amsr_matchup is not None:
            amsr_matchup = add_modis_06(amsr_matchup, avhrr_file, config_options)
        if synop_matchup is not None:
            synop_matchup = add_modis_06(synop_matchup, avhrr_file, config_options)  
    #add additional vars to cloudsat and calipso objects and print them to file:
    cl_matchup, ca_matchup = add_additional_clousat_calipso_index_vars(cl_matchup, ca_matchup)
    cl_matchup, ca_matchup, iss_matchup = add_elevation_corrected_imager_ctth(cl_matchup, ca_matchup, iss_matchup)

    #imager_name
    avhrr_obj_name = 'pps'
    if config.CCI_CLOUD_VALIDATION:
        avhrr_obj_name = 'cci'
    if config.MAIA_CLOUD_VALIDATION:
        avhrr_obj_name = 'maia'

    # Write cloudsat matchup 
    if cl_matchup is not None:
        cl_match_file = rematched_file_base.replace(
            'atrain_datatype', 'cloudsat')
        writeCloudsatAvhrrMatchObj(cl_match_file, cl_matchup, 
                                   avhrr_obj_name = avhrr_obj_name)
    else:
        logger.debug('No CloudSat Match File created')

    # Write iss matchup   
    if iss_matchup is not None:
        is_match_file = rematched_file_base.replace(
            'atrain_datatype', 'iss')
        writeIssAvhrrMatchObj(is_match_file, iss_matchup, 
                                   avhrr_obj_name = avhrr_obj_name)
    else:
        logger.debug('No Iss Match File created')

    # Write amsr matchup   
    if amsr_matchup is not None:
        am_match_file = rematched_file_base.replace(
            'atrain_datatype', 'amsr')
        writeAmsrAvhrrMatchObj(am_match_file, amsr_matchup, 
                               avhrr_obj_name = avhrr_obj_name)
    else:
        logger.debug('No Amsr Match File created')

    # Write synop matchup   
    if synop_matchup is not None:
        am_match_file = rematched_file_base.replace(
            'atrain_datatype', 'synop')
        writeSynopAvhrrMatchObj(am_match_file, synop_matchup, 
                               avhrr_obj_name = avhrr_obj_name)
    else:
        logger.debug('No Synop Match File created')
     
    # Write calipso matchup
    if ca_matchup is not None:
        ca_match_file = rematched_file_base.replace('atrain_datatype', 'caliop')
        writeCaliopAvhrrMatchObj(ca_match_file, ca_matchup, 
                                 avhrr_obj_name = avhrr_obj_name) 
    nwp_obj = None
    return {'cloudsat': cl_matchup, 'calipso': ca_matchup, 
            'iss': iss_matchup, 'amsr': amsr_matchup,
            'synop': synop_matchup, 
            'basename': basename, 'values':values}


def get_matchups(cross, options, reprocess=False):
    """
    Retrieve Cloudsat- and Calipso-AVHRR matchups. If *reprocess* is False, and if
    matchup files exist, get matchups directly from the processed files.
    Otherwise process Cloudsat, Calipso, and PPS files first.
    """
    caObj = None
    clObj = None
    isObj = None
    amObj = None
    syObj = None
    values = {}
    try:
        values["satellite"] = cross.satellite1.lower()
    except AttributeError:
        raise ValueError('Need satellite1 and time (cross: %s)' % cross)
    
    if reprocess is False or USE_EXISTING_RESHAPED_FILES:
        diff_avhrr_seconds=None
        avhrr_file=None
        if not USE_EXISTING_RESHAPED_FILES:
            if (PPS_VALIDATION ):
                avhrr_file, tobj = find_radiance_file(cross, options)
                values_avhrr = get_satid_datetime_orbit_from_fname(avhrr_file)
            if (CCI_CLOUD_VALIDATION):
                avhrr_file, tobj = find_cci_cloud_file(cross, options)
            if avhrr_file is not None:
                values_avhrr = get_satid_datetime_orbit_from_fname(avhrr_file)
                date_time_avhrr = values_avhrr["date_time"]
                td = date_time_avhrr-cross.time
                diff_avhrr_seconds=abs(td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6


        #CLOUDSAT
        if not config.CLOUDSAT_MATCHING:
            logger.info("CloudSat matching turned off config.ISS_MATCHING.")
        else:    
            values["atrain_sat"] = "cloudsat" 
            values["atrain_datatype"] = "cloudsat" 
            cl_match_file,  date_time = find_avhrr_file(cross, options['reshape_dir'], 
                                                  options['reshape_file'], values=values)
            if not cl_match_file:
                logger.info("No processed CloudSat match files found. "  
                            "Generating from source data if required.")
                date_time=cross.time
            else:
                clObj = readCloudsatAvhrrMatchObj(cl_match_file) 
                basename = '_'.join(os.path.basename(cl_match_file).split('_')[1:5])

        #ISS
        if not config.ISS_MATCHING:
            logger.info("ISS matching turned off config.ISS_MATCHING.")
        else:    
            values["atrain_sat"] = "iss"
            values["atrain_datatype"] = "iss"
            is_match_file, date_time = find_avhrr_file(cross, options['reshape_dir'], 
                                                  options['reshape_file'], values=values)
            if not is_match_file:
                logger.info("No processed Iss match files found."  
                            " Generating from source data if required.")
                date_time=cross.time
            else:
                isObj = readIssAvhrrMatchObj(is_match_file) 
                basename = '_'.join(os.path.basename(is_match_file).split('_')[1:5])

        #AMSR
        if not config.AMSR_MATCHING:
            logger.info("AMSR matching turned off config.AMSR_MATCHING.")
        else:    
            values["atrain_sat"] = "amsr"
            values["atrain_datatype"] = "amsr"
            am_match_file, date_time = find_avhrr_file(cross, options['reshape_dir'], 
                                                  options['reshape_file'], values=values)
            if not am_match_file:
                logger.info("No processed Amsr match files found."  
                            " Generating from source data if required.")
                date_time=cross.time
            else:              
                amObj = readAmsrAvhrrMatchObj(am_match_file) 
                basename = '_'.join(os.path.basename(am_match_file).split('_')[1:5])


        #SYNOP
        if not config.SYNOP_MATCHING:
            logger.info("SYNOP matching turned off config.SYNOP_MATCHING.")
        else:    
            values["atrain_sat"] = "synop"
            values["atrain_datatype"] = "synop"
            synop_match_file, date_time = find_avhrr_file(cross, options['reshape_dir'], 
                                                       options['reshape_file'], values=values)
            if not synop_match_file:
                logger.info("No processed Synop match files found."  
                            " Generating from source data if required.")
                date_time=cross.time
            else:              
                syObj = readSynopAvhrrMatchObj(synop_match_file) 
                basename = '_'.join(os.path.basename(synop_match_file).split('_')[1:5])

        #CALIPSO               
        if not config.CALIPSO_MATCHING:
            logger.info("Calipso matching turned off config.AMSR_MATCHING.")
        else:        
            values["atrain_sat"] = "caliop"
            values["atrain_datatype"] = "caliop"
            ca_match_file, date_time = find_avhrr_file(cross, options['reshape_dir'], 
                                                  options['reshape_file'], values=values)
            if not  ca_match_file:
                logger.info( 
                    ("No processed CALIPSO match files found. "
                     "Generating from source data."))
                date_time=cross.time
            else:            
                caObj = readCaliopAvhrrMatchObj(ca_match_file)
                basename = '_'.join(os.path.basename(ca_match_file).split('_')[1:5])


    if (caObj is None and clObj is None and isObj is None and 
        amObj is None and syObj is None):
        pass
    else:
        values['date_time'] = date_time 
        values['year'] = date_time.year      
        values['basename'] = basename
        values['month']="%02d"%(date_time.month)

    if config.USE_EXISTING_RESHAPED_FILES:
        if caObj is None and config.CALIPSO_REQUIRED:
            raise MatchupError(
                "Couldn't find calipso already processed matchup file, "
                "USE_EXISTING_RESHAPED_FILES = True!") 
        if clObj is None and config.CLOUDSAT_REQUIRED:
            raise MatchupError(
                "Couldn't find cloudsat already processed matchup file," 
                "USE_EXISTING_RESHAPED_FILES = True!") 
        if isObj is None and config.ISS_REQUIRED:
            raise MatchupError(
                "Couldn't find iss already processed matchup file," 
                "USE_EXISTING_RESHAPED_FILES = True!") 
        if amObj is None and config.AMSR_REQUIRED:
            raise MatchupError(
                "Couldn't find amsr already processed matchup file," 
                "USE_EXISTING_RESHAPED_FILES = True!") 
        if syObj is None and config.SYNOP_REQUIRED:
            raise MatchupError(
                "Couldn't find synop already processed matchup file," 
                "USE_EXISTING_RESHAPED_FILES = True!") 

    if (caObj is None and clObj is None and isObj is None and 
        amObj is None  and syObj is None):
        out_dict = get_matchups_from_data(cross, options) 
    elif caObj is None and config.CALIPSO_REQUIRED:
        out_dict = get_matchups_from_data(cross, options)
    elif clObj is None and config.CLOUDSAT_REQUIRED:
        out_dict =  get_matchups_from_data(cross, options)
    elif isObj is None and config.ISS_REQUIRED:
        out_dict =  get_matchups_from_data(cross, options)
    elif amObj is None and config.AMSR_REQUIRED:
        out_dict =  get_matchups_from_data(cross, options)
    elif syObj is None and config.SYNOP_REQUIRED:
        out_dict =  get_matchups_from_data(cross, options)
    else:
        out_dict = {'calipso': caObj, 'cloudsat': clObj,  
                    'iss': isObj, 'amsr': amObj,
                    'synop': syObj,
                    'basename': basename,'values':values}

    if out_dict['cloudsat'] is None and config.CLOUDSAT_REQUIRED:
        raise MatchupError("Couldn't find cloudsat matchup and "
                           "CLOUDSAT_REQUIRED is True!")
    if out_dict['calipso'] is None and config.CALIPSO_REQUIRED:
        raise MatchupError("Couldn't find calipso matchup!"
                           "CALIPSO_REQUIRED is True!")
    if out_dict['iss'] is None and config.ISS_REQUIRED:
        raise MatchupError("Couldn't find iss matchup!"
                           "ISS_REQUIRED is True!")   
    if out_dict['amsr'] is None and config.AMSR_REQUIRED:
        raise MatchupError("Couldn't find amsr matchup!"
                           "AMSR_REQUIRED is True!") 
    if out_dict['synop'] is None and config.SYNOP_REQUIRED:
        raise MatchupError("Couldn't find synop matchup!"
                           "SYNOP_REQUIRED is True!") 
    return out_dict

def plot_some_figures(clsatObj, caObj, values, basename, process_mode, 
                      config_options, amObj = None, synopObj = None):

    logger.info("Plotting")
    file_type = ['eps', 'png', 'pdf']
    file_type = PLOT_TYPES 
        
    plotpath = insert_info_in_filename_or_path(config_options['plot_dir'], values,
                                               datetime_obj=values['date_time'])  

    ##TRAJECTORY
    if caObj is not None:
        imlon = caObj.avhrr.longitude.copy()
        imlat = caObj.avhrr.latitude.copy()
        trajectorypath = os.path.join(plotpath, "trajectory_plot")
        if not os.path.exists(trajectorypath):
            os.makedirs(trajectorypath)
        trajectoryname = os.path.join(trajectorypath, 
                                      "%skm_%s_trajectory" % (int(config.RESOLUTION),
                                                              values['basename']))
        plotSatelliteTrajectory(imlon, 
                                imlat,
                                trajectoryname, 
                                config.AREA_CONFIG_FILE, 
                                file_type,
                                **config_options)

    if (clsatObj is not None and 
        "CPR_Cloud_mask" in clsatObj.cloudsat.all_arrays.keys() and 
        caObj is not None):
        #HEIGHT
        drawCalClsatGEOPROFAvhrrPlot(clsatObj, 
                                     caObj, 
                                     caObj.avhrr.imager_ctth_m_above_seasurface, 
                                     plotpath,
                                     basename, 
                                     process_mode, 
                                     file_type,
                                     instrument=IMAGER_INSTRUMENT)
        #TIME DIFF SATZ 
        drawCalClsatAvhrrPlotTimeDiff(clsatObj, 
                                      caObj,
                                      plotpath, basename, 
                                      config.RESOLUTION,
                                      instrument=IMAGER_INSTRUMENT)
        drawCalClsatAvhrrPlotSATZ(clsatObj, 
                                  caObj,
                                  plotpath, basename, 
                                  config.RESOLUTION, file_type,
                                  instrument=IMAGER_INSTRUMENT)

    if (clsatObj is not None and 
        'RVOD_liq_water_path' in clsatObj.cloudsat.all_arrays.keys()):

        elevation = np.where(np.less_equal(clsatObj.cloudsat.elevation,0),
                             -9, clsatObj.cloudsat.elevation)
        data_ok = np.ones(clsatObj.cloudsat.elevation.shape,'b')                

        phase='LW'  
        drawCalClsatCWCAvhrrPlot(clsatObj, 
                                 elevation, 
                                 data_ok, 
                                 plotpath, basename, 
                                 phase,
                                 instrument=IMAGER_INSTRUMENT)
        phase='IW'  
        drawCalClsatCWCAvhrrPlot(clsatObj, 
                                 elevation, 
                                 data_ok, 
                                 plotpath, basename, phase,
                                 instrument=IMAGER_INSTRUMENT)

def split_process_mode_and_dnt_part(process_mode_dnt):        
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
    return process_mode, dnt_flag
 

def process_one_mode(process_mode_dnt, caObj, clsatObj, issObj, amObj,syObj,
                     min_optical_depth, values, config_options, basename):
    
    #Get result filename
    #=============================================================
    process_mode, dnt_flag = split_process_mode_and_dnt_part(process_mode_dnt)
    min_depth_to_file_name = ""
    if process_mode == 'OPTICAL_DEPTH':
        min_depth_to_file_name="-%.2f"%(min_optical_depth)
    values['mode']= process_mode_dnt + min_depth_to_file_name
    result_path = insert_info_in_filename_or_path(config_options['result_dir'], 
                                                  values, 
                                                  datetime_obj=values['date_time'])
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    result_file = config_options['result_file'].format(
        resolution=str(config.RESOLUTION),
        basename=values['basename'],
        truth_sat = "xxx")
    statfilename = os.path.join(result_path, result_file)                           
    #=============================================================
    # Draw plot
    logger.debug("Plotting")
    if process_mode_dnt in config.PLOT_MODES:
        plot_some_figures(clsatObj, caObj, values, basename, process_mode, 
                          config_options, amObj=amObj)
    #==============================================================
    #Calculate Statistics
    logger.debug("Calculating statistics")
    CalculateStatistics(process_mode, statfilename, caObj, clsatObj, 
                        issObj, amObj, syObj, dnt_flag)
    #=============================================================


                               
def run(cross, run_modes, config_options, reprocess=False):
    """
    The main work horse.    
    """    
    logger.info("Case: %s", str(cross))
    sensor = INSTRUMENT.get(cross.satellite1.lower(), 'avhrr')
    if sensor.lower() != IMAGER_INSTRUMENT.lower() :
        logger.error("Uncertain of sensor: %s or %s?", 
                     sensor.upper(), IMAGER_INSTRUMENT.upper())
    if (not config.USE_CMA_FOR_CFC_STATISTICS and 
        not config.USE_CT_FOR_CFC_STATISTICS and
        not config.USE_CMAPROB_FOR_CFC_STATISTICS):
        logger.error(
            "\n###########################"
            "\n\tSet one of USE_*_FOR_CFC_STATISTICS=True in config.py!"
            "\n###########################")
        raise MatchupError("Configure problems, see messages above.")

    #Get the data that we need:
    matchup_results = get_matchups(cross, config_options, reprocess)
    caObj = matchup_results['calipso']
    issObj = matchup_results['iss']
    amObj = matchup_results['amsr']
    syObj = matchup_results['synop']
    clsatObj = matchup_results['cloudsat']
    values = matchup_results['values']
    basename = matchup_results['basename']
    if caObj is not None and caObj.calipso.cloudsat_index is None:
        logger.info("Adding stuff missing in old reshaped files")
        clsatObj, caObj = add_additional_clousat_calipso_index_vars(clsatObj, caObj)
    logger.info("Adding validation height missing in old reshaped files")
    clsatObj, caObj = add_validation_ctth(clsatObj, caObj)
    #Calculate hight from sea surface 
    clsatObj, caObj, issObj = add_elevation_corrected_imager_ctth(clsatObj, caObj, issObj)
    calipso_original = CalipsoObject()
    #Save data orignal data that we might edit for some modes
    if caObj is not None:
        calipso_original.layer_top_altitude = caObj.calipso.layer_top_altitude.copy()
        calipso_original.layer_base_altitude = caObj.calipso.layer_base_altitude.copy()
        calipso_original.cloud_fraction = caObj.calipso.cloud_fraction.copy()
        calipso_original.feature_classification_flags = caObj.calipso.feature_classification_flags.copy()
        calipso_original.validation_height = caObj.calipso.validation_height.copy()
        calipso_original.layer_top_pressure = caObj.calipso.layer_top_pressure.copy()
        calipso_original.layer_base_pressure = caObj.calipso.layer_base_pressure.copy()

    #For each mode, do the statistics:
    if (caObj is not None and 
        config.COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC and
        (config.ALSO_USE_5KM_FILES or config.RESOLUTION==5) and 
        caObj.calipso.total_optical_depth_5km is None):
        logger.warning("\n\t Rematched_file is missing total_optical_depth_5km field"
                       "\n\t Consider reprocessing with: "
                       "\n\t COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC=True"
                       "\n\t ALSO_USE_5KM_FILES=True or RESOLUTION==5")
        
    for process_mode_dnt in run_modes:
        logger.info("Process mode: %s", process_mode_dnt)
        optical_depths = [None]         #Update this if you always want to do filtering!/Nina
        if process_mode_dnt in ["OPTICAL_DEPTH","OPTICAL_DEPTH_DAY",
                                "OPTICAL_DEPTH_NIGHT","OPTICAL_DEPTH_TWILIGHT"]:
            optical_depths = config.MIN_OPTICAL_DEPTH           
            
        # split process_mode_dnt into two parts. One with process_mode and one dnt_flag
        process_mode, dnt_flag = split_process_mode_and_dnt_part(process_mode_dnt)
        for min_optical_depth in optical_depths:
            #For some modes these are updated, so reset calipso data to original
            if caObj is not None:
                #########################################################################
                caObj.calipso.layer_top_altitude = calipso_original.layer_top_altitude.copy()
                caObj.calipso.layer_base_altitude = calipso_original.layer_base_altitude.copy()
                caObj.calipso.cloud_fraction = calipso_original.cloud_fraction.copy()
                caObj.calipso.feature_classification_flags = calipso_original.feature_classification_flags.copy()
                caObj.calipso.validation_height = calipso_original.validation_height.copy()
                caObj.calipso.layer_top_pressure = calipso_original.layer_top_pressure.copy()
                caObj.calipso.layer_base_pressure = calipso_original.layer_base_pressure.copy()
                #########################################################################
            # If mode = OPTICAL_DEPTH -> Change cloud -top and -base profile
            if caObj is not None and process_mode == 'OPTICAL_DEPTH': 
                use_old_method = config.KG_OLD_METHOD_CLOUD_CENTER_AS_HEIGHT
                retv = CalipsoCloudOpticalDepth_new(
                    caObj.calipso,
                    min_optical_depth,
                    use_old_method=use_old_method)
                caObj.calipso.layer_top_altitude = retv[0]
                caObj.calipso.layer_base_altitude = retv[1]
                caObj.calipso.cloud_fraction = retv[2]
                caObj.calipso.feature_classification_flags = retv[3]
                caObj.calipso.validation_height = retv[4]
                caObj.calipso.layer_top_pressure = retv[5]
                caObj.calipso.layer_base_pressure = retv[6]
           # If mode = STANDARD -> Change cloud -top and -base profile
            if caObj is not None and process_mode == 'STANDARD' and RESOLUTION==5:
                retv = CalipsoCloudOpticalDepth_new(caObj.calipso, 0.0)
                caObj.calipso.layer_top_altitude = retv[0]
                caObj.calipso.layer_base_altitude = retv[1]
                caObj.calipso.cloud_fraction = retv[2]
                caObj.calipso.feature_classification_flags = retv[3]
                caObj.calipso.validation_height = retv[4]
                caObj.calipso.layer_top_pressure = retv[5]
                caObj.calipso.layer_base_pressure = retv[6]
            if caObj is not None:    
                check_total_optical_depth_and_warn(caObj)
                if  caObj is not None and 'STANDARD' in process_mode and RESOLUTION==1:
                    caObj.calipso.validation_height = CalipsoOpticalDepthHeightFiltering1km(caObj)
            if  caObj is not None and process_mode == 'OPTICAL_DEPTH_THIN_IS_CLEAR':
                logger.info("Setting thin clouds to clear"
                            ", using 5km data in mode OPTICAL_DEPTH_THIN_IS_CLEAR")
                retv = CalipsoOpticalDepthSetThinToClearFiltering1km(caObj) 
                caObj.calipso.cloud_fraction = retv[0]
                caObj.calipso.validation_height = retv[1]
            #Time to process results files for one mode:    
            process_one_mode(process_mode_dnt, 
                             caObj, clsatObj, issObj, amObj, syObj,   
                             min_optical_depth, values, 
                             config_options, basename)
    #We are done, free some memory:        
    caObj = None
    clsatObj = None
    issObj = None
    amObj = None
