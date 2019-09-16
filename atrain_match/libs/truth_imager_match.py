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
Program truth_imager_match.py is run via process_master.py or process_master_only_match.py.

This program is used to process matchups for the inter-comparison
of IMAGER PPS results and CloudSat/CALIPSO observations. It may be run
repeatedly and supervised or process_master_only_match.py.

This particular version of Adam's original CloudSat/CALIPSO matchup and
analysis program has been complemented with the following:

 * Program is updated wo be able to use CALIOP-CALIPSO, CPR (CloudSat), AMSR_E 
   or CATS (ISS) as truth. Modules used to handle the truths:
      cloudsat.py
      calipso.py
      amsr.py
      iss.py
      synop.py
      mora.py

 * Program can read satellite data from: PPS, CCI and MAIA. When satllite data 
   comes from PPS-MODIS also modis lvl-2 data can be matched.
   Files to read imager satellite data:
      read_cloudproducts_and_nwp_pps.py  
      read_cloudproducts_maia.py
      read_cloudproducts_cci.py          
      read_modis_products.py
      read_cloudproducts_patmosx.py 

 * Format of matchup files, and reading and writing can be found in matchobject_io.py.

 * The main running program is: process_master_only_match.py or process_master.py.

 * From 20190913 statistcs will be handled by truth_imager_make_statistics.py.

 * Imager cloud top height datasets are recalculated to heights above mean
   sea level using CloudSat and CALIPSO elevation data

 * The MODIS cloud flag has been added to the extracted CALIPSO dataset. This
   enables direct comparisons to the MODIS cloud mask! Consequently,
   corresponding MODIS Cloud Mask statistics are calculated and printed.

Updated 20190913 
Nina


"""

#change log i found in git
import os
import sys
import numpy as np
import logging

from datetime import timedelta, datetime
try: 
    from datetime import timezone as my_tz
except ImportError:
    import pytz as my_tz
from glob import glob
import time

logger = logging.getLogger(__name__)

import config
from config import INSTRUMENT

from utils.common import MatchupError, ProcessingError
from truths.cloudsat import (reshapeCloudsat, 
                             match_cloudsat_imager,
                             merge_cloudsat)
from truths.amsr import (reshape_amsr, match_amsr_imager)
from truths.mora import (reshape_mora, match_mora_imager)
from truths.synop import (reshape_synop, match_synop_imager)
from truths.iss import reshape_iss, match_iss_imager
from truths.calipso import ( #check_total_optical_depth_and_warn,
                            total_and_top_layer_optical_depth_5km,
                            reshape_calipso, 
                            discard_calipso_files_outside_time_range,
                            match_calipso_imager, 
                            find_break_points,
                            add_1km_to_5km,
                            add_singleshot_to5km,
                            add_5km_variables_to_1km,
                            adjust_5km_to_1km_resolution)
from matchobject_io import (write_truth_imager_match_obj, 
                            read_truth_imager_match_obj)

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
    tobj1 = cross_time #end_time
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

def find_closest_nwp_file(cloudproducts, AM_PATHS, values, SETTINGS):
    date_time = datetime.fromtimestamp((cloudproducts.sec1970_end*0.5 + 
                                        cloudproducts.sec1970_start*0.5), my_tz.utc)
    delta_3h = timedelta(hours=SETTINGS['MAX_NWP_TDIFF_HOURS'])
    tlist = get_time_list(date_time,  [delta_3h, delta_3h], 60*60) # time_window +/- 3h
    for tobj in tlist: 
        for prognosis_length in range(15):
            values['plus_hours'] = "{:03d}".format(prognosis_length)
            grib_datetime = tobj - timedelta(hours=prognosis_length)
            grib_dir = insert_info_in_filename_or_path(
                AM_PATHS['grib_dir'],
                values, datetime_obj=grib_datetime)
            grib_file_pattern = insert_info_in_filename_or_path(
                AM_PATHS['grib_file'],
                values, 
                datetime_obj=grib_datetime)
            #print("globbing", os.path.join(grib_dir, grib_file_pattern))
            tmplist = glob(os.path.join(grib_dir, grib_file_pattern))
            if len(tmplist)>0:
                return tmplist[0]
    return None

def find_truth_files_inner(date_time, time_window, AM_PATHS, values, truth='calipso'):
    """Find the matching Calipso file"""
    tlist = get_time_list(date_time, time_window, 600)
    flist = []
    for tobj in tlist:    
        calipso_dir = insert_info_in_filename_or_path(
            AM_PATHS[truth + '_dir'],
            values, datetime_obj=tobj)
        calipso_file_pattern = insert_info_in_filename_or_path(
            AM_PATHS[truth + '_file'],
            values, 
            datetime_obj=tobj)
        tmplist = glob(os.path.join(calipso_dir, calipso_file_pattern))
        logger.debug("globbing" + os.path.join(calipso_dir, calipso_file_pattern))
        flist.extend([ s for s in tmplist if s not in flist ])      
    return flist

def get_satid_datetime_orbit_from_fname(filename, SETTINGS, cross=None, as_oldstyle=False):
    from imager_cloud_products.read_cloudproducts_and_nwp_pps import get_satid_datetime_orbit_from_fname_pps
    from imager_cloud_products.read_cloudproducts_cci import get_satid_datetime_orbit_from_fname_cci
    from imager_cloud_products.read_cloudproducts_oca import get_satid_datetime_orbit_from_fname_oca
    from imager_cloud_products.read_cloudproducts_maia import get_satid_datetime_orbit_from_fname_maia
    from imager_cloud_products.read_cloudproducts_patmosx import get_satid_datetime_orbit_from_fname_patmosx
    #Get satellite name, time, and orbit number from imager_file
    if SETTINGS['PPS_VALIDATION']:
        values = get_satid_datetime_orbit_from_fname_pps(filename, as_oldstyle=as_oldstyle)  
    if SETTINGS['CCI_CLOUD_VALIDATION']:
        values = get_satid_datetime_orbit_from_fname_cci(filename)
    if SETTINGS['MAIA_VALIDATION']:
        values = get_satid_datetime_orbit_from_fname_maia(filename)
    if SETTINGS['PATMOSX_VALIDATION']:
        values = get_satid_datetime_orbit_from_fname_patmosx(filename, SETTINGS, cross)
    if SETTINGS['OCA_VALIDATION']:
        values = get_satid_datetime_orbit_from_fname_oca(filename, SETTINGS, cross)
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
        lines_lines=values.get("lines_lines", "*"),
        val_dir=config._validation_results_dir,
        extrai=values.get('extrai',""),
        plus_hours=values.get("plus_hours",""),
        year=values.get('year',"unknown"),
        month=values.get('month',"unknown"),
        mode=values.get('mode',"unknown"),
        #min_opt_depth=values.get('min_opt_depth',""),
        atrain_datatype=values.get("atrain_datatype","atrain_datatype"))
    if datetime_obj is None:
        return file_or_name_path
    name = datetime_obj.strftime(file_or_name_path)
    return name

def find_truth_files(date_time, AM_PATHS, SETTINGS, values, truth='calipso'):
    my_sec_THR = SETTINGS['sec_timeThr']
    TRUTH_FILE_LENGTH = config.CALIPSO_FILE_LENGTH
    if truth in ['cloudsat']:
        TRUTH_FILE_LENGTH = config.CLOUDSAT_FILE_LENGTH
    if truth in ['iss']:
        TRUTH_FILE_LENGTH = config.ISS_FILE_LENGTH
    if truth in ['amsr']:
        TRUTH_FILE_LENGTH = config.AMSR_FILE_LENGTH
    if truth in ['mora']:
        my_sec_THR = SETTINGS['sec_timeThr_synop']
        TRUTH_FILE_LENGTH = config.MORA_FILE_LENGTH
    if truth in ['synop']:
        TRUTH_FILE_LENGTH = config.SYNOP_FILE_LENGTH
        my_sec_THR = SETTINGS['sec_timeThr_synop']
    #might need to geth this in before looking for matchups
    tdelta_before = timedelta(seconds = (TRUTH_FILE_LENGTH +
                                         my_sec_THR))
    tdelta = timedelta(
        seconds = (SETTINGS['SAT_ORBIT_DURATION'] +  my_sec_THR))

    time_window = (tdelta_before, tdelta)
    t_files = find_truth_files_inner(date_time, time_window, AM_PATHS, 
                                     values, truth=truth)
    if len(t_files) > 1:
        logger.debug("More than one %s file found within time window!", 
                     truth.upper())
    elif len(t_files) == 0:
        logger.info("No %s file found within time window!", 
                    truth.upper())
        return None
    if truth in ['calipso']:    
        t_files = require_h5(t_files, SETTINGS)
    t_files = sorted(t_files)
    if truth in ['iss']:
        t_files.sort(key=lambda x: x.rsplit('.')[-2])
    truth_basenames = [ "\n          " + os.path.basename(s) for s in t_files ]
    logger.info("%s files:  %s",truth.upper(), "\n ".join(truth_basenames))
    return t_files
    

def find_radiance_file(cross, AM_PATHS):
    found_file, tobj= find_imager_file(cross, 
                                      AM_PATHS['radiance_dir'], 
                                      AM_PATHS['radiance_file'])
    if not found_file:
        raise MatchupError("No dir or file found with radiance data!\n" + 
                           "Searching for %s %s" % (AM_PATHS['radiance_dir'],AM_PATHS['radiance_file']))
    return found_file, tobj
def find_cci_cloud_file(cross, AM_PATHS):
    found_file, tobj= find_imager_file(cross, 
                                      AM_PATHS['cci_dir'], 
                                      AM_PATHS['cci_file'])
    if not found_file:
        raise MatchupError("No dir or file found with cci cloud data!\n" + 
                           "Searching under %s" % AM_PATHS['cci_dir'])
    return found_file, tobj
def find_oca_cloud_file(cross, AM_PATHS):
    found_file, tobj= find_imager_file(cross, 
                                      AM_PATHS['oca_dir'], 
                                      AM_PATHS['oca_file'])
    if not found_file:
        raise MatchupError("No dir or file found with oca cloud data!\n" + 
                           "Searching under %s" % AM_PATHS['oca_dir'])
    return found_file, tobj
def find_maia_cloud_file(cross, AM_PATHS):
    found_file, tobj= find_imager_file(cross, 
                                       AM_PATHS['maia_dir'], 
                                       AM_PATHS['maia_file'])
    if not found_file:
        raise MatchupError("No dir or file found with maia cloud data!\n" + 
                           "Searching under %s" % AM_PATHS['maia_dir'])
    return found_file, tobj

def find_patmosx_cloud_file(cross, AM_PATHS):
    found_file, tobj= find_imager_file(cross, 
                                       AM_PATHS['patmosx_dir'], 
                                       AM_PATHS['patmosx_file'])
    if not found_file:
        raise MatchupError("No dir or file found with patmosx cloud data!\n" + 
                           "Searching under %s" % AM_PATHS['patmosx_dir'])
    return found_file, tobj


def find_imager_file(cross, filedir_pattern, filename_pattern, values=None):
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
        imager_dir = insert_info_in_filename_or_path(filedir_pattern,
                                                    values, datetime_obj=tobj)
        #print imager_dir
        if os.path.exists(imager_dir):
            found_dir = imager_dir
            if imager_dir not in checked_dir.keys():
                checked_dir[imager_dir] = 1
                logger.info("Found directory %s", found_dir)
        else:
            if not found_dir and imager_dir not in checked_dir.keys():
                checked_dir[imager_dir] = 1
                logger.info("This directory does not exist, pattern: %s", imager_dir)
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


def require_h5(files, SETTINGS):
    """
    Convert any '.hdf' files in *files* to '.h5'. Returns a list of '.h5' files.
    
    """
    if not SETTINGS["H4H5_EXECUTABLE"]:
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

def get_pps_file(imager_file, AM_PATHS, values, type_of_file, file_dir, 
                 OnlyPrintInDebugMode=False, 
                 FailIfRequestedAndMissing=False):
    if not type_of_file in AM_PATHS and OnlyPrintInDebugMode:
        logger.debug("No %s file in cfg-file.", type_of_file)
        return None  
    elif not type_of_file in AM_PATHS:              
        logger.info("No %s file in cfg-file.", type_of_file)
        return None
    date_time = values["date_time"]
    pps_file_name = insert_info_in_filename_or_path(AM_PATHS[type_of_file], 
                                                     values, datetime_obj=date_time)                        
    path = insert_info_in_filename_or_path(AM_PATHS[file_dir], 
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
            logger.info("No %s file found corresponding to %s.", type_of_file, imager_file)
            return None

def check_cfc_configuration(file_name_dict, SETTINGS):
    if  (file_name_dict['cma'] is None and 
         file_name_dict['cloudtype'] is None and 
         file_name_dict['cmaprob'] is None):
        raise MatchupError("No cma, cloudtype or cmaprob file "
                           "found atrain_match.cfg")
    if (SETTINGS['USE_CT_FOR_CFC_STATISTICS'] and 
        file_name_dict['cloudtype'] is None):   
        logger.error(
            "\n\tError: USE_CT_FOR_CFC_STATISTICS=True, but ..."
            "\n\t... no cloudtype file in atrain_match.cfg")
        raise MatchupError("Configure problems, see messages above.")
    if (SETTINGS['USE_CMA_FOR_CFC_STATISTICS'] and 
        file_name_dict['cma'] is None):   
        logger.error(
            "\n\tError: USE_CMA_FOR_CFC_STATISTICS=True, but..."
            "\n\t... no cma file in atrain_match.cfg")
        raise MatchupError("Configure problems, see messages above.")
    if (SETTINGS['USE_CMAPROB_FOR_CFC_STATISTICS'] and 
        file_name_dict['cmaprob'] is None):   
        logger.error(
            "\n\tError: USE_CMAPROB_FOR_CFC_STATISTICS=True, but..."
            "\n\t... no cmaprob file in atrain_match.cfg")
        raise MatchupError("Configure problems, see messages above.")

def find_files_from_imager(imager_file, AM_PATHS, SETTINGS, as_oldstyle=False):
    """
    Find all files needed to process matchup from source data files.
    """
    # Let's get the satellite and the date-time of the pps radiance
    # (imager/viirs) file:
    logger.debug("IMAGER: %s", imager_file)
    values = get_satid_datetime_orbit_from_fname(imager_file,
                                                 SETTINGS,
                                                 cross=None, 
                                                 as_oldstyle=as_oldstyle)
    date_time = values["date_time"]
    file_name_dict={}
    file_name_dict['cma'] =  get_pps_file(imager_file, AM_PATHS, values, 
                             'cma_file', 'cma_dir', 
                             FailIfRequestedAndMissing=True)
    file_name_dict['cloudtype'] = get_pps_file(imager_file, AM_PATHS, values, 
                                               'cloudtype_file', 'cloudtype_dir', 
                                               FailIfRequestedAndMissing=True)
    file_name_dict['cmaprob'] = get_pps_file(imager_file, AM_PATHS, values, 'cmaprob_file', 'cmaprob_dir', 
                                             FailIfRequestedAndMissing=True) 
    #Check cfc configuration
    check_cfc_configuration(file_name_dict, SETTINGS)
    # For CTTH can have several files:    
    ctth_files = {}
    file_name_dict['ctth'] = ctth_files #If no ctth matching requested
    if 'ctth_file' in AM_PATHS.keys():  
        for ctth_type in SETTINGS['CTTH_TYPES']:
            values['ctth_type'] = ctth_type
            ctth_files[ctth_type] = get_pps_file(imager_file, AM_PATHS, values, 
                                                 'ctth_file', 'ctth_dir', 
                                                 FailIfRequestedAndMissing=True)
        file_name_dict.update({'ctth': ctth_files})

    file_name_dict['cpp'] = get_pps_file(imager_file, AM_PATHS, values, 'cpp_file', 'cpp_dir', 
                                         FailIfRequestedAndMissing=True)                         
    file_name_dict['sunsatangles'] = get_pps_file(imager_file, AM_PATHS, values, 'sunsatangles_file', 
                                                  'sunsatangles_dir', 
                                                  FailIfRequestedAndMissing=True)
    file_name_dict['physiography'] = get_pps_file(imager_file, AM_PATHS, values, 'physiography_file', 
                                                  'physiography_dir', 
                                                  FailIfRequestedAndMissing=True)
    file_name_dict['r37'] = get_pps_file(imager_file, AM_PATHS, values, 'r37_file', 'r37_dir', 
                                         FailIfRequestedAndMissing=True)

    file_name_dict['nwp_tsur'] = get_pps_file(imager_file, AM_PATHS, values, 
                                              'nwp_tsur_file', 'nwp_nwp_dir', 
                                              FailIfRequestedAndMissing=True)
    file_name_dict['emis'] = get_pps_file(imager_file, AM_PATHS, values, 
                                          'emis_file', 'emis_dir', 
                                          FailIfRequestedAndMissing=True)

    file_name_dict['seaice'] = get_pps_file(imager_file, AM_PATHS, values, 
                                            'seaice_file', 'seaice_dir', 
                                            FailIfRequestedAndMissing=True)
    #Textures and Thresholds:
    file_name_dict['text_t11'] = get_pps_file(imager_file, AM_PATHS, values, 
                                              'text_t11_file', 'text_dir', 
                                              FailIfRequestedAndMissing=True)
    file_name_dict['thr_t11ts'] = get_pps_file(imager_file, AM_PATHS, values, 
                                               'thr_t11ts_file', 'thr_dir', 
                                               FailIfRequestedAndMissing=True)
    ####More NWP, textures and thresholds for v2014:
    for nwp_file in ['nwp_t500','nwp_t700',
                     'nwp_t850','nwp_t950', 'nwp_ciwv', 'nwp_ttro']:   
        file_name_dict[nwp_file] = get_pps_file(imager_file, AM_PATHS, values, 
                                                 nwp_file+'_file', 'nwp_nwp_dir', 
                                                 OnlyPrintInDebugMode=True)
    for text_file in ['text_r06',  'text_t37t12', 'text_t37']:
        file_name_dict[text_file] = get_pps_file(imager_file, AM_PATHS, values, 
                                                 text_file+'_file', 'text_dir', 
                                                 OnlyPrintInDebugMode=True)

    for thr_file in ['thr_t11ts_inv', 'thr_t11t37_inv', 
                     'thr_t37t12_inv', 'thr_t11t12_inv', 
                     'thr_t11ts', 'thr_t11t37', 'thr_t37t12', 'thr_t11t12',
                     'thr_r09', 'thr_r06', 'thr_t85t11_inv', 'thr_t85t11']:
        file_name_dict[thr_file] = get_pps_file(imager_file, AM_PATHS, values, 
                                                thr_file+'_file', 'thr_dir', 
                                                OnlyPrintInDebugMode=True)   
    file_name_dict['nwp_segments'] = get_pps_file(imager_file, AM_PATHS,
                                                  values, 
                                                  'segment_file', 'segment_dir') 
    ######################################## 
    ppsfiles = ppsFiles(file_name_dict)
    return  ppsfiles

def get_cloudsat_matchups(cloudsat_files, cloudsat_files_lwp, cloudproducts, SETTINGS):
    """
    Read Cloudsat data and match with the given PPS data.
    """
    cloudsat_lwp = None
    cloudsat = None
    if cloudsat_files is not None:   
        logger.debug("Reading cloudsat for type GEOPROF.")
        cloudsat = reshapeCloudsat(cloudsat_files, cloudproducts, SETTINGS)
    if cloudsat_files_lwp is not None: 
        logger.debug("Reading cloudsat for type CWC-RVOD.")
        cloudsat_lwp = reshapeCloudsat(cloudsat_files_lwp, cloudproducts, SETTINGS)    
    if cloudsat is not None and cloudsat_lwp is not None:
        logger.info("Merging CloudSat GEOPROF and CWC-RVOD data to one object")
        cloudsat = merge_cloudsat(cloudsat, cloudsat_lwp)
    elif cloudsat is None:
        cloudsat = cloudsat_lwp        
    logger.debug("Matching CloudSat with imager")
    cl_matchup = match_cloudsat_imager(cloudsat, cloudproducts,  SETTINGS)
    return cl_matchup

def get_iss_matchups(iss_files, cloudproducts,  SETTINGS):
    """
    Read Iss data and match with the given PPS data.
    """
    iss = reshape_iss(iss_files, cloudproducts,  SETTINGS)
    cl_matchup = match_iss_imager(iss, cloudproducts,  SETTINGS)
    return cl_matchup

def get_amsr_matchups(amsr_files, cloudproducts,  SETTINGS):
    """
    Read Amsr data and match with the given PPS data.
    """
    amsr = reshape_amsr(amsr_files, cloudproducts,  SETTINGS)
    am_matchup = match_amsr_imager(amsr, cloudproducts,  SETTINGS)
    return am_matchup

def get_synop_matchups(synop_files, cloudproducts,  SETTINGS):
    """
    Read Synop data and match with the given PPS data.
    """
    synop = reshape_synop(synop_files, cloudproducts,  SETTINGS)
    synop_matchup = match_synop_imager(synop, cloudproducts,  SETTINGS)
    return synop_matchup

def get_mora_matchups(mora_files, cloudproducts, SETTINGS):
    """
    Read Mora data and match with the given PPS data.
    """
    mora = reshape_mora(mora_files, cloudproducts,  SETTINGS)
    mora_matchup = match_mora_imager(mora, cloudproducts,  SETTINGS)
    return mora_matchup



def get_calipso_matchups(calipso_files, values, 
                         cloudproducts, 
                         AM_PATHS, SETTINGS, 
                         cafiles1km=None, cafiles5km=None, 
                         cafiles5km_aerosol=None):
    """
    Read Calipso data and match with the given PPS data.
    """
    #First remove files clearly outside time limit from the lists
    #When combinating 5km and 1km data some expensive calculations are done
    #before cutting the data that fits the time condition.
    #It is unnessecary to do this for files where no-pixel will match!
    calipso_files = discard_calipso_files_outside_time_range(
        calipso_files, cloudproducts, values, SETTINGS)
    if cafiles1km is not None:
        cafiles1km = discard_calipso_files_outside_time_range(
            cafiles1km, cloudproducts, values, SETTINGS, res=1)
    if cafiles5km is not None:
        cafiles5km = discard_calipso_files_outside_time_range(
            cafiles5km, cloudproducts, values, SETTINGS, res=5)
    if cafiles5km_aerosol is not None:
        cafiles5km_aerosol = discard_calipso_files_outside_time_range(
            cafiles5km_aerosol, cloudproducts, values, SETTINGS, res=5, ALAY=True)
    CALIPSO_version = 4
    if len(calipso_files)>0 and 'V3' in os.path.basename(calipso_files[0]):
        logger.info("CALIPSO version 3 data!")
        CALIPSO_version = 3

    calipso  = reshape_calipso(calipso_files)
    #find time breakpoints, but don't cut the data yet ..
    startBreak, endBreak = find_break_points(calipso,  cloudproducts, SETTINGS)
    if cafiles1km is not None and CALIPSO_version == 3 and config.RESOLUTION == 5:
        #RESOLUTION 5km also have 1km data
        logger.info("Calipso version 3 data used and old 1 km restore method!")
        calipso1km = reshape_calipso(cafiles1km, res=1)
        calipso5km = calipso
        calipso = add_1km_to_5km(calipso1km, calipso5km)
        calipso = time_reshape_calipso(calipso, 
                                       start_break=startBreak, end_break=endBreak) 
        calipso = total_and_top_layer_optical_depth_5km(calipso, SETTINGS, resolution=5)

    elif cafiles5km is not None and CALIPSO_version == 4  and config.RESOLUTION == 1:
        #RESOLUTION 1km also have 5km data calipso version 4
        logger.info("Calipso version 4, single shot fraction and "
                    "old 5km restored optical depth method used!")
        calipso1km = calipso
        calipso5km = reshape_calipso(cafiles5km, res=5)
        calipso5km = add_singleshot_to5km(calipso5km, SETTINGS) 
        calipso5km = total_and_top_layer_optical_depth_5km(calipso5km, SETTINGS, resolution=5)
        calipso1km = add_5km_variables_to_1km(calipso1km, calipso5km, CALIPSO_version)
        calipso = calipso1km.extract_elements(starti=startBreak, 
                                              endi=endBreak) 

    elif cafiles5km is not None and CALIPSO_version == 3 and config.RESOLUTION == 1:
        #RESOLUTION 1km also have 5km data calipso version 3
        logger.info("Calipso version 3 data used and old 5 km restored optical depth method!")
        calipso1km = calipso
        calipso5km = reshape_calipso(cafiles5km, res=5)
        calipso5km = total_and_top_layer_optical_depth_5km(calipso5km, SETTINGS, resolution=5)
        calipso1km = add_5km_variables_to_1km(calipso1km, calipso5km, CALIPSO_version)
        calipso = calipso1km.extract_elements(starti=startBreak, 
                                              endi=endBreak) 
    elif CALIPSO_version == 4 and config.RESOLUTION == 5 and SETTINGS['ALSO_USE_SINGLE_SHOT_CLOUD_CLEARED']:

        #RESOLUTION exclusively 5km data but additional clouds taken from 330 m single shot resolution
        logger.info("Calipso version 4 data used and new single shot restore method!")
        #calipso5km = reshape_calipso(cafiles5km, res=5)
        calipso5km  = reshape_calipso(calipso_files)
        calipso = add_singleshot_to5km(calipso5km, SETTINGS) 

        calipso = calipso.extract_elements(starti=startBreak, 
                                           endi=endBreak) 
        calipso = total_and_top_layer_optical_depth_5km(calipso, SETTINGS, resolution=5)
    elif CALIPSO_version == 4 and config.RESOLUTION == 5 and SETTINGS['ALSO_USE_1KM_FILES']:

        #RESOLUTION exclusively 5km data but additional clouds taken from 1 km data
        logger.info("Calipso version 4 data used but old method combining 1 km and 5 km data!")
        #calipso5km = reshape_calipso(cafiles5km, res=5)
        calipso5km  = reshape_calipso(calipso_files)
        calipso1km = reshape_calipso(cafiles1km, res=1)
        calipso = add_1km_to_5km(calipso1km, calipso5km) 
        calipso = calipso.extract_elements(starti=startBreak, 
                                           endi=endBreak) 
        calipso = total_and_top_layer_optical_depth_5km(calipso, SETTINGS, resolution=5)
    else:
        logger.warning("Old metod, only one resolution used, expect bad results!")
        calipso = calipso.extract_elements(starti=startBreak, 
                                           endi=endBreak) 
        if config.RESOLUTION == 5:
            calipso = total_and_top_layer_optical_depth_5km(calipso, SETTINGS, resolution=5)

    #aerosol-data
    calipso_aerosol = None
    if cafiles5km_aerosol is not None:
        calipso5km_aerosol = reshape_calipso(cafiles5km_aerosol, res=5, ALAY=True)
        if config.RESOLUTION == 1:
            calipso_aerosol = adjust_5km_to_1km_resolution(calipso5km_aerosol)
        elif config.RESOLUTION == 5:
            calipso_aerosol = calipso5km_aerosol
        calipso_aerosol = calipso_aerosol.extract_elements(starti=startBreak, 
                                                           endi=endBreak) 
    # free some memory    
    calipso1km = None
    calipso5km = None
        
    logger.debug("Matching CALIPSO with imager")
    ca_matchup = match_calipso_imager(values, calipso, calipso_aerosol,
                                      cloudproducts, SETTINGS)
    return ca_matchup
def read_cloud_cci(imager_file):
    from imager_cloud_products.read_cloudproducts_cci import cci_read_all
    return cci_read_all(imager_file)
def read_cloud_oca(imager_file):
    from imager_cloud_products.read_cloudproducts_oca import oca_read_all
    return oca_read_all(imager_file)
def read_cloud_maia(imager_file):
    from imager_cloud_products.read_cloudproducts_maia import maia_read_all
    return maia_read_all(imager_file)

def read_cloud_patmosx(imager_file, cross, SETTINGS):
    from imager_cloud_products.read_cloudproducts_patmosx import patmosx_read_all
    return patmosx_read_all(imager_file, cross, SETTINGS)


def read_pps_data(pps_files, imager_file, SETTINGS):
    from imager_cloud_products.read_cloudproducts_and_nwp_pps import pps_read_all
    return pps_read_all(pps_files, imager_file, SETTINGS)

def get_additional_calipso_files_if_requested(calipso_files, SETTINGS):
    import glob
    calipso5km = None
    calipso1km = None
    calipso5km_aerosol=None

    if config.RESOLUTION == 5:
        if SETTINGS['ALSO_USE_1KM_FILES'] == True:
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
            calipso1km = sorted(require_h5(calipso1km, SETTINGS))
            calipso_files = sorted(calipso5km)
            calipso5km = None

            if len(calipso1km) == 0:
                raise MatchupError("Did not find any matching 1km calipso files")

            if len(calipso_files) != len(calipso1km):
                raise MatchupError("Inconsistent number of calipso files...\n" + 
                                   "\tlen(calipso_files) = %d\n" % len(calipso_files) + 
                                   "\tlen(calipso1km) = %d" % len(calipso1km))
    if config.RESOLUTION == 1:
        if SETTINGS['ALSO_USE_5KM_FILES'] == True:
            logger.info("Search for CALIPSO 5km data too")
            calipso5km = []
            calipso1km = []
            for file1km in calipso_files:
                file5km = file1km.replace('/1km/', '/5km/').\
                          replace('01kmCLay', '05kmCLay').\
                          replace('-Standard-V4-10.','*').\
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
            calipso5km = sorted(require_h5(calipso5km, SETTINGS))
            calipso_files = sorted(calipso1km)
            calipso1km = None
            if len(calipso5km) == 0:
                raise MatchupError("Did not find any matching 5km calipso files")
            if len(calipso_files) != len(calipso5km):
                raise MatchupError("Inconsistent number of calipso files...\n" + 
                                   "\tlen(calipso_files) = %d\n" % len(calipso_files) + 
                                   "\tlen(calipso1km) = %d" % len(calipso5km))
    if SETTINGS['MATCH_AEROSOL_CALIPSO']:
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
        logger.debug("found these aerosol files:") 
        logger.debug("\n".join(calipso5km_aerosol))             
        calipso5km_aerosol = sorted(require_h5(calipso5km_aerosol, SETTINGS))
        logger.debug("found these aerosol files:") 
        logger.debug("\n".join(calipso5km_aerosol)) 
    return calipso5km, calipso1km, calipso5km_aerosol              


def add_additional_clousat_calipso_index_vars(match_clsat, match_calipso):
    #add cloudsat modisflag to calipso obj
    
    if match_clsat is not None and match_calipso is not None:
        match_calipso.calipso.cal_modis_cflag = None
        #map cloudsat to calipso and the other way around!
        from utils.match import match_lonlat
        source = (match_clsat.cloudsat.longitude.astype(np.float64).reshape(-1,1), 
                  match_clsat.cloudsat.latitude.astype(np.float64).reshape(-1,1))
        target = (match_calipso.calipso.longitude.astype(np.float64).reshape(-1,1), 
                  match_calipso.calipso.latitude.astype(np.float64).reshape(-1,1))
        mapper, dummy = match_lonlat(source, target, radius_of_influence=1000, n_neighbours=1)
        match_calipso.calipso.cloudsat_index = mapper.rows.filled(config.NODATA).ravel()
        target = (match_clsat.cloudsat.longitude.astype(np.float64).reshape(-1,1), 
                  match_clsat.cloudsat.latitude.astype(np.float64).reshape(-1,1))
        source = (match_calipso.calipso.longitude.astype(np.float64).reshape(-1,1), 
                  match_calipso.calipso.latitude.astype(np.float64).reshape(-1,1))
        mapper, dummy = match_lonlat(source, target, radius_of_influence=1000, n_neighbours=1)
        match_clsat.cloudsat.calipso_index = mapper.rows.filled(config.NODATA).ravel()

        # Transfer CloudSat MODIS cloud flag to CALIPSO representation
        index = match_calipso.calipso.cloudsat_index.copy()
        index[index<0] = 0
        match_calipso.calipso.cal_modis_cflag = np.where(
            match_calipso.calipso.cloudsat_index>=0, 
            match_clsat.cloudsat.MODIS_cloud_flag[index],
            -9)
        
    if match_clsat is not None and match_calipso is not None:
        index = match_clsat.cloudsat.calipso_index.copy()
        index[index<0] = 0
        for var_2d_name in ['feature_classification_flags', 
                           'layer_base_altitude', 
                           'layer_top_altitude',
                           'feature_optical_depth_532',
                           'feature_optical_depth_532_5km']:
            if hasattr(match_calipso.calipso, var_2d_name):
                data_calipso = getattr(match_calipso.calipso, var_2d_name)
                if data_calipso is None:
                    continue
                temp_data = np.array([np.where(
                    match_clsat.cloudsat.calipso_index>=0, 
                    data_calipso[index,i], -9) for i in range(data_calipso.shape[1])])
                setattr(match_clsat.cloudsat,  'calipso_{:s}'.format(var_2d_name), temp_data.transpose()) 

        for var_1d_name in ['column_optical_depth_tropospheric_aerosols_532',
                           'column_optical_depth_tropospheric_aerosols_532_5km',
                           'column_optical_depth_aerosols_532',
                           'column_optical_depth_aerosols_532_5km']:

            if hasattr(match_calipso.calipso, var_1d_name):
                data_calipso = getattr(match_calipso.calipso, var_1d_name)
                if data_calipso is None:
                    continue
                temp_data = np.where(match_clsat.cloudsat.calipso_index>=0, data_calipso[index], -9) 
                setattr(match_clsat.cloudsat,  'calipso_{:s}'.format(var_1d_name), temp_data) 
            
        """
        match_clsat.cloudsat.calipso_feature_classification_flags= np.where(
            match_clsat.cloudsat.calipso_index>=0,
            match_calipso.calipso.feature_classification_flags[index,0],
            -9)
        # first bas layer use height not pressure as cloudsat uses height
        match_clsat.cloudsat.calipso_layer_base_altitude = np.where(
            match_clsat.cloudsat.calipso_index>=0,
            match_calipso.calipso.layer_base_altitude[index,0],
            -9)
        for layer in range(1,10):
            match_clsat.cloudsat.calipso_layer_base_altitude = np.where(
                np.logical_and(np.logical_and(match_clsat.cloudsat.calipso_index>=0,
                                              match_calipso.calipso.layer_base_altitude[index,layer]>0),
                               match_calipso.calipso.layer_base_altitude[index,layer]< 
                               match_calipso.calipso.layer_base_altitude[index,layer-1]),
                
                match_calipso.calipso.layer_base_altitude[index,layer],
            match_clsat.cloudsat.calipso_layer_base_altitude)
        match_clsat.cloudsat.calipso_layer_base_altitude[match_clsat.cloudsat.calipso_layer_base_altitude<-999] = -9.0
        # first bas layer use height not pressure as cloudsat uses height
        match_clsat.cloudsat.calipso_layer_top_altitude = np.where(
            match_clsat.cloudsat.calipso_index>=0,
            match_calipso.calipso.layer_top_altitude[index,0],
            -9)
        for layer in range(1,10):
            match_clsat.cloudsat.calipso_layer_top_altitude = np.where(
                np.logical_and(np.logical_and(match_clsat.cloudsat.calipso_index>=0,
                                              match_calipso.calipso.layer_top_altitude[index,layer]>0),
                               match_calipso.calipso.layer_top_altitude[index,layer]> 
                               match_calipso.calipso.layer_top_altitude[index,layer-1]),
                
                match_calipso.calipso.layer_top_altitude[index,layer],
            match_clsat.cloudsat.calipso_layer_top_altitude)
        match_clsat.cloudsat.calipso_layer_top_altitude[match_clsat.cloudsat.calipso_layer_top_altitude<-999] = -9.0
    """   
    return match_clsat, match_calipso

def add_modis_lvl2_clousat_(match_clsat, match_calipso):
    #add modis lvl2 to cloudsat
    if match_clsat is not None and match_calipso is not None:
        match_clsat.modis.all_arrays["height"] = match_calipso.modis.all_arrays["height"][match_clsat.cloudsat.calipso_index]
        match_clsat.modis.all_arrays["temperature"] = match_calipso.modis.all_arrays["temperature"][match_clsat.cloudsat.calipso_index]
        match_clsat.modis.all_arrays["pressure"] = match_calipso.modis.all_arrays["pressure"][match_clsat.cloudsat.calipso_index]
        match_clsat.modis.all_arrays["cloud_emissivity"] = match_calipso.modis.all_arrays["cloud_emissivity"][match_clsat.cloudsat.calipso_index]
        match_clsat.modis.all_arrays["latitude_5km"] = match_calipso.modis.all_arrays["latitude_5km"][match_clsat.cloudsat.calipso_index]
        match_clsat.modis.all_arrays["longitude_5km"] = match_calipso.modis.all_arrays["longitude_5km"][match_clsat.cloudsat.calipso_index]
    return match_clsat


def add_elevation_corrected_imager_ctth(match_clsat, match_calipso, issObj, SETTINGS):
    ## Cloudsat ##
    if match_clsat is None or match_clsat.imager.ctth_height is None:
        pass
    elif  match_clsat.imager.imager_ctth_m_above_seasurface is None:
        # First make sure that PPS cloud top heights are converted to height
        # above sea level just as CloudSat height are defined. Use
        # corresponding DEM data.
        elevation = np.where(np.less_equal(match_clsat.cloudsat.elevation,0),
                             0,match_clsat.cloudsat.elevation)		
        num_csat_data_ok = len(match_clsat.cloudsat.elevation)
        logger.debug("Length of CLOUDSAT array: %d", num_csat_data_ok )
        imager_ctth_m_above_seasurface = match_clsat.imager.ctth_height.copy()
        #import pdb;pdb.set_trace()
        if SETTINGS["CCI_CLOUD_VALIDATION"] or SETTINGS["PATMOSX_VALIDATION"]: 
            #ctth already relative mean sea level
            imager_ctth_m_above_seasurface = match_calipso.imager.ctth_height
        else: #ctth relative topography
            got_height = imager_ctth_m_above_seasurface>=0                    
            imager_ctth_m_above_seasurface[got_height] += elevation[got_height]*1.0
        match_clsat.imager.imager_ctth_m_above_seasurface = imager_ctth_m_above_seasurface
        if num_csat_data_ok == 0:
            logger.info("Processing stopped: Zero lenght of matching arrays!")
            raise ProcessingError("Zero lenght of matching arrays!")
    ## Calipso ##        
    # First make sure that PPS cloud top heights are converted to height above sea level
    # just as CALIPSO heights are defined. Use corresponding DEM data.
    if match_calipso is None or match_calipso.imager.ctth_height is None:
        pass
    elif  match_calipso.imager.imager_ctth_m_above_seasurface is None:
        cal_elevation = np.where(np.less_equal(match_calipso.calipso.elevation,0),
                                 0,match_calipso.calipso.elevation)
        num_cal_data_ok = len(match_calipso.calipso.elevation)
        logger.debug("Length of CALIOP array: %d", num_cal_data_ok)
        imager_ctth_m_above_seasurface = match_calipso.imager.ctth_height.copy()
        logger.debug("CCI_CLOUD_VALIDATION %s", str(SETTINGS["CCI_CLOUD_VALIDATION"]))
        if SETTINGS["CCI_CLOUD_VALIDATION"]: 
            #ctth relative mean sea level
            imager_ctth_m_above_seasurface = match_calipso.imager.ctth_height
        else: #ctth relative topography
            got_height = imager_ctth_m_above_seasurface>=0                    
            imager_ctth_m_above_seasurface[got_height] += cal_elevation[got_height]*1.0
        match_calipso.imager.imager_ctth_m_above_seasurface = imager_ctth_m_above_seasurface
    if issObj is None or issObj.imager.ctth_height is None:
        pass
    elif  issObj.imager.imager_ctth_m_above_seasurface is None:
        iss_elevation = np.where(np.less_equal(issObj.iss.elevation,0),
                                 0,issObj.iss.elevation)
        num_iss_data_ok = len(issObj.iss.elevation)
        logger.info("Length of ISS array: %d", num_iss_data_ok)
        imager_ctth_m_above_seasurface = issObj.imager.ctth_height.copy()
        if SETTINGS["CCI_CLOUD_VALIDATION"]: 
            #ctth relative mean sea level
            pass
        else: #ctth relative topography
            got_height = imager_ctth_m_above_seasurface>=0                    
            imager_ctth_m_above_seasurface[got_height] += iss_elevation[got_height]*1.0
        issObj.imager.imager_ctth_m_above_seasurface = imager_ctth_m_above_seasurface
    return match_clsat, match_calipso, issObj


def get_matchups_from_data(cross, AM_PATHS, SETTINGS):
    """
    Retrieve Cloudsat- and Calipso-IMAGER matchups from Cloudsat, Calipso, and
    PPS files.
    """
    #STEP 1 get imager files
    if SETTINGS['PPS_VALIDATION']:
        imager_file, tobj = find_radiance_file(cross, AM_PATHS)
        pps_files = find_files_from_imager(imager_file, AM_PATHS, SETTINGS) 
    if SETTINGS['CCI_CLOUD_VALIDATION']:
        imager_file, tobj = find_cci_cloud_file(cross, AM_PATHS)
    if SETTINGS['MAIA_VALIDATION']:
        imager_file, tobj = find_maia_cloud_file(cross, AM_PATHS)
    if SETTINGS['PATMOSX_VALIDATION']:
        imager_file, tobj = find_patmosx_cloud_file(cross, AM_PATHS)
    if SETTINGS['OCA_VALIDATION']:
        imager_file, tobj = find_oca_cloud_file(cross, AM_PATHS)
    if not imager_file:
        raise MatchupError("No imager file found!\ncross = " + str(cross))
    values = get_satid_datetime_orbit_from_fname(imager_file, SETTINGS, cross)
    date_time = values["date_time"]

    #Step 2 get truth satellite files 
    truth_files = {}
    for truth in ['cloudsat', 'amsr', 'iss', 'synop', 'mora', 'cloudsat_lwp', 'calipso']:
        truth_files[truth] = None
        if (SETTINGS[truth.replace("_lwp","").upper()+'_MATCHING'] and truth + '_file' in AM_PATHS.keys()):
            truth_files[truth] = find_truth_files(date_time, AM_PATHS, SETTINGS, values, truth=truth)
        elif not SETTINGS[truth.replace("_lwp","").upper()+'_MATCHING']:
            logger.info("NO {truth} File, {truth} matching not requested "
                        "{truth}_MATCHING=False".format(truth=truth))  
        elif  truth + '_file' not in AM_PATHS.keys():
            logger.info("NO {truth}_file in atrain_match.cfg".format(truth=truth.lower()))
    #CALIPSO get some extra files:
    if truth_files['calipso'] is not None:
        extra_files = get_additional_calipso_files_if_requested( truth_files['calipso'], SETTINGS)
        calipso5km, calipso1km, calipso5km_aerosol  = extra_files
    if (all(truth_files_i is None for truth_files_i in truth_files.values())): 
        raise MatchupError(
            "Couldn't find any matching CALIPSO/CLoudSat/ISS data")

   
    #STEP 3 Read imager data:    
    if (SETTINGS['PPS_VALIDATION'] ):
        cloudproducts =read_pps_data(pps_files, imager_file, SETTINGS)
        if os.path.isfile(SETTINGS['CNN_PCKL_PATH']):
            from utils.pps_prototyping_util import add_cnn_features_full
            imager_obj.cnn_dict = add_cnn_features_full(cloudproducts.imager_channeldata, 
                                                       cloudproducts, 
                                                       SETTINGS)
    if (SETTINGS['CCI_CLOUD_VALIDATION']):
        cloudproducts = read_cloud_cci(imager_file) 
        cloudproducts.satellite = values["satellite"]
    if (SETTINGS['MAIA_VALIDATION']):
        cloudproducts = read_cloud_maia(imager_file)  
        cloudproducts.satellite = values["satellite"]
    if (SETTINGS['PATMOSX_VALIDATION']):
        cloudproducts = read_cloud_patmosx(imager_file, cross, SETTINGS)
        cloudproducts.satellite = values["satellite"]
    if (SETTINGS['OCA_VALIDATION']):
        cloudproducts = read_cloud_oca(imager_file)
        cloudproducts.satellite = values["satellite"]

    #STEP 4 get matchups 
    #ClloudSat
    cloudsat_matchup = None
    if (SETTINGS['PPS_VALIDATION'] and SETTINGS['CLOUDSAT_MATCHING'] and 
        truth_files['cloudsat'] is not None):
        logger.info("Read CLOUDSAT data")
        cloudsat_matchup = get_cloudsat_matchups(truth_files['cloudsat'],
                                                 truth_files['cloudsat_lwp'],
                                                 cloudproducts, SETTINGS) 
    #ISS:  
    iss_matchup = None
    if (SETTINGS['PPS_VALIDATION'] and SETTINGS['ISS_MATCHING'] and  
        truth_files['iss'] is not None):
        logger.info("Read ISS data")
        iss_matchup = get_iss_matchups(truth_files['iss'],
                                       cloudproducts, SETTINGS)
    #AMSR
    amsr_matchup = None
    if (SETTINGS['PPS_VALIDATION'] and SETTINGS['AMSR_MATCHING'] and 
        truth_files['amsr'] is not None):
        logger.info("Read AMSR data")
        amsr_matchup = get_amsr_matchups(truth_files['amsr'],
                                         cloudproducts, SETTINGS)
    #SYNOP
    synop_matchup = None
    if (SETTINGS['PPS_VALIDATION'] and SETTINGS['SYNOP_MATCHING'] and 
        truth_files['synop'] is not None):
        logger.info("Read SYNOP data")
        synop_matchup = get_synop_matchups(truth_files['synop'], 
                                           cloudproducts, SETTINGS)
    #MORA
    mora_matchup = None
    if (SETTINGS['PPS_VALIDATION'] and SETTINGS['MORA_MATCHING'] and truth_files['mora'] is not None):
        logger.info("Read MORA data")
        mora_matchup = get_mora_matchups(truth_files['mora'],
                                         cloudproducts, SETTINGS)
    #CALIPSO:
    calipso_matchup = None
    if SETTINGS['CALIPSO_MATCHING'] and truth_files['calipso'] is not None:
        logger.info("Read CALIPSO data")        
        calipso_matchup= get_calipso_matchups(truth_files['calipso'],
                                              values,
                                              cloudproducts,
                                              AM_PATHS, SETTINGS,
                                              calipso1km, calipso5km, calipso5km_aerosol)
 
    if calipso_matchup is None and SETTINGS['CALIPSO_REQUIRED']:
        raise MatchupError("No matches with CALIPSO.")
    elif cloudsat_matchup is None and SETTINGS['CLOUDSAT_REQUIRED']:
        raise MatchupError("No matches with CLOUSDAT.")
    elif iss_matchup is None and SETTINGS['ISS_REQUIRED']:
        raise MatchupError("No matches with ISS.")
    elif amsr_matchup is None and SETTINGS['AMSR_REQUIRED']:
        raise MatchupError("No matches with AMSR.")
    elif synop_matchup is None and SETTINGS['SYNOP_REQUIRED']:
        raise MatchupError("No matches with SYNOP.")
    elif mora_matchup is None and SETTINGS['MORA_REQUIRED']:
        raise MatchupError("No matches with MORA.")
    elif (calipso_matchup is None and cloudsat_matchup is None and iss_matchup is None 
          and amsr_matchup is None and synop_matchup is None 
          and mora_matchup is None):
        raise MatchupError("No matches with any truth.")


    # Get satellite name, time, and orbit number from imager_file
    date_time = values["date_time"]
    #basename = '_'.join(os.path.basename(imager_file).split('_')[:4])
    values = get_satid_datetime_orbit_from_fname(imager_file, SETTINGS, cross)
    basename = values["basename"]
    rematched_path = date_time.strftime(AM_PATHS['reshape_dir'].format(
        val_dir=config._validation_results_dir,
        satellite=values["satellite"],
        resolution=str(config.RESOLUTION),
    ))
    rematched_file = date_time.strftime(AM_PATHS['reshape_file'].format(
        satellite=values["satellite"],
        orbit=values["orbit"],
        resolution=str(config.RESOLUTION),
        instrument=INSTRUMENT.get(values["satellite"],'avhrr'),
        atrain_datatype="atrain_datatype",
        extrai=values.get("extrai","")
    ))
    rematched_file_base = rematched_path + rematched_file
    
    # Create directories if they don't exist yet
    if not os.path.exists(os.path.dirname(rematched_path)):
        logger.info("Creating dir %s:", rematched_path)
        os.makedirs(os.path.dirname(rematched_path))

    for matchup, name in zip([cloudsat_matchup, iss_matchup, amsr_matchup, 
                              synop_matchup, mora_matchup, calipso_matchup],
                             ['CloudSat', 'ISS', 'AMSR-E',
                              'SYNOP', 'MORA', 'CALIPSO']):
        if matchup is None:
            continue
        #add modis lvl2    
        if SETTINGS['MATCH_MODIS_LVL2']:
            from imager_cloud_products.read_modis_products import add_modis_06  
            if matchup.imager_instrument in ['modis']:
                matchup = add_modis_06(matchup, imager_file, AM_PATHS) 
        if SETTINGS['ADD_NWP']:
            import pps_nwp
            from libs.extract_imager_along_track import _interpolate_height_and_temperature_from_pressure
            nwp_file = find_closest_nwp_file(cloudproducts, AM_PATHS, 
                                             values, SETTINGS)
            logger.debug(nwp_file)
            if pps_nwp is None:
                continue
            gribfile = pps_nwp.GRIBFile(nwp_file, (matchup.imager.longitude,
                                                   matchup.imager.latitude))
            setattr(matchup.imager, "nwp_height", 
                    gribfile.get_gh_vertical()[0,:,:].astype(np.float32).transpose())
            setattr(matchup.imager, "nwp_surface_h",  
                    gribfile.get_gh_surface()[:].astype(np.float32))
            setattr(matchup.imager, "nwp_temperature",  
                    gribfile.get_t_vertical()[0,:,:].astype(np.float32).transpose())
            setattr(matchup.imager, "nwp_h2m",  
                    gribfile.get_h_2meter()[:].astype(np.float32))
            setattr(matchup.imager, "nwp_t2m",  
                    gribfile.get_t_2meter()[:].astype(np.float32))
            setattr(matchup.imager, "nwp_u10m",  
                    gribfile.get_u_10meter()[:].astype(np.float32))
            setattr(matchup.imager, "nwp_v10m",  
                    gribfile.get_v_10meter()[:].astype(np.float32))
            #get pressure variables in hPs
            field = gribfile.get_p_vertical()
            if field.units =='Pa':
                field = 0.01*field[:]
                field.units = 'hPa'
            matchup.imager.nwp_pressure = field[0,:,:].astype(np.float32).transpose()
            field = gribfile.get_p_surface()
            if field.units =='Pa':
                field = 0.01*field[:]
                field.units = 'hPa'
            matchup.imager.nwp_psur = field[:].astype(np.float32).transpose()
            data = _interpolate_height_and_temperature_from_pressure(matchup.imager, 440)
            setattr(matchup.imager, 'nwp_h440', data)
            data = _interpolate_height_and_temperature_from_pressure(matchup.imager, 680)
            setattr(matchup.imager, 'nwp_h680', data)

    #add additional vars to cloudsat and calipso objects and print them to file:
    cloudsat_matchup, calipso_matchup = add_additional_clousat_calipso_index_vars(cloudsat_matchup, calipso_matchup)
    cloudsat_matchup, calipso_matchup, iss_matchup = add_elevation_corrected_imager_ctth(cloudsat_matchup, calipso_matchup, iss_matchup, SETTINGS)

 
    #imager_name
    imager_obj_name = 'pps'
    if SETTINGS['CCI_CLOUD_VALIDATION']:
        imager_obj_name = 'cci'
    if SETTINGS['OCA_VALIDATION']:
        imager_obj_name = 'oca'
    if SETTINGS['MAIA_VALIDATION']:
        imager_obj_name = 'maia'
    if SETTINGS['PATMOSX_VALIDATION']:
        imager_obj_name = 'patmosx'

    #write matchups
    for matchup, name in zip([cloudsat_matchup, iss_matchup, amsr_matchup, 
                              synop_matchup, mora_matchup, calipso_matchup],
                             ['CloudSat', 'ISS', 'AMSR-E',
                              'SYNOP', 'MORA', 'CALIPSO']):
        if matchup is None:
            logger.debug("No {:s} Match File created".format(name))
        else:    
            truth_sat = matchup.truth_sat
            match_file = rematched_file_base.replace(
                'atrain_datatype', truth_sat)
            write_truth_imager_match_obj(match_file, matchup,
                                     SETTINGS, 
                                     imager_obj_name = imager_obj_name)

    #no longer return the matchup data?
    aux_obj = None
    return {'cloudsat': cloudsat_matchup, 
            'calipso': calipso_matchup, 
            'iss': iss_matchup, 
            'amsr': amsr_matchup,
            'synop': synop_matchup, 
            'mora': mora_matchup, 
            'basename': basename, 
            'values':values}

def check_if_got_all_match_files(cross, AM_PATHS, SETTINGS):
    values = {}
    values["satellite"] = cross.satellite1.lower()
    for truth in ['cloudsat', 'amsr', 'iss', 'synop', 'mora', 'calipso']:
        if not SETTINGS[truth.upper() + '_MATCHING']:
            logger.info(
                "{truth} matching turned off {truth}_MATCHING]=False.".format(
                    truth=truth.upper()))
        else:    
            values["atrain_sat"] = truth
            values["atrain_datatype"] = truth
            match_file, date_time = find_imager_file(
                cross, 
                AM_PATHS['reshape_dir'], 
                AM_PATHS['reshape_file'], 
                values=values)
            if match_file is None:
                return False
    # If no reshaped file where missing, return true:            
    return True
                               
def run(cross, AM_PATHS, SETTINGS, reprocess=False):
    """
    The main work horse.    
    """    
    logger.info("Case: %s", str(cross))
    sensor = INSTRUMENT.get(cross.satellite1.lower(), 'imager')
    # Match the data that we need:
    if reprocess or not check_if_got_all_match_files(cross, AM_PATHS, SETTINGS):
        matchup_results = get_matchups_from_data(cross, AM_PATHS, SETTINGS)

 
