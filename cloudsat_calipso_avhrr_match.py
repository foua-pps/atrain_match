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
                    PLOT_ONLY_PNG,
                    CCI_CLOUD_VALIDATION,
                    PPS_VALIDATION,
                    ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS,
                    USE_5KM_FILES_TO_FILTER_CALIPSO_DATA,
                    PPS_FORMAT_2012_OR_EARLIER,
                    RESOLUTION)


from pps_basic_configure import *
from pps_error_messages import write_log

#import config
from common import attach_subdir_from_config, MatchupError

from cloudsat_calipso_avhrr_statistics import *
#from trajectory_plot import plotSatelliteTrajectory
from trajectory_plotting import plotSatelliteTrajectory
from cloudsat_calipso_avhrr_prepare import *

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
from calipso import  (use5km_find_detection_height_and_total_optical_thickness_faster, 
                      add1kmTo5km,
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
              'eos2': 'modis'} 
from datetime import datetime, timedelta
from glob import glob

class ppsFiles(object):
    def __init__(self, file_name_dict):
        self.cloudtype = None
        self.ctth = None
        self.cpp = None
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

class NWPObj(object):
    def __init__(self, array_dict):
        self.surftemp = None
        self.t500 = None
        self.t700 = None
        self.t850 = None
        self.t950 = None
        self.ttro = None
        self.ciwv = None
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
        self.emis1 = None
        self.emis8 = None
        self.emis9 = None
        self.__dict__.update(array_dict) 

test = 0

def readCpp(filename, cpp_type):
    import h5py #@UnresolvedImport
    h5file = h5py.File(filename, 'r')
    if cpp_type in h5file.keys():
        value = h5file[cpp_type].value
        gain = h5file[cpp_type].attrs['gain']
        intercept = h5file[cpp_type].attrs['intercept']
        nodat = h5file[cpp_type].attrs['no_data_value']
        product = np.where(value != nodat,value * gain + intercept, value)   
    h5file.close()
    return product


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

def get_satid_datetime_orbit_from_fname(avhrr_filename,as_oldstyle=False):
    #import runutils
    #satname, _datetime, orbit = runutils.parse_scene(avhrr_filename)
    #returnd orbit as int, loosing leeding zeros, use %05d to get it right.
    # Get satellite name, time, and orbit number from avhrr_file
    if PPS_VALIDATION  and not PPS_FORMAT_2012_OR_EARLIER and not as_oldstyle:
        sl_ = os.path.basename(avhrr_filename).split('_')
        date_time= datetime.strptime(sl_[5], '%Y%m%dT%H%M%S%fZ')
        values= {"satellite": sl_[3],
                 "date_time": date_time,
                 "orbit": sl_[4],
                 #"date":"%04d%02d%02d"%(date_time.year,dat_time.month, date_time.day),
                 "date": date_time.strftime("%Y%m%d"),
                 "year":date_time.year,
                 "month":"%02d"%(date_time.month),  
                 "lines_lines": "*",
                 "time":date_time.strftime("%H%M"),
                 "ppsfilename":avhrr_filename}
        values['basename'] = values["satellite"] + "_" + values["date"] + "_" + values["time"] + "_" + values["orbit"]
    if PPS_VALIDATION and (PPS_FORMAT_2012_OR_EARLIER or as_oldstyle):
        sl_ = os.path.basename(avhrr_filename).split('_')
        date_time= datetime.strptime(sl_[1] + sl_[2], '%Y%m%d%H%M')
        values= {"satellite": sl_[0],
                 "date_time": date_time,
                 "orbit": sl_[3].split('.')[0],
                 "date":sl_[1],
                 "year":date_time.year,
                 "month":"%02d"%(date_time.month),  
                 #"lines_lines": sl_[5] + "_" + sl_[6],
                 "lines_lines": "*",
                 "time":sl_[2],
                 "ppsfilename":avhrr_filename}
        values['basename'] = values["satellite"] + "_" + values["date"] + "_" + values["time"] + "_" + values["orbit"]
    if CCI_CLOUD_VALIDATION:
         #avhrr_file = "20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
        sl_ = os.path.basename(avhrr_filename).split('-')
        date_time = datetime.strptime(sl_[0], '%Y%m%d%H%M%S')

        sat_id = "noaa18"#sl_[5]_lower,
        values= {"satellite": sat_id,
                 "date_time": date_time,
                 "orbit": "99999",
                 "date":date_time.strftime("%Y%m%d"),
                 "year":date_time.year,
                 "month":"%02d"%(date_time.month),    
                 "time":date_time.strftime("%H%M"),
                 #"basename":sat_id + "_" + date_time.strftime("%Y%m%d_%H%M_99999"),#"20080613002200-ESACCI",
                 "ccifilename":avhrr_filename,
                 "ppsfilename":None}
        values['basename'] = values["satellite"] + "_" + values["date"] + "_" + values["time"] + "_" + values["orbit"]


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
        write_log("INFO", "Do not use CALIPSO in CWC-RWOD mode, Continue")
        calipso_files = []
    else:
        #might need to geth this in before looking for matchups
        tdelta = timedelta(seconds = (config.SAT_ORBIT_DURATION + config.sec_timeThr))
        time_window = (tdelta, tdelta)
        calipso_files = find_calipso_files_inner(date_time, time_window, options, values)
        if len(calipso_files) > 1:
            write_log("INFO", "More than one Calipso file found within time window!")
        elif len(calipso_files) == 0:
            raise MatchupError("Couldn't find calipso matchup!")
        calipso_files = sorted(require_h5(calipso_files))
        calipso_basenames = [ os.path.basename(s) for s in calipso_files ]
        write_log("INFO", "Calipso files: " + str(calipso_basenames))
        return calipso_files

def find_cloudsat_files(date_time, options, values):
    #might need to geth this in before looking for matchups
    tdelta = timedelta(seconds = (config.SAT_ORBIT_DURATION + config.sec_timeThr))
    time_window = (tdelta, tdelta)
    cloudsat_files = find_cloudsat_files_inner(date_time, time_window, options, values)
    if len(cloudsat_files) > 1:
        write_log("INFO", "More than one Cloudsat file found within time window!")
    elif len(cloudsat_files) == 0:
        write_log("INFO", "No Cloudsat file found within time window!")
        #raise MatchupError("Couldn't find cloudsat matchup!")
    cloudsat_files = sorted(require_h5(cloudsat_files))
    cloudsat_basenames = [ os.path.basename(s) for s in cloudsat_files ]
    write_log("INFO", "Cloudsat files: " + str(cloudsat_basenames))
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
    if ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS:
        ddt2=timedelta(seconds=0)
    else:
        ddt2=ddt
    time_low = cross_time - ddt
    time_high = cross_time + ddt2
    write_log("INFO", "Searching for avhrr/viirs file with start time  between" 
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

def find_avhrr_file(cross, filedir_pattern, filename_pattern, values={}):
    (tlist, cross_time, cross_satellite) = get_time_list_and_cross_time(cross)
    time_window=cross.time_window
    write_log('INFO',"Time window: ", time_window)
    write_log('INFO',"Cross time: {cross_time}".format(cross_time=cross_time))
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
                write_log('INFO',"Found directory "
                          "{dirctory} ".format(dirctory=found_dir))
        else:
            if not found_dir and avhrr_dir not in checked_dir.keys():
                checked_dir[avhrr_dir] = 1
                write_log('INFO',"This directory does not exist, pattern:"
                          " {directory}".format(directory=avhrr_dir))
            continue

    
        file_pattern = insert_info_in_filename_or_path(filename_pattern,
                                                       values, datetime_obj=tobj)
        files = glob(os.path.join(found_dir, file_pattern))
        if len(files) > 0:
            no_files_found = False
            write_log('INFO',"Found files: " + os.path.basename(str(files[0])))
            return files[0], tobj
    if not found_dir:       
        return None, None
    if no_files_found:       
         write_log('INFO', "Found no files for patterns of type:"
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
                write_log('INFO', "Converting %r to HDF5" % f)
                if os.system(command):
                    raise RuntimeError("Couldn't convert %r to HDF5" % f)
                h5_files.append(h5_file)
        else:
            raise ValueError("File format of %r not recognized" % f)
    
    return h5_files

def get_pps_file(avhrr_file, options, values, type_of_file, file_dir):
    if not type_of_file in options:
        write_log("INFO", "No %s file in cfg-file!"%(type_of_file))
        return None
    date_time = values["date_time"]
    cloudtype_name = insert_info_in_filename_or_path(options[type_of_file], 
                                                     values, datetime_obj=date_time)                        
    path = insert_info_in_filename_or_path(options[file_dir], 
                                           values, datetime_obj=date_time)  
    try:
        file_name = glob(os.path.join(path, cloudtype_name))[0]
        #write_log("DEBUG", type_of_file +": ", file_name)
        return file_name
    except IndexError:
        write_log("INFO","No %s file found corresponding to %s." %( type_of_file, avhrr_file))
        return None

def find_files_from_avhrr(avhrr_file, options, as_oldstyle=False):
    """
    Find all files needed to process matchup from source data files.
    """
    # Let's get the satellite and the date-time of the pps radiance
    # (avhrr/viirs) file:
    write_log("INFO", "Avhrr or viirs file = %s" % avhrr_file)
    values = get_satid_datetime_orbit_from_fname(avhrr_file,
                                                 as_oldstyle=as_oldstyle)
    date_time = values["date_time"]

    cloudtype_name = insert_info_in_filename_or_path(options['cloudtype_file'], 
                                                     values, datetime_obj=date_time)
    path = insert_info_in_filename_or_path(options['cloudtype_dir'], 
                                           values, datetime_obj=date_time)
    try:
        print os.path.join(path, cloudtype_name)
        cloudtype_file = glob(os.path.join(path, cloudtype_name))[0]
    except IndexError:
        raise MatchupError("No cloudtype file found corresponding to %s." % avhrr_file)
    write_log("INFO", "CLOUDTYPE: " + cloudtype_file)
    ctth_name = insert_info_in_filename_or_path(options['ctth_file'], 
                                                values, datetime_obj=date_time)
    path =  insert_info_in_filename_or_path(options['ctth_dir'], 
                                                values, datetime_obj=date_time)  
    try:
        ctth_file = glob(os.path.join(path, ctth_name))[0]
    except IndexError:
        raise MatchupError("No ctth file found corresponding to %s." % avhrr_file)
    write_log("INFO", "CTTH: " + ctth_file)
    if VAL_CPP: 
        cpp_name = insert_info_in_filename_or_path(options['cpp_file'],
                                                   values, datetime_obj=date_time)
        path = insert_info_in_filename_or_path(options['cpp_dir'],
                                               values, datetime_obj=date_time)
        try:
            cpp_file = glob(os.path.join(path, cpp_name))[0]
        except IndexError:
            raise MatchupError("No cpp file found corresponding to %s." % avhrr_file)
        write_log("INFO", "CPP: " + cpp_file)
    else:     
       write_log("INFO", "Not validation of CPP ")
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
    write_log("INFO", "SUNSATANGLES: " + sunsatangles_file)

    if not 'physiography_file' in options:
        write_log("WARNING", "No physiography file searched for!")
        physiography_file = None
    else:
        physiography_name = insert_info_in_filename_or_path(options['physiography_file'],
                                                            values, datetime_obj=date_time)
        path =insert_info_in_filename_or_path(options['physiography_dir'], 
                                              values, datetime_obj=date_time)
        try:
            physiography_file = glob(os.path.join(path, physiography_name))[0]
        except IndexError:
            physiography_file = None
            raise MatchupError("No physiography file found corresponding to %s." % avhrr_file)
        write_log("INFO", "PHYSIOGRAPHY: " + physiography_file)

    if not 'nwp_tsur_file' in options:
        write_log("WARNING", "No surface temperature file searched for!")
        nwp_tsur_file=None
    else:
        nwp_tsur_name = insert_info_in_filename_or_path(options['nwp_tsur_file'],
                                                        values, datetime_obj=date_time)
        path = insert_info_in_filename_or_path(options['nwp_nwp_dir'],
                                               values, datetime_obj=date_time)
        try:
            nwp_tsur_file = glob(os.path.join(path, nwp_tsur_name))[0]
        except IndexError:
            raise MatchupError("No nwp_tsur file found corresponding to %s." % avhrr_file)
        write_log("INFO", "NWP_TSUR: " + nwp_tsur_file)

    file_name_dict={}
    for nwp_file in ['nwp_tsur','nwp_t500','nwp_t700',
                     'nwp_t850','nwp_t950', 'nwp_ciwv', 'nwp_ttro']:   
        file_name_dict[nwp_file] = get_pps_file(avhrr_file, options, values, 
                                                 nwp_file+'_file', 'nwp_nwp_dir')

    emis_file = get_pps_file(avhrr_file, options, values, 
                                           'emis_file', 'emis_dir')
    file_name_dict['emis'] = emis_file

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
                           'ctth': ctth_file,
                           'cpp': cpp_file,
                           'nwp_tsur': nwp_tsur_file,
                           'sunsatangles': sunsatangles_file,
                           'physiography': physiography_file})

    ppsfiles = ppsFiles(file_name_dict)
    return  ppsfiles


def _plot_avhrr_track(match_file, avhrr, track):
    """
    Helper function for calling `cloudsat_calipso_avhrr_plot.map_avhrr_track`.
    
    """
    from cloudsat_calipso_avhrr_plot import map_avhrr_track
    
    avhrr_lonlat = avhrr.longitude, avhrr.latitude
    track_lonlat = track.longitude, track.latitude
    fig = map_avhrr_track(avhrr_lonlat, track_lonlat)
    
    plot_file = match_file.replace('.h5', '.png')
    fig.savefig(plot_file)

def get_cloudsat_matchups(cloudsat_files, cloudtype_file, avhrrGeoObj, avhrrObj,
                          ctype, ctth, nwp_obj, avhrrAngObj, cppLwp, plot_file=None):
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

    if plot_file != None:
        write_log('INFO',"plot cloudsat - avhrr track")
        _plot_avhrr_track(plot_file, avhrrGeoObj, cloudsat)
    else:
        write_log('INFO',("No cloudsat plot-file"))
    cl_matchup, cl_min_diff, cl_max_diff = match_fun(cloudtype_file, cloudsat,
                                                     avhrrGeoObj, avhrrObj, ctype,
                                                     ctth, nwp_obj.surftemp, avhrrAngObj, cppLwp)

    return cl_matchup, (cl_min_diff, cl_max_diff)

def get_total_optical_depth_and_optical_depth_top_layer_from_5km_data(calipso, calipso5km=None, resolution=5):
    write_log('INFO',"Find total optical depth from 5km data")
    if calipso5km is None and resolution == 5:
        calipso5km = calipso
    calipso.feature_optical_depth_532_top_layer5km = -9.0 + 0*calipso.number_layers_found.ravel()
    calipso.total_optical_depth_5km = -9.0 + 0*calipso.number_layers_found.ravel()
    if resolution==5:
        pixels = np.logical_and(calipso5km.number_layers_found.ravel()>0,
                                calipso5km.feature_optical_depth_532[:,0].ravel()>=0)   
        calipso.feature_optical_depth_532_top_layer5km[pixels] = calipso5km.feature_optical_depth_532[pixels, 0]
        calipso.total_optical_depth_5km[pixels] =  calipso5km.feature_optical_depth_532[pixels, 0]       
        for lay in range(1,np.max(calipso5km.number_layers_found[pixels]),1):  
            pixels = np.logical_and(pixels,calipso5km.feature_optical_depth_532[:, lay]>=0)
            calipso.total_optical_depth_5km[pixels] +=  calipso5km.feature_optical_depth_532[pixels, lay]
    else:
        for pixel in range(calipso5km.utc_time.shape[0]): 
            if calipso5km.number_layers_found[pixel]>0 and calipso5km.feature_optical_depth_532[pixel, 0]>=0:
                calipso.feature_optical_depth_532_top_layer5km[pixel*5:pixel*5+5] = calipso5km.feature_optical_depth_532[pixel, 0] 
                calipso.total_optical_depth_5km[pixel*5:pixel*5+5] =  calipso5km.feature_optical_depth_532[pixel, 0]
                for lay in range(1,calipso5km.number_layers_found[pixel],1):  
                    if calipso5km.feature_optical_depth_532[pixel, lay]>=0:
                        calipso.total_optical_depth_5km[pixel*5:pixel*5+5] +=  calipso5km.feature_optical_depth_532[pixel, lay]  
    return calipso 

 

def get_calipso_matchups(calipso_files, values, 
                         avhrrGeoObj, avhrrObj, ctype, ctth, 
                         nwp_obj, avhrrAngObj, options, cppCph=None, 
                         nwp_segments=None, cafiles1km=None, cafiles5km=None, 
                         cafiles5km_aero=None):
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
    if cafiles1km != None:
        cafiles1km = discardCalipsoFilesOutsideTimeRange(
            cafiles1km, avhrrGeoObj, values, res=1)
    if cafiles5km != None:
        cafiles5km = discardCalipsoFilesOutsideTimeRange(
            cafiles5km, avhrrGeoObj, values, res=5)
    if cafiles5km_aero!=None:
        cafiles5km_aero = discardCalipsoFilesOutsideTimeRange(
            cafiles5km_aero, avhrrGeoObj, values, res=5, ALAY=True)
        
    calipso  = reshapeCalipso(calipso_files)
    #find time breakpoints, but don't cut the data yet ..
    startBreak, endBreak = find_break_points(calipso,  avhrrGeoObj, values)
    if cafiles1km != None:
        #RESOLUTION 5km also have 1km data
        calipso1km = reshapeCalipso(cafiles1km, res=1)
        calipso5km = calipso
        #data are also time reshaped in this function (add1km..)
        calipso = add1kmTo5km(calipso1km, calipso5km, startBreak, endBreak) 
        calipso = get_total_optical_depth_and_optical_depth_top_layer_from_5km_data(
            calipso, resolution=5)
    elif cafiles5km !=None:
        #RESOLUTION 1km also have 5km data
        calipso1km = calipso
        calipso5km = reshapeCalipso(cafiles5km, res=5)
        calipso1km = add5kmVariablesTo1kmresolution(calipso1km, calipso5km)
        calipso1km = get_total_optical_depth_and_optical_depth_top_layer_from_5km_data(
            calipso1km, calipso5km=calipso5km, resolution=1)
        if USE_5KM_FILES_TO_FILTER_CALIPSO_DATA:
            write_log('INFO',"Find detection height using 5km data")
            #data are also time reshaped in this function use5km ...
            calipso = use5km_find_detection_height_and_total_optical_thickness_faster(
                calipso1km, 
                calipso5km, 
                startBreak, 
                endBreak)
        else:
            calipso = time_reshape_calipso(calipso1km,  
                                           start_break=startBreak, 
                                           end_break=endBreak) 
    else:
        calipso = time_reshape_calipso(calipso,  
                                       start_break=startBreak, 
                                       end_break=endBreak)
        if RESOLUTION == 5:
            calipso = get_total_optical_depth_and_optical_depth_top_layer_from_5km_data(
                calipso, resolution=RESOLUTION)
    if cafiles5km_aero!=None:
        calipso5km_aero = reshapeCalipso(cafiles5km_aero, res=5, ALAY=True)
        if RESOLUTION == 1:
            calipso5km_aero = adjust5kmTo1kmresolution(calipso5km_aero)
        calipso5km_aero = time_reshape_calipso(calipso5km_aero,  
                                               start_break=startBreak, 
                                               end_break=endBreak)
        #for now insert aerosol flag into calipso_cloud object
        calipso.aerosol_flag = calipso5km_aero.feature_classification_flags[::]
    # free some memory    
    calipso1km = None
    calipso5km = None
        
    write_log('INFO',"Matching with avhrr")
    tup = match_calipso_avhrr(values, calipso,
                              avhrrGeoObj, avhrrObj, ctype,
                              ctth, cppCph, nwp_obj, avhrrAngObj, 
                              nwp_segments, options)
    ca_matchup, ca_min_diff, ca_max_diff = tup
    #import pdb; pdb.set_trace()
    return ca_matchup, (ca_min_diff, ca_max_diff)

def read_nwp(file_name, type_of_nwp):
    import epshdf #@UnresolvedImport
    if file_name is not None : # and config.CLOUDSAT_TYPE == "GEOPROF":
        nwpinst = epshdf.read_nwpdata(file_name)
        write_log("INFO", "Read NWP %s"%(type_of_nwp))
        return nwpinst.gain*nwpinst.data.astype('d') + nwpinst.intercept
    else:
        write_log("INFO","NO NWP %s File, Continue"%(type_of_nwp))
        return None

def read_segment_data(filename):
    import h5py
    product = {}
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        for attribute in list(h5file.attrs):
            product[attribute] = h5file.attrs[attribute]
        for attribute in list(h5file['satdef'].attrs):
            product[attribute] = h5file['satdef'].attrs[attribute]
        product['colidx'] = h5file['satdef']['colidx'].value
        product['rowidx'] = h5file['satdef']['rowidx'].value
        write_log("INFO", "Read segment info moisture")
        for moist_data in ['moist', 'surfaceMoist']:
            data = h5file['segment_area'][moist_data]
            gain = h5file.attrs['m_gain']
            intercept = h5file.attrs['m_intercept']
            nodata = h5file.attrs['m_nodata']
            #data[data!=nodata] = data[data!=nodata] * (gain) + intercept
            product[moist_data] = data
        write_log("INFO", "Read segment info pressure")
        for pressure_data in ['pressure', 'surfacePressure', 'ptro']:
            #pressure is in Pa in segments file
            data = h5file['segment_area'][pressure_data]
            gain = h5file.attrs['p_gain']
            intercept = h5file.attrs['p_intercept']
            nodata = h5file.attrs['p_nodata']
            data[data!=nodata] = data[data!=nodata] * (gain/100) + intercept/100 #Pa => hPa
            product[pressure_data] = data
        write_log("INFO", "Read segment info height")
        for geoheight_data in ['geoheight', 'surfaceGeoHeight']:
            #geo height is in meters in segment file
            data = h5file['segment_area'][geoheight_data]
            gain = h5file.attrs['h_gain']
            intercept = h5file.attrs['h_intercept']
            nodata = h5file.attrs['h_nodata']
            data[data!=nodata] = data[data!=nodata] * gain + intercept
            product[geoheight_data] = data
        write_log("INFO", "Read segment info temperature")
        for temperature_data in ['temp']:
            # temperature are measured in Kelvin in segmentfile
            data = h5file['segment_area'][temperature_data]
            gain = h5file.attrs['t_gain']
            intercept = h5file.attrs['t_intercept']
            nodata = h5file.attrs['t_nodata']
            data[data!=nodata] = data[data!=nodata] * gain + intercept
            product[temperature_data] = data
        for temperature_data in ['t850', 'ttro', 'surfaceLandTemp', 'surfaceSeaTemp']:
            data = h5file['segment_area'][temperature_data]
            gain = h5file.attrs['t_gain']
            intercept = h5file.attrs['t_intercept']
            nodata = h5file.attrs['t_nodata']
            data_float = np.array(data, dtype=np.float)
            data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
            product[temperature_data] = data_float
        for misc_data in ['meanElevation', 'fractionOfLand']:
            product[misc_data] = h5file['segment_area'][misc_data]
        write_log("INFO", "Read segment info brightness temperature")
        try:
            for tb_data in ['tb11clfree_sea',
                            'tb12clfree_sea',
                            'tb11clfree_land',
                            'tb12clfree_land']:
                data = h5file['segment_area'][tb_data]
                gain = h5file.attrs['t_gain']
                intercept = h5file.attrs['tb_intercept']
                nodata = h5file.attrs['t_nodata']
                data_float = np.array(data, dtype=np.float)
                data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
                product[tb_data] = data_float
        except ValueError:
            for tb_data, h5name in zip(['tb11clfree_sea',
                                        'tb12clfree_sea',
                                        'tb11clfree_land',
                                        'tb12clfree_land'],
                                       ['tb4clfree_sea',
                                        'tb5clfree_sea',
                                        'tb4clfree_land',
                                       'tb5clfree_land']):
                data = h5file['segment_area'][h5name]
                gain = h5file.attrs['t_gain']
                intercept = h5file.attrs['tb_intercept']
                nodata = h5file.attrs['t_nodata']
                data_float = np.array(data, dtype=np.float)
                data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
                product[tb_data] = data_float
        h5file.close()
        return product
    else:
        write_log("INFO","NO segment %s File, Continue"%(filename))
        return None


def read_thr(filename, h5_obj_type, thr_type):
    import h5py #@UnresolvedImport
    product = None
    if thr_type in ["emis1", "emis8", "emis9"]:
        if filename is not None: 
            h5file = h5py.File(filename, 'r')
            if 1==1:#h5_obj_type in h5file.keys():
                value = h5file[h5_obj_type].value
                gain = h5file.attrs['gain']
                intercept = h5file.attrs['intercept']
                product = value * gain + intercept
                write_log("INFO", "Read EMIS: %s"%(thr_type))
            else:
                write_log("ERROR","Could not read %s File, Continue"%(thr_type))
            h5file.close()   
        else:
            write_log("INFO","NO EMIS %s File, Continue"%(thr_type))
        return product  
    if filename is not None: 
        h5file = h5py.File(filename, 'r')
        if h5_obj_type in h5file.keys():
            value = h5file[h5_obj_type].value
            gain = h5file[h5_obj_type].attrs['gain']
            intercept = h5file[h5_obj_type].attrs['intercept']
            product = value * gain + intercept
            write_log("INFO", "Read THR: %s"%(thr_type))
        else:
            write_log("ERROR","Could not read %s File, Continue"%(thr_type))
        h5file.close()   
    else:
        write_log("INFO","NO THR %s File, Continue"%(thr_type))
    return product

def readViirsData_h5(filename):
    import h5py
    from pps_io import NewAvhrrData, AvhrrChannelData
    avhrrdata = NewAvhrrData()
    ich=0
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        nof_images = h5file['what'].attrs['sets']
        for num in xrange(1, nof_images+1,1):
            image = "image%d"%(num)
            channel = h5file[image].attrs['channel']
            if channel in ["M05",
                           "M07",
                           "M10",
                           "M12",
                           "M14",
                           "M15",
                           "M16",
                           "M11",
                           "M09"]:
 
                avhrrdata.channel.append(AvhrrChannelData())
                avhrrdata.channel[ich].data = h5file[image]['data'].value
                avhrrdata.channel[ich].des = "VIIRS %s"%(channel)
                avhrrdata.channel[ich].gain = h5file[image]['what'].attrs['gain']
                avhrrdata.channel[ich].intercept = h5file[image]['what'].attrs['offset']
                avhrrdata.nodata = h5file[image]['what'].attrs['nodata']
                avhrrdata.missingdata = h5file[image]['what'].attrs['missingdata']
                avhrrdata.no_data = h5file[image]['what'].attrs['nodata']
                avhrrdata.missing_data = h5file[image]['what'].attrs['missingdata']
                ich = ich + 1
    h5file.close()

    return avhrrdata

def readModisData_h5(filename):
    import h5py
    from pps_io import NewImagerData, ImagerChannelData
    avhrrdata = NewImagerData()
    ich=0
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        nof_images = h5file['what'].attrs['sets']
        for num in xrange(1, nof_images+1,1):
            image = "image%d"%(num)
            channel = h5file[image].attrs['channel']
            if channel in [
                    "1", "2", "3", "4", 
                    "5", "6", "7", "8", 
                    "9", "10", "11", "12", 
                    "13lo", "13hi", "14lo", "14hi", 
                    "15", "16", "17", "18",
                    "19", "20", "21", "22",
                    "23", "24", "25", "26",
                    "27", "28", "29", "30",
                    "31", "32", "33", "34",
                    "35", "36"]:
                avhrrdata.channel.append(ImagerChannelData())
                avhrrdata.channel[ich].data = h5file[image]['data'].value
                avhrrdata.channel[ich].des = "MODIS %s"%(channel)
                avhrrdata.channel[ich].gain = h5file[image]['what'].attrs['gain']
                avhrrdata.channel[ich].intercept = h5file[image]['what'].attrs['offset']
                avhrrdata.nodata = h5file[image]['what'].attrs['nodata']
                avhrrdata.missingdata = h5file[image]['what'].attrs['missingdata']
                avhrrdata.no_data = h5file[image]['what'].attrs['nodata']
                avhrrdata.missing_data = h5file[image]['what'].attrs['missingdata']
                ich = ich + 1
    h5file.close()

    return avhrrdata

def read_pps_data(pps_files, avhrr_file, cross):
    import pps_io #@UnresolvedImport
    import epshdf #@UnresolvedImport
    write_log("INFO","Read AVHRR geolocation data")
    avhrrGeoObj = pps_io.readAvhrrGeoData(avhrr_file)    
    write_log("INFO","Read AVHRR Sun -and Satellites Angles data")
    avhrrAngObj = pps_io.readSunSatAngles(pps_files.sunsatangles) #, withAbsoluteAzimuthAngles=True)
    if 'npp' in [cross.satellite1, cross.satellite2]:
        write_log("INFO","Read VIIRS data")
        avhrrObj = pps_io.readViirsData(avhrr_file)
        #avhrrObj = readViirsData_h5(avhrr_file)
    elif (cross.satellite1 in ['eos1','eos2'] or 
        cross.satellite2 in ['eos1', 'eos2']):
        #avhrrObj = pps_io.readModisData(avhrr_file)
        avhrrObj = readModisData_h5(avhrr_file)
    else:
        write_log("INFO","Read AVHRR data")
        avhrrObj = pps_io.readAvhrrData(avhrr_file)    
    cppLwp = None
    cppCph = None
    if VAL_CPP:    
        write_log("INFO","Read CPP data")
        from ppshdf_cloudproducts import CppProducts #@UnresolvedImport
        if PPS_FORMAT_2012_OR_EARLIER:
            try:
                cpp = CppProducts.from_h5(pps_files.cpp,
                                          product_names=['cph','lwp'])
                cppLwp = cpp.products['lwp'].array
                cppCph = cpp.products['cph'].array
                write_log("INFO", "CPP chp and lwp data read")
            except KeyError:
                #import traceback
                #traceback.print_exc()
                cppLwp = readCpp(pps_files.cpp, 'lwp')
                cppCph = readCpp(pps_files.cpp, 'cph')
                write_log("INFO", "CPP lwp, cph data read")
        else:
            cpp = CppProducts.from_h5(pps_files.cpp,
                                      product_names=['cpp_phase','cpp_lwp'],
                                      scale_up=True)
            # LWP convert from kg/m2 to g/m2
            cppLwp = 1000. * cpp.products['cpp_lwp'].array
            cppCph = cpp.products['cpp_phase'].array
            
    write_log("INFO","Read PPS Cloud Type")
    try:
        ctype = epshdf.read_cloudtype(pps_files.cloudtype, 1, 1, 1)  
    except:
        write_log("INFO","Could not use pps program to read, use mpop instead")    
        #read with mpop instead
        from mpop.satin.nwcsaf_pps import NwcSafPpsChannel
        ctype_mpop = NwcSafPpsChannel()
        print ctype_mpop
        ctype_mpop.read(pps_files.cloudtype)
        write_log("INFO","Done reading cloudtype") 
        print vars(ctype_mpop).keys()
        for varname in vars(ctype_mpop).keys():
            write_log("INFO",varname) 
        write_log("INFO","Done reading cloudtype") 
        #need to make empty ctype object here
        ctype_mpop.ct_quality.data
        #ctype=ctype_mpop
        ctype.cloudtype = ctype_mpop.cloudtype.data
        ctype.quality_flag = ctype_mpop.ct_quality.data
        ctype.conditions_flag = ctype_mpop.ct_conditions.data
        ctype.status_flag = ctype_mpop.ct_statusflag.data
        #print ctype_mpop.ct_quality

    write_log("INFO","Read PPS CTTH")
    try:
        ctth = epshdf.read_cloudtop(pps_files.ctth, 1, 1, 1, 0, 1)
    except:
        write_log("INFO","No CTTH")
        ctth = None

    write_log("INFO","Read PPS NWP data")   
    nwp_dict={}
    nwp_dict['surftemp'] = read_nwp(pps_files.nwp_tsur, "surface temperature")
    nwp_dict['t500'] = read_nwp(pps_files.nwp_t500, "temperature 500hPa")
    nwp_dict['t700'] = read_nwp(pps_files.nwp_t700, "temperature 700HPa")
    nwp_dict['t850'] = read_nwp(pps_files.nwp_t850, "temperature 850hPa")
    nwp_dict['t950'] = read_nwp(pps_files.nwp_t950, "temperature 950hPa")
    nwp_dict['ttro'] = read_nwp(pps_files.nwp_ttro, "tropopause temperature")
    nwp_dict['ciwv'] = read_nwp(pps_files.nwp_ciwv, 
                                "atmosphere_water_vapor_content")

    write_log("INFO","Read PPS threshold data")  
    for ttype in ['r06', 't11', 't37t12', 't37']:
        h5_obj_type = ttype +'_text'
        text_type = 'text_' + ttype
        nwp_dict[text_type] = read_thr(getattr(pps_files,text_type), 
                                       h5_obj_type,text_type)
    for h5_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv', 
                        't11ts', 't11t37', 't37t12', 't11t12',
                        'r09', 'r06', 't85t11_inv', 't85t11']:
        thr_type = 'thr_' + h5_obj_type
        nwp_dict[thr_type] = read_thr(getattr(pps_files,thr_type), 
                                      h5_obj_type, thr_type)
    write_log("INFO","Read PPS Emissivity data")  
    for h5_obj_type in ['emis1', 'emis8','emis9']:
        emis_type = h5_obj_type
        nwp_dict[emis_type] = read_thr(getattr(pps_files,"emis"), 
                                       h5_obj_type, emis_type)
    nwp_obj = NWPObj(nwp_dict)
    write_log("INFO","Read PPS NWP segment resolution data") 
    segment_data_object = read_segment_data(getattr(pps_files,'nwp_segments'))
    return avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, nwp_obj, cppLwp, cppCph, segment_data_object 

def read_cloud_cci(avhrr_file):
    from read_cloudproducts_cci import cci_read_all
    return cci_read_all(avhrr_file)

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
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, nwp_obj, cppLwp, cppCph, nwp_segment = retv
        date_time = values["date_time"]
    if (CCI_CLOUD_VALIDATION):
        avhrr_file, tobj = find_cci_cloud_file(cross, config_options)
        #avhrr_file = "20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
        values = get_satid_datetime_orbit_from_fname(avhrr_file)
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surftemp, cppLwp, cppCph =read_cloud_cci(avhrr_file)
        nwp_segment = None
        nwp_obj = NWPObj({'surftemp':surftemp})        
        avhrrGeoObj.satellite = values["satellite"];
        date_time = values["date_time"]

    calipso_files = find_calipso_files(date_time, config_options, values)

    if (PPS_VALIDATION):
        cloudsat_files = find_cloudsat_files(date_time, config_options, values)
        print cloudsat_files
        if (isinstance(cloudsat_files, str) == True or 
            (isinstance(cloudsat_files, list) and len(cloudsat_files) != 0)):
            write_log("INFO","Read CLOUDSAT %s data" % config.CLOUDSAT_TYPE)
            cl_matchup, cl_time_diff = get_cloudsat_matchups(cloudsat_files, 
                                                             pps_files.cloudtype,
                                                             avhrrGeoObj, avhrrObj, ctype,
                                                             ctth, nwp_obj, avhrrAngObj, cppLwp, config_options)
        else:
            write_log("INFO", "NO CLOUDSAT File, Continue")
    else:
        write_log("INFO", "NO CLOUDSAT File,"
                  "CCI-cloud validation only for calipso, Continue")

    if (isinstance(calipso_files, str) == True or 
        (isinstance(calipso_files, list) and len(calipso_files) != 0)):
        import glob
        calipso5km = None
        calipso1km = None
        
        if config.RESOLUTION == 5:
            if config.ALSO_USE_1KM_FILES == True:
                write_log("INFO", "Search for CALIPSO 1km data too")
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
            else:
                calipso1km = None

        if config.RESOLUTION == 1:
            if config.ALSO_USE_5KM_FILES == True:
                write_log("INFO", "Search for CALIPSO 5km data too")
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
            else:
                calipso5km = None

        calipso5km_aero=None
        if config.MATCH_AEROSOL_CALIPSO:
            calipso5km_aero=[]
            for cfile in calipso_files:
                file5km_aero = cfile.replace('/CLAY/', '/ALAY/').\
                               replace('CLay', 'ALay').\
                               replace('/1km/', '/5km/').\
                               replace('01km', '05km').\
                               replace('-ValStage1-V3-30.', '*').\
                               replace('-ValStage1-V3-01.', '*').\
                               replace('-ValStage1-V3-02.', '*').\
                               replace('-Prov-V3-01.', '*').\
                               replace('-Prov-V3-02.', '*').\
                               replace('-Prov-V3-30.', '*')
                files_found_aero = glob.glob(file5km_aero)
                if len(files_found_aero)==0:
                    #didn't find h5 file, might be hdf file instead
                    file5km_aero = file1km.replace('.h5','.hdf')
                    files_found_aero = glob.glob(file5km_aero)
                if len(files_found_aero)>0: 
                        calipso5km_aero.append(files_found_aero[0])                   
            calipso5km_aero = sorted(require_h5(calipso5km_aero))
            print "found these aerosol files", calipso5km_aero
                
                
        write_log("INFO", "Read CALIPSO data")        
        ca_matchup, ca_time_diff = get_calipso_matchups(calipso_files, 
                                                        values,
                                                        avhrrGeoObj, avhrrObj, ctype,
                                                        ctth, nwp_obj, avhrrAngObj, 
                                                        config_options, cppCph,
                                                        nwp_segment,
                                                        calipso1km, calipso5km, calipso5km_aero)

    else:
        write_log("INFO", "NO CALIPSO File, Continue")
    #import pdb; pdb.set_trace()

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
        write_log('INFO', "Creating dir %s:"%(rematched_path))
        os.makedirs(os.path.dirname(rematched_path))
    
    # Write cloudsat matchup
    if config.CLOUDSAT_TYPE == "GEOPROF":
        try:
            cl_match_file = rematched_file_base.replace(
                'atrain_datatype', 'cloudsat-%s' % config.CLOUDSAT_TYPE)
            writeCloudsatAvhrrMatchObj(cl_match_file, cl_matchup)
        except NameError:
            cl_matchup = None
            cl_time_diff = (NaN, NaN)
            print('CloudSat is not defined. No CloudSat Match File created')

    else:
        cl_match_file = rematched_file_base.replace(
            'atrain_datatype', 'cloudsat-%s' % config.CLOUDSAT_TYPE)
        writeCloudsatAvhrrMatchObj(cl_match_file, cl_matchup)

    # Write calipso matchup
    if config.CLOUDSAT_TYPE == "CWC-RVOD":
        ca_matchup = None
        ca_time_diff = (NaN, NaN)
    else:
        ca_match_file = rematched_file_base.replace('atrain_datatype', 'caliop')
        writeCaliopAvhrrMatchObj(ca_match_file, ca_matchup) 
    
    
    return {'cloudsat': cl_matchup, 'cloudsat_time_diff': cl_time_diff,
            'calipso': ca_matchup, 'calipso_time_diff': ca_time_diff,
            'basename': basename, 'values':values}


def get_matchups(cross, options, reprocess=False):
    """
    Retrieve Cloudsat- and Calipso-AVHRR matchups. If *reprocess* is False, and if
    matchup files exist, get matchups directly from the processed files.
    Otherwise process Cloudsat, Calipso, and PPS files first.
    """
    caObj = None
    clObj = None
    calipso_min_and_max_timediffs = NaN, NaN
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
        values["atrain_sat"] = "cloudsat-%s" % config.CLOUDSAT_TYPE[0]
        values["atrain_datatype"] = "cloudsat-%s" % config.CLOUDSAT_TYPE[0]
        cl_match_file, tobj = find_avhrr_file(cross, options['reshape_dir'], options['reshape_file'], values=values)
        if not cl_match_file:
            write_log('INFO', "No processed CloudSat match files found." + 
                      " Generating from source data.")
            clObj = None
            date_time=date_time_cross
        else:
            date_time=tobj
            td = tobj- date_time_cross
            matchup_diff_seconds=abs(td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6
            clObj = readCloudsatAvhrrMatchObj(cl_match_file) 
            basename = '_'.join(os.path.basename(cl_match_file).split('_')[1:5])
            if (diff_avhrr is None or 
                matchup_diff_seconds<=diff_avhrr_seconds or 
                abs(matchup_diff_seconds-diff_avhrr_seconds)<300 or
                abs(matchup_diff_seconds-diff_avhrr_seconds) <300):
                write_log('INFO', "CloudSat Matchups read from previously " + 
                          "processed data.")
                date_time=tobj
            else:
                write_log('INFO', "CloudSat Matchups will be processed for better match" + 
                          " %s."%values_avhrr["basename"]) 
                write_log('INFO', "CloudSat Matchups not read from previously " + 
                          "processed data %s."%basename)  
                clObj = None
                
        values["atrain_sat"] = "caliop"
        values["atrain_datatype"] = "caliop"
        ca_match_file, tobj = find_avhrr_file(cross, options['reshape_dir'], options['reshape_file'], values=values)
        if not  ca_match_file:
            write_log('INFO', 
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
                write_log('INFO', 
                          "CALIPSO Matchups read from previously processed data.")
                write_log('INFO', 'Filename: ' + ca_match_file)
                calipso_min_and_max_timediffs = (caObj.diff_sec_1970.min(), 
                                                 caObj.diff_sec_1970.max())
            else:
                write_log('INFO', "Calipso Matchups will be processed for better match" + 
                           " %s."%values_avhrr["basename"]) 
                write_log('INFO', "Calipso Matchups not read from previously " + 
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
            'basename': basename,
            'calipso_time_diff': calipso_min_and_max_timediffs,
            'values':values
            }

def run(cross, process_mode_dnt, config_options, min_optical_depth, reprocess=False):

    """
    The main work horse.
    
    """
    
    write_log('INFO', "Case: %s" % str(cross))
    write_log('INFO', "Process mode: %s" % process_mode_dnt)
    cloudsat_type = config.CLOUDSAT_TYPE
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
    write_log("INFO", "Sensor = " + sensor)

    #pdb.set_trace()


    # Now fetch all the datasets for the section of the AREA where all
    # three datasets match. Also get maximum and minimum time differences to AVHRR (in seconds)
    matchup_results = get_matchups(cross, config_options, reprocess)
    caObj = matchup_results['calipso']
    clsatObj = matchup_results['cloudsat']
    values = matchup_results['values']
    #import pdb;pdb.set_trace()

    clsat_min_diff, clsat_max_diff = matchup_results.get('cloudsat_time_diff', (NaN, NaN))
    ca_min_diff, ca_max_diff = matchup_results.get('calipso_time_diff', (NaN, NaN))
    
    basename = matchup_results['basename']
    base_sat = basename.split('_')[0]
    base_year = basename.split('_')[1][:4]
    base_month = basename.split('_')[1][4:6]
    
    if config.RESOLUTION not in [1, 5]:
        write_log("INFO","Define resolution")
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9) 

    ## Cloudsat ##
    if clsatObj is not None:
        if cloudsat_type == 'GEOPROF':
            cllon = clsatObj.cloudsat.longitude.copy()
            cllat = clsatObj.cloudsat.latitude.copy()
        elif cloudsat_type == 'CWC-RVOD':
            cllon = clsatObj.cloudsat.longitude.copy()
            cllat = clsatObj.cloudsat.latitude.copy()
            return

        # Issue a warning if startpoint or endpoint latitude of CloudSat and
        # CALIPSO differ by more than 0.1 degrees Furthermore, if startpoint
        # differs it means that CALIPSO data is not available for the first
        # part of the matchup cross section. This means that we must find the
        # first corresponding CloudSat point to this CALIPSO start point before 
        # we can do the plotting (statistics calculations are not
        # affected). Consequently, set the CALIPSO_DISPLACED flag and find
        # correct startpoint just before starting the plotting of CALIPSO data!

        CALIPSO_DISPLACED = 0
        latdiff = abs(clsatObj.cloudsat.latitude[0] - caObj.calipso.latitude[0])
        write_log('INFO',"latdiff: ", latdiff)
        if latdiff > 0.1:
            write_log('INFO', "CloudSat/CALIPSO startpoint differ by %f degrees." % latdiff)
            write_log('INFO', "Cloudsat start lon, lat: %f, %f" % (clsatObj.cloudsat.longitude[0], clsatObj.cloudsat.latitude[0]))
            write_log('INFO', "CALIPSO start lon, lat: %f, %f" % (caObj.calipso.longitude[0], caObj.calipso.latitude[0]))
            CALIPSO_DISPLACED = 1
            for j in range(clsatObj.cloudsat.latitude.shape[0]):
                if (abs(clsatObj.cloudsat.latitude[j] - caObj.calipso.latitude[0]) < 0.05) and\
                       (abs(clsatObj.cloudsat.longitude[j] - caObj.calipso.longitude[0]) < 0.1):
                    calipso_displacement=int(j*config.CLOUDSAT_TRACK_RESOLUTION)
                    write_log('INFO', "CALIPSO_DISPLACEMENT: %d" % calipso_displacement)
                    break
        # First make sure that PPS cloud top heights are converted to height
        # above sea level just as CloudSat and CALIPSO heights are defined. Use
        # corresponding DEM data.
        elevation = np.where(np.less_equal(clsatObj.cloudsat.elevation,0),\
                            -9,clsatObj.cloudsat.elevation)			# If clsatObj.cloudsat.elevation is <= 0 elevation(i,j)=-9, else the value = clsatObj.cloudsat.elevation(i,j)
        data_ok = np.ones(clsatObj.cloudsat.elevation.shape,'b')
        write_log('INFO', "Length of CLOUDSAT array: ", len(data_ok))
        avhrr_ctth_csat_ok = np.repeat(clsatObj.avhrr.ctth_height[::],data_ok)

        if CCI_CLOUD_VALIDATION: #ctth relative mean sea level
            avhrr_ctth_csat_ok = np.where(np.greater(avhrr_ctth_csat_ok,0.0),avhrr_ctth_csat_ok[::],avhrr_ctth_csat_ok)
        else: #ctth relative topography
            avhrr_ctth_csat_ok = np.where(np.greater(avhrr_ctth_csat_ok,0.0),avhrr_ctth_csat_ok[::]+elevation*1.0,avhrr_ctth_csat_ok)

        if len(data_ok) == 0:
            write_log('INFO',"Processing stopped: Zero lenght of matching arrays!")
            print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit()
    else:
        data_ok = None
        avhrr_ctth_csat_ok = None

    ## Calipso ##        
    #import pdb;pdb.set_trace()
    calon = caObj.calipso.longitude.copy()
    calat = caObj.calipso.latitude.copy()

    # First make sure that PPS cloud top heights are converted to height above sea level
    # just as CloudSat and CALIPSO heights are defined. Use corresponding DEM data.
    cal_elevation = np.where(np.less_equal(caObj.calipso.elevation,0),
                                -9,caObj.calipso.elevation)
    cal_data_ok = np.ones(caObj.calipso.elevation.shape,'b')
    write_log('INFO', "Length of CALIOP array: ", len(cal_data_ok))

    avhrr_ctth_cal_ok = np.repeat(caObj.avhrr.ctth_height[::],cal_data_ok)

    if CCI_CLOUD_VALIDATION: #ctth relative mean sea level
        avhrr_ctth_cal_ok = np.where(np.greater(avhrr_ctth_cal_ok,0.0),avhrr_ctth_cal_ok[::],avhrr_ctth_cal_ok)
    else: #ctth relative topography
        avhrr_ctth_cal_ok = np.where(np.greater(avhrr_ctth_cal_ok,0.0),avhrr_ctth_cal_ok[::]+cal_elevation,avhrr_ctth_cal_ok)
        
    if (len(cal_data_ok) == 0):
        write_log('INFO', "Processing stopped: Zero lenght of matching arrays!")
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit()
        
    # If everything is OK, now create filename for statistics output file and
    # open it for writing.  Notice that more than one file (but maximum 2) can
    # be created for one particular noaa orbit.
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
    statfile = open(statfilename,"w")
    if process_mode == "BASIC":
        if clsatObj is not None:
            statfile.write("CloudSat min and max time diff: %f %f \n" %(clsat_min_diff,clsat_max_diff))
        else:
            statfile.write('No CloudSat \n')
        statfile.write("CALIPSO min and max time diff: %f %f \n" %(ca_min_diff,ca_max_diff))
    else:
        if clsatObj is not None:
            statfile.write("CloudSat min and max time diff: See results for BASIC! \n")
        else:
            statfile.write('No CloudSat \n')
        statfile.write("CALIPSO min and max time diff: See results for BASIC! \n")
    if clsatObj is not None:
        statfile.write("Start-Stop-Length Cloudsat: %f %f %f %f %s \n" %(clsatObj.cloudsat.latitude[0],clsatObj.cloudsat.longitude[0],clsatObj.cloudsat.latitude[len(clsatObj.cloudsat.latitude)-1],clsatObj.cloudsat.longitude[len(clsatObj.cloudsat.latitude)-1],len(data_ok)))
    else:
        statfile.write('No CloudSat \n')
    statfile.write("Start-Stop-Length CALIPSO: %f %f %f %f %s \n" %(caObj.calipso.latitude[0],caObj.calipso.longitude[0],caObj.calipso.latitude[len(caObj.calipso.latitude)-1],caObj.calipso.longitude[len(caObj.calipso.latitude)-1],len(cal_data_ok)))


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
        print "WARNING:", "rematched_file is missing total_optical_depth_5km field"
        print "INFO:", "consider reprocessing with "
        print "COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC=True"
        print "ALSO_USE_5KM_FILES=True or RESOLUTION==5"

    check_total_optical_depth_and_warn(caObj)

    if process_mode != 'BASIC' and RESOLUTION==1:
        caObj = CalipsoOpticalDepthHeightFiltering1km(caObj)

    if process_mode == 'OPTICAL_DEPTH_THIN_IS_CLEAR' and RESOLUTION==1:
        write_log('INFO',"Setting thin clouds to clear"
                  ", using 5km data in mode OPTICAL_DEPTH_THIN_IS_CLEAR")
        caObj = CalipsoOpticalDepthSetThinToClearFiltering1km(caObj)      

    # Extract CALIOP Vertical Feature Classification (bits 10-12) from 16 bit
    # representation for topmost cloud layer #Nina 20140120 this is cloud type nit vert feature !!
    cal_vert_feature = np.ones(caObj.calipso.layer_top_altitude[::,0].shape)*-9
    feature_array = 4*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],11),1) + 2*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],10),1) + np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],9),1)
    cal_vert_feature = np.where(np.not_equal(caObj.calipso.feature_classification_flags[::,0],1),feature_array[::],cal_vert_feature[::])  
    #feature_array = 2*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],4),1) + np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],3),1)

    # Prepare for plotting, cloud emissivity and statistics calculations
    
    caliop_toplay_thickness = np.ones(caObj.calipso.layer_top_altitude[::,0].shape)*-9

    thickness = (caObj.calipso.layer_top_altitude[::,0]-caObj.calipso.layer_base_altitude[::,0])*1000.
    caliop_toplay_thickness = np.where(np.greater(caObj.calipso.layer_top_altitude[::,0],-9),thickness,caliop_toplay_thickness)
    
    caliop_height = []
    caliop_base = []
    caliop_max_height = np.ones(caObj.calipso.layer_top_altitude[::,0].shape)*-9
    
    for i in range(10):
        hh = np.where(np.greater(caObj.calipso.layer_top_altitude[::,i],-9),
                            caObj.calipso.layer_top_altitude[::,i] * 1000.,-9)
                                        
        caliop_max_height = np.maximum(caliop_max_height,
                                       caObj.calipso.layer_top_altitude[::,i] * 1000.)
        # This is actually unnecessary - we know that layer 1 is always the
        # highest layer!!  However, arrays caliop_height and caliop_base are
        # needed later for plotting/ KG

        bb = np.where(np.greater(caObj.calipso.layer_base_altitude[::,i],-9),
                            caObj.calipso.layer_base_altitude[::,i] * 1000.,-9)
        #if (hh>bb):
        caliop_height.append(hh)
        caliop_base.append(bb)
        thickness = hh - bb
        
    x = np.repeat(caObj.calipso.number_layers_found.ravel(),
                        np.greater(caObj.calipso.number_layers_found.ravel(),0))
    #print "Number of points with more than 0 layers: ",x.shape[0]
    cal_data_ok = np.greater(caliop_max_height,-9.)

    if np.size(caObj.avhrr.surftemp)>1 or caObj.avhrr.surftemp != None:
        print caObj.avhrr.surftemp.shape, cal_data_ok.shape
        cal_surftemp_ok = np.repeat(caObj.avhrr.surftemp[::], cal_data_ok)

    if clsatObj is not None:
        # Transfer CloudSat MODIS cloud flag to CALIPSO representation
        cal_MODIS_cflag = np.zeros(len(cal_data_ok),'b')
        for i in range(len(cal_data_ok)):
            if not CALIPSO_DISPLACED: 
                cloudsat_index = int(i/config.CLOUDSAT_TRACK_RESOLUTION + 0.5)
                #print len(clsatObj.cloudsat.MODIS_cloud_flag),cloudsat_index
                if len(clsatObj.cloudsat.MODIS_cloud_flag) > (cloudsat_index):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
            else:
                cloudsat_index = int((i+calipso_displacement)/config.CLOUDSAT_TRACK_RESOLUTION + 0.5)
                if len(clsatObj.cloudsat.MODIS_cloud_flag) > (cloudsat_index-1):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
    else:
        cal_MODIS_cflag = None
                
    ##########################
    ### 1 KM DATA CWC-RVOD ###                       
    if config.RESOLUTION == 1 and cloudsat_type == 'CWC-RVOD' and clsatObj is not None:
        elevationcwc = np.where(np.less_equal(clsatObj.cloudsatcwc.elevation,0),
                            -9, clsatObj.cloudsatcwc.elevation)
        data_okcwc = np.ones(clsatObj.cloudsatcwc.elevation.shape,'b')
                
    ### 5 KM DATA CWC-RVOD ###                       
    elif config.RESOLUTION == 5 and cloudsat_type == 'CWC-RVOD' and clsatObj is not None:
        elevationcwc = np.where(np.less_equal(clsatObj.cloudsat5kmcwc.elevation,0),
                            -9, clsatObj.cloudsat5kmcwc.elevation,-9)
        data_okcwc = np.ones(clsatObj.cloudsat5kmcwc.elevation.shape,'b')
   
    #======================================================================
    # Draw plot
    if process_mode_dnt in config.PLOT_MODES:
        write_log('INFO', "Plotting")

        plotpath = insert_info_in_filename_or_path(config_options['plot_dir'], values,
                                                   datetime_obj=values['date_time'])      
        trajectorypath = os.path.join(plotpath, "trajectory_plot")
        if not os.path.exists(trajectorypath):
            os.makedirs(trajectorypath)
        trajectoryname = os.path.join(trajectorypath, 
                                      "%skm_%s_trajectory" % (int(config.RESOLUTION),
                                                              values['basename']))
        # To make it possible to use the same function call to drawCalClsatGEOPROFAvhrr*kmPlot
        # in any processing mode:
        #This is always None!
        emissfilt_calipso_ok = None
        if clsatObj is None:
            file_type = ['eps', 'png']
            if PLOT_ONLY_PNG==True:
                file_type = ['png']
            if 'trajectory_plot_area' in config_options:
                plotSatelliteTrajectory(calon, calat,
                                        trajectoryname, 
                                        config.AREA_CONFIG_FILE, 
                                        file_type,
                                        area_id=config_options['trajectory_plot_area'])
            else:
                plotSatelliteTrajectory(calon, calat,
                                        trajectoryname, 
                                        config.AREA_CONFIG_FILE,
                                        file_type)

            drawCalClsatAvhrrPlotTimeDiff(calat, None, 
                                          caObj.diff_sec_1970, 
                                          plotpath, basename, 
                                          config.RESOLUTION, file_type,
                                          instrument=sensor)
            drawCalClsatAvhrrPlotSATZ(calat, None, 
                                      caObj.avhrr.satz, 
                                      plotpath, basename, 
                                      config.RESOLUTION, file_type,
                                      instrument=sensor)
            drawCalClsatGEOPROFAvhrrPlot(None, 
                                         caObj.calipso, None, data_ok,
                                         None, caliop_base,
                                         caliop_height, cal_data_ok,
                                         avhrr_ctth_cal_ok, plotpath,
                                         basename, process_mode, 
                                         emissfilt_calipso_ok, file_type,
                                         instrument=sensor)
        else:                    
            if cloudsat_type=='GEOPROF':
                file_type = ['eps', 'png']
                drawCalClsatAvhrrPlotSATZ(cllat, 
                                          clsatObj.avhrr.satz, 
                                          caObj.avhrr.satz, 
                                          plotpath, basename, 
                                          config.RESOLUTION, file_type,
                                          instrument=sensor)
                drawCalClsatGEOPROFAvhrrPlot(clsatObj.cloudsat, 
                                             clsatObj.avhrr, caObj.calipso, 
                                             caObj.avhrr, elevation, data_ok,
                                             CALIPSO_DISPLACED, caliop_base,
                                             caliop_height, cal_data_ok,
                                             avhrr_ctth_cal_ok, plotpath,
                                             basename, process_mode, 
                                             emissfilt_calipso_ok, file_type,
                                             instrument=sensor)
                drawCalClsatAvhrrPlotTimeDiff(cllat, 
                                              clsatObj.diff_sec_1970, 
                                              caObj.diff_sec_1970, 
                                              plotpath, basename, 
                                              config.RESOLUTION, file_type,
                                              instrument=sensor)
                if 'trajectory_plot_area' in config_options:
                    plotSatelliteTrajectory(cllon, cllat, trajectoryname,
                                            config.AREA_CONFIG_FILE,
                                            file_type,
                                            area_id=config_options['trajectory_plot_area'])
                else:
                    plotSatelliteTrajectory(cllon, cllat, 
                                            trajectoryname,
                                            config.AREA_CONFIG_FILE,
                                            file_type)
                
            elif cloudsat_type=='CWC-RVOD':
                drawCalClsatAvhrrPlotTimeDiff(cllat, 
                                              clsatObj.diff_sec_1970, 
                                              caObj.diff_sec_1970, 
                                              plotpath, basename, 
                                              config.RESOLUTION,
                                              instrument=sensor)
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
                
    #==============================================================
    #Calculate Statistics
    if cloudsat_type == 'GEOPROF':
        process_calipso_ok = 0
        
    write_log('INFO', "Calculating statistics")
    CalculateStatistics(process_mode, clsatObj, statfile, caObj, cal_MODIS_cflag,
                        cal_vert_feature, avhrr_ctth_csat_ok, data_ok,
                        cal_data_ok, avhrr_ctth_cal_ok, caliop_max_height,
                        process_calipso_ok, dnt_flag)
