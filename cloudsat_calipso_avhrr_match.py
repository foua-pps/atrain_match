"""
TODO: The following description needs updating. atrain_match is now run via
process_master.py.

Program cloudsat_calipso_avhrr_match.py

This program is used to process and output statistics for the inter-comparison
of AVHRR PPS results and CloudSat/CALIPSO observations. It may be run
repeatedly and supervised by program cloudsat_calipso_process_master.py.

This particular version of Adam's original CloudSat/CALIPSO matchup and
analysis program has been complemented with the following:

 * A method to calculate cloud emissivities for the uppermost CALIPSO cloud
   layer. With the use of parameters EMISS_FILTERING, EMISS_MIN_HEIGHT and
   EMISS_LIMIT the thinnest uppermost CALIPSO cloud layers can be analysed and
   the entire column can be disregarded if the emissivity falls below the
   EMISS_LIMIT value.  Cloud emissivities Ec are calculated as follows:

                    Ec = (I-Iclear)/(B(Tc)-Iclear)
   where

       I = Measured radiance in AVHRR channel 4 (11 micron) To be calculated as
           the Planck radiance for the associated brightness temperature

       Iclear = Estimated radiance in cloud free situations To be calculate as
                the Planck radiance for the NWP-analysed surface temperature
                (i.e., neglecting further atmospheric contributions)

       B(Tc) = Planck radiance for the uppermost cloud layer using CALIPSO
       mid-layer temperatures


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

The program is capable of running in a wide range of modes (according to
description above). These various modes are selected by enabling (disabling)
the following parameters:

PLOT_OPTION, EMISS_FILTERING, ICE_COVER_SEA, ICE_FREE_SEA, SNOW_COVER_LAND,
SNOW_FREE_LAND, COASTAL_ZONE

However, notice that only one mode can be chosen for each run. The only
exception is PLOT_OPTION (i.e., the generation of a PNG plot) which can be
combined with EMISS_FILTERING. This also means that processing of individual
surface categories only generates statistics and not any plots.

Input data has to be supplied at directories defined by MAIN_DIR and SUB_DIR
parameters below.
 
Every exection of the program prints out statistics for PPS Cloud Mask, Cloud
Type and Cloud Top Height directly on the screen.

Don't forget to set all parameters needed for standard ACPG/AHAMAP execution
since the matchup software uses parts of the ACPG/AHAMAP software.

Dependencies: For a successful run of the program the following supporting
              python modules must be available in the default run directory:

              cloudsat.py
              calipso.py
              calipso_avhrr_matchup.py
              cloudsat_avhrr_matchup.py
              radiance_tb_tables_kgtest.py

              For full consistency make sure that MAIN_DIR and SUB_DIR
              parameters are the same also in modules cloudsat.py and
              calipso.py.

Output data files: Main results are generally written in the directory
                   MAIN_DIR/SUB_DIR but plotting results are stored at ./Plot
                   and temporary results at ./Data directories. Thus, make sure
                   that these directories exist as subdirectories at the
                   default run directory.  This is now made automatic /Erik

Finally, notice that the matching of the PPS, CloudSat and CALIPSO datasets
have been calculated using the fix area arctic_super_5010 defined over the
Arctic region in Lambert Azimuthal Equal Area projection. Thus, for matching
data to other regions please modify modules calipso.py and cloudsat.py and
replace area arctic_super_5010 with the desired area.  This is now made below
/Erik

/KG March 2010

"""
import os
import sys

import numpy as np
from numpy import NaN
from config import (VAL_CPP,
                    PLOT_ONLY_PNG,
                    CCI_CLOUD_VALIDATION,
                    PPS_VALIDATION,
                    ALWAYS_USE_AVHRR_ORBIT_THAT_STARTS_BEFORE_CROSS)

from radiance_tb_tables_kgtest import get_central_wavenumber
# Just use the brightness temperature to
# radiance conversion/KG


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
from calipso import reshapeCalipso, match_calipso_avhrr
from calipso import (writeCaliopAvhrrMatchObj, 
                     readCaliopAvhrrMatchObj,
                     use5km_remove_thin_clouds_from_1km,
                     add1kmTo5km)
import inspect
import numpy as np
from cloudsat_calipso_avhrr_plot import (drawCalClsatAvhrrPlotTimeDiff,
                                         drawCalClsatGEOPROFAvhrrPlot, 
                                         drawCalClsatAvhrrPlotSATZ,
                                         drawCalClsatCWCAvhrrPlot)
from config import CALIPSO_CLOUD_FRACTION

INSTRUMENT = {'npp': 'viirs',
              'noaa18': 'avhrr'}
from datetime import datetime, timedelta
from glob import glob

class ppsFiles(object):
    def __init__(self, **kwargs):
        if 'cloudtype' in kwargs:
            self.cloudtype = kwargs['cloudtype']
        else:
            self.cloudtype = None
        if 'ctth' in kwargs:
            self.ctth = kwargs['ctth']
        else:
            self.ctth = None
        if 'cpp' in kwargs:
            self.cpp = kwargs['cpp']
        else:
            self.cpp = None
        if 'sunsatangles' in kwargs:
            self.sunsatangles = kwargs['sunsatangles']
        else:
            self.sunsatangles = None
        if 'nwp_tsur' in kwargs:
            self.nwp_tsur = kwargs['nwp_tsur']
        else:
            self.nwp_tsur = None
    

test = 0

def readCpp(filename, type):
    import h5py #@UnresolvedImport
    h5file = h5py.File(filename, 'r')
    if type in h5file.keys():
        value = h5file[type].value
        gain = h5file[type].attrs['gain']
        intersec = h5file[type].attrs['intercept']
        nodat = h5file[type].attrs['no_data_value']
        product = np.where(value != nodat,value * gain + intersec, value)   
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
        flist.extend([ s for s in tmplist if s not in flist ])      
    return flist

def get_satid_datetime_orbit_from_fname(avhrr_filename):
    import runutils
    #satname, _datetime, orbit = runutils.parse_scene(avhrr_filename)
    #returnd orbit as int, loosing leeding zeros, use %05d to get it right.
    # Get satellite name, time, and orbit number from avhrr_file
    sl_ = os.path.basename(avhrr_filename).split('_')
    date_time= datetime.strptime(sl_[1] + sl_[2], '%Y%m%d%H%M')
    values= {"satellite": sl_[0],
             "date_time": date_time,
             "orbit": sl_[3],
             "date":sl_[1],
             "year":date_time.year,
             "month":"%02d"%(date_time.month),    
             "time":sl_[2],
             "basename":sl_[1] + sl_[2]+ sl_[3],
             "ppsfilename":avhrr_filename}
    return values

def insert_info_in_filename_or_path(file_or_name_path, values, datetime_obj=None):
    #file_or_name_path = file_or_name_path.format(**values)
    file_or_name_path = file_or_name_path.format(
        satellite=values["satellite"],
        orbit=values.get("orbit","*"),
        instrument = INSTRUMENT.get(values["satellite"],"avhrr"),
        resolution=config.RESOLUTION,
        area=config.AREA,
        val_dir=config._validation_results_dir,
        year=values.get('year',"unknown"),
        month=values.get('month',"unknown"),
        mode=values.get('mode',"unknown"),
        min_opt_depth=('min_opt_depth',"unknown"),
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
    time_low = cross_time-ddt
    time_high = cross_time -ddt2
    write_log("INFO", "Searching for avhrr/viirs file with start time  between" 
              ": %s and %s  "%(time_low.strftime("%d %H:%M"),time_high.strftime("%d %H:%M"))) 
    time_window = (ddt, ddt2)
    # Make list of file times to search from:
    tlist = get_time_list(cross_time, time_window, delta_t_in_seconds=60)
    return tlist, cross_time, cross_satellite



def find_avhrr_file_old(cross, options):
    """
    Find the *satellite* avhrr file closest to *datetime*.
    """
    (tlist, cross_time, cross_satellite) = get_time_list_and_cross_time(cross)
    found_dir = None
    for tobj in tlist:
        rad_dir = insert_info_in_filename_or_path(options['radiance_dir'],
                                                  values, datetime_obj=tobj)
        if os.path.exists(rad_dir):
            found_dir = rad_dir
            break
    if not found_dir:
        raise MatchupError("No dir found with radiance data!\n" + 
                           "Searching under %s" % options['radiance_dir'])
    for tobj in tlist:
        file_pattern = insert_info_in_filename_or_path(options['radiance_file'],
                                                       values,datetime_obj=tobj)
        radiance_files = glob(os.path.join(found_dir, file_pattern))
        if len(radiance_files) > 0:
            write_log('INFO',"Found files: " + str(radiance_files))
            return radiance_files[0]
    return None

def find_radiance_file(cross, options):
    found_file, tobj= find_avhrr_file(cross, 
                                      options['radiance_dir'], 
                                      options['radiance_file'])
    if not found_file:
        raise MatchupError("No dir or file found with radiance data!\n" + 
                           "Searching under %s" % options['radiance_dir'])
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
    for tobj in tlist:
        #print values
        avhrr_dir = insert_info_in_filename_or_path(filedir_pattern,
                                                  values, datetime_obj=tobj)
        if os.path.exists(avhrr_dir):
            found_dir = avhrr_dir
            write_log('INFO',"Found directory "
                      "{dirctory} ".format(dirctory=found_dir))
            break
    if not found_dir:
        write_log('INFO',"This directory does not exist, pattern:"
                  " {directory}".format(directory=filedir_pattern))
        return None, None
    no_files_found = True
    for tobj in tlist:
        file_pattern = insert_info_in_filename_or_path(filename_pattern,
                                                       values, datetime_obj=tobj)  
        files = glob(os.path.join(found_dir, file_pattern))
        if len(files) > 0:
            no_files_found = False
            write_log('INFO',"Found files: " + os.path.basename(str(files[0])))
            return files[0], tobj
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

def find_cloudsat_files(date_time, options):
    """
    Find calipso and cloudsat files closest to date_time.
    """

    #TOTDO%%%%%%%%%%%%
    import file_finders
    cloudsat_finder = file_finders.CloudsatFileFinder(config.CLOUDSAT_DIR,
                                                      config.RESOLUTION,
                                                      config.CLOUDSAT_TYPE)
    attach_subdir_from_config(cloudsat_finder)
    cloudsat_finder.set_time_window(config.SAT_ORBIT_DURATION + config.sec_timeThr)

    cloudsat_files = sorted(require_h5(cloudsat_finder.find(date_time)))
    return cloudsat_files
    #%%%%%%%%%%%%5

def find_files_from_avhrr(avhrr_file, options):
    """
    Find all files needed to process matchup from source data files.
    """
    # Let's get the satellite and the date-time of the pps radiance
    # (avhrr/viirs) file:
    write_log("INFO", "Avhrr or viirs file = %s" % avhrr_file)
    values = get_satid_datetime_orbit_from_fname(avhrr_file)
    date_time = values["date_time"]
    cloudtype_name = insert_info_in_filename_or_path(options['cloudtype_file'], 
                                                     values, datetime_obj=date_time)                        
    path = insert_info_in_filename_or_path(options['cloudtype_dir'], 
                                           values, datetime_obj=date_time)  
    try:
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
    if not 'nwp_tsur_file' in options:
        write_log("WARNING", "No t-surface file searched for!")
        nwp_tsur_file=None
    else:
        nwp_tsur_name = insert_info_in_filename_or_path(options['nwp_tsur_file'],
                                                        values, datetime_obj=date_time)
        path = insert_info_in_filename_or_path(options['nwp_tsur_dir'],
                                               values, datetime_obj=date_time)
        try:
            nwp_tsur_file = glob(os.path.join(path, nwp_tsur_name))[0]
        except IndexError:
            raise MatchupError("No nwp_tsur file found corresponding to %s." % avhrr_file)
        write_log("INFO", "NWP_TSUR: " + nwp_tsur_file)
        
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

    ppsfiles = ppsFiles(cloudtype=cloudtype_file,
                        ctth=ctth_file,
                        cpp=cpp_file,
                        nwp_tsur=nwp_tsur_file,
                        sunsatangles=sunsatangles_file)

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
                          ctype, ctth, surft, avhrrAngObj, cppLwp, plot_file=None):
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
                                                     ctth, surft, avhrrAngObj, cppLwp)

    return cl_matchup, (cl_min_diff, cl_max_diff)


def get_calipso_matchups(calipso_files, values, 
                         avhrrGeoObj, avhrrObj, ctype, ctth, 
                         surft, avhrrAngObj, options, cppCph=None, cafiles1km=None, cafiles5km=None):
    """
    Read Calipso data and match with the given PPS data.
    """
    if cafiles1km != None:
        #pdb.set_trace()
        calipso1km = reshapeCalipso(cafiles1km, avhrrGeoObj, values, False, 1)[0]
        calipso5km, startBreak, endBreak = reshapeCalipso(calipso_files,avhrrGeoObj, values, False)
        calipso = add1kmTo5km(calipso1km, calipso5km, startBreak, endBreak)
        calipso1km = None
        calipso5km = None
    elif cafiles5km !=None:
        calipso1km, startBreak, endBreak  = reshapeCalipso(calipso_files, avhrrGeoObj, values, False, 1)
        calipso5km = reshapeCalipso(cafiles5km,avhrrGeoObj, values, False, 5)[0]
        write_log('INFO',"Cut optically thin clouds at selected optical depth, using 5km data")
        calipso = use5km_remove_thin_clouds_from_1km(calipso1km, calipso5km, startBreak, endBreak)
        calipso1km = None
        calipso5km = None
    else:
        calipso = reshapeCalipso(calipso_files, 
                                 avhrrGeoObj, values, True)[0]
    write_log('INFO',"Matching with avhrr")
    tup = match_calipso_avhrr(values, calipso,
                              avhrrGeoObj, avhrrObj, ctype,
                              ctth, cppCph, surft, avhrrAngObj, options)
    ca_matchup, ca_min_diff, ca_max_diff = tup
    #import pdb; pdb.set_trace()
    return ca_matchup, (ca_min_diff, ca_max_diff)

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
    else:
        write_log("INFO","Read AVHRR data")
        avhrrObj = pps_io.readAvhrrData(avhrr_file)    
    cppLwp = None
    cppCph = None
    if VAL_CPP:    
        write_log("INFO","Read CPP data")
        try:
            from ppshdf_cloudproducts import CppProducts #@UnresolvedImport
            cpp = CppProducts.from_h5(pps_files.cpp, product_names=['cph','lwp'])        
            cppLwp = cpp.products['lwp'].array
            cppCph = cpp.products['cph'].array
            write_log("INFO", "CPP chp and lwp data read")
        except KeyError:
        #import traceback
        #traceback.print_exc()
            cppLwp = readCpp(pps_files.cpp, 'lwp')
            cppCph = readCpp(pps_files.cpp, 'cph')
            write_log("INFO", "CPP lwp, cph data read")
    write_log("INFO","Read PPS Cloud Type")
    ctype = epshdf.read_cloudtype(pps_files.cloudtype, 1, 1, 0)
    write_log("INFO","Read PPS CTTH")
    try:
        ctth = epshdf.read_cloudtop(pps_files.ctth, 1, 1, 1, 0, 1)
    except:
        ctth = None  
    surft = None
    if pps_files.nwp_tsur is not None : # and config.CLOUDSAT_TYPE == "GEOPROF":
        if 1:
            nwpinst = epshdf.read_nwpdata(pps_files.nwp_tsur)
            write_log("INFO", "Read NWP surface temperature")
            surft = nwpinst.gain*nwpinst.data.astype('d') + nwpinst.intercept
        #except:
        #    write_log("INFO", 
        #              "Corrupted NWP surface temperature File, Continue")
        #    surft = None
        else:
            write_log("INFO","NO NWP surface temperature File, Continue")
            surft = None
    return avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surft, cppLwp, cppCph 

def read_cloud_cci(avhrr_file):
    from MakeCloudTopPPSObjFromCLoudCSINetcdf4file import cci_read_ctth
    return cci_read_ctth(avhrr_file)

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
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surft, cppLwp, cppCph =read_pps_data(pps_files, avhrr_file, cross)
        date_time = values["date_time"]
    if (CCI_CLOUD_VALIDATION):
        avhrr_file, tobj = find_cci_cloud_file(cross, config_options)
        #avhrr_file = "20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
        avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surft, cppLwp, cppCph =read_cloud_cci(avhrr_file)
        date_time=datetime.strptime("200806130022", '%Y%m%d%H%M')
        values= {"satellite": "noaa18",
                 "date_time": date_time,
                 "orbit": "9999",
                 "date":"20080613",
                 "year":date_time.year,
                 "month":"%02d"%(date_time.month),    
                 "time":"0022",
                 "basename":"20080613002200-ESACCI",
                 "ccifilename":avhrr_file,
                 "ppsfilename":"noaa18_20080613_0022_99999_satproj_00000_13793_cloudtype.h5"}

        avhrrGeoObj.satellite = "noaa18";
    calipso_files = find_calipso_files(date_time, config_options, values)

    if (PPS_VALIDATION):
        cloudsat_files = find_cloudsat_files(date_time, config_options)
        if (isinstance(cloudsat_files, str) == True or 
            (isinstance(cloudsat_files, list) and len(cloudsat_files) != 0)):
            write_log("INFO","Read CLOUDSAT %s data" % config.CLOUDSAT_TYPE)
            cl_matchup, cl_time_diff = get_cloudsat_matchups(cloudsat_files, 
                                                             pps_files.cloudtype,
                                                             avhrrGeoObj, avhrrObj, ctype,
                                                             ctth, surft, avhrrAngObj, cppLwp, config_options)
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

                if len(calipso_files) == 0:
                    raise MatchupError("Did not find any matching 5km and 1km calipso files")
                if len(calipso_files) != len(calipso5km):
                    raise MatchupError("Inconsistent number of calipso files...\n" + 
                                       "\tlen(calipso_files) = %d\n" % len(calipso_files) + 
                                       "\tlen(calipso1km) = %d" % len(calipso5km))
            else:
                calipso5km = None
                
        write_log("INFO", "Read CALIPSO data")        
        ca_matchup, ca_time_diff = get_calipso_matchups(calipso_files, 
                                                        values,
                                                        avhrrGeoObj, avhrrObj, ctype,
                                                        ctth, surft, avhrrAngObj, 
                                                        config_options, cppCph,  
                                                        calipso1km, calipso5km)

    else:
        write_log("INFO", "NO CALIPSO File, Continue")
    #import pdb; pdb.set_trace()

    # Get satellite name, time, and orbit number from avhrr_file
    date_time = values["date_time"]
    basename = '_'.join(os.path.basename(avhrr_file).split('_')[:4])
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
        date_time = cross.time1
    except AttributeError:
        raise ValueError('cross is not a valid SNO cross. (cross: %s)' % cross)
    
    if values["satellite"] in ['calipso', 'cloudsat']:
        values["satellite"] = cross.satellite2.lower()
        date_time = cross.time2
        
    if reprocess is False:
        #need to pyt in the info res, atrain data type before go inte find avhrr??
        # or change read of files
        values["atrain_sat"] = "cloudsat-%s" % config.CLOUDSAT_TYPE[0]
        values["atrain_datatype"] = "cloudsat-%s" % config.CLOUDSAT_TYPE[0]
        cl_match_file, tobj = find_avhrr_file(cross, options['reshape_dir'], options['reshape_file'], values=values)
        if not cl_match_file:
            write_log('INFO', "No processed CloudSat match files found." + 
                      " Generating from source data.")
            clObj = None
        else:
            date_time=tobj
            clObj = readCloudsatAvhrrMatchObj(cl_match_file) 
            basename = '_'.join(os.path.basename(cl_match_file).split('_')[1:5])
            write_log('INFO', "CloudSat Matchups read from previously " + 
                      "processed data.")
        
        values["atrain_sat"] = "caliop"
        values["atrain_datatype"] = "caliop"
        ca_match_file, tobj = find_avhrr_file(cross, options['reshape_dir'], options['reshape_file'], values=values)
        if not  ca_match_file:
            write_log('INFO', 
                      ("No processed CALIPSO match files found. "+
                       "Generating from source data."))
            caObj = None
        else:
            #print ca_match_file
            date_time=tobj
            caObj = readCaliopAvhrrMatchObj(ca_match_file)
            basename = '_'.join(os.path.basename(ca_match_file).split('_')[1:5])
            write_log('INFO', 
                      "CALIPSO Matchups read from previously processed data.")
            write_log('INFO', 'Filename: ' + ca_match_file)
            calipso_min_and_max_timediffs = (caObj.diff_sec_1970.min(), 
                                             caObj.diff_sec_1970.max())
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


def get_cloud_emissivity(satellite, calipsoObj, calipso_okay, radtb_table_path=None):
    """Derive the cloud emissivity for each Calipso matchup"""
    if satellite in ['noaa17', 'noaa18', 'noaa19']:
        platform = "noaa"
        satnumber = int(satellite[4:6])
    elif satellite in ['metop02', 'metop01', 'metop03']:
        platform = "metop"
        satnumber = int(satellite[5:7])
    elif satellite in ['npp']:
        platform = 'npp'
        satnumber = 1
    else:
        raise NotImplementedError("Support for satellite %s is not yet implemented." % satellite)

    cloud_e = np.zeros(calipso_okay.shape[0], 'd')
    radtbObj = None
    if radtb_table_path:
        from rad_tb_tables import radtbTable
        radtbObj = radtbTable(platform, satnumber, 
                              INSTRUMENT.get(platform, 'avhrr'), 
                              path=radtb_table_path, channel='M15',
                              detector_number=1)
        radtbObj.read()

    elif platform not in ["noaa"]:
        raise NotImplementedError("Support for platform " + 
                                  "%s is not yet implemented." % platform)

    midlayer_temp_kelvin = calipsoObj.calipso.cloud_mid_temperature + 273.15
    caliop_max_height_midlaytemp = np.where(np.greater(calipsoObj.calipso.cloud_top_profile[0, ::], -9), 
                                            midlayer_temp_kelvin[0, ::], -9)

    if radtbObj:
        for i in range(calipso_okay.shape[0]):
            if calipso_okay[i] and (calipsoObj.avhrr.surftemp[i] < 360. 
                                    and calipsoObj.avhrr.surftemp[i] > 180.):                
                radiance = radtbObj.get_radiance(calipsoObj.avhrr.bt11micron[i])
                rad_clear = radtbObj.get_radiance(calipsoObj.avhrr.surftemp[i])
                if radiance > rad_clear:
                    rad_clear = radiance   # Just avoiding too much mismatch
                                           # between forecasted and real
                                           # surface temps
                btc = radtbObj.get_radiance(caliop_max_height_midlaytemp[i])
                
                if btc < rad_clear:
                    cloud_e[i] = (radiance - rad_clear)/(btc - rad_clear)
                else:
                    cloud_e[i]=-9.0 # Give up on all temperature inversion cases!
            else:
                cloud_e[i]=-9.0

        return cloud_e

    dum1, cwnum, dum2 = get_central_wavenumber(satnumber, 273.15)
    if calipsoObj.avhrr.surftemp != None:
        for i in range(calipso_okay.shape[0]):
            if calipso_okay[i]:
                radiance = tb2radiance_using_central_wavenumber_klm(cwnum, 
                                                                    calipsoObj.avhrr.bt11micron[i], 
                                                                    satnumber,'4')
                rad_clear = tb2radiance_using_central_wavenumber_klm(cwnum,
                                                                     calipsoObj.avhrr.surftemp[i],
                                                                     satnumber,'4')
                if radiance > rad_clear:
                    rad_clear = radiance   # Just avoiding too much mismatch
                                           # between forecasted and real
                                           # surface temps
                btc = tb2radiance_using_central_wavenumber_klm(cwnum,
                                                               caliop_max_height_midlaytemp[i],
                                                               satnumber,'4')
                if btc < rad_clear:
                    cloud_e[i] = (radiance - rad_clear)/(btc - rad_clear)
                else:
                    cloud_e[i]=-9.0 # Give up on all temperature inversion cases!
            else:
                cloud_e[i]=-9.0

    return cloud_e



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
    print values
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
            clsatObj = CloudsatAvhrrSatz(clsatObj)
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
        avhrr_ctth_csat_ok = np.where(np.greater(avhrr_ctth_csat_ok,0.0),avhrr_ctth_csat_ok[::]+elevation*1.0,avhrr_ctth_csat_ok)
        if len(data_ok) == 0:
            write_log('INFO',"Processing stopped: Zero lenght of matching arrays!")
            print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit()
    else:
        data_ok = None
        avhrr_ctth_csat_ok = None

    ## Calipso ##        
    caObj = CalipsoAvhrrSatz(caObj)
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
    avhrr_ctth_cal_ok = np.where(np.greater(avhrr_ctth_cal_ok,0.0),avhrr_ctth_cal_ok[::]+cal_elevation,avhrr_ctth_cal_ok)                    

    if (len(cal_data_ok) == 0):
        write_log('INFO', "Processing stopped: Zero lenght of matching arrays!")
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit()
        
    # If everything is OK, now create filename for statistics output file and
    # open it for writing.  Notice that more than one file (but maximum 2) can
    # be created for one particular noaa orbit.

    #resultpath = "%s/%s/%ikm/%s/%s/%s/%s" % (config.RESULT_DIR, 
    #                                         base_sat,
    #                                         int(config.RESOLUTION), 
    #                                         base_year, base_month, 
    #                                         config.AREA, process_mode_dnt)
    min_depth_to_file_name=""
    if process_mode == 'OPTICAL_DEPTH':
        min_depth_to_file_name="-%.1f"%(min_optical_depth)
        #resultpath = "%s/%s/%ikm/%s/%s/%s/%s-%.1f" % (config.RESULT_DIR, 
        #                                              base_sat,
        #                                              int(config.RESOLUTION), 
        #                                              base_year, base_month,
        #                                              config.AREA, process_mode_dnt, 
        #                                              min_optical_depth)      
    #result_dir = {val_dir}/Results/{satellite}/{resolution}km/{year}/{month}/{area}/{mode}-{min_opt_depth}/
    values['mode']= process_mode_dnt
    values['min_opt_depth']=min_depth_to_file_name
    result_path = insert_info_in_filename_or_path(config_options['result_dir'], values, datetime_obj=values['date_time'])
    #result_path = config_options['result_dir'].format(
    #        val_dir=config._validation_results_dir,
    #        satellite=values["satellite"],
    #        resolution=str(config.RESOLUTION),
    #        area=config.AREA,
    #        month=values["month"],
    #        year=values["year"],
    #        mode=process_mode_dnt,
    #        min_opt_depth=min_depth_to_file_name
    #        )
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    #statname = "%ikm_%s_cloudsat_calipso_avhrr_stat.dat" % (int(config.RESOLUTION),
    #                                                        basename)
    result_file = config_options['result_file'].format(
            resolution=str(config.RESOLUTION),
            basename=values['basename'] )

    statfilename = os.path.join(result_path, result_file)

    if CALIPSO_CLOUD_FRACTION == True:
        statfilename = "%s/CCF_%ikm_%s_cloudsat_calipso_avhrr_stat.dat" % (resultpath,int(config.RESOLUTION),basename)        
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

##    if CALIPSO_CLOUD_FRACTION == True: ==> Let's skip the use of ssccf for calculating cloud fraction, done in another way/KG
####        (new_cloud_top, new_cloud_base, new_optical_depth, new_cloud_fraction, new_fcf, new_ssccf) = \
##            CalipsoCloudFraction(caObj.calipso.cloud_top_profile, \
##                                 caObj.calipso.cloud_base_profile, \
##                                 caObj.calipso.optical_depth, \
##                                 caObj.calipso.cloud_fraction, \
##                                 caObj.calipso.feature_classification_flags, \
##                                 caObj.calipso.single_shot_cloud_cleared_fraction)
##        caObj.calipso.cloud_top_profile = new_cloud_top
##        caObj.calipso.cloud_base_profile = new_cloud_base
##        caObj.calipso.optical_depth = new_optical_depth
##        caObj.calipso.cloud_fraction = new_cloud_fraction
##        caObj.calipso.feature_classification_flags = new_fcf
##        caObj.calipso.single_shot_cloud_cleared_fraction = new_ssccf



    # If mode = OPTICAL_DEPTH -> Change cloud -top and -base profile
    if process_mode == 'OPTICAL_DEPTH':#Remove this if-statement if you always want to do filtering!/KG
        (new_cloud_top, new_cloud_base, new_cloud_fraction, new_fcf) = \
                        CloudsatCloudOpticalDepth(caObj.calipso.cloud_top_profile, caObj.calipso.cloud_base_profile, \
                                                  caObj.calipso.optical_depth, caObj.calipso.cloud_fraction, caObj.calipso.feature_classification_flags, min_optical_depth)
        caObj.calipso.cloud_top_profile = new_cloud_top
        caObj.calipso.cloud_base_profile = new_cloud_base
        caObj.calipso.cloud_fraction = new_cloud_fraction
        caObj.calipso.feature_classification_flags = new_fcf
        
    # Extract CALIOP Vertical Feature Classification (bits 10-12) from 16 bit representation
    # for topmost cloud layer
    cal_vert_feature = np.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9
    feature_array = 4*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[0,::],11),1) + 2*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[0,::],10),1) + np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[0,::],9),1)
    NinaTestar = False    
    if NinaTestar :   
        new_cloud_top = NinaTestarMedelCloudBaseAndTop(caObj.calipso.cloud_top_profile, caObj.calipso.cloud_base_profile)
        caObj.calipso.cloud_top_profile = new_cloud_top


    # Extract CALIOP Vertical Feature Classification (bits 10-12) from 16 bit
    # representation for topmost cloud layer
    cal_vert_feature = np.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9
    feature_array = 4*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[0,::],11),1) + 2*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[0,::],10),1) + np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[0,::],9),1)

    cal_vert_feature = np.where(np.not_equal(caObj.calipso.feature_classification_flags[0,::],1),feature_array[::],cal_vert_feature[::])   
    # Prepare for plotting, cloud emissivity and statistics calculations
    
    caliop_toplay_thickness = np.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9

    thickness = (caObj.calipso.cloud_top_profile[0,::]-caObj.calipso.cloud_base_profile[0,::])*1000.
    caliop_toplay_thickness = np.where(np.greater(caObj.calipso.cloud_top_profile[0,::],-9),thickness,caliop_toplay_thickness)

    
    caliop_height = []
    caliop_base = []
    caliop_max_height = np.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9
    
    for i in range(10):
        hh = np.where(np.greater(caObj.calipso.cloud_top_profile[i,::],-9),
                            caObj.calipso.cloud_top_profile[i,::] * 1000.,-9)
                                        
        caliop_max_height = np.maximum(caliop_max_height,
                                       caObj.calipso.cloud_top_profile[i, ::] * 1000.)
        # This is actually unnecessary - we know that layer 1 is always the
        # highest layer!!  However, arrays caliop_height and caliop_base are
        # needed later for plotting/ KG

        bb = np.where(np.greater(caObj.calipso.cloud_base_profile[i,::],-9),
                            caObj.calipso.cloud_base_profile[i,::] * 1000.,-9)
        #if (hh>bb):
        caliop_height.append(hh)
        caliop_base.append(bb)
        thickness = hh - bb
        
    x = np.repeat(caObj.calipso.number_of_layers_found.ravel(),
                        np.greater(caObj.calipso.number_of_layers_found.ravel(),0))
    #print "Number of points with more than 0 layers: ",x.shape[0]
    cal_data_ok = np.greater(caliop_max_height,-9.)

    if caObj.avhrr.surftemp != None:
        cal_surftemp_ok = np.repeat(caObj.avhrr.surftemp[::],cal_data_ok)

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
                
    if process_mode == 'EMISSFILT': # Apply only when cloud heights are above EMISS_MIN_HEIGHT
        # Now calculate cloud emissivity emiss_cloud for topmost CALIOP and CloudSat layer
        config_path = os.environ.get('IMAGERRAD_HOME', './etc') + "/Tables"
        emiss_cloud = get_cloud_emissivity(sno_satname, caObj, cal_data_ok, config_path)
        write_log('INFO', "Emissivity filtering applied!")

        caliop_min_height_ok = np.greater(caliop_max_height, config.EMISS_MIN_HEIGHT)
        emissfilt_calipso_ok = np.logical_or(np.logical_and(np.greater(emiss_cloud, config.EMISS_LIMIT),caliop_min_height_ok),np.logical_or(np.equal(caliop_max_height,-9.),np.less_equal(caliop_max_height, config.EMISS_MIN_HEIGHT)))

    ##########################
    ### 1 KM DATA CWC-RVOD ###                       
    elif config.RESOLUTION == 1 and cloudsat_type == 'CWC-RVOD' and clsatObj is not None:
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

        #plot_dir = {val_dir}/Plot/{satellite}/{resolution}km/%Y/%m/{area}/  
        plotpath = insert_info_in_filename_or_path(config_options['plot_dir'], values, datetime_obj=values['date_time'])      
        #print plotpath
        #plotpath = os.path.join(config.PLOT_DIR,
        #                        base_sat, "%ikm" % config.RESOLUTION, 
        #                        base_year, base_month, config.AREA)
            
        trajectorypath = os.path.join(plotpath, "trajectory_plot")
        if not os.path.exists(trajectorypath):
            os.makedirs(trajectorypath)

        trajectoryname = os.path.join(trajectorypath, 
                                      "%skm_%s_trajectory" % (int(config.RESOLUTION),
                                                              values['basename']))
        # To make it possible to use the same function call to drawCalClsatGEOPROFAvhrr*kmPlot
        # in any processing mode:
        if 'emissfilt_calipso_ok' not in locals():
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
        if process_mode == 'EMISSFILT':
            process_calipso_ok = emissfilt_calipso_ok
        else:
            process_calipso_ok = 0
        
    write_log('INFO', "Calculating statistics")
    CalculateStatistics(process_mode, clsatObj, statfile, caObj, cal_MODIS_cflag,
                        cal_vert_feature, avhrr_ctth_csat_ok, data_ok,
                        cal_data_ok, avhrr_ctth_cal_ok, caliop_max_height,
                        process_calipso_ok, dnt_flag)
