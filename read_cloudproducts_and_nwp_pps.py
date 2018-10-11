import numpy as np
import os
import netCDF4
import h5py
import logging
import time
import calendar
from datetime import datetime
logger = logging.getLogger(__name__)
from config import (NODATA, PPS_FORMAT_2012_OR_EARLIER,  CTTH_TYPES, 
                    RESOLUTION)
from runutils import do_some_geo_obj_logging
from common import InputError
ATRAIN_MATCH_NODATA = NODATA
#logger.debug('Just so you know: this module has a logger...')

DSEC_PER_AVHRR_SCALINE = 1.0/6. * 4 # A "work for the time being" solution.
if RESOLUTION == 1:
    DSEC_PER_AVHRR_SCALINE = 1.0/6. # Full scan period, i.e. the time
                                    # interval between two consecutive
                                    # lines (sec)
def get_satid_datetime_orbit_from_fname_pps(avhrr_filename,as_oldstyle=False):
    #import runutils
    #satname, _datetime, orbit = runutils.parse_scene(avhrr_filename)
    #returnd orbit as int, loosing leeding zeros, use %05d to get it right.
    # Get satellite name, time, and orbit number from avhrr_file
    if PPS_FORMAT_2012_OR_EARLIER or as_oldstyle:
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
        values['basename'] = (values["satellite"] + 
                              "_" + values["date"] + 
                              "_" + values["time"] + 
                              "_" + values["orbit"])
    else: #PPS v2014-filenames
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
        values['basename'] = (values["satellite"] + 
                              "_" + values["date"] + 
                              "_" + values["time"] + 
                              "_" + values["orbit"])
    values["jday"]=date_time.timetuple().tm_yday

    return values    
        

def createAvhrrTime(Obt, values=None, Trust_sec_1970=False):
    """ Function to make crate a matrix with time for each pixel 
    from objects start adn end time """
    from config import IMAGER_INSTRUMENT 
    #filename = os.path.basename(filename)
    # Ex.: npp_20120827_2236_04321_satproj_00000_04607_cloudtype.h5
    if IMAGER_INSTRUMENT == 'viirs':
    #if filename.split('_')[0] == 'npp':
        if Obt.sec1970_start < 0: #10800
            logger.warning( 
                      "NPP start time negative! " + str(Obt.sec1970_start))
            datetime=values["date_time"]
            Obt.sec1970_start = calendar.timegm(datetime.timetuple())
            #Obt.sec1970_start = calendar.timegm((year, mon, day, hour, mins, sec)) + hundredSec
        num_of_scan = np.int(Obt.num_of_lines / 16.)      
        linetime = np.linspace(Obt.sec1970_start, Obt.sec1970_end, num_of_scan)
        Obt.time = np.repeat(linetime, 16)
        logger.info("NPP start time :  %s", time.gmtime(Obt.sec1970_start))
        logger.info("NPP end time : %s", time.gmtime(Obt.sec1970_end))
    elif Trust_sec_1970  :
        Obt.time = np.linspace(Obt.sec1970_start, Obt.sec1970_end, Obt.num_of_lines)
    else:
        if Obt.sec1970_end < Obt.sec1970_start:
            """
            In some GAC edition the end time is negative. If so this if statement 
            tries to estimate the endtime depending on the start time plus number of 
            scanlines multiplied with the estimate scan time for the instrument. 
            This estimation is not that correct but what to do?
            """
            logger.warning("I need DSEC_PER_AVHRR_SCALINE, to continue ")
            Obt.sec1970_end = int(DSEC_PER_AVHRR_SCALINE * Obt.num_of_lines + Obt.sec1970_start)
        datetime=values["date_time"]
        sec1970_start_filename = calendar.timegm(datetime.timetuple())
        diff_filename_infile_time = sec1970_start_filename-Obt.sec1970_start
        diff_hours= abs( diff_filename_infile_time/3600.0  )
        if (diff_hours<13):
            logger.debug("Time in file and filename do agree. Difference  %d hours.", diff_hours)
        if (diff_hours>13):
            """
            This if statement takes care of a bug in start and end time, 
            that occurs when a file is cut at midnight
            Former condition needed line number in file name:
            if ((values["ppsfilename"].split('_')[-3] != '00000' and PPS_FORMAT_2012_OR_EARLIER) or
            (values["ppsfilename"].split('_')[-2] != '00000' and not PPS_FORMAT_2012_OR_EARLIER)):
            Now instead check if we aer more than 13 hours off. 
            If we are this is probably the problem, do the check and make sure results are fine afterwards.
            """
            logger.warning("Time in file and filename do not agree! Difference  %d hours.", diff_hours)
            timediff = Obt.sec1970_end - Obt.sec1970_start
            old_start = time.gmtime(Obt.sec1970_start + (24 * 3600)) # Adds 24 h to get the next day in new start
            new_start = calendar.timegm(time.strptime('%i %i %i' %(old_start.tm_year, \
                                                                   old_start.tm_mon, \
                                                                   old_start.tm_mday), \
                                                                   '%Y %m %d'))
            Obt.sec1970_start = new_start
            Obt.sec1970_end = new_start + timediff
            diff_filename_infile_time = sec1970_start_filename-Obt.sec1970_start
            diff_hours= abs( diff_filename_infile_time/3600.0)
            if (diff_hours>20):
                raise InputError("Time in file and filename do not agree.",
                                 " %d s diff"%(diff_hours*60.0))        
        Obt.time = np.linspace(Obt.sec1970_start, Obt.sec1970_end, Obt.num_of_lines)
    return Obt
    

class NWPObj(object):
    def __init__(self, array_dict):
        self.surftemp = None
        self.t500 = None
        self.t700 = None
        self.t850 = None
        self.t950 = None
        self.t1000 = None
        self.t900 = None
        self.t250 = None
        self.t800 = None
        self.psur = None
        self.ptro = None
        self.t2m = None
        self.ciwv = None
        self.text_r06 = None
        self.text_t11 = None
        self.text_t37t12 = None
        self.text_t11t12 = None
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
        self.emis6 = None
        self.emis8 = None
        self.emis9 = None
        self.snowa = None
        self.snowd = None
        self.seaice = None
        self.r37 = None
        self.landuse = None
        self.fractionofland = None
        self.elevation = None
        self.__dict__.update(array_dict) 

class smallDataObject(object):
    def __init__(self):
        self.data=None
        self.gain = 1.0
        self.intercept = 0.0
        self.no_data = 0.0
        self.missing_data = 0.0
class imagerAngObj(object):
    def __init__(self):        
        self.satz = smallDataObject()
        self.sunz = smallDataObject()
        self.satazimuth = smallDataObject()
        self.sunazimuth = smallDataObject()
        self.azidiff = smallDataObject()

class imagerGeoObj(object):
    def __init__(self):
        self.latitude = None
        self.longitude = None
        self.azidiff = None
        self.num_of_lines = None

class NewImagerData:
    def __init__(self):
        self.des = ""
        self.no_data = -999
        self.missing_data = -999
        self.nodata = -999
        self.missingdata = -999
        self.channel = []

class ImagerChannelData:
    def __init__(self):
        self.channel = ""
        self.des = ""
        self.gain = 1.0
        self.intercept = 0.0
        self.SZA_corr_done = False
        self.data = None
class CmaObj:
    #skeleton container for v2014 cma
    def __init__(self):
        self.cma_ext = None
        self.cma_prob = None
        self.cma_aflag = None
        self.cma_testlist0 = None #new in v2018
        self.cma_testlist1 = None #new in v2018
        self.cma_testlist2 = None #new in v2018
        self.cma_testlist3 = None #new in v2018
        self.cma_testlist4 = None #new in v2018
        self.cma_testlist5 = None #new in v2018

class CtypeObj:
    #skeleton container for v2014 cloudtype
    def __init__(self):
        self.cloudtype = None
        self.ct_statusflag = None
        self.ct_quality = None
        self.ct_conditions = None
        self.phaseflag = None #v2012
class CtthObj:
    #skeleton container for v2014 cloudtype
    def __init__(self):
        self.height = None
        self.height_corr = None
        self.temperature = None
        self.pressure = None
        self.ctth_statusflag = None
class CppObj:
    #skeleton container for v2018 cpp
    def __init__(self):
        self.cpp_cot = None
        self.cpp_cwp = None
        self.cpp_dcot = None
        self.cpp_dcwp = None
        self.cpp_lwp = None
        self.cpp_iwp = None
        self.cpp_phase = None
        self.cpp_phase_extended = None
        self.cpp_reff = None

#def read_ctth_v2012 might be needed?
def read_ctth_h5(filename):
    h5file = h5py.File(filename, 'r')
    ctth = CtthObj()
    ctth.height = h5file['ctth_alti'].value.astype(np.float)
    ctth.temperature = h5file['ctth_tempe'].value.astype(np.float)
    ctth.pressure = h5file['ctth_pres'].value.astype(np.float)
    ctth.ctth_statusflag = h5file['ctth_status_flag'].value
    ctth.h_gain = h5file['ctth_alti'].attrs['scale_factor']
    ctth.h_intercept = h5file['ctth_alti'].attrs['add_offset']
    ctth.t_gain = h5file['ctth_tempe'].attrs['scale_factor']
    ctth.t_intercept = h5file['ctth_tempe'].attrs['add_offset']
    ctth.p_gain = h5file['ctth_pres'].attrs['scale_factor']
    ctth.p_intercept = h5file['ctth_pres'].attrs['add_offset']
    ctth.h_nodata = h5file['ctth_alti'].attrs['_FillValue']
    ctth.t_nodata = h5file['ctth_tempe'].attrs['_FillValue']
    ctth.p_nodata = h5file['ctth_pres'].attrs['_FillValue']
    hmask = ctth.height == ctth.h_nodata
    tmask = ctth.temperature == ctth.t_nodata
    pmask = ctth.pressure == ctth.p_nodata
    ctth.height[~hmask] = ctth.height[~hmask]*ctth.h_gain + ctth.h_intercept
    ctth.pressure[~pmask] = ctth.pressure[~pmask]*ctth.p_gain + ctth.p_intercept
    ctth.temperature[~tmask] = ctth.temperature[~tmask]*ctth.t_gain + ctth.t_intercept
    ctth.height[hmask] = ATRAIN_MATCH_NODATA
    ctth.pressure[pmask] = ATRAIN_MATCH_NODATA    
    ctth.temperature[tmask] = ATRAIN_MATCH_NODATA 

    logger.debug("min-h: %d, max-h: %d, h_nodata: %d",
                 np.min(ctth.height), np.max(ctth.height), ctth.h_nodata)
    h5file.close()
    return ctth

def read_ctth_nc(filename):
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    ctth = CtthObj()
    ctth.height = pps_nc.variables['ctth_alti'][0,:,:].astype(np.float)
    ctth.temperature = pps_nc.variables['ctth_tempe'][0,:,:].astype(np.float)
    ctth.pressure = pps_nc.variables['ctth_pres'][0,:,:].astype(np.float)
    ctth.ctth_statusflag = pps_nc.variables['ctth_status_flag'][0,:,:]
    #Currently unpacked arrays later in calipso.py
    ctth.h_gain=1.0 
    ctth.h_intercept=0.0
    ctth.p_gain=1.0 #
    ctth.p_intercept=0.0
    ctth.t_gain=1 #
    ctth.t_intercept=0.0
    ctth.h_nodata = pps_nc.variables['ctth_alti']._FillValue
    ctth.t_nodata = pps_nc.variables['ctth_tempe']._FillValue
    ctth.p_nodata = pps_nc.variables['ctth_pres']._FillValue
    logger.debug("min-h: %d, max-h: %d, h_nodata: %d",
                 np.min(ctth.height), np.max(ctth.height.data), ctth.h_nodata)
    #already scaled
    if np.ma.is_masked(ctth.height):     
        ctth.height.data[ctth.height.mask]  = ATRAIN_MATCH_NODATA
        ctth.height = ctth.height.data
    if np.ma.is_masked(ctth.pressure):       
        ctth.pressure.data[ctth.pressure.mask]  = ATRAIN_MATCH_NODATA 
        ctth.pressure = ctth.pressure.data
    if np.ma.is_masked(ctth.temperature):        
        ctth.temperature.data[ctth.temperature.mask]  = ATRAIN_MATCH_NODATA
        ctth.temperature = ctth.temperature.data
    pps_nc.close()
    return ctth

def read_cloudtype_h5(filename):
    h5file = h5py.File(filename, 'r')
    ctype = CtypeObj()
    ctype.cloudtype = h5file['ct'].value
    ctype.ct_conditions = h5file['ct_conditions'].value
    ctype.ct_statusflag = h5file['ct_status_flag'].value
    ctype.ct_quality = h5file['ct_quality'].value
    h5file.close()
    return ctype

def read_cloudtype_nc(filename):
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    ctype = CtypeObj()
    ctype.cloudtype = pps_nc.variables['ct'][0,:,:]
    ctype.ct_conditions = pps_nc.variables['ct_conditions'][0,:,:]
    ctype.ct_statusflag = pps_nc.variables['ct_status_flag'][0,:,:]
    ctype.ct_quality = pps_nc.variables['ct_quality'][0,:,:]
    pps_nc.close()
    return ctype

def read_cma_h5(filename):
    h5file = h5py.File(filename, 'r')
    cma = CmaObj()
    if 'cma_extended' not in h5file.keys():
        if 'cloud_probability' in h5file.keys():
            logger.error("This CMA-file seem lika a CMAPROB file!")
    cma.cma_ext = h5file['cma_extended'].value
    cma.cma_bin = 0*cma.cma_ext.copy()
    cma.cma_bin[cma.cma_ext==1] = 1.0
    cma.cma_bin[cma.cma_ext==2] = 1.0
    cma.cma_bin[cma.cma_ext<0] =  ATRAIN_MATCH_NODATA
    cma.cma_bin[cma.cma_ext>10] =  ATRAIN_MATCH_NODATA

    #try KeyError 'cma'
    h5file.close()
    return cma

def read_cmaprob_h5(filename, cma):
    h5file = h5py.File(filename, 'r')
    if cma is None:
        cma = CmaObj()
    if 'cma_extended' in h5file.keys():
        if 'cloud_probability' not in h5file.keys():
            logger.error("\n This file looks like a normal CMA-file (not CMAPROB) %s", filename)
    name = 'cloud_probability'
    if name not in h5file.keys():
        logger.info("This CMA-file is old.")
        name = "cmaprob"
    cma.cma_prob = h5file[name].value
    h5file.close()
    return cma

def read_cmaprob_nc(filename, cma):
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    if cma is None:
        cma = CmaObj()
    if 'cma_extended' in pps_nc.variables.keys():
        if 'cmaprob' not in pps_nc.variables.keys():
            logger.error("\n This file looks not like a CMAPROB-file like a normal CMA-file %s", filename)
    cma.cma_prob = pps_nc.variables['cmaprob'][0,:,:]
    pps_nc.close()
    return cma

def read_cma_nc(filename):
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    cma = CmaObj()
    if 'cma_extended' not in pps_nc.variables.keys():
        if 'cloud_probability' in pps_nc.variables.keys():
            logger.error("Probably you shold set CMAP_PROB_VALIDATION=True!")
    cma.cma_ext = pps_nc.variables['cma_extended'][0,:,:]
    cma.cma_bin = 0*cma.cma_ext.copy()
    cma.cma_bin[cma.cma_ext==1] = 1.0
    cma.cma_bin[cma.cma_ext==2] = 1.0
    cma.cma_bin[cma.cma_ext<0] =  ATRAIN_MATCH_NODATA
    cma.cma_bin[cma.cma_ext>10] =  ATRAIN_MATCH_NODATA

    #cma.cma_testlist0 = pps_nc.variables['cma_testlist0'][0,:,:]
    #cma.cma_testlist1 = pps_nc.variables['cma_testlist1'][0,:,:]
    #cma.cma_testlist2 = pps_nc.variables['cma_testlist2'][0,:,:]
    #cma.cma_testlist3 = pps_nc.variables['cma_testlist3'][0,:,:]
    #cma.cma_testlist4 = pps_nc.variables['cma_testlist4'][0,:,:]
    #cma.cma_testlist5 = pps_nc.variables['cma_testlist5'][0,:,:]
    for var_name in [
            'cma_aerosol', #new updated name
            'cma_aerosolflag',
            'cma_dust',
            'cma_testlist0',
            'cma_testlist1',
            'cma_testlist2',
            'cma_testlist3',
            'cma_testlist4',
            'cma_testlist5']:
        if var_name in pps_nc.variables.keys():
            array = pps_nc.variables[var_name][0,:,:]
            atrain_name = var_name
            if var_name == 'cma_aerosol':
                atrain_name =    'cma_aerosolflag'            
            if np.ma.is_masked(array):
                mask = array.mask
                data = array.data
                data[mask]= 0
                setattr(cma, atrain_name, data)
            else:
                setattr(cma, atrain_name, array)
        else:
            logger.info("No %s in cma file", var_name)
    pps_nc.close()    
    return cma


def readImagerData_nc(pps_nc):
    imager_data = NewImagerData()
    for var in pps_nc.variables.keys():
        if 'image' in var:
            image = pps_nc.variables[var]
            if image.id_tag in ['satzenith', 'sunzenith', 
                                'azimuthdiff', 
                                'sunazimuth', 'satazimuth']:
                continue
                
            logger.debug("reading channel %s", image.description)
            one_channel = ImagerChannelData()
            #channel = image.channel                     
            data_temporary = image[0,:,:]
            if np.ma.is_masked(one_channel.data):
                one_channel.data = data_temporary.data
                one_channel.data[data_temporary.mask] = ATRAIN_MATCH_NODATA
            else:
                one_channel.data = data_temporary
            one_channel.des = image.description
            if hasattr(pps_nc.variables[var], "sun_zenith_angle_correction_applied"):
                corr_done_attr = pps_nc.variables[var].sun_zenith_angle_correction_applied
                if corr_done_attr.upper() in ["TRUE"]:
                    one_channel.SZA_corr_done = True
            one_channel.gain = 1.0 #data already scaled
            one_channel.intercept = 0.0 #data already scaled
            imager_data.channel.append(one_channel) 
            fill_value = pps_nc.variables[var]._FillValue
            imager_data.nodata = fill_value
            imager_data.missingdata = fill_value
            imager_data.no_data = fill_value
            imager_data.missing_data = fill_value
    return imager_data

def read_pps_angobj_nc(pps_nc):
    """Read angles info from filename
    """
    AngObj = imagerAngObj()
    for varname in pps_nc.variables.keys():
        this_is = "non_angle_variable"
        if varname in ['satzenith', 
                       'sunzenith', 
                       'azimuthdiff', 
                       'sunazimuth', 
                       'satazimuth']:
            this_is = varname
        else:    
            #some times we got angels in imager file
            #they are then named imageX as varname
            if hasattr(pps_nc.variables[varname], 'id_tag'):
                this_is = pps_nc.variables[varname].id_tag 
 
        if this_is in['satzenith']:
            AngObj.satz.data = pps_nc.variables[varname][0,:,:].astype(np.float)
            AngObj.satz.no_data = pps_nc.variables[varname]._FillValue
            AngObj.satz.intercept = pps_nc.variables[varname].add_offset
            AngObj.satz.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['sunzenith']:
            AngObj.sunz.data = pps_nc.variables[varname][0,:,:].astype(np.float)
            AngObj.sunz.no_data = pps_nc.variables[varname]._FillValue
            AngObj.sunz.intercept = pps_nc.variables[varname].add_offset
            AngObj.sunz.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['azimuthdiff']:
            AngObj.azidiff.data = pps_nc.variables[varname][0,:,:].astype(np.float)
            AngObj.azidiff.no_data = pps_nc.variables[varname]._FillValue
            AngObj.azidiff.intercept = pps_nc.variables[varname].add_offset
            AngObj.azidiff.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['sunazimuth']:
            AngObj.sunazimuth.data = pps_nc.variables[varname][0,:,:].astype(np.float)
            AngObj.sunazimuth.no_data = pps_nc.variables[varname]._FillValue
            AngObj.sunazimuth.intercept = pps_nc.variables[varname].add_offset
            AngObj.sunazimuth.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['satazimuth']:
            AngObj.satazimuth.data = pps_nc.variables[varname][0,:,:].astype(np.float)
            AngObj.satazimuth.no_data = pps_nc.variables[varname]._FillValue
            AngObj.satazimuth.intercept = pps_nc.variables[varname].add_offset
            AngObj.satazimuth.gain = pps_nc.variables[varname].scale_factor
 
    #already scaled
    if np.ma.is_masked(AngObj.sunz.data):        
        AngObj.sunz.data[AngObj.sunz.data.mask] = ATRAIN_MATCH_NODATA
        AngObj.sunz.data = AngObj.sunz.data.data
    if np.ma.is_masked(AngObj.satz.data): 
        AngObj.satz.data[AngObj.satz.data.mask] = ATRAIN_MATCH_NODATA
        AngObj.satz.data = AngObj.satz.data.data
    if np.ma.is_masked(AngObj.azidiff.data): 
        AngObj.azidiff.data[AngObj.azidiff.data.mask] = ATRAIN_MATCH_NODATA
        AngObj.azidiff.data = AngObj.azidiff.data.data
    if np.ma.is_masked(AngObj.sunazimuth.data): 
        AngObj.sunazimuth.data[AngObj.sunazimuth.data.mask] = ATRAIN_MATCH_NODATA
        AngObj.sunazimuth.data = AngObj.sunazimuth.data.data
    if np.ma.is_masked(AngObj.satazimuth.data): 
        AngObj.satazimuth.data[AngObj.satazimuth.data.mask] = ATRAIN_MATCH_NODATA
        AngObj.satazimuth.data = AngObj.satazimuth.data.data
    return AngObj

def read_pps_angobj_h5(filename):
    """Read angles info from filename
    """
    h5file = h5py.File(filename, 'r')
    AngObj = imagerAngObj() 
  
    for var in h5file.keys():
        if 'image' in var:
            image = h5file[var]     
            if (image.attrs['description'] == "sun zenith angle" or
                image.attrs['description'] == "Solar zenith angle"):
                print "reading sunz"
                AngObj.sunz.data = image['data'].value.astype(np.float)
                AngObj.sunz.gain = image['what'].attrs['gain']
                AngObj.sunz.intercept = image['what'].attrs['offset']
                AngObj.sunz.no_data = image['what'].attrs['nodata']
                AngObj.sunz.missing_data = image['what'].attrs['missingdata']
            elif (image.attrs['description'] == "satellite zenith angle" or
                  image.attrs['description'] == "Satellite zenith angle"):
                AngObj.satz.data = image['data'].value.astype(np.float)
                AngObj.satz.gain = image['what'].attrs['gain']
                AngObj.satz.intercept = image['what'].attrs['offset']
                AngObj.satz.no_data = image['what'].attrs['nodata']
                AngObj.satz.missing_data = image['what'].attrs['missingdata']
            elif (image.attrs['description'] == 
                  "relative sun-satellite azimuth difference angle" or 
                  image.attrs['description'] == 
                  "Relative satellite-sun azimuth angle"):
                AngObj.azidiff.data = image['data'].value.astype(np.float)
                AngObj.azidiff.gain = image['what'].attrs['gain']
                AngObj.azidiff.intercept = image['what'].attrs['offset']
                AngObj.azidiff.no_data = image['what'].attrs['nodata']
                AngObj.azidiff.missing_data = image['what'].attrs['missingdata']
    sunzmask = np.logical_or(AngObj.sunz.data ==AngObj.sunz.no_data,
                             AngObj.sunz.data ==AngObj.sunz.missing_data)
    satzmask = np.logical_or(AngObj.satz.data ==AngObj.satz.no_data,
                             AngObj.satz.data ==AngObj.satz.missing_data)
    diffmask = np.logical_or(AngObj.azidiff.data ==AngObj.azidiff.no_data,
                             AngObj.azidiff.data ==AngObj.azidiff.missing_data)
    AngObj.sunz.data[~sunzmask] = AngObj.sunz.data[~sunzmask]*AngObj.sunz.gain + AngObj.sunz.intercept
    AngObj.satz.data[~satzmask] = AngObj.satz.data[~satzmask]*AngObj.satz.gain + AngObj.satz.intercept
    AngObj.azidiff.data[~diffmask] = AngObj.azidiff.data[~diffmask]*AngObj.azidiff.gain + AngObj.azidiff.intercept
    AngObj.sunz.data[sunzmask] = ATRAIN_MATCH_NODATA
    AngObj.satz.data[satzmask] = ATRAIN_MATCH_NODATA
    AngObj.azidiff.data[diffmask] = ATRAIN_MATCH_NODATA
    h5file.close()
    return AngObj

def read_pps_geoobj_nc(pps_nc):
    """Read geolocation and time info from filename
    """
    GeoObj = imagerGeoObj()


    longitude = pps_nc.variables['lon'][::]
    if np.ma.is_masked(longitude):
        GeoObj.longitude = pps_nc.variables['lon'][::].data
        GeoObj.longitude[longitude.mask] = -999.0
    else: 
        GeoObj.longitude = longitude
    latitude = pps_nc.variables['lat'][::]
    if np.ma.is_masked(latitude):
        GeoObj.latitude = pps_nc.variables['lat'][::].data
        GeoObj.latitude[latitude.mask] = -999.0
    else: 
        GeoObj.latitude = latitude
    GeoObj.nodata = pps_nc.variables['lon']._FillValue
    GeoObj.num_of_lines = GeoObj.latitude.shape[0]
    #import pdb
    #pdb.set_trace()
    time_temp = pps_nc.variables['time'].units #to 1970 s
    seconds = np.float64(pps_nc.variables['time'][::]) #from PPS often 0                        
    if 'T' in time_temp:
        time_obj = time.strptime(time_temp,'seconds since %Y-%m-%dT%H:%M:%S+00:00')
    else:
        time_obj = time.strptime(time_temp,'seconds since %Y-%m-%d %H:%M:%S.%f +00:00')
 
    sec_since_1970 = calendar.timegm(time_obj)

    GeoObj.sec1970_start = (sec_since_1970 +
                            np.float64(np.min(pps_nc.variables['time_bnds'][::])) + 
                            seconds)
    GeoObj.sec1970_end = (sec_since_1970 + 
                          np.float64(np.max(pps_nc.variables['time_bnds'][::])) + 
                          seconds)
    #print type(GeoObj.sec1970_start)
    GeoObj.sec1970_start = np.float64(GeoObj.sec1970_start)
    GeoObj.sec1970_end = np.float64(GeoObj.sec1970_end)
    do_some_geo_obj_logging(GeoObj)
    return  GeoObj

def read_pps_geoobj_h5(filename):
    """Read geolocation and time info from filename
    """
    h5file = h5py.File(filename, 'r')
    GeoObj = imagerGeoObj()
    in_fillvalue1 = h5file['where/lon/what'].attrs['nodata']
    in_fillvalue2 = h5file['where/lon/what'].attrs['missingdata']
    GeoObj.nodata = -999.0
    gain = h5file['where/lon/what'].attrs['gain']
    intercept = h5file['where/lon/what'].attrs['offset']
    GeoObj.longitude = h5file['where/lon']['data'].value*gain + intercept
    GeoObj.latitude = h5file['where/lat']['data'].value*gain + intercept

    GeoObj.longitude[h5file['where/lon']['data'].value==in_fillvalue1] = GeoObj.nodata
    GeoObj.latitude[h5file['where/lat']['data'].value==in_fillvalue1] = GeoObj.nodata
    GeoObj.longitude[h5file['where/lon']['data'].value==in_fillvalue2] = GeoObj.nodata
    GeoObj.latitude[h5file['where/lat']['data'].value==in_fillvalue2] = GeoObj.nodata

    GeoObj.num_of_lines = GeoObj.latitude.shape[0]
    GeoObj.sec1970_start = h5file['how'].attrs['startepochs']
    GeoObj.sec1970_end =  h5file['how'].attrs['endepochs']
    do_some_geo_obj_logging(GeoObj)

    return  GeoObj

def read_cpp_h5(filename):
    density = 1e3
    h5file = h5py.File(filename, 'r')
    cpp_obj = CppObj()    
    for cpp_key in cpp_obj.__dict__.keys():
        data = read_cpp_h5_one_var(h5file, cpp_key)
        if cpp_key in ["cpp_lwp"]:
            logger.debug("Convert from CPP-lwp from kg/m-2 to g/m-2")
            data[data>0] = density * data[data>0] 
        setattr(cpp_obj, cpp_key, data)
    h5file.close()    
    return cpp_obj 

def read_cpp_h5_one_var(h5file, cpp_key):
    if cpp_key in h5file.keys():
        logger.debug("Read %s", cpp_key)
        cpp_var_value = h5file[cpp_key].value
        nodata = h5file[cpp_key].attrs['_FillValue']  
        if cpp_key in ["cpp_phase", "cpp_phase_extended"]:
            gain = 1.0
            intercept = 0.0
        else:    
            gain = h5file[cpp_key].attrs['scale_factor']
            intercept = h5file[cpp_key].attrs['add_offset']

        cpp_data = np.where(cpp_var_value != nodata, 
                           cpp_var_value * gain + intercept, 
                           ATRAIN_MATCH_NODATA) 
        return  cpp_data
    else:
        logger.debug("NO %s field, Continue ", cpp_key)
        return None

def read_cpp_nc_one_var(ncFile, cpp_key):
    density = 1e3
    if cpp_key in ncFile.variables.keys():
        logger.debug("Read %s", cpp_key)
        cpp_var = ncFile.variables[cpp_key][0,:,:]
        if np.ma.is_masked(cpp_var):
            cpp_data = cpp_var.data.astype(np.float)
            cpp_data[cpp_var.mask] = ATRAIN_MATCH_NODATA
        else:
            cpp_data = cpp_var
        if cpp_key in ["cpp_lwp"]:
            logger.debug("Convert from CPP-lwp from to g/m-2")
            cpp_data[cpp_data>0] = density * cpp_data[cpp_data>0] 
        return  cpp_data
    else:
        logger.debug("NO %s field, Continue ", cpp_key)
        return None


def read_cpp_nc(filename):
    pps_nc_cpp = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    cpp_obj = CppObj()    
    for cpp_key in cpp_obj.__dict__.keys():
        data = read_cpp_nc_one_var(pps_nc_cpp, cpp_key)
        setattr(cpp_obj, cpp_key, data)
    pps_nc_cpp.close()
    return cpp_obj    


def read_nwp_h5(filename, nwp_key):

    import h5py 
    h5file = h5py.File(filename, 'r')
    if nwp_key in h5file.keys():
        logger.debug("Read NWP %s", nwp_key)
        value = h5file[nwp_key].value
        gain = h5file[nwp_key].attrs['gain']
        intercept = h5file[nwp_key].attrs['intercept']
        nodat = h5file[nwp_key].attrs['nodata']
        data = np.where(value != nodat,value * gain + intercept, value)
        h5file.close()
        return data    
    else:
        logger.debug("NO NWP %s File, Continue", nwp_key)
        return None

def read_etc_nc(ncFile, etc_key):
    if etc_key in ncFile.variables.keys():
        logger.debug("Read %s", etc_key)
        nwp_var = ncFile.variables[etc_key][0,:,:]
        if np.ma.is_masked(nwp_var):
            if 'emis' in etc_key:
                #set emissivity 1.0 where we miss data
                nwp_var.data[nwp_var.mask] = 1.0
            nwp_var = nwp_var.data
        return  nwp_var
    else:
        logger.debug("NO %s field, Continue ", etc_key)
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
        logger.debug("Read segment info moisture")
        for moist_data in ['moist', 'surfaceMoist']:
            data = h5file['segment_area'][moist_data]
            gain = h5file.attrs['m_gain']
            intercept = h5file.attrs['m_intercept']
            nodata = h5file.attrs['m_nodata']
            #data[data!=nodata] = data[data!=nodata] * (gain) + intercept
            product[moist_data] = data
        logger.debug("Read segment info pressure")
        for pressure_data in ['pressure', 'surfacePressure', 'ptro']:
            #pressure is in Pa in segments file
            data = h5file['segment_area'][pressure_data]
            gain = h5file.attrs['p_gain']
            intercept = h5file.attrs['p_intercept']
            nodata = h5file.attrs['p_nodata']
            #data[data!=nodata] = data[data!=nodata] * (gain/100) + intercept/100 #Pa => hPa
            data_float = np.array(data, dtype=np.float)
            data_float[data!=nodata] = data_float[data!=nodata] * (gain/100) + intercept/100 #Pa => hPa
            product[pressure_data] = data_float
        logger.debug("Read segment info height")
        for geoheight_data in ['geoheight', 'surfaceGeoHeight']:
            #geo height is in meters in segment file
            data = h5file['segment_area'][geoheight_data]
            gain = h5file.attrs['h_gain']
            intercept = h5file.attrs['h_intercept']
            nodata = h5file.attrs['h_nodata']
            data[data!=nodata] = data[data!=nodata] * gain + intercept
            product[geoheight_data] = data
        logger.debug("Read segment info temperature")
        for temperature_data in ['temp']:
            # temperature are measured in Kelvin in segmentfile
            data = h5file['segment_area'][temperature_data]
            gain = h5file.attrs['t_gain']
            intercept = h5file.attrs['t_intercept']
            nodata = h5file.attrs['t_nodata']
            data_float = np.array(data, dtype=np.float)
            data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
            product[temperature_data] = data_float
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
        logger.debug("Read segment info brightness temperature")
        
        for tb_data in ['tb11clfree_sea',
                        'tb12clfree_sea',
                        'tb11clfree_land',
                        'tb12clfree_land',
                        'tb4clfree_sea',
                        'tb5clfree_sea',
                        'tb4clfree_land',
                        'tb5clfree_land'
                        'tb11lowcloud_sea',
                        'tb12lowcloud_sea',
                        'tb11lowcloud_land',
                        'tb12lowcloud_land']:
            try:
                
                data = h5file['segment_area'][tb_data]
                gain = h5file.attrs['t_gain']
                intercept = h5file.attrs['tb_intercept']
                nodata = h5file.attrs['t_nodata']
                data_float = np.array(data, dtype=np.float)
                data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
                name= tb_data
                name = name.replace('4','11')
                name = name.replace('5','12')
                product[name] = data_float                
            except ValueError:
                pass
        h5file.close()
        return product
    else:
        logger.info("NO segment %s File, Continue", filename)
        return None


def read_thr_h5(filename, h5_obj_type, thr_type):
    import h5py 
    product = None
    if thr_type in ["emis1","emis6", "emis8", "emis9"]:
        if filename is not None: 
            h5file = h5py.File(filename, 'r')
            if 1==1:#h5_obj_type in h5file.keys():
                value = h5file[h5_obj_type].value
                gain = h5file.attrs['gain']
                intercept = h5file.attrs['intercept']
                product = value * gain + intercept
                product[product<0] = 1.0
                product[product>1.0] = 1.0
                logger.debug("Read EMIS: %s", thr_type)
            else:
                logger.info("ERROR","Couldn't read %s file, Continue", thr_type)
            h5file.close()   
        else:
            logger.debug("NO EMIS %s File, Continue", thr_type)
        return product  
    if filename is not None: 
        h5file = h5py.File(filename, 'r')
        if h5_obj_type in h5file.keys():
            value = h5file[h5_obj_type].value
            gain = h5file[h5_obj_type].attrs['gain']
            intercept = h5file[h5_obj_type].attrs['intercept']
            product = value * gain + intercept
            logger.debug("Read THR: %s"%(thr_type))
        else:
            logger.error("Could not read %s File, Continue", thr_type)
        h5file.close()   
    else:
        logger.debug("NO THR %s File, Continue", thr_type)
    return product


def readImagerData_h5(filename):
    h5file = h5py.File(filename, 'r')
    imager_data = NewImagerData()
    for var in h5file.keys():
        if 'image' in var:
            image = h5file[var]            
            if 'description' not in image.attrs.keys() and "modis" in filename:
                my_description = "MODIS %s"%(image.attrs['channel'])
            else:
                my_description = image.attrs['description']
            logger.debug("reading channel %s", my_description )
            one_channel = ImagerChannelData()                   
            one_channel.data = image['data'].value
            one_channel.des = my_description
            one_channel.gain = 1.0
            one_channel.intercept = 0.0
            gain = image['what'].attrs['gain']
            intercept = image['what'].attrs['offset']
            imager_data.channel.append(one_channel) 
            imager_data.nodata = image['what'].attrs['nodata']
            imager_data.missingdata = image['what'].attrs['missingdata']
            imager_data.no_data = imager_data.nodata
            imager_data.missing_data = imager_data.missingdata
            mask = np.logical_or(one_channel.data == imager_data.no_data,
                                 one_channel.data == imager_data.missing_data)
            one_channel.data[~mask] = one_channel.data[~mask]*gain + intercept
            one_channel.data[mask] = ATRAIN_MATCH_NODATA
    h5file.close()
    return imager_data


def read_all_intermediate_files(pps_files):
    nwp_dict={}
    if pps_files.seaice is None:
        pass
    elif '.nc' in pps_files.seaice:
        pps_nc_seaice = netCDF4.Dataset(pps_files.seaice, 'r', format='NETCDF4')
        nwp_dict["seaice"] = read_etc_nc(pps_nc_seaice, "seaice")
        pps_nc_seaice.close()
    else:
        logger.info("Not reading PPS seaice data")  
    if pps_files.physiography is None:
        logger.info("Not reading PPS physiography data") 
    elif '.nc' in pps_files.physiography:
        pps_nc_physiography = netCDF4.Dataset(pps_files.physiography, 'r', format='NETCDF4')
        nwp_dict["landuse"] = read_etc_nc(pps_nc_physiography, "landuse")
        nwp_dict["fractionofland"] = read_etc_nc(pps_nc_physiography, "fractionofland")
        nwp_dict["elevation"] = read_etc_nc(pps_nc_physiography, "elevation")
        pps_nc_physiography.close()
    else:
        logger.info("Not reading PPS physiography data") 
    if pps_files.r37 is None:
        pass
    else:
        pps_nc_r37 = netCDF4.Dataset(pps_files.r37, 'r', format='NETCDF4')
        nwp_dict["r37_sza_correction_done"] = read_etc_nc(pps_nc_r37, "r37")
        pps_nc_r37.close()        
    if pps_files.nwp_tsur is None:
        pass
    elif '.nc' in pps_files.nwp_tsur:
        pps_nc_nwp = netCDF4.Dataset(pps_files.nwp_tsur, 'r', format='NETCDF4')
        nwp_dict['surftemp'] = read_etc_nc(pps_nc_nwp, "tsur")
        nwp_dict['t500'] = read_etc_nc(pps_nc_nwp, "t500")
        nwp_dict['t700'] = read_etc_nc(pps_nc_nwp, "t700")
        nwp_dict['t850'] = read_etc_nc(pps_nc_nwp, "t850")
        nwp_dict['t950'] = read_etc_nc(pps_nc_nwp, "t950")
        nwp_dict['ttro'] = read_etc_nc(pps_nc_nwp, "ttro")
        nwp_dict['ciwv'] = read_etc_nc(pps_nc_nwp, "ciwv")
        nwp_dict['t1000'] = read_etc_nc(pps_nc_nwp, "t1000")
        nwp_dict['t900'] = read_etc_nc(pps_nc_nwp, "t900")
        nwp_dict['t800'] = read_etc_nc(pps_nc_nwp, "t800")
        nwp_dict['t250'] = read_etc_nc(pps_nc_nwp, "t250")
        nwp_dict['ptro'] = read_etc_nc(pps_nc_nwp, "ptro")
        nwp_dict['psur'] = read_etc_nc(pps_nc_nwp, "psur")
        nwp_dict['t2m'] = read_etc_nc(pps_nc_nwp, "t2m")
        nwp_dict['snowa'] = read_etc_nc(pps_nc_nwp, "snowa")
        nwp_dict['snowd'] = read_etc_nc(pps_nc_nwp, "snowd")
    else:   
         nwp_dict['surftemp'] = read_nwp_h5(pps_files.nwp_tsur,"tsur")
         nwp_dict['t500'] = read_nwp_h5(pps_files.nwp_t500, "t500")
         nwp_dict['t700'] = read_nwp_h5(pps_files.nwp_t700, "t700")
         nwp_dict['t850'] = read_nwp_h5(pps_files.nwp_t850, "t850")
         nwp_dict['t950'] = read_nwp_h5(pps_files.nwp_t950, "t950")
         nwp_dict['ttro'] = read_nwp_h5(pps_files.nwp_ttro, "ttro")
         nwp_dict['ciwv'] = read_nwp_h5(pps_files.nwp_ciwv, "ciwv")
    if pps_files.text_t11 is None:
        pass
        logger.info("Not reading PPS texture data")  
    elif '.nc' in pps_files.text_t11:
        pps_nc_txt = netCDF4.Dataset(pps_files.text_t11, 'r', format='NETCDF4')
        for ttype in ['r06', 't11', 't37t12', 't37', 't11t12']:
            text_type = 'text_' + ttype
            nwp_dict[text_type] = read_etc_nc(pps_nc_txt, ttype)
        pps_nc_txt.close()
    else:    
        for ttype in ['r06', 't11', 't37t12', 't37']:
            h5_obj_type = ttype +'_text'
            text_type = 'text_' + ttype
            nwp_dict[text_type] = read_thr_h5(getattr(pps_files,text_type), 
                                              h5_obj_type,text_type)
    if pps_files.thr_t11ts is None:
        pass
        logger.info("Not reading PPS threshold data") 
    elif '.nc' in pps_files.thr_t11ts:
        pps_nc_thr = netCDF4.Dataset(pps_files.thr_t11ts, 'r', format='NETCDF4')
        for nc_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv', 
                            't11ts', 't11t37', 't37t12', 't11t12',
                            'r09', 'r06', 't85t11_inv', 't85t11']:
            thr_type = 'thr_' + nc_obj_type
            nwp_dict[thr_type] = read_etc_nc(pps_nc_thr,nc_obj_type)
        pps_nc_thr.close()    
    else:    
        for h5_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv', 
                            't11ts', 't11t37', 't37t12', 't11t12',
                            'r09', 'r06', 't85t11_inv', 't85t11']:
            thr_type = 'thr_' + h5_obj_type
            nwp_dict[thr_type] = read_thr_h5(getattr(pps_files,thr_type), 
                                             h5_obj_type, thr_type)
    if pps_files.emis is None:
        pass
        logger.info("Not reading PPS Emissivity data") 
    elif '.nc' in pps_files.emis:
        pps_nc_thr = netCDF4.Dataset(pps_files.emis, 'r', format='NETCDF4')
        for emis_type in ['emis1',"emis6", 'emis8','emis9']:
            nwp_dict[emis_type] = read_etc_nc(pps_nc_thr, emis_type)
        pps_nc_thr.close()    
    else:
        for h5_obj_type in ['emis1',"emis6", 'emis8','emis9']:
            emis_type = h5_obj_type
            nwp_dict[emis_type] = read_thr_h5(getattr(pps_files,"emis"), 
                                              h5_obj_type, emis_type)
    if len(CTTH_TYPES)>1:        
        for ctth_type in CTTH_TYPES[1:]: #already read first
            nwp_dict[ctth_type] = read_ctth_nc(pps_files.ctth[ctth_type]) 
    nwp_obj = NWPObj(nwp_dict)
    return nwp_obj

def pps_read_all(pps_files, avhrr_file):
    logger.info("Read Imager geolocation data")
    if '.nc' in avhrr_file:
        pps_nc = netCDF4.Dataset(avhrr_file, 'r', format='NETCDF4')
        imagerGeoObj = read_pps_geoobj_nc(pps_nc)
    else:    
        #use mpop?
        imagerGeoObj = read_pps_geoobj_h5(avhrr_file)  
    #create time info for each pixel  
    values = get_satid_datetime_orbit_from_fname_pps(avhrr_file)  
    imagerGeoObj = createAvhrrTime(imagerGeoObj, values)
    logger.info("Read sun and satellites angles data")
    if '.nc' in pps_files.sunsatangles:
        pps_nc_ang = netCDF4.Dataset(pps_files.sunsatangles, 'r', format='NETCDF4')
        avhrrAngObj = read_pps_angobj_nc(pps_nc_ang)
        pps_nc_ang.close()
    else:
        #use mpop?
        avhrrAngObj = read_pps_angobj_h5(pps_files.sunsatangles)
    logger.info("Read Imager data")
    if '.nc' in avhrr_file:
        avhrrObj = readImagerData_nc(pps_nc)
        pps_nc.close()
    else:
        avhrrObj = readImagerData_h5(avhrr_file)

    logger.debug("%s, %s, %s", pps_files.cloudtype, pps_files.ctth, pps_files.cma)

    #CPP
    cpp = None
    if pps_files.cpp is not None:
        logger.info("Read CPP data")
        if '.nc' in pps_files.cpp:          
            cpp = read_cpp_nc(pps_files.cpp)
        else:            
            cpp = read_cpp_h5(pps_files.cpp)    
    #CMA
    cma = None
    if pps_files.cma is not None:
        logger.info("Read PPS Cloud mask")
        logger.debug(pps_files.cma )
        if '.nc' in pps_files.cma:
            cma = read_cma_nc(pps_files.cma)
        else:
            cma = read_cma_h5(pps_files.cma)  
    #CMAPROB
    if pps_files.cmaprob is not None:
        logger.info("Read PPS Cloud mask prob")
        if '.nc' in pps_files.cmaprob:
            cma = read_cmaprob_nc(pps_files.cmaprob, cma)
        else:
            cma = read_cmaprob_h5(pps_files.cmaprob, cma)
    #CTYPE
    ctype = None
    if pps_files.cloudtype is not None:
        logger.info("Read PPS Cloud type")        
        if '.nc' in pps_files.cloudtype:
            ctype = read_cloudtype_nc(pps_files.cloudtype)
        else:
            ctype = read_cloudtype_h5(pps_files.cloudtype)
    #CTTH        
    ctth = None
    if len(pps_files.ctth.keys())>=1:
        logger.info("Read PPS CTTH")
        if '.nc' in pps_files.ctth[CTTH_TYPES[0]]:
            #read first ctth as primary one
            ctth = read_ctth_nc(pps_files.ctth[CTTH_TYPES[0]])
        else:            
            ctth = read_ctth_h5(pps_files.ctth[CTTH_TYPES[0]])

    logger.info("Read PPS full resolution intermediate files")
    nwp_obj = read_all_intermediate_files(pps_files)
  
    logger.info("Read PPS NWP segment resolution data") 
    segment_data_object = read_segment_data(getattr(pps_files,'nwp_segments'))
    return avhrrAngObj, ctth, imagerGeoObj, ctype, avhrrObj, nwp_obj, cpp, segment_data_object, cma 

if __name__ == "__main__":
    pass
