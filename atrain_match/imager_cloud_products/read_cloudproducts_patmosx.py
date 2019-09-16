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
  Use this module to read patmosx cloudproducts
  2018 SMHI, N.Hakansson 
"""

import time
import os
import netCDF4
import numpy as np
import calendar
from datetime import datetime
from datetime import timedelta
import logging
logger = logging.getLogger(__name__)

from imager_cloud_products.read_cloudproducts_and_nwp_pps import (
    AllImagerData, 
    CtypeObj, CtthObj, CmaObj,
    create_imager_time,
    ImagerAngObj)
from utils.runutils import do_some_geo_obj_logging
import config 
ATRAIN_MATCH_NODATA = config.NODATA
#from utils.get_flag_info import get_patmosx_ct_flag, get_day_night_twilight_info_patmosx

def get_satid_datetime_orbit_from_fname_patmosx(imager_filename, SETTINGS, cross):
    # Get satellite name, time, and orbit number from imager_file
    # patmosx_v05r03_NOAA-18_asc_d20090102_c20140317.nc 
    #patmosx_noaa-18_des_2009_122.level2b.hdf

    sl_ = os.path.basename(imager_filename).split('_')
    if ".nc" in imager_filename:
        date_time = datetime.strptime(sl_[4], 'd%Y%m%d')
    else:
        date_time = datetime.strptime(sl_[3]+sl_[4], '%Y%j.level2b.hdf')
    asc_or_des = "asc"
    if "_des_" in imager_filename:
        asc_or_des = "des"

    date_time = cross.time
    #date_time_end = cross.time + timedelta(seconds=SETTINGS['SAT_ORBIT_DURATION'])

    sat_id = sl_[2].replace('-','').lower()
    values = {"satellite": sat_id,
              "date_time": date_time,
              "orbit": "99999",
              "date": date_time.strftime("%Y%m%d"),
              "year": date_time.year,
              "month": "%02d" % (date_time.month),
              "time": date_time.strftime("%H%M"),
              "extrai": asc_or_des,
              #"basename":sat_id + "_" + date_time.strftime("%Y%m%d_%H%M_99999"),#"20080613002200-ESACCI",
              "ccifilename": imager_filename,
              "ppsfilename": None}
    values['basename'] = values["satellite"] + "_" + \
        values["date"] + "_" + values["time"] + "_" + values["orbit"] + "_" + asc_or_des
    return values

def patmosx_read_all(filename, cross, SETTINGS):
    if ".nc" in filename:
        return patmosx_read_all_nc(filename, cross, SETTINGS)
    else:   
        return patmosx_read_all_hdf(filename, cross, SETTINGS)

def patmosx_read_all_nc(filename, cross, SETTINGS):
    """Read geolocation, angles info, ctth, and cma
    """
    patmosx_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    logger.info("Opening file %s", filename)
    logger.info("Reading longitude, latitude and time ...")
    cloudproducts = read_patmosx_geoobj(patmosx_nc, filename, cross, SETTINGS)
    logger.info("Reading angles ...")
    cloudproducts.imager_angles = read_patmosx_angobj(patmosx_nc)
    logger.info("Reading cloud type ...")
    # , angle_obj)
    ctype, cma, ctth = read_patmosx_ctype_cmask_ctth(patmosx_nc)
    cloudproducts.cma = cma
    cloudproducts.ctth = ctth
    cloudproducts.ctype = cype

    logger.info("Not reading surface temperature")
    logger.info("Not reading cloud microphysical properties")
    logger.info("Not reading channel data")
    return cloudproducts

def read_patmosx_ctype_cmask_ctth(patmosx_nc):
    """Read cloudtype and flag info from filename
    """
    ctth = CtthObj()
    ctype = CtypeObj()
    cma = CmaObj()
    ctth.height = patmosx_nc.variables['cld_height_acha'][0,:,:].astype(np.float)
    ctth.h_nodata = patmosx_nc.variables['cld_height_acha']._FillValue
    if np.ma.is_masked(ctth.height):     
        ctth.height.data[ctth.height.mask]  = ATRAIN_MATCH_NODATA
        ctth.height = ctth.height.data
    ctth.height = 1000*ctth.height
    cf = patmosx_nc.variables['cloud_fraction'][0,:,:].astype(np.float)
    cma.cma_ext = np.where(cf>=0.5, 1, 0)
    if np.ma.is_masked(cf):
       cma.cma_ext[cf.mask] = 255 
    ctype.phaseflag = None
    ctype.ct_conditions = None
    return ctype, cma, ctth


def read_patmosx_angobj(patmosx_nc):
    """Read angles info from filename
    """
    angle_obj = ImagerAngObj()
    angle_obj.satz.data = patmosx_nc.variables['sensor_zenith_angle'][0,:,:].astype(np.float)
    angle_obj.sunz.data = patmosx_nc.variables['solar_zenith_angle'][0,:,:].astype(np.float)
    angle_obj.azidiff.data = None
    return angle_obj



def read_patmosx_geoobj(patmosx_nc, filename, cross, SETTINGS):
    """Read geolocation and time info from filename
    """
    cloudproducts = AllImagerData()
    latitude_v = patmosx_nc.variables['latitude'][:].astype(np.float)
    longitude_v = patmosx_nc.variables['longitude'][:].astype(np.float)
    cloudproducts.latitude = np.repeat(latitude_v[:,np.newaxis], len(longitude_v), axis=1)
    cloudproducts.longitude = np.repeat(longitude_v[np.newaxis,:], len(latitude_v), axis=0)
    cloudproducts.nodata=-999
    date_time_start = cross.time
    date_time_end = cross.time + timedelta(seconds=SETTINGS['SAT_ORBIT_DURATION'])
    cloudproducts.sec1970_start = calendar.timegm(date_time_start.timetuple())
    cloudproducts.sec1970_end = calendar.timegm(date_time_end.timetuple())
    frac_hour = patmosx_nc.variables['scan_line_time'][0,:,:].astype(np.float)
    if np.ma.is_masked(frac_hour):
        frac_hour = frac_hour.data
    seconds = frac_hour*60*60.0
    cloudproducts.time = seconds +  patmosx_nc.variables['time']
    do_some_geo_obj_logging(cloudproducts)
    return cloudproducts


def patmosx_read_all_hdf(filename, cross, SETTINGS):
    """Read geolocation, angles info, ctth, and cma
    """
    from pyhdf.SD import SD, SDC
    from pyhdf.HDF import HDF, HC
    import pyhdf.VS 
    patmosx_hdf = SD(filename, SDC.READ)
    logger.info("Opening file %s", filename)
    logger.info("Reading longitude, latitude and time ...")
    cloudproducts = read_patmosx_geoobj_hdf(patmosx_hdf, filename, cross, SETTINGS)
    logger.info("Reading angles ...")
    cloudproducts.imager_angles = read_patmosx_angobj_hdf(patmosx_hdf)
    logger.info("Reading cloud type ...")
    # , angle_obj)
    ctype, cma, ctth = read_patmosx_ctype_cmask_ctth_hdf(patmosx_hdf)
    cloudproducts.ctype = ctype
    cloudproducts.cma = cma
    cloudproducts.ctth = ctth
    logger.info("Not reading surface temperature")
    logger.info("Not reading cloud microphysical properties")
    logger.info("Not reading channel data")
    return imagerAngObj, ctth, imagercloudproducts, ctype,  imager_obj, surft, cpp, cma

def read_patmosx_ctype_cmask_ctth_hdf(patmosx_hdf):
    """Read cloudtype and flag info from filename
    """
    ctth = CtthObj()
    ctype = CtypeObj()
    cma = CmaObj()
    name = 'cld_height_acha'
    height_unscaled =  np.array(patmosx_hdf.select(name).get())
    offset = patmosx_hdf.select(name).attributes()['add_offset']
    gain = patmosx_hdf.select(name).attributes()['scale_factor']
    ctth.height = height_unscaled * gain + offset
    ctth.h_nodata = patmosx_hdf.select(name)._FillValue
    #ctth.height = 1000*ctth.height
    ctth.height[height_unscaled==ctth.h_nodata] = ATRAIN_MATCH_NODATA
    name = 'cloud_fraction'
    cf_unscaled =  np.array(patmosx_hdf.select(name).get())
    offset = patmosx_hdf.select(name).attributes()['add_offset']
    gain = patmosx_hdf.select(name).attributes()['scale_factor']
    cf = cf_unscaled * gain + offset
    cf[cf_unscaled== patmosx_hdf.select(name)._FillValue] =  ATRAIN_MATCH_NODATA
    cma.cma_ext = np.where(cf>=0.5, 1, 0)
    ctype.phaseflag = None
    ctype.ct_conditions = None
    return ctype, cma, ctth


def read_patmosx_angobj_hdf(patmosx_hdf):
    """Read angles info from filename
    """
    angle_obj = ImagerAngObj()
    name = 'sensor_zenith_angle'
    temp = patmosx_hdf.select(name).get().astype(np.float)
    offset = patmosx_hdf.select(name).attributes()['add_offset']
    gain = patmosx_hdf.select(name).attributes()['scale_factor']
    angle_obj.satz.data =  temp * gain + offset
    
    name = 'solar_zenith_angle' 
    temp = patmosx_hdf.select(name).get().astype(np.float)
    offset = patmosx_hdf.select(name).attributes()['add_offset']
    gain = patmosx_hdf.select(name).attributes()['scale_factor']
    angle_obj.sunz.data = temp * gain + offset
    angle_obj.azidiff.data = None
    return angle_obj



def read_patmosx_geoobj_hdf(patmosx_hdf, filename, cross, SETTINGS):
    """Read geolocation and time info from filename
    """
    cloudproducts = AllImagerData()
    name = 'latitude'
    temp = patmosx_hdf.select(name).get().astype(np.float)
    offset = patmosx_hdf.select(name).attributes()['add_offset']
    gain = patmosx_hdf.select(name).attributes()['scale_factor']
    latitude_v = temp * gain + offset
    name = 'longitude'
    temp = patmosx_hdf.select(name).get().astype(np.float)
    offset = patmosx_hdf.select(name).attributes()['add_offset']
    gain = patmosx_hdf.select(name).attributes()['scale_factor']
    longitude_v = temp * gain + offset

    cloudproducts.latitude = np.repeat(latitude_v[:,np.newaxis], len(longitude_v), axis=1)
    cloudproducts.longitude = np.repeat(longitude_v[np.newaxis,:], len(latitude_v), axis=0)
    cloudproducts.nodata=-999
    date_time_start = cross.time
    date_time_end = cross.time + timedelta(seconds=SETTINGS['SAT_ORBIT_DURATION'])
    cloudproducts.sec1970_start = calendar.timegm(date_time_start.timetuple())
    cloudproducts.sec1970_end = calendar.timegm(date_time_end.timetuple())
    frac_hour = patmosx_hdf.select('scan_line_time').get().astype(np.float)
    if np.ma.is_masked(frac_hour):
        frac_hour = frac_hour.data
    seconds = frac_hour*60*60.0
    time_s = patmosx_hdf.attributes()['time_coverage_start']
    dt_obj = datetime.strptime(time_s, "%Y-%m-%dT%H:%M:%SZ")
    time_sec_1970 =  calendar.timegm(dt_obj.timetuple())
    cloudproducts.time = seconds +  time_sec_1970
    do_some_geo_obj_logging(cloudproducts)
    return cloudproducts


if __name__ == "__main__":
    pass
