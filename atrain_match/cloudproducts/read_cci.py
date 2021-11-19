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
  Use this module to read cci cloudproducts
  2013 SMHI, N.Hakansson
"""

from atrain_match.cloudproducts.read_pps import (
    AllImagerData, CppObj,
    CtthObj, CmaObj,
    ImagerAngObj)
from atrain_match.utils.runutils import do_some_geo_obj_logging
import atrain_match.config as config
import os
import netCDF4
import numpy as np
import calendar
import datetime
import logging
logger = logging.getLogger(__name__)
ATRAIN_MATCH_NODATA = config.NODATA


def get_satid_datetime_orbit_from_fname_cci(imager_filename):
    # Get satellite name, time, and orbit number from imager_file
    # imager_file = "20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-IMAGERGAC-NOAA18-fv1.0.nc"
    sl_ = os.path.basename(imager_filename).split('-')
    date_time = datetime.datetime.strptime(sl_[0], '%Y%m%d%H%M%S')

    sat_id = sl_[5].lower()
    values = {"satellite": sat_id,
              "date_time": date_time,
              "orbit": "99999",
              "date": date_time.strftime("%Y%m%d"),
              "year": date_time.year,
              "month": "%02d" % (date_time.month),
              "time": date_time.strftime("%H%M"),
              # "basename":sat_id + "_" + date_time.strftime("%Y%m%d_%H%M_99999"), # "20080613002200-ESACCI",
              "ccifilename": imager_filename,
              "ppsfilename": None}
    values['basename'] = values["satellite"] + "_" + values["date"] + "_" + values["time"] + "_" + values["orbit"]
    return values


def daysafter4713bc_to_sec1970(bcdate_array):
    """Translating days after 4713bc 12:00 to sec since 1970
    """
    # import pps_time_util # @UnresolvedImport
    # first date if several in the file, swath at midnight
    bcdate = np.min(bcdate_array)
    # days since 4713 bc 12:00:00
    ddays = np.floor(bcdate)
    # seconds after 12:00:00
    dseconds = np.floor(24*60*60*(bcdate-ddays))
    mseconds = np.floor(10**6*(24*60*60*(bcdate-ddays)-dseconds))
    # to avoid too large numbers count from 1950-01-01 12:00
    # time datetime.datetime(1950, 1, 1, 12, 0, 0, 0) == julian day 2433283
    ddays = ddays-2433283

    time_delta = datetime.timedelta(days=ddays, seconds=dseconds,
                                    microseconds=mseconds)
    the_time = datetime.datetime(1950, 1, 1, 12, 0, 0, 0) + time_delta
    sec_1970 = calendar.timegm(the_time.timetuple()) + the_time.microsecond/(1.0*10**6)
    sec_1970_array = 24*60*60*(bcdate_array - bcdate) + sec_1970
    return sec_1970_array


def cci_read_all(filename):
    """Read geolocation, angles info, ctth, and cloudtype
    """
    logger.info("Opening file %s", filename)
    cci_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    cloudproducts = read_cci_geoobj(cci_nc)
    logger.debug("Reading ctth ...")
    cloudproducts.ctth = read_cci_ctth(cci_nc)
    logger.debug("Reading angles ...")
    cloudproducts.imager_angles = read_cci_angobj(cci_nc)
    logger.debug("Reading cloud type ...")
    cloudproducts.cma = read_cci_cma(cci_nc)

    logger.debug("Reading longitude, latitude and time ...")
    cloudproducts.instrument = "imager"
    for imager in ["avhrr", "viirs", "modis", "seviri"]:
        if imager in os.path.basename(filename).lower():
            cloudproducts.instrument = imager
    logger.debug("Not reading surface temperature")
    logger.debug("Reading cloud phase")
    cloudproducts.cpp = read_cci_phase(cci_nc)
    logger.debug("Reading LWP")
    cloudproducts.cpp = read_cci_lwp(cci_nc, cloudproducts.cpp)
    logger.debug("Reading LSM")
    cloudproducts.aux = read_cci_lsm(cci_nc)
    logger.debug("Not reading channel data")
    if cci_nc:
        cci_nc.close()
    return cloudproducts
  
  
def read_cci_lsm(cci_nc):
    """ Read land-sea mask from CCI file.
    """
    fractionofland = cci_nc.variables['lsflag'][::]
    aux = AuxiliaryObj({'fractionofland': fractionofland})
    return aux


def read_cci_cma(cci_nc):
    """Read cloudtype and flag info from filename
    """
    cma = CmaObj()
    # cci_nc.variables['cc_total'][:])
    # cci_nc.variables['ls_flag'][:])
    # ctype = CtypeObj()
    # pps 0 cloudfree 3 cloudfree snow 1: cloudy, 2:cloudy
    # cci lsflag 0:sea 1:land
    cma.cma_ext = 0*cci_nc.variables['lsflag'][::]
    cma.cma_ext[cci_nc.variables['cc_total'][::] > 0.5] = 1
    # ctype.landseaflag = cci_nc.variables['lsflag'][::]
    return cma


def read_cci_angobj(cci_nc):
    """Read angles info from filename
    """
    my_angle_obj = ImagerAngObj()
    my_angle_obj.satz.data = cci_nc.variables['satellite_zenith_view_no1'][::]
    my_angle_obj.sunz.data = cci_nc.variables['solar_zenith_view_no1'][::]
    my_angle_obj.azidiff = None  # cci_nc.variables['rel_azimuth_view_no1']??
    return my_angle_obj


def read_cci_phase(cci_nc):
    """Read angles info from filename
    """
    cpp_obj = CppObj()
    data = cci_nc.variables['ann_phase'][::]
    data = np.squeeze(data)
    data = np.where(data < 1, ATRAIN_MATCH_NODATA, data)
    setattr(cpp_obj, 'cpp_phase', data)
    
    data = cci_nc.variables['qcflag'][::]
    data = np.squeeze(data)
    setattr(cpp_obj, 'cpp_quality', data)
    
    data = cci_nc.variables['cot'][::]
    data = np.squeeze(data)
    data = np.where(data < 1, ATRAIN_MATCH_NODATA, data)
    setattr(cpp_obj, 'cpp_cot', data)
    
    return cpp_obj
  
  
  def read_cci_lwp(cci_nc, cpp_obj):
    """ Read LWP from CCI file 
    """
    data = cci_nc.variables['cwp'][::]
    data = np.squeeze(data)
    # split LWP from CWP
    data = np.where(cpp_obj.cpp_phase == 2, ATRAIN_MATCH_NODATA, data)
    data = np.where(data <= 0, ATRAIN_MATCH_NODATA, data)
    setattr(cpp_obj, 'cpp_lwp', data)
    return cpp_obj


def read_cci_geoobj(cci_nc):
    """Read geolocation and time info from filename
    """
    cloudproducts = AllImagerData()
    logger.debug("Min lon: %s, max lon: %d",
                 np.min(cci_nc.variables['lon'][::]),
                 np.max(cci_nc.variables['lon'][::]))
    # cci_nc.variables['lon'].add_offset
    # cloudproducts.longitude = cci_nc.variables['lon'][::]
    cloudproducts.nodata = -999.0
    in_fillvalue = cci_nc.variables['lon']._FillValue
    cloudproducts.longitude = cci_nc.variables['lon'][::]
    cloudproducts.longitude[cci_nc.variables['lon'][::] == in_fillvalue] = cloudproducts.nodata
    cloudproducts.latitude = cci_nc.variables['lon'][::]
    cloudproducts.latitude[cci_nc.variables['lon'][::] == in_fillvalue] = cloudproducts.nodata
    np.where(
        np.logical_and(
            np.greater_equal(cci_nc.variables['lon'][::],
                             cci_nc.variables['lon'].valid_min),
            np.less_equal(cci_nc.variables['lon'][::],
                          cci_nc.variables['lon'].valid_max)),
        cci_nc.variables['lon'][::],
        cloudproducts.nodata)
    # cloudproducts.latitude = cci_nc.variables['lat'][::]
    cloudproducts.latitude = np.where(
        np.logical_and(
            np.greater_equal(cci_nc.variables['lat'][::],
                             cci_nc.variables['lat'].valid_min),
            np.less_equal(cci_nc.variables['lat'][::],
                          cci_nc.variables['lat'].valid_max)),
        cci_nc.variables['lat'][::],
        cloudproducts.nodata)

    #  For pps these are calculated in calipso.py, but better read them
    # from file because time are already available on arrays in the netcdf files
    # dsec = calendar.timegm((1993, 1, 1, 0, 0, 0, 0, 0, 0)) # TAI to UTC
    time_temp = daysafter4713bc_to_sec1970(cci_nc.variables['time'][::])
    cloudproducts.time = time_temp[::]  # [:, 0]  # time_temp[:, 0]

    cloudproducts.sec1970_start = np.min(cloudproducts.time)
    cloudproducts.sec1970_end = np.max(cloudproducts.time)
    do_some_geo_obj_logging(cloudproducts)
    return cloudproducts


def read_cci_ctth(cci_nc):
    """Read cloud top: temperature, height and pressure from filename
    """
    ctth = CtthObj()
    cth = cci_nc.variables['cth'][::]
    if hasattr(cth, 'mask'):
        cth_data = np.where(cth.mask, ATRAIN_MATCH_NODATA, cth.data)
    else:
        cth_data = cth.data  # already scaled!

    logger.info("Setting ctth for non cloudy pixels do nodata ...")
    cth_data[cci_nc.variables['cc_total'][::] < 0.5] = ATRAIN_MATCH_NODATA

    cth_corr = cci_nc.variables['cth_corrected'][::]
    if hasattr(cth_corr, 'mask'):
        cth_data_corr = np.where(cth_corr.mask, ATRAIN_MATCH_NODATA, cth_corr.data)
    else:
        cth_data_corr = cth_corr.data  # already scaled!
    logger.debug("Setting ctth for non cloudy pixels do nodata ...")
    cth_data_corr[cci_nc.variables['cc_total'][::] < 0.5] = ATRAIN_MATCH_NODATA

    ctth.h_gain = 1.0
    ctth.h_intercept = 0.0
    ctth.h_nodata = ATRAIN_MATCH_NODATA
    ctth.height = 1000*cth_data
    ctth.height_corr = 1000*cth_data_corr

    ctt = cci_nc.variables['ctt'][::]
    if hasattr(ctt, 'mask'):
        ctt_data = np.where(ctt.mask, ATRAIN_MATCH_NODATA, ctt.data)
    else:
        ctt_data = ctt.data
    logger.debug("Setting ctth for non cloudy pixels do nodata ...")
    ctt_data[cci_nc.variables['cc_total'][::] < 0.5] = ATRAIN_MATCH_NODATA

    ctth.t_gain = 1.0
    ctth.t_intercept = 0.0
    ctth.t_nodata = ATRAIN_MATCH_NODATA
    ctth.temperature = ctt_data

    ctp = cci_nc.variables['ctp'][::]
    if hasattr(ctp, 'mask'):
        ctp_data = np.where(ctp.mask, ATRAIN_MATCH_NODATA, ctp.data)  # Upscaled
    else:
        ctp_data = ctp.data
    ctth.p_gain = 1.0
    ctth.p_intercept = 0.0
    ctth.p_nodata = ATRAIN_MATCH_NODATA
    logger.debug("Setting ctth for non cloudy pixels do nodata ...")
    ctp_data[cci_nc.variables['cc_total'][::] < 0.5] = ATRAIN_MATCH_NODATA
    ctth.pressure = ctp_data

    return ctth


if __name__ == "__main__":
    my_filename = "20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-IMAGERGAC-NOAA18-fv1.0.nc"
    PPS_OBJECTS = read_cci_ctth(my_filename)
    print(PPS_OBJECTS)
