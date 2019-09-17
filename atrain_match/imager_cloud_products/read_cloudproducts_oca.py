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
  Use this module to read OCA cloudproducts
  2019 SMHI, N.Hakansson 
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

from atrain_match.imager_cloud_products.read_cloudproducts_and_nwp_pps import (
    AllImagerData, AuxiliaryObj, 
    CtypeObj, CtthObj, CmaObj, CppObj,
    create_imager_time,
    ImagerAngObj)
from atrain_match.utils.runutils import do_some_geo_obj_logging
import atrain_match.config
ATRAIN_MATCH_NODATA = config.NODATA
#from atrain_match.utils.get_flag_info import get_oca_ct_flag, get_day_night_twilight_info_oca

def get_satid_datetime_orbit_from_fname_oca(imager_filename, SETTINGS, cross):
    # Get satellite name, time, and orbit number from imager_file
    #W_de-airbusDS-friedrichshafen,SAT,SGA1-VII-02-OCA_C_EUM_20190811211148_G_D_20070912114935_20070912115438_T_X____.nc
    sl_ = os.path.basename(imager_filename).split('_')
    date_time = datetime.strptime(sl_[6], '%Y%m%d%H%M%S')
    #date_time = cross.time
    #date_time_end = cross.time + timedelta(seconds=SETTINGS['SAT_ORBIT_DURATION'])

    sat_id = sl_[0].split('-')[0].lower()
    values = {"satellite": sat_id,
              "date_time": date_time,
              "orbit": "99999",
              "date": date_time.strftime("%Y%m%d"),
              "year": date_time.year,
              "month": "%02d" % (date_time.month),
              "time": date_time.strftime("%H%M"),
              "extrai": "", #asc_or_des
              #"basename":sat_id + "_" + date_time.strftime("%Y%m%d_%H%M_99999"),#"20080613002200-ESACCI",
              "imagerfilename": imager_filename}
    values['basename'] = values["satellite"] + "_" + \
        values["date"] + "_" + values["time"] + "_" + values["orbit"] + "_" + ""
    return values

def oca_read_all(filename):
    return oca_read_all_nc(filename)

def oca_read_all_nc(filename):
    """Read geolocation, angles info, ctth, and cma
    """
    oca_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    logger.info("Opening file %s", filename)
    logger.info("Reading longitude, latitude and time ...")
    cloudproducts = read_oca_geoobj(oca_nc, filename)
    logger.info("Reading angles ...")
    cloudproducts.imager_angles = read_oca_angobj(oca_nc)
    logger.info("Reading cloud pressure ...")
    # , angle_obj)
    ctype, cma, ctth = read_oca_ctype_cmask_ctth(oca_nc)
    cloudproducts.cma = cma
    cloudproducts.ctth = ctth
    cloudproducts.ctype = ctype
    cloudproducts.aux = read_oca_secondlayer_info(oca_nc)
    logger.info("Not reading cloud microphysical properties")
    #cloudproducts.cpp = read_oca_cpp(oca_nc)
    logger.info("Not reading surface temperature")
    logger.info("Not reading channel data")
    cloudproducts.longitude[cloudproducts.longitude>180] = cloudproducts.longitude[cloudproducts.longitude>180]-360
    return cloudproducts

def scale_oca_var(oca_var):
    likely_unscaled = oca_var[:,:].data
    #Variables are scaled if inside valid_min/valid_max)
    #gain = getattr(oca_var,'scale_factor')
    #intercept = getattr(oca_var,'add_offset')
    nodata = getattr(oca_var,'missing_value')
    #out = likely_unscaled.astype(np.float)*gain + intercept
    out = likely_unscaled
    out[likely_unscaled == nodata] = ATRAIN_MATCH_NODATA
    return out.astype(np.float), nodata

def read_oca_ctype_cmask_ctth(oca_nc):
    """Read ctth pressure file
    """
    ctth = CtthObj()
    ctype = CtypeObj()
    cma = CmaObj()
    ctth.pressure, ctth.p_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctp'])
    ctth.temperature, ctth.t_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctt'])
    cma.cma_ext = np.where(ctth.pressure>0, 1, 0)
    ctype.phaseflag = None
    ctype.ct_conditions = None
    return ctype, cma, ctth

def read_oca_cpp(oca_nc):
    """Read ctype and flag info from filename
    """
    cpp = CppObj()
    cpp.lwp, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['lwp'])
    cpp.phase, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['phase'])
    return cpp

def read_oca_secondlayer_info(oca_nc):
    """Read ctth pressure file, for second layer
    """
    ctth2 = CtthObj()
    ctth2.pressure, ctth2.p_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctp'])
    ctth2.temperature, ctth2.t_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctt'])
    aux_dict = {"CTTH": ctth2}
    aux_obj = AuxiliaryObj(aux_dict)
    return aux_obj


def read_oca_angobj(oca_nc):
    """Read angles info from filename
    """
    angle_obj = ImagerAngObj()
    angle_obj.satz.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['observation_zenith'])
    angle_obj.sunz.data , nodata= scale_oca_var(
        oca_nc['data']['measurement_data']['solar_zenith'])
    angle_obj.satazimuth.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['observation_azimuth'])
    angle_obj.sunazimuth.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['solar_azimuth'])
    angle_obj.azidiff.data = None
    return angle_obj

#filling on, default _FillValue of 4294967295 used


def read_oca_geoobj(oca_nc, filename):
    """Read geolocation and time info from filename
    """
    cloudproducts = AllImagerData()
    #import pdb;pdb.set_trace()
    cloudproducts.longitude, cloudproducts.nodata = scale_oca_var(oca_nc['data']['measurement_data']['longitude'])
    cloudproducts.latitude, cloudproducts.nodata = scale_oca_var(oca_nc['data']['measurement_data']['latitude'])
    cloudproducts.num_of_lines = cloudproducts.longitude.shape[0]
    cloudproducts.nodata=-999
    #sensing_start_time_utc and sensing_end_time_utc, 
    stime = getattr(oca_nc, "sensing_start_time_utc")
    etime = getattr(oca_nc, "sensing_end_time_utc")
    cloudproducts.sec1970_start = calendar.timegm(
        time.strptime(stime, '%Y%m%d%H%M%S.%f')) 
    cloudproducts.sec1970_end = calendar.timegm(
        time.strptime(etime, '%Y%m%d%H%M%S.%f')) 
    cloudproducts = create_imager_time(cloudproducts, values={})
    do_some_geo_obj_logging(cloudproducts)
    return cloudproducts



if __name__ == "__main__":
    pass
