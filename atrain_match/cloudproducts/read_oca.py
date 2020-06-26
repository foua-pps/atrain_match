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
"""Use this module to read OCA cloudproducts."""

import atrain_match.config as config
from atrain_match.utils.runutils import do_some_geo_obj_logging
from atrain_match.cloudproducts.read_pps import (
    AllImagerData, AuxiliaryObj,
    CtypeObj, CtthObj, CmaObj, CppObj,
    create_imager_time,
    ImagerAngObj)
import time
import os
import h5netcdf
import netCDF4
import numpy as np
import calendar
from datetime import datetime
from datetime import timedelta
import logging
logger = logging.getLogger(__name__)

ATRAIN_MATCH_NODATA = config.NODATA
# from atrain_match.utils.get_flag_info import get_oca_ct_flag, get_day_night_twilight_info_oca
OCA_READ_EXTRA = ['ctp2', 'ctt2', 'moca_model_final', 'err_log_cot', 'err_ctp', 'err_cre', 'err_log_cot2', 'err_ctp2', 'err_cre2']

def get_satid_datetime_orbit_from_fname_oca(imager_filename):
    # Get satellite name, time, and orbit number from imager_file
    # W_de-airbusDS-friedrichshafen, SAT, SGA1-VII-02-OCA_C_EUM_20190811211148_G_D_20070912114935_20070912115438_T_X____.nc
    sl_ = os.path.basename(imager_filename).split('_')
    if 'MYD03' in os.path.basename(imager_filename):
        # MYD03.A2008054.0955_OCA_20080223_095500.nc
        date_time = datetime.strptime(sl_[2] + sl_[3], '%Y%m%d%H%M%S.nc')
        sat_id = 'eos2' #
    else:
        #W_xx-eumetsat-darmstadt,SAT,SGA1-VII-02-CTP_C_EUM_20191028185827_G_D_20070912115940_20070912120443_T_N____.nc
        # W_de-airbusDS-friedrichshafen, SAT, SGA1-VII-02-OCA_C_EUM_20190811211148_G_D_20070912114935_20070912115438_T_X____.nc  
        date_time = datetime.strptime(sl_[6], '%Y%m%d%H%M%S')
        sat_id = sl_[0].split('-')[0].lower()
    # date_time = cross.time
    # date_time_end = cross.time + timedelta(seconds=SETTINGS['SAT_ORBIT_DURATION'])


    values = {"satellite": sat_id,
              "date_time": date_time,
              "orbit": "99999",
              "date": date_time.strftime("%Y%m%d"),
              "year": date_time.year,
              "month": "%02d" % (date_time.month),
              "time": date_time.strftime("%H%M"),
              "extrai": "",  # asc_or_des
              # "basename":sat_id + "_" + date_time.strftime("%Y%m%d_%H%M_99999"), # "20080613002200-ESACCI",
              "imagerfilename": imager_filename}
    values['basename'] = values["satellite"] + "_" + \
        values["date"] + "_" + values["time"] + "_" + values["orbit"] + "_" + ""
    return values


def oca_read_all(filename):
    if 'MYD' in os.path.basename(filename):
        cloudproducts = oca_read_all_nc_modis(filename)
        logger.info("No timeinfo in OCA-MODIS file calculate from filename + 5min")
        time_info = get_satid_datetime_orbit_from_fname_oca(filename)
        cloudproducts.sec1970_start = calendar.timegm(time_info["date_time"].timetuple())
        cloudproducts.sec1970_end = cloudproducts.sec1970_start + 60*5
        cloudproducts = create_imager_time(cloudproducts, values={})
        cloudproducts.instrument = 'modis'
        do_some_geo_obj_logging(cloudproducts)
        return(cloudproducts)
    return oca_read_all_nc(filename)


def oca_read_all_nc(filename):
    """Read geolocation, angles info, ctth, and cma."""
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
    # cloudproducts.cpp = read_oca_cpp(oca_nc)
    logger.info("Not reading surface temperature")
    logger.info("Not reading channel data")
    cloudproducts.longitude[cloudproducts.longitude > 180] = cloudproducts.longitude[cloudproducts.longitude > 180]-360
    return cloudproducts


def scale_oca_var(oca_var):
    likely_unscaled = oca_var[:, :].data
    # Variables are scaled if inside valid_min/valid_max)
    # gain = getattr(oca_var, 'scale_factor')
    # intercept = getattr(oca_var, 'add_offset')
    nodata = getattr(oca_var, 'missing_value')
    # out = likely_unscaled.astype(np.float)*gain + intercept
    out = likely_unscaled
    out[likely_unscaled == nodata] = ATRAIN_MATCH_NODATA
    return out.astype(np.float), nodata


def read_oca_ctype_cmask_ctth(oca_nc):
    """Read ctth pressure file."""
    ctth = CtthObj()
    ctype = CtypeObj()
    cma = CmaObj()
    ctth.pressure, ctth.p_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctp'])
    ctth.temperature, ctth.t_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctt'])
    cma.cma_ext = np.where(ctth.pressure > 0, 1, 0)
    ctype.phaseflag = None
    ctype.ct_conditions = None
    return ctype, cma, ctth


def read_oca_cpp(oca_nc):
    """Read ctype and flag info from filename."""
    cpp = CppObj()
    cpp.lwp, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['lwp'])
    cpp.phase, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['phase'])
    return cpp


def read_oca_secondlayer_info(oca_nc):
    """Read ctth pressure file, for second layer."""
    ctth2 = CtthObj()
    ctth2.pressure, ctth2.p_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctp'])
    ctth2.temperature, ctth2.t_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['ctt'])
    aux_dict = {"CTTH2": ctth2}
    aux_obj = AuxiliaryObj(aux_dict)
    return aux_obj


def read_oca_angobj(oca_nc):
    """Read angles info from filename."""
    angle_obj = ImagerAngObj()
    angle_obj.satz.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['observation_zenith'])
    angle_obj.sunz.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['solar_zenith'])
    angle_obj.satazimuth.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['observation_azimuth'])
    angle_obj.sunazimuth.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['solar_azimuth'])
    angle_obj.azidiff.data = None
    return angle_obj

# filling on, default _FillValue of 4294967295 used


def read_oca_geoobj(oca_nc, filename):
    """Read geolocation and time info from filename."""
    cloudproducts = AllImagerData()
    # import pdb;pdb.set_trace()
    cloudproducts.longitude, cloudproducts.nodata = scale_oca_var(oca_nc['data']['measurement_data']['longitude'])
    cloudproducts.latitude, cloudproducts.nodata = scale_oca_var(oca_nc['data']['measurement_data']['latitude'])
    cloudproducts.num_of_lines = cloudproducts.longitude.shape[0]
    cloudproducts.nodata = -999
    # sensing_start_time_utc and sensing_end_time_utc,
    stime = getattr(oca_nc, "sensing_start_time_utc")
    etime = getattr(oca_nc, "sensing_end_time_utc")
    cloudproducts.sec1970_start = calendar.timegm(
        time.strptime(stime, '%Y%m%d%H%M%S.%f'))
    cloudproducts.sec1970_end = calendar.timegm(
        time.strptime(etime, '%Y%m%d%H%M%S.%f'))
    cloudproducts = create_imager_time(cloudproducts, values={})
    do_some_geo_obj_logging(cloudproducts)
    return cloudproducts



# OCA MODIS completely different format!!!!!!!!!!!!!!!!!!!!!11


def oca_read_all_nc_modis(filename):
    """Read geolocation, angles info, ctth, and cma."""
    oca_nc = h5netcdf.File(filename, 'r')
    my_dir = os.path.dirname(filename)
    my_file = os.path.basename(filename).replace('OCA', 'CLM')
    oca_nc_clm = h5netcdf.File(os.path.join(my_dir,my_file))
    logger.info("Opening file %s", filename)
    logger.info("Reading longitude, latitude and time ...")
    cloudproducts = read_oca_geoobj_modis(oca_nc, filename)
    logger.info("Reading angles ...")
    cloudproducts.imager_angles = read_oca_angobj_modis(oca_nc)
    logger.info("Reading cloud pressure ...")
    # , angle_obj)
    cma = read_oca_cmask_modis(oca_nc_clm)
    ctype, cma, ctth = read_oca_ctype_cmask_ctth_modis(oca_nc)
    cloudproducts.cma = cma
    cloudproducts.ctth = ctth
    cloudproducts.ctype = ctype
    cloudproducts.aux = read_oca_secondlayer_etc_info_modis(oca_nc)
    logger.info("Not reading cloud microphysical properties")
    # cloudproducts.cpp = read_oca_cpp(oca_nc)
    logger.info("Not reading surface temperature")
    logger.info("Not reading channel data")
    cloudproducts.longitude[cloudproducts.longitude > 180] = cloudproducts.longitude[cloudproducts.longitude > 180]-360
    return cloudproducts


def scale_oca_var_modis(oca_var):
    
    likely_unscaled = oca_var[:, :]
    nodata = -999.0
    out = likely_unscaled
    out[likely_unscaled == nodata] = ATRAIN_MATCH_NODATA
    return out.astype(np.float), nodata


def read_oca_ctype_cmask_ctth_modis(oca_nc):
    """Read ctth pressure file."""
    ctth = CtthObj()
    ctype = CtypeObj()
    cma = CmaObj()
    ctth.pressure, ctth.p_nodata = scale_oca_var_modis(
        oca_nc['ctp'])
    ctth.temperature, ctth.t_nodata = scale_oca_var_modis(
        oca_nc['ctt'])
    
    cma.cma_ext = np.where(oca_nc['moca_model_final'][:] > 0, 1, 0)
    ctype.phaseflag = None
    ctype.ct_conditions = None
    return ctype, cma, ctth

def read_oca_cmask_modis(oca_nc):
    """Read ctth pressure file."""
    cma = CmaObj()
    cma.cma_ext = np.where(oca_nc['flag_cm'][:] > 1, 1, 0)
    cma.cma_bin = np.int64(0*cma.cma_ext.copy())
    cma.cma_bin[oca_nc['flag_cm'][:] == 3] = 1.0
    cma.cma_bin[oca_nc['flag_cm'][:] == 4] = 1.0
    return cma

def read_oca_cpp_modis(oca_nc):
    """Read ctype and flag info from filename."""
    cpp = CppObj()
    cpp.lwp, nodata = scale_oca_var_modis(
        oca_nc['data']['measurement_data']['lwp'])
    cpp.phase, nodata = scale_oca_var_modis(
        oca_nc['moca_model_final'])
    return cpp

def read_oca_secondlayer_etc_info_modis(oca_nc):
    """Read ctth pressure file, for second layer."""
    aux_dict = {}
    for param in OCA_READ_EXTRA:
        aux_dict[param], _ = scale_oca_var_modis(
            oca_nc[param])
    aux_obj = AuxiliaryObj(aux_dict)
    return aux_obj


def read_oca_angobj_modis(oca_nc):
    """Read angles info from filename."""
    angle_obj = ImagerAngObj()
    angle_obj.satz.data, _ = scale_oca_var_modis(
        oca_nc['satzen'])
    angle_obj.sunz.data, _ = scale_oca_var_modis(
        oca_nc['sunzen'])
    angle_obj.satazimuth.data = None
    angle_obj.sunazimuth.data = None
    angle_obj.azidiff.data, _ = scale_oca_var_modis(
        oca_nc['sun_sat_azi'])
    return angle_obj

def read_oca_geoobj_modis(oca_nc, filename):
    """Read geolocation and time info from filename."""
    cloudproducts = AllImagerData()
    # import pdb;pdb.set_trace()
    cloudproducts.longitude = oca_nc['longitude'][:]
    cloudproducts.latitude = oca_nc['latitude'][:]
    cloudproducts.num_of_lines = cloudproducts.longitude.shape[0]
    cloudproducts.nodata = -999
    return cloudproducts


if __name__ == "__main__":
    pass
