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
OCA_READ_EXTRA_ONLY_OCA = ['ctp2', 'ctt2', 'moca_model_final', 'err_log_cot', 'err_ctp', 'err_cre', 'err_log_cot2', 'err_ctp2', 'err_cre2']

OCA_READ_EXTRA = OCA_READ_EXTRA_ONLY_OCA + ['flag_cm']

def get_satid_datetime_orbit_from_fname_oca(imager_filename):
    # Get satellite name, time, and orbit number from imager_file
    # W_de-airbusDS-friedrichshafen, SAT, SGA1-VII-02-OCA_C_EUM_20190811211148_G_D_20070912114935_20070912115438_T_X____.nc
    sl_ = os.path.basename(imager_filename).split('_')
    if "OCA_SEV_MET" in os.path.basename(imager_filename):
        # MYD03.A2008054.0955_OCA_20080223_095500.nc
        date_time = datetime.strptime(sl_[-1], '%Y%m%d%H%M.nc')
        sat_id = 'meteosat9' #
    elif 'MYD03' in os.path.basename(imager_filename):
        # MYD03.A2008054.0955_OCA_20080223_095500.nc
        date_time = datetime.strptime(sl_[2] + sl_[3], '%Y%m%d%H%M%S.nc')
        sat_id = 'eos2' #
    elif 'SGA1-VII' in os.path.basename(imager_filename):
        # MYD03.A2008054.0955_OCA_20080223_095500.nc
        date_time = datetime.strptime(sl_[7], '%Y%m%d%H%M%S')
        sat_id = 'metopsga1'
    elif 'MET09+SEVIRI' in os.path.basename(imager_filename):
        date_time = datetime.strptime(sl_[4], '%Y%m%d%H%M%S')
        sat_id = 'meteosat9'
    else:
        # W_xx-eumetsat-darmstadt,SAT,SGA1-VII-02-CTP_C_EUM_20191028185827_G_D_20070912115940_20070912120443_T_N____.nc
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


def oca_read_all(filename, extra_files):
    if 'MYD' in os.path.basename(filename):
        cloudproducts, aux_dict = oca_read_all_nc_modis(filename)
        logger.info("No timeinfo in OCA-MODIS file calculate from filename + 5min")
        seconds = 5*60
        cloudproducts.instrument = 'modis'
    elif 'MET09+SEVIRI' in os.path.basename(filename):
        cloudproducts, aux_dict = oca_read_all_nc_cdr(filename)
        logger.info("Calculate time from filename + 12min")
        seconds = 12*60
    elif "OCA_SEV_MET" in  os.path.basename(filename):
        cloudproducts, aux_dict = oca_read_all_nc_modis(filename)
        logger.info("No timeinfo in OCA-SEVIRI file calculate from filename +15min")
        seconds = 15*60
        cloudproducts.instrument = 'SEVIRI-OCA'
    else:
        cloudproducts, aux_dict = oca_read_all_nc(filename, extra_files)

    if  cloudproducts.sec1970_start is None:   
        time_info = get_satid_datetime_orbit_from_fname_oca(filename)
        cloudproducts.sec1970_start = calendar.timegm(time_info["date_time"].timetuple())
        cloudproducts.sec1970_end = cloudproducts.sec1970_start + seconds
        cloudproducts = create_imager_time(cloudproducts, values={})
        do_some_geo_obj_logging(cloudproducts, extra_files)

    aux_dict = add_claas3(cloudproducts, extra_files, aux_dict)
    cloudproducts.aux = AuxiliaryObj(aux_dict)
    return cloudproducts

def add_claas3(cloudproducts, extra_files, aux_dict):
    """Add claas3 data."""
    if "cma_claas3" in extra_files:
        from atrain_match.utils.get_flag_info import get_day_night_twilight_info_pps2014
        claas3_nc = netCDF4.Dataset(extra_files["cma_claas3"], 'r', format='NETCDF4')
        aux_dict["claas3_cma"] = np.squeeze(claas3_nc['cma'][:, -1::-1, -1::-1])
        aux_dict["claas3_cma_prob"] = np.squeeze(claas3_nc['cma_prob'][:, -1::-1, -1::-1])
        aux_dict["claas3_conditions"] = np.squeeze(claas3_nc['conditions'][:, -1::-1, -1::-1])
        claas3_nc.close()
        claas3_nc = netCDF4.Dataset(extra_files["cth_claas3"], 'r', format='NETCDF4')
        aux_dict["claas3_cth"] = np.squeeze(claas3_nc['cth'][:, -1::-1, -1::-1])
        aux_dict["claas3_ctp"] = np.squeeze(claas3_nc['ctp'][:, -1::-1, -1::-1])
        aux_dict["claas3_ctt"] = np.squeeze(claas3_nc['ctt'][:, -1::-1, -1::-1])
        aux_dict["claas3_cth_unc"] = np.squeeze(claas3_nc['cth_unc'][:, -1::-1, -1::-1])
        aux_dict["claas3_ctp_unc"] = np.squeeze(claas3_nc['ctp_unc'][:, -1::-1, -1::-1])
        aux_dict["claas3_ctt_unc"] = np.squeeze(claas3_nc['ctt_unc'][:, -1::-1, -1::-1])
        #aux_dict["claas3_conditions"] = np.squeeze(claas3_nc['conditions'][:, -1::-1, -1::-1])
        claas3_nc.close()
        claas3_nc = netCDF4.Dataset(extra_files["cpp_claas3"], 'r', format='NETCDF4')
        aux_dict["claas3_cph"] = np.squeeze(claas3_nc['cph'][:, -1::-1, -1::-1])
        aux_dict["claas3_cwp"] = np.squeeze(claas3_nc['cwp'][:, -1::-1, -1::-1])
        claas3_nc.close()
        no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag = get_day_night_twilight_info_pps2014(
            aux_dict["claas3_conditions"])
        cloudproducts.imager_angles.sunz = np.where(night_flag, 100, 90)
        cloudproducts.imager_angles.sunz = np.where(day_flag, 10,  cloudproducts.imager_angles.sunz)
    for key in aux_dict:
        if "claas3" in key:
            aux_dict[key][aux_dict[key].mask] = ATRAIN_MATCH_NODATA
    return aux_dict


def oca_read_all_nc_cdr(filename):
    """Read geolocation, angles info, ctth, and cma."""
    oca_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    logger.info("Opening file %s", filename)
    logger.info("Reading longitude, latitude and time ...")
    cloudproducts = read_oca_geoobj_modis(oca_nc, filename)
    cloudproducts.imager_angles = ImagerAngObj()
    cloudproducts.imager_angles.sunz = np.where(cloudproducts.longitude > -999999, 0, 100)
    #cloudproducts.imager_angles.satz = np.where(cloudproducts.longitude > -999999, 0, 100)

    logger.info("Reading angles ...")
    # No angles!
    logger.info("Reading cloud pressure ...")
    # , angle_obj)
    ctype, cma, ctth, aux_dict = read_oca_ctype_cmask_ctth_cdr(oca_nc)

    cloudproducts.cma = cma
    cloudproducts.ctth = ctth
    cloudproducts.ctype = ctype
    aux_dict = read_oca_secondlayer_info(oca_nc, aux_dict=aux_dict)
    logger.info("Not reading cloud microphysical properties")
    # cloudproducts.cpp = read_oca_cpp(oca_nc)
    logger.info("Not reading surface temperature")
    logger.info("Not reading channel data")
    cloudproducts.longitude[cloudproducts.longitude > 180] = cloudproducts.longitude[cloudproducts.longitude > 180]-360
    oca_nc.close()
    return cloudproducts, aux_dict


def oca_read_all_nc(filename, extra_files):
    """Read geolocation, angles info, ctth, and cma."""
    oca_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    logger.info("Opening file %s", filename)
    logger.info("Reading longitude, latitude and time normal nc...")
    cloudproducts = read_oca_geoobj(oca_nc, filename)
    logger.info("Reading angles ...")
    cloudproducts.imager_angles = read_oca_angobj(oca_nc)
    logger.info("Reading cloud pressure ...")
    # , angle_obj)
    ctype, cma, ctth = read_oca_ctype_cmask_ctth(oca_nc)
    if "oca_cma" in extra_files:
        cma = read_oca_cmask_epssg(extra_files["oca_cma"])
    cloudproducts.cma = cma
    cloudproducts.ctth = ctth
    cloudproducts.ctype = ctype
    aux_dict = read_oca_secondlayer_info(oca_nc)
    logger.info("Not reading cloud microphysical properties")
    # cloudproducts.cpp = read_oca_cpp(oca_nc)
    logger.info("Not reading surface temperature")
    logger.info("Not reading channel data")
    cloudproducts.longitude[cloudproducts.longitude > 180] = cloudproducts.longitude[cloudproducts.longitude > 180]-360
    return cloudproducts, aux_dict


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
    ctth.height, ctth.h_nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['cth'])
    ctth.height = ctth.height * 1000  # km =>m
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


def read_oca_secondlayer_info(oca_nc, aux_dict={}):
    """Read ctth pressure file, for second layer."""
    try:
        ctp2, ____ = scale_oca_var(
            oca_nc['cloud_top_pressure_lower_layer'])
        aux_dict['oca_cloud_top_pressure_lower_layer'] = ctp2
    except (KeyError, IndexError):
        pass
        
    return aux_dict


def read_oca_angobj(oca_nc):
    """Read angles info from filename."""
    angle_obj = ImagerAngObj()
    angle_obj.satz.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['observation_zenith'])
    angle_obj.sunz.data, nodata = scale_oca_var(
        oca_nc['data']['measurement_data']['solar_zenith'])
    #angle_obj.satazimuth.data, nodata = scale_oca_var(
    #    oca_nc['data']['measurement_data']['observation_azimuth'])
    #angle_obj.sunazimuth.data, nodata = scale_oca_var(
    #    oca_nc['data']['measurement_data']['solar_azimuth'])
    angle_obj.azidiff.data = None

    from geotiepoints.viiinterpolator import tie_points_interpolation
    import xarray as xr
    if angle_obj.satz.data.shape[0] <200:        
        print("interpolating!")
        angle_obj.satz.data, angle_obj.sunz.data = np.array(tie_points_interpolation([xr.DataArray(angle_obj.satz.data), xr.DataArray(angle_obj.sunz.data)], 4,8))


    return angle_obj

# filling on, default _FillValue of 4294967295 used


def read_oca_geoobj(oca_nc, filename):
    """Read geolocation and time info from filename."""
    cloudproducts = AllImagerData()
    # import pdb;pdb.set_trace()
    cloudproducts.longitude, cloudproducts.nodata = scale_oca_var(oca_nc['data']['measurement_data']['longitude'])
    cloudproducts.latitude, cloudproducts.nodata = scale_oca_var(oca_nc['data']['measurement_data']['latitude'])

    from geotiepoints.viiinterpolator import tie_points_geo_interpolation
    import xarray as xr
    if cloudproducts.longitude.shape[0] < 200:
        print("interpolating!", cloudproducts.longitude[0::4,0::3].shape)
        lon_interp, lat_interp = tie_points_geo_interpolation(
            
            xr.DataArray(cloudproducts.longitude), xr.DataArray(cloudproducts.latitude),  4, 8)
        cloudproducts.longitude = np.array(lon_interp)
        cloudproducts.latitude = np.array(lat_interp)
        
    cloudproducts.num_of_lines = cloudproducts.longitude.shape[0]
    cloudproducts.nodata = -999
    # sensing_start_time_utc and sensing_end_time_utc,
    stime = getattr(oca_nc, "sensing_start_time_utc")
    etime = getattr(oca_nc, "sensing_end_time_utc")
    try:
        cloudproducts.sec1970_start = calendar.timegm(
            time.strptime(stime, '%Y%m%d%H%M%S.%f'))
        cloudproducts.sec1970_end = calendar.timegm(
            time.strptime(etime, '%Y%m%d%H%M%S.%f'))
        cloudproducts = create_imager_time(cloudproducts, values={})
        do_some_geo_obj_logging(cloudproducts)
    except:
        cloudproducts.sec1970_start = calendar.timegm(
            time.strptime(stime, '%Y-%m-%d %H:%M:%S.%f'))
        cloudproducts.sec1970_end = calendar.timegm(
            time.strptime(etime, '%Y-%m-%d %H:%M:%S.%f'))
        cloudproducts = create_imager_time(cloudproducts, values={})
        do_some_geo_obj_logging(cloudproducts)

    return cloudproducts



# OCA MODIS completely different format!!!!!!!!!!!!!!!!!!!!!11


def oca_read_all_nc_modis(filename, clm_separate_file=False):
    """Read geolocation, angles info, ctth, and cma."""
    oca_nc = h5netcdf.File(filename, 'r')
    my_dir = os.path.dirname(filename)
    if clm_separate_file:
        my_file = os.path.basename(filename).replace('OCA', 'CLM')
        oca_nc_clm = h5netcdf.File(os.path.join(my_dir, my_file))
        logger.info("Opening file %s", filename)
        cma = read_oca_cmask_modis(oca_nc_clm)
    logger.info("Reading longitude, latitude and time ...")
    cloudproducts = read_oca_geoobj_modis(oca_nc, filename)
    logger.info("Reading angles ...")
    cloudproducts.imager_angles = read_oca_angobj_modis(oca_nc)
    logger.info("Reading cloud pressure ...")
    # , angle_obj)
    ctype, cma_dummy, ctth = read_oca_ctype_cmask_ctth_modis(oca_nc)
    cloudproducts.cma = cma
    cloudproducts.ctth = ctth
    cloudproducts.ctype = ctype
    aux_dict = read_oca_secondlayer_etc_info_modis(oca_nc)
    #? aux_dict['flag_cm'] = read_oca_secondlayer_etc_info_modis(oca_nc_clm, params=['flag_cm'])['flag_cm']
    logger.info("Not reading cloud microphysical properties")
    # cloudproducts.cpp = read_oca_cpp(oca_nc)
    logger.info("Not reading surface temperature")
    logger.info("Not reading channel data")
    cloudproducts.longitude[cloudproducts.longitude > 180] = cloudproducts.longitude[cloudproducts.longitude > 180]-360
    return cloudproducts, aux_dict


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
    if 'ctt' in oca_nc.keys():
        ctth.temperature, ctth.t_nodata = scale_oca_var_modis(
            oca_nc['ctt'])
    
    cma.cma_ext = np.where(oca_nc['moca_model_final'][:] > 0, 1, 0)
    ctype.phaseflag = None
    ctype.ct_conditions = None
    return ctype, cma, ctth

def read_oca_ctype_cmask_ctth_cdr(oca_nc):
    """Read ctth pressure file."""
    ctth = CtthObj()
    ctype = CtypeObj()
    cma = CmaObj()
    aux_dict = {}
    pressure, ___ = scale_oca_var(
        oca_nc['cloud_top_pressure'])
    pressure_error, ___ = scale_oca_var(
        oca_nc['cloud_top_pressure_error_log'])
    ctth_pressure_error = np.power(10, pressure_error) * 0.01 #hPa
    aux_dict["oca_cloud_top_pressure"] = pressure * 0.01
    aux_dict["oca_cloud_top_pressure_error"] = ctth_pressure_error
    
    cma_prob, ___ = scale_oca_var(oca_nc['cloud_probability'])
    cma.cma_ext = np.where(cma_prob > 50, 1, 0)
    aux_dict["oca_cloud_probability"] = cma_prob    
    ctype.phaseflag = None
    ctype.ct_conditions = None
    return ctype, cma, ctth, aux_dict

def read_oca_cmask_epssg(filename):
    """Read cmask file."""
    oca_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    cma = CmaObj()
    try:
        cm_flag = oca_nc['data']['measurement_data']['flag_cm'][:]
        cma.cma_ext = np.where(cm_flag > 2, 1, 0)
        cma.cma_ext[cm_flag>200] = -9
        cma.cma_ext[cm_flag==5] = -9
        cma.cma_bin = np.int64(0*cma.cma_ext.copy())
        cma.cma_bin[cm_flag == 3] = 1.0
        cma.cma_bin[cm_flag == 4] = 1.0
        
    except:
        import pdb;pdb.set_trace()
        cma.cma_ext = np.where(oca_nc['moca_model_final'][:] >= 1, 1, 0)
        cma.cma_bin = np.int64(0*cma.cma_ext.copy())
    oca_nc.close()    
    return cma

def read_oca_cmask_modis(oca_nc):
    """Read cmask file."""
    cma = CmaObj()
    if 'flag_cm' in oca_nc.keys():
        cma.cma_ext = np.where(oca_nc['flag_cm'][:] > 1, 1, 0)
        cma.cma_bin = np.int64(0*cma.cma_ext.copy())
        cma.cma_bin[oca_nc['flag_cm'][:] == 3] = 1.0
        cma.cma_bin[oca_nc['flag_cm'][:] == 4] = 1.0
    else:
        import pdb;pdb.set_trace()
        cma.cma_ext = np.where(oca_nc['moca_model_final'][:] >= 1, 1, 0)
        cma.cma_bin = np.int64(0*cma.cma_ext.copy())
    return cma

def read_oca_cpp_modis(oca_nc):
    """Read ctype and flag info from filename."""
    cpp = CppObj()
    cpp.lwp, nodata = scale_oca_var_modis(
        oca_nc['data']['measurement_data']['lwp'])
    cpp.phase, nodata = scale_oca_var_modis(
        oca_nc['moca_model_final'])
    return cpp

def read_oca_secondlayer_etc_info_modis(oca_nc, params=OCA_READ_EXTRA_ONLY_OCA):
    """Read ctth pressure file, for second layer."""
    aux_dict = {}
    for param in params:
        if param in oca_nc.keys():
            aux_dict[param], _ = scale_oca_var_modis(
                oca_nc[param])
    return aux_dict


def read_oca_angobj_modis(oca_nc):
    """Read angles info from filename."""
    angle_obj = ImagerAngObj()
    angle_obj.satz.data, _ = scale_oca_var_modis(
        oca_nc['satzen'])
    angle_obj.sunz.data, _ = scale_oca_var_modis(
        oca_nc['sunzen'])
    angle_obj.satazimuth.data = None
    angle_obj.sunazimuth.data = None
    if 'sun_sat_azi' in oca_nc.keys():
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
