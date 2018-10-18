import numpy as np
import os
import netCDF4
import h5py
import logging
import time
import calendar
from datetime import datetime
logger = logging.getLogger(__name__)
from config import NODATA
from libs.extract_imager_along_track import get_data_from_array
ATRAIN_MATCH_NODATA = NODATA
class MOD06Obj:
    #skeleton container for v2014 cloudtype
    def __init__(self):
        self.height = None
        self.temperature = None
        self.pressure = None
        self.cloud_emissivity = None
        self.cloud_phase = None
        self.lwp = None

def add_modis_06(ca_matchup, pps_imager_file, options):

    mfile = find_modis_lvl2_file_from_pps(pps_imager_file, options)
    modis_06 = read_modis_h5(mfile)
    ca_matchup_truth_sat = getattr(ca_matchup, ca_matchup.truth_sat)
    row_matched = ca_matchup_truth_sat.imager_linnum
    col_matched = ca_matchup_truth_sat.imager_pixnum

    index = {'row': row_matched,
             'col': col_matched}
    index_5km = {'row': np.floor(row_matched/5).astype(np.int),
                 'col': np.floor(col_matched/5).astype(np.int)}

    ca_matchup.modis.height = get_data_from_array(modis_06.height, index)
    ca_matchup.modis.temperature = get_data_from_array(modis_06.temperature, index)
    ca_matchup.modis.pressure = get_data_from_array( modis_06.pressure, index)
    ca_matchup.modis.lwp = get_data_from_array( modis_06.lwp, index)
    ca_matchup.modis.cloud_emissivity = get_data_from_array( modis_06.cloud_emissivity, index)
    ca_matchup.modis.cloud_phase = get_data_from_array( modis_06.cloud_phase, index)
    ca_matchup.modis.latitude_5km = get_data_from_array( modis_06.latitude, index_5km)
    ca_matchup.modis.longitude_5km = get_data_from_array( modis_06.longitude, index_5km)
    return ca_matchup

def find_modis_lvl2_file_from_pps(pps_imager_file, options):
    from libs.truth_imager_match import get_pps_file
    from imager_cloud_products.read_cloudproducts_and_nwp_pps import get_satid_datetime_orbit_from_fname_pps
    values = get_satid_datetime_orbit_from_fname_pps(pps_imager_file)
    modis_06_filename = get_pps_file(pps_imager_file, options, values, 'modis_06_file', 'modis_06_dir')
    logger.info("MODIS-C6 file:  %s", os.path.basename(modis_06_filename))
    return modis_06_filename
    
def read_modis_h5(filename):
    h5file = h5py.File(filename, 'r')
    modis_06 = MOD06Obj()
    modis_06.height = h5file['mod06']['Data Fields']['cloud_top_height_1km'].value
    modis_06.temperature = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].value
    modis_06.pressure = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].value
    modis_06.lwp = h5file['mod06']['Data Fields']['Cloud_Water_Path'].value


    h_gain = h5file['mod06']['Data Fields']['cloud_top_height_1km'].attrs['scale_factor']
    h_intercept = h5file['mod06']['Data Fields']['cloud_top_height_1km'].attrs['add_offset']
    t_gain = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].attrs['scale_factor']
    t_intercept = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].attrs['add_offset']
    p_gain = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].attrs['scale_factor']
    p_intercept = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].attrs['add_offset']
    h_nodata = h5file['mod06']['Data Fields']['cloud_top_height_1km'].attrs['_FillValue']
    t_nodata = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].attrs['_FillValue']
    p_nodata = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].attrs['_FillValue']

    lwp_gain = h5file['mod06']['Data Fields']['Cloud_Water_Path'].attrs['scale_factor']
    lwp_intercept = h5file['mod06']['Data Fields']['Cloud_Water_Path'].attrs['add_offset']
    lwp_nodata = h5file['mod06']['Data Fields']['Cloud_Water_Path'].attrs['_FillValue']

    hmask = modis_06.height == h_nodata
    tmask = modis_06.temperature == t_nodata
    pmask = modis_06.pressure == p_nodata
    lwpmask = modis_06.lwp == lwp_nodata
    modis_06.height = modis_06.height.astype(np.float)
    modis_06.pressure = modis_06.pressure.astype(np.float)
    modis_06.temperature = modis_06.temperature.astype(np.float)
    modis_06.lwp = modis_06.lwp.astype(np.float)
    modis_06.height[~hmask] = modis_06.height[~hmask]*h_gain + h_intercept
    modis_06.pressure[~pmask] = modis_06.pressure[~pmask]*p_gain + p_intercept
    modis_06.temperature[~tmask] = modis_06.temperature[~tmask]*t_gain - t_intercept*t_gain
    modis_06.lwp[~lwpmask] = modis_06.lwp[~lwpmask]*lwp_gain + lwp_intercept
    modis_06.height[hmask] = ATRAIN_MATCH_NODATA
    modis_06.pressure[pmask] = ATRAIN_MATCH_NODATA    
    modis_06.temperature[tmask] = ATRAIN_MATCH_NODATA 
    modis_06.lwp[lwpmask] = ATRAIN_MATCH_NODATA 
    modis_06.cloud_emissivity = h5file['mod06']['Data Fields']['cloud_emissivity_1km'].value.astype(np.float)
    modis_06.cloud_phase = h5file['mod06']['Data Fields']['Cloud_Phase_Infrared_1km'].value.astype(np.int)
    modis_06.latitude = h5file['mod06']['Geolocation Fields']['Latitude'].value.astype(np.float)
    modis_06.longitude = h5file['mod06']['Geolocation Fields']['Longitude'].value.astype(np.float)
    h5file.close()
    return modis_06
