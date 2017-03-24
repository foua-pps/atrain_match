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
ATRAIN_MATCH_NODATA = NODATA
class MOD06Obj:
    #skeleton container for v2014 cloudtype
    def __init__(self):
        self.height = None
        self.temperature = None
        self.pressure = None
        self.cloud_emissivity = None

def add_modis_06(ca_matchup, pps_imager_file, options):

    mfile = find_modis_lvl2_file_from_pps(pps_imager_file, options)
    modis_06 = read_modis_h5(mfile)
    row_matched = ca_matchup.calipso.avhrr_linnum
    col_matched = ca_matchup.calipso.avhrr_pixnum
    index = zip(row_matched,col_matched)
    index_5km = zip(np.floor(row_matched/5).astype(np.int),np.floor(col_matched/5).astype(np.int))
    npix = len(row_matched)
    #import time
    #tic= time.time()
    ca_matchup.modis.height = np.array([ modis_06.height[rid, cid] for rid, cid in index])
    #print time.time()-tic; tic=time.time()
    ca_matchup.modis.temperature = np.array([ modis_06.temperature[rid, cid] for rid, cid in index])
    #print time.time()-tic; tic=time.time()
    ca_matchup.modis.pressure = np.array([ modis_06.pressure[rid, cid] for rid, cid in index])
    ca_matchup.modis.cloud_emissivity = np.array([ modis_06.cloud_emissivity[rid, cid] for rid, cid in index])
    ca_matchup.modis.latitude_5km = np.array([ modis_06.latitude[rid, cid] for rid, cid in index_5km])
    ca_matchup.modis.longitude_5km = np.array([ modis_06.longitude[rid, cid] for rid, cid in index_5km])
    return ca_matchup




def find_modis_lvl2_file_from_pps(pps_imager_file, options):
    from glob import glob
    from cloudsat_calipso_avhrr_match import get_pps_file, insert_info_in_filename_or_path
    from read_cloudproducts_and_nwp_pps import get_satid_datetime_orbit_from_fname_pps
    values = get_satid_datetime_orbit_from_fname_pps(pps_imager_file)
    date_time = values["date_time"]
    modis_06_filename = insert_info_in_filename_or_path(options['modis_06_file'], 
                                                     values, datetime_obj=date_time)
    path = insert_info_in_filename_or_path(options['modis_06_dir'], 
                                           values, datetime_obj=date_time)
    file_names = glob(os.path.join(path, modis_06_filename ))[0]
    print path, modis_06_filename
    print file_names                  
    modis_06_filename = get_pps_file(pps_imager_file, options, values, 'modis_06_file', 'modis_06_dir')
    print modis_06_filename
    return modis_06_filename
    
def read_modis_h5(filename):
    h5file = h5py.File(filename, 'r')
    modis_06 = MOD06Obj()
    modis_06.height = h5file['mod06']['Data Fields']['cloud_top_height_1km'].value.astype(np.float)
    modis_06.temperature = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].value.astype(np.float)
    modis_06.pressure = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].value.astype(np.float)

    #Currently unpacked arrays later in calipso.py
    #TODO: move this here also for h5!
    modis_06.h_gain = h5file['mod06']['Data Fields']['cloud_top_height_1km'].attrs['scale_factor']
    modis_06.h_intercept = h5file['mod06']['Data Fields']['cloud_top_height_1km'].attrs['add_offset']
    modis_06.t_gain = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].attrs['scale_factor']
    modis_06.t_intercept = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].attrs['add_offset']
    modis_06.p_gain = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].attrs['scale_factor']
    modis_06.p_intercept = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].attrs['add_offset']
    modis_06.h_nodata = h5file['mod06']['Data Fields']['cloud_top_height_1km'].attrs['_FillValue']
    modis_06.t_nodata = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'].attrs['_FillValue']
    modis_06.p_nodata = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'].attrs['_FillValue']
    hmask = modis_06.height == modis_06.h_nodata
    tmask = modis_06.temperature == modis_06.t_nodata
    pmask = modis_06.pressure == modis_06.p_nodata
    modis_06.height[~hmask] = modis_06.height[~hmask]*modis_06.h_gain + modis_06.h_intercept
    modis_06.pressure[~pmask] = modis_06.pressure[~pmask]*modis_06.p_gain + modis_06.p_intercept
    modis_06.temperature[~tmask] = modis_06.temperature[~tmask]*modis_06.t_gain - modis_06.t_intercept*modis_06.t_gain
    modis_06.height[hmask] = ATRAIN_MATCH_NODATA
    modis_06.pressure[pmask] = ATRAIN_MATCH_NODATA    
    modis_06.temperature[tmask] = ATRAIN_MATCH_NODATA 
    modis_06.cloud_emissivity = h5file['mod06']['Data Fields']['cloud_emissivity_1km'].value.astype(np.float)
    modis_06.latitude = h5file['mod06']['Geolocation Fields']['Latitude'].value.astype(np.float)
    modis_06.longitude = h5file['mod06']['Geolocation Fields']['Longitude'].value.astype(np.float)
    return modis_06
