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
from atrain_match.libs.extract_imager_along_track import get_data_from_array
from atrain_match.config import NODATA
import numpy as np
import os
import h5py
import logging
logger = logging.getLogger(__name__)
ATRAIN_MATCH_NODATA = NODATA


class MOD06Obj:
    # skeleton container for MODIS Level2 data
    def __init__(self):
        self.height = None
        self.temperature = None
        self.pressure = None
        self.cloud_emissivity = None
        self.cloud_phase = None
        self.lwp = None
        self.multilayer = None
        self.optical_depth = None

def add_modis_06(ca_matchup, AM_PATHS, cross):

    mfile = find_modis_lvl2_file(AM_PATHS, cross)
    if mfile.endswith('.h5'):
        modis_06 = read_modis_h5(mfile)
    else:
        modis_06 = read_modis_hdf(mfile)
    ca_matchup_truth_sat = getattr(ca_matchup, ca_matchup.truth_sat)
    row_matched = ca_matchup_truth_sat.imager_linnum
    col_matched = ca_matchup_truth_sat.imager_pixnum

    index = {'row': row_matched,
             'col': col_matched}
    index_5km = {'row': np.floor(row_matched/5).astype(np.int),
                 'col': np.floor(col_matched/5).astype(np.int)}

    ca_matchup.modis_lvl2.height = get_data_from_array(modis_06.height, index)
    ca_matchup.modis_lvl2.temperature = get_data_from_array(modis_06.temperature, index)
    ca_matchup.modis_lvl2.pressure = get_data_from_array(modis_06.pressure, index)
    ca_matchup.modis_lvl2.lwp = get_data_from_array(modis_06.lwp, index)
    ca_matchup.modis_lvl2.cloud_emissivity = get_data_from_array(modis_06.cloud_emissivity, index)
    ca_matchup.modis_lvl2.multilayer = get_data_from_array(modis_06.multilayer, index)
    ca_matchup.modis_lvl2.optical_depth = get_data_from_array(modis_06.optical_depth, index)
    ca_matchup.modis_lvl2.cloud_phase = get_data_from_array(modis_06.cloud_phase, index)
    ca_matchup.modis_lvl2.latitude_5km = get_data_from_array(modis_06.latitude, index_5km)
    ca_matchup.modis_lvl2.longitude_5km = get_data_from_array(modis_06.longitude, index_5km)
    return ca_matchup


def find_modis_lvl2_file(AM_PATHS, cross):
    from atrain_match.libs.truth_imager_match import find_main_cloudproduct_file
    modis_06_filename, _ = find_main_cloudproduct_file(cross,
                                                    AM_PATHS['modis_06_dir'],
                                                    AM_PATHS['modis_06_file'])
    logger.info("MODIS-C6 file:  {:s}".format(os.path.basename(modis_06_filename)))
    return modis_06_filename


def read_modis_h5(filename):
    h5file = h5py.File(filename, 'r')
    modis_06 = MOD06Obj()
    modis_06.height = h5file['mod06']['Data Fields']['cloud_top_height_1km'][:]
    modis_06.temperature = h5file['mod06']['Data Fields']['cloud_top_temperature_1km'][:]
    modis_06.pressure = h5file['mod06']['Data Fields']['cloud_top_pressure_1km'][:]
    modis_06.lwp = h5file['mod06']['Data Fields']['Cloud_Water_Path'][:]
    modis_06.multilayer = h5file['mod06']['Data Fields']['Cloud_Multi_Layer_Flag'][:]
    modis_06.optical_depth = h5file['mod06']['Data Fields']['Cloud_Optical_Thickness'][:]
    modis_06.cloud_emissivity = h5file['mod06']['Data Fields']['cloud_emissivity_1km'][:].astype(np.float)
    modis_06.cloud_phase = h5file['mod06']['Data Fields']['Cloud_Phase_Infrared_1km'][:].astype(np.int)
    modis_06.latitude = h5file['mod06']['Geolocation Fields']['Latitude'][:].astype(np.float)
    modis_06.longitude = h5file['mod06']['Geolocation Fields']['Longitude'][:].astype(np.float)

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

    od_gain = h5file['mod06']['Data Fields']['Cloud_Optical_Thickness'].attrs['scale_factor']
    od_intercept = h5file['mod06']['Data Fields']['Cloud_Optical_Thickness'].attrs['add_offset']
    od_nodata = h5file['mod06']['Data Fields']['Cloud_Optical_Thickness'].attrs['_FillValue']

    hmask = modis_06.height == h_nodata
    tmask = modis_06.temperature == t_nodata
    pmask = modis_06.pressure == p_nodata
    lwpmask = modis_06.lwp == lwp_nodata
    odmask = modis_06.optical_depth == od_nodata

    modis_06.height = modis_06.height.astype(np.float)
    modis_06.pressure = modis_06.pressure.astype(np.float)
    modis_06.temperature = modis_06.temperature.astype(np.float)
    modis_06.lwp = modis_06.lwp.astype(np.float)
    modis_06.optical_depth = modis_06.optical_depth.astype(np.float)
    modis_06.height[~hmask] = modis_06.height[~hmask] * h_gain + h_intercept
    modis_06.pressure[~pmask] = modis_06.pressure[~pmask] * p_gain + p_intercept
    modis_06.temperature[~tmask] = modis_06.temperature[~tmask] * t_gain - t_intercept * t_gain
    modis_06.lwp[~lwpmask] = modis_06.lwp[~lwpmask] * lwp_gain + lwp_intercept
    modis_06.optical_depth[~odmask] = modis_06.optical_depth[~odmask] * od_gain + od_intercept
    modis_06.height[hmask] = ATRAIN_MATCH_NODATA
    modis_06.pressure[pmask] = ATRAIN_MATCH_NODATA
    modis_06.temperature[tmask] = ATRAIN_MATCH_NODATA
    modis_06.lwp[lwpmask] = ATRAIN_MATCH_NODATA
    modis_06.optical_depth[odmask] = ATRAIN_MATCH_NODATA
    h5file.close()
    return modis_06


def read_modis_hdf(filename):
    from pyhdf.SD import SD, SDC
    h4file = SD(filename, SDC.READ)
    datasets = h4file.datasets()
    modis_06 = MOD06Obj()
    modis_06.height = h4file.select('cloud_top_height_1km').get()
    modis_06.temperature = h4file.select('cloud_top_temperature_1km').get()
    modis_06.pressure = h4file.select('cloud_top_pressure_1km').get()
    modis_06.lwp = h4file.select('Cloud_Water_Path').get()
    modis_06.multilayer = h4file.select('Cloud_Multi_Layer_Flag').get()
    modis_06.optical_depth = h4file.select('Cloud_Optical_Thickness').get()
    modis_06.cloud_emissivity = h4file.select('cloud_emissivity_1km').get().astype(np.float)
    modis_06.cloud_phase = h4file.select('Cloud_Phase_Infrared_1km').get().astype(np.int)
    modis_06.latitude = h4file.select('Latitude').get().astype(np.float)
    modis_06.longitude = h4file.select('Longitude').get().astype(np.float)

    h_gain = h4file.select('cloud_top_height_1km').attributes()['scale_factor']
    h_intercept = h4file.select('cloud_top_height_1km').attributes()['add_offset']
    t_gain = h4file.select('cloud_top_temperature_1km').attributes()['scale_factor']
    t_intercept = h4file.select('cloud_top_temperature_1km').attributes()['add_offset']
    p_gain = h4file.select('cloud_top_pressure_1km').attributes()['scale_factor']
    p_intercept = h4file.select('cloud_top_pressure_1km').attributes()['add_offset']
    h_nodata = h4file.select('cloud_top_height_1km').attributes()['_FillValue']
    t_nodata = h4file.select('cloud_top_temperature_1km').attributes()['_FillValue']
    p_nodata = h4file.select('cloud_top_pressure_1km').attributes()['_FillValue']

    lwp_gain = h4file.select('Cloud_Water_Path').attributes()['scale_factor']
    lwp_intercept = h4file.select('Cloud_Water_Path').attributes()['add_offset']
    lwp_nodata = h4file.select('Cloud_Water_Path').attributes()['_FillValue']

    od_gain = h4file.select('Cloud_Optical_Thickness').attributes()['scale_factor']
    od_intercept = h4file.select('Cloud_Optical_Thickness').attributes()['add_offset']
    od_nodata = h4file.select('Cloud_Optical_Thickness').attributes()['_FillValue']

    hmask = modis_06.height == h_nodata
    tmask = modis_06.temperature == t_nodata
    pmask = modis_06.pressure == p_nodata
    lwpmask = modis_06.lwp == lwp_nodata
    odmask = modis_06.optical_depth == od_nodata

    modis_06.height = modis_06.height.astype(np.float)
    modis_06.pressure = modis_06.pressure.astype(np.float)
    modis_06.temperature = modis_06.temperature.astype(np.float)
    modis_06.lwp = modis_06.lwp.astype(np.float)
    modis_06.optical_depth = modis_06.optical_depth.astype(np.float)
    modis_06.height[~hmask] = modis_06.height[~hmask] * h_gain + h_intercept
    modis_06.pressure[~pmask] = modis_06.pressure[~pmask] * p_gain + p_intercept
    modis_06.temperature[~tmask] = modis_06.temperature[~tmask] * t_gain - t_intercept * t_gain
    modis_06.lwp[~lwpmask] = modis_06.lwp[~lwpmask] * lwp_gain + lwp_intercept
    modis_06.optical_depth[~odmask] = modis_06.optical_depth[~odmask] * od_gain + od_intercept
    modis_06.height[hmask] = ATRAIN_MATCH_NODATA
    modis_06.pressure[pmask] = ATRAIN_MATCH_NODATA
    modis_06.temperature[tmask] = ATRAIN_MATCH_NODATA
    modis_06.lwp[lwpmask] = ATRAIN_MATCH_NODATA
    modis_06.optical_depth[odmask] = ATRAIN_MATCH_NODATA

    return modis_06
