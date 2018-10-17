#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Read Calipso/VIIRS/IMAGER matchup data object from hdf5 file
"""

import numpy as np
import h5py   
from common import write_match_objects

class DataObject(object):
    """
    Class to handle data objects with several arrays.
    
    """
    def __getattr__(self, name):
        try:
            return self.all_arrays[name]
        except KeyError:
            raise AttributeError("%s instance has no attribute '%s'" % (
                self.__class__.__name__, name))
    
    def __setattr__(self, name, value):
        if name == 'all_arrays':
            object.__setattr__(self, name, value)
        else:
            self.all_arrays[name] = value

    def __add__(self, other):
        """Adding two objects together"""
        for key in self.all_arrays:
            try:
                if self.all_arrays[key].ndim != self.all_arrays[key].ndim:
                    raise ValueError("Can't concatenate arrays " + 
                                     "of different dimensions!")
            except AttributeError, e:
                #print "Don't concatenate member " + key + "... " + str(e)
                self.all_arrays[key] = other.all_arrays[key]
                continue
            try:
                if self.all_arrays[key].ndim == 1:  
                    self.all_arrays[key] = np.concatenate(
                        [self.all_arrays[key],
                         other.all_arrays[key]])
                elif key in ['segment_nwp_geoheight',
                             'segment_nwp_moist',
                             'segment_nwp_pressure',
                             'segment_nwp_temp']:
                     self.all_arrays[key] = np.concatenate(
                         [self.all_arrays[key],
                          other.all_arrays[key]], 0)
                elif self.all_arrays[key].ndim == 2: 
                    self.all_arrays[key] = np.concatenate(
                        [self.all_arrays[key],
                         other.all_arrays[key]], 0)
                        
            except ValueError, e:
                #print "Don't concatenate member " + key + "... " + str(e)
                self.all_arrays[key] = other.all_arrays[key]
        return self

    def mask_nodata(self, nodata):
        for key in self.all_arrays:
            if key in ['latitude']:
                pass
            else:
                try:
                    self.all_arrays[key] = np.ma.array(
                        self.all_arrays[key], 
                        mask = self.all_arrays[key]<=nodata)
                except:
                    print "cloud not mask %s"%(key)
            
            
class ppsImagerObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'imager_ctth_m_above_seasurface': None,
            'longitude': None,
            'latitude': None,
            'sec_1970': None,
            'ctth_height': None,
            'ctth_pressure': None,
            'ctth_temperature': None,
            'ctthold_height': None,
            'ctthold_pressure': None,
            'ctthold_temperature': None,
            'ctthnna1_height': None,
            'ctthnna1_pressure': None,
            'ctthnna1_temperature': None,
            'ctthnna1nt_height': None,
            'ctthnna1nt_pressure': None,
            'ctthnna1nt_temperature': None,
            'ctthnnant_height': None,
            'ctthnnant_pressure': None,
            'ctthnnant_temperature': None,
            'ctthnnvnt_height': None,
            'ctthnnvnt_pressure': None,
            'ctthnnvnt_temperature': None,
            'ctthnnm2nt_height': None,
            'ctthnnm2nt_pressure': None,
            'ctthnnm2nt_temperature': None,
            'ctthnnmint_height': None,
            'ctthnnmint_pressure': None,
            'ctthnnmint_temperature': None,
            'ctthnnmintnco2_height': None,
            'ctthnnmintnco2_pressure': None,
            'ctthnnmintnco2_temperature': None,
            'ctth_opaque': None,  # True if opaque retrieval was applied
            'cloudtype': None,
            'cloudmask': None,
            'cfc_mean': None,
            'cma_aerosolflag':None,
            'cma_aerosol':None,
            'cma_dust': None,
            'cma_prob': None,
            'cma_prob_mean': None,
            'cma_testlist0':None,
            'cma_testlist1':None,
            'cma_testlist2':None,
            'cma_testlist3':None,
            'cma_testlist4':None,
            'cma_testlist5':None,
            'cpp_cot':None,
            'cpp_cwp':None,
            'cpp_dcot':None,
            'cpp_iwp':None,
            'cpp_lwp':None,
            'cpp_phase':None,
            'cpp_reff':None,
            'seaice':None,
            'snowd':None,
            'snowa':None,
            'cloudtype_qflag': None,
            'cloudtype_phaseflag': None,
            'cloudtype_quality': None,
            'cloudtype_conditions': None,
            'cloudtype_status': None,
            'ctth_status': None,
            'surftemp': None,
            't250': None,
            't500': None,
            't700': None,
            't800': None,
            't850': None,
            't900': None,
            't950': None,
            't1000': None,
            't2m': None,
            'h2m': None,
            'u10m': None,
            'v10m': None,
            'ttro': None,
            'ptro': None,
            'psur': None,
            'ciwv': None,
            #Imager channels currently used by pps
            'r06micron':  None,  #modis_1
            'r09micron':  None,  #modis_2
            'bt37micron': None,  #modis_20
            'bt11micron': None,  #modis 31
            'bt12micron': None,  #modis 32
            'bt86micron': None,  #modis 29
            'r22micron': None,
            'r13micron': None,   #modis 26
            'r16micron':  None,  #moodis 6
            'r06micron_sza_correction_done':  None,  #modis_1
            'r09micron_sza_correction_done':  None,  #modis_2
            'r13micron_sza_correction_done':  None,  #modis 26
            'r16micron_sza_correction_done':  None,  #modis 6
            'r22micron_sza_correction_done':  None,  #VIIRS
            'r37_sza_correction_done':  None,  #VIIRS
            'r37': None,
            #developing channels from modis change to pps_id_name i.e r12micron 
            #when introduced into pps
            'modis_3_sza_correction_done':  None,
            'modis_4_sza_correction_done':  None,
            'modis_5_sza_correction_done':  None,
            'modis_7_sza_correction_done':  None,
            'modis_8_sza_correction_done': None,
            'modis_3':  None,
            'modis_4':  None,
            'modis_5':  None,
            'modis_7':  None,
            'modis_8':  None,
            'modis_9':  None,
            'modis_10': None,
            'modis_11': None,
            'modis 12': None,
            'modis_13lo': None,
            'modis_13hi': None,
            'modis_14lo': None,
            'modis_14hi': None,
            'modis_15': None,
            'modis 16': None,
            'modis_17': None,
            'modis_18': None,
            'modis_19': None,
            'modis 21': None,
            'modis_22': None,
            'modis_23': None,
            'modis_24': None,
            'modis_25': None,
            'modis_27': None,
            'modis_28': None,
            'modis_30': None,
            'modis_33': None,
            'modis_34': None,
            'modis_35': None,
            'modis_36': None,
            'seviri_bt134': None,
            'seviri_bt73': None,
            'seviri_bt134': None,
            'seviri_bt97': None,
            'seviri_bt62': None,

            'fraction_of_land': None,
            'fractionofland': None,
            'landuse': None,
            'text_r06': None,
            'text_t11': None,
            'text_t37t12': None,
            'text_t11t12': None,
            'text_t37': None,
            'thr_t11ts_inv': None,
            'thr_t11t37_inv': None,
            'thr_t37t12_inv': None,
            'thr_t11t12_inv': None,
            'thr_t85t11_inv': None,
            'thr_t11ts': None,
            'thr_t11t37': None,
            'thr_t37t12': None,
            'thr_t11t12': None,  
            'thr_t85t11': None, 
            'thr_r06': None, 
            'thr_r09': None, 
            'satz': None,
            'sunz': None,
            'azidiff': None,
            'satazimuth': None,
            'sunazimuth': None,
            'lwp': None,
            'emis1': None,
            'emis6': None,
            'emis8': None,
            'emis9': None,
            'warmest_t11': None,
            'warmest_t12': None,
            'warmest_t37': None,
            'warmest_r09': None,
            'warmest_r06': None,
            'warmest_r16': None,
            'darkest_t11': None,
            'darkest_t12': None,
            'darkest_t37': None,
            'darkest_r09': None,
            'darkest_r06': None,
            'darkest_r16': None,
            'coldest_t11': None,
            'coldest_t12': None,
            'coldest_t37': None,
            'coldest_r09': None,
            'coldest_r06': None,
            'coldest_r16': None,
            'warmest_r06_correction_done': None, 
            'warmest_r09_correction_done': None, 
            'warmest_r16_correction_done': None,
            'darkest_r06_correction_done': None, 
            'darkest_r09_correction_done': None, 
            'darkest_r16_correction_done': None,
            'coldest_r06_correction_done': None, 
            'coldest_r09_correction_done': None, 
            'coldest_r16_correction_done': None,
            #nwp data on segment resolution
            #'segment_nwp_geoheight': None,
            #'segment_nwp_moist': None,
            #'segment_nwp_pressure': None,
            #'segment_nwp_temp': None,
            'segment_nwp_h440': None,
            'segment_nwp_h680': None,
            'segment_nwp_surfaceLandTemp': None,
            'segment_nwp_surfaceSeaTemp': None,
            'segment_nwp_surfaceGeoHeight': None,
            'segment_nwp_surfaceMoist': None,
            'segment_nwp_surfacePressure': None,
            'segment_nwp_fractionOfLand': None,
            'segment_nwp_meanElevation': None,
            #'segment_nwp_ptro': None,
            #'segment_nwp_ttro': None,
            #'segment_nwp_t850': None,
            #'segment_nwp_tb11clfree_sea': None,
            #'segment_nwp_tb12clfree_sea': None,
            #'segment_nwp_tb11clfree_land': None,
            #'segment_nwp_tb12clfree_land': None,
            #'segment_nwp_tb11cloudy_surface': None,
            #'segment_nwp_tb12cloudy_surface': None,
        }
        
class ModisObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'height': None,
            'temperature': None,
            'pressure': None,
            'cloud_emissivity': None,
            'cloud_phase':None,
            'lwp': None}

        
class CalipsoObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'cal_MODIS_cflag': None,
            'cloudsat_index': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'elevation': None,
            'longitude': None,
            'latitude': None,
            'cloud_fraction': None,
            'validation_height': None,
            'layer_top_altitude': None,
            'layer_top_temperature': None,
            'layer_top_pressure': None,
            'midlayer_temperature': None,
            'layer_base_altitude': None,
            'layer_base_pressure': None,
            'number_layers_found': None,
            'igbp_surface_type': None,
            'nsidc_surface_type': None,
            'snow_ice_surface_type': None,
            #'nsidc_surface_type_texture': None,
            'profile_utc_time': None, 
            'sec_1970': None,
            'profile_time_tai': None,
            'feature_classification_flags': None,
            'day_night_flag': None,
            'feature_optical_depth_532': None,
            'single_shot_cloud_cleared_fraction': None,
            'cfc_single_shots_1km_from_5km_file': None,
            'column_optical_depth_tropospheric_aerosols_532': None,
            'column_optical_depth_tropospheric_aerosols_532_5km': None,
            'column_optical_depth_aerosols_532': None,
            'column_optical_depth_aerosols_532_5km': None,
            'profile_id':None,
            #If a combination of 5 and 1km data are used for RESOLUTION=1
            #A vector with the corresponding optical thickness for 5km data
            # is stored also for 1km data. Because of that I put the 5km in the name
            #/2013-08-17/Nina
            'feature_optical_depth_532_top_layer_5km': None,
            'detection_height_5km': None,
            'total_optical_depth_5km': None,
            #'layer_top_altitude_5km': None,
            #'layer_base_altitude_5km': None,
            'tropopause_height': None
            }

class CloudsatObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'clsat_max_height':None,
            'longitude': None,
            'latitude': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'cloud_fraction': None,
            'validation_height': None,
            'validation_height_base': None,
            'elevation': None,
            'Profile_time': None,
            'sec_1970': None,
            'TAI_start': None,
            #'Temp_min_mixph_K': None,
            #'Temp_max_mixph_K': None,
            # The data:
            'CPR_Cloud_mask': None,
            #'CPR_Echo_Top': None,
            #'Clutter_reduction_flag': None,
            #'Data_quality': None,
            #'Data_targetID': None,
            #'Gaseous_Attenuation': None,
            'MODIS_Cloud_Fraction': None,
            'MODIS_cloud_flag': None,
            #'Radar_Reflectivity': None,
            'Height': None,
            #'SigmaZero': None,
            #'SurfaceHeightBin': None,
            #'SurfaceHeightBin_fraction': None,
            #'sem_NoiseFloor': None,
            #'sem_NoiseFloorVar': None,
            #'sem_NoiseGate': None,
            'RVOD_liq_water_path': None,
            'RVOD_ice_water_path': None,
            'RVOD_CWC_status': None,
            'LO_RVOD_liquid_water_path': None,
            'IO_RVOD_ice_water_path': None,
            'RVOD_liq_water_content': None,
            'RVOD_ice_water_content': None,
            #'RVOD_liq_water_path_uncertainty': None,
            #'RVOD_ice_water_path_uncertainty': None,
            #'LO_RVOD_liquid_water_path': None,
            #'LO_RVOD_liquid_water_path_uncertainty': None,
            #'IO_RVOD_ice_water_path': None,
            #'IO_RVOD_ice_water_path_uncertainty': None,
            #'RVOD_liq_water_content': None,
            #'RVOD_liq_water_content_uncertainty': None,
            #'RVOD_ice_water_content': None,
            #'RVOD_ice_water_content_uncertainty': None,
            #'LO_RVOD_liquid_water_content': None,
            #'LO_RVOD_liquid_water_content_uncertainty': None,
            #'IO_RVOD_ice_water_content': None,
            #'IO_RVOD_ice_water_content_uncertainty': None,
            #'RVOD_CWC_status': None
            'calipso_layer_base_altitude': None,
            'calipso_layer_top_altitude': None,
            'calipso_feature_classification_flags': None
                           }

class IssObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'sec_1970': None,
            'elevation': None,
            'cloud_fraction': None,
            'validation_height': None,
            "cats_fore_fov_angle": None, #(2319,)
            "cats_fore_fov_latitude": None, #(2319, 3)
            "cats_fore_fov_longitude": None, #(2319, 3)
            "index_top_bin_fore_fov": None, #(2319,)
            "solar_azimuth_angle": None, #(2319,)
            "solar_zenith_angle": None, #(2319,)
            "aerosol_type_fore_fov": None, #(2319, 10)
            "cloud_phase_fore_fov": None, #(2319, 10)
            "cloud_phase_score_fore_fov": None, #(2319, 10)
            "constrained_lidar_ratio_flag": None, #(3, 2319, 10, 2)
            "dem_surface_altitude_fore_fov": None, #(2319,)
            "day_night_flag": None, #(2319,)
            "extinction_qc_flag_1064_fore_fov": None, #(2319, 10)
            "feature_type_fore_fov": None, #(2319, 10)
            "feature_type_score_fore_fov": None, #(2319, 10)
            "layer_base_altitude_fore_fov": None, #(2319, 10)
            "layer_base_bin_fore_fov": None, #(2319, 10)
            "layer_base_pressure_fore_fov": None, #(2319, 10)
            "layer_base_temperature_fore_fov": None, #(2319, 10)
            "layer_effective_multiple_scattering_factor_1064_fore_fov": None, #(2319, 10)
            "layer_top_altitude_fore_fov": None, #(2319, 10)
            "layer_top_bin_fore_fov": None, #(2319, 10)
            "layer_top_pressure_fore_fov": None, #(2319, 10)
            "layer_top_temperature_fore_fov": None, #(2319, 10)
            "lidar_ratio_selection_method_1064_fore_fov": None, #(2319, 10)
            "lidar_surface_altitude_fore_fov": None, #(2319,)
            "number_layers_fore_fov": None, #(2319,)
            "opacity_fore_fov": None, #(2319, 10)
            "profile_id": None, #(2319,)
            "profile_utc_date": None, #(2319,)
            "profile_utc_time": None, #(2319, 3)
            "sky_condition_fore_fov": None, #(2319,)
            "surface_type_fore_fov": None, #(2319,)
            "bin_altitude_array": None, #(533,)
            "feature_optical_depth_1064_fore_fov": None, #(2319, 10)
            "feature_optical_depth_uncertainty_1064_fore_fov": None, #(2319, 10)
            "ice_water_path_1064_fore_fov": None, #(2319, 10)
            "ice_water_path_1064_uncertainty_fore_fov": None, #(2319, 10)
            "integrated_attenauted_backscatter_1064_fore_fov": None, #(2319, 10)
            "integrated_attenauted_backscatter_532_fore_fov": None, #(2319, 10)
            "integrated_attenauted_backscatter_uncertainty_1064_fore_fov": None, #(2319, 10)
            "integrated_attenauted_backscatter_uncertainty_532_fore_fov": None, #(2319, 10)
            "integrated_attenuated_total_color_ratio_fore_fov": None, #(2319, 10)
            "integrated_attenuated_total_color_ratio_uncertainty_fore_fov": None, #(2319, 10)
            "integrated_volume_depolarization_ratio_1064_fore_fov": None, #(2319, 10)
            "integrated_volume_depolarization_ratio_uncertainty_1064_fore_fov": None, #(2319, 10)
            "lidar_ratio_1064_fore_fov": None, #(2319, 10)
            "measured_two_way_transmittance_1064_fore_fov": None, #(2319, 10)
            "measured_two_way_transmittance_uncertainty_1064_fore_fov": None, #(2319, 10)
            "two_way_transmittance_measurement_region_fore_fov": None, #(2319, 10)            
            #"attenuated_backscatter_statistics_1064_fore_fov": None, #(2319, 10, 4)
            #"attenuated_backscatter_statistics_532_fore_fov": None, #(2319, 10, 4)
            #"attenuated_total_color_ratio_statistics_fore_fov": None, #(2319, 10, 4)
            #"volume_depolarization_ratio_statistics_1064_fore_fov": None, #(2319, 10, 4)
        }

class AmsrObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'imager_linnum_nneigh': None,
            'imager_pixnum_nneigh': None,
            'sec_1970': None,
            'lwp': None}

class MoraObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'cloud_base_height': None,
            'sec_1970': None}

class SynopObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'imager_linnum_nneigh': None,
            'imager_pixnum_nneigh': None,
            'total_cloud_cover': None,
            'cloud_fraction': None,
            'nh':None,
            'cl': None,
            'cm': None,
            'ch': None,
            'vvv': None,
            'ww': None,
            'temp': None,
            'dtemp': None,
            'sec_1970': None,
            'pressure': None}

class SynopImagerTrackObject:
    def __init__(self):
        self.imager = ppsImagerObject()
        self.modis = ModisObject()
        self.synop = SynopObject()
        self.diff_sec_1970 = None
        self.truth_sat = 'synop' # Well not a satelite, but ...
    def __add__(self, other):
        """Concatenating two objects together"""
        self.imager = self.imager + other.imager
        self.synop = self.synop + other.synop
        #self.modis = self.modis + other.modis
        try:
            self.diff_sec_1970 = np.concatenate([self.diff_sec_1970,
                                                 other.diff_sec_1970])
        except ValueError, e:
            #print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = other.diff_sec_1970
        return self

class AmsrImagerTrackObject:
    def __init__(self):
        self.imager = ppsImagerObject()
        self.modis = ModisObject()
        self.amsr = AmsrObject()
        self.diff_sec_1970 = None
        self.truth_sat = 'amsr' #Satellite is EOS-Aqua or EOS-Terra
                                #As we also use MODIS data from those
                                #Name the truth_sat after the instrument
    def __add__(self, other):
        """Concatenating two objects together"""
        self.imager = self.imager + other.imager
        self.amsr = self.amsr + other.amsr
        #self.modis = self.modis + other.modis
        try:
            self.diff_sec_1970 = np.concatenate([self.diff_sec_1970,
                                                 other.diff_sec_1970])
        except ValueError, e:
            #print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = other.diff_sec_1970

        return self


class MoraImagerTrackObject:
    def __init__(self):
        self.imager = ppsImagerObject()
        self.modis = ModisObject()
        self.mora = MoraObject()
        self.diff_sec_1970 = None
        self.truth_sat = 'mora' #Satellite is EOS-Aqua or EOS-Terra
                                #As we also use MODIS data from those
                                #Name the truth_sat after the instrument
    def __add__(self, other):
        """Concatenating two objects together"""
        self.imager = self.imager + other.imager
        self.mora = self.mora + other.mora
        #self.modis = self.modis + other.modis
        try:
            self.diff_sec_1970 = np.concatenate([self.diff_sec_1970,
                                                 other.diff_sec_1970])
        except ValueError, e:
            #print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = other.diff_sec_1970

        return self

class IssImagerTrackObject:
    def __init__(self):
        self.imager=ppsImagerObject()
        self.modis = ModisObject()
        self.iss=IssObject()
        self.diff_sec_1970=None
        self.truth_sat = 'iss'
    def __add__(self, other):
        """Concatenating two objects together"""
        self.imager = self.imager + other.imager
        self.iss = self.iss + other.iss
        self.modis = self.modis + other.modis
        try:
            self.diff_sec_1970 = np.concatenate([self.diff_sec_1970,
                                                 other.diff_sec_1970])
        except ValueError, e:
            #print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = other.diff_sec_1970

        return self

class CloudsatImagerTrackObject:
    def __init__(self):
        self.imager=ppsImagerObject()
        self.modis = ModisObject()
        self.cloudsat=CloudsatObject()
        self.diff_sec_1970=None
        self.truth_sat = 'cloudsat'
    def __add__(self, other):
        """Concatenating two objects together"""
        self.imager = self.imager + other.imager
        self.cloudsat = self.cloudsat + other.cloudsat
        self.modis = self.modis + other.modis
        try:
            self.diff_sec_1970 = np.concatenate([self.diff_sec_1970,
                                                 other.diff_sec_1970])
        except ValueError, e:
            #print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = other.diff_sec_1970

        return self

class CalipsoImagerTrackObject:
    def __init__(self):
        self.imager = ppsImagerObject()
        self.modis = ModisObject()
        self.calipso = CalipsoObject()
        self.calipso_aerosol = CalipsoObject()
        self.diff_sec_1970 = None
        self.truth_sat = 'calipso'

    def make_nsidc_surface_type_texture(self, kernel_sz = 51):
        """Derive the stdv of the ice dataset"""

        if self.calipso.all_arrays['nsidc_surface_type'] is not None:
            self.calipso.all_arrays['nsidc_surface_type_texture'] = sliding_std(
                self.calipso.all_arrays['nsidc_surface_type'], kernel_sz)
    
    def __add__(self, other):
        """Concatenating two objects together"""
        self.imager = self.imager + other.imager
        self.calipso = self.calipso + other.calipso
        self.calipso_aerosol = self.calipso_aerosol + other.calipso_aerosol
        self.modis = self.modis + other.modis
        try:
            self.diff_sec_1970 = np.concatenate([self.diff_sec_1970,
                                                 other.diff_sec_1970])
        except ValueError, e:
            #print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = other.diff_sec_1970

        return self

"""
These variables belonging to calipso object now have new names.
They now keep their name from the calipso file.
Except for the varaible profile_time which is called profile_time_tai.
So not to forget it is tai time in it.
Here we remember what names they used to have to be able to 
reprocess old reshaped-files.
When reprocessing old reshaped-files,
which we might want to do, these are needed.
"""
traditional_atrain_match_to_new_names ={
    #"time":                      "profile_time_tai",
    #"utc_time":                  "profile_utc_time",
    #these were never used:
    #"optical_depth_uncertainty": "feature_optical_depth_uncertainty_532",

    "cloud_mid_temperature": "midlayer_temperature",
    #"ice_water_path5km": "ice_water_path",
    #"ice_water_path_uncertainty5km": "ice_water_path_uncertainty",
    #"horizontal_averaging5km",: "horizontal_averaging",
    #"opacity5km: "opacity_flag",

    "cloud_top_profile_pressure": "layer_top_pressure",
    "cloud_top_profile":          "layer_top_altitude",
    "cloud_base_profile":         "layer_base_altitude",
    "number_of_layers_found":     "number_layers_found",
    "igbp":                       "igbp_surface_type",
    "nsidc":                      "nsidc_surface_type",
    "optical_depth":              "feature_optical_depth_532",
    "optical_depth_top_layer5km": "feature_optical_depth_532_top_layer_5km"
    }
  
        
def readCaliopImagerMatchObjOldFormat(h5file, retv, var_to_read=None):
    #print "OLD FORMAT"
    for group, data_obj in [(h5file['/calipso'], retv.calipso),
                            (h5file['/imager'], retv.imager)]:
        for dataset in group.keys():  
            atrain_match_name = dataset
            if (dataset in traditional_atrain_match_to_new_names.keys()):
                atrain_match_name = traditional_atrain_match_to_new_names[dataset]  
            if atrain_match_name in data_obj.all_arrays.keys():
                the_data = group[dataset].value
                if dataset in ["cloud_top_profile",
                               "cloud_base_profile",
                               "cloud_base_profile_pressure",
                               "cloud_mid_temperature",
                               #"horizontal_averaging5km",
                               #"ice_water_path5km",
                               #"ice_water_path_uncertainty5km",
                               #"opacity5km",
                               #"optical_depth_uncertainty",
                               "optical_depth",
                               "single_shot_cloud_cleared_fraction",
                               "lidar_surface_elevation",
                               "feature_classification_flags"]:
                    the_data = the_data.transpose()
                data_obj.all_arrays[atrain_match_name] = the_data
    return retv

def get_stuff_to_read_from_a_reshaped_file(h5file, retv):
    h5_groups = []
    data_objects = []
    if 'calipso' in h5file.keys():
        h5_groups.append(h5file['/calipso'])
        data_objects.append(retv.calipso)
    if 'calipso_aerosol' in h5file.keys():
        h5_groups.append(h5file['/calipso_aerosol'])
        data_objects.append(retv.calipso_aerosol)
    if 'pps' in h5file.keys():
        h5_groups.append(h5file['/pps'])
        data_objects.append(retv.imager)
    if 'cci' in h5file.keys():
        h5_groups.append(h5file['/cci'])
        data_objects.append(retv.imager)
    if 'maia' in h5file.keys():
        h5_groups.append(h5file['/maia'])
        data_objects.append(retv.imager)
    if 'modis_lvl2' in h5file.keys():
        h5_groups.append(h5file['/modis_lvl2'])
        data_objects.append(retv.modis)        
    if 'cloudsat' in  h5file.keys():
        h5_groups.append(h5file['/cloudsat'])
        data_objects.append(retv.cloudsat)
    if 'iss' in  h5file.keys():
        h5_groups.append(h5file['/iss'])
        data_objects.append(retv.iss)
    if 'amsr' in  h5file.keys():
        h5_groups.append(h5file['/amsr'])
        data_objects.append(retv.amsr)
    if 'mora' in  h5file.keys():
        h5_groups.append(h5file['/mora'])
        data_objects.append(retv.mora)
    if 'synop' in  h5file.keys():
        h5_groups.append(h5file['/synop'])
        data_objects.append(retv.synop)
    return (h5_groups, data_objects)
    
def readCaliopImagerMatchObjNewFormat(h5file, retv, var_to_read=None, var_to_skip=None):
    (h5_groups, data_objects) =  get_stuff_to_read_from_a_reshaped_file(h5file, retv)
    for group, data_obj in zip(h5_groups, data_objects):
        for dataset in group.keys():  
            atrain_match_name = dataset
            if atrain_match_name in data_obj.all_arrays.keys():
                if var_to_read is not None:
                  if atrain_match_name not in var_to_read:  
                      continue
                if var_to_skip is not None:
                    if var_to_skip in atrain_match_name:
                        #print "skipping",atrain_match_name 
                        continue  
                if atrain_match_name in ["snow_ice_surface_type"]:
                    atrain_match_name = "nsidc_surface_type"
                data_obj.all_arrays[atrain_match_name] = group[dataset].value
    return retv            

def readCaliopImagerMatchObj(filename, var_to_read=None, var_to_skip=None):
    retv = CalipsoImagerTrackObject()    
    h5file = h5py.File(filename, 'r')
    if "cloud_top_profile" in h5file['/calipso'].keys():
        retv = readCaliopImagerMatchObjOldFormat(h5file, retv, var_to_read=None)
        #print "OLD FORMAT"
    else:
        retv = readCaliopImagerMatchObjNewFormat(h5file, retv, var_to_read=None, var_to_skip=var_to_skip)
    retv.diff_sec_1970 = h5file['diff_sec_1970'].value
    h5file.close()
    #retv.make_nsidc_surface_type_texture()
    return retv

def readTruthImagerMatchObj(filename, retv):
    h5file = h5py.File(filename, 'r')
    (h5_groups, data_objects) =  get_stuff_to_read_from_a_reshaped_file(h5file, retv)
    for group, data_obj in zip(h5_groups, data_objects):
        for dataset in group.keys():        
            if dataset in data_obj.all_arrays.keys():
                data_obj.all_arrays[dataset] = group[dataset].value 
    retv.diff_sec_1970 = h5file['diff_sec_1970'].value
    h5file.close()
    return retv

def readCloudsatImagerMatchObj(filename):
    retv = CloudsatImagerTrackObject()  
    return readTruthImagerMatchObj(filename, retv)
def readIssImagerMatchObj(filename): 
    retv = IssImagerTrackObject()   
    return readTruthImagerMatchObj(filename, retv)
def readAmsrImagerMatchObj(filename): 
    retv = AmsrImagerTrackObject()   
    return readTruthImagerMatchObj(filename, retv)
def readMoraImagerMatchObj(filename): 
    retv = MoraImagerTrackObject()   
    return readTruthImagerMatchObj(filename, retv)
def readSynopImagerMatchObj(filename): 
    retv = SynopImagerTrackObject()   
    return readTruthImagerMatchObj(filename, retv)

def read_files(files, truth='calipso'):
    if 'cali' in truth:
        tObj = CalipsoImagerTrackObject()  
    if 'cloudsat' in  truth:
        tObj = CloudsatImagerTrackObject()  
    if 'iss' in  truth:
        tObj = IssImagerTrackObject()  
    if 'amsr' in  truth:
        tObj = AmsrImagerTrackObject()  
    if 'mora' in  truth:
        tObj = MoraImagerTrackObject() 
    if 'synop' in  truth:
        tObj = SynopImagerTrackObject()  
    for filename in files:
        if 'cali' in truth:
            tObj += readCaliopImagerMatchObj(filename)  
        if 'cloudsat' in  truth:
            tObj += readCloudsatImagerMatchObj(filename)  
        if 'iss' in  truth:
            tObj += readIssImagerMatchObj(filename)  
        if 'amsr' in  truth:
            tObj += readAmsrImagerMatchObj(filename)
        if 'mora' in  truth:
            tObj += readMoraImagerMatchObj(filename)
        if 'synop' in  truth:
            tObj += readSynopImagerMatchObj(filename)
    return tObj 
# write matchup files
def writeCaliopImagerMatchObj(filename, ca_obj, imager_obj_name = 'pps'):
    """
    Write *ca_obj* to *filename*.    
    """
    groups = {'calipso': ca_obj.calipso.all_arrays,
              'calipso_aerosol': ca_obj.calipso_aerosol.all_arrays,
              imager_obj_name: ca_obj.imager.all_arrays,
              'modis_lvl2': ca_obj.modis.all_arrays }
    write_match_objects(filename, ca_obj.diff_sec_1970, groups)    
    status = 1
    return status

def writeCloudsatImagerMatchObj(filename,cl_obj, imager_obj_name = 'pps'):    
    groups = {'cloudsat': cl_obj.cloudsat.all_arrays,
              'modis_lvl2': cl_obj.modis.all_arrays,
              imager_obj_name: cl_obj.imager.all_arrays}
    write_match_objects(filename, cl_obj.diff_sec_1970, groups)    
    status = 1
    return status

def writeIssImagerMatchObj(filename,iss_obj, imager_obj_name = 'pps'):
    groups = {'iss': iss_obj.iss.all_arrays,
              'modis_lvl2': iss_obj.modis.all_arrays,
              imager_obj_name: iss_obj.imager.all_arrays}
    write_match_objects(filename, iss_obj.diff_sec_1970, groups)    
    status = 1
    return status

def writeAmsrImagerMatchObj(filename,amsr_obj, imager_obj_name = 'pps'):
    groups = {'amsr': amsr_obj.amsr.all_arrays,
              'modis_lvl2': amsr_obj.modis.all_arrays,
              imager_obj_name: amsr_obj.imager.all_arrays}
    write_match_objects(filename, amsr_obj.diff_sec_1970, groups)    
    status = 1
    return status

def writeMoraImagerMatchObj(filename,mora_obj, imager_obj_name = 'pps'):
    groups = {'mora': mora_obj.mora.all_arrays,
              'modis_lvl2': mora_obj.modis.all_arrays,
              imager_obj_name: mora_obj.imager.all_arrays}
    write_match_objects(filename, mora_obj.diff_sec_1970, groups)    
    status = 1
    return status

def writeSynopImagerMatchObj(filename,synop_obj, imager_obj_name = 'pps'):
    groups = {'synop': synop_obj.synop.all_arrays,
              'modis_lvl2': synop_obj.modis.all_arrays,
              imager_obj_name: synop_obj.imager.all_arrays}
    write_match_objects(filename, synop_obj.diff_sec_1970, groups)    
    status = 1
    return status

def sliding_std(x, size=5):
    """derive a sliding standard deviation of a data array"""
    from scipy.ndimage.filters import uniform_filter
    c1 = uniform_filter(x.astype('float'), size=size)
    c2 = uniform_filter(x.astype('float')*x.astype('float'), size=size)
    return abs(c2 - c1*c1)**.5


the_used_variables = [
    'imager_ctth_m_above_seasurface',
    'longitude',
    'latitude',
    'sec_1970',
    'ctth_height',
    'ctth_pressure',
    'ctth_temperature',
    'ctthold_height',
    'ctthold_pressure',
    'ctthold_temperature',
    'ctthnna1_height',
    'ctthnna1_pressure',
    'ctthnna1_temperature',
    'ctthnna1nt_height',
    'ctthnna1nt_pressure',
    'ctthnna1nt_temperature',
    'ctthnnant_height',
    'ctthnnant_pressure',
    'ctthnnant_temperature',
    'ctthnnvnt_height',
    'ctthnnvnt_pressure',
    'ctthnnvnt_temperature',
    'ctthnnm2nt_height',
    'ctthnnm2nt_pressure',
    'ctthnnm2nt_temperature',
    'ctthnnmint_height',
    'ctthnnmint_pressure',
    'ctthnnmint_temperature',
    'ctthnnmintnco2_height',
    'ctthnnmintnco2_pressure',
    'ctthnnmintnco2_temperature',
    'ctth_opaque',  # True if opaque retrieval was applied
    'cloudtype',
    'cloudmask',
    'cma_prob',
    'cfc_mean',
    'cma_prob_mean',
    'cma_aerosolflag',
    'cma_aerosol',
    'cma_dust',
    'cma_testlist0',
    'cma_testlist1',
    'cma_testlist2',
    'cma_testlist3',
    'cma_testlist4',
    'cma_testlist5',
    'cma_prob',
    'cpp_cot',
    'cpp_cwp',
    'cpp_dcot',
    'cpp_iwp',
    'cpp_lwp',
    'cpp_phase',
    'cpp_reff',
    'seaice',
    'snowd',
    'snowa',
    'cloudtype_qflag',
    'cloudtype_phaseflag',
    'cloudtype_quality',
    'cloudtype_conditions',
    'cloudtype_status',
    'ctth_status',
    'surftemp',
    't250',
    't500',
    't700',
    't800',
    't850',
    't900',
    't950',
    't1000',
    't2m',
    'h2m',
    'u10m',
    'v10m',
    'ttro',
    'ptro',
    'psur',
    'ciwv',
    #Imager channels currently used by pps
    'r06micron',  #modis_1
    'r09micron',  #modis_2
    'bt37micron',  #modis_20
    'bt11micron',  #modis 31
    'bt12micron',  #modis 32
    'bt86micron',  #modis 29
    'r22micron',
    'r13micron',   #modis 26
    'r16micron',  #moodis 6
    'r06micron_sza_correction_done',  #modis_1
    'r09micron_sza_correction_done',  #modis_2
    'r13micron_sza_correction_done',  #modis 26
    'r16micron_sza_correction_done',  #modis 6
    'r22micron_sza_correction_done',  #VIIRS
    'r37_sza_correction_done',  #VIIRS
    'r37',
    #developing channels from modis change to pps_id_name i.e r12micron 
    #when introduced into pps
    'modis_3_sza_correction_done',
    'modis_4_sza_correction_done',
    'modis_5_sza_correction_done',
    'modis_7_sza_correction_done',
    'modis_8_sza_correction_done',
    'modis_3',
    'modis_4',
    'modis_5',
    'modis_7',
    'modis_8',
    'modis_9',
    'modis_10',
    'modis_11',
    'modis 12',
    'modis_13lo',
    'modis_13hi',
    'modis_14lo',
    'modis_14hi',
    'modis_15',
    'modis 16',
    'modis_17',
    'modis_18',
    'modis_19',
    'modis 21',
    'modis_22',
    'modis_23',
    'modis_24',
    'modis_25',
    'modis_27',
    'modis_28',
    'modis_30',
    'modis_33',
    'modis_34',
    'modis_35',
    'modis_36',
    'seviri_bt134',
    'seviri_bt73',
    'seviri_bt134',
    'seviri_bt97',
    'seviri_bt62',
    'landuse',
    'fraction_of_land',
    'fractionofland',
    'text_r06',
    'text_t11',
    'text_t37t12',
    'text_t11t12',
    'text_t37',
    'thr_t11ts_inv',
    'thr_t11t37_inv',
    'thr_t37t12_inv',
    'thr_t11t12_inv',
    'thr_t85t11_inv',
    'thr_t11ts',
    'thr_t11t37',
    'thr_t37t12',
    'thr_t11t12',  
    'thr_t85t11', 
    'thr_r06', 
    'thr_r09', 
    'satz',
    'sunz',
    'azidiff',
    'satazimuth',
    'sunazimuth',
    'lwp',
    'emis1',
    'emis6',
    'emis8',
    'emis9',
    'warmest_t11',
    'warmest_t12',
    'warmest_t37',
    'warmest_r09',
    'warmest_r06',
    'warmest_r16',
    'darkest_t11',
    'darkest_t12',
    'darkest_t37',
    'darkest_r09',
    'darkest_r06',
    'darkest_r16',
    'coldest_t11',
    'coldest_t12',
    'coldest_t37',
    'coldest_r09',
    'coldest_r06',
    'coldest_r16',
    'warmest_r06_correction_done', 
    'warmest_r09_correction_done', 
    'warmest_r16_correction_done',
    'darkest_r06_correction_done', 
    'darkest_r09_correction_done', 
    'darkest_r16_correction_done',
    'coldest_r06_correction_done', 
    'coldest_r09_correction_done', 
    'coldest_r16_correction_done',
    'segment_nwp_h440',
    'segment_nwp_h680',
    'longitude',
    'latitude',
    'imager_linnum',
    'imager_pixnum',
    'imager_linnum_nneigh',
    'imager_pixnum_nneigh',
    'sec_1970',
    'elevation',
    #MODIS LVL2
    'height',
    'temperature',
    'pressure',
    'cloud_emissivity',
    'cloud_phase',
    'lwp',
    #AMSR
    'lwp',   
    'imager_amsr_dist',
    #MORA
    'cloud_base_height',   
    #Cloudsat 
    'cloud_fraction',
    'validation_height',
    'RVOD_liq_water_path',
    'RVOD_ice_water_path',
    'LO_RVOD_liquid_water_path',
    'IO_RVOD_ice_water_path',
    'RVOD_CWC_status',
    #CALIPSO write to not combine
    'cal_MODIS_cflag',
    'cloudsat_index',
    #CALIPSO only (ISS?)
    'profile_id',
    'layer_top_altitude',
    'layer_top_altitude_fore_fov',
    'layer_top_temperature',
    'layer_top_pressure',
    'layer_base_altitude',
    'layer_base_pressure',
    'midlayer_temperature',
    'number_layers_found',
    'igbp_surface_type',
    'nsidc_surface_type',
    'snow_ice_surface_type',
    'surface_type_fore_fov',
    'feature_classification_flags',
    'feature_optical_depth_532',
    'single_shot_cloud_cleared_fraction',
    'cfc_single_shots_1km_from_5km_file',
    'feature_optical_depth_532_top_layer_5km',
    'total_optical_depth_5km',
    'detection_height_5km',         
    'column_optical_depth_cloud_532',
    'column_optical_depth_cloud_uncertainty_532',
    'column_optical_depth_tropospheric_aerosols_532_5km',
    'column_optical_depth_tropospheric_aerosols_532',
    'column_optical_depth_aerosols_532_5km',
    'column_optical_depth_aerosols_532',
    'Number_cloudy_single_shots',
    #CLOUDSAT only
    'clsat_max_height',
    'validation_height_base',
    'MODIS_Cloud_Fraction',
    'MODIS_cloud_flag',
    'calipso_layer_base_altitude',
    'calipso_layer_top_altitude',
    'calipso_feature_classification_flags']



# ----------------------------------------
if __name__ == "__main__":

    import os.path
    TESTDIR = ("/local_disk/laptop/NowcastingSaf/FA/cloud_week_2013may" + 
               "/atrain_matchdata/2012/10/arctic_europe_1km")
    TESTFILE = os.path.join(TESTDIR, 
                            "1km_npp_20121012_1246_04968_caliop_viirs_match.h5")
    TESTFILE2 = os.path.join(TESTDIR,
                             "1km_npp_20121004_0700_04851_caliop_viirs_match.h5")
    caObj = readCaliopImagerMatchObj(TESTFILE)
    caObj2 = readCaliopImagerMatchObj(TESTFILE2)

    caObj = caObj + caObj2
