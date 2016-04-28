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

"""Read Calipso/VIIRS/AVHRR matchup data object from hdf5 file
"""

import numpy as np

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
            
            
class ppsAvhrrObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'sec_1970': None,
            'ctth_height': None,
            'ctth_pressure': None,
            'ctth_temperature': None,
            'ctth_opaque': None,  # True if opaque retrieval was applied
            'cloudtype': None,
            'cloudtype_qflag': None,
            'cloudtype_phaseflag': None,
            'cloudtype_quality': None,
            'cloudtype_conditions': None,
            'cloudtype_status': None,
            'ctth_status': None,
            'surftemp': None,
            't500': None,
            't700': None,
            't850': None,
            't950': None,
            'ttro': None,
            'ciwv': None,
            #Imager channels currently used by pps
            'r06micron':  None,  #modis_1
            'r09micron':  None,  #modis_2
            'bt37micron': None, #modis_20
            'bt11micron': None, #modis 31
            'bt12micron': None, #modis 32
            'bt86micron': None, #modis 29
            'r22micron': None,
            'r13micron': None, #modis 26
            'r16micron':  None,  #moodis 6
            #developing channels from modis change to pps_id_name i.e r12micron 
            #when introduced into pps
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
            'fraction_of_land': None,
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
            'lwp': None,
            'emis1': None,
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
            #nwp data on segment resolution
            'segment_nwp_geoheight': None,
            'segment_nwp_moist': None,
            'segment_nwp_pressure': None,
            'segment_nwp_temp': None,
            'segment_nwp_surfaceLandTemp': None,
            'segment_nwp_surfaceSeaTemp': None,
            'segment_nwp_surfaceGeoHeight': None,
            'segment_nwp_surfaceMoist': None,
            'segment_nwp_surfacePressure': None,
            'segment_nwp_fractionOfLand': None,
            'segment_nwp_meanElevation': None,
            'segment_nwp_ptro': None,
            'segment_nwp_ttro': None,
            'segment_nwp_t850': None,
            'segment_tb11clfree_sea': None,
            'segment_tb12clfree_sea': None,
            'segment_tb11clfree_land': None,
            'segment_tb12clfree_land': None,
            'segment_tb11cloudy_surface': None,
            'segment_tb12cloudy_surface': None,
        }
        
        
class CalipsoObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'avhrr_linnum': None,
            'avhrr_pixnum': None,
            'elevation': None,
            'longitude': None,
            'latitude': None,
            'cloud_fraction': None,
            'layer_top_altitude': None,
            'midlayer_temperature': None,
            'layer_base_altitude': None,
            'number_layers_found': None,
            'igbp_surface_type': None,
            'nsidc_surface_type': None,
            'nsidc_surface_type_texture': None,
            'profile_utc_time': None, 
            'sec_1970': None,
            'profile_time_tai': None,
            'feature_classification_flags': None,
            'day_night_flag': None,
            'feature_optical_depth_532': None,
            'single_shot_cloud_cleared_fraction': None,
            #If a combination of 5 and 1km data are used for RESOLUTION=1
            #A vector with the corresponding optical thickness for 5km data
            # is stored also for 1km data. Because of that I put the 5km in the name
            #/2013-08-17/Nina
            'feature_optical_depth_532_top_layer_5km': None,
            'detection_height_5km': None,
            'total_optical_depth_5km': None
            }



class CalipsoAvhrrTrackObject:
    def __init__(self):
        self.avhrr = ppsAvhrrObject()
        self.calipso = CalipsoObject()
        self.calipso_aerosol = CalipsoObject()
        self.diff_sec_1970 = None

    def make_nsidc_surface_type_texture(self, kernel_sz = 51):
        """Derive the stdv of the ice dataset"""

        if self.calipso.all_arrays['nsidc_surface_type'] is not None:
            self.calipso.all_arrays['nsidc_surface_type_texture'] = sliding_std(
                self.calipso.all_arrays['nsidc_surface_type'], kernel_sz)
    
    def __add__(self, other):
        """Concatenating two objects together"""
        self.avhrr = self.avhrr + other.avhrr
        self.calipso = self.calipso + other.calipso
        self.calipso_aerosol = self.calipso_aerosol + other.calipso_aerosol
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
  
        
def readCaliopAvhrrMatchObjOldFormat(h5file, retv):
    #print "OLD FORMAT"
    for group, data_obj in [(h5file['/calipso'], retv.calipso),
                            (h5file['/avhrr'], retv.avhrr)]:
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

def readCaliopAvhrrMatchObjNewFormat(h5file, retv, var_to_read=None):
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
        data_objects.append(retv.avhrr)
    if 'cci' in h5file.keys():
        h5_groups.append(h5file['/cci'])
        data_objects.append(retv.avhrr)

    for group, data_obj in zip(h5_groups, data_objects):
        for dataset in group.keys():  
            atrain_match_name = dataset
            if atrain_match_name in data_obj.all_arrays.keys():
                if var_to_read is not None:
                  if atrain_match_name not in var_to_read:  
                      pass
                data_obj.all_arrays[atrain_match_name] = group[dataset].value
    return retv            
# ----------------------------------------

def readCaliopAvhrrMatchObj(filename, var_to_read=None):
    import h5py          
    retv = CalipsoAvhrrTrackObject()    
    h5file = h5py.File(filename, 'r')
    if "cloud_top_profile" in h5file['/calipso'].keys():
        retv = readCaliopAvhrrMatchObjOldFormat(h5file, retv, var_to_read=None)
        #print "OLD FORMAT"
    else:
        retv = readCaliopAvhrrMatchObjNewFormat(h5file, retv, var_to_read=None)
    retv.diff_sec_1970 = h5file['diff_sec_1970'].value
    h5file.close()
    retv.make_nsidc_surface_type_texture()
    return retv

# ----------------------------------------
def writeCaliopAvhrrMatchObj(filename, ca_obj, avhrr_obj_name = 'pps'):
    """
    Write *ca_obj* to *filename*.
    
    """
    from common import write_match_objects    
    groups = {'calipso': ca_obj.calipso.all_arrays,
              'calipso_aerosol': ca_obj.calipso_aerosol.all_arrays,
              avhrr_obj_name: ca_obj.avhrr.all_arrays}
    write_match_objects(filename, ca_obj.diff_sec_1970, groups)    

    status = 1
    return status


def sliding_std(x, size=5):
    """derive a sliding standard deviation of a data array"""
    from scipy.ndimage.filters import uniform_filter
    c1 = uniform_filter(x.astype('float'), size=size)
    c2 = uniform_filter(x.astype('float')*x.astype('float'), size=size)
    return abs(c2 - c1*c1)**.5



# ----------------------------------------
if __name__ == "__main__":

    import os.path
    TESTDIR = ("/local_disk/laptop/NowcastingSaf/FA/cloud_week_2013may" + 
               "/atrain_matchdata/2012/10/arctic_europe_1km")
    TESTFILE = os.path.join(TESTDIR, 
                            "1km_npp_20121012_1246_04968_caliop_viirs_match.h5")
    TESTFILE2 = os.path.join(TESTDIR,
                             "1km_npp_20121004_0700_04851_caliop_viirs_match.h5")
    caObj = readCaliopAvhrrMatchObj(TESTFILE)
    caObj2 = readCaliopAvhrrMatchObj(TESTFILE2)

    caObj = caObj + caObj2
