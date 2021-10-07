#!/usr/bin/env python
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
# -*- coding: utf-8 -*-

# Author(s):

#   Adam.Dybbroe,
#   N.Hakansson

"""Read Calipso/VIIRS/IMAGER matchup data object from hdf5 file
"""

import numpy as np
import h5py
from atrain_match.utils.common import write_match_objects


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
        # Check if we have an empty object
        # modis objects does not have longitude attribute
        is_empty_self = True
        is_empty_other = True
        for key in self.all_arrays.keys():
            if self.all_arrays[key] is not None and len(self.all_arrays[key]) > 0:
                is_empty_self = False
        for key in other.all_arrays.keys():
            if other.all_arrays[key] is not None and len(other.all_arrays[key]) > 0:
                is_empty_other = False
        if is_empty_self:
            # print("First object is None!, returning second object")
            return other
        if is_empty_other:
            # print("Second object is None!, returning first object")
            return self
        for key in self.all_arrays:
            try:
                if self.all_arrays[key].ndim != self.all_arrays[key].ndim:
                    raise ValueError("Can't concatenate arrays " +
                                     "of different dimensions!")
            except AttributeError:
                # print "Don't concatenate member " + key + "... " + str(e)
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
            except ValueError:
                # print "Don't concatenate member " + key + "... " + str(e)
                self.all_arrays[key] = other.all_arrays[key]
        return self

    def extract_elements(self, idx=None, starti=0, endi=0):
        """Extract elements with index idx"""
        # to replace calipso_track_from_matched

        for key, value in self.all_arrays.items():
            if key in ["TAI_start"]:
                continue
            if value is None:
                self.all_arrays[key] = None
            elif value.size == 1:
                pass
            elif idx is not None:
                self.all_arrays[key] = value[idx.ravel(), ...]
            else:
                self.all_arrays[key] = value[starti:endi, ...]
            if value is not None and len(value.shape) > 1 and value.shape[1] == 1:
                self.all_arrays[key] = self.all_arrays[key].ravel()

        return self

    def mask_nodata(self, nodata):
        for key in self.all_arrays:
            if key in ['latitude']:
                pass
            else:
                try:
                    self.all_arrays[key] = np.ma.array(
                        self.all_arrays[key],
                        mask=self.all_arrays[key] <= nodata)
                except:
                    print("cloud not mask %s" % (key))


class ExtractedImagerObject(DataObject):
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
            'cloudtype': None,
            'cloudmask': None,
            'cfc_mean': None,
            'cma_prob': None,
            'cma_prob_mean': None,

            'cpp_lwp': None,
            'cpp_phase': None,
            # Quality flags
            'cloudtype_qflag': None,
            'cloudtype_phaseflag': None,
            'cloudtype_quality': None,
            'cloudtype_conditions': None,
            'cloudtype_status': None,
            'ctth_status': None,
            # Angles
            'satz': None,
            'sunz': None,

        }


class ModisObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)
        self.all_arrays = {
            'height': None,
            'temperature': None,
            'pressure': None,
            'cloud_emissivity': None,
            'cloud_phase': None,
            'lwp': None}


class ExtraObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)
        self.all_arrays = {
            'name': None}


class CalipsoObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)
        self.all_arrays = {
            # Normal name = calipso.name.lower()

            # Imager matching needed for all truths:
            'longitude': None,
            'latitude': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'elevation': None,  # DEM_elevation => elevation in (m)"
            'cloud_fraction': None,
            'validation_height': None,
            'sec_1970': None,
            'minimum_laser_energy_532': None,
            'layer_top_altitude': None,
            'layer_top_temperature': None,
            'layer_top_pressure': None,
            'midlayer_temperature': None,
            'layer_base_altitude': None,
            'layer_base_pressure': None,
            'number_layers_found': None,
            'igbp_surface_type': None,
            'nsidc_surface_type': None,  # V4 renamed from 'snow_ice_surface_type'
            'snow_ice_surface_type': None,
            # 'nsidc_surface_type_texture': None,
            'profile_time_tai': None,  # renamed from "Profile_Time"
            'feature_classification_flags': None,
            'day_night_flag': None,
            'feature_optical_depth_532': None,
            'tropopause_height': None,
            'profile_id': None,

            # If a combination of 5 and 1km data are used for RESOLUTION=1
            # "column_optical_depth_tropospheric_aerosols_1064_5km": None,
            # "column_optical_depth_tropospheric_aerosols_1064": None,
            "column_optical_depth_tropospheric_aerosols_532_5km": None,
            "column_optical_depth_tropospheric_aerosols_532": None,
            "column_optical_depth_aerosols_532_5km": None,
            "column_optical_depth_aerosols_532": None,
            # "column_optical_depth_tropospheric_aerosols_uncertainty_1064_5km": None,
            # "column_optical_depth_tropospheric_aerosols_uncertainty_532_5km": None,
            "column_optical_depth_cloud_532_5km": None,
            # "column_optical_depth_cloud_uncertainty_532_5km": None,
            "feature_optical_depth_532_5km": None,
            "layer_top_altitude_5km": None,
            "layer_top_pressure_5km": None,
            "number_layers_found_5km": None,
            # Variables derived for 5km data
            # Also included if a combination of 5 and 1km data are used for RESOLUTION=1
            'detection_height_5km': None,
            'total_optical_depth_5km': None,
            "feature_optical_depth_532_top_layer_5km": None,
            'cfc_single_shots_1km_from_5km_file': None,
            "average_cloud_top_pressure_single_shots": None,
            "average_cloud_top_pressure_single_shots_5km": None,
            "average_cloud_top_single_shots": None,
            "average_cloud_top_single_shots_5km": None,
            "average_cloud_base_single_shots": None,
            "average_cloud_base_single_shots_5km": None,
            "single_shot_data": None,
            # Variables derived from 5km file to 1kmresolution_
            'cfc_single_shots_1km_from_5km_file': None,

            # From cloudsat:
            'cal_modis_cflag': None,
            'cloudsat_index': None,
        }


class CloudsatObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)
        self.all_arrays = {
            'clsat_max_height': None,
            'longitude': None,
            'latitude': None,
            'imager_linnum': None,
            'imager_pixnum': None,
            'cloud_fraction': None,
            'validation_height': None,
            'validation_height_base': None,
            'elevation': None,
            'sec_1970': None,
            'CPR_Cloud_mask': None,
            'MODIS_Cloud_Fraction': None,
            'MODIS_cloud_flag': None,
            'Height': None,
            'LO_RVOD_liquid_water_path': None,
            'IO_RVOD_ice_water_path': None,
            'LO_RO_liquid_water_path': None,
            'IO_RO_ice_water_path': None,
            #'liq_water_path': None,  # kg!/m2 R05
            #'ice_water_path': None,  # kg!/m2 R05
            'RVOD_liq_water_path': None,  # g/m2 R04
            'RVOD_ice_water_path': None,  # g/m2 R04
            'RO_liq_water_path': None,  # g/m2 R05
            'RO_ice_water_path': None,  # g/m2 R05
            'precip_liq_water_path_gm2': None,  # g/m2
            'cloud_liq_water_path_gm2': None,  # g/m2 
            'precip_ice_water_path_gm2': None,  # g/m2 
            'cloud_ice_water_path_gm2': None,  # g/m2
            'liq_water_path_gm2': None,  # g/m2
            'ice_water_path_gm2': None,  # g/m2
            'precip_liq_water_path': None,  # kg/m2 R05
            'cloud_liq_water_path': None,  # kg/m2 R05
            'precip_ice_water_path': None,  # kg/m2 R05
            'cloud_ice_water_path': None,  # kg/m2 R05
            'liq_water_path': None,  # kg/m2 R05
            'ice_water_path': None,  # kg/m2 R05
            #'ice_water_content': None,
            #'liq_water_content': None,
            "RVOD_CWC_status": None,
            "RO_CWC_status": None,
            'Phase': None,
            'Profile_time': None,
            'TAI_start': None,
            # From calipso
            'calipso_layer_base_altitude': None,
            'calipso_layer_top_altitude': None,
            'calipso_feature_classification_flags': None
        }


class IssObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)
        self.all_arrays = {
            # Name: Iss name .lower()
            'longitude': None,
            'latitude': None,
            # Derived:
            'imager_linnum': None,
            'imager_pixnum': None,
            'sec_1970': None,
            'elevation': None,
            'cloud_fraction': None,
            'validation_height': None,
            'total_optical_depth_5km': None,
            # Used
            'cloud_phase_fore_fov': None,
            'feature_type_fore_fov': None,
            'extinction_qc_flag_1064_fore_fov': None,
            'layer_top_altitude_fore_fov': None,
            'sky_condition_fore_fov': None,
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
            'lwp': None,
            'pixel_status': None,
            'quality': None,
            'surface_type': None}


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
            'nh': None,
            'cl': None,
            'cm': None,
            'ch': None,
            'vvv': None,
            'ww': None,
            'temp': None,
            'dtemp': None,
            'sec_1970': None,
            'pressure': None}


class TruthImagerTrackObject:
    def __init__(self, truth='calipso'):
        self.imager = ExtractedImagerObject()
        self.modis_lvl2 = ModisObject()
        if truth in 'calipso':
            self.calipso = CalipsoObject()
            self.calipso_aerosol = CalipsoObject()
        elif truth in 'cloudsat':
            self.cloudsat = CloudsatObject()
        elif truth in 'amsr':
            self.amsr = AmsrObject()
        elif truth in 'synop':
            self.synop = SynopObject()
        elif truth in 'mora':
            self.mora = MoraObject()
        elif truth in 'iss':
            self.iss = IssObject()
        self.extra = ExtraObject()
        self.diff_sec_1970 = None
        self.truth_sat = truth
        self.imager_instrument = 'imager'

    def make_nsidc_surface_type_texture(self, kernel_sz=51):
        """Derive the stdv of the ice dataset"""

        if self.calipso.all_arrays['nsidc_surface_type'] is not None:
            self.calipso.all_arrays['nsidc_surface_type_texture'] = sliding_std(
                self.calipso.all_arrays['nsidc_surface_type'], kernel_sz)

    def __add__(self, other):
        """Concatenating two objects together"""
        for object_name in ['imager', 'calipso', 'calipso_aerosol', 'amsr',
                            'cloudsat', 'iss', 'mora', 'synop', 'modis_lvl2', 'modis', 'extra']:
            if hasattr(self, object_name):
                setattr(self, object_name,
                        getattr(self, object_name) +
                        getattr(other, object_name))
        try:
            self.diff_sec_1970 = np.concatenate([self.diff_sec_1970,
                                                 other.diff_sec_1970])
        except ValueError:
            # print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = other.diff_sec_1970
        return self

    def extract_elements(self, idx=None, starti=None, endi=None):
        for object_name in ['imager', 'calipso', 'calipso_aerosol', 'amsr',
                            'cloudsat', 'iss', 'mora', 'synop', 'modis', 'extra']:
            if hasattr(self, object_name):
                obj = getattr(self, object_name)
                setattr(self, object_name, obj.extract_elements(idx=idx, starti=starti, endi=endi))
        try:
            if idx is not None:
                self.diff_sec_1970 = self.diff_sec_1970[idx]
            else:
                self.diff_sec_1970 = self.diff_sec_1970[starti:endi]
        except ValueError:
            # print "Don't concatenate member diff_sec_1970... " + str(e)
            self.diff_sec_1970 = self.diff_sec_1970
        return self


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
    if 'oca' in h5file.keys():
        h5_groups.append(h5file['/oca'])
        data_objects.append(retv.imager)
    if 'patmosx' in h5file.keys():
        h5_groups.append(h5file['/patmosx'])
        data_objects.append(retv.imager)
    if 'modis_lvl2' in h5file.keys():
        h5_groups.append(h5file['/modis_lvl2'])
        data_objects.append(retv.modis_lvl2)
    if 'cloudsat' in h5file.keys():
        h5_groups.append(h5file['/cloudsat'])
        data_objects.append(retv.cloudsat)
    if 'iss' in h5file.keys():
        h5_groups.append(h5file['/iss'])
        data_objects.append(retv.iss)
    if 'amsr' in h5file.keys():
        h5_groups.append(h5file['/amsr'])
        data_objects.append(retv.amsr)
    if 'mora' in h5file.keys():
        h5_groups.append(h5file['/mora'])
        data_objects.append(retv.mora)
    if 'synop' in h5file.keys():
        h5_groups.append(h5file['/synop'])
        data_objects.append(retv.synop)
    if 'cmaprob_cots' in h5file:
        h5_groups.append(h5file['/cmaprob_cots'])
        data_objects.append(retv.extra)
    if 'extra' in h5file:
        h5_groups.append(h5file['/extra'])
        data_objects.append(retv.extra)
    return (h5_groups, data_objects)


def read_truth_imager_match_obj(filename, truth='calipso',
                                read_all=True,
                                read_var=[],
                                skip_var=[]):
    retv = TruthImagerTrackObject(truth=truth)
    h5file = h5py.File(filename, 'r')
    (h5_groups, data_objects) = get_stuff_to_read_from_a_reshaped_file(h5file, retv)
    for group, data_obj in zip(h5_groups, data_objects):
        for dataset in group.keys():
            if dataset in skip_var:
                continue
            if (read_all or dataset in read_var or
                    (len(read_var) == 0 and dataset.data_obj.all_arrays.keys())):
                atrain_match_name = dataset
                if atrain_match_name in ["snow_ice_surface_type"]:
                    atrain_match_name = "nsidc_surface_type"
                setattr(data_obj, atrain_match_name, group[dataset][...])
    retv.diff_sec_1970 = h5file['diff_sec_1970'][...]
    h5file.close()
    return retv


def read_files(files, truth='calipso', read_all=True, read_var=[], skip_var=[]):
    my_files = files.copy()
    tObj = read_truth_imager_match_obj(my_files.pop(), truth=truth, read_all=read_all, read_var=read_var, skip_var=skip_var)
    if len(my_files) > 0:
        for filename in my_files:
            tObj += read_truth_imager_match_obj(filename, truth=truth, read_all=read_all, read_var=read_var, skip_var=skip_var)
    return tObj


# write matchup files


def write_truth_imager_match_obj(filename, match_obj, SETTINGS=None, imager_obj_name='pps'):
    """Write *match_obj* to *filename*."""
    datasets = {'diff_sec_1970': match_obj.diff_sec_1970}
    groups = {imager_obj_name: match_obj.imager.all_arrays}
    imager_attrs = {'imager_instrument': match_obj.imager_instrument}
    groups_attrs = {imager_obj_name: imager_attrs}
    for name in ['calipso', 'calipso_aerosol', 'iss', 'modis_lvl2',
                 'amsr', 'synop', 'mora', 'cloudsat', 'extra']:
        if hasattr(match_obj, name):
            groups[name] = getattr(match_obj, name).all_arrays
    write_match_objects(filename, datasets, groups, groups_attrs, SETTINGS=SETTINGS)
    return 1


def sliding_std(x, size=5):
    """derive a sliding standard deviation of a data array"""
    from scipy.ndimage.filters import uniform_filter
    c1 = uniform_filter(x.astype('float'), size=size)
    c2 = uniform_filter(x.astype('float')*x.astype('float'), size=size)
    return abs(c2 - c1 * c1)**.5


the_used_variables = [
    'longitude',
    'latitude',
    'sec_1970',
    'imager_linnum',
    'imager_pixnum',
    'imager_linnum_nneigh',
    'imager_pixnum_nneigh',
    'sec_1970',
    'elevation',
    # MODIS LVL2
    'height',
    'temperature',
    'pressure',
    'cloud_emissivity',
    'cloud_phase',
    'lwp',
    # AMSR
    'lwp',
    'imager_amsr_dist',
    'pixel_status',
    'quality',
    'surface_type',
    # MORA
    'cloud_base_height',
    # Cloudsat
    'cloud_fraction',
    'validation_height',
    'LO_RVOD_liquid_water_path',
    'IO_RVOD_ice_water_path',
    'LO_RO_liquid_water_path',
    'IO_RO_ice_water_path',
    'RVOD_liq_water_path',  # g/m2 R04
    'RVOD_ice_water_path',  # g/m2 R04
    'RO_liq_water_path',  # g/m2 R05
    'RO_ice_water_path',  # g/m2 R05
    'precip_liq_water_path_gm2',  # g/m2
    'cloud_liq_water_path_gm2',  # g/m2 
    'precip_ice_water_path_gm2',  # g/m2 
    'cloud_ice_water_path_gm2',  # g/m2
    'liq_water_path_gm2',  # g/m2
    'ice_water_path_gm2',  # g/m2
    'precip_liq_water_path',  # kg/m2 R05
    'cloud_liq_water_path',  # kg/m2 R05
    'precip_ice_water_path',  # kg/m2 R05
    'cloud_ice_water_path',  # kg/m2 R05
    'liq_water_path',  # kg/m2 R05
    'ice_water_path',  # kg/m2 R05
    #'ice_water_content',
    #'liq_water_content',
    "RVOD_CWC_status",
    # CALIPSO write do not combine
    'cal_modis_cflag',
    'cloudsat_index',

    # CALIPSO only (ISS?)
    'minimum_laser_energy_532',

    "average_cloud_top_pressure_single_shots",
    "average_cloud_top_pressure_single_shots_1km",
    "average_cloud_top_single_shots",
    "average_cloud_top_single_shots_1km",
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
    'single_shot_data',
    'cfc_single_shots_1km_from_5km_file',
    'feature_optical_depth_532_top_layer_5km',
    'feature_optical_depth_532_5km',
    'total_optical_depth_5km',
    'detection_height_5km',
    'column_optical_depth_cloud_532',
    'column_optical_depth_cloud_uncertainty_532',
    'column_optical_depth_tropospheric_aerosols_532_5km',
    'column_optical_depth_tropospheric_aerosols_532',
    'column_optical_depth_aerosols_532_5km',
    'column_optical_depth_aerosols_532',
    "layer_top_altitude_5km",
    "layer_top_pressure_5km",
    'number_cloudy_single_shots',
    "average_cloud_base_single_shots_5km",
    "average_cloud_top_pressure_single_shots_5km",
    "average_cloud_top_single_shots_5km",
    # CLOUDSAT only
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
    match_calipso = read_truth_imager_match_obj(TESTFILE)
    match_calipso2 = read_truth_imager_match_obj(TESTFILE2)

    match_calipso = match_calipso + match_calipso2
