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
import numpy as np
import copy
import h5py
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import get_neighbour_info
from pyresample.kd_tree import get_sample_from_neighbour_info
from pyresample.kd_tree import resample_nearest
from pyresample import geometry
from pyresample.geometry import AreaDefinition
import cartopy.crs as ccrs
from atrain_match.config import AREA_CONFIG_FILE_PLOTS_ON_AREA
from pyresample import load_area
# import pyresample as pr
import os
import matplotlib
from scipy.interpolate import griddata
from scipy import ndimage
import matplotlib.pyplot as plt
# matplotlib.use("TkAgg")
from atrain_match.matchobject_io import DataObject
from atrain_match.utils.get_flag_info import (get_calipso_low_clouds,
                                              get_calipso_cad_score,
                                              get_calipso_clouds_of_type_i,
                                              get_calipso_high_clouds)
matplotlib.rcParams.update({'font.size': 30})

COMPRESS_LVL = 6
cots = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
        0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 2.0, 3.0, 4.0, 5.0]
cots_mid = []
for ind, filtcot in enumerate(cots):
    if ind == len(cots)-1:
        filtcot_upper = 999.999
    else:
        filtcot_upper = cots[ind+1]
    cots_mid += ["{:0.3f}".format(0.5 * (filtcot + filtcot_upper))]



class MatchImagerCalipso(DataObject):
    def __init__(self):
        DataObject.__init__(self)
        self.all_arrays = {
            'detected_clouds': None,
            'undetected_clouds': None,
            'false_clouds': None,
            'detected_clear': None,
            'new_detected_clouds': None,
            'new_false_clouds': None,
            'lats': None,
            'lons': None}

    def get_some_info_from_caobj(self, match_calipso, PROCES_FOR_ART=False, PROCES_FOR_PRESSENTATIONS=False):
        self.set_false_and_missed_cloudy_and_clear(match_calipso=match_calipso, PROCES_FOR_ART=PROCES_FOR_ART,
                                                   PROCES_FOR_PRESSENTATIONS=PROCES_FOR_PRESSENTATIONS)
        self.get_detection_for_cots(match_calipso)
        self.set_r13_extratest(match_calipso=match_calipso)
        self.get_thr_offset(match_calipso=match_calipso)
        self.get_lapse_rate(match_calipso=match_calipso)
        self.get_ctth_bias(match_calipso=match_calipso)
        self.get_ctth_bias_low(match_calipso=match_calipso)
        self.get_ctth_bias_high(match_calipso=match_calipso)
        self.get_ctth_bias_low_temperature(match_calipso=match_calipso)
        self.height_bias_type = {}
        self.detected_height_type = {}
        for cc_type in range(8):
            self.get_ctth_bias_type(match_calipso=match_calipso, calipso_cloudtype=cc_type)

    def set_false_and_missed_cloudy_and_clear(self, match_calipso, PROCES_FOR_ART, PROCES_FOR_PRESSENTATIONS):
        lat = match_calipso.imager.all_arrays['latitude']
        lon = match_calipso.imager.all_arrays['longitude']
        if match_calipso.imager.all_arrays['cloudmask'] is not None:
            cloudy_imager = np.logical_or(match_calipso.imager.all_arrays['cloudmask'] == 1,
                                          match_calipso.imager.all_arrays['cloudmask'] == 2)
            clear_imager = np.logical_or(match_calipso.imager.all_arrays['cloudmask'] == 0,
                                         match_calipso.imager.all_arrays['cloudmask'] == 3)
        else:
            cloudy_imager = np.logical_and(match_calipso.imager.all_arrays['cloudtype'] > 4,
                                           match_calipso.imager.all_arrays['cloudtype'] < 21)
            clear_imager = np.logical_and(match_calipso.imager.all_arrays['cloudtype'] > 0,
                                          match_calipso.imager.all_arrays['cloudtype'] < 5)

        # cloudy_imager = np.logical_and(match_calipso.imager.all_arrays['cloudtype'] > 4,
        #                             match_calipso.imager.all_arrays['cloudtype'] < 21)
        # clear_imager = np.logical_and(match_calipso.imager.all_arrays['cloudtype'] > 0,
        #                            match_calipso.imager.all_arrays['cloudtype'] < 5)
        nlay = np.where(match_calipso.calipso.all_arrays['number_layers_found'] > 0, 1, 0)
        meancl = ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
        if self.cc_method == 'BASICOLD' and self.isGAC:
            # BASIC before 20181009
            # Buggy for 1km data as cloud_fraction was calcualted with mean filter
            # From neighbours.
            calipso_cloudy = match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5
            calipso_clear = np.not_equal(calipso_cloudy, True)
        elif self.cc_method == 'BASIC':
            calipso_cloudy = nlay > 0  # With ALL 300m data for 5km! Updated
            calipso_clear = nlay == 0  # np.not_equal(calipso_cloudy, True)
        elif self.cc_method == 'COTFILT':
            filtcot = self.cotfilt_value
            # Let's try to correct for the "artificial" jump in CFC at cot value 1.0 due to
            # added 1 km CALIPSO clouds with cot set to 1.0. These clouds should be moved to
            # the other extreme at cot = 5.0 (they are definitely bright clouds!)
            is_artificial_cloudy = np.logical_and(
                match_calipso.calipso.all_arrays['cloud_fraction'] < 1.0,
                match_calipso.calipso.all_arrays['total_optical_depth_5km'] == 1.0)
            is_artificial_cloudy = np.logical_and(
                is_artificial_cloudy,
                match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5)
            if filtcot > 1.0:
                # ("Artificial" jump in CFC due to added 1 km CALIPSO clouds with tau=1.0 ==>
                # remove jump and add these cases to category 5.0!
                # Note that 300m data have total_optical_thickness=1.0, should ahve been 5.0
                # We must add them here!
                calipso_cloudy = np.logical_and(
                    match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5,
                    match_calipso.calipso.all_arrays['total_optical_depth_5km'] >= filtcot)
                calipso_cloudy = np.logical_or(calipso_cloudy, is_artificial_cloudy)
                calipso_clear = np.not_equal(calipso_cloudy, True)
            else:
                # Note that 300m data have total_optical_thickness=1.0, should ahve been 5.0
                # They are still here as 1.0 is larger than filtcot!
                calipso_cloudy = np.logical_and(
                    match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5,
                    # match_calipso.calipso.all_arrays['total_optical_depth_5km']>filtcot)  #  serious error!!!!
                    match_calipso.calipso.all_arrays['total_optical_depth_5km'] >= filtcot)
                # match_calipso.calipso.all_arrays['total_optical_depth_5km']>0.225)
                # Changed COT threshold from 0.15 to 0.225 when going to CALIPSO version 4.10/KG 2017-04-20
                calipso_clear = np.not_equal(calipso_cloudy, True)
        elif self.cc_method == 'KG' and self.isGAC:
            calipso_cloudy = np.logical_and(
                match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5,
                match_calipso.calipso.all_arrays['total_optical_depth_5km'] > 0.15)
            calipso_clear = np.not_equal(calipso_cloudy, True)
        elif self.cc_method == 'KG':
            calipso_cloudy = nlay > 0
            calipso_clear = np.not_equal(calipso_cloudy, True)
        elif self.cc_method == 'Nina' and self.isGAC:
            calipso_cloudy = np.logical_and(
                match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5,
                match_calipso.calipso.all_arrays['total_optical_depth_5km'] > 0.15)
            # exclude pixels that might be cloud contaminated
            calipso_clear = np.logical_and(nlay == 0, meancl < 0.01)
            calipso_clear = np.logical_and(
                calipso_clear,
                match_calipso.calipso.all_arrays['total_optical_depth_5km'] < 0)
        elif self.cc_method == 'Nina':
            # calipso_cloudy = np.logical_or(
            #   match_calipso.calipso.all_arrays['total_optical_depth_5km']>0.15,
            #   np.logical_and(match_calipso.calipso.all_arrays['total_optical_depth_5km']<0,
            #                  nlay>0))
            calipso_cloudy = nlay > 0
            # exclude pixels that might be cloud contaminated
            calipso_clear = np.logical_and(nlay == 0, meancl < 0.01)
            # and not or! dangerous for broken clouds!
            calipso_clear = np.logical_and(
                calipso_clear,
                match_calipso.calipso.all_arrays['total_optical_depth_5km'] < 0)
        elif self.cc_method == 'Abhay' and self.isGAC:
            # Excude pixels with low_cad_score, compare with mode BASIC
            conf_medium_or_high, conf_no, conf_low = get_calipso_cad_score(match_calipso)
            calipso_cloudy = np.logical_or(
                # High confidence cloudy
                np.logical_and(nlay > 0, conf_medium_or_high),
                # clouds clouds from 300m data
                np.logical_and(match_calipso.calipso.all_arrays['cloud_fraction'] > 0.0,
                               match_calipso.calipso.all_arrays['cloud_fraction'] < 1.0))
            calipso_clear = nlay == 0
        elif self.cc_method == 'Abhay':
            # Excude pixels with low_cad_score
            conf_medium_or_high, conf_no, conf_low = get_calipso_cad_score(match_calipso)
            calipso_cloudy = np.logical_and(nlay > 0, conf_medium_or_high)
            calipso_clear = nlay == 0
        use_ok_lat_lon = np.logical_and(np.logical_and(lat >= -90, lat <= 90),
                                        np.logical_and(lon >= -180, lat <= 180))
        cloudy_imager = np.logical_and(cloudy_imager, use_ok_lat_lon)
        clear_imager = np.logical_and(clear_imager, use_ok_lat_lon)
        sunz = match_calipso.imager.all_arrays['sunz']
        if self.filter_method == "satz":
            satz = match_calipso.imager.all_arrays['satz']
            cloudy_imager = np.logical_and(cloudy_imager, satz <= 30)
            clear_imager = np.logical_and(clear_imager, satz <= 30)
        if self.DNT in ["day"]:
            cloudy_imager = np.logical_and(cloudy_imager, sunz <= 80)
            clear_imager = np.logical_and(clear_imager, sunz <= 80)
        if self.DNT in ["night"]:
            cloudy_imager = np.logical_and(cloudy_imager, sunz >= 95)
            clear_imager = np.logical_and(clear_imager, sunz >= 95)
        if self.DNT in ["twilight"]:
            cloudy_imager = np.logical_and(cloudy_imager, sunz > 80)
            clear_imager = np.logical_and(clear_imager, sunz > 80)
            cloudy_imager = np.logical_and(cloudy_imager, sunz < 95)
            clear_imager = np.logical_and(clear_imager, sunz < 95)
        if self.DNT in ["all"]:
            pass
        undetected_clouds = np.logical_and(calipso_cloudy, clear_imager)
        false_clouds = np.logical_and(calipso_clear, cloudy_imager)
        detected_clouds = np.logical_and(calipso_cloudy, cloudy_imager)
        detected_clear = np.logical_and(calipso_clear, clear_imager)
        use = np.logical_or(np.logical_or(detected_clouds, detected_clear),
                            np.logical_or(false_clouds, undetected_clouds))
        detected_height = np.logical_and(detected_clouds,
                                         match_calipso.imager.all_arrays['ctth_height'] > -9)
        detected_height = np.logical_and(detected_height,
                                         match_calipso.imager.all_arrays['ctth_height'] < 45000)

        # settings for article eos_modis:
        if PROCES_FOR_ART:
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['warmest_t12'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['coldest_t12'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['psur'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['surftemp'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t950'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t850'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t700'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t500'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t250'] > -9)
            detected_height = np.logical_and(
                detected_height, match_calipso.imager.all_arrays['text_t11'] > -9)  # without this 1793146 pixels
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['text_t11t12'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['warmest_t11'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['coldest_t11'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['modis_27'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['modis_28'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['modis_33'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['text_t37'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['warmest_t37'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['coldest_t37'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt11micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt12micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt37micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt86micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['ciwv'] > -9)
            detected_height = np.logical_and(
                detected_height, match_calipso.imager.all_arrays["ctthold_height"] > -9)  # only included in art
            detected_height = np.logical_and(
                detected_height, match_calipso.imager.all_arrays["ctthnnant_height"] > -9)  # only included in art
            detected_height = np.logical_and(detected_height, match_calipso.modis.all_arrays["height"] > -9)

        if PROCES_FOR_PRESSENTATIONS:
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['warmest_t12'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['coldest_t12'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['psur'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['surftemp'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t950'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t850'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t700'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t500'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['t250'] > -9)
            detected_height = np.logical_and(
                detected_height, match_calipso.imager.all_arrays['text_t11'] > -9)  # without this 1793146 pixels
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['text_t11t12'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['warmest_t11'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['coldest_t11'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['text_t37'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['warmest_t37'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['coldest_t37'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt11micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt12micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt37micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['bt86micron'] > -1)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays['ciwv'] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.modis.all_arrays["height"] > -9)
            detected_height = np.logical_and(detected_height, match_calipso.imager.all_arrays["ctth_height"] > -9)

        detected_height = np.logical_and(detected_height,
                                         match_calipso.calipso.all_arrays['layer_top_altitude'][:, 0] > -1)
        detected_temperature = np.logical_and(detected_clouds,
                                              match_calipso.imager.all_arrays['ctth_temperature'] > -9)

        self.false_clouds = false_clouds[use]
        self.detected_clouds = detected_clouds[use]
        self.undetected_clouds = undetected_clouds[use]
        self.detected_clear = detected_clear[use]
        self.latitude = match_calipso.imager.latitude[use]
        self.longitude = match_calipso.imager.longitude[use]
        self.use = use
        self.detected_height = detected_height[use]
        self.detected_temperature = detected_temperature[use]

    def get_detection_for_cots(self, match_calipso):
        # ccm method sould be  BASIC/BASICOLD/COTFILT (not one removing thin clouds)
        if match_calipso.imager.all_arrays['cloudmask'] is not None:
            cloudy_imager = np.logical_or(match_calipso.imager.all_arrays['cloudmask'] == 1,
                                          match_calipso.imager.all_arrays['cloudmask'] == 2)
            clear_imager = np.logical_or(match_calipso.imager.all_arrays['cloudmask'] == 0,
                                         match_calipso.imager.all_arrays['cloudmask'] == 3)
        else:
            cloudy_imager = np.logical_and(match_calipso.imager.all_arrays['cloudtype'] > 4,
                                           match_calipso.imager.all_arrays['cloudtype'] < 21)
            clear_imager = np.logical_and(match_calipso.imager.all_arrays['cloudtype'] > 0,
                                          match_calipso.imager.all_arrays['cloudtype'] < 5)
        is_artificial_cloudy = np.logical_and(
            match_calipso.calipso.all_arrays['cloud_fraction'] < 1.0,
            match_calipso.calipso.all_arrays['total_optical_depth_5km'] == 1.0)
        is_artificial_cloudy = np.logical_and(
            is_artificial_cloudy,
            match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5)
        setattr(self, "detected_clouds_filtcot", {})
        setattr(self, "undetected_clouds_filtcot", {})
        self.detected_clouds_filtcot = {}
        self.undetected_clouds_filtcot = {}
        calipso_cloudy_all = match_calipso.calipso.all_arrays['cloud_fraction'] > 0.5
        for ind, filtcot in enumerate(cots):
            if ind == len(cots)-1:
                filtcot_upper = 99999999.999
            else:
                filtcot_upper = cots[ind+1]
            calipso_cloudy = np.logical_and(
                match_calipso.calipso.all_arrays['total_optical_depth_5km'] < filtcot_upper,
                match_calipso.calipso.all_arrays['total_optical_depth_5km'] >= filtcot)
            calipso_cloudy = np.logical_and(calipso_cloudy, calipso_cloudy_all)
            if filtcot >= 5.0:
                # Artifically cloudy should have optical thickness 5.0
                calipso_cloudy = np.logical_or(calipso_cloudy, is_artificial_cloudy)
            else:
                calipso_cloudy = np.logical_and(calipso_cloudy, ~is_artificial_cloudy)
            calipso_clear = np.not_equal(calipso_cloudy, True)
            detected_clouds = np.logical_and(calipso_cloudy, cloudy_imager)
            undetected_clouds = np.logical_and(calipso_cloudy, clear_imager)
            # self.false_clouds = false_clouds[use]
            filtcot_mid = cots_mid[ind]
            self.detected_clouds_filtcot[filtcot_mid] = detected_clouds[self.use]
            self.undetected_clouds_filtcot[filtcot_mid] = undetected_clouds[self.use]
            # self.detected_clear = detected_clear[use]

    def set_r13_extratest(self, match_calipso):
        if ('r13micron' not in match_calipso.imager.all_arrays.keys() or
                match_calipso.imager.all_arrays['r13micron'] is None):
            self.new_false_clouds = np.zeros(self.false_clouds.shape)
            self.new_detected_clouds = np.zeros(self.false_clouds.shape)
            return
        r13 = match_calipso.imager.all_arrays['r13micron']
        sunz = match_calipso.imager.all_arrays['sunz']
        sunz_cos = sunz.copy()
        sunz_cos[sunz > 87] = 87
        r13[sunz < 90] = r13[sunz < 90]/np.cos(np.radians(sunz_cos[sunz < 90]))
        isCloud_r13 = np.logical_and(r13 > 2.0, match_calipso.imager.all_arrays['ciwv'] > 3)
        new_detected_clouds = np.logical_and(self.detected_clouds,
                                             isCloud_r13[self.use])
        new_false_clouds = np.logical_and(self.detected_clear,
                                          isCloud_r13[self.use])
        self.new_false_clouds = new_false_clouds
        self.new_detected_clouds = new_detected_clouds

    def get_thr_offset(self, match_calipso):
        if ('surftemp' not in match_calipso.imager.all_arrays.keys() or
                match_calipso.imager.all_arrays['surftemp'] is None):
            self.t11ts_offset = np.zeros(self.false_clouds.shape)
            self.t11t12_offset = np.zeros(self.false_clouds.shape)
            self.t11t37_offset = np.zeros(self.false_clouds.shape)
            self.t37t12_offset = np.zeros(self.false_clouds.shape)
            return
        t11ts_offset = (match_calipso.imager.all_arrays['bt11micron'] -
                        match_calipso.imager.all_arrays['surftemp'] -
                        match_calipso.imager.all_arrays['thr_t11ts_inv'])
        t11t12_offset = (match_calipso.imager.all_arrays['bt11micron'] -
                         match_calipso.imager.all_arrays['bt12micron'] -
                         match_calipso.imager.all_arrays['thr_t11t12'])
        t37t12_offset = (match_calipso.imager.all_arrays['bt37micron'] -
                         match_calipso.imager.all_arrays['bt12micron'] -
                         match_calipso.imager.all_arrays['thr_t37t12'])
        t11t37_offset = (match_calipso.imager.all_arrays['bt11micron'] -
                         match_calipso.imager.all_arrays['bt37micron'] -
                         match_calipso.imager.all_arrays['thr_t11t37'])
        t11ts_offset = t11ts_offset[self.use]
        t11ts_offset[self.detected_clouds] = 99999
        t11ts_offset[self.undetected_clouds] = 99999
        self.t11ts_offset = t11ts_offset
        t11t12_offset = t11t12_offset[self.use]
        t11t12_offset[self.detected_clouds] = -99999
        t11t12_offset[self.undetected_clouds] = -99999
        self.t11t12_offset = t11t12_offset
        t37t12_offset = t37t12_offset[self.use]
        t37t12_offset[self.detected_clouds] = -99999
        t37t12_offset[self.undetected_clouds] = -99999
        self.t37t12_offset = t37t12_offset
        t11t37_offset = t11t37_offset[self.use]
        t11t37_offset[self.detected_clouds] = -99999
        t11t37_offset[self.undetected_clouds] = -99999
        self.t11t37_offset = t11t37_offset

    def get_lapse_rate(self, match_calipso):
        if ('surftemp' not in match_calipso.imager.all_arrays.keys() or
                match_calipso.imager.all_arrays['surftemp'] is None):
            self.lapse_rate = np.zeros(self.false_clouds.shape)
            return
        low_clouds = get_calipso_low_clouds(match_calipso)
        delta_h = (match_calipso.calipso.all_arrays['layer_top_altitude'][:, 0]
                   - 0.001*match_calipso.calipso.all_arrays['elevation'][:])
        delta_t = (273.15 + match_calipso.calipso.all_arrays['layer_top_temperature']
                   [:, 0] - match_calipso.imager.all_arrays['surftemp'][:])
        lapse_rate = delta_t/delta_h
        lapse_rate[match_calipso.calipso.all_arrays['layer_top_temperature'][:, 0] < -500] = 0
        lapse_rate[match_calipso.calipso.all_arrays['layer_top_altitude'][:, 0] > 35.0] = 0
        lapse_rate[~low_clouds] = 0.0
        self.lapse_rate = lapse_rate[self.use]

    def get_ctth_bias(self, match_calipso):
        height_c = (1000*(
            match_calipso.calipso.all_arrays['layer_top_altitude'][self.use, 0]) -
            match_calipso.calipso.all_arrays['elevation'][self.use])
        if "cci" in self.satellites:
            height_c = 1000*match_calipso.calipso.all_arrays[
                'layer_top_altitude'][self.use, 0]
        if "modis_lvl2" in self.satellites:
            height_c = 1000*match_calipso.calipso.all_arrays[
                'layer_top_altitude'][self.use, 0]

        height_pps = match_calipso.imager.all_arrays['ctth_height'][self.use]
        delta_h = height_pps - height_c
        self.height_bias = delta_h
        self.height_bias[~self.detected_height] = 0

        self.detected_height_both = np.where(self.detected_height, False, self.detected_height)
        self.height_mae_diff = 0*self.height_bias
        if match_calipso.modis.all_arrays["height"] is not None:
            height_modis = match_calipso.modis.all_arrays["height"][self.use]
            delta_h_modis = height_modis - 1000*match_calipso.calipso.all_arrays[
                'layer_top_altitude'][self.use, 0]
            self.detected_height_both = np.logical_and(self.detected_height, height_modis > 0)
            mae_pps = np.abs(self.height_bias.copy())
            mae_modis = np.abs(delta_h_modis)
            diff_mae = mae_modis - mae_pps
            diff_mae[~self.detected_height_both] = 0

            self.height_mae_diff = diff_mae

        try:
            # tsur = match_calipso.imager.all_arrays['surftemp']
            tsur = match_calipso.imager.all_arrays['segment_nwp_temp'][:, 0]
            # temperature_pps = match_calipso.imager.all_arrays['ctth_temperature']
            temperature_pps = match_calipso.imager.all_arrays['bt12micron']
            temp_diff = temperature_pps - tsur
            rate_neg = -1.0/6.5
            rate_pos = + 1.0/3.0
            new_pps_h = rate_neg*temp_diff*1000
            new_pps_h[temp_diff > 0] = rate_pos*temp_diff[temp_diff > 0]*1000
            new_pps_h[new_pps_h < 100] = 100
            keep = (
                (tsur - match_calipso.imager.all_arrays['ctth_temperature']) /
                (tsur - match_calipso.imager.all_arrays['ttro'])) > 0.33
            keep = new_pps_h > 3000
            new_pps_h[keep] = match_calipso.imager.all_arrays['ctth_height'][keep]
            self.lapse_bias = new_pps_h[self.use] - height_c
            self.lapse_bias[~self.detected_height] = 0
            self.lapse_bias[temperature_pps < 0] = 0
        except (TypeError, KeyError):
            # raise
            self.lapse_bias = 0 * height_c

    def get_ctth_bias_low(self, match_calipso):
        low_clouds = get_calipso_low_clouds(match_calipso)
        detected_low = np.logical_and(self.detected_height,
                                      low_clouds[self.use])
        delta_h = self.height_bias.copy()
        delta_h[~detected_low] = 0
        self.height_bias_low = delta_h
        self.detected_height_low = detected_low
        delta_h = self.lapse_bias.copy()
        delta_h[~detected_low] = 0
        self.lapse_bias_low = delta_h

    def get_ctth_bias_high(self, match_calipso):
        high_clouds = get_calipso_high_clouds(match_calipso)
        detected_high = np.logical_and(self.detected_height, high_clouds[self.use])
        delta_h = self.height_bias.copy()
        delta_h[~detected_high] = 0
        self.height_bias_high = delta_h.copy()
        self.detected_height_high = detected_high
        delta_h = self.lapse_bias.copy()
        delta_h[~detected_high] = 0
        self.lapse_bias_high = delta_h

    def get_ctth_bias_low_temperature(self, match_calipso):
        low_clouds = get_calipso_low_clouds(match_calipso)
        detected_low = np.logical_and(self.detected_height, low_clouds[self.use])
        temperature_pps = match_calipso.imager.all_arrays['ctth_temperature'][self.use]
        try:
            temperature_c = 273.15 + match_calipso.calipso.all_arrays['midlayer_temperature'][self.use, 0]

            delta_t = temperature_pps - temperature_c
            delta_t[~detected_low] = 0
            delta_t[temperature_c < 0] = 0
            delta_t[temperature_pps < 0] = 0
        except (TypeError, KeyError):
            delta_t = 0*temperature_pps
        try:
            temperature_pps = match_calipso.imager.all_arrays['bt11micron'][self.use]
            delta_t_t11 = temperature_pps - temperature_c
            delta_t_t11[~detected_low] = 0
            delta_t_t11[temperature_c < 0] = 0
            delta_t_t11[temperature_pps < 0] = 0
        except (TypeError, KeyError):
            delta_t_t11 = 0*delta_t
        self.temperature_bias_low = delta_t
        # self.temperature_bias_low[~detected_low]=0
        self.temperature_bias_low_t11 = delta_t_t11

    def get_ctth_bias_type(self, match_calipso, calipso_cloudtype=0):
        wanted_clouds = get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=calipso_cloudtype)
        detected_typei = np.logical_and(self.detected_height, wanted_clouds[self.use])
        delta_h = self.height_bias.copy()
        delta_h[~detected_typei] = 0
        self.height_bias_type[calipso_cloudtype] = delta_h
        self.detected_height_type[calipso_cloudtype] = detected_typei




class StatsOnFibonacciLattice(DataObject):
    def __init__(self):
        DataObject.__init__(self)
        self.all_arrays = {
            'definition': None,
            'lats': None,
            'lons': None,
            'num_false_clouds': None,
            'Min_lapse_rate': None,
            'num_detected_clouds': None,
            'num_new_false_clouds': None,
            'num_new_detected_clouds': None,
            'num_undetected_clouds': None,
            'num_detected_clear': None,
            'num_clear': None,
            'num_clouds': None,
            'N': None,
            'Kuipers': None}

    def write_fibonacci_grid(self, filename=None):
        """Write match objects to HDF5 file *filename*."""
        # pplot_obj.flattice.set_flattice(radius_km=radius_km)
        # pplot_obj.flattice.PLOT_DIR = BASE_PLOT_DIR
        # pplot_obj.flattice.DNT = DNT
        # pplot_obj.flattice.satellites = satellites
        # pplot_obj.flattice.filter_method = filter_method
        # pplot_obj.flattice.cc_method = method
        # pplot_obj.flattice.cotfilt_value = chosen_cot # Testing adding extra parameter
        # plot_obj.flattice.isGAC = isGAC
        if filename is None:
            filename = "%s/h5files/f_lattice_%s_ccm_%s_%sfilter_dnt_%s_r%skm.h5" % (
                self.PLOT_DIR, self.satellites, self.cc_method, self.filter_method,
                self.DNT,  self.radius_km)
        self.h5_file_file = filename
        with h5py.File(filename, 'w') as f:
            for key in self.all_arrays.keys():
                print(key)
                data = self.all_arrays[key]
                if data is not None:
                    if hasattr(data, 'keys'):
                        # Some lattice data stored as dictionaries
                        for key_inner in data.keys():
                            f.create_dataset(key + str(key_inner), data=data[key_inner], compression=COMPRESS_LVL)
                    elif np.size(data) > 1:
                        f.create_dataset(key, data=data, compression=COMPRESS_LVL)
                    else:
                        f.create_dataset(key, data=data)
            # f.create_dataset('fibonacci_radius_km', self.radius_km)
            f.close()

    def set_flattice(self, radius_km=200):
        self.radius_km = radius_km
        self.lons, self.lats = get_fibonacci_spread_points_on_earth(radius_km=radius_km)
        fig = plt.figure(figsize=(36, 18))
        #ax = fig.add_subplot(111)
        #plt.plot(self.lons, self.lats, 'b*')
        # plt.show()
        self.sum_ctth_bias_low = np.zeros(self.lats.shape)
        self.sum_lapse_bias_low = np.zeros(self.lats.shape)
        self.sum_ctth_bias_high = np.zeros(self.lats.shape)
        self.sum_ctth_mae_low = np.zeros(self.lats.shape)
        self.sum_ctth_mae_high = np.zeros(self.lats.shape)
        self.sum_ctth_mae = np.zeros(self.lats.shape)
        self.sum_ctth_mae_diff = np.zeros(self.lats.shape)
        self.sum_lapse_bias_high = np.zeros(self.lats.shape)
        self.sum_ctth_bias_temperature_low = np.zeros(self.lats.shape)
        self.sum_ctth_bias_temperature_low_t11 = np.zeros(self.lats.shape)
        self.Min_lapse_rate = np.zeros(self.lats.shape)
        self.num_ctth_error_above_1km = np.zeros(self.lats.shape)
        self.num_new_false_clouds = np.zeros(self.lats.shape)
        self.num_new_detected_clouds = np.zeros(self.lats.shape)
        self.num_false_clouds = np.zeros(self.lats.shape)
        self.num_detected_clouds = np.zeros(self.lats.shape)
        self.num_undetected_clouds = np.zeros(self.lats.shape)
        self.num_detected_clouds_filtcot = {}
        self.num_undetected_clouds_filtcot = {}
        for cot in cots_mid:
            self.num_detected_clouds_filtcot[cot] = np.zeros(self.lats.shape)
            self.num_undetected_clouds_filtcot[cot] = np.zeros(self.lats.shape)
        self.num_detected_clear = np.zeros(self.lats.shape)
        self.num_detected_height_low = np.zeros(self.lats.shape)
        self.num_detected_height_high = np.zeros(self.lats.shape)
        self.num_detected_height = np.zeros(self.lats.shape)
        self.num_detected_height_both = np.zeros(self.lats.shape)
        self.Min_t11ts_offset = np.zeros(self.lats.shape)
        self.Max_t11t12_offset = np.zeros(self.lats.shape)
        self.Max_t37t12_offset = np.zeros(self.lats.shape)
        self.Max_t11t37_offset = np.zeros(self.lats.shape)
        self.sum_height_bias_type = {}
        self.num_detected_height_type = {}
        for cc_type in range(8):
            self.sum_height_bias_type[cc_type] = 1.0*np.zeros(self.lats.shape)
            self.num_detected_height_type[cc_type] = 1.0*np.zeros(self.lats.shape)

    def np_float_array(self):

        self.sum_ctth_mae_low = 1.0*np.array(self.sum_ctth_mae_low)
        self.sum_ctth_mae_high = 1.0*np.array(self.sum_ctth_mae_high)
        self.sum_ctth_mae = 1.0*np.array(self.sum_ctth_mae)
        self.sum_ctth_mae_diff = 1.0*np.array(self.sum_ctth_mae_diff)

        self.sum_ctth_bias_low = 1.0*np.array(self.sum_ctth_bias_low)
        self.sum_lapse_bias_low = 1.0*np.array(self.sum_lapse_bias_low)
        self.sum_ctth_bias_temperature_low = 1.0*np.array(self.sum_ctth_bias_temperature_low)
        self.sum_ctth_bias_temperature_low_t11 = 1.0*np.array(self.sum_ctth_bias_temperature_low_t11)
        self.Min_lapse_rate = 1.0*np.array(self.Min_lapse_rate)
        self.num_detected_clouds = 1.0*np.array(self.num_detected_clouds)
        self.num_undetected_clouds = 1.0*np.array(self.num_undetected_clouds)
        self.num_false_clouds = 1.0*np.array(self.num_false_clouds)
        self.num_detected_clear = 1.0*np.array(self.num_detected_clear)
        self.num_new_detected_clouds = 1.0*np.array(self.num_new_detected_clouds)
        self.num_new_false_clouds = 1.0*np.array(self.num_new_false_clouds)
        self.num_detected_height_low = 1.0*np.array(self.num_detected_height_low)

    def find_number_of_clouds_clear(self):
        self.np_float_array()
        self.num_clear = self.num_detected_clear + self.num_false_clouds
        self.num_clouds = self.num_detected_clouds + self.num_undetected_clouds
        self.N = self.num_clear + self.num_clouds
        self.Number_of = np.ma.masked_array(self.N, mask=self.N < 0)

    def _remap_a_score_on_an_area(self, plot_area_name='antarctica', 
                                  vmin=0.0, vmax=1.0,
                                  score='Kuipers', screen_out_valid=False):
        if 'robin' in plot_area_name:
            # Need crs from cartopy to get map boarder, the elipse. 
            crs = ccrs.Robinson()
            # {'a': 6378137.0, 'proj': 'robin', 'lon_0': 0}
            # crs.x_limits crs.y_limits 
            # (-17005833.33052523, 17005833.33052523)
            # (-8625155.12857459, 8625155.12857459)
            # Bounds needed for imshow or picture will be white
            crs.bounds = (crs.x_limits[0], crs.x_limits[1], 
                          crs.y_limits[0], crs.y_limits[1])
            # create pyresample area_def object
            area_def = AreaDefinition('robinson',
                                      'robinson',
                                      'robinson',
                                      projection=crs.proj4_params,
                                      width=1000, height=500, 
                                      area_extent=(crs.x_limits[0],  # not same as crs.bounds
                                                   crs.y_limits[0],
                                                   crs.x_limits[1],
                                                   crs.y_limits[1])
            )
        else:
            # For rectangular areas use the area_def from pyresample
            try:
                area_def = load_area(AREA_CONFIG_FILE_PLOTS_ON_AREA, plot_area_name)
            except AttributeError as e:
                print("Have you set AREA_CONFIG_FILE_PLOTS_ON_AREA?")
                print(e)
                raise(ValueError)
            crs = area_def.to_cartopy_crs()

        plot_label = score.replace('_', '-')
        # if "mae" in score:
        #     plot_label = ""

        data = getattr(self, score)
        data = data.copy()
        self.truncate_data(data, vmin, vmax)
        lons = self.lons
        lats = self.lats
        swath_def = geometry.SwathDefinition(lons=lons, lats=lats)
        result = resample_nearest(
            swath_def, data, area_def,
            radius_of_influence=self.radius_km*1000*2.5, fill_value=None)
        my_cmap, vmin, vmax = self._get_colormap(score, vmin, vmax, 
                                                 screen_out_valid=screen_out_valid)

        fig = plt.figure(figsize=(16, 9))
        ax = fig.subplots(nrows=1, subplot_kw={'projection': crs})  # Creates the area projection
        ax.coastlines()
        ax.set_global()
        ax.gridlines(edgecolor='black')
        if np.ma.is_masked(result):
            result.data[result.mask] = np.nan
            result.mask = result.data > vmax  # Empty mask, important
        else:
            result = np.ma.masked_array(result, mask=result > vmax)
        result.data[result.mask] = np.nan

        if 'robin' in plot_area_name:
            # imshow and meshgrid needs 2D data
            # hack to make data in corners white (data in corners caused by bug in PROJ)
            lon, lat = area_def.get_lonlats()
            yshape, xshape = result.data.shape
            result.data[:,0:xshape//2][lon[:,0:xshape//2]>0] = 2*vmax
            result.data[:,xshape//2:][lon[:,xshape//2:]<0] = 2*vmax
            
        # yshape, xshape = result.data.shape
        # xi = list(range(xshape))
        # yi = list(range(yshape))
        # lonsi, latsi = np.meshgrid(xi, yi)
        # print(lonsi.shape)
        # plt.pcolormesh(lon, lat, result.data, transform=crs,
        #               vmin=vmin, vmax=vmax, cmap=my_cmap)
        # plt.pcolormesh(lonsi, latsi, result.data, transform=crs,
        #               vmin=vmin, vmax=vmax, cmap=my_cmap)

        plt.imshow(result,  transform=crs, extent=crs.bounds,
                   vmin=vmin, vmax=vmax, label=plot_label, cmap=my_cmap)
        ax.set_global()
        plt.colorbar()
        plt.title(plot_label)
        filename = self.PLOT_DIR_SCORE + self.PLOT_FILENAME_START + plot_area_name + '.png'
        plt.savefig(filename, bbox_inches='tight')
        # import pyresample as pr
        # pr.plot.save_quicklook(self.PLOT_DIR_SCORE + self.PLOT_FILENAME_START +
        #                       plot_area_name + extra + 'quick.png',
        #                       area_def, result,
        #                       vmin=vmin, vmax=vmax, label=plot_label, cmap=my_cmap)


    def _get_colormap(self, score, vmin, vmax, screen_out_valid=False):
        """Find colormap for the score."""
        my_cmap = copy.copy(matplotlib.cm.BrBG)
        if "FAR" in score:
            # matplotlib.rcParams['image.cmap'] = "BrBG"
            my_cmap = copy.copy(matplotlib.cm.BrBG)
        elif "diff" in score:
            my_cmap = copy.copy(matplotlib.cm.BrBG)
        elif "mae" in score:
            my_cmap = copy.copy(matplotlib.cm.Reds)
        elif "ctth_pe" in score:
            my_cmap = copy.copy(matplotlib.cm.BrBG_r)
        else:
            my_cmap = copy.copy(matplotlib.cm.BrBG)
        if "cds" in score:
            from atrain_match.reshaped_files_scr.cds_colormap import get_cds_colormap
            my_cmap = get_cds_colormap()
        if score in "Bias" and screen_out_valid:
            # This screens out values between -5 and +5%
            vmax = 25
            vmin = -25
            my_cmap = copy.copy(matplotlib.cm.get_cmap("BrBG", lut=100))
            cmap_vals = my_cmap(np.arange(100))  # extractvalues as an array
            cmap_vals[39:61] = [0.9, 0.9, 0.9, 1]  # change the first value
            my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "newBrBG", cmap_vals)
        if score in "RMS" and screen_out_valid:
            # This screens out values beteen 0 and 20%. 41/100=20%
            vmax = 50
            vmin = 0
            my_cmap = copy.copy(matplotlib.cm.get_cmap("BrBG", lut=100))
            cmap_vals = my_cmap(np.arange(100))  # extract values as an array
            cmap_vals[0:41] = [0.9, 0.9, 0.9, 1]  # change the first value
            my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "newBrBG", cmap_vals)
        my_cmap.set_over(color='white', alpha=1.0)
        my_cmap.set_under(color='0.5', alpha=0.2)
        my_cmap.set_bad(color='0.5', alpha=0.2)
        return my_cmap, vmin, vmax

    def truncate_data(self, data, vmin, vmax):
        """Truncate data to [vmin, vmax]."""
        if np.ma.is_masked(data):
            the_mask = data.mask
            data[np.logical_and(data > vmax, ~the_mask)] = vmax
            data[np.logical_and(data < vmin, ~the_mask)] = vmin
            # do not wan't low scores hitrates set to nodata!
        else:
            data[data > vmax] = vmax
            data[data < vmin] = vmin  
 

    def remap_and_plot_score_on_several_areas(self, vmin=0.0, vmax=1.0,
                                              score='Kuipers', screen_out_valid=False):
        self.PLOT_DIR_SCORE = self.PLOT_DIR + "/%s/%s/" % (score, self.satellites)
        self.PLOT_FILENAME_START = "fig_%s_ccm_%s_%sfilter_dnt_%s_%s_r%skm_" % (
            self.satellites, self.cc_method, self.filter_method,
            self.DNT, score, self.radius_km)
 
        if not os.path.exists(self.PLOT_DIR_SCORE):
            os.makedirs(self.PLOT_DIR_SCORE)
        plt.close('all')    
        for plot_area_name in [
                # 'antarctica',
                # 'arctica',
                'robinson',
                'ease_nh',
                'ease_sh']:
            if plot_area_name.startswith('robin') and 'metop' in self.satellites:
                print("Not making robinson projection for metops.")
                continue
            self._remap_a_score_on_an_area(plot_area_name=plot_area_name,
                                           vmin=vmin, vmax=vmax, score=score, 
                                           screen_out_valid=screen_out_valid)

    def calculate_kuipers(self):
        """Calculate kuipers."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        num_clear = self.num_clear
        num_clouds = self.num_clouds
        num_detected_clouds = self.num_detected_clouds
        num_detected_clear = self.num_detected_clear
        # Typically we have num_clear/num_clouds = 30/70
        # In areas with only clouds or only clears the Kuipers will be ==0
        # Even if all clouds/clears are classified correctly!
        # Do something for these set to none or update
        Kuipers_devider = (num_clouds)*(num_clear)
        Kuipers_devider[Kuipers_devider == 0] = 1.0
        Kuipers = (num_detected_clouds*num_detected_clear -
                   self.num_false_clouds*self.num_undetected_clouds)/Kuipers_devider
        the_mask = np.logical_or(self.num_clear < 20, self.num_clouds < 20)
        the_mask = np.logical_or(the_mask, self.num_clouds < 0.01*self.num_clear)
        the_mask = np.logical_or(the_mask, self.num_clear < 0.01*self.num_clouds)
        use = np.logical_or(self.lats >= 70, ~the_mask)
        self.Kuipers_total_mean_polar = (
            (np.sum(num_detected_clouds[use])*np.sum(num_detected_clear[use]) -
             np.sum(self.num_false_clouds[use])*np.sum(self.num_undetected_clouds[use]))
            /((np.sum(num_clouds[use]))*(np.sum(num_clear[use]))))
        Kuipers = np.ma.masked_array(Kuipers, mask=the_mask)
        self.Kuipers = Kuipers

    def calculate_lapse_rate(self):
        """Calculate lapse_rate."""
        self.np_float_array()
        the_mask = self.Min_lapse_rate > - 0.001
        lapse_rate = np.ma.masked_array(self.Min_lapse_rate, mask=the_mask)
        self.lapse_rate = lapse_rate

    def calculate_t11ts_offset(self):
        """Calculate t11ts_offset."""
        self.np_float_array()
        the_mask = self.num_clear < 10
        t11ts_offset = np.ma.masked_array(self.Min_t11ts_offset, mask=the_mask)
        self.t11ts_offset = t11ts_offset

    def calculate_t11t12_offset(self):
        """Calculate t11t12_offset."""
        self.np_float_array()
        the_mask = self.num_clear < 10
        t11t12_offset = np.ma.masked_array(self.Max_t11t12_offset, mask=the_mask)
        self.t11t12_offset = t11t12_offset

    def calculate_t37t12_offset(self):
        """Calculate t37t12_offset."""
        self.np_float_array()
        the_mask = self.num_clear < 10
        t37t12_offset = np.ma.masked_array(self.Max_t37t12_offset, mask=the_mask)
        self.t37t12_offset = t37t12_offset

    def calculate_t11t37_offset(self):
        """Calculate t11t37_offset."""
        self.np_float_array()
        the_mask = self.num_clear < 10
        t11t37_offset = np.ma.masked_array(self.Max_t11t37_offset, mask=the_mask)
        self.t11t37_offset = t11t37_offset

    def calculate_temperature_bias(self):
        """Calculate temperature_bias."""
        self.np_float_array()
        the_mask = self.num_detected_height_low < 10
        ctth_bias_temperature_low = self.sum_ctth_bias_temperature_low*1.0/self.num_detected_height_low
        ctth_bias_temperature_low = np.ma.masked_array(ctth_bias_temperature_low, mask=the_mask)
        self.ctth_bias_temperature_low = ctth_bias_temperature_low

    def calculate_temperature_bias_t11(self):
        """Calculate temperature_bias_t11."""
        self.np_float_array()
        the_mask = self.num_detected_height_low < 10
        ctth_bias_temperature_low_t11 = self.sum_ctth_bias_temperature_low_t11*1.0/self.num_detected_height_low
        ctth_bias_temperature_low_t11 = np.ma.masked_array(ctth_bias_temperature_low_t11, mask=the_mask)
        self.ctth_bias_temperature_low_t11 = ctth_bias_temperature_low_t11

    def calculate_height_bias(self):
        """Calculate height_bias."""
        self.np_float_array()
        the_mask = self.num_detected_height_low < 10
        ctth_bias_low = self.sum_ctth_bias_low*1.0/self.num_detected_height_low
        ctth_bias_low = np.ma.masked_array(ctth_bias_low, mask=the_mask)
        self.ctth_bias_low = ctth_bias_low
        the_mask = self.num_detected_height_high < 10
        ctth_bias_high = self.sum_ctth_bias_high*1.0/self.num_detected_height_high
        ctth_bias_high = np.ma.masked_array(ctth_bias_high, mask=the_mask)
        self.ctth_bias_high = ctth_bias_high

    def calculate_height_mae(self):
        """Calculate height_mae."""
        self.np_float_array()
        the_mask = self.num_detected_height_low < 10
        ctth_mae_low = self.sum_ctth_mae_low*1.0/self.num_detected_height_low
        ctth_mae_low = np.ma.masked_array(ctth_mae_low, mask=the_mask)
        self.ctth_mae_low = ctth_mae_low
        the_mask = self.num_detected_height_high < 10
        ctth_mae_high = self.sum_ctth_mae_high*1.0/self.num_detected_height_high
        ctth_mae_high = np.ma.masked_array(ctth_mae_high, mask=the_mask)
        self.ctth_mae_high = ctth_mae_high
        the_mask = self.num_detected_height < 10
        ctth_mae = self.sum_ctth_mae*1.0/self.num_detected_height
        ctth_mae = np.ma.masked_array(ctth_mae, mask=the_mask)
        self.ctth_mae = ctth_mae
        the_mask = self.num_detected_height_both < 10
        ctth_mae_diff = self.sum_ctth_mae_diff*1.0/self.num_detected_height_both
        ctth_mae_diff = np.ma.masked_array(ctth_mae_diff, mask=the_mask)
        self.ctth_mae_diff = ctth_mae_diff
        the_mask = self.num_detected_height_both < 0
        self.num_detected_height_both = np.ma.masked_array(self.num_detected_height_both, mask=the_mask)

    def calculate_ctth_pe1(self):
        """Calculate ctth_pe1."""
        self.np_float_array()
        the_mask = self.num_detected_height < 10
        ctth_pe1 = self.num_ctth_error_above_1km*100.0/self.num_detected_height
        ctth_pe1 = np.ma.masked_array(ctth_pe1, mask=the_mask)
        self.ctth_pe1 = ctth_pe1

    def calculate_height_bias_lapse(self):
        """Calculate height_bias_lapse."""
        self.np_float_array()
        the_mask = self.num_detected_height_low < 10
        lapse_bias_low = self.sum_lapse_bias_low*1.0/self.num_detected_height_low
        lapse_bias_low = np.ma.masked_array(lapse_bias_low, mask=the_mask)
        self.lapse_bias_low = lapse_bias_low
        the_mask = self.num_detected_height_high < 10
        lapse_bias_high = self.sum_lapse_bias_high*1.0/self.num_detected_height_high
        lapse_bias_high = np.ma.masked_array(lapse_bias_high, mask=the_mask)
        self.lapse_bias_high = lapse_bias_high

    def calculate_height_bias_type(self):
        """Calculate height_bias_type."""
        self.np_float_array()
        for cc_type in range(8):
            the_mask = self.num_detected_height_type[cc_type] < 10
            ctth_bias_type_i = self.sum_height_bias_type[cc_type]*1.0/self.num_detected_height_type[cc_type]
            ctth_bias_type_i = np.ma.masked_array(ctth_bias_type_i, mask=the_mask)
            setattr(self, "ctth_bias_type_{:d}".format(cc_type), ctth_bias_type_i)

    def calculate_cds(self):
        """Calculate cds."""
        self.np_float_array()
        the_mask = self.num_detected_clouds + self.num_undetected_clouds < 100
        cds = 0 * self.num_detected_clouds.copy()
        cots_largest_first = cots_mid[::-1]
        for cot in cots_largest_first:
            update = np.logical_and(
                self.num_detected_clouds_filtcot[cot] > 0,
                self.num_detected_clouds_filtcot[cot] >= self.num_undetected_clouds_filtcot[cot])
            cds[update] = cot
        cds = np.ma.masked_array(cds, mask=the_mask)
        setattr(self, "cds", cds)

    def calculate_hitrate(self):
        """Calculate hitrate."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        Hitrate = (self.num_detected_clouds + self.num_detected_clear)*1.0/(
            self.num_clear + self.num_clouds)
        the_mask = self.N < 20
        Hitrate = np.ma.masked_array(Hitrate, mask=the_mask)
        use = np.logical_or(self.lats >= 70, ~the_mask)
        self.Hitrate_total_mean_polar = (
            np.sum(self.num_detected_clouds[use]) +
            np.sum(self.num_detected_clear[use]))*1.0/(
                np.sum(self.N[use]))
        self.Hitrate = Hitrate

    def calculate_increased_hitrate(self):
        """Calculate increased_hitrate."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        Hitrate = (
            self.num_detected_clouds + self.num_detected_clear)*1.0/(
                self.num_clear + self.num_clouds)
        new_Hitrate = (
            self.num_detected_clouds + self.num_new_detected_clouds -
            self.num_new_false_clouds + self.num_detected_clear)*1.0/(
                self.num_clear + self.num_clouds)
        the_mask = self.N < 20
        increased_Hitrate = np.ma.masked_array(new_Hitrate - Hitrate, mask=the_mask)
        self.increased_Hitrate = increased_Hitrate

    def calculate_threat_score(self):
        """Calculate threat_score."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        ThreatScore = (
            self.num_detected_clouds)*1.0/(self.num_clouds + self.num_false_clouds)
        the_mask = self.num_clouds < 20
        ThreatScore = np.ma.masked_array(ThreatScore, mask=the_mask)
        self.Threat_Score = ThreatScore

    def calculate_threat_score_clear(self):
        """Calculate threat_score_clear."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        ThreatScoreClear = (
            self.num_detected_clear)*1.0/(self.num_clear +
                                          self.num_undetected_clouds)
        the_mask = self.num_clear < 20
        ThreatScoreClear = np.ma.masked_array(ThreatScoreClear,
                                              mask=the_mask)
        self.Threat_Score_Clear = ThreatScoreClear

    def calculate_pod_clear(self):
        """Calculate pod_clear."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        PODclear = (
            self.num_detected_clear)*1.0/(self.num_clear)
        the_mask = self.num_clear < 20
        PODclear = np.ma.masked_array(PODclear, mask=the_mask)
        self.PODclear = PODclear
        use = np.logical_or(self.lats >= 70, ~the_mask)
        self.PODclear_total_mean_polar = (
            np.sum(self.num_detected_clear[use]))*1.0/(
                np.sum(self.num_clear[use]))

    def calculate_pod_cloudy(self):
        """Calculate pod_cloudy."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        PODcloudy = (
            self.num_detected_clouds)*1.0/(self.num_clouds)
        the_mask = self.num_clouds < 20
        PODcloudy = np.ma.masked_array(PODcloudy, mask=the_mask)
        self.PODcloudy = PODcloudy
        use = np.logical_or(self.lats >= 70, ~the_mask)
        self.PODcloudy_total_mean_polar = (
            np.sum(self.num_detected_clouds[use]))*1.0/(
                np.sum(self.num_clouds[use]))

    def calculate_far_clear(self):
        """Calculate far_clear."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        FARclear = (
            self.num_undetected_clouds)*1.0/(self.num_detected_clear +
                                             self.num_undetected_clouds)
        the_mask = self.num_clear < 20
        FARclear = np.ma.masked_array(FARclear, mask=the_mask)
        self.FARclear = FARclear

    def calculate_far_cloudy(self):
        """Calculate far_cloudy."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        FARcloudy = (
            self.num_false_clouds)*1.0/(self.num_detected_clouds +
                                        self.num_false_clouds)
        the_mask = self.num_clouds < 20
        FARcloudy = np.ma.masked_array(FARcloudy, mask=the_mask)
        self.FARcloudy = FARcloudy
        use = np.logical_or(self.lats >= 70, ~the_mask)
        self.FARcloudy_total_mean_polar = (
            np.sum(self.num_false_clouds[use]))*1.0/(
                np.sum(self.num_detected_clouds[use]) +
                np.sum(self.num_false_clouds[use]))

    def calculate_calipso_cfc(self):
        """Calculate calipso_cfc."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        calipso_cfc = 100*(
            self.num_detected_clouds + self.num_undetected_clouds)*1.0/(self.N)
        the_mask = self.N < 20
        calipso_cfc = np.ma.masked_array(calipso_cfc, mask=the_mask)
        self.calipso_cfc = calipso_cfc

    def calculate_pps_cfc(self):
        """Calculate pps_cfc."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        pps_cfc = 100*(
            self.num_detected_clouds + self.num_false_clouds)*1.0/(self.N)
        the_mask = self.N < 20
        pps_cfc = np.ma.masked_array(pps_cfc, mask=the_mask)
        self.pps_cfc = pps_cfc

    def calculate_bias(self):
        """Calculate bias."""
        self.np_float_array()
        self.find_number_of_clouds_clear()
        self.calculate_calipso_cfc()
        self.calculate_pps_cfc()
        Bias = self.pps_cfc - self.calipso_cfc
        the_mask = self.N < 20
        Bias = np.ma.masked_array(Bias, mask=the_mask)
        use = np.logical_or(self.lats >= 70, ~the_mask)
        self.Bias_total_mean_polar = 100*(
            - np.sum(self.num_undetected_clouds[use]) +
            np.sum(self.num_false_clouds[use]))*1.0/(np.sum(self.N[use]))
        self.Bias = Bias

    def calculate_rms(self):
        """Calculate rms."""
        self.np_float_array()
        self.calculate_calipso_cfc()
        self.find_number_of_clouds_clear()
        self.calculate_bias()
        RMS = np.sqrt((self.num_false_clouds*(100.0 - 0.0 - self.Bias)**2 +
                       self.num_undetected_clouds*(0.0 - 100.0 - self.Bias)**2 +
                       self.num_detected_clear*self.Bias**2 +
                       self.num_detected_clouds*self.Bias**2)/(
            self.N))
        the_mask = self.N < 20
        RMS = np.ma.masked_array(RMS, mask=the_mask)
        self.RMS = RMS


class PerformancePlottingObject:
    def __init__(self):
        self.flattice = StatsOnFibonacciLattice()

    def add_detection_stats_on_fib_lattice(self, my_obj):
        # Start with the area and get lat and lon to calculate the stats:
        if len(my_obj.longitude) == 0:
            print("Skipping file, no matches !")
            return
        lats = self.flattice.lats[:]
        max_distance = self.flattice.radius_km*1000*2.5
        area_def = SwathDefinition(*(self.flattice.lons,
                                     self.flattice.lats))
        target_def = SwathDefinition(*(my_obj.longitude,
                                       my_obj.latitude))
        valid_in, valid_out, indices, distances = get_neighbour_info(
            area_def, target_def, radius_of_influence=max_distance,
            epsilon=100, neighbours=1)
        cols = get_sample_from_neighbour_info('nn', target_def.shape,
                                              np.array(range(0, len(lats))),
                                              valid_in, valid_out,
                                              indices)
        cols = cols[valid_out]
        detected_clouds = my_obj.detected_clouds[valid_out]
        detected_clear = my_obj.detected_clear[valid_out]
        detected_height_low = my_obj.detected_height_low[valid_out]
        detected_height_high = my_obj.detected_height_high[valid_out]
        detected_height = my_obj.detected_height[valid_out]
        detected_height_both = my_obj.detected_height_both[valid_out]
        false_clouds = my_obj.false_clouds[valid_out]
        undetected_clouds = my_obj.undetected_clouds[valid_out]
        new_detected_clouds = my_obj.new_detected_clouds[valid_out]
        undetected_clouds_filtcot = {}
        detected_clouds_filtcot = {}
        for cot in cots_mid:
            undetected_clouds_filtcot[cot] = my_obj.undetected_clouds_filtcot[cot][valid_out]
            detected_clouds_filtcot[cot] = my_obj.detected_clouds_filtcot[cot][valid_out]

        new_false_clouds = my_obj.new_false_clouds[valid_out]
        lapse_rate = my_obj.lapse_rate[valid_out]
        t11ts_offset = my_obj.t11ts_offset[valid_out]
        t11t12_offset = my_obj.t11t12_offset[valid_out]
        t37t12_offset = my_obj.t37t12_offset[valid_out]
        t11t37_offset = my_obj.t11t37_offset[valid_out]
        height_bias_low = my_obj.height_bias_low[valid_out]
        height_bias = my_obj.height_bias[valid_out]
        height_mae_diff = my_obj.height_mae_diff[valid_out]
        temperature_bias_low = my_obj.temperature_bias_low[valid_out]
        temperature_bias_low_t11 = my_obj.temperature_bias_low_t11[valid_out]
        lapse_bias_low = my_obj.lapse_bias_low[valid_out]
        height_bias_high = my_obj.height_bias_high[valid_out]
        lapse_bias_high = my_obj.lapse_bias_high[valid_out]
        is_clear = np.logical_or(detected_clear, false_clouds)
        # lets make things faster, I'm tired of waiting!
        cols[distances > max_distance] = -9  # don't use pixles matched too far away!
        import time
        tic = time.time()
        arr, counts = np.unique(cols, return_index=False, return_counts=True)
        for d in arr[arr > 0]:
            use = cols == d
            ind = np.where(use)[0]
            # if ind.any():
            self.flattice.num_false_clouds[d] += np.sum(false_clouds[ind])
            self.flattice.num_detected_clouds[d] += np.sum(detected_clouds[ind])
            self.flattice.num_detected_clear[d] += np.sum(detected_clear[ind])
            self.flattice.num_undetected_clouds[d] += np.sum(undetected_clouds[ind])
            self.flattice.num_new_false_clouds[d] += np.sum(new_false_clouds[ind])
            self.flattice.num_new_detected_clouds[d] += np.sum(new_detected_clouds[ind])
            self.flattice.num_detected_height_low[d] += np.sum(detected_height_low[ind])
            self.flattice.num_detected_height_high[d] += np.sum(detected_height_high[ind])
            self.flattice.num_detected_height[d] += np.sum(detected_height[ind])
            self.flattice.num_detected_height_both[d] += np.sum(detected_height_both[ind])
            self.flattice.sum_ctth_bias_low[d] += np.sum(height_bias_low[ind])
            self.flattice.sum_ctth_mae_low[d] += np.sum(np.abs(height_bias_low[ind]))
            self.flattice.sum_ctth_mae[d] += np.sum(np.abs(height_bias[ind]))
            self.flattice.num_ctth_error_above_1km[d] += np.sum(np.abs(height_bias[ind]) > 1000)
            self.flattice.sum_ctth_mae_diff[d] += np.sum(height_mae_diff[ind])
            self.flattice.sum_lapse_bias_low[d] += np.sum(lapse_bias_low[ind])
            self.flattice.sum_ctth_bias_high[d] += np.sum(height_bias_high[ind])
            self.flattice.sum_ctth_mae_high[d] += np.sum(np.abs(height_bias_high[ind]))
            self.flattice.sum_lapse_bias_high[d] += np.sum(lapse_bias_high[ind])
            self.flattice.sum_ctth_bias_temperature_low[d] += np.sum(temperature_bias_low[ind])
            self.flattice.sum_ctth_bias_temperature_low_t11[d] += np.sum(temperature_bias_low_t11[ind])
            self.flattice.Min_lapse_rate[d] = np.min([self.flattice.Min_lapse_rate[d],
                                                      np.min(lapse_rate[ind])])
            if np.sum(is_clear[ind]) > 0:
                self.flattice.Min_t11ts_offset[d] = np.min([self.flattice.Min_t11ts_offset[d],
                                                            np.percentile(t11ts_offset[ind][is_clear[ind]], 5)])
                self.flattice.Max_t11t12_offset[d] = np.max([self.flattice.Max_t11t12_offset[d],
                                                             np.percentile(t11t12_offset[ind][is_clear[ind]], 95)])
                self.flattice.Max_t37t12_offset[d] = np.max([self.flattice.Max_t37t12_offset[d],
                                                             np.percentile(t37t12_offset[ind][is_clear[ind]], 95)])
                self.flattice.Max_t11t37_offset[d] = np.max([self.flattice.Max_t11t37_offset[d],
                                                             np.percentile(t11t37_offset[ind][is_clear[ind]], 95)])
            cc_type = 0
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 1
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 2
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 3
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 4
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 5
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 6
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 7
            self.flattice.sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.num_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])

            for cot in cots_mid:
                self.flattice.num_detected_clouds_filtcot[cot][d] += np.sum(detected_clouds_filtcot[cot][ind])
                self.flattice.num_undetected_clouds_filtcot[cot][d] += np.sum(undetected_clouds_filtcot[cot][ind])

        print("mapping took %1.4f seconds" % (time.time()-tic))


def get_fibonacci_spread_points_on_earth(radius_km):
    # Earth area = 510072000km2
    # 4000 point with radius~200km
    # 1000 point with radium~100km
    # 25000 radius 80km
    # 64000 radius 5km
    EARTH_AREA = 510072000
    POINT_AREA = radius_km * radius_km * 3.14
    n = int(EARTH_AREA / POINT_AREA)
    # http://arxiv.org/pdf/0912.4540.pdf
    # Alvaro Gonzalez: Measurement of areas on sphere usig Fibonacci and latitude-longitude grid.
    # import math
    lin_space = np.array(range(- n//2, n//2))
    pi = 3.14
    theta = (1 + np.sqrt(5))*0.5
    longitude = (lin_space % theta)*360/theta
    temp = (2.0 * lin_space) / (n)
    temp[temp > 1.0] = 0.999
    temp[temp < - 1.0] = -0.999
    latitude = np.arcsin(temp)*180/pi
    longitude[longitude > 180] = longitude[longitude > 180] - 360
    longitude[longitude < - 180] = longitude[longitude < - 180] + 360
    # latitude[latitude>90]=180 - latitude[latitude>90]
    # latitude[latitude<-90]=-180 -latitude[latitude<-90]
    longitude = longitude[latitude < 90]
    latitude = latitude[latitude < 90]
    longitude = longitude[latitude > - 90]
    latitude = latitude[latitude > - 90]

    if np.isnan(np.max(latitude)):
        raise ValueError
    return longitude, latitude
