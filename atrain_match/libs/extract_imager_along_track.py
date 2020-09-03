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
import logging
import os
logger = logging.getLogger(__name__)
CHANNEL_MICRON_DESCRIPTIONS = {'11': ["avhrr channel 4 - 11um",
                                      "Avhrr channel channel4.",
                                      "AVHRR ch4",
                                      "AVHRR ch 4",
                                      "AVHRR ch_4",
                                      "SEVIRI IR_108",
                                      "channel4",
                                      "AVHRR 4",
                                      "AVHRR-2 4",
                                      "AVHRR-3 4",
                                      "MODIS 31",
                                      "VIIRS M15",
                                      "Avhrr channel channel4."],
                               '12': ["avhrr channel 5 - 12um",
                                      "Avhrr channel channel5.",
                                      "AVHRR ch5",
                                      "AVHRR ch 5",
                                      "AVHRR ch_5",
                                      "SEVIRI IR_120",
                                      "channel5",
                                      "AVHRR 5",
                                      "AVHRR-2 5",
                                      "AVHRR-3 5",
                                      "MODIS 32",
                                      "VIIRS M16",
                                      "Avhrr channel channel5."],
                               '06': ["VIIRS M05",
                                      "AVHRR ch 1",
                                      "SEVIRI VIS006",
                                      "AVHRR ch1",
                                      "AVHRR ch_1",
                                      "AVHRR 1",
                                      "AVHRR-2 1",
                                      "AVHRR-3 1",
                                      "MODIS 1"],
                               '09': ["VIIRS M07",
                                      "AVHRR ch 2",
                                      "AVHRR ch2",
                                      "AVHRR ch_2",
                                      "AVHRR 2",
                                      "AVHRR-2 2",
                                      "AVHRR-3 2",
                                      "SEVIRI VIS008",
                                      "MODIS 2"],
                               '16': ["VIIRS M10",
                                      "AVHRR ch 3a",
                                      "3a",
                                      "SEVIRI IR_016",
                                      "AVHRR ch3a",
                                      "AVHRR ch_3a",
                                      "AVHRR 3A",
                                      "AVHRR-2 3A",
                                      "AVHRR-3 3A",
                                      "MODIS 6"],
                               '37': ["VIIRS M12",
                                      "AVHRR ch 3b",
                                      "SEVIRI IR_039",
                                      "AVHRR ch_3b",
                                      "AVHRR ch3b",
                                      "AVHRR-2 3B",
                                      "AVHRR-3 3B",
                                      "AVHRR 3B",
                                      "MODIS 20",
                                      "3b",
                                      "AVHRR 3B"],
                               '22': ["VIIRS M11"],
                               '13': ["VIIRS M09",
                                      "MODIS 26"],
                               '86': ["VIIRS M14",
                                      "SEVIRI IR_087",
                                      "MODIS 29"],
                               'seviri_bt73': ['SEVIRI WV_073'],
                               'seviri_bt134': ['SEVIRI IR_134'],
                               'seviri_bt97': ['SEVIRI IR_097'],
                               'seviri_bt62': ['SEVIRI WV_062'],

                               'modis_3': ['MODIS 3'],
                               'modis_4': ['MODIS 4'],
                               'modis_5': ['MODIS 5'],
                               'modis_7': ['MODIS 7'],
                               'modis_8': ['MODIS 8'],
                               'modis_9': ['MODIS 9'],
                               'modis_10': ['MODIS 10'],
                               'modis_11': ['MODIS 11'],
                               'modis_12': ['MODIS 12'],
                               'modis_13lo': ['MODIS 13lo'],
                               'modis_13hi': ['MODIS 13hi'],
                               'modis_14lo': ['MODIS 14lo'],
                               'modis_14hi': ['MODIS 14hi'],
                               'modis_15': ['MODIS 15'],
                               'modis_16': ['MODIS 16'],
                               'modis_17': ['MODIS 17'],
                               'modis_18': ['MODIS 18'],
                               'modis_19': ['MODIS 19'],
                               'modis_21': ['MODIS 21'],
                               'modis_22': ['MODIS 22'],
                               'modis_23': ['MODIS 23'],
                               'modis_24': ['MODIS 24'],
                               'modis_25': ['MODIS 25'],
                               'modis_27': ['MODIS 27'],
                               'modis_28': ['MODIS 28'],
                               'modis_30': ['MODIS 30'],
                               'modis_33': ['MODIS 33'],
                               'modis_34': ['MODIS 34'],
                               'modis_35': ['MODIS 35'],
                               'modis_36': ['MODIS 36']
                               }
CHANNEL_MICRON_IMAGER_PPS = {'11': 3,
                             '12': 4,
                             '06': 0,
                             '09': 1,
                             '37': 2,
                             '86': -1,
                             '16': 5,
                             '22': -1,
                             '13': -1}
CURRENTLY_UNUSED_SEVIRI_CHANNELS = ['seviri_bt73',
                                    'seviri_bt134',
                                    'seviri_bt97',
                                    'seviri_bt62']
CURRENTLY_UNUSED_MODIS_CHANNELS = ['modis_3',
                                   'modis_4',
                                   'modis_5',
                                   'modis_7',
                                   'modis_8',
                                   'modis_9',
                                   'modis_10',
                                   'modis_11',
                                   'modis_12',
                                   'modis_13lo',
                                   'modis_13hi',
                                   'modis_14lo',
                                   'modis_14hi',
                                   'modis_15',
                                   'modis_16',
                                   'modis_17',
                                   'modis_18',
                                   'modis_19',
                                   'modis_21',
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
                                   'modis_36']


def get_data_from_array(array, matched):
    if array is None:
        return None
    return np.array([array[matched['row'][idx], matched['col'][idx]]
                     for idx in range(matched['row'].shape[0])])


def get_data_from_array_nneigh(array, matched):
    if array is None:
        return None
    out = np.zeros(matched['row'].shape)
    for i in range(matched['row'].shape[1]):
        nodata = matched['row'][:, i] < 0
        matched_i = {'row': matched['row'][:, i].copy(),
                     'col': matched['col'][:, i].copy()}
        matched_i['row'][nodata] = 0
        matched_i['col'][nodata] = 0
        out[:, i] = get_data_from_array(array, matched_i)
        out[nodata, i] = -9
    return out


def get_mean_data_from_array_nneigh(array, matched):
    if array is None:
        return None
    out_all = get_data_from_array_nneigh(array, matched)
    n_ok = np.sum(out_all >= 0, axis=-1)  # before
    out_all[out_all < 0] = 0
    sum_out_all = np.sum(out_all, axis=-1)
    return sum_out_all * 1.0 / n_ok


def get_channel_data_from_object(imager_obj, chn_des, matched, nodata=-9):
    """Get the IMAGER/VIIRS channel data on the track

    matched: dict of matched indices (row, col)

    """
    channels = imager_obj.channel
    numOfChannels = len(channels)
    # for ich in range(numOfChannels):
    #    if channels[ich].des in CHANNEL_MICRON_DESCRIPTIONS[chn_des]:
    #        chnum = ich
    chnum = [ich for ich in range(numOfChannels)
             if channels[ich].des in CHANNEL_MICRON_DESCRIPTIONS[chn_des]]

    if len(chnum) == 0:
        # chnum = CHANNEL_MICRON_IMAGER_PPS[chn_des]
        logger.debug("Did not find pps channel number for channel "
                     "{:s}".format(chn_des))
        return None, ""
    else:
        chnum = chnum[0]
    if matched is None:
        return channels[chnum].data, ""

    chdata_on_track = get_data_from_array(channels[chnum].data, matched)
    # np.array([channels[chnum].data[matched['row'][idx], matched['col'][idx]]
    #                   for idx in range(matched['row'].shape[0])])
    extra_info = ""
    if channels[chnum].SZA_corr_done:
        extra_info = "_sza_correction_done"
    return np.array(chdata_on_track), extra_info


def _interpolate_height_and_temperature_from_pressure(imager_obj,
                                                      level, list_of_levels=None):
    """ Function to find height att pressure level (level)
    from segment_nwp, pressure and height vectors.
    High means high in pressure. The level closest to ground i hi, and lo is at lower
    pressure further up in atmosphere.
    """
    if hasattr(imager_obj, "nwp_height") and imager_obj.nwp_height is not None:
        values_h = imager_obj.nwp_height
        pressure_v = imager_obj.nwp_pressure
        surface_h = imager_obj.nwp_surface_h
        psur = imager_obj.nwp_psur
    elif hasattr(imager_obj, "segment_nwp_geoheight") and imager_obj.segment_nwp_geoheight is not None:
        values_h = imager_obj.segment_nwp_geoheight
        pressure_v = imager_obj.segment_nwp_pressure
        surface_h = imager_obj.segment_nwp_surfaceGeoHeight
        psur = imager_obj.segment_nwp_surfacePressure
    else:
        return None
    # import pdb
    # pdb.set_trace()
    nlev = pressure_v.shape[1]
    npix = pressure_v.shape[0]
    k = np.arange(npix)
    if list_of_levels is None:
        higher_index = np.array([nlev - 1 - np.searchsorted(pressure_v[ind, :], level, side='right',
                                                            sorter=range(nlev - 1, -1, -1))
                                 for ind in range(npix)])
    else:
        higher_index = np.array([nlev - 1 - np.searchsorted(pressure_v[ind, :], list_of_levels[ind], side='right',
                                                            sorter=range(nlev - 1, -1, -1))
                                 for ind in range(npix)])
        level = list_of_levels
    higher_index[higher_index >= (nlev - 1)] = nlev - 2
    lower_index = higher_index + 1
    # update "lo" where level is between surface and first level in array
    below_level_1 = level > pressure_v[:, 0]
    lower_index[below_level_1] = 0
    # get pressure and height for layer below and above level
    hi = pressure_v[k, higher_index]
    lo = pressure_v[k, lower_index]
    height_hi_ = values_h[k, higher_index]*1.0
    height_lo_ = values_h[k, lower_index]*1.0
    # update "hi" where level is between surface and first level in array
    hi[below_level_1] = psur[below_level_1]
    height_hi_[below_level_1] = surface_h[below_level_1]
    # log pressures
    hi = np.log(hi)
    lo = np.log(lo)
    level = np.log(level)
    # interpolate
    out_h = height_hi_ - (hi - level) * (height_hi_ - height_lo_) / (hi - lo)
    return out_h


def insert_nwp_segments_data(nwp_segments, row_matched, col_matched, obt):
    npix = row_matched.shape[0]
    """
        # obt.imager.segment_nwgeoheight
        obt.imager.segment_nwp_moist
        obt.imager.segment_nwp_pressure
        obt.imager.segment_nwp_temp
        obt.imager.segment_surfaceLandTemp
        obt.imager.segment_surfaceSeaTemp
        obt.imager.segment_surfaceGeoHeight
        obt.imager.segment_surfaceMoist
        obt.imager.segment_surfacePressure
        obt.imager.segment_fractionOfLand
        obt.imager.segment_meanElevation
        obt.imager.segment_ptro
        obt.imager.segment_ttro
        # obt.imager.segment_t850
        obt.imager.segment_tb11clfree_sea
        obt.imager.segment_tb12clfree_sea
        obt.imager.segment_tb11clfree_land
        obt.imager.segment_tb12clfree_land
        obt.imager.segment_tb11cloudy_surface
        obt.imager.segment_tb12cloudy_surface
        """
    def get_segment_row_col_idx(nwp_segments, row_matched, col_matched):
        segment_colidx = nwp_segments['colidx']
        segment_rowidx = nwp_segments['rowidx']
        seg_row = np.zeros(np.size(row_matched)) - 9
        seg_col = np.zeros(np.size(col_matched)) - 9
        for s_col in range(nwp_segments['norows']):
            for s_row in range(nwp_segments['nocols']):
                within_segment = np.logical_and(
                    np.logical_and(
                        row_matched >= (segment_rowidx[s_row, s_col]
                                        - nwp_segments['segSizeX']/2),
                        row_matched < (segment_rowidx[s_row, s_col]
                                       + nwp_segments['segSizeX']/2)),
                    np.logical_and(
                        col_matched >= (segment_colidx[s_row, s_col]
                                        - nwp_segments['segSizeY']/2),
                        col_matched < (segment_colidx[s_row, s_col]
                                       + nwp_segments['segSizeY']/2)))
                seg_row[within_segment] = s_row
                seg_col[within_segment] = s_col
        return seg_row.astype(np.int16), seg_col.astype(np.int16)
    seg_row, seg_col = get_segment_row_col_idx(nwp_segments, row_matched, col_matched)
    for data_set in ['surfaceLandTemp',
                     'surfaceSeaTemp',
                     'surfaceGeoHeight',
                     'surfaceMoist',
                     'surfacePressure',
                     'fractionOfLand',
                     'meanElevation',
                     'ptro',
                     'ttro',
                     't850',
                     'tb11clfree_sea',
                     'tb12clfree_sea',
                     'tb11clfree_land',
                     'tb12clfree_land',
                     'tb11lowcloud_sea',
                     'tb12lowcloud_sea',
                     'tb11lowcloud_land',
                     'tb12lowcloud_land']:
        if data_set in nwp_segments.keys():
            # 'tb11cloudy_surface',
            # 'tb12cloudy_surface ',
            setattr(obt.imager, 'segment_nwp_' + data_set,
                    np.array([nwp_segments[data_set][seg_row[idx], seg_col[idx]]
                              for idx in range(npix)]))
        elif 'clfree' in data_set or 'lowcloud' in data_set:
            # these are nor always present
            pass

    for data_set in ['moist', 'pressure', 'geoheight', 'temp']:
        setattr(obt.imager, 'segment_nwp_' + data_set,
                np.array([nwp_segments[data_set][seg_row[idx], seg_col[idx]]
                          for idx in range(npix)]))
    # Remove nodata and not used upper part of atmosphere
    N = obt.imager.segment_nwp_pressure.shape[1]
    pressure_n_to_keep = np.sum(np.max(obt.imager.segment_nwp_pressure, axis=0) > 50)
    logger.debug("Not saving upper %d levels of 3-D nwp from segment file" % (N-pressure_n_to_keep))
    logger.debug("Keeping %d lower levels of 3-D nwp from segment file" % (pressure_n_to_keep))
    for data_set in ['segment_nwp_moist', 'segment_nwp_pressure',
                     'segment_nwp_geoheight', 'segment_nwp_temp']:
        data = getattr(obt.imager, data_set)
        setattr(obt.imager, data_set, data[:, 0:pressure_n_to_keep])
    return obt


def insert_nwp_h440_h680_data(obt):
    data = _interpolate_height_and_temperature_from_pressure(obt.imager, 440)
    setattr(obt.imager, 'segment_nwp_h440', data)
    data = _interpolate_height_and_temperature_from_pressure(obt.imager, 680)
    setattr(obt.imager, 'segment_nwp_h680', data)
    return obt


# ---------------------------------------------------------------------------

from atrain_match.cloudproducts.read_oca import OCA_READ_EXTRA
def imager_track_from_matched(obt, SETTINGS, cloudproducts,
                              extract_radiances=True,
                              extract_cma=True,
                              extract_ctth=True,
                              extract_ctype=True,
                              extract_cpp=True,
                              extract_aux_segments=True,
                              extract_aux=True,
                              aux_params=None,
                              extract_some_data_for_x_neighbours=False,
                              find_mean_data_for_x_neighbours=False):
    aux_params_all = ["surftemp",
                      "t500", "t700", "t850", "t950", "ttro",
                      "ciwv",
                      "t900", "t1000", "t800", "t250", "t2m",
                      "ptro", "psur",
                      "h2m", "u10m", "v10m", "t2m",
                      "snowa", "snowd", "seaice",
                      "landuse", "fractionofland", "elevation",
                      "r37_sza_correction_done",

    ] + OCA_READ_EXTRA # And NN-extra!
    
    aux_obj = cloudproducts.aux
    if aux_params is None:
        #aux_params = aux_params_all
        aux_params = dir(aux_obj)
        # For amsr-E matching (many neighbors) use only the needed nwp data
        aux_params = [param for param in aux_params if '__' not in param]
        aux_params = [param for param in aux_params if 'CTTH' not in param]
    
    imager_obj = cloudproducts.imager_channeldata
    angle_obj = cloudproducts.imager_angles
    ctth = cloudproducts.ctth
    cma = cloudproducts.cma
    ctype = cloudproducts.ctype
    cpp = cloudproducts.cpp
    nwp_segments = cloudproducts.nwp_segments

    truth = getattr(obt, obt.truth_sat)
    row_matched = truth.imager_linnum
    col_matched = truth.imager_pixnum
    row_col = {'row': row_matched, 'col': col_matched}
    if extract_some_data_for_x_neighbours or find_mean_data_for_x_neighbours:
        row_matched_nneigh = truth.imager_linnum_nneigh
        col_matched_nneigh = truth.imager_pixnum_nneigh
        row_col_nneigh = {'row': row_matched_nneigh, 'col': col_matched_nneigh}

    obt.imager.latitude = get_data_from_array(cloudproducts.latitude, row_col)
    obt.imager.longitude = get_data_from_array(cloudproducts.longitude, row_col)
    if extract_ctype and ctype is not None:
        obt.imager.cloudtype = get_data_from_array(ctype.cloudtype, row_col)
    if extract_cma and cma is not None:
        obt.imager.cloudmask = get_data_from_array(cma.cma_ext, row_col)
        obt.imager.cloudmask_bin = get_data_from_array(cma.cma_bin, row_col)
        if find_mean_data_for_x_neighbours:
            obt.imager.cfc_mean = get_mean_data_from_array_nneigh(cma.cma_bin, row_col_nneigh)
    for varname in ['cma_testlist0', 'cma_testlist1', 'cma_testlist2',
                    'cma_testlist3', 'cma_testlist4', 'cma_testlist5',
                    'cma_prob', 'cma_aerosolflag', 'cma_dust', 'cma_quality']:
        if extract_cma and hasattr(cma, varname):
            setattr(obt.imager, varname,
                    get_data_from_array(getattr(cma, varname), row_col))
            if find_mean_data_for_x_neighbours and varname == 'cma_prob':
                obt.imager.cma_prob_mean = get_mean_data_from_array_nneigh(cma.cma_prob, row_col_nneigh)

    # cloud-type flags
    if extract_ctype and SETTINGS["PPS_VALIDATION"]:
        for (variable, outname) in zip(
                ['ct_quality', 'ct_conditions', 'ct_statusflag',
                 'qualityflag', 'phaseflag'],
                ['cloudtype_quality', 'cloudtype_conditions', 'cloudtype_status',
                 'cloudtype_qflag', 'cloudtype_pflag']):
            if hasattr(ctype, variable):
                setattr(obt.imager, outname,
                        get_data_from_array(getattr(ctype, variable), row_col))
    for nwp_info in aux_params:
        if extract_aux and hasattr(aux_obj, nwp_info):
            data = getattr(aux_obj, nwp_info)
            if data is not None:
                setattr(obt.imager, nwp_info, get_data_from_array(data, row_col))
        else:
            logger.debug("missing {:s}".format(nwp_info))
    CTTH_TYPES = SETTINGS["CTTH_TYPES"]
    if len(CTTH_TYPES) > 1 and SETTINGS["PPS_VALIDATION"]:
        for ctth_type in CTTH_TYPES[1:]:
            if hasattr(aux_obj, ctth_type):
                ctth_obj = getattr(aux_obj, ctth_type)
            else:
                continue
            for data_set in ["pressure", "temperature", "height"]:
                data = getattr(ctth_obj, data_set)
                name = "%s_%s" % (ctth_type.lower(), data_set)
                setattr(obt.imager, name, get_data_from_array(data, row_col))
    from atrain_match.utils.pps_prototyping_util import (get_coldest_values,
                                                         get_darkest_values,
                                                         get_warmest_values)
    if imager_obj is not None:
        pass
        # aux_obj = get_t11t12_texture_data_from_object(imager_obj, aux_obj, '11', '12',
        #                                              'text_t11t12_square??')
    for texture in ["text_r06", "text_t11", "text_t37", "text_t37t12",
                    "text_t37t12_square", "text_t11t12_square", "text_t11t12"]:
        if hasattr(aux_obj, texture):
            data = getattr(aux_obj, texture)
            setattr(obt.imager, texture, get_data_from_array(data, row_col))
    if imager_obj is not None and SETTINGS["SAVE_NEIGHBOUR_INFO"]:
        neighbour_obj = get_warmest_values(imager_obj, row_col)
        for key in ["warmest_r06", "warmest_r09", "warmest_r16"]:
            setattr(obt.imager, key + neighbour_obj.extra_info_sza_corr,
                    getattr(neighbour_obj, key))
        for key in ["warmest_t11", "warmest_t12", "warmest_t37"]:
            setattr(obt.imager, key,
                    getattr(neighbour_obj, key))
        neighbour_obj = get_darkest_values(imager_obj, row_col)
        for key in ["darkest_r06", "darkest_r09", "darkest_r16"]:
            setattr(obt.imager, key + neighbour_obj.extra_info_sza_corr,
                    getattr(neighbour_obj, key))
        for key in ["darkest_t11", "darkest_t12", "darkest_t37"]:
            setattr(obt.imager, key,
                    getattr(neighbour_obj, key))
        neighbour_obj = get_coldest_values(imager_obj, row_col)
        for key in ["coldest_r06", "coldest_r09", "coldest_r16"]:
            setattr(obt.imager, key + neighbour_obj.extra_info_sza_corr,
                    getattr(neighbour_obj, key))
        for key in ["coldest_t11", "coldest_t12", "coldest_t37"]:
            setattr(obt.imager, key,
                    getattr(neighbour_obj, key))
    # Thresholds:
    for thr in ["thr_t11ts_inv",
                "thr_t85t11_inv",
                "thr_t11t37_inv",
                "thr_t37t12_inv",
                "thr_t11t12_inv",
                "thr_t85t11",
                "thr_t11ts",
                "thr_t11t37",
                "thr_t37t12",
                "thr_t11t12",
                "thr_r06",
                "thr_r09"]:
        if hasattr(aux_obj, thr):
            data = getattr(aux_obj, thr)
            setattr(obt.imager, thr, get_data_from_array(data, row_col))
    for emis in ["emis1", "emis6", "emis8", "emis9"]:
        if hasattr(aux_obj, emis):
            data = getattr(aux_obj, emis)
            setattr(obt.imager, emis, get_data_from_array(data, row_col))
    if imager_obj is not None:
        temp_data, info = get_channel_data_from_object(imager_obj, '06', row_col)
        setattr(obt.imager, "r06micron" + info, temp_data)
        # r09
        temp_data, info = get_channel_data_from_object(imager_obj, '09', row_col)
        setattr(obt.imager, "r09micron" + info, temp_data)
        # bt37
        temp_data, info = get_channel_data_from_object(imager_obj, '37', row_col)
        setattr(obt.imager, "bt37micron" + info, temp_data)
        # b11
        temp_data, info = get_channel_data_from_object(imager_obj, '11', row_col)
        setattr(obt.imager, "bt11micron" + info, temp_data)
        # b12
        temp_data, info = get_channel_data_from_object(imager_obj, '12', row_col)
        setattr(obt.imager, "bt12micron" + info, temp_data)
        # b86
        temp_data, info = get_channel_data_from_object(imager_obj, '86', row_col)
        setattr(obt.imager, "bt86micron" + info, temp_data)
        # b16
        temp_data, info = get_channel_data_from_object(imager_obj, '16', row_col)
        setattr(obt.imager, "r16micron" + info, temp_data)
        # b22
        temp_data, info = get_channel_data_from_object(imager_obj, '22', row_col)
        setattr(obt.imager, "r22micron" + info, temp_data)
        # b13
        temp_data, info = get_channel_data_from_object(imager_obj, '13', row_col)
        setattr(obt.imager, "r13micron" + info, temp_data)
        if obt.imager_instrument.lower() in ['modis']:
            for modis_channel in CURRENTLY_UNUSED_MODIS_CHANNELS:
                modis_track, info = get_channel_data_from_object(imager_obj,
                                                                 modis_channel, row_col)
                setattr(obt.imager, modis_channel + info, modis_track)
        if obt.imager_instrument.lower() in ['seviri']:
            for seviri_channel in CURRENTLY_UNUSED_SEVIRI_CHANNELS:
                seviri_track, info = get_channel_data_from_object(imager_obj,
                                                                  seviri_channel, row_col)
                setattr(obt.imager, seviri_channel + info, seviri_track)
    # Angles, scale with gain and intercept when reading
    for angle in ['satz', 'sunz', 'azidiff', 'sunazimuth', 'satazimuth']:
        data = getattr(angle_obj, angle)
        if data is not None:
            setattr(obt.imager, angle, get_data_from_array(data.data, row_col))
    if ctth is None:
        logger.info("Not extracting ctth")
    else:
        logger.debug("Extracting ctth along track ")
        if hasattr(ctth, 'ctth_statusflag') and SETTINGS["PPS_VALIDATION"]:
            obt.imager.ctth_status = get_data_from_array(ctth.ctth_statusflag, row_col)
        for ctth_product in ['height', 'temperature', 'pressure', 'height_corr']:
            data = getattr(ctth, ctth_product)
            if data is None:
                continue
            setattr(obt.imager, "ctth_" + ctth_product, get_data_from_array(data, row_col))
        if (SETTINGS["PPS_VALIDATION"] and hasattr(ctth, 'processingflag')):
            is_opaque = np.bitwise_and(np.right_shift(ctth.processingflag, 2), 1)
            obt.imager.ctth_opaque = get_data_from_array(is_opaque, row_col)
    # NWP on ctth resolution
    if nwp_segments is not None:
        obt = insert_nwp_segments_data(nwp_segments, row_matched, col_matched, obt)
    if cpp is None:
        logger.debug("Not extracting cpp")
    elif extract_some_data_for_x_neighbours:
        for data_set_name in cpp.__dict__.keys():
            data = getattr(cpp, data_set_name)
            if data is not None:
                setattr(obt.imager, data_set_name,
                        get_data_from_array_nneigh(data, row_col_nneigh))
        for nwp_info in ["landuse", "fractionofland"]:
            data = getattr(aux_obj, nwp_info)
            setattr(obt.imager, nwp_info, get_data_from_array_nneigh(data, row_col_nneigh))
    else:
        logger.debug("Extracting cpp along track ")
        for data_set_name in cpp.__dict__.keys():
            data = getattr(cpp, data_set_name)
            if data is not None:
                setattr(obt.imager, data_set_name,
                        get_data_from_array(data, row_col))

    obt = insert_nwp_h440_h680_data(obt)

    # ADD CNN features

    from atrain_match.utils.pps_prototyping_util import add_cnn_features
    if os.path.isfile(SETTINGS['CNN_PCKL_PATH']):
        filters_dict = add_cnn_features(cloudproducts.cnn_dict, row_col,
                                        obt.imager.latitude, obt.imager.longitude,
                                        SETTINGS)
        for filter_name in filters_dict.keys():
            setattr(obt.imager, filter_name, filters_dict[filter_name])

    return obt
