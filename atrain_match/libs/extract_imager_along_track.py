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


def get_atrain_name(ch):
    out = ch.id_tag.replace('tb', 'bt').replace('ch_', '') + 'micron'
    if ch.SZA_corr_done:
        out += "_sza_correction_done"
    return out

def get_channel_data_from_object(imager_obj, chn_des, matched, nodata=-9):
    """Get the IMAGER/VIIRS channel data on the track

    matched: dict of matched indices (row, col)

    """
    if matched is None:
        return imager_obj.channel[chn_des].data
    else:
        chdata_on_track = get_data_from_array(imager_obj.channel[chn_des].data, matched)
        return np.array(chdata_on_track)

def _interpolate_pressure_from_height(imager_obj, list_of_levels):
    """ Function to find height att pressure level (level)
    from segment_nwp, pressure and height vectors.le
    High means high in pressure. The level closest to ground i hi, and lo is at lower
    pressure further up in atmosphere.
    """
    if hasattr(imager_obj, "nwp_height") and imager_obj.nwp_height is not None:
        height_h = imager_obj.nwp_height
        pressure_v = imager_obj.nwp_pressure
        surface_h = imager_obj.nwp_surface_h
        psur = imager_obj.nwp_psur
    else:
        return None
    
    nlev = pressure_v.shape[1]
    npix = pressure_v.shape[0]
    k = np.arange(npix)

    higher_index = np.array([np.searchsorted(height_h[ind, :], list_of_levels[ind], side='right')
                             for ind in range(npix)])

    
    higher_index[higher_index >= (nlev - 1)] = nlev - 2
    lower_index = higher_index - 1

        
    # update "lo" where level is between surface and first level in array
    below_level_1 = list_of_levels < height_h[:, 0]
    lower_index[below_level_1] = 0


    # get pressure and height for layer below and above level
    hi = pressure_v[k, higher_index]
    lo = pressure_v[k, lower_index]
    height_hi_ = height_h[k, higher_index]*1.0
    height_lo_ = height_h[k, lower_index]*1.0
    # update "hi" where level is between surface and first level in array
    lo[below_level_1] = psur[below_level_1]
    height_lo_[below_level_1] = surface_h[below_level_1]
    # log pressures
    hi = np.log(hi)
    lo = np.log(lo)
    #level = np.log(level)
    # interpolate
    #out_h = height_hi_ - (hi - level) * (height_hi_ - height_lo_) / (hi - lo)
    #(hi - level) * (height_hi_ - height_lo_) / (hi - lo) = height_hi_ - out_h
    out_pressure = np.exp(hi - (height_hi_ - list_of_levels) *  (hi - lo) /  (height_hi_ - height_lo_))
    out_pressure[list_of_levels < 0] = -9

    return out_pressure


def _interpolate_height_or_temperature_from_pressure(imager_obj,
                                                     level,
                                                     list_of_levels=None,
                                                     temperature=False,):
    """ Function to find height or temperature at pressure level (level)
    from nwp, pressure and height vectors.
    High means high in pressure. The level closest to ground i hi, and lo is at lower
    pressure further up in atmosphere.
    """
    if hasattr(imager_obj, "nwp_height") and imager_obj.nwp_height is not None:
        pressure_v = imager_obj.nwp_pressure
        psur = imager_obj.nwp_psur
        if temperature:
            values_h = imager_obj.nwp_temperature
            surface_height = imager_obj.nwp_tsur
        else:
            values_h = imager_obj.nwp_height
            surface_height = imager_obj.nwp_surface_h
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
    height_hi_[below_level_1] = surface_height[below_level_1]
    # log pressures
    hi = np.log(hi)
    lo = np.log(lo)
    level = np.log(level)
    # interpolate
    out_h = height_hi_ - (hi - level) * (height_hi_ - height_lo_) / (hi - lo)
    return out_h


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
                              extract_unc=True,
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
     # For amsr-E matching (many neighbors) use only the needed nwp data (aux_params != None)
    if aux_params is None:
        # aux_params = aux_params_all
        aux_params = dir(aux_obj)
        aux_params = [param for param in aux_params if '__' not in param]

    
    imager_obj = cloudproducts.imager_channeldata
    angle_obj = cloudproducts.imager_angles
    ctth = cloudproducts.ctth
    cma = cloudproducts.cma
    ctype = cloudproducts.ctype
    cpp = cloudproducts.cpp
    unc_obj = cloudproducts.unc
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
                
    # Emiss, threhsolds, texture, nwp etc
    for key in aux_params:
        if extract_aux and hasattr(aux_obj, key):
            if 'CTTH' in key:
                ctth_type = key
                ctth_obj = getattr(aux_obj, ctth_type)
                for data_set in ["pressure", "temperature", "height"]:
                    data = getattr(ctth_obj, data_set)
                    name = "%s_%s" % (ctth_type.lower(), data_set)
                    setattr(obt.imager, name, get_data_from_array(data, row_col))
            else:    
                data = getattr(aux_obj, key)
                if data is not None:
                    setattr(obt.imager, key, get_data_from_array(data, row_col))
        else:
            logger.debug("missing {:s}".format(key))
  
    from atrain_match.utils.pps_prototyping_util import (get_coldest_values,
                                                         get_darkest_values,
                                                         get_warmest_values)

    if imager_obj is not None:
        imager_channels = [key for key in imager_obj.channel if 'ch_' in key]

    if imager_obj is not None and SETTINGS["SAVE_NEIGHBOUR_INFO"]  and extract_radiances:
        warm_row_col = get_warmest_values(imager_obj, row_col)
        for key in imager_channels:
            atrain_name = get_atrain_name(imager_obj.channel[key])
            atrain_name = 'warmest_'  +atrain_name.replace('micron',''). replace('b','')
            data = get_channel_data_from_object(imager_obj, key, warm_row_col)
            setattr(obt.imager, atrain_name, data)            
        cold_row_col = get_coldest_values(imager_obj, row_col)
        for key in imager_channels:
            atrain_name = get_atrain_name(imager_obj.channel[key])
            atrain_name = 'coldest_'  + atrain_name.replace('micron',''). replace('b','')
            data = get_channel_data_from_object(imager_obj, key, cold_row_col)
            setattr(obt.imager, atrain_name, data)
        dark_row_col = get_darkest_values(imager_obj, row_col)
        for key in imager_channels:
            atrain_name = get_atrain_name(imager_obj.channel[key])
            atrain_name = 'darkest_'  + atrain_name.replace('micron',''). replace('b','')
            data = get_channel_data_from_object(imager_obj, key, dark_row_col)
            setattr(obt.imager, atrain_name, data)

        for delta_r in range(-2,3):
            for delta_c in range(-2,3):
                if delta_r == 0 and delta_c == 0:
                    # Actual pixel
                    continue
                for key in imager_channels:
                    max_row, max_col =imager_obj.channel[key].data.shape
                    atrain_name = key.replace('micron','')
                    atrain_name = atrain_name + '_neighbour_r{:d}_{:d}'.format(delta_r, delta_c)
                    new_r = row_col['row'] + delta_r
                    new_r[new_r < 0] = 0
                    new_c = row_col['col'] + delta_c
                    new_c[new_c < 0] = 0
                    new_r[new_r >= max_row] = max_row-1
                    new_c[new_c >= max_col] = max_col-1
                    new_row_col = {'row': new_r,
                                   'col': new_c}
                    new_row_col['row']
                    data = get_channel_data_from_object(imager_obj, key, new_row_col)
                    setattr(obt.imager, atrain_name, data)

    if imager_obj is not None:            
        if 'qual_flags' in imager_obj.channel:
            qual = imager_obj.channel['qual_flags'].data
            setattr(obt.imager, 'qual_flags', qual[truth.imager_linnum,:])
        
    # Imager data        
    if imager_obj is not None and extract_radiances:
        for key in imager_channels:
            atrain_name = get_atrain_name(imager_obj.channel[key])
            data = get_data_from_array(imager_obj.channel[key].data, row_col)
            setattr(obt.imager, atrain_name, data)  
      
    # Angles, scale with gain and intercept when reading
    for angle in ['satz', 'sunz', 'azidiff', 'sunazimuth', 'satazimuth']:
        data = getattr(angle_obj, angle)
        if data is not None:
            setattr(obt.imager, angle, get_data_from_array(data.data, row_col))
    if ctth is None:
        logger.info("Not extracting ctth")
    elif extract_ctth:
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

    if cpp is None or not extract_cpp:
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

    # extract cci uncertainties if requested
    if unc_obj is not None and extract_unc:
        if extract_some_data_for_x_neighbours:
            extractor = row_col_nneigh
        else:
            extractor = row_col
            
        for data_set_name in unc_obj.__dict__.keys():
            data = getattr(unc_obj, data_set_name)
            if data is not None:
                setattr(obt.imager, data_set_name,
                        get_data_from_array(data, extractor))

    # ADD CNN features

    from atrain_match.utils.pps_prototyping_util import add_cnn_features
    if os.path.isfile(SETTINGS['CNN_PCKL_PATH']):
        filters_dict = add_cnn_features(cloudproducts.cnn_dict, row_col,
                                        obt.imager.latitude, obt.imager.longitude,
                                        SETTINGS)
        for filter_name in filters_dict.keys():
            setattr(obt.imager, filter_name, filters_dict[filter_name])

    return obt

def imager_track_from_matched_hrit(obt, SETTINGS, cloudproducts,
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


    truth = getattr(obt, obt.truth_sat)
    row_matched = truth.imager_linnum
    col_matched = truth.imager_pixnum
    row_col = {'row': row_matched, 'col': col_matched}

    obt.imager.latitude = get_data_from_array(cloudproducts.latitude, row_col)
    obt.imager.longitude = get_data_from_array(cloudproducts.longitude, row_col)
    obt.imager.vis006 = get_data_from_array(cloudproducts.vis006, row_col)
    obt.imager.vis008 = get_data_from_array(cloudproducts.vis008, row_col)
    obt.imager.ir_016 = get_data_from_array(cloudproducts.ir_016, row_col)
    obt.imager.ir_039 = get_data_from_array(cloudproducts.ir_039, row_col)
    obt.imager.ir_062 = get_data_from_array(cloudproducts.ir_062, row_col)
    obt.imager.ir_073 = get_data_from_array(cloudproducts.ir_073, row_col)
    obt.imager.ir_087 = get_data_from_array(cloudproducts.ir_087, row_col)
    obt.imager.ir_097 = get_data_from_array(cloudproducts.ir_097, row_col)
    obt.imager.ir_108 = get_data_from_array(cloudproducts.ir_108, row_col)
    obt.imager.ir_120 = get_data_from_array(cloudproducts.ir_120, row_col)
    obt.imager.ir_134 = get_data_from_array(cloudproducts.ir_134, row_col)
    obt.imager.satzen = get_data_from_array(cloudproducts.satzen, row_col)
    obt.imager.solzen = get_data_from_array(cloudproducts.solzen, row_col)
    obt.imager.satazi = get_data_from_array(cloudproducts.satazi, row_col)
    obt.imager.time = get_data_from_array(cloudproducts.time, row_col)

    return obt
