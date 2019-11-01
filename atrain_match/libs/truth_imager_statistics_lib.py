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
"""Library that perfrom statistics on one matchup file."""

from atrain_match.utils.stat_util import my_iqr
from atrain_match.utils.get_flag_info import (
    get_calipso_low_medium_high_classification,
    get_cloudsat_low_medium_high_classification,
    get_day_night_twilight_info_pps2014,
    get_day_night_twilight_info_pps2012,
    get_day_night_twilight_info_cci2014)
from atrain_match.utils.common import ProcessingError
import numpy as np
import atrain_match.config as config
import logging
logger = logging.getLogger(__name__)


def calculate_ctth_stats(val_subset, imager_ctth_m_above_seasurface, truth_sat_validation_height, imager_is_cloudy):
    """Calculate CTTH statistics."""

    imager_have_hight_for_selection = np.logical_and(
        val_subset,
        np.greater_equal(imager_ctth_m_above_seasurface, 0))
    truth_have_hight_for_selection = np.logical_and(
        val_subset,
        np.greater_equal(truth_sat_validation_height, 0))
    # validate where both have height:
    val_subset = np.logical_and(imager_have_hight_for_selection,
                                truth_have_hight_for_selection)
    # print "debug", np.sum(val_subset)
    # note how many true clouds (iss/cloudsat/caliop) had no hight for imager:
    only_truth_had_height = np.logical_and(~imager_have_hight_for_selection,
                                           truth_have_hight_for_selection)
    # print "debug", np.sum(only_truth_had_height)
    n_only_truth_had_height = np.sum(only_truth_had_height)
    n_only_truth_had_height_both_had_cloud = np.sum(np.logical_and(only_truth_had_height,
                                                                   imager_is_cloudy))
    # print "debug", np.sum(n_only_truth_had_height_both_had_cloud)
    imager_height_work = imager_ctth_m_above_seasurface[val_subset.ravel(), ...]
    truth_sat_validation_height_work = truth_sat_validation_height[val_subset.ravel(), ...]
    # import pdb; pdb.set_trace()
    corr_caliop_imager = -9.0
    bias = -9.0
    RMS_difference = -9.0
    #        RMS_difference_biascorr = -9.0
    diff_squared_biascorr = np.array([-9.0])
    MAE = -9
    n_above_250 = 0
    n_above_500 = 0
    n_above_1000 = 0
    n_above_2500 = 0
    if len(truth_sat_validation_height_work) > 0:
        if len(imager_height_work) > 20:
            corr_caliop_imager = np.corrcoef(truth_sat_validation_height_work,
                                             imager_height_work)[0, 1]
        else:
            corr_caliop_imager = -99.0
        diff = imager_height_work - truth_sat_validation_height_work
        bias = np.mean(diff)
        MAE = np.mean(np.abs(diff))
        diff_squared = diff * diff
        RMS_difference = np.sqrt(np.mean(diff_squared))
        diff_squared_biascorr = (diff - bias) * (diff - bias)
        n_above_250 = np.sum(np.abs(diff) > 250)
        n_above_500 = np.sum(np.abs(diff) > 500)
        n_above_1000 = np.sum(np.abs(diff) > 1000)
        n_above_2500 = np.sum(np.abs(diff) > 2500)


#        RMS_difference_biascorr = np.sqrt(np.mean(diff_squared_biascorr))

    # return (corr_caliop_imager, bias, RMS_difference, imager_height_work, diff_squared_biascorr)
    return "%3.2f %3.2f %3.2f %d %d %d %3.2f %d %d %d %d" % (
        corr_caliop_imager,
        bias,
        RMS_difference,
        len(imager_height_work),
        n_only_truth_had_height,
        n_only_truth_had_height_both_had_cloud,
        MAE,
        n_above_250,
        n_above_500,
        n_above_1000,
        n_above_2500,
    )


def get_subset_for_mode(match_obj, mode):
    """Find the selection of pixels for one mode (like POLAR_SNOW_FREE_ICE)."""
    match_obj_truth_sat = getattr(match_obj, match_obj.truth_sat)
    latitude_abs = np.abs(getattr(match_obj_truth_sat, 'latitude'))
    if match_obj.truth_sat.lower() in ['calipso']:
        nsidc_st = getattr(match_obj_truth_sat, 'nsidc_surface_type')
        igbp_st = getattr(match_obj_truth_sat, 'igbp_surface_type')
    else:
        nsidc_st = None
        igbp_st = None

    # First prepare possible subsetting of CALIOP/CLOUDSAT/ISS without NSIDC
    # and IGBP surface types
    if mode == 'BASIC':
        cal_subset = np.bool_(np.ones(latitude_abs.shape))
    elif mode == 'SATZ_LOW':
        cal_subset = match_obj.imager.all_arrays['satz'] < 20
    elif mode == 'SATZ_HIGH':
        cal_subset = match_obj.imager.all_arrays['satz'] >= 20
    elif mode == 'OPTICAL_DEPTH':
        cal_subset = np.bool_(np.ones(latitude_abs.shape))
    elif mode == 'STANDARD':
        cal_subset = np.bool_(np.ones(latitude_abs.shape))
    elif mode == 'OPTICAL_DEPTH_THIN_IS_CLEAR':
        cal_subset = np.bool_(np.ones(latitude_abs.shape))
    elif mode == 'TROPIC_ZONE':
        cal_subset = latitude_abs <= 10
    elif mode == 'SUB_TROPIC_ZONE':
        cal_subset = np.logical_and((latitude_abs > 10),
                                    (latitude_abs <= 45))
    elif mode == 'HIGH-LATITUDES':
        cal_subset = np.logical_and((latitude_abs > 45),
                                    (latitude_abs <= 75))
    elif mode == 'POLAR':
        cal_subset = latitude_abs > 75
    elif nsidc_st is None and igbp_st is None:
        cal_subset = np.bool_(np.zeros(latitude_abs.shape))
        logger.warning("Will not run igbp/nsidc dependent mode: %s for %s",
                       mode, match_obj.truth_sat)
        return None
    # Then prepare possible subsetting of CALIOP datasets according to NSIDC
    # and IGBP surface types  if we have them
    elif mode == 'ICE_COVER_SEA':
        cal_subset = np.logical_and(
            np.logical_and(np.less_equal(nsidc_st, 100), np.greater(nsidc_st, 10)),
            np.equal(igbp_st, 17))
    elif mode == 'ICE_FREE_SEA':
        cal_subset = np.logical_and(np.equal(nsidc_st, 0), np.equal(igbp_st, 17))
    elif mode == 'SNOW_COVER_LAND':
        cal_subset = np.logical_and(
            np.logical_and(np.less(nsidc_st, 104), np.greater(nsidc_st, 10)),
            np.not_equal(igbp_st, 17))
        #  Notice that some uncertainty remains about the meaning of IGBP
        #  category 15 = "snow and ice". Can this possibly include also
        #  the Arctic ice sheet? We hope that it is not!
        #  However, if it is, the whole classification here might be wrong
        #  since this will affect also the definition of IGBP category 17./KG
    elif mode == 'SNOW_FREE_LAND':
        cal_subset = np.logical_and(np.equal(nsidc_st, 0),
                                    np.not_equal(igbp_st, 17))
    elif mode == 'COASTAL_ZONE':
        cal_subset = np.equal(nsidc_st, 255)
    elif mode == 'TROPIC_ZONE_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and(np.equal(nsidc_st, 0),
                                        np.not_equal(igbp_st, 17))
        cal_subset_area = latitude_abs <= 10
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'TROPIC_ZONE_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and(np.equal(nsidc_st, 0),
                                        np.equal(igbp_st, 17))
        cal_subset_area = latitude_abs <= 10
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'SUB_TROPIC_ZONE_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and((latitude_abs > 10),
                                        (latitude_abs <= 45))
        cal_subset_area = np.logical_and(np.equal(nsidc_st, 0),
                                         np.not_equal(igbp_st, 17))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'SUB_TROPIC_ZONE_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and((latitude_abs > 10),
                                        (latitude_abs <= 45))
        cal_subset_area = np.logical_and(np.equal(nsidc_st, 0),
                                         np.equal(igbp_st, 17))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)

    elif mode == 'HIGH-LATITUDES_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and(np.equal(nsidc_st, 0),
                                        np.not_equal(igbp_st, 17))
        cal_subset_area = np.logical_and((latitude_abs > 45),
                                         (latitude_abs <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'HIGH-LATITUDES_SNOW_COVER_LAND':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less(nsidc_st, 104),
                           np.greater(nsidc_st, 10)),
            np.not_equal(igbp_st, 17))
        cal_subset_area = np.logical_and((latitude_abs > 45), (latitude_abs <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'HIGH-LATITUDES_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and(np.equal(nsidc_st, 0),
                                        np.equal(igbp_st, 17))
        cal_subset_area = np.logical_and((latitude_abs > 45),
                                         (latitude_abs <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'HIGH-LATITUDES_ICE_COVER_SEA':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less_equal(nsidc_st, 100),
                           np.greater(nsidc_st, 10)),
            np.equal(igbp_st, 17))
        cal_subset_area = np.logical_and((latitude_abs > 45),
                                         (latitude_abs <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'POLAR_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and(np.equal(nsidc_st, 0),
                                        np.not_equal(igbp_st, 17))
        cal_subset_area = latitude_abs > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'POLAR_SNOW_COVER_LAND':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less(nsidc_st, 104),
                           np.greater(nsidc_st, 10)),
            np.not_equal(igbp_st, 17))
        cal_subset_area = latitude_abs > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'POLAR_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and(
            np.equal(nsidc_st, 0), np.equal(igbp_st, 17))
        cal_subset_area = latitude_abs > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'POLAR_ICE_COVER_SEA':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less_equal(nsidc_st, 100),
                           np.greater(nsidc_st, 10)),
            np.equal(igbp_st, 17))
        cal_subset_area = latitude_abs > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    else:
        raise ProcessingError('Unknown mode')

    return cal_subset


def get_day_night_info(match_obj, SETTINGS):
    daynight_flags = None
    match_obj_imager = getattr(match_obj, 'imager')  # Same as match_obj.imager
    match_obj_truth_sat = getattr(match_obj, match_obj.truth_sat)  # match_obj.calipso or match_obj.iss
    if not SETTINGS["PPS_VALIDATION"]:
        daynight_flags = get_day_night_twilight_info_cci2014(
            match_obj_imager.sunz)
    if SETTINGS["PPS_VALIDATION"] and hasattr(match_obj_imager, 'cloudtype_qflag'):
        if match_obj_imager.cloudtype_qflag is not None:
            daynight_flags = get_day_night_twilight_info_pps2012(
                match_obj_imager.cloudtype_qflag)
    if SETTINGS["PPS_VALIDATION"] and hasattr(match_obj_imager, 'cloudtype_conditions'):
        if match_obj_imager.cloudtype_conditions is not None:
            daynight_flags = get_day_night_twilight_info_pps2014(
                match_obj_imager.cloudtype_conditions)
    if SETTINGS["PPS_VALIDATION"] and daynight_flags is None:
        daynight_flags = get_day_night_twilight_info_cci2014(
            match_obj_imager.sunz)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag) = daynight_flags
    if (no_qflag.sum() + night_flag.sum() + twilight_flag.sum() +
            day_flag.sum()) != match_obj_truth_sat.longitude.size:
        raise ProcessingError("Something wrong with quality flags. It does not sum up.")
    return daynight_flags


def find_imager_clear_cloudy(match_obj, SETTINGS):
    """Find imager clear and cloudy pixels."""
    if 'SYNOP' in match_obj.truth_sat.upper():
        imager_clear = match_obj.imager.cfc_mean < SETTINGS["PPS_SYNOP_CLEAR_MAX_CFC"]
        imager_cloudy = match_obj.imager.cfc_mean >= SETTINGS["PPS_SYNOP_CLOUDY_MIN_CFC"]
    elif SETTINGS["USE_CMA_FOR_CFC_STATISTICS"]:
        imager_clear = np.logical_or(np.equal(match_obj.imager.cloudmask, 3),
                                     np.equal(match_obj.imager.cloudmask, 0))
        imager_cloudy = np.logical_or(np.equal(match_obj.imager.cloudmask, 1),
                                      np.equal(match_obj.imager.cloudmask, 2))
    elif SETTINGS["USE_CT_FOR_CFC_STATISTICS"]:
        imager_clear = np.logical_and(np.less_equal(match_obj.imager.cloudtype, 4),
                                      np.greater(match_obj.imager.cloudtype, 0))
        imager_cloudy = np.logical_and(np.greater(match_obj.imager.cloudtype, 4),
                                       np.less(match_obj.imager.cloudtype, 20))
    elif SETTINGS["USE_CMAPROB_FOR_CFC_STATISTICS"]:
        CMA_PROB_CLOUDY_LIMIT = SETTINGS["CMA_PROB_CLOUDY_LIMIT"]
        imager_clear = np.logical_and(np.less(match_obj.imager.cma_prob, CMA_PROB_CLOUDY_LIMIT),
                                      np.greater_equal(match_obj.imager.cma_prob, 0))
        imager_cloudy = np.logical_and(np.greater_equal(match_obj.imager.cma_prob, CMA_PROB_CLOUDY_LIMIT),
                                       np.less_equal(match_obj.imager.cma_prob, 100.0))

    else:
        raise ProcessingError("You need at least one of USE_*_FOR_CFC_STATISICS"
                              " in atrain_match.cfg to True.")
    return imager_clear, imager_cloudy


def find_truth_clear_cloudy(match_obj, val_subset, SETTINGS):
    """Find clear and cloudy pixels for the truth."""
    # For the combined 1km + 5km dataset cloud_fraction can only have values
    # (0.0, 0.2, 0.4, 0.6, 0.8, 1.0). So the threshold should
    # really be set to 0.4, i.e., at least two 1 km columns should be cloudy!.
    # Imager cloudy clear
    match_obj_truth_sat = getattr(match_obj, match_obj.truth_sat)
    if 'CALIPSO' in match_obj.truth_sat.upper():
        truth_clear = np.logical_and(
            np.less(match_obj_truth_sat.cloud_fraction, SETTINGS["CALIPSO_CLEAR_MAX_CFC"]), val_subset)
        truth_cloudy = np.logical_and(
            np.greater_equal(match_obj_truth_sat.cloud_fraction, SETTINGS["CALIPSO_CLOUDY_MIN_CFC"]), val_subset)
    elif 'SYNOP' in match_obj.truth_sat.upper():
        truth_clear = np.logical_and(
            np.less(match_obj_truth_sat.cloud_fraction, SETTINGS["SYNOP_CLEAR_MAX_CFC"]), val_subset)
        truth_cloudy = np.logical_and(
            np.greater_equal(match_obj_truth_sat.cloud_fraction, SETTINGS["SYNOP_CLOUDY_MIN_CFC"]), val_subset)
    else:
        truth_clear = np.logical_and(
            np.less_equal(match_obj_truth_sat.cloud_fraction, 0.5), val_subset)
        truth_cloudy = np.logical_and(
            np.greater(match_obj_truth_sat.cloud_fraction, 0.5), val_subset)
    return truth_clear, truth_cloudy


def get_lwp_diff_inner_cloudsat(match_obj, val_subset, wide_selection=False):
    """Find lwp difference between imager and CPR (CloudSat)."""
    selection = np.logical_and(match_obj.imager.cpp_lwp >= 0,
                               match_obj.cloudsat.RVOD_liq_water_path >= 0)
    selection = np.logical_and(selection, match_obj.imager.cpp_phase == 1)

    if wide_selection:
        pass
    else:
        # exclude risk for precipitation contamination
        selection = np.logical_and(selection, np.bitwise_and(
            np.right_shift(match_obj.cloudsat.RVOD_CWC_status, 2), 1) == 0)
        # exclude not seen in GEOPROF
        selection = np.logical_and(selection, np.bitwise_and(np.right_shift(
            match_obj.cloudsat.RVOD_CWC_status, 0), 10) == 0)  # clear geoprof
        # exclude cloudsat ice water path
        selection = np.logical_and(selection, match_obj.cloudsat.RVOD_ice_water_path <= 0)

    selection = np.logical_and(val_subset, selection)
    lwp_diff = match_obj.imager.cpp_lwp - match_obj.cloudsat.RVOD_liq_water_path
    lwp_diff = lwp_diff[selection]

    selection1 = np.logical_and(match_obj.imager.cpp_lwp >= 0,
                                match_obj.cloudsat.LO_RVOD_liquid_water_path >= 0)
    selection1 = np.logical_and(selection1, match_obj.imager.cpp_phase == 1)
    selection1 = np.logical_and(selection1, match_obj.cloudsat.cloud_fraction > 0)
    # selection = np.logical_and(selection, match_obj.imager.fractionofland <= 0)
    selection1 = np.logical_and(val_subset, selection1)
    lwp_diff_lo = match_obj.imager.cpp_lwp - match_obj.cloudsat.LO_RVOD_liquid_water_path
    lwp_diff_lo = lwp_diff_lo[selection1]
    return lwp_diff, lwp_diff_lo, match_obj.imager.cpp_lwp, match_obj.cloudsat.RVOD_liq_water_path, selection


def print_cpp_lwp_stats(match_obj, statfile, val_subset):
    """Print lwp statistics for AMSR-E or CPR (CloudSat) matchup file."""
    # CLOUD LWP EVALUATION
    # =======================
    # AMSR-E - IMAGER
    # OR CLOUDSAT CWC-RVOD- IMAGER
    if match_obj.imager.cpp_lwp is None:
        logger.warning("There are no cpp data.")
        return
    if "amsr" in match_obj.truth_sat:
        from atrain_match.utils.validate_lwp_util import get_lwp_diff
        lwp_diff = get_lwp_diff(match_obj, val_subset)
    elif "cloudsat" in match_obj.truth_sat:
        lwp_diff, lwp_diff_lo, dummy, dummy2, dummy3 = get_lwp_diff_inner_cloudsat(match_obj, val_subset)

    if len(lwp_diff) > 0:
        bias = np.mean(lwp_diff)
        diff_squared = lwp_diff*lwp_diff
        RMS_difference = np.sqrt(np.mean(diff_squared))
        N = len(lwp_diff)
        median = np.median(lwp_diff)
        iqr = my_iqr(lwp_diff)
    else:
        bias = -9
        diff_squared = -9
        RMS_difference = -9
        N = 0
        median = -9
        iqr = -9
    if "cloudsat" in match_obj.truth_sat and len(lwp_diff_lo) > 0:
        bias_lo = np.mean(lwp_diff_lo)
        diff_squared_lo = lwp_diff_lo*lwp_diff_lo
        RMS_difference_lo = np.sqrt(np.mean(diff_squared_lo))
        N_lo = len(lwp_diff_lo)
        median_lo = np.median(lwp_diff_lo)
        iqr_lo = my_iqr(lwp_diff_lo)
    else:
        bias_lo = -9
        diff_squared_lo = -9
        RMS_difference_lo = -9
        N_lo = 0
        iqr_lo = -9
        median_lo = -9

    statfile.write("CLOUD LWP %s-IMAGER TABLE: %3.2f %3.2f %d\n" % (
        match_obj.truth_sat.upper(), bias, RMS_difference, N))
    statfile.write("CLOUD LWP %s-IMAGER TABLE lo: %3.2f %3.2f %d\n" % (
        match_obj.truth_sat.upper(), bias_lo, RMS_difference_lo, N_lo))
    statfile.write("CLOUD LWP %s-IMAGER bias: %3.2f \n" % (
        match_obj.truth_sat.upper(), bias))
    statfile.write("CLOUD LWP %s-IMAGER median: %3.2f \n" % (
        match_obj.truth_sat.upper(), median))
    statfile.write("CLOUD LWP %s-IMAGER IQR: %3.2f \n" % (
        match_obj.truth_sat.upper(), iqr))
    statfile.write("CLOUD LWP %s-IMAGER std: %3.2f \n" % (
        match_obj.truth_sat.upper(), RMS_difference))
    statfile.write("CLOUD LWP %s-IMAGER bias lo: %3.2f \n" % (
        match_obj.truth_sat.upper(), bias_lo))
    statfile.write("CLOUD LWP %s-IMAGER median lo: %3.2f \n" % (
        match_obj.truth_sat.upper(), median_lo))
    statfile.write("CLOUD LWP %s-IMAGER IQR lo: %3.2f \n" % (
        match_obj.truth_sat.upper(), iqr_lo))
    statfile.write("CLOUD LWP %s-IMAGER std lo: %3.2f \n" % (
        match_obj.truth_sat.upper(), RMS_difference_lo))


def print_cpp_stats(match_obj, statfile, val_subset, SETTINGS):
    """Print Cpp statistics for CALIOPSO."""
    # CLOUD PHASE EVALUATION
    # =======================
    # CLOUD PHASE: CALIOP/ISS - IMAGER
    if match_obj.imager.cpp_phase is None:
        logger.warning("There are no cpp data.")
        return
    from atrain_match.utils.validate_cph_util import get_calipso_phase_inner, CALIPSO_PHASE_VALUES
    val_subset = np.logical_and(
        val_subset,
        match_obj.calipso.cloud_fraction >= SETTINGS["CALIPSO_CLOUDY_MIN_CFC"])
    cal_phase = get_calipso_phase_inner(
        match_obj.calipso.feature_classification_flags,
        max_layers=10,
        same_phase_in_top_three_lay=True)
    truth_water = np.equal(cal_phase, CALIPSO_PHASE_VALUES['water'])
    truth_ice = np.logical_or(
        np.equal(cal_phase, CALIPSO_PHASE_VALUES['ice']),
        np.equal(cal_phase, CALIPSO_PHASE_VALUES['horizontal_oriented_ice']))
    truth_water = np.logical_and(truth_water.data, ~cal_phase.mask)
    truth_ice = np.logical_and(truth_ice.data, ~cal_phase.mask)
    pps_water = np.equal(match_obj.imager.cpp_phase, 1)
    pps_ice = np.equal(match_obj.imager.cpp_phase, 2)
    pps_ice = np.logical_and(pps_ice, val_subset)
    pps_water = np.logical_and(pps_water, val_subset)

    n_ice_ice = np.sum(
        np.logical_and(truth_ice, pps_ice))
    n_water_water = np.sum(
        np.logical_and(truth_water, pps_water))
    n_ice_water = np.sum(
        np.logical_and(truth_ice, pps_water))
    n_water_ice = np.sum(
        np.logical_and(truth_water, pps_ice))

    nice = n_ice_ice + n_ice_water
    nwater = n_water_water + n_water_ice
    # nwater_pps = n_water_water+n_ice_water
    # nice_pps = n_water_ice+n_ice_ice

    pod_water = -9.0
    pod_ice = -9.0
    hitrate = -9
    if nwater > 0:
        pod_water = 100*float(n_water_water)/nwater
    if nice > 0:
        pod_ice = 100*float(n_ice_ice)/nice
    if nice + nwater > 0:
        hitrate = (n_ice_ice + n_water_water)*1.0/(nice+nwater)
    statfile.write("CLOUD PHASE %s-IMAGER TABLE: %s %s %s %s \n" %
                   (match_obj.truth_sat.upper(), n_ice_ice, n_ice_water, n_water_ice, n_water_water))
    statfile.write("CLOUD PHASE %s-IMAGER POD-WATER: %3.2f \n" % (match_obj.truth_sat.upper(), pod_water))
    statfile.write("CLOUD PHASE %s-IMAGER POD-ICE: %3.2f \n" % (match_obj.truth_sat.upper(), pod_ice))
    statfile.write("CLOUD PHASE %s-IMAGER Hitrate: %3.2f \n" % (match_obj.truth_sat.upper(), hitrate))


def print_cmask_stats(match_obj, statfile, val_subset, SETTINGS):
    """Print cloudmask statistics."""
    # CLOUD MASK EVALUATION
    # =======================
    # CORRELATION CLOUD MASK: CALIOP/ISS - IMAGER
    truth_clear, truth_cloudy = find_truth_clear_cloudy(match_obj, val_subset, SETTINGS)
    pps_clear, pps_cloudy = find_imager_clear_cloudy(match_obj, SETTINGS)
    pps_clear = np.logical_and(pps_clear, val_subset)
    pps_cloudy = np.logical_and(pps_cloudy, val_subset)
    n_clear_clear = np.repeat(
        pps_clear, np.logical_and(truth_clear, pps_clear)).shape[0]
    n_cloudy_cloudy = np.repeat(
        pps_cloudy, np.logical_and(truth_cloudy, pps_cloudy)).shape[0]
    n_clear_cloudy = np.repeat(
        pps_cloudy, np.logical_and(truth_clear, pps_cloudy)).shape[0]
    n_cloudy_clear = np.repeat(
        pps_clear, np.logical_and(truth_cloudy, pps_clear)).shape[0]
    nclear = n_clear_clear+n_clear_cloudy  # np.repeat(truth_clear, truth_clear).shape[0]
    ncloudy = n_cloudy_cloudy+n_cloudy_clear  # np.repeat(truth_cloudy, truth_cloudy).shape[0]
    ncloudy_pps = n_cloudy_cloudy+n_clear_cloudy
    nclear_pps = n_cloudy_clear+n_clear_clear

    pod_cloudy = -9.0*0.01
    far_cloudy = -9.0*0.01
    pod_clear = -9.0*0.01
    far_clear = -9.0*0.01
    if ncloudy > 0:
        pod_cloudy = float(n_cloudy_cloudy)/ncloudy
    if ncloudy_pps > 0:
        far_cloudy = float(n_clear_cloudy)/ncloudy_pps
    if nclear > 0:
        pod_clear = float(n_clear_clear)/nclear
    if nclear_pps > 0:
        far_clear = float(n_cloudy_clear)/nclear_pps

    if (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy) > 0:
        mean_caliop = ((n_clear_clear+n_clear_cloudy)*0.0 + (n_cloudy_clear+n_cloudy_cloudy)*1.0) / \
            (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        mean_pps = ((n_clear_clear+n_cloudy_clear)*0.0 + (n_cloudy_cloudy+n_clear_cloudy)*1.0) / \
            (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        bias = mean_pps-mean_caliop
    else:
        bias = -9.0*0.01

    statfile.write("CLOUD MASK %s-IMAGER TABLE: %s %s %s %s \n" %
                   (match_obj.truth_sat.upper(), n_clear_clear, n_clear_cloudy, n_cloudy_clear, n_cloudy_cloudy))
    statfile.write("CLOUD MASK %s-IMAGER POD-CLOUDY: %3.2f \n" % (match_obj.truth_sat.upper(), pod_cloudy*100))
    statfile.write("CLOUD MASK %s-IMAGER POD-CLEAR:  %3.2f \n" % (match_obj.truth_sat.upper(), pod_clear*100))
    statfile.write("CLOUD MASK %s-IMAGER FAR-CLOUDY: %3.2f \n" % (match_obj.truth_sat.upper(), far_cloudy*100))
    statfile.write("CLOUD MASK %s-IMAGER FAR-CLEAR:  %3.2f \n" % (match_obj.truth_sat.upper(), far_clear*100))
    statfile.write("CLOUD MASK %s-IMAGER BIAS percent: %3.2f \n" % (match_obj.truth_sat.upper(), bias*100))


def print_cmask_prob_stats(match_obj, statfile, val_subset, SETTINGS):
    """Print cloud probability statistics."""
    # CLOUD MASK PROB EVALUATION
    # =======================
    # CORRELATION CLOUD MASK: CALIOP/ISS - IMAGER
    if match_obj.imager.cma_prob is None:
        return
    if 'SYNOP' in match_obj.truth_sat.upper():
        cma_prob = match_obj.imager.cma_prob_mean
    else:
        cma_prob = match_obj.imager.cma_prob
    truth_clear, truth_cloudy = find_truth_clear_cloudy(match_obj, val_subset, SETTINGS)

    # selection:
    truth_clear = np.logical_and(truth_clear, val_subset)
    truth_cloudy = np.logical_and(truth_cloudy, val_subset)
    step = 5  # percents
    clear_string = ""
    cloudy_string = ""
    for lower in range(0, 100, step):
        upper = lower + step
        if upper == 100:
            upper = 101
        pps_in_interval = np.logical_and(cma_prob >= lower,
                                         cma_prob < upper)
        n_clear = np.sum(np.logical_and(truth_clear, pps_in_interval))
        n_cloudy = np.sum(np.logical_and(truth_cloudy, pps_in_interval))
        clear_string += "%s " % (n_clear)
        cloudy_string += "%s " % (n_cloudy)
    statfile.write("CLOUD MASK PROB %s-IMAGER TABLE STEP: %s \n" % (match_obj.truth_sat.upper(), step))
    statfile.write("CLOUD MASK PROB %s-IMAGER TABLE CLEAR: %s \n" % (match_obj.truth_sat.upper(), clear_string))
    statfile.write("CLOUD MASK PROB %s-IMAGER TABLE CLOUDY: %s \n" % (match_obj.truth_sat.upper(), cloudy_string))


def print_modis_stats(match_obj, statfile, val_subset, cal_modis_cflag, SETTINGS):
    """Print CPR (CloudSat) MODIS flag statistics compared to CALIOP or CPR (CloudSat)."""
    # CORRELATION CLOUD MASK: CALIOP - MODIS
    truth_clear, truth_cloudy = find_truth_clear_cloudy(match_obj, val_subset, SETTINGS)
    if cal_modis_cflag is None:
        return
    if len(val_subset) != len(cal_modis_cflag):
        logger.error("Lenght mismatch error for cal_MODIS_cflag")
        return

    modis_clear = np.logical_and(
        np.logical_or(np.equal(cal_modis_cflag, 1),
                      np.equal(cal_modis_cflag, 0)), val_subset)
    modis_cloudy = np.logical_and(
        np.logical_or(np.equal(cal_modis_cflag, 3),
                      np.equal(cal_modis_cflag, 2)), val_subset)

    n_clear_clear = np.repeat(
        modis_clear,
        np.logical_and(truth_clear, modis_clear)).shape[0]
    n_cloudy_cloudy = np.repeat(
        modis_cloudy,
        np.logical_and(truth_cloudy, modis_cloudy)).shape[0]
    n_clear_cloudy = np.repeat(
        modis_cloudy,
        np.logical_and(truth_clear, modis_cloudy)).shape[0]
    n_cloudy_clear = np.repeat(
        modis_clear,
        np.logical_and(truth_cloudy, modis_clear)).shape[0]
    nclear = np.repeat(truth_clear, truth_clear).shape[0]
    ncloudy = np.repeat(truth_cloudy, truth_cloudy).shape[0]
    ncloudy_modis = n_cloudy_cloudy+n_clear_cloudy
    nclear_modis = n_cloudy_clear+n_clear_clear

    if ncloudy > 0:
        pod_cloudy = float(n_cloudy_cloudy)/ncloudy
    else:
        pod_cloudy = -9.0
    if ncloudy_modis > 0:
        far_cloudy = float(n_clear_cloudy)/ncloudy_modis
    else:
        far_cloudy = -9.0
    if nclear > 0:
        pod_clear = float(n_clear_clear)/nclear
    else:
        pod_clear = -9.0
    if nclear_modis > 0:
        far_clear = float(n_cloudy_clear)/nclear_modis
    else:
        far_clear = -9.0

    if (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy) > 0:
        mean_caliop = ((n_clear_clear+n_clear_cloudy)*0.0 + (n_cloudy_clear+n_cloudy_cloudy)*1.0) / \
            (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        mean_modis = ((n_clear_clear+n_cloudy_clear)*0.0 + (n_cloudy_cloudy+n_clear_cloudy)*1.0) / \
            (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        bias = mean_modis-mean_caliop
    else:
        bias = -9.0
    statfile.write("CLOUD MASK %s-MODIS TABLE: %s %s %s %s \n" %
                   (match_obj.truth_sat.upper(), n_clear_clear, n_clear_cloudy, n_cloudy_clear, n_cloudy_cloudy))

    statfile.write("CLOUD MASK %s-MODIS FROM CLOUDSAT FLAG POD-CLOUDY:  %3.2f \n" %
                   (match_obj.truth_sat.upper(), pod_cloudy*100))
    statfile.write("CLOUD MASK %s-MODIS FROM CLOUDSAT FLAG POD-CLEAR:   %3.2f \n" %
                   (match_obj.truth_sat.upper(), pod_clear*100))
    statfile.write("CLOUD MASK %s-MODIS FROM CLOUDSAT FLAG FAR-CLOUDY:  %3.2f \n" %
                   (match_obj.truth_sat.upper(), far_cloudy*100))
    statfile.write("CLOUD MASK %s-MODIS FROM CLOUDSAT FLAG FAR-CLEAR:   %3.2f \n" %
                   (match_obj.truth_sat.upper(), far_clear*100))
    statfile.write("CLOUD MASK %s-MODIS FROM CLOUDSAT FLAG BIAS percent: %3.2f \n" %
                   (match_obj.truth_sat.upper(), bias*100))


def print_calipso_stats_ctype(match_calipso, statfile, val_subset, low_medium_high_class, SETTINGS):
    """Print cloud type statistics (CALIOP)."""
    if SETTINGS["CCI_CLOUD_VALIDATION"]:
        logger.info("Cloudtype validation not useful for CCI validation")
        return
    if match_calipso.imager.cloudtype is None:
        logger.warning("There are no cloudtype data.")
        return
    # CLOUD TYPE EVALUATION - Based exclusively on CALIPSO data (Vertical Feature Mask)
    # =======================
    calipso_low = np.logical_and(low_medium_high_class['low_clouds'],
                                 val_subset)
    calipso_medium = np.logical_and(low_medium_high_class['medium_clouds'],
                                    val_subset)
    calipso_high = np.logical_and(low_medium_high_class['high_clouds'],
                                  val_subset)
    calipso_medium_tp = np.logical_and(low_medium_high_class['medium_clouds_tp'],
                                       val_subset)
    calipso_high_tp = np.logical_and(low_medium_high_class['high_clouds_tp'],
                                     val_subset)
    calipso_medium_op = np.logical_and(low_medium_high_class['medium_clouds_op'],
                                       val_subset)
    calipso_high_op = np.logical_and(low_medium_high_class['high_clouds_op'],
                                     val_subset)

    if match_calipso.imager.cloudtype_conditions is not None:
        logger.debug("Assuming cloudtype structure from pps v2014")
        imager_low = np.logical_and(
            np.logical_and(np.greater_equal(match_calipso.imager.cloudtype, 5),
                           np.less_equal(match_calipso.imager.cloudtype, 6)),
            val_subset)
        imager_medium = np.logical_and(
            np.equal(match_calipso.imager.cloudtype, 7), val_subset)
        imager_high_op = np.logical_and(
            np.logical_and(np.greater_equal(match_calipso.imager.cloudtype, 8),
                           np.less_equal(match_calipso.imager.cloudtype, 9)),
            val_subset)
        imager_cirrus = np.logical_and(
            np.logical_and(np.greater_equal(match_calipso.imager.cloudtype, 11),
                           np.less_equal(match_calipso.imager.cloudtype, 15)),
            val_subset)
        imager_high = imager_high_op  # np.logical_or(imager_high_op, imager_high_semi)
        imager_frac = np.logical_and(np.equal(match_calipso.imager.cloudtype, 10),
                                     val_subset)
        imager_low = np.logical_or(imager_low, imager_frac)

    else:
        logger.error("Assuming cloudtype structure from pps v2012")

    calipso_clear = np.logical_and(
        np.less(match_calipso.calipso.cloud_fraction, SETTINGS["CALIPSO_CLEAR_MAX_CFC"]), val_subset)
    calipso_cloudy = np.logical_and(
        np.greater_equal(match_calipso.calipso.cloud_fraction, SETTINGS["CALIPSO_CLOUDY_MIN_CFC"]), val_subset)
    imager_clear = np.logical_and(
        np.logical_and(np.less_equal(match_calipso.imager.cloudtype, 4),
                       np.greater(match_calipso.imager.cloudtype, 0)),
        val_subset)

    # Notice that we have unfortunately changed order in notation compared to cloud mask
    # Here the PPS category is mentioned first and then the CALIOP category

    n_low_low = np.repeat(
        imager_low,
        np.logical_and(calipso_low, imager_low)).shape[0]
    n_low_medium = np.repeat(
        imager_low,
        np.logical_and(calipso_medium, imager_low)).shape[0]
    n_low_high = np.repeat(
        imager_low,
        np.logical_and(calipso_high, imager_low)).shape[0]
    n_medium_low = np.repeat(
        imager_medium,
        np.logical_and(calipso_low, imager_medium)).shape[0]
    n_medium_medium = np.repeat(
        imager_medium,
        np.logical_and(calipso_medium, imager_medium)).shape[0]
    n_medium_high = np.repeat(
        imager_medium,
        np.logical_and(calipso_high, imager_medium)).shape[0]
    n_high_low = np.repeat(
        imager_high,
        np.logical_and(calipso_low, imager_high)).shape[0]
    n_high_medium = np.repeat(
        imager_high,
        np.logical_and(calipso_medium, imager_high)).shape[0]
    n_high_high = np.repeat(
        imager_high,
        np.logical_and(calipso_high, imager_high)).shape[0]
    n_cirrus_low = np.repeat(
        imager_cirrus,
        np.logical_and(calipso_low, imager_cirrus)).shape[0]
    n_cirrus_medium_tp = np.repeat(
        imager_cirrus,
        np.logical_and(calipso_medium_tp, imager_cirrus)).shape[0]
    n_cirrus_high_tp = np.repeat(
        imager_cirrus,
        np.logical_and(calipso_high_tp, imager_cirrus)).shape[0]
    n_cirrus_medium_op = np.repeat(
        imager_cirrus,
        np.logical_and(calipso_medium_op, imager_cirrus)).shape[0]
    n_cirrus_high_op = np.repeat(
        imager_cirrus,
        np.logical_and(calipso_high_op, imager_cirrus)).shape[0]
    n_clear_low = np.repeat(
        imager_clear,
        np.logical_and(calipso_low, imager_clear)).shape[0]
    n_clear_medium = np.repeat(
        imager_clear,
        np.logical_and(calipso_medium, imager_clear)).shape[0]
    n_clear_high = np.repeat(
        imager_clear,
        np.logical_and(calipso_high, imager_clear)).shape[0]
    n_low_clear = np.repeat(
        imager_low,
        np.logical_and(calipso_clear, imager_low)).shape[0]
    n_medium_clear = np.repeat(
        imager_medium,
        np.logical_and(calipso_clear, imager_medium)).shape[0]
    n_high_clear = np.repeat(
        imager_high,
        np.logical_and(calipso_clear, imager_high)).shape[0]
    n_cirrus_clear = np.repeat(
        imager_cirrus,
        np.logical_and(calipso_clear, imager_cirrus)).shape[0]

    pod_low = -9.0
    far_low = -9.0
    pod_medium = -9.0
    far_medium = -9.0
    pod_high = -9.0
    far_high = -9.0
    far_cirrus = -9.0
    N_cal_low = n_low_low + n_medium_low + n_high_low + n_cirrus_low
    N_pps_low = (n_low_low + n_low_medium+n_low_high)
    N_cal_medium = n_low_medium+n_medium_medium+n_high_medium + n_cirrus_medium_tp + n_cirrus_medium_op
    N_pps_medium = n_medium_low + n_medium_medium+n_medium_high
    N_cal_high = n_low_high+n_medium_high+n_high_high + n_cirrus_high_op + n_cirrus_high_tp
    N_pps_high = n_high_low + n_high_medium + n_high_high
    N_pps_cirrus = n_cirrus_high_op + n_cirrus_high_tp + n_cirrus_low + n_cirrus_medium_op + n_cirrus_medium_tp
    if (N_cal_low) > 0:
        pod_low = float(n_low_low) / N_cal_low
    if (N_pps_low) > 0:
        far_low = float(N_pps_low - n_low_low) / N_pps_low

    if N_cal_medium > 0:
        pod_medium = float(n_medium_medium + n_cirrus_medium_tp)/N_cal_medium
    if N_pps_medium > 0:
        far_medium = float(n_medium_low+n_medium_high)/N_pps_medium

    if N_cal_high > 0:
        pod_high = float(n_high_high+n_cirrus_high_tp)/N_cal_high
    if N_pps_high > 0:
        far_high = float(n_high_low+n_high_medium)/N_pps_high
    if N_pps_cirrus > 0:
        far_cirrus = float(n_cirrus_low+n_cirrus_medium_op + n_cirrus_high_op)/N_pps_cirrus

    statfile.write("CLOUD TYPE %s-IMAGER TABLE: %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % (
        match_calipso.truth_sat.upper(),
        n_low_low, n_low_medium, n_low_high,
        n_medium_low, n_medium_medium, n_medium_high,
        n_high_low, n_high_medium, n_high_high,
        n_cirrus_low,
        n_cirrus_medium_tp, n_cirrus_high_tp,
        n_cirrus_medium_op, n_cirrus_high_op))
    statfile.write("CLOUD TYPE %s-IMAGER PROB: %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f \n" % (
        match_calipso.truth_sat.upper(),
        pod_low, pod_medium, pod_high, far_low, far_medium, far_high, far_cirrus))
    statfile.write("CLOUD TYPE %s-IMAGER TABLE MISSED: %s %s %s %s %s %s %s \n" % (
        match_calipso.truth_sat.upper(),
        n_clear_low, n_clear_medium, n_clear_high,
        n_low_clear, n_medium_clear, n_high_clear, n_cirrus_clear))


def print_height_all_low_medium_high(NAME, val_subset, statfile,
                                     low_medium_high_class, imager_ctth_m_above_seasurface,
                                     truth_sat_validation_height, imager_is_cloudy):
    """Print all CTTH statistics for one case (thin/thick/geostyle etc.)."""
    out_stats = calculate_ctth_stats(val_subset, imager_ctth_m_above_seasurface,
                                     truth_sat_validation_height, imager_is_cloudy)
    statfile.write("CLOUD HEIGHT %s ALL: %s\n" % (NAME, out_stats))
    if low_medium_high_class is None:
        # Nothing more can be done!
        return
    cal_low_ok = np.logical_and(low_medium_high_class['low_clouds'],
                                val_subset)
    out_stats = calculate_ctth_stats(cal_low_ok, imager_ctth_m_above_seasurface,
                                     truth_sat_validation_height, imager_is_cloudy)
    statfile.write("CLOUD HEIGHT %s LOW: %s \n" % (NAME, out_stats))
    cal_mid_ok = np.logical_and(low_medium_high_class['medium_clouds'],
                                val_subset)
    out_stats = calculate_ctth_stats(cal_mid_ok, imager_ctth_m_above_seasurface,
                                     truth_sat_validation_height, imager_is_cloudy)
    statfile.write("CLOUD HEIGHT %s MEDIUM: %s \n" % (NAME, out_stats))
    cal_high_ok = np.logical_and(low_medium_high_class['high_clouds'],
                                 val_subset)
    out_stats = calculate_ctth_stats(cal_high_ok, imager_ctth_m_above_seasurface,
                                     truth_sat_validation_height, imager_is_cloudy)
    statfile.write("CLOUD HEIGHT %s HIGH: %s \n" % (NAME, out_stats))


def print_stats_ctop(match_obj, statfile, val_subset, low_medium_high_class, SETTINGS):
    """Print all CTTH statistics."""
    if match_obj.imager.ctth_height is None:
        logger.warning("There are no ctth height data.")
        return

    # CORRELATION: CALIOP - IMAGER HEIGHT
    # FIRST TOTAL FIGURES

    match_obj_imager = getattr(match_obj, 'imager')  # Same as match_obj.imager
    match_obj_truth_sat = getattr(match_obj, match_obj.truth_sat)  # match_obj.calipso or match_obj.iss
    imager_ctth_m_above_seasurface = match_obj_imager.imager_ctth_m_above_seasurface
    logger.warning("WARNING Only validating CTTH for cloudy pixels!")
    (dummy, imager_is_cloudy) = find_imager_clear_cloudy(match_obj, SETTINGS)
    imager_ctth_m_above_seasurface[~imager_is_cloudy] = -9
    # imager_ctth_m_above_seasurface[match_obj_imager.cloudmask==0] = -9
    # imager_ctth_m_above_seasurface[match_obj_imager.cloudmask==3] = -9
    truth_sat_validation_height = match_obj_truth_sat.validation_height
    truth_sat_validation_height = np.array(match_obj_truth_sat.validation_height, dtype=np.float64)

    val_subset = np.logical_and(
        val_subset,
        match_obj_truth_sat.cloud_fraction >= SETTINGS["CALIPSO_CLOUDY_MIN_CFC"])

    # print "ALL CLOUDS:"
    print_height_all_low_medium_high(match_obj.truth_sat.upper(),
                                     val_subset,
                                     statfile, low_medium_high_class,
                                     imager_ctth_m_above_seasurface,
                                     truth_sat_validation_height,
                                     imager_is_cloudy)

    if match_obj.truth_sat.upper() not in ["CALIPSO"]:
        if match_obj.truth_sat.upper() in ["ISS"]:
            logger.warning("WARNING WARNING WARNING only printing over all statistics "
                           "for cloudtop for ISS")
        return
    if (SETTINGS["COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC"] and
            match_obj.imager.cloudtype is not None):
        statfile.write("CLOUD HEIGHT GEO-STYLE\n")
        from scipy import ndimage
        # GEO uses pixels with homogene CT in 9x9 pixels
        # Let us use 9 homogene CT pixels. Not fully the same, but similar.
        # And variation (definition?) caliop pressure <200hPa
        # And variation CPR height less than 3km
        # And CALIPO clouds thinner than 0.2 removed. however pixel kept
        # For 1km data we have to either keep or fully remove the pixel
        maxct = ndimage.filters.maximum_filter1d(match_obj.imager.cloudtype, size=9)
        minct = ndimage.filters.minimum_filter1d(match_obj.imager.cloudtype, size=9)
        val_geo = np.logical_and(
            val_subset,
            np.equal(maxct, minct))
        if hasattr(match_obj, 'calipso'):
            var_pressure = (ndimage.filters.maximum_filter1d(match_obj.calipso.layer_top_pressure[:, 0], size=9) -
                            ndimage.filters.minimum_filter1d(match_obj.calipso.layer_top_pressure[:, 0], size=9))
            val_geo = np.logical_and(
                val_geo,
                var_pressure < 200)  # Pressure variation less than 200hPa
        if hasattr(match_obj, 'cloudsat'):
            var_height = (ndimage.filters.maximum_filter1d(truth_sat_validation_height, size=9) -
                          ndimage.filters.minimum_filter1d(truth_sat_validation_height, size=9))
            val_geo = np.logical_and(
                val_geo,
                var_height < 3000)  # Height variation less than 3km
        average_height_truth = ndimage.filters.uniform_filter1d(truth_sat_validation_height*1.0, size=9)
        average_height_truth[truth_sat_validation_height < 0] = -9
        average_height_imager = ndimage.filters.uniform_filter1d(imager_ctth_m_above_seasurface*1.0, size=9)
        average_height_imager[imager_ctth_m_above_seasurface < 0] = -9
        print_height_all_low_medium_high("CALIOP-GEO-STYLE",
                                         val_geo,
                                         statfile, low_medium_high_class,
                                         average_height_imager,
                                         average_height_truth,
                                         imager_is_cloudy)
        val_geo = np.logical_and(
            val_geo,
            np.greater_equal(match_obj.calipso.feature_optical_depth_532_top_layer_5km, 0.2))
        statfile.write("CLOUD HEIGHT GEO-STYLE-EXCLUDE-THIN-PIXELS\n")
        print_height_all_low_medium_high("CALIOP-GEO-STYLE-EXCLUDE-THIN-PIXELS",
                                         val_geo,
                                         statfile, low_medium_high_class,
                                         average_height_imager,
                                         average_height_truth,
                                         imager_is_cloudy)

    if SETTINGS["COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC"]:
        statfile.write("CLOUD HEIGHT SINGLE-LAYER\n")
        val_subset_single = np.logical_and(
            val_subset,
            np.equal(match_obj.calipso.number_layers_found, 1))
        print_height_all_low_medium_high("CALIOP-SINGLE-LAYER",
                                         val_subset_single,
                                         statfile, low_medium_high_class,
                                         imager_ctth_m_above_seasurface,
                                         truth_sat_validation_height,
                                         imager_is_cloudy)

    if (SETTINGS["COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC"] and
            (SETTINGS["ALSO_USE_5KM_FILES"] or config.RESOLUTION == 5)):
        statfile.write("CLOUD HEIGHT SINGLE-LAYER, NOT THIN\n")
        lim = 2*SETTINGS["OPTICAL_DETECTION_LIMIT"]
        val_subset_single_not_thinnest = np.logical_and(
            val_subset_single,
            np.greater_equal(match_obj.calipso.feature_optical_depth_532_top_layer_5km, lim))
        print_height_all_low_medium_high("CALIOP-SINGLE-LAYER>%f" % (lim),
                                         val_subset_single_not_thinnest,
                                         statfile, low_medium_high_class,
                                         imager_ctth_m_above_seasurface,
                                         truth_sat_validation_height,
                                         imager_is_cloudy)

        statfile.write("CLOUD HEIGHT NOT VERY THIN TOP LAYER\n")
        lim = SETTINGS["OPTICAL_DETECTION_LIMIT"]
        val_subset_not_thinnest_top_layer = np.logical_and(
            val_subset,
            np.greater_equal(match_obj.calipso.feature_optical_depth_532_top_layer_5km, lim))
        print_height_all_low_medium_high("CALIOP-TOP-LAYER>%f" % (lim),
                                         val_subset_not_thinnest_top_layer,
                                         statfile, low_medium_high_class,
                                         imager_ctth_m_above_seasurface,
                                         truth_sat_validation_height,
                                         imager_is_cloudy)

        lim = SETTINGS["OPTICAL_DETECTION_LIMIT"]
        statfile.write("CLOUD HEIGHT VERY THIN TOP LAYER\n")
        val_subset_thinnest_top_layer = np.logical_and(
            val_subset,
            np.less_equal(match_obj.calipso.feature_optical_depth_532_top_layer_5km, lim))
        print_height_all_low_medium_high("CALIOP-TOP-LAYER<=%f" % (lim),
                                         val_subset_thinnest_top_layer,
                                         statfile, low_medium_high_class,
                                         imager_ctth_m_above_seasurface,
                                         truth_sat_validation_height,
                                         imager_is_cloudy)


def print_main_stats(match_obj, statfile):
    """Print some geolocation and time information about the match."""
    val_object = getattr(match_obj, match_obj.truth_sat)
    num_val_data_ok = len(match_obj.diff_sec_1970)
    statfile.write("%s min and max time diff: %3.2f %3.2f \n" % (
        match_obj.truth_sat.upper(),
        match_obj.diff_sec_1970.min(),
        match_obj.diff_sec_1970.max()))
    statfile.write("%s start and stop Latitude: %3.2f %3.2f \n" % (
        match_obj.truth_sat.upper(),
        val_object.latitude[0],
        val_object.latitude[-1]))
    statfile.write("%s start and stop Longitude: %3.2f %3.2f \n" % (
        match_obj.truth_sat.upper(),
        val_object.longitude[0],
        val_object.longitude[-1]))
    statfile.write("%s-IMAGER number of matches: %d\n" % (
        match_obj.truth_sat.upper(),
        num_val_data_ok))


def calculate_statistics(mode, statfilename, match_calipso, match_clsat, match_iss, match_amsr, match_synop,
                         SETTINGS,
                         dnt_flag=None):
    """Calculate all requested statistics, for all matches with one imager cloudproduct file (main function)."""

    def get_day_night_subset(match_obj, val_subset, SETTINGS):
        (no_qflag, night_flag, twilight_flag,
         day_flag, all_dnt_flag) = get_day_night_info(match_obj, SETTINGS)

        if dnt_flag is None:
            logger.debug('dnt_flag = %s', 'ALL PIXELS')
            dnt_subset = np.logical_and(val_subset, all_dnt_flag)
        elif dnt_flag.upper() == 'DAY':
            logger.debug('dnt_flag = %s', dnt_flag.upper())
            dnt_subset = np.logical_and(val_subset, day_flag)
        elif dnt_flag.upper() == 'NIGHT':
            logger.debug('dnt_flag = %s', dnt_flag.upper())
            dnt_subset = np.logical_and(val_subset, night_flag)
        elif dnt_flag.upper() == 'TWILIGHT':
            logger.debug('dnt_flag = %s', dnt_flag.upper())
            dnt_subset = np.logical_and(val_subset, twilight_flag)
        else:
            raise ProcessingError("Unknown DNT-flag %s" % (dnt_flag.upper()))
        return dnt_subset

    if match_clsat is not None:
        logger.info("Cloudsat Statistics")
        val_subset = get_subset_for_mode(match_clsat, mode)
        if val_subset is not None:
            val_subset = get_day_night_subset(match_clsat, val_subset, SETTINGS)
            statfile = open(statfilename.replace('xxx', 'cloudsat'), "w")
            if match_clsat.cloudsat.all_arrays['cloud_fraction'] is not None:
                low_medium_high_class = get_cloudsat_low_medium_high_classification(match_clsat)
                print_main_stats(match_clsat, statfile)
                print_cmask_stats(match_clsat, statfile, val_subset, SETTINGS)
                print_cmask_prob_stats(match_clsat, statfile, val_subset, SETTINGS)
                print_modis_stats(match_clsat, statfile, val_subset, match_clsat.cloudsat.MODIS_cloud_flag, SETTINGS)
                print_stats_ctop(match_clsat, statfile, val_subset, low_medium_high_class, SETTINGS)
            if match_clsat.cloudsat.all_arrays['RVOD_liq_water_path'] is not None:
                print_cpp_lwp_stats(match_clsat, statfile, val_subset)
            statfile.close()

    if match_calipso is not None:
        logger.info("Calipo Statistics")
        val_subset = get_subset_for_mode(match_calipso, mode)
        if val_subset is not None:
            statfile = open(statfilename.replace('xxx', 'calipso'), "w")
            low_medium_high_class = get_calipso_low_medium_high_classification(match_calipso)
            # semi_flag, opaque_flag = get_semi_opaque_info(match_calipso)
            val_subset = get_day_night_subset(match_calipso, val_subset, SETTINGS)
            print_main_stats(match_calipso, statfile)
            print_cmask_stats(match_calipso, statfile, val_subset, SETTINGS)
            print_cmask_prob_stats(match_calipso, statfile, val_subset, SETTINGS)
            print_modis_stats(match_calipso, statfile, val_subset, match_calipso.calipso.cal_modis_cflag, SETTINGS)
            print_calipso_stats_ctype(match_calipso, statfile, val_subset, low_medium_high_class, SETTINGS)
            print_stats_ctop(match_calipso, statfile, val_subset, low_medium_high_class, SETTINGS)
            print_cpp_stats(match_calipso, statfile, val_subset, SETTINGS)
            statfile.close()

    if match_iss is not None:
        val_subset = get_subset_for_mode(match_iss, mode)
        if val_subset is not None:
            statfile = open(statfilename.replace('xxx', 'iss'), "w")
            val_subset = get_day_night_subset(match_iss, val_subset, SETTINGS)
            print_main_stats(match_iss, statfile)
            print_cmask_stats(match_iss, statfile, val_subset, SETTINGS)
            print_cmask_prob_stats(match_iss, statfile, val_subset, SETTINGS)
            # print_calipso_stats_ctype(match_iss, statfile, val_subset, cal_vert_feature)
            print_stats_ctop(match_iss, statfile, val_subset, None, SETTINGS)
            statfile.close()

    if match_amsr is not None:
        logger.info("AMSR-E Statistics")
        val_subset = get_subset_for_mode(match_amsr, mode)
        if val_subset is not None:
            val_subset = get_day_night_subset(match_amsr, val_subset, SETTINGS)
            statfile = open(statfilename.replace('xxx', 'amsr'), "w")
            print_main_stats(match_amsr, statfile)
            print_cpp_lwp_stats(match_amsr, statfile, val_subset)
            statfile.close()

    if match_synop is not None:
        logger.info("SYNOP Statistics")
        val_subset = get_subset_for_mode(match_synop, mode)
        if val_subset is not None:
            val_subset = get_day_night_subset(match_synop, val_subset, SETTINGS)
            statfile = open(statfilename.replace('xxx', 'synop'), "w")
            print_main_stats(match_synop, statfile)
            print_cmask_stats(match_synop, statfile, val_subset, SETTINGS)
            print_cmask_prob_stats(match_synop, statfile, val_subset, SETTINGS)
            statfile.close()
