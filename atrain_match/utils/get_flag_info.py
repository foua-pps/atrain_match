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
logger = logging.getLogger(__name__)

# 0 = unknown / not determined
# 1 = randomly oriented ice
# 2 = water
# 3 = horizontally oriented ic


def get_calipso_phase_cloud(match, phase='water'):
    """Get CALIOP pixels determined to be water clouds."""
    # bits 6-7, start at 1 counting
    cflag = match.calipso.feature_classification_flags[::, 0]
    cal_vert_feature = np.zeros(cflag.shape) - 9.0
    feature_array = (2 * np.bitwise_and(np.right_shift(cflag, 6), 1) +
                     1 * np.bitwise_and(np.right_shift(cflag, 5), 1))
    cal_vert_feature = np.where(
        np.not_equal(cflag, 1), feature_array, cal_vert_feature)
    if phase == 'water':
        selected = cal_vert_feature == 2
    else:
        selected = np.logical_or(cal_vert_feature == 1,
                                 cal_vert_feature == 3)
    return selected


def get_semi_opaque_info_pps2014(ctth_status):
    """Get opaque/semitransparent flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (ctth_status >> 7 & 1)
    semi_flag = temp_val == 1
    opaque_flag = temp_val == 0
    return semi_flag, opaque_flag


def get_semi_opaque_info_pps2012(ctth_opaque):
    """Get opaque/semitransparent flag from PPS-v2012."""
    logger.info("Assuming  pps v2012")
    semi_flag = ctth_opaque == 0
    opaque_flag = ctth_opaque == 1
    return semi_flag, opaque_flag


def get_sunglint_info_pps2014(cloudtype_conditions):
    """Get sunglint flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (cloudtype_conditions >> 3 & 1)
    sunglint_flag = temp_val == 1
    return sunglint_flag


def get_high_terrain_info_pps2014(cloudtype_conditions):
    """Get high terrain flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (cloudtype_conditions >> 6 & 1)
    mountin_flag = temp_val == 1
    logger.debug("Number of mountain {:d}".format(len(cloudtype_conditions[mountin_flag])))
    return mountin_flag


def get_mountin_info_pps2014(cloudtype_conditions):
    """Get mountain flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (cloudtype_conditions >> 7 & 1)
    mountin_flag = temp_val == 1
    logger.debug("Number of mountain {:d}".format(len(cloudtype_conditions[mountin_flag])))
    return mountin_flag


def get_inversion_info_pps2014(cloudtype_status):
    """Get inversion flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    flag_temp = (cloudtype_status >> 0 & 1) + 0
    inversion_flag = flag_temp == 1
    return inversion_flag


def get_land_coast_sea_info_pps2014(cloudtype_conditions):
    """Get land/sea/coast flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    sealand_val = (cloudtype_conditions >> 4 & 1) + (cloudtype_conditions >> 5 & 1) * 2
    no_qflag = sealand_val == 0
    land_flag = sealand_val == 1
    sea_flag = sealand_val == 2
    coast_flag = sealand_val == 3
    land_or_sea = np.logical_or(land_flag, sea_flag)
    true_coast = np.logical_and(coast_flag, np.equal(land_or_sea, False))
    logger.debug("Number coast {:d}".format(len(cloudtype_conditions[coast_flag])))
    logger.debug("Number true coast {:d}".format(len(cloudtype_conditions[true_coast])))
    all_lsc_flag = np.bool_(np.ones(cloudtype_conditions.shape))
    return (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag)


def get_land_coast_sea_info_pps2012(cloudtype_qflag):
    """Get land/sea/coast flag from PPS-v2012."""
    logger.info("Assuming cloudtype flags structure from pps v2012")
    land_sea_val = (cloudtype_qflag >> 0 & 1)
    coast_val = (cloudtype_qflag >> 1 & 1)
    no_qflag = None
    land_flag = np.logical_and(land_sea_val == 1, coast_val == 0)
    sea_flag = np.logical_and(land_sea_val == 0, coast_val == 0)
    coast_flag = coast_val == 1
    land_or_sea = np.logical_or(land_flag, sea_flag)
    true_coast = np.logical_and(coast_flag, np.equal(land_or_sea, False))
    logger.info("Number coast {:d}".format(len(cloudtype_qflag[coast_flag])))
    logger.info("Number true coast {:d}".format(len(cloudtype_qflag[true_coast])))
    all_lsc_flag = np.bool_(np.ones(cloudtype_qflag.shape))
    return (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag)


def get_ice_info_pps2014(cloudtype_status):
    """Get ice flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    ice_flag_temp = (cloudtype_status >> 3 & 1) + 0
    ice_flag = ice_flag_temp == 1
    return ice_flag


def get_ice_info_pps2012(cloudtype_qflag):
    """Get ice flag from PPS-v2012."""
    logger.info("Assuming cloudtype flags structure from pps v2012")
    ice_flag_temp = (cloudtype_qflag >> 15 & 1) + 0
    ice_flag = ice_flag_temp == 1
    return ice_flag


def get_day_night_twilight_info_pps2014(cloudtype_conditions):
    """Get day/night/twilight flag from PPS-v2014 or later."""
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    daynight_val = (cloudtype_conditions >> 1 & 1) + (cloudtype_conditions >> 2 & 1) * 2
    no_qflag = daynight_val == 0
    night_flag = daynight_val == 1
    day_flag = daynight_val == 2
    twilight_flag = daynight_val == 3
    logger.debug("Number of day {:d}".format(len(cloudtype_conditions[day_flag])))
    logger.debug("Number of night {:d}".format(len(cloudtype_conditions[night_flag])))
    logger.debug("Number of twilight {:d}".format(len(cloudtype_conditions[twilight_flag])))
    all_dnt_flag = np.bool_(np.ones(cloudtype_conditions.shape))
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)


def get_day_night_twilight_info_pps2012(cloudtype_qflag):
    """Get day/night/twilight flag from PPS-v2012."""
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    night_flag = (((cloudtype_qflag >> 2) & 1) == 1) & ~no_qflag
    twilight_flag = (((cloudtype_qflag >> 3) & 1) == 1) & ~no_qflag
    day_flag = (((cloudtype_qflag >> 2) & 1) == 0) & (((cloudtype_qflag >> 3) & 1) == 0) & ~no_qflag
    logger.info("Number of day {:d}".format(len(cloudtype_qflag[day_flag])))
    logger.info("Number of night {:d}".format(len(cloudtype_qflag[night_flag])))
    logger.info("Number of twilight {:d}".format(len(cloudtype_qflag[twilight_flag])))
    all_dnt_flag = np.bool_(np.ones(cloudtype_qflag.shape))
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)


def get_sunglint_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    sunglint_flag = (((cloudtype_qflag >> 4) & 1) == 1) & ~no_qflag
    return sunglint_flag


def get_mountin_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    mountin_flag = (((cloudtype_qflag >> 5) & 1) == 1) & ~no_qflag
    return mountin_flag


def get_inversion_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    inversion_flag = (((cloudtype_qflag >> 6) & 1) == 1) & ~no_qflag
    return inversion_flag

# high confidence = |CAD score| >= 70
# medium confidence = 50 <= |CAD score| < 70
# low confidence = 20 <= |CAD score| < 50
# no confidence = |CAD score| < 20


def get_calipso_cad_score(match_calipso):
    """Get quality (CAD score) from CALIPSO (not for clear pixels)."""
    # bits 4-5, start at 1 counting
    cflag = match_calipso.calipso.feature_classification_flags[::, 0]
    cal_vert_feature = np.zeros(cflag.shape) - 9.0
    feature_array = (2 * np.bitwise_and(np.right_shift(cflag, 4), 1) +
                     1 * np.bitwise_and(np.right_shift(cflag, 3), 1))
    cal_vert_feature = np.where(
        np.not_equal(cflag, 1), feature_array, cal_vert_feature)
    is_medium_or_high = cal_vert_feature >= 2
    is_no_confidence = cal_vert_feature == 0
    is_no_low = cal_vert_feature == 1
    return is_medium_or_high, is_no_confidence, is_no_low


def get_calipso_cad_score_also_nodata(match_calipso):
    # bits 4-5, start at 1 counting
    cflag = match_calipso.calipso.feature_classification_flags[::, 0]
    cal_vert_feature = np.zeros(cflag.shape) - 9.0
    feature_array = (2 * np.bitwise_and(np.right_shift(cflag, 4), 1) +
                     1 * np.bitwise_and(np.right_shift(cflag, 3), 1))
    # Include also nodata retrievals
    cal_vert_feature[:] = feature_array[:]
    is_medium_or_high = cal_vert_feature >= 2
    is_no_confidence = cal_vert_feature == 0
    is_no_low = cal_vert_feature == 1
    return is_medium_or_high, is_no_confidence, is_no_low

# 0 = not determined
# 1 = clean marine
# 2 = dust
# 3 = polluted continental
# 4 = clean continental
# 5 = polluted dust
# 6 = smoke
# 7 = other


def get_calipso_aerosol_of_type_i(match_calipso, atype=0):
    """Get CALIPSO aerosols of type i."""
    # bits 10-12, start at 1 counting
    cflag = match_calipso.calipso_aerosol.feature_classification_flags[::, 0]
    cal_vert_feature = np.zeros(cflag.shape) - 9.0
    feature_array = (4 * np.bitwise_and(np.right_shift(cflag, 11), 1) +
                     2 * np.bitwise_and(np.right_shift(cflag, 10), 1) +
                     1 * np.bitwise_and(np.right_shift(cflag, 9), 1))
    cal_vert_feature = np.where(
        np.not_equal(cflag, 1), feature_array, cal_vert_feature)
    is_requested_type = cal_vert_feature == atype
    return is_requested_type


# If feature type = cloud, bits 10-12 will specify the cloud type.
# 0 = low overcast, transparent
# 1 = low overcast, opaque
# 2 = transition stratocumulus
# 3 = low, broken cumulus
# 4 = altocumulus (transparent)
# 5 = altostratus (opaque)
# 6 = cirrus (transparent)
# 7 = deep convective (opaque)


def get_calipso_clouds_of_type_i_feature_classification_flags_one_layer(cflag, calipso_cloudtype=0):
    """Get CALIPSO clouds of type i from one layer."""
    # bits 10-12, start at 1 counting
    cal_vert_feature = np.zeros(cflag.shape) - 9.0
    feature_array = (4 * np.bitwise_and(np.right_shift(cflag, 11), 1) +
                     2 * np.bitwise_and(np.right_shift(cflag, 10), 1) +
                     1 * np.bitwise_and(np.right_shift(cflag, 9), 1))
    cal_vert_feature = np.where(
        np.not_equal(cflag, 1), feature_array, cal_vert_feature)
    is_requested_type = cal_vert_feature == calipso_cloudtype
    return is_requested_type


def get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0):
    """Get CALIPSO clouds of type i from top layer."""
    # bits 10-12, start at 1 counting
    cflag = match_calipso.calipso.feature_classification_flags[::, 0]
    return get_calipso_clouds_of_type_i_feature_classification_flags_one_layer(
        cflag,
        calipso_cloudtype=calipso_cloudtype)


def get_calipso_op(match_calipso):
    """Get CALIPSO opaque clouds."""
    # type 1, 2, 5, 7 are opaque cloudtypes
    calipso_low = np.logical_or(
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=1),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=2)),
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=5),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=7)))
    return calipso_low


def get_calipso_tp(match_calipso):
    """Get CALIPSO semi-transparent clouds."""
    calipso_low = np.logical_or(
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=3)),
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=4),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=6)))
    return calipso_low


def get_calipso_low_clouds(match_calipso):
    """Get CALIPSO low clouds."""
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_low = np.logical_or(
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=1)),
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=2),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=3)))
    return calipso_low


def get_calipso_low_clouds_tp(match_calipso):
    """Get CALIPSO low and transparent clouds."""
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_low = np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=3))
    return calipso_low


def get_calipso_low_clouds_op(match_calipso):
    """Get CALIPSO low and opaque clouds."""
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_low = np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=1),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=2))
    return calipso_low


def get_calipso_high_clouds(match_calipso):
    """Get CALIPSO high clouds."""
    # type 6, 7 are high cloudtypes
    calipso_high = np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=6),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=7))
    return calipso_high


def get_calipso_medium_clouds(match_calipso):
    """Get CALIPSO medium clouds."""
    # type 6, 7 are high cloudtypes
    calipso_high = np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=4),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=5))
    return calipso_high


def get_calipso_medium_and_high_clouds_tp(match_calipso):
    """Get CALIPSO medium transparent and high transparent clouds."""
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_transp = np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=4),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=6))
    return calipso_transp


def get_calipso_low_medium_high_classification(match_calipso):
    """Get CALIPSO clouds that are transparent/opaque and/or low/medium/high."""
    mlh_class = {}
    mlh_class['clouds_op'] = get_calipso_op(match_calipso)
    mlh_class['clouds_tp'] = get_calipso_tp(match_calipso)
    mlh_class['low_clouds'] = get_calipso_low_clouds(match_calipso)
    mlh_class['medium_clouds'] = get_calipso_medium_clouds(match_calipso)
    mlh_class['high_clouds'] = get_calipso_high_clouds(match_calipso)
    mlh_class['low_clouds_op'] = get_calipso_low_clouds_op(match_calipso)
    mlh_class['low_clouds_tp'] = get_calipso_low_clouds_tp(match_calipso)
    mlh_class['medium_clouds_tp'] = get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=4)
    mlh_class['medium_clouds_op'] = get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=5)
    mlh_class['high_clouds_tp'] = get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=6)
    mlh_class['high_clouds_op'] = get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=7)
    return mlh_class


def get_cloudsat_low_medium_high_classification(match_clsat):
    """Get low, medium and high clouds from CPR (CloudSat)."""
    mlh_class = {}
    h440 = None
    h680 = None
    if (hasattr(match_clsat.imager, 'segment_nwp_h440') and
            match_clsat.imager.all_arrays['segment_nwp_h440'] is not None):
        h680 = match_clsat.imager.all_arrays['segment_nwp_h680']
        h440 = match_clsat.imager.all_arrays['segment_nwp_h440']
    if (hasattr(match_clsat.imager, 'nwp_h440') and
            match_clsat.imager.all_arrays['nwp_h440'] is not None):
        h680 = match_clsat.imager.all_arrays['nwp_h680']
        h440 = match_clsat.imager.all_arrays['nwp_h440']

    clsat_h = match_clsat.cloudsat.validation_height
    if h680 is not None and h440 is not None:
        mlh_class['low_clouds'] = np.less_equal(clsat_h, h680)
        mlh_class['medium_clouds'] = np.logical_and(np.greater(clsat_h, h680),
                                                    np.less(clsat_h, h440))
        mlh_class['high_clouds'] = np.greater_equal(clsat_h, h440)
    return mlh_class


def get_land_coast_sea_info_cci2014(lsflag):
    """Get land/sea/cost flag from CCI lsflag."""
    logger.info("Assuming cloudtype flags structure from CCI v3?")
    land_flag = lsflag
    sea_flag = 1 - lsflag
    coast_flag = None
    all_lsc_flag = np.bool_(np.ones(lsflag))
    no_qflag = np.bool_(np.zeros(lsflag))
    return (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag)


def get_day_night_twilight_info_cci2014(sunz):
    """Get day/night/twilight flag from sun zenith angle."""
    sunz = np.array(sunz)
    logger.info("Getting day/night info from sunz")
    no_qflag = np.zeros(sunz.shape)
    day_flag = np.where(np.less_equal(sunz, 80), True, False)
    night_flag = np.where(np.greater_equal(sunz, 95), True, False)
    twilight_flag = np.where(
        np.logical_and(np.greater(sunz, 80),
                       np.less(sunz, 95)),
        True, False)
    no_qflag = np.where(np.isnan(sunz), True, False)
    logger.debug("number of day {:d}".format(np.sum(day_flag)))
    logger.debug("number of night {:d}".format(np.sum(night_flag)))
    logger.debug("number of twilight {:d}".format(np.sum(twilight_flag)))
    all_dnt_flag = np.bool_(np.ones(sunz.shape))
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)


def get_maia_ct_flag(ct_flag):
    """Get maia cloud type flag."""
    # bits 4-8, start at 0 counting
    maia_ct_flag = (16 * np.bitwise_and(np.right_shift(ct_flag, 8), 1) +
                    8 * np.bitwise_and(np.right_shift(ct_flag, 7), 1) +
                    4 * np.bitwise_and(np.right_shift(ct_flag, 6), 1) +
                    2 * np.bitwise_and(np.right_shift(ct_flag, 5), 1) +
                    1 * np.bitwise_and(np.right_shift(ct_flag, 4), 1))
    return maia_ct_flag


def get_day_night_twilight_info_maia(cm_flag):
    """Get day/night/twilight flag from MAIA cm_flag."""
    # bit 13-14
    maia_cm_flag = (2 * np.bitwise_and(np.right_shift(cm_flag, 14), 1) +
                    1 * np.bitwise_and(np.right_shift(cm_flag, 13), 1))
    day_flag = np.where(maia_cm_flag == 2, True, False)  # include also sunglint in day ==3
    night_flag = np.where(maia_cm_flag == 0, True, False)
    twilight_flag = np.where(maia_cm_flag == 1, True, False)
    logger.debug("number of day {:d}".format(np.sum(day_flag)))
    logger.debug("number of night {:d}".format(np.sum(night_flag)))
    logger.debug("number of twilight {:d}".format(np.sum(twilight_flag)))
    all_dnt_flag = np.ones(maia_cm_flag.shape)
    no_qflag = 0 * all_dnt_flag  # dummy flag
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)


def get_pixels_where_test_is_passed(cma_testlist, bit_nr=0):
    """Get testflag for one CMA test."""
    # print bit_nr,
    # print
    test_was_fullfilled = (cma_testlist.astype(np.uint16) >> bit_nr & 1)
    return test_was_fullfilled
