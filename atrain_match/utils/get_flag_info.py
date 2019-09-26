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



def get_semi_opaque_info_pps2014(ctth_status):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (ctth_status >> 7 & 1)
    semi_flag = temp_val == 1
    opaque_flag = temp_val == 0
    return  semi_flag, opaque_flag

def get_semi_opaque_info_pps2012(ctth_opaque):
    logger.info("Assuming  pps v2012")
    semi_flag = ctth_opaque == 0
    opaque_flag = ctth_opaque == 1
    return  semi_flag, opaque_flag

def get_sunglint_info_pps2014(cloudtype_conditions):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (cloudtype_conditions >> 3 & 1)
    sunglint_flag = temp_val == 1
    return  sunglint_flag

def get_high_terrain_info_pps2014(cloudtype_conditions):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (cloudtype_conditions >> 6 & 1)
    mountin_flag = temp_val == 1
    logger.debug("Number of mountain %d"%(len(cloudtype_conditions[mountin_flag==True])))
    return  mountin_flag

def get_mountin_info_pps2014(cloudtype_conditions):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    temp_val = (cloudtype_conditions >> 7 & 1)
    mountin_flag = temp_val == 1
    logger.debug("Number of mountain %d"%(len(cloudtype_conditions[mountin_flag==True])))
    return  mountin_flag

def get_inversion_info_pps2014(cloudtype_status):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    flag_temp = (cloudtype_status >> 0 & 1) + 0
    inversion_flag = flag_temp == 1
    return inversion_flag

def get_land_coast_sea_info_pps2014(cloudtype_conditions):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    sealand_val = (cloudtype_conditions >> 4 & 1) + (cloudtype_conditions >> 5 & 1) * 2
    no_qflag = sealand_val == 0
    land_flag =  sealand_val == 1
    sea_flag =  sealand_val == 2
    coast_flag =   sealand_val == 3
    land_or_sea=np.logical_or(land_flag, sea_flag)
    true_coast=np.logical_and(coast_flag, np.equal(land_or_sea, False))
    logger.debug("Number coast %d"%(len(cloudtype_conditions[coast_flag==True])))
    logger.debug("Number true coast %d"%(len(cloudtype_conditions[true_coast==True])))
    all_lsc_flag =  np.bool_(np.ones(cloudtype_conditions.shape))
    return (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag)

def get_land_coast_sea_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    land_sea_val = (cloudtype_qflag >> 0 & 1)
    coast_val = (cloudtype_qflag >> 1 & 1)
    no_qflag = None
    land_flag =  np.logical_and(land_sea_val == 1, coast_val == 0)
    sea_flag =  np.logical_and(land_sea_val == 0, coast_val == 0)
    coast_flag = coast_val == 1
    land_or_sea=np.logical_or(land_flag, sea_flag)
    true_coast=np.logical_and(coast_flag, np.equal(land_or_sea, False))
    logger.info("Number coast %d"%( len(cloudtype_qflag[coast_flag==True])))
    logger.info("Number true coast %d"%( len(cloudtype_qflag[true_coast==True])))
    all_lsc_flag =  np.bool_(np.ones(cloudtype_qflag.shape))
    return (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag)




def get_ice_info_pps2014(cloudtype_status):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    ice_flag_temp = (cloudtype_status >> 3 & 1) + 0
    ice_flag = ice_flag_temp == 1
    return ice_flag

def get_ice_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    ice_flag_temp = (cloudtype_qflag >> 15 & 1) + 0
    ice_flag = ice_flag_temp == 1
    return ice_flag

def get_day_night_twilight_info_pps2014(cloudtype_conditions):
    logger.debug("Assuming cloudtype flags structure from pps v2014")
    daynight_val = (cloudtype_conditions >> 1 & 1) + (cloudtype_conditions >> 2 & 1) * 2
    no_qflag = daynight_val == 0
    night_flag =  daynight_val == 1
    day_flag =   daynight_val == 2
    twilight_flag =  daynight_val == 3
    logger.debug("Number of day %d"%(len(cloudtype_conditions[day_flag==True])))
    logger.debug("Number of night %d"%(len(cloudtype_conditions[night_flag==True])))
    logger.debug("Number of twilight %d"%(len(cloudtype_conditions[twilight_flag==True])))
    all_dnt_flag =  np.bool_(np.ones(cloudtype_conditions.shape))
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)

def get_day_night_twilight_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    night_flag = (((cloudtype_qflag >> 2) & 1) == 1) & ~no_qflag
    twilight_flag = (((cloudtype_qflag >> 3) & 1) == 1) & ~no_qflag
    day_flag =  (((cloudtype_qflag >> 2) & 1) == 0) & (((cloudtype_qflag >> 3) & 1) == 0) & ~no_qflag
    logger.info("Number of day %d"%(len(cloudtype_qflag[day_flag==True])))
    logger.info("Number of night %d"%(len(cloudtype_qflag[night_flag==True])))
    logger.info("Number of twilight %d"%(len(cloudtype_qflag[twilight_flag==True])))
    all_dnt_flag =  np.bool_(np.ones(cloudtype_qflag.shape))
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)

def get_sunglint_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    sunglint_flag = (((cloudtype_qflag >> 4) & 1) == 1) & ~no_qflag
    return  sunglint_flag

def get_mountin_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    mountin_flag = (((cloudtype_qflag >> 5) & 1) == 1) & ~no_qflag
    return  mountin_flag

def get_inversion_info_pps2012(cloudtype_qflag):
    logger.info("Assuming cloudtype flags structure from pps v2012")
    no_qflag = cloudtype_qflag == 0
    inversion_flag = (((cloudtype_qflag >> 6) & 1) == 1) & ~no_qflag
    return inversion_flag

# 0 = not determined
# 1 = clean marine
# 2 = dust
# 3 = polluted continental
# 4 = clean continental
# 5 = polluted dust
# 6 = smoke
# 7 = other


def get_calipso_aerosol_of_type_i(match_calipso, atype=0):
    # bits 10-12, start at 1 counting
    cflag = match_calipso.calipso_aerosol.feature_classification_flags[::, 0]
    cal_vert_feature = np.zeros(cflag.shape) - 9.0
    feature_array = (4 * np.bitwise_and(np.right_shift(cflag, 11), 1) +
                     2 * np.bitwise_and(np.right_shift(cflag, 10), 1) +
                     1 * np.bitwise_and(np.right_shift(cflag, 9), 1))
    cal_vert_feature = np.where(
        np.not_equal(cflag, 1), feature_array, cal_vert_feature)

    is_requested_type =  cal_vert_feature == atype
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
    # bits 10-12, start at 1 counting
    cal_vert_feature = np.zeros(cflag.shape) - 9.0
    feature_array = (4 * np.bitwise_and(np.right_shift(cflag, 11), 1) +
                     2 * np.bitwise_and(np.right_shift(cflag, 10), 1) +
                     1 * np.bitwise_and(np.right_shift(cflag, 9), 1))
    cal_vert_feature = np.where(
        np.not_equal(cflag, 1), feature_array, cal_vert_feature)

    is_requested_type =  cal_vert_feature == calipso_cloudtype
    return is_requested_type

def get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0):
    # bits 10-12, start at 1 counting
    cflag = match_calipso.calipso.feature_classification_flags[::, 0]
    return get_calipso_clouds_of_type_i_feature_classification_flags_one_layer(cflag, calipso_cloudtype=calipso_cloudtype)

def get_calipso_op(match_calipso):
    # type 1, 2, 5, 7 are low cloudtypes
    calipso_low =  np.logical_or(
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=1),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=2)),
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=5),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=7)))
    return calipso_low
    # type 0, 1, 2, 3 are low cloudtypes
def get_calipso_tp(match_calipso):
    calipso_low =  np.logical_or(
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=3)),
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=4),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=6)))
    return calipso_low




def get_calipso_low_clouds(match_calipso):
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_low =  np.logical_or(
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=1)),
        np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=2),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=3)))
    return calipso_low

def get_calipso_low_clouds_tp(match_calipso):
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_low =  np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=0),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=3))
    return calipso_low
def get_calipso_low_clouds_op(match_calipso):
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_low =  np.logical_or(
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=1),
            get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=2))
    return calipso_low

def get_calipso_high_clouds(match_calipso):
    # type 6, 7 are high cloudtypes
    calipso_high =   np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=6),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=7))
    return calipso_high
def get_calipso_medium_clouds(match_calipso):
    # type 6, 7 are high cloudtypes
    calipso_high =   np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=4),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=5))
    return calipso_high

def get_calipso_medium_and_high_clouds_tp(match_calipso):
    # type 0, 1, 2, 3 are low cloudtypes
    calipso_transp =   np.logical_or(
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=4),
        get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=6))
    return calipso_transp





def get_calipso_low_medium_high_classification(match_calipso):
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
    mlh_class = {}
    if not hasattr( match_clsat.imager, 'segment_nwp_h440') :
        return None
    if match_clsat.imager.all_arrays['segment_nwp_h440'] is None:
        return  None
    clsat_h = match_clsat.cloudsat.validation_height
    mlh_class['low_clouds'] = np.less_equal(clsat_h, match_clsat.imager.all_arrays['segment_nwp_h680'])
    mlh_class['medium_clouds'] = np.logical_and(np.greater(clsat_h, match_clsat.imager.all_arrays['segment_nwp_h680']),
                                                np.less(clsat_h, match_clsat.imager.all_arrays['segment_nwp_h440']))
    mlh_class['high_clouds'] = np.greater_equal(clsat_h, match_clsat.imager.all_arrays['segment_nwp_h440'])
    return mlh_class

# cci FLAGS
def get_land_coast_sea_info_cci2014(lsflag):
    logger.info("Assuming cloudtype flags structure from CCI v2014?")
    land_flag =  lsflag
    sea_flag =  1 - lsflag
    coast_flag = None
    all_lsc_flag =  np.bool_(np.ones(lsflag))
    return (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag)
def get_day_night_twilight_info_cci2014(sunz):
    sunz = np.array(sunz)
    logger.info("Getting day/night info from sunz")
    no_qflag = np.zeros(sunz.shape)
    day_flag = np.where(np.less_equal(sunz, 80), 1, 0)
    night_flag =  np.where(np.greater_equal(sunz, 95), 1, 0)
    twilight_flag =  np.where(
        np.logical_and(np.greater(sunz, 80),
                       np.less(sunz, 95)),
        1, 0)
    no_qflag = np.where(np.isnan(sunz), 1, 0)
    logger.debug("number of day {:d}".format(len(sunz[day_flag==True])))
    logger.debug("number of night {:d}".format(len(sunz[night_flag==True])))
    logger.debug("number of twilight {:d}".format(len(sunz[twilight_flag==True])))
    all_dnt_flag =  np.bool_(np.ones(sunz.shape))
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)

def get_maia_ct_flag(ct_flag):
    # bits 4-8, start at 0 counting
    maia_ct_flag = (16 * np.bitwise_and(np.right_shift(ct_flag, 8), 1) +
                    8 * np.bitwise_and(np.right_shift(ct_flag, 7), 1) +
                    4 * np.bitwise_and(np.right_shift(ct_flag, 6), 1) +
                    2 * np.bitwise_and(np.right_shift(ct_flag, 5), 1) +
                    1 * np.bitwise_and(np.right_shift(ct_flag, 4), 1))
    return  maia_ct_flag
def get_day_night_twilight_info_maia(cm_flag):
    # bit 13-14
    maia_cm_flag = (2 * np.bitwise_and(np.right_shift(cm_flag, 14), 1) +
                    1 * np.bitwise_and(np.right_shift(cm_flag, 13), 1))
    day_flag = np.where(maia_cm_flag == 2, 1, 0)  # include also sunglint in day ==3
    night_flag =  np.where(maia_cm_flag == 0, 1, 0)
    twilight_flag =  np.where(maia_cm_flag == 1, 1, 0)
    logger.debug("number of day {:d}".format(len(maia_cm_flag[day_flag==True])))
    logger.debug("number of night {:d}".format(len(maia_cm_flag[night_flag==True])))
    logger.debug("number of twilight {:d}".format(len(maia_cm_flag[twilight_flag==True])))
    all_dnt_flag =  np.ones(maia_cm_flag.shape)
    no_qflag = 0 * all_dnt_flag  # dummy flag
    return (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag)

def get_pixels_where_test_is_passed(cma_testlist, bit_nr=0):
    # print bit_nr,
    # print
    test_was_fullfilled = (cma_testlist.astype(np.uint16) >> bit_nr & 1)
    return test_was_fullfilled
