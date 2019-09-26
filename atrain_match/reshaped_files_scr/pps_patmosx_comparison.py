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

"""Read all matched data and make tables
"""
import os
from glob import glob
import numpy as np
from scipy import ndimage
from matchobject_io import (read_files)


from utils.stat_util import (my_hist,
                       my_iqr,
                       my_rms,
                       my_mae,
                       half_sample_mode,
                       half_sample_mode,
                       my_pe250m,
                       my_pe500m,
                       my_pe1000m,
                       my_pe2000m,
                       my_pe2500m,
                       my_pe5000m)
from my_dir import ADIR
def crop_object(match_obj, use_in=None):
    y =match_obj.imager.all_arrays['ctth_height']
    if 'ctth_height_corr' in match_obj.imager.all_arrays.keys() and   match_obj.imager.all_arrays['ctth_height_corr'] is not None:
        y =match_obj.imager.all_arrays['ctth_height_corr']
    x = match_obj.calipso.all_arrays['validation_height']
    pps_profile_id = match_obj.calipso.sec_1970#profile_id[:, 0]
    use = np.logical_and(y >= 0, x >= 0)
    if use_in is not None:
        use = np.logical_and(use, use_in)
    for arnameca, valueca in match_obj.calipso.all_arrays.items():
        if match_obj.calipso.all_arrays[arnameca] is None:
            pass
        elif arnameca in ['ctth_height', 'validation_height', 'sec_1970', 'imager_ctth_m_above_seasurface', 'ctth_height_corr']:
            match_obj.calipso.all_arrays[arnameca] = match_obj.calipso.all_arrays[arnameca][use]
        else:
            match_obj.calipso.all_arrays[arnameca] = None
    for arnameca, valueca in match_obj.imager.all_arrays.items():
        if match_obj.imager.all_arrays[arnameca] is None:
            pass
        elif arnameca in ['ctth_height', 'validation_height', 'sec_1970', 'imager_ctth_m_above_seasurface', 'ctth_height_corr']:
            match_obj.imager.all_arrays[arnameca] = match_obj.imager.all_arrays[arnameca][use]
        else:
            match_obj.imager.all_arrays[arnameca] = None
    return match_obj


def remove_missing(match_objPPS, match_objPATMOSX, common_index):

    patmosx_profile_id = match_objPATMOSX.calipso.sec_1970#profile_id[:, 0]
    pps_profile_id = match_objPPS.calipso.sec_1970#profile_id[:, 0]

    use_patmosx_same_profile = np.array([p_id in common_index for p_id in patmosx_profile_id])
    use_pps_same_profile = np.array([p_id in  common_index for p_id in pps_profile_id])
    use_patmosx = use_patmosx_same_profile
    use_pps = use_pps_same_profile
    """
    pps_profile_id[~use_pps]= -111
    # remove doubles
    unique, index = np.unique(pps_profile_id, return_index=True)
    pps_profile_id_new = -111 + 0*pps_profile_id.copy()
    pps_profile_id_new[index] = unique
    pps_profile_id = pps_profile_id_new

    unique, index = np.unique(patmosx_profile_id, return_index=True)
    patmosx_profile_id_new = -333 + 0*patmosx_profile_id.copy()
    patmosx_profile_id_new[index] = unique
    patmosx_profile_id = patmosx_profile_id_new

    use_patmosx_same_profile = np.array([p_id in pps_profile_id[use_pps] for p_id in patmosx_profile_id])
    use_pps_same_profile = np.array([p_id in patmosx_profile_id[use_patmosx] for p_id in pps_profile_id])
    use_patmosx = np.logical_and(use_patmosx, use_patmosx_same_profile)
    use_pps = np.logical_and(use_pps, use_pps_same_profile)
    """
    return use_pps, use_patmosx

def print_stats(match_objPPS, match_objPATMOSX, use_pps, use_patmosx):

    # print(sorted(patmosx_profile_id[use_patmosx])[-10:])
    # print(sorted(pps_profile_id[use_pps])[-10:])
    x = match_objPPS.calipso.all_arrays['validation_height']
    x_patmosx = match_objPATMOSX.calipso.all_arrays['validation_height']
    y_pps = match_objPPS.imager.all_arrays['imager_ctth_m_above_seasurface']
    if 'ctth_height_corr' in match_objPPS.imager.all_arrays.keys() and   match_objPPS.imager.all_arrays['ctth_height_corr'] is not None:
        y_pps =match_objPPS.imager.all_arrays['ctth_height_corr']
    y_patmosx =match_objPATMOSX.imager.all_arrays['ctth_height']

    print(np.sum(use_pps), np.sum(use_patmosx))

    use_patmosx = np.logical_and(use_patmosx, x_patmosx > 0)
    use_pps = np.logical_and(use_pps, x > 0)

    bias_pps = y_pps[use_pps] - x[use_pps]
    bias_patmosx = y_patmosx[use_patmosx] - x_patmosx[use_patmosx]
    abias_pps = np.abs(bias_pps)
    abias_patmosx = np.abs(bias_patmosx)

    print("MAE: {:3.1f}, ({:3.1f}), PE05: {:3.1f}, ({:3.1f}), Median: {:3.1f}, ({:3.1f}) IQR: {:3.1f}, ({:3.1f}) N: {:d}, ({:d}) BIAS: {:3.1f}, ({:3.1f}) STD:{:3.1f}, ({:3.1f})".format(
        np.mean(abias_pps),
        np.mean(abias_patmosx),
        my_pe500m(abias_pps),
        my_pe500m(abias_patmosx),
        np.median(bias_pps),
        np.median(bias_patmosx),
        my_iqr(bias_pps),
        my_iqr(bias_patmosx),
        len(bias_pps),
        len(bias_patmosx),
        np.mean(bias_pps),
        np.mean(bias_patmosx),
        np.std(bias_pps),
        np.std(bias_patmosx),                                       ))

if __name__ == "__main__":

    PATMOSX_ROOT_DIR = (ADIR + "/VALIDATION_PATMOSX/Reshaped_Files/noaa18/5km/2009/*/*h5")

    PPS_ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/2009/5km_noaa18_2009*cali*h5")
    PPS14_ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/2009/5km_noaa18_2009*cali*h5")
    CCI_ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files_cci_noaa18_2009/V2/*2009*h5")


    patmosx_files = glob(PATMOSX_ROOT_DIR)

    match_objPATMOSX = read_files(glob(PATMOSX_ROOT_DIR))
    match_objPATMOSX = crop_object(match_objPATMOSX, use_in=None)
    match_objPPS = read_files(glob(PPS_ROOT_DIR))
    match_objPPS = crop_object(match_objPPS, use_in=None)
    match_objPPS14 = read_files(glob(PPS14_ROOT_DIR))
    match_objPPS14 = crop_object(match_objPPS14, use_in=None)
    match_objCCI = read_files(glob(CCI_ROOT_DIR))
    match_objCCI = crop_object(match_objCCI, use_in=None)

    patmosx_profile_id = match_objPATMOSX.calipso.sec_1970#profile_id[:, 0]
    pps_profile_id = match_objPPS.calipso.sec_1970#profile_id[:, 0]
    pps14_profile_id = match_objPPS14.calipso.sec_1970#profile_id[:, 0]
    cci_profile_id = match_objCCI.calipso.sec_1970#profile_id[:, 0]
    common_index1 = np.intersect1d(patmosx_profile_id, pps_profile_id)
    common_index2 = np.intersect1d(cci_profile_id, pps14_profile_id)
    common_index = np.intersect1d(common_index1, common_index2)

    use_pps_v14, use_patmosx = remove_missing(match_objPPS14, match_objPATMOSX, common_index)
    use_pps, use_cci = remove_missing(match_objPPS, match_objCCI, common_index)


    print("PPS-v2014")
    print_stats(match_objPPS14, match_objPATMOSX, use_pps_v14, use_patmosx)
    print("PPS-vCCI")
    print_stats(match_objCCI, match_objPATMOSX, use_cci, use_patmosx)
    print("PPS-v2018")
    print_stats(match_objPPS, match_objPATMOSX, use_pps, use_patmosx)
    print("Dummy")
    match_objPPS.imager.all_arrays['imager_ctth_m_above_seasurface'][:] = np.mean(match_objPPS.calipso.all_arrays['validation_height'][use_pps])
    print(np.mean(match_objPPS.calipso.all_arrays['validation_height'][use_pps]))
    print_stats(match_objPPS, match_objPATMOSX, use_pps, use_patmosx)
