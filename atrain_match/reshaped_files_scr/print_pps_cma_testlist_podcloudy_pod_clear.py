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
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)
import warnings
warnings.filterwarnings("ignore")
from utils.get_flag_info import get_pixels_where_test_is_passed

v2014 = False
from my_dir import ADIR


def print_common_stats(match_calipso, use, name_dict, mints, maxts, surface_type, illumination, basename_outfile="basename_outfile"):
    if v2014:
        sunz = 100*match_calipso.imager.all_arrays['sunz']
    else:
        sunz = match_calipso.imager.all_arrays['sunz']
    print np.min(sunz), np.max(sunz)

    outfile=(ADIR + "/DATA_MISC/cma_log_files_can_be_removed_after_some_time"
             "/log_cma_thin_%s_%s_%s_temp%d_%d_.txt"%(basename_outfile, surface_type, illumination, mints, maxts))
    print outfile
    outfile_h = open(outfile, 'w')
    all_out_text = ""
    nlay =np.where(match_calipso.calipso.all_arrays['number_layers_found']>0, 1, 0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
    """
    isCalipsoCloudy = match_calipso.calipso.all_arrays['cloud_fraction']>=0.5

    isCalipsoCloudy = np.logical_and(
        nlay >= 0,
        match_calipso.calipso.all_arrays['cloud_fraction']>=0.5)
    isCalipsoCloudy = np.logical_and(
        isCalipsoCloudy,
        match_calipso.calipso.all_arrays['total_optical_depth_5km']>-0.1)
    isCalipsoClear = match_calipso.calipso.all_arrays['cloud_fraction']<0.5
    isCalipsoClear = nlay == 0
    isCalipsoClear = np.logical_and(isCalipsoClear, meancl < 0.01)
    isCalipsoClear = np.logical_and(
        isCalipsoClear,
        match_calipso.calipso.all_arrays['total_optical_depth_5km']<0)
    """
    isCalipsoCloudy = match_calipso.calipso.all_arrays['cloud_fraction']>=0.5
    isCalipsoCloudy = np.logical_and(
        isCalipsoCloudy,
        np.logical_or(match_calipso.calipso.all_arrays['total_optical_depth_5km']>0.1,
                       match_calipso.calipso.all_arrays['total_optical_depth_5km']<0))
    isCalipsoClear = match_calipso.calipso.all_arrays['cloud_fraction']<0.0001
    isCalipsoSnowIce = np.logical_and(
        isCalipsoClear,
        match_calipso.calipso.all_arrays['nsidc_surface_type']>50)
    isCalipsoSnowIceWithCloudAbove = np.logical_and(
        isCalipsoCloudy,
        match_calipso.calipso.all_arrays['nsidc_surface_type']>5)
    # print np.sum(isCalipsoSnowIce)
    isCalipsoNotSnowIce = np.logical_and(
        isCalipsoClear,
        match_calipso.calipso.all_arrays['nsidc_surface_type'] <= 0)
    isClearPPS = np.logical_or(np.equal(match_calipso.imager.cloudmask, 3),
                                   np.equal(match_calipso.imager.cloudmask, 0))
    isCloudyPPS = np.logical_or(np.equal(match_calipso.imager.cloudmask, 1),
                                   np.equal(match_calipso.imager.cloudmask, 2))

    gotLight = sunz < 95
    nodata = np.sum(match_calipso.imager.all_arrays['cloudtype'][use]>200)
    use = np.logical_and(use, np.logical_or(isCloudyPPS, isClearPPS))
    use = np.logical_and(use, np.logical_or(isCalipsoCloudy, isCalipsoClear))
    use = np.logical_and(use, match_calipso.imager.all_arrays['surftemp']<=maxts)
    use = np.logical_and(use, match_calipso.imager.all_arrays['surftemp']>mints)
    from utils.get_flag_info import  get_land_coast_sea_info_pps2014
    (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag
    ) = get_land_coast_sea_info_pps2014(match_calipso.imager.all_arrays['cloudtype_conditions'])
    if surface_type in ["land"]:
        # use = np.logical_and(use, np.not_equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, land_flag)
    elif surface_type in ["emiss"]:
        # use = np.logical_and(use, np.not_equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, np.logical_or(land_flag, coast_flag))
        use = np.logical_and(use, match_calipso.imager.all_arrays['emis1']<1.0)
    elif surface_type in ["emiss_coast"]:
        # use = np.logical_and(use, np.not_equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, coast_flag)
        use = np.logical_and(use, match_calipso.imager.all_arrays['emis1']<1.0)
    elif surface_type in ["emiss_land"]:
        # use = np.logical_and(use, np.not_equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, land_flag)
        use = np.logical_and(use, match_calipso.imager.all_arrays['emis1']<1.0)
    elif surface_type in ["noemiss_coast"]:
        # use = np.logical_and(use, np.not_equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, coast_flag)
        use = np.logical_and(use, match_calipso.imager.all_arrays['emis1']==1.0)
    elif surface_type in ["noemiss_land"]:
        # use = np.logical_and(use, np.not_equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, land_flag)
        use = np.logical_and(use, match_calipso.imager.all_arrays['emis1']==1.0)
    elif surface_type in ["noemiss"]:
        # use = np.logical_and(use, np.not_equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, np.logical_or(land_flag, coast_flag))
        use = np.logical_and(use, match_calipso.imager.all_arrays['emis1']==1.0)
        import matplotlib.pyplot as plt
        plt.plot(match_calipso.imager.all_arrays['longitude'][use], match_calipso.imager.all_arrays['latitude'][use], 'b.')
        plt.savefig('emiss_missing.png')
        # plt.show()
    elif surface_type in ["sea"]:
        # use = np.logical_and(use, np.equal(match_calipso.calipso.igbp_surface_type, 17))
        use = np.logical_and(use, sea_flag)
    elif surface_type in ["coast"]:
        use = np.logical_and(use, coast_flag)
    elif surface_type in ["all"]:
        pass
    else:
        sys.exit()
    if illumination in ["day"]:
        use = np.logical_and(use, sunz <= 80)
    elif illumination in ["night"]:
        use = np.logical_and(use, sunz >= 95)
    elif illumination in ["twilight"]:
        use = np.logical_and(use, np.logical_and(
            sunz < 95,
            sunz > 80))
    elif illumination in ["dnt"]:
        pass
    else:
        sys.exit()

    N_detected_clouds = np.sum(np.logical_and(isCalipsoCloudy[use],
                                              isCloudyPPS[use]))
    N_detected_clear = np.sum(np.logical_and(isCalipsoClear[use],
                                             isClearPPS[use]))
    N_false_clouds = np.sum(np.logical_and(isCalipsoClear[use],
                                           isCloudyPPS[use]))
    N_undetected_clouds = np.sum(np.logical_and(isCalipsoCloudy[use],
                                                isClearPPS[use]))
    N_clouds = N_detected_clouds + N_undetected_clouds
    N_clear = N_detected_clear + N_false_clouds

    Kuipers_devider = 1.0*(N_clouds)*(N_clear)
    if Kuipers_devider == 0:
        Kuipers_devider = 1.0
    Kuipers = (N_detected_clouds*N_detected_clear -
               N_false_clouds*N_undetected_clouds)/Kuipers_devider

    NALL = N_clouds + N_clear
    PODcloudy = N_detected_clouds*1.0/N_clouds
    FARcloudy = N_false_clouds *1.0/np.sum(isCloudyPPS[use])
    PODclear = N_detected_clear*1.0/N_clear
    Hitrate = (N_detected_clear + N_detected_clouds)*1.0/NALL

    Num = np.sum(use)
    cfc = np.sum(isCalipsoCloudy[use])*1.0/Num
    part_nodata = nodata*1.0/(nodata + Num)
    all_out_text += "N: %d POD-cloudy: %3.1f FAR-cloudy: %3.1f POD-clear %3.1f %3.3f Hitrate:  %3.1f Kuipers:  %3.3f\n"%(
        Num, 100.0*PODcloudy, 100.0*FARcloudy , 100.0*PODclear, cfc, 100*Hitrate, Kuipers )
    print all_out_text
    all_pix = use.copy()
    if v2014:
          outfile_h.write(all_out_text+"\n")
          Num = np.sum(use)
          part_nodata = nodata*1.0/(nodata + Num)
          return
    for var in ['cma_testlist0',
                'cma_testlist1',
                'cma_testlist2',
                'cma_testlist3',
                'cma_testlist4',
                'cma_testlist5',
            ]:
        # print var
        for bit_nr in range(0, 16):
            if bit_nr not in name_dict[var].keys():
                print "not using", var, bit_nr
                continue
            test_is_on = get_pixels_where_test_is_passed(
                match_calipso.imager.all_arrays[var], bit_nr=bit_nr)
            all_pix = np.logical_and(all_pix, np.equal(test_is_on, False))
            use_this_test = np.logical_and(use, test_is_on)
            # match_calipso.imager.all_arrays[var].
            # print np.sum(test_is_on), np.sum(all_pix), np.sum(use)
            # print np.sum(match_calipso.imager.all_arrays[var]==4)
            # print np.sum(np.logical_and(~test_is_on,
            #                           match_calipso.imager.all_arrays[var]==4))
            PODcloudy = (np.sum(np.logical_and(isCalipsoCloudy[use_this_test],
                                               isCloudyPPS[use_this_test]))*1.0
                         /np.sum(isCalipsoCloudy[use]))
            FARcloudy = (np.sum(np.logical_and(isCalipsoClear[use_this_test],
                                           isCloudyPPS[use_this_test]))*1.0
                         /np.sum(isCloudyPPS[use_this_test]))
            PODclear = (np.sum(np.logical_and(isCalipsoClear[use_this_test],
                                              isClearPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoClear[use]))
            PODsnow = (np.sum(np.logical_and(isCalipsoSnowIce[use_this_test],
                                            isClearPPS[use_this_test]))*1.0
                       /np.sum(np.logical_and(isCalipsoSnowIce[use], gotLight[use])))
            FARsnow_not_clouds = (np.sum(np.logical_and(isCalipsoNotSnowIce[use_this_test],
                                                        isClearPPS[use_this_test]))*1.0
                                  /np.sum(isClearPPS[use_this_test]))
            FARclear = (np.sum(np.logical_and(isCalipsoCloudy[use_this_test],
                                              isClearPPS[use_this_test]))*1.0
                        /np.sum(isClearPPS[use_this_test]))
            MISSclear = (np.sum(np.logical_and(isCalipsoClear[use_this_test],
                                              isCloudyPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoClear[use]))
            MISScloudy = (np.sum(np.logical_and(isCalipsoCloudy[use_this_test],
                                              isClearPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoCloudy[use]))
            MISScloudyOverSnow = (np.sum(np.logical_and(
                isCalipsoSnowIceWithCloudAbove[use_this_test],
                isClearPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoCloudy[use]))
            MISSsnow = (np.sum(np.logical_and(
                np.logical_and(
                    isCalipsoSnowIce[use_this_test],
                    gotLight[use_this_test]),
                isCloudyPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoClear[use]))
            LOSTsnow = (np.sum(np.logical_and(
                np.logical_and(
                    isCalipsoSnowIce[use_this_test],
                    gotLight[use_this_test]),
                isCloudyPPS[use_this_test]))*1.0
                        /np.sum(np.logical_and(isCalipsoSnowIce[use], gotLight[use])))

            PRINT_FOR_OUTPUT_ON_SCREEN = False
            test_name = name_dict[var][bit_nr]
            if PRINT_FOR_OUTPUT_ON_SCREEN:
                all_out_text += "%s test_bit: %s"%(var, str(bit_nr).rjust(2, ' ') )

            clear_test=False
            if (var == 'cma_testlist5' and bit_nr<9) or (var == 'cma_testlist4' and bit_nr>1) or (var == 'cma_testlist5' and bit_nr==12)  or (var == 'cma_testlist5' and bit_nr==13):
                clear_test=True
            num_passed = np.sum(use_this_test)
            print num_passed, name_dict[var][bit_nr]
            if clear_test and 'now' in name_dict[var][bit_nr] and num_passed>0:
                all_out_text += "N: %s POD-clear: %s FAR-clear: %s POD-snow: %s FAR_not_clouds %s LOSTsnow %s MISScloudy: %s (%s) %s \n"%(
                    str(num_passed).rjust(8, ' '),
                    ("%3.2f"%(PODclear*100)).rjust(5, ' '),
                    ("%3.2f"%(FARclear*100)).rjust(5, ' '),
                    ("%3.2f"%(PODsnow*100)).rjust(5, ' '),
                    ("%3.2f"%(FARsnow_not_clouds*100)).rjust(5, ' '),
                    ("%3.2f"%(LOSTsnow*100)).rjust(5, ' '),
                    ("%3.2f"%(MISScloudy*100)).rjust(5, ' '),
                    ("%3.2f"%(MISScloudyOverSnow*100)).rjust(5, ' '),
                    test_name)
            elif clear_test and num_passed > 0:
                all_out_text += "N: %s POD-clear: %s FAR-clear: %s MISS-cloudy: %s %s \n"%(
                    str(num_passed).rjust(8, ' '),
                    ("%3.2f"%(PODclear*100)).rjust(5, ' '),
                    ("%3.2f"%(FARclear*100)).rjust(5, ' '),
                    ("%3.2f"%(MISScloudy*100)).rjust(5, ' '),
                    test_name)
            elif num_passed > 0:
                all_out_text_i = "N: %s POD-cloudy: %s FAR-cloudy: %s MISS-clear: %s MISS-snow: %s %s \n"%(
                    str(num_passed).rjust(8, ' '),
                    ("%3.2f"%(PODcloudy*100)).rjust(5, ' '),
                    ("%3.2f"%(FARcloudy*100)).rjust(5, ' '),
                    ("%3.2f"%(MISSclear*100)).rjust(5, ' '),
                    ("%3.2f"%(MISSsnow*100)).rjust(5, ' '),
                    test_name)
                all_out_text += "N: %s PCy: %s Fy: %s mPCl: %s : %s (%s) \n"%(
                    str(num_passed).rjust(8, ' '),
                    ("%3.2f"%(PODcloudy*100)).rjust(5, ' '),
                    ("%3.2f"%(FARcloudy*100)).rjust(5, ' '),
                    ("%3.2f"%(MISSclear*100)).rjust(5, ' '),
                    ("%3.2f"%(MISSsnow*100)).rjust(5, ' '),
                    test_name)

            # print "%s test_bit:%d N: %d POD-cloudy: %3.2f FAR-cloudy: %3.2f POD-clear %3.2f FAR-clear %3.2f"%(
            #   var, bit_nr, Num, PODcloudy, FARcloudy , PODclear, FARclear)

    # print "should be zero", np.sum(np.logical_and(isCalipsoClear[all_pix],
    #                                             isCloudyPPS[all_pix])),
    # print np.sum(np.logical_and(isCalipsoClear[use],
    #                           isCloudyPPS[use]))

    outfile_h.write(all_out_text)
    Num = np.sum(use)
    part_nodata = nodata*1.0/(nodata + Num)
    # print "N: %d POD-cloudy: %3.2f FAR-cloudy: %3.2f POD-clear %3.2f"%(
    #   Num, PODcloudy, FARcloudy , PODclear)


ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files/"
            "global_modis_14th_created20170519/")
ROOT_DIR_GAC = (ADIR + "/DATA_MISC/reshaped_files/"
            "ATRAIN_RESULTS_GAC_oldctth/Reshaped_Files/noaa18/")
ROOT_DIR_GAC = (ADIR + "/DATA_MISC/reshaped_files/"
            "ATRAIN_RESULTS_GAC_oldctth_20161201/Reshaped_Files/noaa18/")
ROOT_DIR_GAC = (ADIR + "/NN_CTTH/data/validation_data/imager_gac_noaa19/")
ROOT_DIR_NPP1 = (ADIR + "/DATA_MISC/reshaped_files/"
                "ATRAIN_RESULTS_NPP_new_thr_alos_85/Reshaped_Files"
                "/npp/1km/2015/07/no_area/")
ROOT_DIR_NPP2 = (ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
                "ATRAIN_RESULTS_NPP_test_with_lower_37_emiss/Reshaped_Files"
                "/npp/1km/2015/07/")
ROOT_DIR_NPP2 = (ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
                "ATRAIN_RESULTS_NPP_watercloud/Reshaped_Files"
                "/npp/1km/2015/07/")
ROOT_DIR_NPP2 = (ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
                "ATRAIN_RESULTS_NPP_emiss37/Reshaped_Files"
                "/npp/1km/2015/07/")
ROOT_DIR_NPP2 = (ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
                "ATRAIN_RESULTS_NPP_fine_snow/Reshaped_Files"
                "/npp/1km/2015/07/")
ROOT_DIR_NPP2 = (ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
                "ATRAIN_RESULTS_NPP_C4_build15/Reshaped_Files"
                "/npp/1km/2015/07/")
ROOT_DIR_GAC9 = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_tuned_nnIMAGER/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_thr_wvp/Reshaped_Files/noaa18/5km/2009/*/"
ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_rttov12/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_percentile_thr/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_rolles_method/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_rolle_percentile/Reshaped_Files/noaa18/5km/2009/*/"
#ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_larger_bins/Reshaped_Files/noaa18/5km/2009/*/"
#ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_larger_bins_and_noise/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_rM_tP_lB_N_mM/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_defaults_5/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_5_cfg_adjusted/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_adjusted_2of/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_ts_noise/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_more_clouds/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_21/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_23/Reshaped_Files/noaa18/5km/2009/*/"
ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files_jenkins_gac/ATRAIN_RESULTS_GAC_t11t12_higher_after_problem_cases/Reshaped_Files/noaa18/5km/2009/*/"
ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files_jenkins_gac/ATRAIN_RESULTS_GAC_one_larger_refl/Reshaped_Files/noaa18/5km/2009/*/"
# ROOT_DIR_GAC = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC_v20142/Reshaped_Files/noaa18/5km/2009/*/"

# files = glob(ROOT_DIR + "Reshaped_Files_merged/eos2/1km/2010/*/*h5")
# files = glob(ROOT_DIR_GAC + "5km/20??/*/*/*h5")
# files = glob(ROOT_DIR_GAC  + "/*/*h5")

files = glob(ROOT_DIR_NPP2 + "*/*h5")
test_list_file = ADIR + "/git/acpg_develop/acpg/pges/cloudmask/pps_pge01_cmasktests.h"
# test_list_file = "pps_pge01_cmasktests.h"
TEST_NAMEFILE = open(test_list_file, 'r')
name_dict = {'cma_testlist0': {},
             'cma_testlist1': {},
             'cma_testlist2': {},
             'cma_testlist3': {},
             'cma_testlist4': {},
             'cma_testlist5': {}}
for line in TEST_NAMEFILE:
    if 'define SM_ACMG_' in line:
        list_of_line = line.split(' ')
        (name, list_nr, bit) = (list_of_line[1].replace('SM_ACMG_', ''), list_of_line[2].replace(', ', ''), list_of_line[3])
        var = 'cma_testlist' + list_nr
        name_dict[var][int(bit)] = name

match_calipsoPPS = CalipsoImagerTrackObject()
for filename in files:
    print filename
    match_calipsoPPS += readCaliopImagerMatchObj(filename)
    print "Scene %s"%(os.path.basename(filename))
use = match_calipsoPPS.imager.all_arrays['bt11micron']>-9
# use = match_calipsoPPS.imager.all_arrays['sunz']>95
basename_outfile = filename.split("/Reshaped_Files")[0]
basename_outfile = basename_outfile.split("/")[-1]
print basename_outfile
for illumination in ["dnt", "day", "night", "twilight"]:
    for surface_type in [ "noemiss_coast", "noemiss_land", "land", "sea", "all", "coast", "emiss_land", "emiss_coast"]:
        # for mints, maxts in zip([190, 190, 240, 260, 260, 270, 280, 290, 280, 300, 320, 190, 270],
        #                       [380, 240, 260, 280, 270, 280, 290, 300, 300, 320, 380, 270, 380 ]):
        for mints, maxts in zip([190, 190, 270], # 190, 260, 280, 290, 300],
                                [380, 270, 380]):#, 260, 280, 290, 300, 380 ]):
            pass
            print_common_stats(match_calipsoPPS, use, name_dict, mints, maxts, surface_type, illumination, basename_outfile=basename_outfile)
# for surface_type in [ "all", "land", "sea", ]:
#   for mints, maxts in zip([ 190],
#                           [ 380 ]):
#       print_common_stats(match_calipsoPPS, use, name_dict, mints, maxts, surface_type, basename_outfile=basename_outfile)
