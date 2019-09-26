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
"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       ppsMatch_Imager_CalipsoObject)
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
from utils.get_flag_info import get_calipso_clouds_of_type_i
from utils.get_flag_info import (get_semi_opaque_info_pps2014,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_low_clouds)
from my_dir import ADIR
def make_boxplot(match_calipso, name):
    low_clouds = get_calipso_low_clouds(match_calipso)
    high_clouds = get_calipso_high_clouds(match_calipso)
    medium_clouds = get_calipso_medium_clouds(match_calipso)
    height_c = (1000*match_calipso.calipso.all_arrays['layer_top_altitude'][:, 0] -
                match_calipso.calipso.all_arrays['elevation'])
    height_c1 = (1000*match_calipso.calipso.all_arrays['layer_top_altitude'][:, 0] -
                 match_calipso.calipso.all_arrays['elevation'])
    height_c2 = (1000*match_calipso.calipso.all_arrays['layer_top_altitude'][:, 1] -
                 match_calipso.calipso.all_arrays['elevation'])
    height_c3 = (1000*match_calipso.calipso.all_arrays['layer_top_altitude'][:, 2] -
                 match_calipso.calipso.all_arrays['elevation'])
    height_c4 = (1000*match_calipso.calipso.all_arrays['layer_top_altitude'][:, 3] -
                 match_calipso.calipso.all_arrays['elevation'])
    height_pps = match_calipso.imager.all_arrays['ctth_height']
    bias_1 = height_pps - height_c1
    bias_2 = height_pps - height_c2
    bias_3 = height_pps - height_c3
    bias_4 = height_pps - height_c4
    thin = np.logical_and(match_calipso.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.30,
                          match_calipso.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0)
    very_thin = np.logical_and(match_calipso.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.10,
                          match_calipso.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0)
    thin_top = np.logical_and(match_calipso.calipso.all_arrays['number_layers_found']>1, thin)
    thin_1_lay = np.logical_and(match_calipso.calipso.all_arrays['number_layers_found']==1, thin)
    # height_c[thin_top] = height_c2[thin_top]
    # height_c[np.abs(bias_1)<np.abs(bias_2)] = height_c1[np.abs(bias_1)<np.abs(bias_2)]
    # height_c[np.abs(bias_2)<np.abs(bias_1)] = height_c2[np.abs(bias_2)<np.abs(bias_1)]
    # bias = height_pps - height_c
    # height_c[np.abs(bias_3)<np.abs(bias)] = height_c3[np.abs(bias_3)<np.abs(bias)]
    # height_c[~thin_top] = height_c1[~thin_top]
    # height_c[thin_top] = height_c2[thin_top]



    use = np.logical_and(height_pps > - 1,
                         match_calipso.calipso.all_arrays['layer_top_altitude'][:, 0]>=0)
    use = np.logical_and(height_pps < 45000, use)

    low = np.logical_and(low_clouds, use)
    medium = np.logical_and(medium_clouds, use)
    high = np.logical_and(high_clouds, use)
    c_all = np.logical_or(high, np.logical_or(low, medium))
    high_very_thin = np.logical_and(high, very_thin)
    high_thin = np.logical_and(high, np.logical_and(~very_thin, thin))
    high_thick = np.logical_and(high, ~thin)
    # print "thin, thick high", np.sum(high_thin), np.sum(high_thick)
    bias = height_pps - height_c
    limit = np.percentile(bias[use], 5)
    # print limit
    abias = np.abs(bias)
    MAE = np.mean(abias[c_all])
    # abias[abias>2000]=2000
    print name.ljust(30, " "), "%3.1f"%(np.mean(abias[c_all])), "%3.1f"%(np.mean(abias[low])), "%3.1f"%(np.mean(abias[medium])), "%3.1f"%(np.mean(abias[high])), "%3.1f"%(limit)

    c_all = np.logical_or(np.logical_and(~very_thin, high), np.logical_or(low, medium))
    number_of = np.sum(c_all)

    # print name.ljust(30, " "), "%3.1f"%(np.sum(abias[c_all]<250)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<500)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<1000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<1500)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<2000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<3000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<4000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<5000)*100.0/number_of)
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    fig = plt.figure(figsize = (6, 9))
    ax = fig.add_subplot(111)
    plt.xticks(rotation=70)
    # plt.tight_layout()
    # plt.subplots_adjust(left=0.2)
    # plt.subplots_adjust(left=10, bottom=10, right=10, top=10, wspace=0, hspace=0)

    ax.fill_between(np.arange(0, 8), -500, 500, facecolor='green', alpha=0.6)
    ax.fill_between(np.arange(0, 8), -1000, 1000, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(0, 8), -1500, 1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0, 8), 2000, 15000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0, 8), -2000, -15000, facecolor='red', alpha=0.2)
    for y_val in [-5, -4, -3, -2, 2, 3, 4, 5]:
        plt.plot(np.arange(0, 8), y_val*1000 + 0*np.arange(0, 8), ':k')
        plt.plot(np.arange(0, 8), -10*1000 + 0*np.arange(0, 8), ':k')
    plt.plot(np.arange(0, 8), 0 + 0*np.arange(0, 8), 'k')
    # plt.boxplot([bias[low], bias[medium], bias[high], bias[high_thick], bias[high_thin], bias[high_very_thin]], whis=[5, 95], sym='',
    #           labels=["low", "medium", "high-all", "high-thick\n od>0.3", "high-thin \n 0.1<od<0.3", "high-vthin\n od<0.1"], showmeans=True)
    plt.boxplot([bias[low], bias[medium], bias[high]], whis=[5, 95], sym='',
                labels=["low", "medium", "high"], showmeans=True)
    ax.set_ylim(-14000, 8000)
    title_name = "CTTH-2018  "
    if "2014" in name:
        title_name = "CTTH-2014  "

    plt.title("%s MAE = %3.0f"%(title_name, MAE))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_LAPSE_RATE_INVESTIGATION/ctth_box_plot_%s_5_95_filt.png"%(name))
    # plt.show()



def investigate_nn_ctth():
    ROOT_DIR_GAC_nnNina = (ADIR + "/DATA_MISC/reshaped_files/"
                       "ATRAIN_RESULTS_GAC_nnNina/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn = (ADIR + "/DATA_MISC/reshaped_files/"
                       "ATRAIN_RESULTS_GAC_nn21/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_old = (ADIR + "/DATA_MISC/reshaped_files/"
                        "ATRAIN_RESULTS_GAC_old/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_oldCTTH_12x12 = (ADIR + "/DATA_MISC/reshaped_files/"
                        "ATRAIN_RESULTS_GAC_oldCTTH_12x12/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_v2014 = (ADIR + "/DATA_MISC/reshaped_files/"
                        "ATRAIN_RESULTS_GAC_v2014/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_v2014_12x12 = (ADIR + "/DATA_MISC/reshaped_files/"
                                "ATRAIN_RESULTS_GAC_v2014_12x12/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn_new = (ADIR + "/DATA_MISC/reshaped_files/"
                           "ATRAIN_RESULTS_GAC_nn20161125/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn_imager = (ADIR + "/DATA_MISC/reshaped_files/"
                           "ATRAIN_RESULTS_GAC_nn20161130/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn_imager_tuned = (ADIR + "/DATA_MISC/reshaped_files/"
                           "ATRAIN_RESULTS_GAC_tuned_nnIMAGER/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn_imager1_tuned = (ADIR + "/DATA_MISC/reshaped_files/"
                           "ATRAIN_RESULTS_GAC_nnIMAGER1/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn_imager_with_gac = (ADIR + "/DATA_MISC/reshaped_files/"
                           "ATRAIN_RESULTS_GAC_IMAGER_with_gac/Reshaped_Files/noaa18/")
    re_name = re.compile("_RESULTS_GAC_(\w+)\/")
    caobj_dict = {}
    for ROOT_DIR, name in zip([ROOT_DIR_GAC_nnNina, ROOT_DIR_GAC_nn, ROOT_DIR_GAC_old, ROOT_DIR_GAC_v2014,
                               ROOT_DIR_GAC_nn_new, ROOT_DIR_GAC_nn_imager,
                               ROOT_DIR_GAC_nn_imager_tuned, ROOT_DIR_GAC_nn_imager1_tuned,
                               ROOT_DIR_GAC_v2014_12x12,
                               ROOT_DIR_GAC_oldCTTH_12x12,
                               ROOT_DIR_GAC_nn_imager_with_gac],
                              ["gac_nnLessIsMore", "gac_nn21", "gac_CTTHold", "gac_2014",
                               "gac_nn20161125", "gac_nnIMAGER",
                               "gac_nnIMAGER_tuned", "gac_nnIMAGER1_tuned",
                               "gac_v2014_12x12", "gac_CTTHold_12x12", "gac_17var_modis_noaa19"]):
        files = glob(ROOT_DIR + "5km/2009/*/*/*h5")
        match_calipso = CalipsoImagerTrackObject()
        for filename in files:
            match_calipso += readCaliopImagerMatchObj(filename)
        caobj_dict[name] = match_calipso
        make_boxplot(match_calipso, name)
    # make_compare(caobj_dict['old'], caobj_dict['nn20161125'], 'test')
    # make_compare(caobj_dict['nn20161130'], caobj_dict['nn20161125'], 'test2')

def investigate_nn_ctth_viirs():
    ROOT_DIR_v2014 = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_C4_2014/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_v2018 = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_C4/Reshaped_Files/npp/1km/2015/07/*/")
    # ROOT_DIR_14bug_maia = (
    #   ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
    #   "NPP_FULL_ORBIT_2014/Reshaped_Files/")
    ROOT_DIR_14bug = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_v2014_before_ctthbug_correction/"
        "Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_v2014_old = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_v2014_bug_corrected_20170313/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_imager = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_nnIMAGER_20170313/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_imager_tuned = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_nnIMAGER_20170313_tuned/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_imager1 = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_nnIMAGER1_20170313/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_viirs = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_nnVIIRS_20170310/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_viirs_new = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_nnVIIRS_20170313_new/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_viirs_CLAY4 = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_nnVIIRS_20170310_CLAY4/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_viirs_lm = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_viirs_lm/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_imager_lm = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_imager_lm/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_nn_imager_wg = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_IMAGER_with_gac/Reshaped_Files/npp/1km/2015/07/*/")
    caobj_dict = {}
    for ROOT_DIR, name in zip(
            [ROOT_DIR_v2014, ROOT_DIR_v2018, ROOT_DIR_nn_imager_wg, ROOT_DIR_nn_imager1, ROOT_DIR_nn_imager, ROOT_DIR_nn_imager_tuned, ROOT_DIR_v2014_old, ROOT_DIR_14bug, ROOT_DIR_nn_viirs, ROOT_DIR_nn_viirs_CLAY4, ROOT_DIR_nn_viirs_new, ROOT_DIR_nn_viirs_lm, ROOT_DIR_nn_imager_lm],
            ["npp_CTTH-2014-C4", "npp_CTTH-2018-C4", "npp_CTTHnn_IMAGER_with_gac", "npp_CTTHnn_IMAGER1", "npp_CTTHnn_IMAGER", "npp_CTTHnn_IMAGER_tuned", "npp_CTTHv2014", "npp_CTTHv2014_buggy", "npp_CTTHnn_VIIRS", "npp_CTTHnn_VIIRS_C4", "npp_CTTHnn_VIIRS_tuned", "npp_nnVIIRS_LessIsMore", "npp_nnIMAGER_LessIsMore"]):
        print ROOT_DIR
        files = glob(ROOT_DIR + "*.h5")
        print files
        match_calipso = CalipsoImagerTrackObject()
        for filename in files:
            # print filename
            match_calipso += readCaliopImagerMatchObj(filename)
        caobj_dict[name] = match_calipso
        make_boxplot(match_calipso, name)

def investigate_nn_ctth_modis_may():
   # november
    # may
    """
    ROOT_DIR_MODIS_nn = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_MAY/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_viirs = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_MAY_nnviirs_20161205/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_mersi2 = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_MAY_nnmersi2_20161206/Reshaped_Files/merged/")
    """
    ROOT_DIR_MODIS_old = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "global_modis_14th_created20161108/Reshaped_Files/merged/")

    caobj_dict = {}
    for ROOT_DIR, name in zip(
            [ROOT_DIR_MODIS_nn, ROOT_DIR_MODIS_nn_viirs,
             ROOT_DIR_MODIS_nn_mersi2, ROOT_DIR_MODIS_old],
            ["modis_november_CTTHnn_IMAGER",
             "modis_november_CTTHnn_VIIRS",
             "modis_november_CTTHnn_MERSI2",
             "modis_november_CTTHold"]):
        print name
        files = glob(ROOT_DIR + "/*11*.h5")
        match_calipso = CalipsoImagerTrackObject()
        for filename in files:
            # print filename
            match_calipso += readCaliopImagerMatchObj(filename)
        caobj_dict[name] = match_calipso
        make_boxplot(match_calipso, name)
    # make_compare(caobj_dict["modis_nn18var"],
    #            caobj_dict["modis_CTTHold"],
    #            'compare_modis')


def investigate_nn_ctth_modis_november():
    # november
    ROOT_DIR_MODIS_nn_imager = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_viirs = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER_nnviirs_20161205/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_mersi2 = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER_AEROSOL/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_viirs_tuned = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER_nnVIIRS_20170315/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_mersi2_tuned = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER_nnMERSI2/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_imager_tuned = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER_nnIMAGER_20170315/Reshaped_Files/merged/")

    ROOT_DIR_MODIS_old = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "global_modis_14th_created20161108/Reshaped_Files/merged/")

    caobj_dict = {}
    for ROOT_DIR, name in zip(
            [ROOT_DIR_MODIS_nn_imager,
             ROOT_DIR_MODIS_nn_viirs,
             ROOT_DIR_MODIS_nn_mersi2,
             ROOT_DIR_MODIS_old,
             ROOT_DIR_MODIS_nn_viirs_tuned,
             ROOT_DIR_MODIS_nn_mersi2_tuned,
             ROOT_DIR_MODIS_nn_imager_tuned],
            ["modis_nov_nnIMAGER",
             "modis_nov_nnVIIRS",
             "modis_nov_nnMERSI2",
             "modis_nov_CTTHold",

             "modis_nov_nnVIIRS_tuned",
             "modis_nov_nnMERSI2_tuned",
             "modis_nov_nnIMAGER_tuned",
]):
        files = glob(ROOT_DIR + "/*11*.h5")
        match_calipso = CalipsoImagerTrackObject()
        for filename in files:
            # print filename
            match_calipso += readCaliopImagerMatchObj(filename)
        caobj_dict[name] = match_calipso
        make_boxplot(match_calipso, name)
    # make_compare(caobj_dict["modis_nn18var"],
    #            caobj_dict["modis_CTTHold"],
    #            'compare_modis')

if __name__ == "__main__":
    # investigate_nn_ctth()
    # investigate_nn_ctth_modis_november()
    investigate_nn_ctth_viirs()

