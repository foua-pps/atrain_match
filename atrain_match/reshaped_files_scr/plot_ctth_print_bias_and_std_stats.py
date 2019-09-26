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
                            DataObject,
                            CloudsatImagerTrackObject,
                            readCloudsatImagerMatchObj,
                            CalipsoImagerTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       ppsMatch_Imager_CalipsoObject)
from utils.get_flag_info import get_calipso_clouds_of_type_i_feature_classification_flags_one_layer
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import kurtosis, skewtest, skew, mode, kurtosis
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 18})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from my_dir import ADIR
from plot_ctth_bias_distributions import PlotAndDataObject, extract_data

tag_dict = {"old": "(a) PPS-v2014",
            "mlvl2": "(b) MODIS-C6",
            "nnant": "(c) NN-AVHRR",
            "pps": "(c) NN-AVHRR",
            "nna1nt": "(d) NN-AVHRR1",
            "nnvnt": "(e) NN-VIIRS",
            "nnm2nt": "(f) NN-MERSI-2",
            "nnmint": "(h) NN-MetImage",
            "nnmintnco2": "(g) NN-MetImage-$NoCO_2$",

}
# https://stats.stackexchange.com/questions/278237/half-sample-mode-estimate-of-sample-of-weighted-data
def half_sample_mode(x, already_sorted=False):
    if len(x) < 3:
        return np.mean(x)
    if already_sorted:
        sorted_x = x # No need to sort
    else:
        sorted_x = np.sort(x)
    half_idx = int((len(x) + 1) / 2) # Round up to include the middle value, in the case of an odd-length array

    # Calculate all interesting ranges that span half of all data points
    ranges = sorted_x[ - half_idx:] - sorted_x[:half_idx]
    smallest_range_idx = np.argmin(ranges)

    # Now repeat the procedure on the half that spans the smallest range
    x_subset = sorted_x[smallest_range_idx : (smallest_range_idx + half_idx)]
    return half_sample_mode(x_subset, already_sorted=True)

def my_iqr(data):
    return np.percentile(data, 75) - np.percentile(data, 25)

def my_mode(bias):
    bmin = -40
    bmax = 40
    delta_h = 100.0
    bins = np.arange(bmin*1000, bmax*1000, delta_h)
    hist_heights, bins = np.histogram(bias, bins=bins)
    n_pix = len(bias)
    hist_heights = hist_heights*100.0/n_pix
    maxind = np.argmax(hist_heights)
    maxind2 = len(hist_heights) - 1 - np.argmax(hist_heights[::-1])
    if maxind != maxind2:
        print maxind, maxind2
        raise ValueError
    mode =bins[maxind] + delta_h*0.5
    return mode


def print_for_one(plt_obj, compare, truth='height_c'):
    x = getattr(plt_obj, truth)
    y = getattr(plt_obj, compare)
    use = np.logical_and(x>=0, getattr(plt_obj, 'use_all'))
    use_low = np.logical_and(use, plt_obj.low_clouds)
    use_medium = np.logical_and(use, plt_obj.medium_clouds)
    use_high = np.logical_and(use, plt_obj.high_clouds)

    name_conversion = {'old': 'PPS-v2014',
                       'mlvl2': 'MODIS-C6',
                       'pps': 'NN-AVHRR',
                       'nnant': 'NN-AVHRR',
                       'nna1nt': 'NN-AVHRR1',
                       'nnvnt': 'NN-VIIRS',
                       'nnm2nt': 'NN-MERSI-2',
                       'nnmintnco2': 'NN-MetImage-NoCO$_{\\text 2}$',
                       'nnmint': 'NN-MetImage'}

    compare_name = name_conversion[compare.split('_')[1]]

    bias = y - x#+1465
    AE = np.abs(bias)
    std = np.std(bias[use])

    use_i = use
    print "%s & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %3.1f  \\\\"%(
        compare_name,
        # np.sum(AE[use_i]<=1000)*100.0/len(AE[use_i]),
        np.mean(AE[use_i]),
        np.percentile(bias[use_i], 75) - np.percentile(bias[use_i], 25),
        np.sqrt(np.mean((bias[use_i]**2))),
        # np.sqrt(np.mean((bias[np.logical_and(use_i, AE<=2000)]**2))),
        # np.mean((AE[np.logical_and(use_i, AE>2000)])),
        np.std(bias[use_i]),
        np.sum(AE[use_i] > 250)*100.0/len(AE[use_i]),
        np.sum(AE[use_i] > 500)*100.0/len(AE[use_i]),
        np.sum(AE[use_i] > 1000)*100.0/len(AE[use_i]),
        np.sum(AE[use_i] > 2000)*100.0/len(AE[use_i]),
        np.sum(AE[use_i] > 2500)*100.0/len(AE[use_i]),
        np.median(bias[use_i]),
        # my_mode(bias[use_i]),
        half_sample_mode(bias[use_i]),
        np.mean(bias[use_i]),
        skew(bias[use_i])#,
        # kurtosis(bias[use_i])
        )


    """
    print "%s & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d &  %d & %d & %d & %d & %d &   \\\\"%(
        compare_name,
        np.mean(AE[use_low]),
        np.mean(AE[use_medium]),
        np.mean(AE[use_high]),
        np.mean(bias[use_low]),
        np.mean(bias[use_medium]),
        np.mean(bias[use_high]),
        np.median(bias[use_low]),
        np.median(bias[use_medium]),
        np.median(bias[use_high]),
        np.percentile(bias[use_low], 75) - np.percentile(bias[use_low], 25),
        np.percentile(bias[use_medium], 75) - np.percentile(bias[use_medium], 25),
        np.percentile(bias[use_high], 75) - np.percentile(bias[use_high], 25),
        np.sqrt(np.mean((bias[use_low]**2))),
        np.sqrt(np.mean((bias[use_medium]**2))),
        np.sqrt(np.mean((bias[use_high]**2))))
    """

    if plt_obj.cflag is not None:
        pass
        """
        print " %s & %d &"%(compare_name, np.sum(AE[use]>500)*100.0/len(AE[use])),
        for ind in range(8):
            use_i = np.logical_and(
                use,
                get_calipso_clouds_of_type_i_feature_classification_flags_one_layer(plt_obj.cflag, ind))
            if np.sum(use_i) > 0:
                print "%d &"%(np.sum(AE[use_i]>500)*100.0/len(AE[use_i])),
                # print " %d &"%(np.mean(AE[use_i ])),
            else:
               print "- &",
        print "\\\\"
        """
        """
        print " %s %d &"%(compare_name, np.mean(AE[use])),
        for ind in range(8):
            use_i = np.logical_and(
                use,
                get_calipso_clouds_of_type_i_feature_classification_flags_one_layer(plt_obj.cflag, ind))
            if np.sum(use_i) > 0:
                print "%d &"%(100.0*np.sum(use_i)/np.sum(use)),
                print " %d &"%(np.mean(AE[use_i ])),
            else:
               print "- &",
        print "\\\\"
        """
        """
        print " %s %d &"%(compare_name, np.median(bias[use])),
        for ind in range(8):
            use_i = np.logical_and(
                use,
                get_calipso_clouds_of_type_i_feature_classification_flags_one_layer(plt_obj.cflag, ind))
            if np.sum(use_i) > 0:
                print " %d &"%(np.median(bias[use_i ])),
            else:
                print "- &",
        print "\\\\"
        """
    """
    else:
        print " %s & %d& %d& %d &  %d& %d& %d \\\\ "%(compare_name, #np.mean(AE[use]),
                                                      np.mean(AE[use_low]),
                                                      np.mean(AE[use_medium]),
                                                      np.mean(AE[use_high]),
                                                      np.median(bias[use_low]),
                                                      np.median(bias[use_medium]),
                                                      np.median(bias[use_high]))


    print "N   : %s %d %d %d  %d "%(compare_name,
                                    len(AE[use]),
                                    len(AE[use_low]),
                                    len(AE[use_medium]),
                                    len(AE[use_high]))


    print "MAE : %s & %3.1f& %3.1f& %3.1f & %3.1f\\\\ "%(compare_name, np.mean(AE[use]),
                                                np.mean(AE[use_low]),
                                                np.mean(AE[use_medium]),
                                                np.mean(AE[use_high]))
    print "Bias: %s & %3.1f & %3.1f & %3.1f &  %3.1f\\\\ "%(compare_name, np.mean(bias[use]),
                                                np.mean(bias[use_low]),
                                                np.mean(bias[use_medium]),
                                                np.mean(bias[use_high]))

    print "Std :%s & %3.1f & %3.1f & %3.1f  & %3.1f\\\\ "%(compare_name, np.std(bias[use]),
                                               np.std(bias[use_low]),
                                               np.std(bias[use_medium]),
                                               np.std(bias[use_high]))

    """
    """
    print "%s & %d & %d & %d &  %d"%(compare_name, np.mean(bias[use]),
                                                np.mean(bias[use_low]),
                                                np.mean(bias[use_medium]),
                                                np.mean(bias[use_high])),

    print " & %d & %d & %d  & %d "%( np.std(bias[use]),
                                               np.std(bias[use_low]),
                                               np.std(bias[use_medium]),
                                               np.std(bias[use_high])),


    print " & %3.1f & %3.1f & %3.1f &  %3.1f \\\\"%(skew(bias[use]),
                                     skew(bias[use_low]),
                                     skew(bias[use_medium]),
                                     skew(bias[use_high]))
    """
    # print " & %3.1f & %3.1f & %3.1f &  %3.1f \\\\"%(skewtest(bias[use])[1],
    #                                skewtest(bias[use_low])[1],
    #                                skewtest(bias[use_medium])[1],
    #                                skewtest(bias[use_high])[1])
    #
    """
    print "mode %s & %d & %d & %d &  %d"%(compare_name, my_mode(bias[use]),
                                     my_mode(bias[use_low]),
                                     my_mode(bias[use_medium]),
                                     my_mode(bias[use_high]))
    print "%s & %d & %d & %d &  %d \\\\"%(compare_name, np.median(bias[use]),
                                                        np.median(bias[use_low]),
                                                     np.median(bias[use_medium]),
                                                        np.median(bias[use_high]))
   """
    """
    lim = 4500
    print "%s & %3.1f & %3.1f & %3.1f &  %3.1f"%(
        compare_name,
        np.sum(np.logical_and(use, np.abs(bias) > lim))*100.0/np.sum(use),
        np.sum(np.logical_and(use_low, np.abs(bias) > lim))*100.0/np.sum(use_low),
        np.sum(np.logical_and(use_medium, np.abs(bias) > lim))*100.0/np.sum(use_medium),
        np.sum(np.logical_and(use_high, np.abs(bias) > lim))*100.0/np.sum(use_high))

    """
    # print " precentile %3.4f & %d & %d & %d  & %d\\\\ "%(ind,
    #                                                    np.percentile(bias[use], ind),
    #                                                    np.percentile(bias[use_low], ind),
    #                                                    np.percentile(bias[use_medium], ind),
    #                                                    np.percentile(bias[use_high], ind))

    # IQR
    """
    print "%s & %d & %d& %d  & %d\\\\ "%(compare_name,
        np.percentile(bias[use], 75) -
        np.percentile(bias[use], 25),
        np.percentile(bias[use_low], 75) -
        np.percentile(bias[use_low], 25),
        np.percentile(bias[use_medium], 75) -
        np.percentile(bias[use_medium], 25),
        np.percentile(bias[use_high], 75) -
        np.percentile(bias[use_high], 25)
    )
    """
    # MAD
    # print "%s & %d & %d& %d  & %d\\\\ "%(
    #   compare_name,
    #   np.median(np.abs(bias[use] - np.median(bias[use]))),
    #   np.median(np.abs(bias[use_low] - np.median(bias[use_low]))),
    #   np.median(np.abs(bias[use_medium] - np.median(bias[use_medium]))),
    #   np.median(np.abs(bias[use_high] - np.median(bias[use_high])))
    # )

    """
        lim = np.percentile(bias[use_low], ind)
        print np.mean(bias[np.logical_and(use_low, bias>lim)]),
        print np.std(bias[np.logical_and(use_low, bias<lim)])
    print "satz %3.2f, %3.2f"%(np.min(plt_obj.satz[use]),
                               np.max(plt_obj.satz[use]) )
    for lim in [1500, 4000]:
        print "lim", lim
        print np.mean(bias[np.logical_and(use_low, bias>lim)]),
        print np.mean(bias[np.logical_and(use_low, bias<=lim)]),
        print np.std(bias[np.logical_and(use_low, bias<=lim)])
        print " %d, %d, %3.2f"%(
            np.sum(np.logical_and(use_low, bias > lim)),
            np.sum(use_low),
            np.sum(np.logical_and(use_low, bias > lim))*100.0/np.sum(use_low))
        print " %d, %d, %3.2f"%(
            np.sum(np.logical_and(use_low, np.abs(bias) > lim)),
            np.sum(use_low),
            np.sum(np.logical_and(use_low, np.abs(bias) > lim))*100.0/np.sum(use_low))
    print "part of all data with error above 4km %3.2f %3.2f %3.2f %3.2f"%(
        np.sum(np.logical_and(use, np.abs(bias) > lim))*100.0/np.sum(use),
        np.sum(np.logical_and(use_low, np.abs(bias) > lim))*100.0/np.sum(use_low),
        np.sum(np.logical_and(use_medium, np.abs(bias) > lim))*100.0/np.sum(use_medium),
        np.sum(np.logical_and(use_high, np.abs(bias) > lim))*100.0/np.sum(use_high))
    """
def print_all(plt_obj_cali_new, plt_obj_csat_new, month):

    print "CALIOP height", month
    print_for_one(plt_obj_cali_new, 'height_old')
    print_for_one(plt_obj_cali_new, 'height_mlvl2')
    print_for_one(plt_obj_cali_new, 'height_pps')
    print_for_one(plt_obj_cali_new, 'height_nnvnt')
    print_for_one(plt_obj_cali_new, 'height_nnm2nt')
    print_for_one(plt_obj_cali_new, 'height_nnmintnco2')
    print_for_one(plt_obj_cali_new, 'height_nnmint')
    print_for_one(plt_obj_cali_new, 'height_nna1nt')

    # print "CALIOP pressure", month
    # print_for_one(plt_obj_cali_new, 'pressure_old', truth='pressure_c')
    # print_for_one(plt_obj_cali_new, 'pressure_mlvl2', truth='pressure_c')
    # print_for_one(plt_obj_cali_new, 'pressure_nnant', truth='pressure_c')
    # print_for_one(plt_obj_cali_new, 'pressure_nnvnt', truth='pressure_c')
    # print_for_one(plt_obj_cali_new, 'pressure_nnm2nt', truth='pressure_c')
    # print_for_one(plt_obj_cali_new, 'pressure_nnmintnco2', truth='pressure_c')
    # print_for_one(plt_obj_cali_new, 'pressure_nnmint', truth='pressure_c')
    # print_for_one(plt_obj_cali_new, 'pressure_nna1nt', truth='pressure_c')

    print "CloudSat height", month
    print_for_one(plt_obj_csat_new, 'height_old')
    print_for_one(plt_obj_csat_new, 'height_mlvl2')
    print_for_one(plt_obj_csat_new, 'height_pps')
    print_for_one(plt_obj_csat_new, 'height_nnvnt')
    print_for_one(plt_obj_csat_new, 'height_nnm2nt')
    print_for_one(plt_obj_csat_new, 'height_nnmintnco2')
    print_for_one(plt_obj_csat_new, 'height_nnmint')
    print_for_one(plt_obj_csat_new, 'height_nna1nt')

def get_plot_object_nn_ctth_modis_lvl2_cloudsat(month):
    day_str="01st"
    ROOT_DIR = (
        ADIR + "/DATA_MISC/reshaped_files/"
        # "global_modis_%s_created20170519/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
        "global_modis_%s_created20180316/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
        # "global_modis_%s_created20170330/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
    # match_clsat = CloudsatImagerTrackObject()
    plt_obj = PlotAndDataObject()
    print ROOT_DIR%(day_str, month)
    files = glob(ROOT_DIR%(day_str, month))
    for filename in files:
        print filename
        match_clsat_new = readCloudsatImagerMatchObj(filename)
        plt_obj += extract_data(match_clsat_new, sat='cloudsat')
    return plt_obj

def get_plot_object_nn_ctth_modis_lvl2(month):
    day_str="01st"
    ROOT_DIR = (
        ADIR + "/DATA_MISC/reshaped_files/"
        # "global_modis_%s_created20170504/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")
        "global_modis_%s_created20180316/Reshaped_Files_merged_calipso_cbase/eos2/1km/2010/%s/*h5")
        # "global_modis_%s_created20170519/Reshaped_Files_merged_calipso_cbase/eos2/1km/2010/%s/*h5")
    plt_obj = PlotAndDataObject()
    print ROOT_DIR%(day_str, month)
    files = glob(ROOT_DIR%(day_str, month))
    for filename in files:
        print filename
        match_calipso_new = readCaliopImagerMatchObj(filename)
        plt_obj += extract_data(match_calipso_new, sat='calipso')
    return plt_obj

def do_the_printing():

    plt_obj_cali = PlotAndDataObject()
    plt_obj_csat = PlotAndDataObject()

    merged_months = ""
    for month in ["02", "04", "06", "08", "10", "12"]:
    # for month in ["08"]:
        merged_months += month
        plt_obj_cali_new = get_plot_object_nn_ctth_modis_lvl2(month)
        plt_obj_csat_new = get_plot_object_nn_ctth_modis_lvl2_cloudsat(month)
        print_all(plt_obj_cali_new, plt_obj_csat_new, month)
        plt_obj_cali += plt_obj_cali_new
        plt_obj_csat += plt_obj_csat_new

    print_all(plt_obj_cali, plt_obj_csat, merged_months)

if __name__ == "__main__":
    do_the_printing()

