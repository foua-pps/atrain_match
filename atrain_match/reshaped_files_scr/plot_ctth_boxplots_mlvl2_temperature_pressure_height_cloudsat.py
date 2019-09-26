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
from matchobject_io import (readCloudsatImagerMatchObj,
                            CloudsatImagerTrackObject)

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})

from utils.get_flag_info import (get_semi_opaque_info_pps2014)
from my_dir import ADIR
def make_boxplot(match_clsat, name, month):

    height_c = match_clsat.cloudsat.all_arrays['clsat_max_height']
    low_clouds = np.logical_and(height_c < 2500, height_c > - 9)
    medium_clouds = np.logical_and(height_c >= 2500, height_c <= 5000)
    high_clouds = np.logical_and(height_c > 5000, height_c > - 9)
    height_mlvl2 = match_clsat.modis.all_arrays['height']
    height_pps = match_clsat.imager.all_arrays['imager_ctth_m_above_seasurface']
    sunz = match_clsat.imager.all_arrays['sunz']
    use = np.logical_and(height_pps > - 1,
                         height_c >= 0)
    use = np.logical_and(height_pps < 45000, use)
    use = np.logical_and(use, height_mlvl2 > - 1)
    use = np.logical_and(use, height_mlvl2 < 45000)
    low = np.logical_and(low_clouds, use)
    medium = np.logical_and(medium_clouds, use)
    high = np.logical_and(high_clouds, use)
    over_high_ground = np.logical_and(use, match_clsat.cloudsat.all_arrays['elevation']>5000)
    c_all = use #np.logical_or(high, np.logical_or(low, medium))
    c_all_day = np.logical_and(use, sunz < 80)
    c_all_twilight = np.logical_and(use, np.logical_and(sunz <= 95, sunz >= 80))
    c_all_night = np.logical_and(use, sunz > 95)

    pps_bias = height_pps - height_c
    pps_abias = np.abs(pps_bias)
    pps_MAE = np.mean(pps_abias[c_all])
    pps_MAEd = np.mean(pps_abias[c_all_day])
    pps_MAEn = np.mean(pps_abias[c_all_night])
    pps_MAEt = np.mean(pps_abias[c_all_twilight])
    mlvl2_bias = height_mlvl2 - height_c
    mlvl2_abias = np.abs(mlvl2_bias)
    mlvl2_MAE = np.mean(mlvl2_abias[c_all])
    mlvl2_MAEd = np.mean(mlvl2_abias[c_all_day])
    mlvl2_MAEn = np.mean(mlvl2_abias[c_all_night])
    mlvl2_MAEt = np.mean(mlvl2_abias[c_all_twilight])
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    fig = plt.figure(figsize = (6, 9))
    ax = fig.add_subplot(111)
    plt.xticks(rotation=50)
    ax.fill_between(np.arange(0, 5), -500, 500, facecolor='green', alpha=0.6)
    ax.fill_between(np.arange(0, 5), -1000, 1000, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(0, 5), -1500, 1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0, 5), 2000, 15000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0, 5), -2000, -15000, facecolor='red', alpha=0.2)
    for y_val in [-5, -4, -3, -2, 2, 3, 4, 5]:
        plt.plot(np.arange(0, 5), y_val*1000 + 0*np.arange(0, 5), ':k', alpha=0.4)
        plt.plot(np.arange(0, 5), -10*1000 + 0*np.arange(0, 5), ':k', alpha=0.4)
    plt.plot(np.arange(0, 5), 0 + 0*np.arange(0, 5), ':k', alpha=0.4)
    bplot = ax.boxplot([pps_bias[c_all], mlvl2_bias[c_all]], whis=[5, 95], sym='',
                       labels=["%s \n MAE=%3.0f"%(name, pps_MAE), "modis_lvl2 \nMAE=%3.0f"%(mlvl2_MAE)],
                       showmeans=True, patch_artist=True)
    ax.set_ylim(-8000, 8000)
    for box in bplot['boxes']:
        box.set_facecolor('0.9')
    plt.title("Cloudsat PPS/MODIS-LVL2 \nHeight bias comparison %s"%(month))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_BOX_cloudsat/"
                "ctth_box_plot_csat_pps_and_lvl2modis_%s_5_95_filt.png"%(month))

    # LowMediumHigh
    fig = plt.figure(figsize = (6, 9))
    ax = fig.add_subplot(111)
    plt.xticks(rotation=50)
    ax.fill_between(np.arange(0, 8), -500, 500, facecolor='green', alpha=0.6)
    ax.fill_between(np.arange(0, 8), -1000, 1000, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(0, 8), -1500, 1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0, 8), 2000, 15000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0, 8), -2000, -15000, facecolor='red', alpha=0.2)
    for y_val in [-5, -4, -3, -2, 2, 3, 4, 5]:
        plt.plot(np.arange(0, 8), y_val*1000 + 0*np.arange(0, 8), ':k', alpha=0.4)
        plt.plot(np.arange(0, 8), -10*1000 + 0*np.arange(0, 8), ':k', alpha=0.4)
    plt.plot(np.arange(0, 8), 0 + 0*np.arange(0, 8), ':k', alpha=0.4)
    bplot = ax.boxplot([pps_bias[low], pps_bias[medium], pps_bias[high], pps_bias[over_high_ground]],
                       whis=[5, 95], sym='',
                       labels=["low <2.5km", "medium", "high>5km", "ground>5km"],
                       showmeans=True, patch_artist=True)
    ax.set_ylim(-8000, 8000)
    for box in bplot['boxes']:
        box.set_facecolor('0.9')
    plt.title("Cloudsat PPS-%s %s \nHeight bias comparison MAE= %3.0f"%(name, month, pps_MAE))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_BOX_cloudsat/"
                "ctth_box_plot_csat_modis_%s_%s_5_95_filt.png"%(name, month))

    # LowMediumHigh
    fig = plt.figure(figsize = (6, 9))
    ax = fig.add_subplot(111)
    plt.xticks(rotation=50)
    ax.fill_between(np.arange(0, 8), -500, 500, facecolor='green', alpha=0.6)
    ax.fill_between(np.arange(0, 8), -1000, 1000, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(0, 8), -1500, 1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0, 8), 2000, 15000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0, 8), -2000, -15000, facecolor='red', alpha=0.2)
    for y_val in [-5, -4, -3, -2, 2, 3, 4, 5]:
        plt.plot(np.arange(0, 8), y_val*1000 + 0*np.arange(0, 8), ':k', alpha=0.4)
        plt.plot(np.arange(0, 8), -10*1000 + 0*np.arange(0, 8), ':k', alpha=0.4)
    plt.plot(np.arange(0, 8), 0 + 0*np.arange(0, 8), ':k', alpha=0.4)
    bplot = ax.boxplot([mlvl2_bias[low], mlvl2_bias[medium], mlvl2_bias[high], mlvl2_bias[over_high_ground]],
                       whis=[5, 95], sym='',
                       labels=["low <2.5km", "medium", "high>5km", "ground>5km"],
                       showmeans=True, patch_artist=True      )
    ax.set_ylim(-8000, 8000)
    for box in bplot['boxes']:
        box.set_facecolor('0.9')
    plt.title("Cloudsat MODIS-LVL2 %s \nHeight bias comparison MAE= %3.0f"%(month, mlvl2_MAE))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_BOX_cloudsat/"
                "ctth_box_plot_csat_modis_lvl2_C6_%s_5_95_filt.png"%(month))


    # DNT
    fig = plt.figure(figsize = (16, 9))
    ax = fig.add_subplot(111)
    plt.xticks(rotation=50)
    ax.fill_between(np.arange(0, 8), -500, 500, facecolor='green', alpha=0.6)
    ax.fill_between(np.arange(0, 8), -1000, 1000, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(0, 8), -1500, 1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0, 8), 2000, 15000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0, 8), -2000, -15000, facecolor='red', alpha=0.2)
    for y_val in [-5, -4, -3, -2, 2, 3, 4, 5]:
        plt.plot(np.arange(0, 8), y_val*1000 + 0*np.arange(0, 8), ':k', alpha=0.4)
        plt.plot(np.arange(0, 8), -10*1000 + 0*np.arange(0, 8), ':k', alpha=0.4)
    plt.plot(np.arange(0, 8), 0 + 0*np.arange(0, 8), ':k', alpha=0.4)
    bplot = ax.boxplot([pps_bias[c_all_day], pps_bias[c_all_night], pps_bias[c_all_twilight],
                        mlvl2_bias[c_all_day], mlvl2_bias[c_all_night], mlvl2_bias[c_all_twilight]], whis=[5, 95], sym='',
                       labels=["day %s \n MAE=%3.0f"%(name, pps_MAEd),
                               "night %s \n MAE=%3.0f"%(name, pps_MAEn),
                               "twilight %s \n MAE=%3.0f"%(name, pps_MAEt),
                               "day modis_lvl2 \nMAE=%3.0f"%(mlvl2_MAEd),
                               "night modis_lvl2 \nMAE=%3.0f"%(mlvl2_MAEn),
                               "twilight modis_lvl2 \nMAE=%3.0f"%(mlvl2_MAEt)], showmeans=True, patch_artist=True)
    ax.set_ylim(-8000, 8000)
    for box in bplot['boxes']:
        box.set_facecolor('0.9')
    plt.title("Cloudsat PPS/MODIS-LVL2 \nHeight bias comparison %s"%(month))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_BOX_cloudsat/"
                "ctth_box_plot_dnt_csat_pps_and_lvl2modis_%s_5_95_filt.png"%(month))

    # VIOLIN
    fig = plt.figure(figsize = (6, 9))
    ax = fig.add_subplot(111)
    plt.xticks(rotation=50)
    # ax.fill_between(np.arange(0, 5), -500, 500, facecolor='green', alpha=0.6)
    # ax.fill_between(np.arange(0, 5), -1000, 1000, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(0, 5), -1500, 1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0, 5), 2000, 45000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0, 5), -2000, -45000, facecolor='red', alpha=0.2)
    ax.x_tick_list = [1, 3]
    ax.x_tick_labels = labels=["%s \n MAE=%3.0f"%(name, pps_MAE), "modis_lvl2 \nMAE=%3.0f"%(mlvl2_MAE)]

    for y_val in [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]:
        plt.plot(np.arange(0, 5), y_val*2500 + 0*np.arange(0, 5), ':k', alpha=0.4)
    bplot = ax.violinplot([pps_bias[c_all], mlvl2_bias[c_all]], positions=[1, 3],
                          widths=1.2, showextrema=False, showmedians=True)
    plt.setp(ax, xticks=[1, 3],
             xticklabels=["%s \n MAE=%3.0f"%(name, pps_MAE), "modis_lvl2 \nMAE=%3.0f"%(mlvl2_MAE)])
    # labels=["%s \n MAE=%3.0f"%(name, pps_MAE), "modis_lvl2 \nMAE=%3.0f"%(mlvl2_MAE)]
    ax.set_ylim(-10000, 10000)
    plt.title("Cloudsat PPS/MODIS-LVL2 \nHeight bias comparison %s"%(month))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_BOX_cloudsat/"
                "ctth_violin_csat_pps_and_lvl2modis_%s_5_95_filt.png"%(month))



def investigate_nn_ctth_modis_lvl2():
    # november

    ROOT_DIR_MODIS_nn_imager = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "global_modis_14th_created20170324/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")

#   ROOT_DIR_MODIS_old = (
#       ADIR + "/DATA_MISC/reshaped_files/"
#       "global_modis_14th_created20161108/Reshaped_Files/merged/*%s*h5")

    for month in [ "06", "09"]:#, "02", "03", "04", "05", "07", "08", "10", "11", "12", "01"]:
        for ROOT_DIR, name in zip(
                [ROOT_DIR_MODIS_nn_imager],
                 # ROOT_DIR_MODIS_old],
                ["nnIMAGER"]):
            # name = "%s_%s"%(name, month)
            print ROOT_DIR
            files = glob(ROOT_DIR%(month))
            match_clsat = CloudsatImagerTrackObject()
            for filename in files:
                # print filename
                match_clsat +=  readCloudsatImagerMatchObj(filename)
            make_boxplot(match_clsat, name, month )


if __name__ == "__main__":
    investigate_nn_ctth_modis_lvl2()

