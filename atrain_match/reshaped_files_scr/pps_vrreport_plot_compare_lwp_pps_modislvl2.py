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

import matplotlib.pyplot as plt
import matplotlib

from utils.get_flag_info import (get_land_coast_sea_info_pps2014,
                           get_calipso_low_medium_high_classification,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_low_clouds)
from utils.stat_util import (my_iqr, my_hist, my_rms, my_mae)

from scipy.stats import kurtosis, skewtest, skew, mode, kurtosis
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 16})
from my_dir import ADIR
def get_land_sea(match_obj):

    use = np.logical_and(match_obj.modis.all_arrays["lwp"]>=0,
                         match_obj.imager.all_arrays["cpp_lwp"]>=0)
    use = np.logical_and(use,
                         match_obj.imager.all_arrays["cpp_phase"]==1)
    use = np.logical_and(match_obj.modis.all_arrays["cloud_phase"]==1,
                         use)
    retv = get_land_coast_sea_info_pps2014(
        match_obj.imager.all_arrays['cloudtype_conditions'])
    (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag) =  retv
    use_land = np.logical_and(use, land_flag)
    use_coast = np.logical_and(use, coast_flag)
    use_sea = np.logical_and(use, sea_flag)
    return use_land, use_sea, use_coast, use


def my_label(data, use):
    label = (#"{:s}\n"
             "bias = {:3.1f}\n"
             "RMS = {:3.1f}\n"
             "median = {:3.1f}\n"
             "MAE = {:3.1f}\n"
             "IQR = {:3.1f}\n"
             "PE>10 = {:3.1f}\n"
             "N = {:d}\n".format(
                 #text
                 np.mean(data[use]),
                 my_rms(data[use]),
                 np.median(data[use]),
                 my_mae(data[use]),
                 my_iqr(data[use]),
                 len(data[use][np.abs(data[use]) > 10])*100.0/len(data[use]),
                 len(data[use])
                            ))
    return label


def plot_all(match_obj, density, my_str=""):
    use_land, use_sea, use_coast, use  = get_land_sea(match_obj)
    x = match_obj.modis.all_arrays["lwp"]
    y = match_obj.imager.all_arrays["cpp_lwp"]*density
    error = y - x
    print x
    print np.max(y)
    fig = plt.figure(figsize=(15, 11))
    ax = fig.add_subplot(221)
    hist_heights, x_, gaussian = my_hist(error, use, return_also_corresponding_gaussian=True)
    ax.fill(x_, gaussian, color='silver')
    plt.plot(x_, hist_heights, 'b-', label=my_label(error, use))
    plt.title("All")
    ax.set_xlim([-500, 600])
    ax.set_ylim([00, 14])
    ax.set_ylabel("percent")
    plt.legend()
    ax = fig.add_subplot(222)
    hist_heights, x_, gaussian = my_hist(error, use_land, return_also_corresponding_gaussian=True)
    ax.fill(x_, gaussian, color='silver')
    plt.plot(x_, hist_heights, 'b-', label=my_label(error, use_land))
    plt.title("Land")
    ax.set_xlim([-500, 600])
    ax.set_ylim([00, 14])
    plt.legend()
    ax = fig.add_subplot(223)
    hist_heights, x_ , gaussian= my_hist(error, use_sea, return_also_corresponding_gaussian=True)
    ax.fill(x_, gaussian, color='silver')
    plt.plot(x_, hist_heights, 'b-', label=my_label(error, use_sea))
    plt.title("Sea")
    ax.set_xlim([-500, 600])
    ax.set_ylim([00, 14.00])
    ax.set_ylabel("percent")
    ax.set_xlabel("PPS LWP - MODIS-C6 LWP $g/m^2$")
    plt.legend()
    ax = fig.add_subplot(224)
    hist_heights, x_, gaussian = my_hist(error, use_coast, return_also_corresponding_gaussian=True)
    ax.fill(x_, gaussian, color='silver')
    plt.plot(x_, hist_heights, 'b-', label=my_label(error, use_coast))
    plt.title("Coast")
    ax.set_xlim([-500, 600])
    ax.set_ylim([00, 14.00])
    ax.set_xlabel("PPS LWP - MODIS-C6 LWP $g/m^2$")
    plt.legend()
    #plt.plot(x[use], y[use], 'b.', alpha=0.002)
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_comp_lwp_modis.png", bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_comp_lwp_modis.pdf", bbox_inches='tight')
    plt.show()
    fig = plt.figure(figsize=(15, 11))
    ax = fig.add_subplot(221)
    from trajectory_plotting import plot_satellite_trajectory
    import config
    plot_satellite_trajectory(match_obj.calipso.all_arrays["longitude"][use],
                            match_obj.calipso.all_arrays["latitude"][use],
                            ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_lwp_modis_lvl2_dist",
                            config.AREA_CONFIG_FILE,
                            fig_type=['png'])
    from histogram_plotting import distribution_map
    distribution_map(match_obj.calipso.all_arrays["longitude"][use],
                     match_obj.calipso.all_arrays["latitude"][use])
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_lwp_modis_lvl2_dist", bbox_inches='tight')

def get_ca_object_nn_ctth_modis_lvl2():
    ROOT_DIR = (
        ADIR + "/DATA_MISC/"
        "reshaped_files_validation_2018/"
        "global_modis_v2018_created20180920/"
        "Reshaped_Files_merged_caliop/eos2/1km/2010/*/*%s*cali*h5")
    files = glob(ROOT_DIR%("20100201"))
    files = files + glob(ROOT_DIR%("20100401"))
    files = files + glob(ROOT_DIR%("20100601"))
    files = files + glob(ROOT_DIR%("20100801"))
    files = files + glob(ROOT_DIR%("20101001"))
    files = files + glob(ROOT_DIR%("20101201"))
    #ROOT_DIR = (
    #    ADIR + "/DATA_MISC/"
    #    "reshaped_files_validation_2018/"
    #    "global_gac_v2018_created20180918/"
    #    "Reshaped_Files/noaa18/5km/*/*cali*h5")
    density = 1.0
    if "v2014" in ROOT_DIR:
        density = 1e3

    match_obj = CalipsoImagerTrackObject()
    print ROOT_DIR

    for filename in files:
        print filename
        match_obj += readCaliopImagerMatchObj(filename)
    return match_obj, density

def do_the_printing():

    match_obj_cali, density = get_ca_object_nn_ctth_modis_lvl2()


    plot_all(match_obj_cali, density, "")

if __name__ == "__main__":
    do_the_printing()
