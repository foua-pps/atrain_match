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
from glob import glob
import numpy as np
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)

import matplotlib.pyplot as plt
from utils.get_flag_info import (get_semi_opaque_info_pps2014,
                           get_day_night_twilight_info_pps2014,
                           get_land_coast_sea_info_pps2014,
                           get_mountin_info_pps2014,
                           get_inversion_info_pps2014,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_medium_and_high_clouds_tp,
                           get_calipso_clouds_of_type_i,
                           get_calipso_low_clouds)
from utils.stat_util import (HR_cma, K_cma,
                       PODcy, FARcy,
                       PODcl, FARcl)
from my_dir import ADIR
out_filename = ADIR + "/Documents/A_PPS_v2017/Validation_2018/results_cma_snow.txt"
out_file_h = open(out_filename,'a')

def my_measures(calipso_cfc, pps_cfc, thin, use):
    thin = thin[use]
    calipso_cfc = calipso_cfc[use]
    pps_cfc = pps_cfc[use]

    indict = {}
    indict["det_cloudy"] = np.sum(np.logical_and(pps_cfc==1, calipso_cfc==1))
    indict["det_clear"] = np.sum(np.logical_and(pps_cfc==0, calipso_cfc==0))
    indict["false_cloudy"] = np.sum(np.logical_and(pps_cfc==1, calipso_cfc==0))
    indict["undet_cloudy"] = np.sum(np.logical_and(pps_cfc==0, calipso_cfc==1))
    indict["N"] = len(pps_cfc)

    undet_cloudy_th = np.sum(np.logical_and(pps_cfc==0, np.logical_and(~thin,calipso_cfc==1)))
    det_cloudy_th = np.sum(np.logical_and(pps_cfc==1, np.logical_and(~thin,calipso_cfc==1)))

    measures = {}
    measures["bias"] = np.mean(pps_cfc-calipso_cfc)*100.0
    measures["HR"] = HR_cma(indict)
    measures["K"] = K_cma(indict)
    measures["PODcy"] = PODcy(indict)
    measures["PODcy02"] =  det_cloudy_th*100.0/(det_cloudy_th + undet_cloudy_th)
    measures["Farcy"] =  FARcy(indict)
    measures["PODcl"] = PODcl(indict)
    measures["Farcl"] = FARcl(indict)
    measures["N"] =  indict["N"]
    return measures

def plot_cfc_table(match_calipso,cfc_limit=0.9,sat="modis"):

    from utils.get_flag_info import get_calipso_clouds_of_type_i

    cfc = match_calipso.calipso.all_arrays['cloud_fraction']
    od = match_calipso.calipso.all_arrays['total_optical_depth_5km']
    cma = np.logical_or(match_calipso.imager.all_arrays['cloudmask']==1,
                         match_calipso.imager.all_arrays['cloudmask']==2)
    cl = np.logical_or(match_calipso.imager.all_arrays['cloudmask']==0,
                         match_calipso.imager.all_arrays['cloudmask']==3)
    pps_snowi =  match_calipso.imager.all_arrays['cloudmask']==3
    if match_calipso.imager.all_arrays['cloudmask'] is None:
        cl =np.logical_and(np.less_equal(match_calipso.imager.cloudtype,4),np.greater(match_calipso.imager.cloudtype,0))
        cma = np.logical_and(np.greater(match_calipso.imager.cloudtype,4),np.less(match_calipso.imager.cloudtype,20))
        pps_snowi =np.logical_and(np.less_equal(match_calipso.imager.cloudtype,4),np.greater(match_calipso.imager.cloudtype,2))


    calipso_cfc = cfc.copy()
    calipso_cfc[cfc<=cfc_limit] = 0
    calipso_cfc[cfc>cfc_limit] = 1
    calipso_cfc[cfc<0] = -1
    pps_cfc = 0*cfc.copy()
    pps_cfc[cma] = 1
    use = np.logical_or(cl,cma)



    europe = np.logical_and(match_calipso.imager.all_arrays['longitude']<60,
                            match_calipso.imager.all_arrays['longitude']>-25)
    europe = np.logical_and(europe,match_calipso.imager.all_arrays['latitude']<72)
    europe = np.logical_and(europe,match_calipso.imager.all_arrays['latitude']>35)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag
    ) = get_day_night_twilight_info_pps2014(match_calipso.imager.all_arrays['cloudtype_conditions'])
    non_polar_night = np.logical_or(~night_flag, np.abs(match_calipso.imager.all_arrays['latitude'])<75)
    measures = my_measures(calipso_cfc, pps_cfc, thin, use)
    measures_d = my_measures(calipso_cfc, pps_cfc, thin, np.logical_and(use,day_flag))
    measures_n = my_measures(calipso_cfc, pps_cfc, thin, np.logical_and(use, night_flag))
    measures_t = my_measures(calipso_cfc, pps_cfc, thin, np.logical_and(use,twilight_flag))
    measures_europe = my_measures(calipso_cfc, pps_cfc, thin, np.logical_and(europe,use))
    measures_nonpolar_night =  my_measures(calipso_cfc, pps_cfc, thin, np.logical_and(non_polar_night,use))
    def print_one_line(measures):
        info = ("{:3.1f} {:3.2f} {:3.2f} {:3.1f} {:3.1f} {:3.1f} {:3.1f} {:3.1f} {:d}\n"
                .format(
                    measures["bias"],
                    measures["HR"],
                    measures["K"],
                    measures["PODcy"],
                    measures["PODcy02"],
                    measures["Farcy"],
                    measures["PODcl"],
                    measures["Farcl"],
                    measures["N"]
                ))
        print(info)
        return info


    out_file_h.write("---------------------\n" )
    out_file_h.write(sat + ": \n" )
    info = print_one_line(measures)
    out_file_h.write("CALIOP ALL:" + info)
    info = print_one_line(measures_d)
    out_file_h.write("CALIOP DAY:" + info)
    info = print_one_line(measures_n)
    out_file_h.write("CALIOP NIGHT:" + info)
    info = print_one_line(measures_t)
    out_file_h.write("CALIOP TWILIGHT:" + info)
    info = print_one_line(measures_europe)
    out_file_h.write("CALIOP EUROPE:" + info)
    info = print_one_line(measures_nonpolar_night)
    out_file_h.write("CALIOP NON-POLAR-NIGHT:" + info)

    from trajectory_plotting import plot_satellite_trajectory
    import config
    plot_satellite_trajectory(match_calipso.calipso.all_arrays["longitude"][day_flag],
                            match_calipso.calipso.all_arrays["latitude"][day_flag],
                            ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_cfc_%s_dist"%(sat),
                            config.AREA_CONFIG_FILE,
                            fig_type=['png'])
    try:
        plot_satellite_trajectory(match_calipso.calipso.all_arrays["longitude"][np.logical_and(europe,use)],
                                match_calipso.calipso.all_arrays["latitude"][np.logical_and(europe,use)],
                                ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_cfc_%s_dist_europe"%(sat),
                                config.AREA_CONFIG_FILE,
                                fig_type=['png'])
    except:
        pass
    from histogram_plotting import distribution_map
    distribution_map(match_calipso.calipso.all_arrays["longitude"][use],
                     match_calipso.calipso.all_arrays["latitude"][use])
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_cfc_%s_dist.png"%(sat), bbox_inches='tight')
    try:
        distribution_map(match_calipso.calipso.all_arrays["longitude"][np.logical_and(europe,use)],
                         match_calipso.calipso.all_arrays["latitude"][np.logical_and(europe,use)])
        plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_cfc_%s_dist_europe.png"%(sat), bbox_inches='tight')
    except:
        pass
    plt.close('all')
# ----------------------------------------

if __name__ == "__main__":
    from pps_vrreport_ctth_stats import read_files

    ROOT_DIR_v2014 = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_C4_2014/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_v2018 = (
        ADIR + "/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_C4/Reshaped_Files/npp/1km/2015/07/*/")


    ROOT_DIR_v2014_npp = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_viirs_v2014_created20180914/Reshaped_Files_merged_caliop/npp/1km/2015/*/")
    ROOT_DIR_v2018_npp = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_viirs_v2018_created20180907/Reshaped_Files_merged_caliop/npp/1km/2015/*/")
    ROOT_DIR_v2014_modis = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_modis_v2014_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/")
    ROOT_DIR_v2018_modis = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_modis_v2018_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/")
    ROOT_DIR_v2018_modis_ice = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_modis_v2018_created20181001_cmap_osiice_dust/Reshaped_Files_merged_caliop/eos2/1km/2010/*/")
    ROOT_DIR_v2014_gac = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/*/")
    ROOT_DIR_v2018_gac = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/*/")

    def process_one_case(ROOT_DIR_INNER, exclude_2009=False):
        files = glob(ROOT_DIR_INNER + "*cali*h5")
        if exclude_2009:
            files = [filename for filename in files if "/2009/" not in filename]
        #print files
        sat = ROOT_DIR_INNER.split("global_")[1].split("/Reshaped")[0]
        limit = 0.9
        if "gac" in sat:
            limit = 0.5
        match_calipso = read_files(files)
        plot_cfc_table(match_calipso, limit, sat=sat)
    for ROOT_DIR in [ROOT_DIR_v2014_gac,
                     ROOT_DIR_v2018_gac]:
        process_one_case(ROOT_DIR, True)

    for ROOT_DIR in [ROOT_DIR_v2014_npp,
                     ROOT_DIR_v2018_npp,
                     ROOT_DIR_v2014_modis,
                     ROOT_DIR_v2018_modis,
                     ROOT_DIR_v2018_modis_ice]:
    #,
    #                 ROOT_DIR_v2014_gac,
    #                 ROOT_DIR_v2018_gac]:
        process_one_case(ROOT_DIR)
