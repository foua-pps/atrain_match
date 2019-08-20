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
out_filename = ADIR + "/Documents/A_PPS_v2017/Validation_2018/results_cma_osisaf.txt"
out_file_h = open(out_filename,'w')

def my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, use):
    if use is not None:
        calipso_cfc = calipso_cfc[use]
        pps_cfc = pps_cfc[use]
        pps_snowi = pps_snowi[use]
        calipso_snowi = calipso_snowi[use]

    c_cloudy_w = np.logical_and(calipso_snowi==0, calipso_cfc==1)
    c_cloudy_i = np.logical_and(calipso_snowi>0, calipso_cfc==1)
    c_water = np.logical_and(calipso_snowi==0, calipso_cfc==0)
    c_ice = np.logical_and(calipso_snowi>0, calipso_cfc==0)
    pps_ice = pps_snowi
    pps_water = np.logical_and(pps_cfc==0, pps_snowi==0)
    pps_cloud = pps_cfc==1

    indict = {}
    indict["pps_cloudy_c_cloudy_w"] = np.sum(np.logical_and(pps_cloud, c_cloudy_w))
    indict["pps_cloudy_c_cloudy_i"] = np.sum(np.logical_and(pps_cloud, c_cloudy_i))
    indict["pps_cloudy_c_water"] = np.sum(np.logical_and(pps_cloud, c_water))
    indict["pps_cloudy_c_ice"] = np.sum(np.logical_and(pps_cloud, c_ice))
    indict["pps_water_c_cloudy_w"] = np.sum(np.logical_and(pps_water, c_cloudy_w))
    indict["pps_water_c_cloudy_i"] = np.sum(np.logical_and(pps_water, c_cloudy_i))
    indict["pps_water_c_water"] = np.sum(np.logical_and(pps_water, c_water))
    indict["pps_water_c_ice"] = np.sum(np.logical_and(pps_water, c_ice))
    indict["pps_ice_c_cloudy_w"] = np.sum(np.logical_and(pps_ice, c_cloudy_w))
    indict["pps_ice_c_cloudy_i"] = np.sum(np.logical_and(pps_ice, c_cloudy_i))
    indict["pps_ice_c_water"] = np.sum(np.logical_and(pps_ice, c_water))
    indict["pps_ice_c_ice"] = np.sum(np.logical_and(pps_ice, c_ice))
    indict["N"] = len(pps_cfc)
    indict["pod_cloudy"] = (indict["pps_cloudy_c_cloudy_w"]+indict["pps_cloudy_c_cloudy_i"])*100.0/(np.sum(c_cloudy_w) + np.sum(c_cloudy_i))
    indict["pod_water"] = (indict["pps_water_c_water"])*100.0/(np.sum(c_water))
    indict["pod_ice"] = (indict["pps_ice_c_ice"])*100.0/(np.sum(c_ice))
    indict["far_cloudy"] = (indict["pps_cloudy_c_ice"]+indict["pps_cloudy_c_water"])*100.0/(indict["pps_cloudy_c_cloudy_w"]+indict["pps_cloudy_c_cloudy_i"] + indict["pps_cloudy_c_ice"]+indict["pps_cloudy_c_water"])
    indict["far_water"] = (indict["pps_water_c_cloudy_i"]+indict["pps_water_c_cloudy_w"])*100.0/(indict["pps_water_c_cloudy_w"]+indict["pps_water_c_cloudy_i"] + indict["pps_water_c_ice"]+indict["pps_water_c_water"])
    indict["far_ice"] = (indict["pps_ice_c_cloudy_i"]+indict["pps_ice_c_cloudy_w"])*100.0/(indict["pps_ice_c_cloudy_w"]+indict["pps_ice_c_cloudy_i"] + indict["pps_ice_c_ice"]+indict["pps_ice_c_water"])





    return indict


def plot_cfc_table(caObj,cfc_limit=0.9,sat="modis"):

    from utils.get_flag_info import get_calipso_clouds_of_type_i
    #cal_subset = np.logical_and(np.equal(nsidc_st,0),np.equal(igbp_st,17))
    cfc = caObj.calipso.all_arrays['cloud_fraction']
    calipso_snowi = caObj.calipso.all_arrays['nsidc_surface_type']
    od = caObj.calipso.all_arrays['total_optical_depth_5km']
    cma = np.logical_or(caObj.imager.all_arrays['cloudmask']==1,
                         caObj.imager.all_arrays['cloudmask']==2)
    cl = np.logical_or(caObj.imager.all_arrays['cloudmask']==0,
                         caObj.imager.all_arrays['cloudmask']==3)
    pps_snowi =  caObj.imager.all_arrays['cloudmask']==3
    if caObj.imager.all_arrays['cloudmask'] is None:
        cl =np.logical_and(np.less_equal(caObj.imager.cloudtype,4),np.greater(caObj.imager.cloudtype,0))
        cma = np.logical_and(np.greater(caObj.imager.cloudtype,4),np.less(caObj.imager.cloudtype,20))
        pps_snowi =np.logical_and(np.less_equal(caObj.imager.cloudtype,4),np.greater(caObj.imager.cloudtype,2))

    (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag) = get_land_coast_sea_info_pps2014(caObj.imager.cloudtype_conditions)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag) = get_day_night_twilight_info_pps2014(caObj.imager.cloudtype_conditions)
    calipso_cfc = cfc.copy()
    calipso_cfc[cfc<=cfc_limit] = 0
    calipso_cfc[cfc>cfc_limit] = 1
    calipso_cfc[cfc<0] = -1
    pps_cfc = 0*cfc.copy()
    pps_cfc[cma] = 1
    use = np.logical_or(cl,cma)
    use = np.logical_and(use, sea_flag)
    thin = None
    measures = my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, 
                           use)
    measures_d = my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, 
                              np.logical_and(use, day_flag))
    measures_t = my_measures(calipso_cfc, pps_cfc,  calipso_snowi, pps_snowi, thin, 
                              np.logical_and(use,twilight_flag))
    measures_n = my_measures(calipso_cfc, pps_cfc,  calipso_snowi, pps_snowi, thin, 
                              np.logical_and(use,night_flag))
    def print_one_line(measures):
        info = "CALIOP-CLOUD(W) CALIOP-CLOUD(I) CALIOP-WATER CALIOP-ICE\n" 
        info += ("*********\nPPS-CLOUD: {:3.1f} {:3.1f} {:3.1f} {:3.1f} \nPPS-WATER {:3.1f} {:3.1f} {:3.1f} {:3.1f} \nPPS-ICE {:3.1f} {:3.1f} {:3.1f} {:3.1f} \n"
                .format(
                    measures["pps_cloudy_c_cloudy_w"] ,
                    measures["pps_cloudy_c_cloudy_i"] ,
                    measures["pps_cloudy_c_water"] ,
                    measures["pps_cloudy_c_ice"] ,
                    measures["pps_water_c_cloudy_w"] ,
                    measures["pps_water_c_cloudy_i"] ,
                    measures["pps_water_c_water"] ,
                    measures["pps_water_c_ice"] ,
                    measures["pps_ice_c_cloudy_w"], 
                    measures["pps_ice_c_cloudy_i"] ,
                    measures["pps_ice_c_water"] ,
                    measures["pps_ice_c_ice"] 
                ))
        info += " POD-C:{:3.1f} POD-W:{:3.1f} POD-I:{:3.1f} \n FAR-C:{:3.1f} FAR-W(cloud):{:3.1f} FAR-I(cloud):{:3.1f} \n".format(
            measures["pod_cloudy"], measures["pod_water"], measures["pod_ice"], measures["far_cloudy"], measures["far_water"],measures["far_ice"],)
        print(info)
        return info



    out_file_h.write("---------------------\n" )
    out_file_h.write(sat + ": \n" )
    info = print_one_line(measures)
    out_file_h.write("CALIOP-ALL: " + info)
    info = print_one_line(measures_d)
    out_file_h.write("CALIOP-DAY: " + info)
    info = print_one_line(measures_n)
    out_file_h.write("CALIOP-NIGHT: " + info)
    info = print_one_line(measures_t)
    out_file_h.write("CALIOP-TWILIGHT: " + info)




# ----------------------------------------

if __name__ == "__main__":
    from matchobject_io import read_files

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
    ROOT_DIR_v2018_npp_synop = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_viirs_v2018_created20181010_synop/Reshaped_Files_merged_synop/npp/1km/2015/*/")

    def process_one_case(ROOT_DIR_INNER, exclude_2009=False):
        files = glob(ROOT_DIR_INNER + "*cali*h5")
        if exclude_2009:
            files = [filename for filename in files if "/2009/" not in filename]
        #print files    
        sat = ROOT_DIR_INNER.split("global_")[1].split("/Reshaped")[0]
        limit = 0.9
        if "gac" in sat:
            limit = 0.5
        caObj = read_files(files, 'calipso')
        plot_cfc_table(caObj, limit, sat=sat)


    for ROOT_DIR in [ROOT_DIR_v2014_gac,
                     ROOT_DIR_v2018_gac]:
        process_one_case(ROOT_DIR, True)
                        
    for ROOT_DIR in [ROOT_DIR_v2014_npp,
                     ROOT_DIR_v2018_npp,
                     ROOT_DIR_v2014_modis,
                     ROOT_DIR_v2018_modis,
                     ROOT_DIR_v2018_modis_ice]:
                     #ROOT_DIR_v2014_gac,
                     #ROOT_DIR_v2018_gac]:
        process_one_case(ROOT_DIR)
