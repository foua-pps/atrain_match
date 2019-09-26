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
out_filename = ADIR + "/Documents/A_PPS_v2017/Validation_2018/results_cma.txt"
out_file_h = open(out_filename,'a')

def my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, use):
    if use is not None:
        calipso_cfc = calipso_cfc[use]
        pps_cfc = pps_cfc[use]


    indict = {}
    indict["det_cloudy"] = np.sum(np.logical_and(pps_cfc==1, calipso_cfc==1))
    indict["det_clear"] = np.sum(np.logical_and(pps_cfc==0, calipso_cfc==0))
    indict["false_cloudy"] = np.sum(np.logical_and(pps_cfc==1, calipso_cfc==0))
    indict["undet_cloudy"] = np.sum(np.logical_and(pps_cfc==0, calipso_cfc==1))
    indict["N"] = len(pps_cfc)
    measures = {}
    measures["bias"] = np.mean(pps_cfc-calipso_cfc)*100.0
    measures["HR"] = HR_cma(indict)
    measures["K"] = K_cma(indict)
    measures["PODcy"] = PODcy(indict)
    measures["Farcy"] =  FARcy(indict)
    measures["PODcl"] = PODcl(indict)
    measures["Farcl"] = FARcl(indict)
    measures["N"] =  indict["N"]

    if calipso_snowi is not None:
        if use is not None:
            pps_snowi = pps_snowi[use]
            calipso_snowi =  calipso_snowi[use]
        indict["det_snow"] = np.sum(np.logical_and(pps_snowi, np.logical_and(calipso_snowi>0, calipso_cfc==0)))
        indict["det_snow_much"] = np.sum(np.logical_and(pps_snowi, np.logical_and(calipso_snowi>=10, calipso_cfc==0)))
        indict["undet_snow_much"] = np.sum(np.logical_and(~pps_snowi, np.logical_and(calipso_snowi>=10, calipso_cfc==0)))
        indict["undet_snow"] = np.sum(np.logical_and(~pps_snowi, np.logical_and(calipso_snowi>0, calipso_cfc==0)))
        indict["false_snow"] = np.sum(np.logical_and(pps_snowi,  np.logical_or(calipso_snowi==0, calipso_cfc==1)))
        indict["false_snow_more_than_cloud_issue"] = np.sum(np.logical_and(pps_snowi, ~calipso_snowi==0 ))
        indict["N_pps_snow"] = np.sum(pps_snowi)

        measures["PODsnow"] =  indict["det_snow"]*100.0/(indict["det_snow"] + indict["undet_snow"])
        measures["PODsnow_much"] =  indict["det_snow_much"]*100.0/(indict["det_snow_much"] + indict["undet_snow_much"])
        measures["FARsnow"] =  indict["false_snow"]*100.0/(indict["N_pps_snow"])
        measures["FARsnow_bad"] =  indict["false_snow_more_than_cloud_issue"]*100.0/(indict["N_pps_snow"])
    if thin is not None:
        if use is not None:
            thin = thin[use]
        undet_cloudy_th = np.sum(np.logical_and(pps_cfc==0, np.logical_and(~thin,calipso_cfc==1)))
        det_cloudy_th = np.sum(np.logical_and(pps_cfc==1, np.logical_and(~thin,calipso_cfc==1)))
        measures["PODcy02"] =  det_cloudy_th*100.0/(det_cloudy_th + undet_cloudy_th)


    return measures


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

    thin = np.logical_and(od>0,
                          od<0.2)
    """
    my_fcf_values = set(cfc)

    for limit in sorted(my_fcf_values):
        use = cfc==limit
        use_thick = np.logical_and(cfc==limit,np.logical_or(od<0, od>0.225))
        ok_cma = np.logical_and(use, cma)
        ok_cl = np.logical_and(use, cl)
        ok_cma_thick = np.logical_and(use_thick, cma)
        #cloudy = np.logcial_and(
        print "%3.1f"%(limit) , len(cfc[use]),
        print "(%d)"%(len(cfc[use_thick]))
        print "cloudy" , "%3.1f"%(len(cfc[ok_cma])*100.0/len(cfc[use])),
        print "(%3.1f)"%(len(cfc[ok_cma_thick])*100.0/np.max([1,len(cfc[use_thick])]))
        print "cear" , "%3.1f"%(len(cfc[ok_cl])*100.0/len(cfc[use])),
        print ""
    """

    calipso_cfc = cfc.copy()
    calipso_cfc[cfc<=cfc_limit] = 0
    calipso_cfc[cfc>cfc_limit] = 1
    calipso_cfc[cfc<0] = -1
    pps_cfc = 0*cfc.copy()
    pps_cfc[cma] = 1
    use = np.logical_or(cl,cma)


    europe = np.logical_and(caObj.imager.all_arrays['longitude']<60,
                            caObj.imager.all_arrays['longitude']>-25)
    europe = np.logical_and(europe,caObj.imager.all_arrays['latitude']<72)
    europe = np.logical_and(europe,caObj.imager.all_arrays['latitude']>35)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag
    ) = get_day_night_twilight_info_pps2014(caObj.imager.all_arrays['cloudtype_conditions'])
    (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag
    ) = get_land_coast_sea_info_pps2014(caObj.imager.all_arrays['cloudtype_conditions'])
    non_polar_night = np.logical_or(~night_flag, np.abs(caObj.imager.all_arrays['latitude'])<75)
    measures = my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi, thin, use)
    measures_d = my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, np.logical_and(use,day_flag))
    measures_n = my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, np.logical_and(use, night_flag))
    measures_t = my_measures(calipso_cfc, pps_cfc,  calipso_snowi, pps_snowi, thin, np.logical_and(use,twilight_flag))
    measures_ds = my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, np.logical_and(sea_flag,np.logical_and(use, day_flag)))
    measures_ts = my_measures(calipso_cfc, pps_cfc,  calipso_snowi, pps_snowi, thin, np.logical_and(sea_flag,np.logical_and(use,twilight_flag)))
    measures_europe = my_measures(calipso_cfc, pps_cfc, calipso_snowi, pps_snowi,  thin, np.logical_and(europe,use))
    measures_nonpolar_night =  my_measures(calipso_cfc, pps_cfc,  calipso_snowi, pps_snowi, thin, np.logical_and(non_polar_night,use))
    def print_one_line(measures):
        info = ("{:3.1f} {:3.2f} {:3.2f} {:3.1f} {:3.1f} {:3.1f} {:3.1f} {:3.1f} {:d} {:3.1f} {:3.1f}  {:3.1f} {:3.1f} \n"
                .format(
                    measures["bias"],
                    measures["HR"],
                    measures["K"],
                    measures["PODcy"],
                    measures["PODcy02"],
                    measures["Farcy"],
                    measures["PODcl"],
                    measures["Farcl"],
                    measures["N"],
                    measures["PODsnow"],
                    measures["PODsnow_much"],
                    measures["FARsnow"],
                    measures["FARsnow_bad"],

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

    info = print_one_line(measures_ds)
    out_file_h.write("CALIOP DAY SEA:" + info)
    info = print_one_line(measures_ts)
    out_file_h.write("CALIOP TWILIGHT SEA:" + info)
    info = print_one_line(measures_europe)
    out_file_h.write("CALIOP EUROPE:" + info)
    info = print_one_line(measures_nonpolar_night)
    out_file_h.write("CALIOP NON-POLAR-NIGHT:" + info)

    from trajectory_plotting import plotSatelliteTrajectory
    import config
    plotSatelliteTrajectory(caObj.calipso.all_arrays["longitude"][use],
                            caObj.calipso.all_arrays["latitude"][use],
                            ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_cfc_%s_dist"%(sat),
                            config.AREA_CONFIG_FILE,
                            fig_type=['png'])
    try:
        plotSatelliteTrajectory(caObj.calipso.all_arrays["longitude"][np.logical_and(europe,use)],
                                caObj.calipso.all_arrays["latitude"][np.logical_and(europe,use)],
                                ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_cfc_%s_dist_europe"%(sat),
                                config.AREA_CONFIG_FILE,
                                fig_type=['png'])
    except:
        pass
    from histogram_plotting import distribution_map
    distribution_map(caObj.calipso.all_arrays["longitude"][use],
                     caObj.calipso.all_arrays["latitude"][use])
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_cfc_%s_dist.png"%(sat), bbox_inches='tight')
    try:
        distribution_map(caObj.calipso.all_arrays["longitude"][np.logical_and(europe,use)],
                         caObj.calipso.all_arrays["latitude"][np.logical_and(europe,use)])
        plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_cfc_%s_dist_europe.png"%(sat), bbox_inches='tight')
    except:
        pass
    plt.close('all')

def val_synop_cfc(caObj, sat):
    europe = np.logical_and(caObj.imager.all_arrays['longitude']<60,
                            caObj.imager.all_arrays['longitude']>-25)
    europe = np.logical_and(europe,caObj.imager.all_arrays['latitude']<72)
    europe = np.logical_and(europe,caObj.imager.all_arrays['latitude']>35)
    use_pps = np.logical_or(caObj.imager.cfc_mean<0.25, caObj.imager.cfc_mean>=0.75)
    use_synop = np.logical_or(caObj.synop.cloud_fraction<0.25, caObj.synop.cloud_fraction>=0.75)
    use_all = np.logical_and(use_pps, use_synop)
    use_all_europe = np.logical_and(use_all, europe)
    synop_cfc = 0*caObj.synop.cloud_fraction.copy()
    synop_cfc[caObj.synop.cloud_fraction>=0.75] = 1.0
    pps_cfc = 0*caObj.imager.cfc_mean.copy()
    pps_cfc[caObj.imager.cfc_mean>=0.75] = 1.0
    def print_one_line(measures):
        info = ("{:3.1f} {:3.2f} {:3.1f} {:3.1f} {:3.1f} {:3.1f} {:d} \n"
                .format(
                    measures["bias"],
                    measures["HR"],
                    #measures["K"],
                    measures["PODcy"],
                    #measures["PODcy02"],
                    measures["Farcy"],
                    measures["PODcl"],
                    measures["Farcl"],
                    measures["N"],
                    #measures["PODsnow"],
                    #measures["PODsnow_much"],
                    #measures["FARsnow"],
                    #measures["FARsnow_bad"],

                ))
        print(info)
        return info
    measures = my_measures(synop_cfc, pps_cfc, None, None,  None, use_all)
    measures_europe = my_measures(synop_cfc, pps_cfc, None, None,  None, use_all_europe)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag
    ) = get_day_night_twilight_info_pps2014(caObj.imager.all_arrays['cloudtype_conditions'])
    measures_d = my_measures(synop_cfc, pps_cfc, None, None,  None, np.logical_and(use_all,day_flag))
    measures_n = my_measures(synop_cfc, pps_cfc, None, None,  None, np.logical_and(use_all, night_flag))
    measures_t = my_measures(synop_cfc, pps_cfc, None, None,  None, np.logical_and(use_all,twilight_flag))

    out_file_h.write("---------------------\n" )
    out_file_h.write(sat + ": \n" )
    info = print_one_line(measures)
    out_file_h.write("SYNOP ALL:" + info)
    info = print_one_line(measures_d)
    out_file_h.write("SYNOP DAY:" + info)
    info = print_one_line(measures_n)
    out_file_h.write("SYNOP NIGHT:" + info)
    info = print_one_line(measures_t)
    out_file_h.write("SYNOP TWILIGHT:" + info)

    info = print_one_line(measures_europe)
    out_file_h.write("SYNOP EUROPE:" + info)


def plot_synop_cfc(caObj):
    print caObj.synop.all_arrays["longitude"]
    from histogram_plotting import distribution_map
    distribution_map(caObj.synop.all_arrays["longitude"],
                     caObj.synop.all_arrays["latitude"])
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_cfc_synop_%s_dist.png"%(sat), bbox_inches='tight')

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

    for ROOT_DIR in [ROOT_DIR_v2018_npp_synop]:
        files = glob(ROOT_DIR + "*syno*h5")
        sat = ROOT_DIR.split("global_")[1].split("/Reshaped")[0]
        caObj = read_files(files,'synop')
        plot_synop_cfc(caObj)
        val_synop_cfc(caObj, sat=sat)

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
