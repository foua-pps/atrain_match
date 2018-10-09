
"""Read all matched data and make some plotting
"""
import os
from glob import glob
import numpy as np
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)

import matplotlib.pyplot as plt
from get_flag_info import (get_semi_opaque_info_pps2014,
                           get_day_night_twilight_info_pps2014,
                           get_land_coast_sea_info_pps2014,
                           get_mountin_info_pps2014,
                           get_inversion_info_pps2014,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_medium_and_high_clouds_tp,
                           get_calipso_clouds_of_type_i,
                           get_calipso_low_clouds)
from stat_util import (HR_cma, K_cma, 
                       PODcy, FARcy, 
                       PODcl, FARcl)

out_filename = "/home/a001865/Documents/A_PPS_v2017/Validation_2018/results_cma_snow.txt"
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

def plot_cfc_table(caObj,cfc_limit=0.9,sat="modis"):

    from get_flag_info import get_calipso_clouds_of_type_i
    
    cfc = caObj.calipso.all_arrays['cloud_fraction']
    od = caObj.calipso.all_arrays['total_optical_depth_5km']
    cma = np.logical_or(caObj.avhrr.all_arrays['cloudmask']==1,
                         caObj.avhrr.all_arrays['cloudmask']==2)
    cl = np.logical_or(caObj.avhrr.all_arrays['cloudmask']==0,
                         caObj.avhrr.all_arrays['cloudmask']==3)
    pps_snowi =  caObj.avhrr.all_arrays['cloudmask']==3
    if caObj.avhrr.all_arrays['cloudmask'] is None:
        cl =np.logical_and(np.less_equal(caObj.avhrr.cloudtype,4),np.greater(caObj.avhrr.cloudtype,0))
        cma = np.logical_and(np.greater(caObj.avhrr.cloudtype,4),np.less(caObj.avhrr.cloudtype,20))
        pps_snowi =np.logical_and(np.less_equal(caObj.avhrr.cloudtype,4),np.greater(caObj.avhrr.cloudtype,2))

        
    calipso_cfc = cfc.copy()
    calipso_cfc[cfc<=cfc_limit] = 0
    calipso_cfc[cfc>cfc_limit] = 1
    calipso_cfc[cfc<0] = -1
    pps_cfc = 0*cfc.copy()
    pps_cfc[cma] = 1
    use = np.logical_or(cl,cma)

    
    
    europe = np.logical_and(caObj.avhrr.all_arrays['longitude']<60,
                            caObj.avhrr.all_arrays['longitude']>-25)
    europe = np.logical_and(europe,caObj.avhrr.all_arrays['latitude']<72)
    europe = np.logical_and(europe,caObj.avhrr.all_arrays['latitude']>35)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag
    ) = get_day_night_twilight_info_pps2014(caObj.avhrr.all_arrays['cloudtype_conditions'])
    non_polar_night = np.logical_or(~night_flag, np.abs(caObj.avhrr.all_arrays['latitude'])<75)
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

    from trajectory_plotting import plotSatelliteTrajectory
    import config
    plotSatelliteTrajectory(caObj.calipso.all_arrays["longitude"][day_flag],
                            caObj.calipso.all_arrays["latitude"][day_flag],
                            "/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_cfc_%s_dist"%(sat), 
                            config.AREA_CONFIG_FILE,
                            fig_type=['png'])
    try:
        plotSatelliteTrajectory(caObj.calipso.all_arrays["longitude"][np.logical_and(europe,use)],
                                caObj.calipso.all_arrays["latitude"][np.logical_and(europe,use)],
                                "/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_cfc_%s_dist_europe"%(sat), 
                                config.AREA_CONFIG_FILE,
                                fig_type=['png'])
    except:
        pass
    from histogram_plotting import distribution_map
    distribution_map(caObj.calipso.all_arrays["longitude"][use], 
                     caObj.calipso.all_arrays["latitude"][use])
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_cfc_%s_dist.png"%(sat), bbox_inches='tight')
    try:
        distribution_map(caObj.calipso.all_arrays["longitude"][np.logical_and(europe,use)], 
                         caObj.calipso.all_arrays["latitude"][np.logical_and(europe,use)])
        plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_cfc_%s_dist_europe.png"%(sat), bbox_inches='tight')
    except:
        pass
    plt.close('all')  
# ----------------------------------------

if __name__ == "__main__":
    from pps_vrreport_ctth_stats import read_files

    ROOT_DIR_v2014 = (
        "/home/a001865/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_C4_2014/Reshaped_Files/npp/1km/2015/07/*/")
    ROOT_DIR_v2018 = (
        "/home/a001865/DATA_MISC/reshaped_files_jenkins_npp_modis/"
        "ATRAIN_RESULTS_NPP_C4/Reshaped_Files/npp/1km/2015/07/*/")


    ROOT_DIR_v2014_npp = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_viirs_v2014_created20180914/Reshaped_Files_merged_caliop/npp/1km/2015/*/")
    ROOT_DIR_v2018_npp = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_viirs_v2018_created20180907/Reshaped_Files_merged_caliop/npp/1km/2015/*/")
    ROOT_DIR_v2014_modis = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_modis_v2014_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/")
    ROOT_DIR_v2018_modis = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_modis_v2018_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/")
    ROOT_DIR_v2018_modis_ice = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_modis_v2018_created20181001_cmap_osiice_dust/Reshaped_Files_merged_caliop/eos2/1km/2010/*/")
    ROOT_DIR_v2014_gac = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/*/")
    ROOT_DIR_v2018_gac = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/*/")

    def process_one_case(ROOT_DIR_INNER, exclude_2009=False):
        files = glob(ROOT_DIR_INNER + "*cali*h5")
        if exclude_2009:
            files = [filename for filename in files if "/2009/" not in filename]
        #print files    
        sat = ROOT_DIR_INNER.split("global_")[1].split("/Reshaped")[0]
        limit = 0.9
        if "gac" in sat:
            limit = 0.5
        caObj = read_files(files)
        plot_cfc_table(caObj, limit, sat=sat)
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
