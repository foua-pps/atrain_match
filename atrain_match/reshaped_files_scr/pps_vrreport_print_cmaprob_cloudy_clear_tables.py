
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

out_filename = "/home/a001865/Documents/A_PPS_v2017/Validation_2018/results_cma_cmap.txt"
out_file_h = open(out_filename,'a')


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

    calipso_cfc = cfc.copy()
    calipso_cfc[cfc<=cfc_limit] = 0
    calipso_cfc[cfc>cfc_limit] = 1
    calipso_cfc[cfc<0] = -1
    pps_cfc = 0*cfc.copy()
    pps_cfc[cma] = 1
    use = np.logical_or(cl,cma)
    pps_cprob = caObj.imager.all_arrays['cma_prob']
    ct_quality = caObj.imager.all_arrays['cloudtype_quality']
    cma_aerosol = caObj.imager.all_arrays['cma_aerosolflag']

    #import pdb
    #pdb.set_trace()
    pps_cprob[np.isnan(pps_cprob)]=-9
    use = np.logical_and(use, pps_cprob>0)
    cl = np.logical_and(use,cl)
    cma = np.logical_and(use,cma)

 
    info = ("CALIPSO-cloudy	CALIPSO-clear	CALIPSO-cloudy	CALISPO-clear\n" + 
            "Cma-clear	Cma-clear	Cma-cloudy	Cma-Cloudy \n")

    for lower in range(0,99,5):
        upper = lower + 5
        if upper ==100:
            upper = 101
        #ppsclear    
        use_i = np.logical_and(cl,pps_cprob>=lower)
        use_i = np.logical_and(use_i,pps_cprob<upper)
        pclear_cloudy_i = np.logical_and(use_i, calipso_cfc==1) 
        pclear_clear_i = np.logical_and(use_i, calipso_cfc==0)
        use_i = np.logical_and(cma,pps_cprob>=lower)
        use_i = np.logical_and(use_i,pps_cprob<upper)
        pcloudy_cloudy_i = np.logical_and(use_i, calipso_cfc==1) 
        pcloudy_clear_i = np.logical_and(use_i, calipso_cfc==0)
        info += "{:d}<=CMAPROB<{:d}: {:d}-({:d}) {:d}-({:d}) {:d}-({:d}) {:d}-({:d})\n".format(
            lower, 
            upper, 
            np.sum(pclear_cloudy_i),
            np.sum(np.logical_and(ct_quality==24,pclear_cloudy_i)),
            np.sum(pclear_clear_i),
            np.sum(np.logical_and(ct_quality==24,pclear_clear_i)),
            np.sum(pcloudy_cloudy_i),
            np.sum(np.logical_and(ct_quality==24,pcloudy_cloudy_i)),
            np.sum(pcloudy_clear_i),
            np.sum(np.logical_and(ct_quality==24,pcloudy_clear_i)))
    print(info)



    out_file_h.write("---------------------\n" )
    out_file_h.write(sat + ": \n" )
    out_file_h.write(info)

    info = ("CALIPSO-cloudy	CALIPSO-clear	CALIPSO-cloudy	CALISPO-clear\n" + 
            "Cma-clear	Cma-clear	Cma-cloudy	Cma-Cloudy \n")

    for lower in range(0,99,5):
        upper = lower + 5
        if upper ==100:
            upper = 101
        #ppsclear    
        use_i = np.logical_and(cl,pps_cprob>=lower)
        use_i = np.logical_and(use_i,pps_cprob<upper)
        pclear_cloudy_i = np.logical_and(use_i, calipso_cfc==1) 
        pclear_clear_i = np.logical_and(use_i, calipso_cfc==0)
        use_i = np.logical_and(cma,pps_cprob>=lower)
        use_i = np.logical_and(use_i,pps_cprob<upper)
        pcloudy_cloudy_i = np.logical_and(use_i, calipso_cfc==1) 
        pcloudy_clear_i = np.logical_and(use_i, calipso_cfc==0)
        info += "{:d}<=CMAPROB<{:d}: {:d}-({:d}) {:d}-({:d}) {:d}-({:d}) {:d}-({:d})\n".format(
            lower, 
            upper, 
            np.sum(pclear_cloudy_i),
            np.sum(np.logical_and(cma_aerosol==1,pclear_cloudy_i)),
            np.sum(pclear_clear_i),
            np.sum(np.logical_and(cma_aerosol==1,pclear_clear_i)),
            np.sum(pcloudy_cloudy_i),
            np.sum(np.logical_and(cma_aerosol==1,pcloudy_cloudy_i)),
            np.sum(pcloudy_clear_i),
            np.sum(np.logical_and(cma_aerosol==1,pcloudy_clear_i)))
    print(info)
    out_file_h.write("---------------------\n" )
    out_file_h.write(sat + ": \n" )
    out_file_h.write(info)

# ----------------------------------------

if __name__ == "__main__":
    from matchobject_io import read_files

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
    ROOT_DIR_v2018_npp_synop = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_viirs_v2018_created20181010_synop/Reshaped_Files_merged_synop/npp/1km/2015/*/")

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

    #for ROOT_DIR in [ROOT_DIR_v2014_gac,
    #                 ROOT_DIR_v2018_gac]:
    #    process_one_case(ROOT_DIR, True)
                        
    for ROOT_DIR in [#ROOT_DIR_v2014_npp,
                     ROOT_DIR_v2018_npp,
                     #ROOT_DIR_v2014_modis,
                     ROOT_DIR_v2018_modis,
                     ROOT_DIR_v2018_modis_ice]:
                     #ROOT_DIR_v2014_gac,
                     #ROOT_DIR_v2018_gac]:
        process_one_case(ROOT_DIR)
