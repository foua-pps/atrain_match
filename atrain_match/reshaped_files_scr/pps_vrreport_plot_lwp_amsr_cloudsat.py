import numpy as np
from glob import glob
import os
from matchobject_io import (readAmsrImagerMatchObj,
                            AmsrImagerTrackObject,
                            readCloudsatImagerMatchObj,
                            CloudsatImagerTrackObject)
import matplotlib.pyplot as plt 
from utils.validate_lwp_util import get_lwp_diff_inner
from libs.truth_imager_statistics import get_lwp_diff_inner_cloudsat

from histogram_plotting import atrain_scatter
from utils.stat_util import (my_hist, my_iqr, my_rms, my_pex, my_mae)

def my_label(data):
    label = (#"{:s}\n"
             "bias = {:3.1f}\n"
             "RMS = {:3.1f}\n"
             "median = {:3.1f}\n"
             "MAE = {:3.1f}\n"
             "IQR = {:3.1f}\n"
             "PE>10 = {:3.1f}\n"
             "N = {:d}\n".format(
                 #text
                 np.mean(data),
                 my_rms(data),
                 np.median(data),
                 my_mae(data),
                 my_iqr(data),
                 my_pex(data,10),
                 len(data)
                            ))
    return label

def do_the_printing(aObj, name):
    conditions = {}
    if "amsr" in  aObj.truth_sat and len(aObj.imager.cpp_lwp.shape)==2:
        xymax=170
        diff, x, y, use = get_lwp_diff_inner(aObj, True)
        lon = aObj.amsr.all_arrays["longitude"]
        lat = aObj.amsr.all_arrays["latitude"]
        use_map = use
        conditions["use"] = use.copy()
        use_surfs = [use]
        surf_names = ["sea"]
        scatter_bin_size=5.0
    elif "amsr" in  aObj.truth_sat:    
        y = aObj.amsr.lwp
        x = aObj.imager.cpp_lwp
        lon = aObj.imager.all_arrays["longitude"]
        lat = aObj.imager.all_arrays["latitude"]
        xymax=170
        diff = x-y
        use = np.ones(diff.shape, dtype=bool)
        use_map = use
        conditions["use"] = use.copy()
        use_surfs = [use]
        surf_names = ["sea"]
        scatter_bin_size=5.0
    elif "cloudsat"  in aObj.truth_sat:
        xymax=300
        diff, lwp_diff_lo, x, y,  use = get_lwp_diff_inner_cloudsat(aObj, True, wide_selection=True)
        lon = aObj.cloudsat.all_arrays["longitude"]
        lat = aObj.cloudsat.all_arrays["latitude"]
        conditions["use"] = use.copy()
        conditions["no_cloudsat_iwp"] =  aObj.cloudsat.RVOD_ice_water_path<=0
        conditions["no_cloudsat_precip"] =  np.bitwise_and(np.right_shift(aObj.cloudsat.RVOD_CWC_status,2),1)==0
        conditions["geoprof_cloudy"] = np.bitwise_and(np.right_shift(aObj.cloudsat.RVOD_CWC_status,0),10)==0
        conditions["without_zeros"] = np.logical_and(aObj.imager.cpp_lwp>0,
                                                    aObj.cloudsat.RVOD_liq_water_path>0)

        use_sea = np.logical_and(use,aObj.imager.fractionofland <=0)
        use_land = np.logical_and(use,aObj.imager.fractionofland >0)
        use_surfs = [use, use_sea, use_land]
        surf_names = ["all", "sea", "land"]
        scatter_bin_size=5.0
    for use_t, surf in zip(use_surfs,surf_names):
        use_i = use_t
        for condition in ["use", "no_cloudsat_iwp", "no_cloudsat_precip", "geoprof_cloudy", "without_zeros"]:
            if condition not in conditions.keys():
                continue
            use_i = np.logical_and(use_i, conditions[condition])
            if surf == "all" and condition == "geoprof_cloudy":
                use_map = use_i
            #print x, y
            diff_i = x[use_i]- y[use_i]
            print name, surf, condition, np.mean(diff_i), my_rms(diff_i), my_iqr(diff_i), np.median(diff_i), np.std(diff_i), len(diff_i)
            
            vmax = len(diff_i)*0.002
            fig = plt.figure()
            ax = fig.add_subplot(111)
            atrain_scatter(fig, ax, x[use_i] , y[use_i],   scatter_bin_size, xymin = 0,  xymax=1000, vmax=vmax, 
                           do_colorbar=True, to_km=1.0, ptype='scatter')
            ax.plot([0,170],[0,170], ':w')
            ax.set_xlabel("PPS LWP g/m^2")
            if "amsr" in  aObj.truth_sat:
                ax.set_ylabel("AMSR-E LWP g/m^2")
            else:
                ax.set_ylabel("CPR (CloudSat) LWP g/m^2")
            plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/scatter_lwp_%s_%s_dist_all_cpp_%s_%s.png"%(aObj.truth_sat, name, surf, condition), bbox_inches='tight')

            fig = plt.figure(figsize=(14, 4.5))
            ax = fig.add_subplot(131)
            hist_heights, x_, hist_heights_gaussian = my_hist(x[use_i]-y[use_i], None, bmin=-3000, bmax=3000, delta_h=5, return_also_corresponding_gaussian=True)

            ax.fill(x_, hist_heights_gaussian, color='silver')
            plt.plot(x_, hist_heights, "r-",  label = my_label(x[use_i]-y[use_i]))
            #ax.set_ylim(0,10)    
            #ax.set_xlim(-4,6)
            ax.set_xlim([-500,600])
            ax.set_ylim([00,14])
            if "amsr" in  aObj.truth_sat:
                ax.set_xlabel("error: PPS - AMSR-E LWP g/m^2")
            else:    
                ax.set_xlabel("error: PPS - CPR (CloudSat) LWP g/m^2")
            ax.set_ylabel("percent of data")
                #plt.yticks(np.arange(0,9,2.0))
            plt.legend()

            ax = fig.add_subplot(132)
            b = atrain_scatter(fig, ax, x[use_i] , y[use_i],   scatter_bin_size, xymin = 0,  xymax=xymax, vmax=vmax, #140 
                           do_colorbar=False, to_km=1.0, ptype='scatter')
            ax.set_xlabel("PPS LWP g/m^2")
            ax.plot([0,170],[0,170], ':w')
            cax = fig.add_axes([0.57, 0.65, 0.03, 0.22])
            cbar = fig.colorbar(b, cax=cax, )

            if "amsr" in  aObj.truth_sat:
                ax.set_ylabel("AMSR-E LWP g/m^2")

            else:
                ax.set_ylabel("CPR (CloudSat) LWP g/m^2")
                
            ax = fig.add_subplot(133)
            hist_heights_x, x_, hist_heights_gaussian = my_hist(x[use_i], None, bmin=0, bmax=300, delta_h=5)
            hist_heights_y, y_, hist_heights_gaussian = my_hist(y[use_i], None, bmin=0, bmax=300, delta_h=5)
            plt.plot(x_, hist_heights_x, "r-",  label = "PPS LWP")
            my_label_s = "CPR (CloudSat) LWP"
            if "amsr" in  aObj.truth_sat:
                my_label_s=("AMSR-E LWP")
            plt.plot(y_, hist_heights_y, "b-",  label = my_label_s)
            plt.legend()
            ax.set_ylim([0,14])
            ax.set_xlabel("LWP g/m^2")
            ax.set_ylabel("percent of data")

            #plt.show()
            plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/scatter_and_error_lwp_%s_%s_dist_all_cpp_170_%s_%s.png"%(aObj.truth_sat, name, surf, condition), bbox_inches='tight')
            #plt.show()
            plt.close('all')
   
    from trajectory_plotting import plotSatelliteTrajectory
    import config
    #plotSatelliteTrajectory(lon, lat,
    #                        "/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_lwp_%s_dist_all"%(aObj.truth_sat), 
    #                        config.AREA_CONFIG_FILE,
    #                        fig_type=['png'])
    plotSatelliteTrajectory(lon[use_map], lat[use_map],
                            "/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_marble_lwp_%s_%s_dist"%(aObj.truth_sat, name), 
                            config.AREA_CONFIG_FILE,
                            fig_type=['png'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    from histogram_plotting import distribution_map
    distribution_map(lon[use_map], 
                     lat[use_map])
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/map_white_lwp_%s_%s_dist_all_cpp.png"%(aObj.truth_sat, name), bbox_inches='tight')

if __name__ == "__main__":
    ROOT_DIR_v2014_amsr = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2014_created20180927/Reshaped_Files_merged_amsr_lwp/noaa18/5km/")
    ROOT_DIR_v2018_amsr = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2018_created20180927/Reshaped_Files_merged_amsr_lwp/noaa18/5km/")
    #ROOT_DIR_v2014_cl = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2014_created20180927/Reshaped_Files_merged/noaa18/5km/*/*cloudsat")
    #ROOT_DIR_v2018_cl = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2018_created20180927/Reshaped_Files_merged/noaa18/5km/*/*cloudsat")

    ROOT_DIR_v2018_modis_cl = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_modis_v2018_created20181001_cmap_osiice_dust/Reshaped_Files_merged_cloudsat_lwp/eos2/1km/*/*/")
    ROOT_DIR_v2018_modis_amsr = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_modis_v2018_created20181001_cmap_osiice_dust/Reshaped_Files_merged_amsr_lwp/eos2/1km/*/*/")
    ROOT_DIR_v2014_modis_amsr = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_modis_v2014_created20180920/Reshaped_Files_merged_amsr_lwp/eos2/1km/*/*/")

    files = glob(ROOT_DIR_v2018_modis_amsr + "*h5")
    print files
    name = files[0].split('global_')[1][0:11]
    
    if "amsr" in  os.path.basename(files[0]):
        aObj = AmsrImagerTrackObject()
        for filename in files:
            print  os.path.basename(filename)
            aObj += readAmsrImagerMatchObj(filename)
    else: 
        aObj = CloudsatImagerTrackObject()
        for filename in files:
            print  os.path.basename(filename)
            if os.path.basename(filename) in [
                    "5km_noaa18_20071202_1600_99999_cloudsat_imager_match.h5"
            ]:

                continue

            aObj_new = readCloudsatImagerMatchObj(filename) 
            if aObj_new.cloudsat.RVOD_liq_water_path is None:
                print  os.path.basename(filename)
            if len(aObj_new.cloudsat.RVOD_liq_water_path) == len(aObj_new.imager.cpp_lwp):
                print "ok",  os.path.basename(filename)
                aObj += aObj_new     
    do_the_printing(aObj, name)
