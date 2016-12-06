"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       ppsMatch_Imager_CalipsoObject)
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
from get_flag_info import get_calipso_clouds_of_type_i
from get_flag_info import (get_semi_opaque_info_pps2014,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_low_clouds)

def make_boxplot(caObj, name):
    if 'old' in name:
        name = "pps_traditional"
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    height_c = (1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] - 
                caObj.calipso.all_arrays['elevation'])
    height_pps = caObj.avhrr.all_arrays['ctth_height']
    use = np.logical_and(height_pps >-1,
                         caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-1)

    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)
    bias = height_pps - height_c
    limit = np.percentile(bias[use],5)
    print limit 
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    fig = plt.figure(figsize = (6,9))        
    ax = fig.add_subplot(111)
    #plt.tight_layout()
    #plt.subplots_adjust(left=0.2)
    #plt.subplots_adjust(left=10, bottom=10, right=10, top=10, wspace=0, hspace=0)
    plt.boxplot([bias[low],bias[medium],bias[high]],whis=[2, 98],sym='',labels=["low","medium","high"],showmeans=True)
    ax.set_ylim(-14000,8000)
    ax.fill_between(np.arange(0,6),-500,500, facecolor='green', alpha=0.6)
    ax.fill_between(np.arange(0,6),-1000,1000, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(0,6),-1500,1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0,6),2000,15000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0,6),-2000,-15000, facecolor='red', alpha=0.2)
    for y_val in [-5,-4,-3,-2,2,3,4,5]:
        plt.plot(np.arange(0,6), y_val*1000 + 0*np.arange(0,6),':k')
        plt.plot(np.arange(0,6), -10*1000 + 0*np.arange(0,6),':k')
    plt.plot(np.arange(0,6), 0 + 0*np.arange(0,6),'k')
    plt.title(name)
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_LAPSE_RATE_INVESTIGATION/ctth_box_plot_%s.png"%(name))
    #plt.show()


def make_compare(caObj, caObj2, name):
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    height_c = (1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] - 
                caObj.calipso.all_arrays['elevation'])
    height_pps1 = caObj.avhrr.all_arrays['ctth_height']
    height_pps2 = caObj2.avhrr.all_arrays['ctth_height']
    use = np.logical_and(height_pps1 >-1,height_pps2>-1)
    use = np.logical_and(use,
                         caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-1)

    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)
    bias1 = height_pps1 - height_c
    bias2 = height_pps2 - height_c
    limit = np.percentile(bias2[use],5)
    print limit 
    fig = plt.figure(figsize = (9,9))        
    ax = fig.add_subplot(111)
    plt.boxplot([bias1[low],bias1[medium],bias1[high], bias2[low],
                 bias2[medium],bias2[high]],whis=[2, 98],sym='',
                labels=["low1","medium1","high1",
                        "low2","medium2","high2"],showmeans=True)
    ax.set_ylim(-15000,15000)
    ax.fill_between(np.arange(0,12),-1500,1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0,12),2000,15000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0,12),-2000,-15000, facecolor='red', alpha=0.2)
    for y_val in xrange(-5000,6000,1000):
        plt.plot(np.arange(0,12), y_val + 0*np.arange(0,12),':k')
    plt.title(name)
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_LAPSE_RATE_INVESTIGATION/ctth_box_plot_%s.png"%(name))
    plt.show()

    """
    fig = plt.figure(figsize = (9,9))        
    ax = fig.add_subplot(111)
    plt.plot(np.sorted(height_pps1 - height_pps2))
    plt.title(name)
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_LAPSE_RATE_INVESTIGATION/ctth_compare_%s.png"%(name))
    plt.show()
    """
def investigate_nn_ctth():
    ROOT_DIR_GAC_nn = ("/home/a001865/DATA_MISC/reshaped_files/"
                       "ATRAIN_RESULTS_GAC_nn21/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_old = ("/home/a001865/DATA_MISC/reshaped_files/"
                        "ATRAIN_RESULTS_GAC_old/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn_new = ("/home/a001865/DATA_MISC/reshaped_files/"
                           "ATRAIN_RESULTS_GAC_nn20161125/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_nn_18 = ("/home/a001865/DATA_MISC/reshaped_files/"
                           "ATRAIN_RESULTS_GAC_nn20161130/Reshaped_Files/noaa18/")

    re_name = re.compile("_RESULTS_GAC_(\w+)\/")
    caobj_dict = {}
    for ROOT_DIR in [ROOT_DIR_GAC_nn, ROOT_DIR_GAC_old, 
                     ROOT_DIR_GAC_nn_new, ROOT_DIR_GAC_nn_18]:
        match = re_name.search(ROOT_DIR)
        name = "no_name"
        if match:
            name = match.group(1)
        files = glob(ROOT_DIR + "5km/2009/*/*/*h5")
        caObj = CalipsoAvhrrTrackObject()
        for filename in files:
            caObj +=  readCaliopAvhrrMatchObj(filename) 
        caobj_dict[name] = caObj    
        make_boxplot(caObj, name) 
    make_compare(caobj_dict['old'],caobj_dict['nn20161125'],'test')
    make_compare(caobj_dict['nn20161130'],caobj_dict['nn20161125'],'test2')


def investigate_nn_ctth_modis():
    #november
    ROOT_DIR_MODIS_nn = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_viirs = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_NOVEMBER_nnviirs_20161205/Reshaped_Files/merged/")
    #may
    
    ROOT_DIR_MODIS_nn = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_MAY/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_viirs = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_MAY_nnviirs_20161205/Reshaped_Files/merged/")
    ROOT_DIR_MODIS_nn_mersi2 = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "ATRAIN_RESULTS_MODIS_MAY_nnmersi2_20161206/Reshaped_Files/merged/")
   
    ROOT_DIR_MODIS_old = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "global_modis_14th_created20161108/Reshaped_Files/merged/")

    caobj_dict = {}
    for ROOT_DIR, name in zip(
            [ROOT_DIR_MODIS_nn, ROOT_DIR_MODIS_nn_viirs, 
             ROOT_DIR_MODIS_nn_mersi2, ROOT_DIR_MODIS_old], 
            ["modis_may_nnAVHRR","modis_may_nnVIIRS","modis_may_nnMERSI2",
             "modis_may_traditional_pps"]):
        print name
        files = glob(ROOT_DIR + "/*05*.h5")
        caObj = CalipsoAvhrrTrackObject()
        for filename in files:
            print filename
            caObj +=  readCaliopAvhrrMatchObj(filename) 
        caobj_dict[name] = caObj    
        make_boxplot(caObj, name) 
    #make_compare(caobj_dict["modis_nn18var"],
    #             caobj_dict["modis_traditional_pps"],
    #             'compare_modis')

if __name__ == "__main__":
    #investigate_nn_ctth()
    investigate_nn_ctth_modis()   
    """
    isModis1km = False
    isNPP_v2014 = False
    isGAC_v2014_morning_sat = False
    isGAC_v2014 = True
    nn = True
    if nn:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_GAC/Reshaped_Files/noaa18/5km/2009/*/cea5km_test/"
        files = glob(ROOT_DIR + "*.h5") 
    elif isModis1km:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20161108/"
        #files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/02/2010??14_07*/*h5")
        files = glob(ROOT_DIR + "Reshaped_Files/merged/*night*01*.h5")
    elif isNPP_v2014:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/sh_reshaped_patch_2014/"
        files = glob(ROOT_DIR + "Reshaped_Files/npp/1km/????/06/arc*/*h5")
    elif isGAC_v2014_morning_sat:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
        files = glob(ROOT_DIR + "noaa17/5km/20??/??/*/*h5")
        files = files + glob(ROOT_DIR + "metop*/5km/20??/??/*/*h5")
        #files = glob(ROOT_DIR + "noaa17/5km/20??/1*/*/*noaa*h5")
    elif isGAC_v2014:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
        files = glob(ROOT_DIR + "merged/noaa18*2009*h5")
#        files = files + glob(ROOT_DIR + "merged/noaa19*2008*h5")
    caObj = CalipsoAvhrrTrackObject()
    for filename in files:
        print  os.path.basename(filename)
        caObj +=readCaliopAvhrrMatchObj(filename)#, var_to_skip='segment')
    make_boxplot(caObj)    
    """
