"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)
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

labels=["low","medium","high-all","high-thick\n od>0.4","high-thin \n 0.1<od<0.4","high-vthin\n od<0.1"],
def make_violinplot(caObj, name, modis_lvl2=False):
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    height_c = (1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] - 
                caObj.calipso.all_arrays['elevation'])
    if modis_lvl2:
        height_pps = caObj.modis.all_arrays['height']
    else:
        height_pps = caObj.imager.all_arrays['ctth_height']
        print "min/max height", np.min( height_pps), np.max(height_pps)
    thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.30, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    very_thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.10, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    thin_top = np.logical_and(caObj.calipso.all_arrays['number_layers_found']>1, thin)
    thin_1_lay = np.logical_and(caObj.calipso.all_arrays['number_layers_found']==1, thin)
    use = np.logical_and(height_pps >-1,
                         caObj.calipso.all_arrays['layer_top_altitude'][:,0]>=0)
    use = np.logical_and(height_pps <45000,use)
    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)
    c_all = np.logical_or(high,np.logical_or(low,medium))
    high_very_thin = np.logical_and(high, very_thin)
    high_thin = np.logical_and(high, np.logical_and(~very_thin,thin))
    high_thick = np.logical_and(high, ~thin)
    #print "thin, thick high", np.sum(high_thin), np.sum(high_thick) 
    bias = height_pps - height_c
    abias = np.abs(bias)
    #abias[abias>2000]=2000
    print name.ljust(30, " "), "%3.1f"%(np.mean(abias[c_all])), "%3.1f"%(np.mean(abias[low])),"%3.1f"%(np.mean(abias[medium])),"%3.1f"%(np.mean(abias[high]))

    c_all = np.logical_or(np.logical_and(~very_thin,high),np.logical_or(low,medium))
    number_of = np.sum(c_all)
     
    #print name.ljust(30, " "), "%3.1f"%(np.sum(abias[c_all]<250)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<500)*100.0/number_of),  "%3.1f"%(np.sum(abias[c_all]<1000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<1500)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<2000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<3000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<4000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<5000)*100.0/number_of)
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    fig = plt.figure(figsize = (6,9))        
    ax = fig.add_subplot(111)
    plt.xticks(rotation=70)
    ax.fill_between(np.arange(0,8),-1500,1500, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0,8),2000,20000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0,8),-2000,-20000, facecolor='red', alpha=0.2)
    for y_val in [-5,-4,-3,-2,2,3,4,5]:
        plt.plot(np.arange(0,8), y_val*1000 + 0*np.arange(0,8),':k', alpha=0.4)
        plt.plot(np.arange(0,8), -10*1000 + 0*np.arange(0,8),':k', alpha=0.4)
    plt.plot(np.arange(0,8), 0 + 0*np.arange(0,8),':k', alpha=0.4)
    bplot = ax.violinplot([bias[low],bias[medium],bias[high],bias[high_thick],bias[high_thin],bias[high_very_thin]], 
                          widths=0.9,showextrema=False, showmedians=True)
    ax.set_ylim(-20000,20000)
    plt.title(name)
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX/ctth_violin_%s_5_95_filt.png"%(name))


def make_violinplot_temperature(caObj, name, modis_lvl2=False):
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    temp_c = caObj.calipso.all_arrays['layer_top_temperature'][:,0] +273.15 
    if modis_lvl2:
        temp_pps = caObj.modis.all_arrays['temperature']
    else:
        temp_pps = caObj.imager.all_arrays['ctth_temperature']  
    if modis_lvl2:
        height_pps = caObj.modis.all_arrays['height']
    else:
        height_pps = caObj.imager.all_arrays['ctth_height']
    thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.30, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    very_thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.10, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    thin_top = np.logical_and(caObj.calipso.all_arrays['number_layers_found']>1, thin)
    thin_1_lay = np.logical_and(caObj.calipso.all_arrays['number_layers_found']==1, thin)
    use = np.logical_and(temp_pps >100,
                         caObj.calipso.all_arrays['layer_top_altitude'][:,0]>=0)
    use = np.logical_and(height_pps <45000,use)
    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)
    c_all = np.logical_or(high,np.logical_or(low,medium))
    high_very_thin = np.logical_and(high, very_thin)
    high_thin = np.logical_and(high, np.logical_and(~very_thin,thin))
    high_thick = np.logical_and(high, ~thin)
    #print "thin, thick high", np.sum(high_thin), np.sum(high_thick) 
    bias = temp_pps - temp_c
    abias = np.abs(bias)
    #abias[abias>2000]=2000
    print name.ljust(30, " "), "%3.1f"%(np.mean(abias[c_all])), "%3.1f"%(np.mean(abias[low])),"%3.1f"%(np.mean(abias[medium])),"%3.1f"%(np.mean(abias[high]))

    c_all = np.logical_or(np.logical_and(~very_thin,high),np.logical_or(low,medium))
    number_of = np.sum(c_all)
     
    #print name.ljust(30, " "), "%3.1f"%(np.sum(abias[c_all]<250)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<500)*100.0/number_of),  "%3.1f"%(np.sum(abias[c_all]<1000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<1500)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<2000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<3000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<4000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<5000)*100.0/number_of)
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    fig = plt.figure(figsize = (6,9))        
    ax = fig.add_subplot(111)
    plt.xticks(rotation=70)
    ax.fill_between(np.arange(0,8),-7.5,7.5, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0,8),10,200, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0,8),-200,-10, facecolor='red', alpha=0.2)
    for y_val in [-5,-4,-3,-2,-1,1,2,3,4,5]:
        plt.plot(np.arange(0,8), y_val*20 + 0*np.arange(0,8),':k', alpha=0.4)
    plt.plot(np.arange(0,8), 0 + 0*np.arange(0,8),':k', alpha=0.4)
    bplot = ax.violinplot([bias[low],bias[medium],bias[high],bias[high_thick],bias[high_thin],bias[high_very_thin]],
                          widths=0.9,showextrema=False, showmedians=True)
    ax.set_ylim(-150,150)
    plt.title(name)
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX/ctth_violin_temperature_%s_5_95_filt.png"%(name))

def make_violinplot_pressure(caObj, name, modis_lvl2=False):
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    pressure_c = caObj.calipso.all_arrays['layer_top_pressure'][:,0]
    if modis_lvl2:
        pressure_pps = caObj.modis.all_arrays['pressure']
    else:
        pressure_pps = 0.01*caObj.imager.all_arrays['ctth_pressure']
    if modis_lvl2:
        height_pps = caObj.modis.all_arrays['height']
    else:
        height_pps = caObj.imager.all_arrays['ctth_height']
    thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.30, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    very_thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.10, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    thin_top = np.logical_and(caObj.calipso.all_arrays['number_layers_found']>1, thin)
    thin_1_lay = np.logical_and(caObj.calipso.all_arrays['number_layers_found']==1, thin)
    use = np.logical_and(pressure_pps >0,
                         caObj.calipso.all_arrays['layer_top_altitude'][:,0]>=0)
    use = np.logical_and(height_pps <45000,use)
    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)
    c_all = np.logical_or(high,np.logical_or(low,medium))
    high_very_thin = np.logical_and(high, very_thin)
    high_thin = np.logical_and(high, np.logical_and(~very_thin,thin))
    high_thick = np.logical_and(high, ~thin)
    #print "thin, thick high", np.sum(high_thin), np.sum(high_thick) 
    bias = pressure_pps - pressure_c
    abias = np.abs(bias)
    #abias[abias>2000]=2000
    print name.ljust(30, " "), "%3.1f"%(np.mean(abias[c_all])), "%3.1f"%(np.mean(abias[low])),"%3.1f"%(np.mean(abias[medium])),"%3.1f"%(np.mean(abias[high]))

    c_all = np.logical_or(np.logical_and(~very_thin,high),np.logical_or(low,medium))
    number_of = np.sum(c_all)
     
    #print name.ljust(30, " "), "%3.1f"%(np.sum(abias[c_all]<250)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<500)*100.0/number_of),  "%3.1f"%(np.sum(abias[c_all]<1000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<1500)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<2000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<3000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<4000)*100.0/number_of), "%3.1f"%(np.sum(abias[c_all]<5000)*100.0/number_of)
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    fig = plt.figure(figsize = (6,9))        
    ax = fig.add_subplot(111)
    plt.xticks(rotation=70)
    ax.fill_between(np.arange(0,8),-150,150, facecolor='green', alpha=0.2)
    ax.fill_between(np.arange(0,8),200,2000, facecolor='red', alpha=0.2)
    ax.fill_between(np.arange(0,8),-2000,-200, facecolor='red', alpha=0.2)
    for y_val in [-6,-4,-2,2,4,6,8,-8]:
        plt.plot(np.arange(0,8), y_val*100 + 0*np.arange(0,8),':k', alpha=0.4)
    plt.plot(np.arange(0,8), 0 + 0*np.arange(0,8),':k', alpha=0.4)
    bplot = ax.violinplot([bias[low],bias[medium],bias[high],bias[high_thick],bias[high_thin],bias[high_very_thin]],
                          widths=0.9,showextrema=False, showmedians=True)
    ax.set_ylim(-1200,1200)
    plt.title(name)
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX/ctth_violin_pressure_%s_5_95_filt.png"%(name))


def investigate_nn_ctth_modis_lvl2():
    #november
 
    ROOT_DIR_MODIS_nn_imager = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "global_modis_14th_created20170324/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")

    ROOT_DIR_MODIS_old = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "global_modis_14th_created20161108/Reshaped_Files/merged/*%s*h5")

    for month in [ "06", "09", "01"]:    
        for ROOT_DIR, name in zip(
                [ROOT_DIR_MODIS_nn_imager,
                 ROOT_DIR_MODIS_nn_imager,
                 ROOT_DIR_MODIS_old], 
                ["modis_nnIMAGER",
                 "modis_lvl2_C6",
                 "modis_CTTHold"]):
            name = "%s_%s"%(name, month)
            print ROOT_DIR
            files = glob(ROOT_DIR%(month))
            caObj = CalipsoImagerTrackObject()
            for filename in files:
                #print filename
                caObj +=  readCaliopImagerMatchObj(filename)
            modis_lvl2 = False
            if "modis_lvl2"  in name:
                modis_lvl2 = True
            make_violinplot(caObj, name, modis_lvl2=modis_lvl2 )
            make_violinplot_pressure(caObj, name, modis_lvl2=modis_lvl2) 
            make_violinplot_temperature(caObj, name, modis_lvl2=modis_lvl2) 
        
if __name__ == "__main__":
    investigate_nn_ctth_modis_lvl2()
  
