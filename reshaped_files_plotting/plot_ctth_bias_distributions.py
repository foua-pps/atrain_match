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
matplotlib.rcParams.update({'font.size': 12})
from get_flag_info import get_calipso_clouds_of_type_i
from get_flag_info import (get_semi_opaque_info_pps2014,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_low_clouds)


def make_profileplot(caObj, month):
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    height_c = 1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0]


    #USE_ONLY_PIXELS_WHERE_PPS_AND_MODIS_C6_HAVE_VALUES
    height_mlvl2 = caObj.modis.all_arrays['height']#+caObj.calipso.all_arrays['elevation']
    height_pps = caObj.avhrr.all_arrays['imager_ctth_m_above_seasurface']
    height_old = caObj.avhrr.all_arrays['ctthold_height']+caObj.calipso.all_arrays['elevation']
    height_nna1 = caObj.avhrr.all_arrays['ctthnna1_height']+caObj.calipso.all_arrays['elevation']
    use = caObj.calipso.all_arrays['layer_top_altitude'][:,0]>0
    use = np.logical_and(use, height_mlvl2>-1)
    use = np.logical_and(use, height_mlvl2<45000)
    use = np.logical_and(use, height_pps>-1)        
    use = np.logical_and(use, height_pps<45000)
    use = np.logical_and(use, height_nna1>-1)        
    use = np.logical_and(use, height_nna1<45000)
    use = np.logical_and(use, height_old>-1)        
    use = np.logical_and(use, height_old<45000)
      
    thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.30, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    very_thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.10, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    thin_top = np.logical_and(caObj.calipso.all_arrays['number_layers_found']>1, thin)
    thin_1_lay = np.logical_and(caObj.calipso.all_arrays['number_layers_found']==1, thin)

    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)

    use_part = True
    part = "all"
    if use_part:
        part = "single_layer_sea_below_60_in5km_as_art"
        use = np.logical_and(use,caObj.calipso.all_arrays['number_layers_found']==1)
        use = np.logical_and(use,np.abs(caObj.calipso.all_arrays['latitude'])<60)
        use = np.logical_and(use,np.equal(caObj.calipso.igbp_surface_type,17))
        #use = np.logical_and(use,height_c<3000+caObj.calipso.all_arrays['elevation'])
        use = np.logical_and(use,caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0)
        use = np.logical_and(use,
                             caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']==
                             caObj.calipso.all_arrays['total_optical_depth_5km'])
        part = "single_layer_sea_below_60_in5km_all_od_top"
        low = np.logical_and(low_clouds,use)
        #low = np.logical_and(height_c<low_clouds,use)
        low = np.logical_and(use,height_c<(3000+caObj.calipso.all_arrays['elevation']))
        medium = np.logical_and(medium_clouds,use)
        high = np.logical_and(height_c>8000,use)

    c_all = np.logical_or(high,np.logical_or(low,medium))
    high_very_thin = np.logical_and(high, very_thin)
    high_thin = np.logical_and(high, np.logical_and(~very_thin,thin))
    #high_cirrus = np.logical_and(height_c>8000,
    high_thick = np.logical_and(high, ~thin)
    #print "thin, thick high", np.sum(high_thin), np.sum(high_thick) 
    pps_bias = height_pps - height_c
    nna1_bias = height_nna1 - height_c
    old_bias = height_old - height_c
    mlvl2_bias = height_mlvl2 - height_c
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    delta_h = 100.0   
    def plot_one(bias_v, selection, bmin, bmax, color, label):
       bins = np.arange(bmin*1000,bmax*1000,delta_h)
       hist_heights,bins = np.histogram(bias_v[selection],bins=bins)
       n_pix = np.sum(selection)
       hist_heights = hist_heights*100.0/n_pix
       plt.plot(0.001*(bins[0:-1]+delta_h*0.5), hist_heights,
                color,label=label)     
       plt.xlabel("Bias imager height - CALIOP height (km) ")
       plt.ylabel("Percent of data")
       ax.set_xlim(bmin,bmax)

    fig = plt.figure(figsize = (5,12))        
    ax = fig.add_subplot(313)
    plt.title("Low")
    plt.xticks(np.arange(-5,5,1.0))
    plt.yticks(np.arange(0,14,1.0))    
    plot_one(pps_bias, low, -14,14, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, low, -14,14, 'r', "MODIS_C6")
    plot_one(old_bias, low, -14,14, 'b', "PPS_old")
    plot_one(nna1_bias, low, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-5,5)
    ax.grid(True)
    ax = fig.add_subplot(312)
    plt.title("Medium")
    plt.xticks(np.arange(-8,8,2.0))
    plt.yticks(np.arange(0,12,1.0))
    plot_one(pps_bias, medium, -12,12, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, medium, -12,12, 'r', "MODIS_C6")
    plot_one(old_bias, medium, -12,12, 'b', "PPS_old")
    plot_one(nna1_bias, medium, -12,12, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-7,7)
    ax.grid(True)
    ax = fig.add_subplot(311)
    plt.title("High")
    plt.xticks(np.arange(-12,12,2.0))
    plt.yticks(np.arange(0,4,0.5))
    plot_one(pps_bias, high, -14,14, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, high, -14,14, 'r', "MODIS_C6")
    plot_one(old_bias, high, -14,14, 'b', "PPS_old")
    plot_one(nna1_bias, high, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)

 
    ax.set_xlim(-12,12)
    ax.set_ylim(0,4.0)
    ax.grid(True)

    #plt.show()
   
    #plt.title("%s MAE = %3.0f"%(name,MAE))
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX/ctth_profile_%s_%s.png"%(month,part))



def investigate_nn_ctth_modis_lvl2():
    #november
 
    ROOT_DIR = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "global_modis_14th_created20170330/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")
    caObj = CalipsoAvhrrTrackObject()
    name=""
    for month in [ "08", "07", "09"]:#"06", "09"]:#, "11", "05", "02","03","04","07","08","10","12"]:#, "01"]:    
        print ROOT_DIR%(month)
        files = glob(ROOT_DIR%(month))
        name+=month 
        #caObj = CalipsoAvhrrTrackObject()
        for filename in files:
            #print filename
            caObj +=  readCaliopAvhrrMatchObj(filename)
    make_profileplot(caObj, month=name)
        
if __name__ == "__main__":
    investigate_nn_ctth_modis_lvl2()
  
