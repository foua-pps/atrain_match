"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CloudsatAvhrrTrackObject,
                            readCloudsatAvhrrMatchObj,
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

def make_profileplot_cloudsat(clsatObj, month):

    height_c = clsatObj.cloudsat.all_arrays['clsat_max_height']
    if clsatObj is not None:
        clsat_max_height = -9 + 0*np.zeros(clsatObj.cloudsat.latitude.shape)
        for i in range(125):
            height = clsatObj.cloudsat.Height[:,i]
            cmask_ok = clsatObj.cloudsat.CPR_Cloud_mask[:,i]
            top_height = height+120
            #top_height[height<240*4] = -9999 #Do not use not sure why these are not used Nina 20170317
            is_cloudy = cmask_ok > 30
            top_height[~is_cloudy] = -9999
            clsat_max_height[clsat_max_height<top_height] =  top_height[clsat_max_height<top_height]
        height_c = clsat_max_height

    low_clouds = np.logical_and(height_c<clsatObj.avhrr.all_arrays['segment_nwp_h680'], height_c>-9)
    medium_clouds = np.logical_and(height_c>=clsatObj.avhrr.all_arrays['segment_nwp_h680'], 
                                   height_c<=clsatObj.avhrr.all_arrays['segment_nwp_h440'])
    high_clouds = np.logical_and(height_c>clsatObj.avhrr.all_arrays['segment_nwp_h440'], height_c>-9)
    #USE_ONLY_PIXELS_WHERE_PPS_AND_MODIS_C6_HAVE_VALUES
    elevation = clsatObj.cloudsat.all_arrays['elevation']
    elevation[elevation<0] = 0
    height_mlvl2 = clsatObj.modis.all_arrays['height']+elevation
    #height_pps = clsatObj.avhrr.all_arrays['imager_ctth_m_above_seasurface']
    height_pps = clsatObj.avhrr.all_arrays['ctthnnant_height']+elevation
    height_old = clsatObj.avhrr.all_arrays['ctthold_height']+elevation
    height_nna1 = clsatObj.avhrr.all_arrays['ctthnna1nt_height']+elevation
    use =  height_c>0
    print len(height_c[use])*1.0/len(height_c)
    use = np.logical_and(use, height_mlvl2>-1)
    use = np.logical_and(use, height_mlvl2<45000)
    use = np.logical_and(use, height_pps>-1)        
    use = np.logical_and(use, height_pps<45000)
    #use = np.logical_and(use, height_nna1>-1)        
    #use = np.logical_and(use, height_nna1<45000)
    #use = np.logical_and(use, height_old>-1)        
    #use = np.logical_and(use, height_old<45000)                                
    #use = np.logical_and(use,clsatObj.avhrr.all_arrays['ctthnna1_height']>-1)
    #use = np.logical_and(use,clsatObj.avhrr.all_arrays['ctthold_height']>-1)
    use = np.logical_and(use,clsatObj.avhrr.all_arrays['imager_ctth_m_above_seasurface']>-1)
    use = np.logical_and(use,clsatObj.avhrr.all_arrays['ctthnna1_height']<45000)
    use = np.logical_and(use,clsatObj.avhrr.all_arrays['ctthold_height']<45000)
    use = np.logical_and(use,clsatObj.avhrr.all_arrays['imager_ctth_m_above_seasurface']<45000)
    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)
    over_land =np.logical_and(elevation>0,use)   
    pps_bias = height_pps - height_c
    nna1_bias = height_nna1 - height_c
    old_bias = height_old - height_c
    mlvl2_bias = height_mlvl2 - height_c
    print "CLOUDSAT"
    for bias_v, name in zip([old_bias, mlvl2_bias, pps_bias,nna1_bias],
                            ["CTTHold", "MODIS-C6", "NN-AVHRR", "NN-AVHRR1" ]):
        print "%s & %3.0f & %3.0f & %3.0f & %3.0f \\\\"%(
            name, 
            np.mean(np.abs(bias_v[use])), 
            np.mean(np.abs(bias_v[low])),
            np.mean(np.abs(bias_v[medium])),
            np.mean(np.abs(bias_v[high])),
            #np.mean(np.abs(bias_v[over_land])),
        )
    n_all = len(bias_v[use])*1.0
    print "%s & %d & %d & %d & %d \\\\"%(name, 
                                         n_all, 
                                         len(bias_v[low]),
                                         len(bias_v[medium]),
                                         len(bias_v[high]))
    print "%s & %d & %3.1f & %3.1f & %3.1f \\\\"%(name, 
                                         n_all,
                                         len(bias_v[low])*100/n_all,
                                         len(bias_v[medium])*100/n_all,
                                         len(bias_v[high])*100/n_all)
    from matplotlib import rcParams
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
       plt.xlabel("Bias imager height - CloudSat height (km) ")
       plt.ylabel("Percent of data")
       ax.set_xlim(bmin,bmax)
    fig = plt.figure(figsize = (5.5,12))        
    ax = fig.add_subplot(313)
    #plt.title("Low")
    plt.text(0.02, 0.90, "c. Low clouds", fontsize=14,
             transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-5,5,1.0))
    plt.yticks(np.arange(0,14,1.0))    
    plot_one(pps_bias, low, -14,14, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, low, -14,14, 'r', "MODIS-C6")
    plot_one(old_bias, low, -14,14, 'b', "PPS-v2014")
    #plot_one(nna1_bias, low, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-5,5)
    ax.set_ylim(0,12.0)
    ax.grid(True)
    ax = fig.add_subplot(312)
    #plt.title("Medium")
    plt. text(0.02, 0.90, "b. Medium clouds", fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-8,8,2.0))
    plt.yticks(np.arange(0,12,1.0))
    plot_one(pps_bias, medium, -12,12, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, medium, -12,12, 'r', "MODIS-C6")
    plot_one(old_bias, medium, -12,12, 'b', "PPS-v2014")
    #plot_one(nna1_bias, medium, -12,12, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-7,7)
    ax.set_ylim(0,6.0)
    ax.grid(True)
    ax = fig.add_subplot(311)
    #plt.title("High")
    plt.text(0.02, 0.90, "a. High clouds", fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-12,12,2.0))
    plt.yticks(np.arange(0,4,0.5))
    plot_one(pps_bias, high, -14,14, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, high, -14,14, 'r', "MODIS-C6")
    plot_one(old_bias, high, -14,14, 'b', "PPS-v2014")
    #plot_one(nna1_bias, high, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)  
    #, bbox_to_anchor=(1.1, 1.05)
    ax.set_xlim(-12,12)
    ax.set_ylim(0,4.0)
    ax.grid(True)
    #plt.show()   
    #plt.title("%s MAE = %3.0f"%(name,MAE))
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX_cloudsat/ctth_profile_1st_%s_cloudsat_all.png"%(month))


def make_profileplot(caObj, month):
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    height_c = 1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0]

    #USE_ONLY_PIXELS_WHERE_PPS_AND_MODIS_C6_HAVE_VALUES
    elevation = caObj.calipso.all_arrays['elevation']
    elevation[elevation<0] = 0
    height_mlvl2 = caObj.modis.all_arrays['height']#+elevation
    #height_pps = caObj.avhrr.all_arrays['imager_ctth_m_above_seasurface']
    height_pps = caObj.avhrr.all_arrays['ctthnnant_height']+elevation
    height_old = caObj.avhrr.all_arrays['ctthold_height']+elevation
    height_nna1 = caObj.avhrr.all_arrays['ctthnna1nt_height']+elevation
    use = caObj.calipso.all_arrays['layer_top_altitude'][:,0]>0
    print len(height_c[use])*1.0/len(height_c)
    use = np.logical_and(use, height_mlvl2>-1)
    use = np.logical_and(use, height_mlvl2<45000)
    use = np.logical_and(use, height_pps>-1)        
    use = np.logical_and(use, height_pps<45000)
    use = np.logical_and(use, height_nna1>-1)        
    use = np.logical_and(use, height_nna1<45000)
    use = np.logical_and(use, height_old>-1)        
    use = np.logical_and(use, height_old<45000)
      
    """
    thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.30, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    very_thin = np.logical_and(caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']<0.10, 
                          caObj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>0) 
    thin_top = np.logical_and(caObj.calipso.all_arrays['number_layers_found']>1, thin)
    thin_1_lay = np.logical_and(caObj.calipso.all_arrays['number_layers_found']==1, thin)
    """

    low = np.logical_and(low_clouds,use)
    medium = np.logical_and(medium_clouds,use)
    high = np.logical_and(high_clouds,use)

    use_part = False
    part = "all"
    if use_part:
        part = "single_layer_sea_below_60_in5km_as_art"
        use = np.logical_and(
            use,caObj.calipso.all_arrays['number_layers_found']==1)
        use = np.logical_and(
            use,np.abs(caObj.calipso.all_arrays['latitude'])<60)
        use = np.logical_and(
            use,np.equal(caObj.calipso.igbp_surface_type,17))
        #use = np.logical_and(use,height_c<3000+caObj.calipso.all_arrays['elevation'])
        use = np.logical_and(use,caObj.calipso.all_arrays[
            'feature_optical_depth_532_top_layer_5km']>0)
        use = np.logical_and(use,
                             caObj.calipso.all_arrays[
                                 'feature_optical_depth_532_top_layer_5km']==
        caObj.calipso.all_arrays['total_optical_depth_5km'])
        part = "single_layer_sea_below_60_in5km_all_od_top"
        low = np.logical_and(low_clouds,use)
        #low = np.logical_and(height_c<low_clouds,use)
        low = np.logical_and(
            use,height_c<(3000+caObj.calipso.all_arrays['elevation']))
        medium = np.logical_and(medium_clouds,use)
        high = np.logical_and(height_c>8000,use)

    ##c_all = np.logical_or(high,np.logical_or(low,medium))
    ##high_very_thin = np.logical_and(high, very_thin)
    ##high_thin = np.logical_and(high, np.logical_and(~very_thin,thin))
    ##high_cirrus = np.logical_and(height_c>8000,
    ##high_thick = np.logical_and(high, ~thin)
    ##print "thin, thick high", np.sum(high_thin), np.sum(high_thick) 
    print len(height_c)
    pps_bias = height_pps - height_c
    nna1_bias = height_nna1 - height_c
    old_bias = height_old - height_c
    mlvl2_bias = height_mlvl2 - height_c
    for bias_v, name in zip([old_bias, mlvl2_bias, pps_bias,nna1_bias],
                            ["CTTHold", "MODIS-C6", "NN-AVHRR", "NN-AVHRR1" ]):
        print "%s & %3.0f & %3.0f & %3.0f & %3.0f \\\\"%(name, 
                           np.mean(np.abs(bias_v[use])), 
                           np.mean(np.abs(bias_v[low])),
                           np.mean(np.abs(bias_v[medium])),
                           np.mean(np.abs(bias_v[high])))
    n_all = len(bias_v[use])*1.0
    print "%s & %d & %d & %d & %d \\\\"%(name, 
                                         n_all, 
                                         len(bias_v[low]),
                                         len(bias_v[medium]),
                                         len(bias_v[high]))
    print "%s & %d & %3.1f & %3.1f & %3.1f \\\\"%(name, 
                                         n_all,
                                         len(bias_v[low])*100/n_all,
                                         len(bias_v[medium])*100/n_all,
                                         len(bias_v[high])*100/n_all)
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
       plt.xlabel("Bias imager height - CALIPSO height (km) ")
       plt.ylabel("Percent of data")
       ax.set_xlim(bmin,bmax)
    fig = plt.figure(figsize = (5.5,12))        
    ax = fig.add_subplot(313)
    #plt.title("Low")
    plt.text(0.02, 0.90, "c. Low clouds", fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-5,5,1.0))
    plt.yticks(np.arange(0,14,1.0))    
    plot_one(pps_bias, low, -14,14, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, low, -14,14, 'r', "MODIS-C6")
    plot_one(old_bias, low, -14,14, 'b', "PPS-v2014")
    #plot_one(nna1_bias, low, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-5,5)
    ax.set_ylim(0,12.0)
    ax.grid(True)
    ax = fig.add_subplot(312)
    #plt.title("Medium")
    plt.text(0.02, 0.90, "b. Medium clouds", fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-8,8,2.0))
    plt.yticks(np.arange(0,12,1.0))
    plot_one(pps_bias, medium, -12,12, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, medium, -12,12, 'r', "MODIS-C6")
    plot_one(old_bias, medium, -12,12, 'b', "PPS-v2014")
    #plot_one(nna1_bias, medium, -12,12, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-7,7)
    ax.set_ylim(0,6.0)
    ax.grid(True)
    ax = fig.add_subplot(311)
    #plt.title("High")
    plt.text(0.02, 0.90, "a. High clouds", fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-12,12,2.0))
    plt.yticks(np.arange(0,4,0.5))
    plot_one(pps_bias, high, -14,14, 'k', "NN-AVHRR")
    plot_one(mlvl2_bias, high, -14,14, 'r', "MODIS-C6")
    plot_one(old_bias, high, -14,14, 'b', "PPS-v2014")
    #plot_one(nna1_bias, high, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)                        
    ax.set_xlim(-12,12)
    ax.set_ylim(0,4.0)
    ax.grid(True)
    #plt.show()   
    #plt.title("%s MAE = %3.0f"%(name,MAE))
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX/ctth_profile_1st_%s_%s.png"%(month,part))



def investigate_nn_ctth_modis_lvl2_cloudsat():
    #november
 
    ROOT_DIR = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "global_modis_01st_created20170504/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
        #"global_modis_14th_created20170330/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")
    clsatObj = CloudsatAvhrrTrackObject()
    name=""
    for month in [ "02", "04","06", "08","10", "12" ]:  #[ "06", "09"]: 
        print ROOT_DIR%(month)
        files = glob(ROOT_DIR%(month))
        name+=month 
        #clsatObj = CloudsatAvhrrTrackObject()
        clsatObj_new = CloudsatAvhrrTrackObject()
        for filename in files:
            print filename
            clsatObj_new +=  readCloudsatAvhrrMatchObj(filename)
        clsatObj +=  clsatObj_new
        make_profileplot_cloudsat(clsatObj_new, month=month)
    make_profileplot_cloudsat(clsatObj, month=name)

def investigate_nn_ctth_modis_lvl2():
    #november
 
    ROOT_DIR = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "global_modis_01st_created20170504/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")
        #"global_modis_14th_created20170330/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")
    caObj = CalipsoAvhrrTrackObject()
    name=""
    for month in [ "02", "04","06", "08","10", "12" ]:  #[ "06", "09"]:    
        print ROOT_DIR%(month)
        files = glob(ROOT_DIR%(month))
        name+=month 
        #caObj = CalipsoAvhrrTrackObject()
        caObj_new = CalipsoAvhrrTrackObject()
        for filename in files:
            #print filename
            caObj_new += readCaliopAvhrrMatchObj(filename)
        make_profileplot(caObj_new, month=month)
        caObj +=  caObj_new

    make_profileplot(caObj, month=name)
        
if __name__ == "__main__":
    #investigate_nn_ctth_modis_lvl2()
    investigate_nn_ctth_modis_lvl2_cloudsat()
    investigate_nn_ctth_modis_lvl2()


