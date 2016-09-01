#!/usr/bin/env python
# -*- coding: utf-8 -*-



"""Read all matched data and make some plotting of thresholds
0 = not determined
1 = clean marine
2 = dust
3 = polluted continental
4 = clean continental
5 = polluted dust
6 = smoke
7 = other
"""

from glob import glob
import os.path
import os
import numpy as np
from scipy import histogram
import re
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
import matplotlib.pyplot as plt

def get_feature_type(caObj, aerosol_type=2):
   cal_vert_feature = np.ones(caObj.calipso_aerosol.layer_top_altitude[::,0].shape)*-9
   cavf = caObj.calipso_aerosol.feature_classification_flags[::,0]
   feature_array = (4*np.bitwise_and(np.right_shift(cavf,11),1) + 
                    2*np.bitwise_and(np.right_shift(cavf,10),1) + 
                    np.bitwise_and(np.right_shift(cavf,9),1))
   cal_vert_feature = np.where(np.not_equal(cavf,1),feature_array[::], cal_vert_feature[::])
   #print cal_vert_feature
   isDust = cal_vert_feature==aerosol_type
   return isDust

def  my_plot_histogram_hexbin(x, y, 
                              bins='log', 
                              gridsize=None,
                              ttitle="", txlabel='x', tylabel='y', ymin=None,
                              file_name=None):
    """Plot a 2d-histogram using hexbins
    """
    xdata = x.ravel()
    ydata = y.ravel()
    if gridsize is not None:
        plt.hexbin(xdata,ydata,bins='log', cmap=plt.cm.gray_r,gridsize=gridsize)
    else:
        plt.hexbin(xdata,ydata,bins='log', cmap=plt.cm.gray_r)    
    plt.colorbar()
    plt.title(ttitle, color='k')
    #if ymin is not None:
   # ax.set_xlim(-1,100)
   # ax.set_ylim(-1,100)
    plt.xlabel(txlabel)
    plt.ylabel(tylabel)
    return 1



def get_calipso_cloudy_and_aerosl(caObj):    
   from scipy import stats, ndimage    
   nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
   meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
   isAerosol = caObj.calipso_aerosol.all_arrays['number_layers_found']>0
   isClear = np.logical_and(meancl<0.1,
                            caObj.calipso.all_arrays['total_optical_depth_5km']<0.0)
   isClear = np.logical_and(isClear, np.not_equal(isAerosol, True))
   isCloudy = np.logical_or( #and in arctic not many small cumulus there!
      #only use those clouds in 1km data
      caObj.calipso.all_arrays['number_layers_found']>0, 
      caObj.calipso.all_arrays['total_optical_depth_5km']>0.5)

   isCloudyAerosolMix = np.logical_and(isCloudy, isAerosol)
   isCloudy = np.logical_and(isCloudy, np.not_equal(isCloudyAerosolMix, True))
   isAerosol = np.logical_and(isAerosol, np.not_equal(isCloudyAerosolMix, True))
   
   return isCloudy, isClear, isAerosol, isCloudyAerosolMix

def do_one_fig(x_vector, y_vector, isClear, isCloudy, isAerosol, isMix,
               x_lab, y_lab,
               filename='temp_plt.png'):
    fig = plt.figure(figsize = (14,12))
    ax = fig.add_subplot(221)
    my_plot_histogram_hexbin(x_vector[isClear], y_vector[isClear], gridsize=30,
                             txlabel=x_lab, tylabel =y_lab)
    ax.set_title('Clear')
    ax = fig.add_subplot(222)
    my_plot_histogram_hexbin(x_vector[isAerosol], y_vector[isAerosol], gridsize=30,
                             txlabel=x_lab, tylabel =y_lab)
    ax.set_title('Aerosol')
    
    ax = fig.add_subplot(223)
    my_plot_histogram_hexbin(x_vector[isCloudy], y_vector[isCloudy], gridsize=30, 
                             txlabel=x_lab, tylabel =y_lab)
    ax.set_title('Cloudy')
    ax = fig.add_subplot(224)
    my_plot_histogram_hexbin(x_vector[isMix], y_vector[isMix], gridsize=30,
                             txlabel=x_lab, tylabel =y_lab)
    ax.set_title('Aerosol & Cloud Mix')
    plt.savefig(filename)

if __name__ == "__main__":
   ROOT_DIR = "/home/a001865/git/atrain_match/modis_merged_??.h5"
   files = glob(ROOT_DIR)

   caObj = CalipsoAvhrrTrackObject()
   for filename in files:
      print os.path.basename(filename)
      newObj = readCaliopAvhrrMatchObj(filename)
      caObj= caObj + newObj

   isCloudy, isClear, isAerosol, isCloudyAerosolMix = get_calipso_cloudy_and_aerosl(caObj)
   from get_pps_flag_info import get_land_coast_sea_info_pps2014
   flgs = get_land_coast_sea_info_pps2014(caObj.avhrr.all_arrays['cloudtype_conditions'])
   (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag) = flgs
   
   sunz = caObj.avhrr.all_arrays['sunz']
   use = sunz<80
   modis_1 = caObj.avhrr.all_arrays['r06micron'][use]
   thr_r06 = caObj.avhrr.all_arrays['thr_r06'][use]
   t11 =  caObj.avhrr.all_arrays['bt11micron'][use]
   t12 =  caObj.avhrr.all_arrays['bt12micron'][use]
   t37 =  caObj.avhrr.all_arrays['bt37micron'][use]
   modis_3 = caObj.avhrr.all_arrays['modis_3'][use]
   modis_8 = caObj.avhrr.all_arrays['modis_8'][use]
   modis_5 = caObj.avhrr.all_arrays['modis_5'][use]
   modis_15 = caObj.avhrr.all_arrays['modis_15'][use]
   modis_13hi = caObj.avhrr.all_arrays['modis_13hi'][use]
   modis_13lo = caObj.avhrr.all_arrays['modis_13lo'][use]
   modis8_5 = modis_8/modis_5
   modis3_5 = modis_3/modis_5
   modis1_5 = modis_1/modis_5
   modis8_1 = modis_8/modis_1
   modis3_1 = modis_3/modis_1
   modis8_15 = modis_8/modis_15
   modis8_13hi = modis_8/modis_13hi
   modis8_13lo = modis_8/modis_13lo
   t11t12 = t11-t12
   t37t12 = t37-t12
   t37t12[t37t12<-10]=-1.0
   isClear = isClear[use]
   isCloudyAerosolMix = isCloudyAerosolMix[use]
   isCloudy = isCloudy[use]
   isAerosol = isAerosol[use]

   day_and_land = land_flag[use]

   isClearLand = np.logical_and(isClear,day_and_land)
   isCloudyAerosolMixLand = np.logical_and(isCloudyAerosolMix, day_and_land)
   isCloudyLand = np.logical_and(isCloudy, day_and_land)
   isAerosolLand = np.logical_and(isAerosol, day_and_land)

   day_and_sea = sea_flag[use]

   isClearSea = np.logical_and(isClear, day_and_sea)
   isCloudyAerosolMixSea = np.logical_and(isCloudyAerosolMix, day_and_sea)
   isCloudySea = np.logical_and(isCloudy, day_and_sea)
   isAerosolSea = np.logical_and(isAerosol, day_and_sea)

   day_and_coast = coast_flag[use]

   isClearCoast = np.logical_and(isClear, day_and_coast)
   isCloudyAerosolMixCoast = np.logical_and(isCloudyAerosolMix, day_and_coast)
   isCloudyCoast = np.logical_and(isCloudy, day_and_coast)
   isAerosolCoast = np.logical_and(isAerosol, day_and_coast)

   import numpy as np
   import matplotlib.pyplot as plt

   do_one_fig(t11, t11t12, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 't11', y_lab= 't11-t12',
           filename ='misc_png/aerosol_ir_modis_t11_t11t12.png')
   do_one_fig(modis_1-100*thr_r06, t37t12, 
              isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
              x_lab = 'r06-thr_r06', y_lab= 't37-t12',
              filename ='misc_png/aerosol_ir_modis_r06_t37t12.png')

   do_one_fig(modis_5, modis8_5, isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
              x_lab = 'modis_5 (1.24)', y_lab= 'modis8/modis5 (0.41/1.24)',
              filename ='misc_png/aerosol_log_modis_m5_qm8m5.png')
   do_one_fig(modis_5, modis3_5, isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
              x_lab = 'modis_5 (1.24)', y_lab= 'modis3/modis5 (0.47/1.24)',
              filename ='misc_png/aerosol_log_modis_m5_qm3m5.png')
   do_one_fig(modis_5, modis8_15, isClear, isCloudy, isAerosol, isCloudyAerosolMix,
              x_lab = 'modis_5 (1.24)', y_lab= 'modis8/modis15 (0.41/0.75)',
           filename ='misc_png/aerosol_log_modis_m5_qm8m15.png')
   do_one_fig(modis_5, modis8_13hi, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 'modis_5 (1.24)', y_lab= 'modis8/modis13hi',
           filename ='misc_png/aerosol_log_modis_m5_qm8m13hi.png')
   do_one_fig(modis_5, modis8_13lo, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 'modis_5 (1.24)', y_lab= 'modis8/modis13lo',
           filename ='misc_png/aerosol_log_modis_m5_qm8m13lo.png')
   do_one_fig(modis8_5, modis8_15, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 'modis_8/modis5', y_lab= 'modis8/modis15',
           filename ='misc_png/aerosol_log_modis_qm8m5_qm8m15.png')
   #sea coast land
   do_one_fig(modis3_5, modis3_1, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 'modis3/modis5 (0.47/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm3m5_qm3m1.png')
   do_one_fig(modis3_5, modis3_1, 
           isClearLand, isCloudyLand, isAerosolLand, isCloudyAerosolMixLand, 
           x_lab = 'modis3/modis5 (0.47/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm3m5_qm3m1_land.png')
   do_one_fig(modis3_5, modis3_1, 
           isClearSea, isCloudySea, isAerosolSea, isCloudyAerosolMixSea, 
           x_lab = 'modis3/modis5 (0.47/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm3m5_qm3m1_sea.png')
   do_one_fig(modis3_5, modis3_1, 
           isClearCoast, isCloudyCoast, 
           isAerosolCoast, isCloudyAerosolMixCoast, 
           x_lab = 'modis3/modis5 (0.47/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm3m5_qm3m1_coast.png')
   #sea coast land
   do_one_fig(modis1_5, modis3_1, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 'modis1/modis5 (0.67/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm1m5_qm3m1.png')
   do_one_fig(modis1_5, modis3_1, 
           isClearLand, isCloudyLand, isAerosolLand, isCloudyAerosolMixLand, 
           x_lab = 'modis1/modis5 (0.67/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm1m5_qm3m1_land.png')
   do_one_fig(modis1_5, modis3_1, 
           isClearSea, isCloudySea, isAerosolSea, isCloudyAerosolMixSea, 
           x_lab = 'modis1/modis5 (0.67/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm1m5_qm3m1_sea.png')
   do_one_fig(modis1_5, modis3_1, 
           isClearCoast, isCloudyCoast, 
           isAerosolCoast, isCloudyAerosolMixCoast, 
           x_lab = 'modis1/modis5 (0.67/1.24)', y_lab= 'modis3/moids1 (0.47/0.67)',
           filename ='misc_png/aerosol_log_modis_qm1m5_qm3m1_coast.png')

   do_one_fig(modis8_5, modis8_1, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 'modis_8/modis5 (0.41/1.24)', y_lab= 'modis8/modis1 (0.41/0.65)',
           filename ='misc_png/aerosol_log_modis_qm8m5_qm8m1.png')
   do_one_fig(modis3_5, modis3_1, 
           isClear, isCloudy, isAerosol, isCloudyAerosolMix, 
           x_lab = 'modis_3/modis5 (0.47/1.24)', y_lab= 'modis3/modis1 (0.47/0.65)',
           filename ='misc_png/aerosol_log_modis_qm3m5_qm3m1.png')

   isDust = get_feature_type(caObj, aerosol_type=2)[use]
   isAerosolDust = np.logical_and(isAerosol, isDust)
   isCloudyAerosolMixDust = np.logical_and(isCloudyAerosolMix, isDust)
   do_one_fig(modis8_5, modis8_1, 
           isClear, isCloudy, isAerosolDust, isCloudyAerosolMixDust, 
           x_lab = 'modis_8/modis5 (0.41/1.24)', y_lab= 'modis8/modis1 (0.41/0.65)',
           filename ='misc_png/aerosol_log_dust_modis_qm8m5_qm8m1.png')
   do_one_fig(modis3_5, modis3_1, 
           isClear, isCloudy, isAerosolDust, isCloudyAerosolMixDust, 
              x_lab = 'modis_3/modis5 (0.47/1.24)', y_lab= 'modis3/modis1 (0.47/0.65)',
              filename ='misc_png/aerosol_log_dust_modis_qm3m5_qm3m1.png')
   do_one_fig(modis_5, modis8_5, isClear, isCloudy, isAerosolDust, isCloudyAerosolMixDust, 
           x_lab = 'modis_5 (1.24)', y_lab= 'modis8/modis5 (0.41/1.24)',
              filename ='misc_png/aerosol_log_dust_modis_m5_qm8m5.png')
   do_one_fig(modis_5, modis3_5, isClear, isCloudy, isAerosolDust, isCloudyAerosolMixDust, 
              x_lab = 'modis_5 (1.24)', y_lab= 'modis3/modis5 (0.47/1.24)',
           filename ='misc_png/aerosol_log_dust_modis_m5_qm3m5.png')
   isSmoke = get_feature_type(caObj, aerosol_type=6)[use]
   isAerosolSmoke = np.logical_and(isAerosol, isSmoke)
   isCloudyAerosolMixSmoke = np.logical_and(isCloudyAerosolMix, isSmoke)
   do_one_fig(modis8_5, modis8_1, 
           isClear, isCloudy, isAerosolSmoke, isCloudyAerosolMixSmoke, 
           x_lab = 'modis_8/modis5 (0.41/1.24)', y_lab= 'modis8/modis1 (0.41/0.65) ',
           filename ='misc_png/aerosol_log_smoke_modis_qm8m5_qm8m1.png')
   do_one_fig(modis3_5, modis3_1, 
           isClear, isCloudy, isAerosolSmoke, isCloudyAerosolMixSmoke, 
           x_lab = 'modis_3/modis5 (0.47/1.24)', y_lab= 'modis3/modis1 (0.47/0.65)',
           filename ='misc_png/aerosol_log_smoke_modis_qm3m5_qm3m1.png')
   do_one_fig(modis_5, modis8_5, isClear, isCloudy, isAerosolSmoke, isCloudyAerosolMixSmoke, 
           x_lab = 'modis_5 (1.24)', y_lab= 'modis8/modis5 (0.41/1.24)',
           filename ='misc_png/aerosol_log_smoke_modis_m5_qm8m5.png')
   do_one_fig(modis_5, modis3_5, isClear, isCloudy, isAerosolSmoke, isCloudyAerosolMixSmoke, 
           x_lab = 'modis_5 (1.24)', y_lab= 'modis3/modis5 (0.47/1.24)',
           filename ='misc_png/aerosol_log_smoke_modis_m5_qm3m5.png')
   isMarine = get_feature_type(caObj, aerosol_type=1)[use]
   isAerosolMarine = np.logical_and(isAerosol, isMarine)
   isCloudyAerosolMixMarine = np.logical_and(isCloudyAerosolMix, isMarine)
   do_one_fig(modis8_5, modis8_1, 
           isClear, isCloudy, isAerosolMarine, isCloudyAerosolMixMarine, 
           x_lab = 'modis_8/modis5 (0.41/1.24)', y_lab= 'modis8/modis1 (0.41/0.65)',
           filename ='misc_png/aerosol_log_marine_modis_qm8m5_qm8m1.png')
   do_one_fig(modis3_5, modis3_1, 
              isClear, isCloudy, isAerosolMarine, isCloudyAerosolMixMarine, 
           x_lab = 'modis_3/modis5 (0.47/1.24)', y_lab= 'modis3/modis1 (0.47/0.65)',
           filename ='misc_png/aerosol_log_marine_modis_qm3m5_qm3m1.png')
   do_one_fig(modis_5, modis8_5, isClear, isCloudy, isAerosolMarine, isCloudyAerosolMixMarine, 
           x_lab = 'modis_5 (1.24)', y_lab= 'modis8/modis5 (0.41/1.24)',
           filename ='misc_png/aerosol_log_marine_modis_m5_qm8m5.png')
   do_one_fig(modis_5, modis3_5, isClear, isCloudy, isAerosolMarine, isCloudyAerosolMixMarine, 
           x_lab = 'modis_5 (1.24)', y_lab= 'modis3/modis5 (0.47/1.24)',
           filename ='misc_png/aerosol_log_marine_modis_m5_qm3m5.png')


