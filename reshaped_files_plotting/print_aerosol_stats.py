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

ROOT_DIR = "/home/a001865/git/atrain_match/modis_merged_05.h5"
files = glob(ROOT_DIR)

from plot_aerosol_histograms import (get_feature_type)

def get_calipso_cloudy_and_aerosl(caObj):    
    from scipy import stats, ndimage    
    nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
    isAerosol = caObj.calipso_aerosol.all_arrays['number_layers_found']>0
    isClear = np.logical_and(meancl<0.1,
                             caObj.calipso.all_arrays['total_optical_depth_5km']<0.0)
    isClear = np.logical_and(isClear, np.not_equal(isAerosol, True))
    isCloudy = np.logical_and( 
       #only use those clouds in 1km data
       caObj.calipso.all_arrays['number_layers_found']>0, 
       caObj.calipso.all_arrays['total_optical_depth_5km']>0.2)

    isCloudyAerosolMix = np.logical_and(isCloudy, isAerosol)
    isCloudy = np.logical_and(isCloudy, np.not_equal(isCloudyAerosolMix, True))
    isAerosol = np.logical_and(isAerosol, np.not_equal(isCloudyAerosolMix, True))
   
    return isCloudy, isClear, isAerosol, isCloudyAerosolMix


def is_pps_aerosol(caObj ,atype=None):
   print atype 
   cf_flag =  caObj.avhrr.all_arrays['cloudtype_status']
   sunz =  caObj.avhrr.all_arrays['sunz']
   isPPSAerosol = cf_flag>=32 #egentligen bit 5 betyder aerosol se upp för fler bitar senare!
   isCloudy, isClear, isAerosol, isMix = get_calipso_cloudy_and_aerosl(caObj)
   use = sunz<=90
   if atype =='heavy':
      isAerosol[caObj.calipso_aerosol.all_arrays['feature_optical_depth_532'][:,0]<2.0]=False
      isMix[caObj.calipso_aerosol.all_arrays['feature_optical_depth_532'][:,0]<2.0]=False
   if atype =='dust':
      isDust = get_feature_type(caObj, aerosol_type=2)[use]
      isAerosol = np.logical_and(isAerosol, isDust)
      isMix = np.logical_and(isMix, isDust)
   if atype =='smoke':
      isSmoke = get_feature_type(caObj, aerosol_type=6)[use]
      isAerosol = np.logical_and(isAerosol, isSmoke)
      isMix = np.logical_and(isMix, isSmoke)
   if atype =='marine':
      isMarine = get_feature_type(caObj, aerosol_type=1)[use]
      isAerosol = np.logical_and(isAerosol, isMarine)
      isMix = np.logical_and(isMix, isMarine)
   num_of_clouds_misclassed_as_aerosol = sum(np.logical_and(isCloudy[use], isPPSAerosol[use]))
   num_of_clear_misclassed_as_aerosol = sum(np.logical_and(isClear[use], isPPSAerosol[use])) 
   num_of_aerosol_detected = sum(np.logical_and(isAerosol[use], isPPSAerosol[use])) 
   num_of_mix_detected = sum(np.logical_and(isMix[use], isPPSAerosol[use])) 
   print (num_of_clouds_misclassed_as_aerosol, 
          num_of_clear_misclassed_as_aerosol, 
          num_of_aerosol_detected,
          num_of_mix_detected)
   POD_aero = (num_of_aerosol_detected +  num_of_mix_detected)*1.0/(
      sum(isAerosol[use])+sum(isMix[use]))
   FAR_aero = (num_of_clouds_misclassed_as_aerosol + num_of_clear_misclassed_as_aerosol)*1.0/(
      sum(isPPSAerosol))
   print "POD aerosl", POD_aero
   print "FAR aerosl", FAR_aero

def make_optical_depth_hist(caObj):
   import matplotlib.pyplot as plt
   import numpy as np
   isCloudy, isClear, isAerosol, isMix = get_calipso_cloudy_and_aerosl(caObj)
   aerosol_optical_depth = caObj.calipso_aerosol.all_arrays['feature_optical_depth_532'][:,0]
   hist, bins = np.histogram(aerosol_optical_depth[aerosol_optical_depth>=0], bins=50)
   width = 0.7 * (bins[1] - bins[0])
   center = (bins[:-1] + bins[1:]) / 2
   plt.bar(center, hist, align='center', width=width)
   plt.show()
   hist, bins = np.histogram(aerosol_optical_depth[aerosol_optical_depth>=0.2], bins=50)
   width = 0.7 * (bins[1] - bins[0])
   center = (bins[:-1] + bins[1:]) / 2
   plt.bar(center, hist, align='center', width=width)
   plt.show()
   hist, bins = np.histogram(aerosol_optical_depth[aerosol_optical_depth>=0.5], bins=50)
   width = 0.7 * (bins[1] - bins[0])
   center = (bins[:-1] + bins[1:]) / 2
   plt.bar(center, hist, align='center', width=width)
   plt.show()

caObj = CalipsoAvhrrTrackObject()
for filename in files:
   print os.path.basename(filename)
   newObj = readCaliopAvhrrMatchObj(filename)
   caObj= caObj + newObj

make_optical_depth_hist(caObj)

is_pps_aerosol(caObj, atype='all')
is_pps_aerosol(caObj, atype='heavy')
is_pps_aerosol(caObj, atype='dust')
is_pps_aerosol(caObj, atype='smoke')
is_pps_aerosol(caObj, atype='marine')


