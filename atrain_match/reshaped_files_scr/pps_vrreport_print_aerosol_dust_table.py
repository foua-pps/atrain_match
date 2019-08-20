#!/usr/bin/env python
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
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)
import matplotlib.pyplot as plt
from my_dir import ADIR
#ROOT_DIR = ADIR + "/git/atrain_match/modis_merged_05.h5"
#ROOT_DIR = ADIR + "/DATA_MISC/reshaped_files/global_modis_14th_created20160615/Reshaped_Files/merged/"
#ROOT_DIR = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_*_AEROSOL_SMOKE/Reshaped_Files/merged/"
#files = glob(ROOT_DIR+"*05*.h5")
#ROOT_DIR = ADIR + "/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_AEROSOL_SMOKE/Reshaped_Files/merged/"
#files = files + glob(ROOT_DIR+"*11*.h5")
#ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files/"
#            "global_modis_14th_created20161108/Reshaped_Files/merged/")
#files = glob(ROOT_DIR+"*day*.h5")
ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_modis_v2018_created20181001_cmap_osiice_dust/Reshaped_Files_merged_caliop/eos2/1km/*/*/*")
files = glob(ROOT_DIR+"*.h5")
from utils.get_flag_info import get_calipso_aerosol_of_type_i

def get_pps_aerosl(caObj):
    sunz = caObj.imager.all_arrays['sunz']
    ts = caObj.imager.all_arrays['surftemp']
    t950 = caObj.imager.all_arrays['t950']
    t850 = caObj.imager.all_arrays['t850']
    t700 = caObj.imager.all_arrays['t700']
    r06 = caObj.imager.all_arrays['r06micron']
    r09 = caObj.imager.all_arrays['r09micron']
    r13 = caObj.imager.all_arrays['r13micron']
    t37 = caObj.imager.all_arrays['bt37micron']
    t11 = caObj.imager.all_arrays['bt11micron']
    t12 = caObj.imager.all_arrays['bt12micron']
    t86 = caObj.imager.all_arrays['bt86micron']
    ctype = caObj.imager.all_arrays['cloudtype']
    thr_t11t12 = caObj.imager.all_arrays['thr_t11t12']
    thr_r06 = 100*caObj.imager.all_arrays['thr_r06']

    thr_t11t12_inv = caObj.imager.all_arrays['thr_t11t12_inv']
    thr_t11ts = caObj.imager.all_arrays['thr_t11ts']
    thr_t11ts_inv = caObj.imager.all_arrays['thr_t11ts_inv']
    feature1 = t11-t12-thr_t11t12+0.25*(thr_t11t12-thr_t11t12_inv)>0
    feature2 = t11-t12-thr_t11t12_inv-0.25*(thr_t11t12-thr_t11t12_inv)<0
    feature3 = t11> ts#+thr_t11ts>0
    feature4 = t11<ts#+thr_t11ts
    feature5 = t950-t700>-100000
    feature6 = t11-t86>0.5 #add if available
    feature7 = t11-t12<0.5 #waterclouds day
    feature8 = np.logical_and(r06<35, r06<350) #waterclouds day (better)
    feature8 = np.logical_or(sunz>95, r06<35) # clean marine, waterclouds day (better)
    feature9 = np.logical_or(sunz<70, t11-t12<0.0)
    #feature10 = ts>260

     
    is_cold_dust = np.logical_and(np.logical_and(feature4,feature5), np.logical_and(
        np.logical_and(feature2,feature9), 
        np.logical_and(feature7, feature8)))
    is_cold_dust_nonday = np.logical_and(feature4, np.logical_and(feature2, np.logical_and(feature6, feature7)))
    is_warm_dust = np.logical_and(feature8, np.logical_and(feature1,feature3))
#    is_warm_dust = np.logical_and(feature3, feature3)
    dust_singal = np.logical_or(is_warm_dust, is_cold_dust)
    safety_modis = np.logical_and(feature6, r13/np.cos(np.radians(sunz))<1.5)
    return np.logical_and(dust_singal,ctype<5)#np.logical_and(dust_singal, feature6)
    #return dust_singal

def get_calipso_cloudy_and_aerosl(caObj):    
    from scipy import stats, ndimage    
    nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
    isAerosol = caObj.calipso_aerosol.all_arrays['number_layers_found']>0
    isClear = nlay==0
    isClear = np.logical_and(isClear, np.not_equal(isAerosol, True))
    isCloudy = nlay>0

    isCloudyAerosolMix = np.logical_and(isCloudy, isAerosol)
    isCloudy = np.logical_and(isCloudy, np.not_equal(isCloudyAerosolMix, True))
    isAerosol = np.logical_and(isAerosol, np.not_equal(isCloudyAerosolMix, True))
    cloud_above = np.logical_and(isCloudyAerosolMix, caObj.calipso_aerosol.all_arrays['layer_top_altitude'][:,0]<
                                 caObj.calipso.all_arrays['layer_top_altitude'][:,0])
    isCloudyAerosolMix = np.logical_and(isCloudyAerosolMix, np.not_equal(cloud_above, True))
    isDust = get_calipso_aerosol_of_type_i(caObj, atype=2)
    isDust = np.logical_and( isDust, isAerosol)
    isCleanMarine = np.logical_and(get_calipso_aerosol_of_type_i(caObj, atype=1), isAerosol)
    isAerosol =  np.logical_and(isAerosol, ~isCleanMarine)
    #isClear = np.logical_or(isClear, isCleanMarine)
    return isCloudy, isClear, isAerosol, isCloudyAerosolMix, isDust, isCleanMarine


def is_pps_aerosol(caObj ,atype=None):
   print atype 
   cf_flag =  caObj.imager.all_arrays['cloudtype_status']
   sunz =  caObj.imager.all_arrays['sunz']
   #isPPSAerosol = cf_flag>=62 #egentligen bit 5 betyder aerosol se upp för fler bitar senare!
   #isPPSAerosol = np.logical_or(isPPSAerosol,get_pps_aerosl(caObj))
   isPPSDust = caObj.imager.cma_dust
   print len(caObj.imager.cma_dust)
   print len(caObj.imager.cma_aerosolflag)

   isPPSAerosol_all = caObj.imager.cma_aerosolflag 

   isCloudy, isClear, isAerosol_ncm, isMix, isDust, isCleanMarine = get_calipso_cloudy_and_aerosl(caObj)

   
   use = caObj.imager.all_arrays['sunz']<20000
   use_e = np.logical_and(use,caObj.imager.all_arrays['latitude']>-15)
   use_e = np.logical_and(use_e,caObj.imager.all_arrays['latitude']<45)
   use_e = np.logical_and(use_e,caObj.imager.all_arrays['longitude']>-30)
   use_e = np.logical_and(use_e,caObj.imager.all_arrays['longitude']<60)
    
   use_d = np.logical_and(sunz<=70, use)
   use_n = np.logical_and(sunz>=95, use)
   use_t = np.logical_and(np.logical_and(sunz>70,sunz<95), use)
   use_de = np.logical_and(sunz<=70, use_e)
   use_ne = np.logical_and(sunz>=90, use_e)
   #use = np.logical_and(use,caObj.imager.all_arrays['surftemp']>273.15)
   print "all days"
   print len(sunz), len(sunz<90)
   print caObj.calipso_aerosol.feature_classification_flags.shape

   i=0
   for isPPSAerosol, isAerosol in zip([isPPSAerosol_all,isPPSDust,isPPSAerosol_all],
                                     [isAerosol_ncm,isDust,isCleanMarine]):
   #["non clean marine", "dust", "clean marine"])
       atype = ["non-marine", "dust", "clean-marine"][i]
       i += 1
       for use_j, name in zip([ use, use_d, use_n, use_t, use_e, use_de, use_ne],
                              ['all','day','night','twilight', 'euro', 'euroday', 'euronight']):
           print atype, name
           use_k = np.logical_and(use_j,use)
           num_a = 1.0*sum(isAerosol[use_k])
           num_cloud = 1.0*sum(isCloudy[use_k])
           num_clear = 1.0*sum(isClear[use_k])
           #select_from = np.logical_and(use_k,isCloudy)
           #the_cloudy = np.random.choice(isPPSAerosol[select_from],np.int(num_a))
           #num_of_clouds_misclassed_as_aerosol = np.sum(the_cloudy)
           #select_from = np.logical_and(use_k,isClear)
           #the_clear = np.random.choice(isPPSAerosol[select_from],np.int(num_a))
           #num_of_clear_misclassed_as_aerosol = np.sum(the_clear)
           

           num_of_clouds_misclassed_as_aerosol2 = sum(np.logical_and(isCloudy, isPPSAerosol)[use_k])*num_a/num_cloud
           num_of_clear_misclassed_as_aerosol2 = sum(np.logical_and(isClear, isPPSAerosol)[use_k])*num_a/num_clear 
           num_of_aerosol_detected = sum(np.logical_and(isAerosol, isPPSAerosol)[use_k]) 


           print (num_a,
                  num_of_clouds_misclassed_as_aerosol2, 
                  num_of_clear_misclassed_as_aerosol2, 

                  num_a-num_of_aerosol_detected,
                  #num_a-num_of_clouds_misclassed_as_aerosol, 
                  #num_a-num_of_clear_misclassed_as_aerosol, 
                  num_of_aerosol_detected,
                  #num_of_clouds_misclassed_as_aerosol, 
                  #num_of_clear_misclassed_as_aerosol, 
                  num_cloud,
                  num_clear)

           POD_aero = (num_of_aerosol_detected)/(num_a)
           #FAR_aero = (num_of_clouds_misclassed_as_aerosol + num_of_clear_misclassed_as_aerosol)*1.0/(num_of_clouds_misclassed_as_aerosol + num_of_clear_misclassed_as_aerosol + num_of_aerosol_detected)
           FAR_aero2 = (num_of_clouds_misclassed_as_aerosol2 + num_of_clear_misclassed_as_aerosol2)*1.0/(num_of_clouds_misclassed_as_aerosol2 + num_of_clear_misclassed_as_aerosol2 + num_of_aerosol_detected)
           print "POD aerosl", POD_aero
           print "FAR aerosl", FAR_aero2
           print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
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

caObj = CalipsoImagerTrackObject()
for filename in files:
    print os.path.basename(filename)
    newObj = readCaliopImagerMatchObj(filename)
    if len(newObj.imager.cma_dust) != len(newObj.imager.longitude):
        print "skipping", os.path.basename(filename)
        continue

    if caObj.imager.cma_dust is None:
        for var_name in ['number_layers_found', 
                         'layer_top_altitude', 'feature_classification_flags']:
            
            caObj.calipso.all_arrays[var_name] = newObj.calipso.all_arrays[var_name]
            caObj.calipso_aerosol.all_arrays[var_name] = newObj.calipso_aerosol.all_arrays[var_name]
        for var_name in ['cma_dust', 'cma_aerosolflag', 'sunz', 'longitude', 'latitude']:
            caObj.imager.all_arrays[var_name] = newObj.imager.all_arrays[var_name]


    else:    
        for var_name in ['number_layers_found', 
                         'layer_top_altitude', 'feature_classification_flags']:
            caObj.calipso.all_arrays[var_name] = np.concatenate([caObj.calipso.all_arrays[var_name], 
                                                                  newObj.calipso.all_arrays[var_name]])
            caObj.calipso_aerosol.all_arrays[var_name] = np.concatenate([caObj.calipso_aerosol.all_arrays[var_name],
                                                                          newObj.calipso_aerosol.all_arrays[var_name]])
        for var_name in ['cma_dust', 'cma_aerosolflag', 'sunz', 'longitude', 'latitude']:
            caObj.imager.all_arrays[var_name] = np.concatenate([caObj.imager.all_arrays[var_name],
                                                                newObj.imager.all_arrays[var_name]])


#make_optical_depth_hist(caObj)

is_pps_aerosol(caObj)



