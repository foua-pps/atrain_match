#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Nina.Hakansson <a001865@c14526.ad.smhi.se>

"""Read all matched data and make some plotting of thresholds
"""

from glob import glob
import os.path
import os
import numpy as np
from scipy import histogram
import re
OFFSET_DIR=os.environ.get('SM_CFG_DIR')
OFFSET_FILE=os.environ.get('SM_THROFFSETS_NAME')   
#OFFSET_FILE="threshold_offsets_gac.cfg"
OFFESET_FILE = os.path.join(OFFSET_DIR,OFFSET_FILE)

isNPP = True
isGAC_v2012 = False
RunWithOutMargins=True
RunEdited= True
PLOT_DIR = '/local_disk/nina_pps/threshold_plots/'
#PLOT_DIR = '/home/rustan/Nina/threshold_plots/NEW2/'
isACPGv2012=False
if isNPP and isACPGv2012:
     isGAC=False
     ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20130911_2012_test2/Reshaped_Files/npp/1km/"

     SATELLITE='npp_2012'
     files = glob(ROOT_DIR + "/????/??/*/*h5")
elif isNPP and  not isACPGv2012:
    isACPGv2012=False
    isGAC=False
    ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20130908_TEST1B_2/Reshaped_Files_all/npp/1km/"

    SATELLITE='npp'
    if RunEdited:
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM/Reshaped_Files/npp/1km/"
        #FewercasesROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_3_removeall_negativtest_correct_some_bugs/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_4_keep_low_quality_snow_abit_more_remove_bad_test_from_night_inv_sheme/Reshaped_Files/npp/1km/"
        
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_6_nightlandinv_landnightmount_2/Reshaped_Files/npp/1km/"
        #ROOT_DIR ="/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_6_nightlandinv_landnightmount_and_night_ice_2/Reshaped_Files/npp/1km/"
        #ROOT_DIR ="/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_6_nightlandinv_landnightmount_and_night_ice_2/Reshaped_Files/npp/1km/"
	#ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_6_nightlandinv_landnightmount_and_night_ice_and_twilight_sea_ice_and_no_ice_and_night_coast_inv/Reshaped_Files/npp/1km/"
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_8_shemes_10/Reshaped_Files/npp/1km/"
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_12_and_half_shemes/Reshaped_Files/npp/1km/"
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_12_and_half_shemes_18v2014coldtestsadded/Reshaped_Files/npp/1km/"
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_12_and_half_shemes_18v2014coldtestsadded_sunz_thr_on_snow_tests_texture_and_ice_thr/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_12_and_half_shemes_18v2014coldtestsadded_sunz_thr_on_snow_tests_texture_and_ice_thr_all/Reshaped_Files/npp/1km/"


    	#ROOT_DIR ="Reshaped_Files_2014/Reshaped_Files/npp6sh85/1km/"
        SATELLITE = 'npp_cm_edited'
    files = glob(ROOT_DIR + "/????/??/*/*h5")

elif isGAC_v2012:
    isGAC=True
    isACPGv2012=True
    ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_GAC_2012_el2010/Reshaped_Files/noaa18/5km/"
    files = glob(ROOT_DIR + "/????/??/cea*/*noaa18*h5")
    SATELLITE='avhrr18'
elif not isNPP and not isACPGv2012:
    isGAC=True
    isACPGv2012=False
    ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_GAC_2014/Reshaped_Files/noaa18/5km/"
    #ROOT_DIR ="Reshaped_Files_2014/Reshaped_Files/noaa18/5km/"
    files = glob(ROOT_DIR + "/????/??/cea*/*noaa18*h5")
    SATELLITE='avhrr18_v2014'

print OFFESET_FILE 


offset_file = open(OFFESET_FILE, 'r')  
re_offset  = re.compile(r"^SM_ACMG_(\w+)[:=]\s*(-*\d+\.*\d*)\D")
OFFSETS={}
OFFSETS['MAX_SUNZEN_DAY'] = 80.0
OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT'] = 86.0 #sm_acmg_sunzen_day_max
for line in offset_file:
    match_offset = re_offset.match(line)
    if match_offset:
        print "read_offset %s=%s"%(match_offset.group(1),match_offset.group(2))
        OFFSETS[match_offset.group(1)] = 1.0*float(match_offset.group(2))



from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from cloudsat_calipso_avhrr_statistics import (
    get_day_night_twilight_info_pps2014,
    get_land_coast_sea_info_pps2014,
    get_ice_info_pps2014,
    get_day_night_twilight_info_pps2012,    
    get_land_coast_sea_info_pps2012,
    get_ice_info_pps2012,
    get_inversion_info_pps2014,
    get_sunglint_info_pps2014,
    get_mountin_info_pps2014)
caObj = CalipsoAvhrrTrackObject()
for filename in files:
    print os.path.basename(filename)
    caObj = caObj + readCaliopAvhrrMatchObj(filename)

isClear = caObj.calipso.all_arrays['number_of_layers_found'] == 0
isCloudy = caObj.calipso.all_arrays['number_of_layers_found'] >0
isSingleLayerCloud = caObj.calipso.all_arrays['number_of_layers_found'] == 1
SeesThrough = caObj.calipso.all_arrays['lidar_surface_elevation'][0] >-999
isHigh = caObj.calipso.all_arrays['cloud_top_profile'][0] >5.0
if isGAC:
    isOpticallyThin = caObj.calipso.all_arrays['optical_depth'][0]<1.0
    notVeryThinTop = caObj.calipso.all_arrays['optical_depth'][0]>0.2
else:
    isOpticallyThin = caObj.calipso.all_arrays['optical_depth_top_layer5km']<5.0
    notVeryThinTop = caObj.calipso.all_arrays['optical_depth_top_layer5km']>0.2
#isTooThin = caObj.calipso.all_arrays['total_optical_depth5km']<0.2
isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                             caObj.avhrr.all_arrays['cloudtype']<21)
isClearPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                            caObj.avhrr.all_arrays['cloudtype']<5)
isPPSCloudyOrClear = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                    caObj.avhrr.all_arrays['cloudtype']<21)


thr_t11t12 = caObj.avhrr.all_arrays['thr_t11t12']
thr_t11t37 = caObj.avhrr.all_arrays['thr_t11t37']
thr_t37t12 = caObj.avhrr.all_arrays['thr_t37t12']
thr_t11ts = caObj.avhrr.all_arrays['thr_t11ts']
thr_t85t11 = caObj.avhrr.all_arrays['thr_t85t11']
thr_t11ts_inv = caObj.avhrr.all_arrays['thr_t11ts_inv']
thr_t11t12_inv = caObj.avhrr.all_arrays['thr_t11t12_inv']
thr_t37t12_inv = caObj.avhrr.all_arrays['thr_t37t12_inv']
thr_t11t37_inv = caObj.avhrr.all_arrays['thr_t11t37_inv']
thr_t85t11_inv = caObj.avhrr.all_arrays['thr_t85t11_inv']
thr_r06 = caObj.avhrr.all_arrays['thr_r06']
surftemp = caObj.avhrr.all_arrays['surftemp']
latitude = caObj.avhrr.all_arrays['latitude']
ciwv = caObj.avhrr.all_arrays['ciwv']
t11 = caObj.avhrr.all_arrays['bt11micron']
t85 = caObj.avhrr.all_arrays['bt86micron']
t12 = caObj.avhrr.all_arrays['bt12micron']
t37 = caObj.avhrr.all_arrays['bt37micron']
r06 = caObj.avhrr.all_arrays['r06micron']
r09 = caObj.avhrr.all_arrays['r09micron']
r16 = caObj.avhrr.all_arrays['r16micron']
t11text = caObj.avhrr.all_arrays['text_t11']
r06text = caObj.avhrr.all_arrays['text_r06']
sunz = caObj.avhrr.all_arrays['sunz']
t37t12text = caObj.avhrr.all_arrays['text_t37t12']
t37text = caObj.avhrr.all_arrays['text_t37']

cloudtype_conditions = caObj.avhrr.all_arrays['cloudtype_conditions']
cloudtype_status = caObj.avhrr.all_arrays['cloudtype_status']
cloudtype = caObj.avhrr.all_arrays['cloudtype']
cloudtype_qflag = caObj.avhrr.all_arrays['cloudtype_qflag']
if isACPGv2012:
    i_flags = get_day_night_twilight_info_pps2012(cloudtype_qflag)
    (no_qflag, isPPSNight, isPPSTwilight, isPPSDay, all_dnt_flag) =  i_flags
    lt_flag = get_land_coast_sea_info_pps2012(cloudtype_qflag)
    (no_qflag, isPPSLand, isPPSSea, isPPSCoast, all_lsc_flag) = lt_flag
    isPPSIce = get_ice_info_pps2012(cloudtype_qflag)
    isPPSInversion = False#= get_inversion_info_pps2012(cloudtype_qflag)
    isPPSSunglint = False#= get_sunglint_info_pps2012(cloudtype_qflag)
    isPPSMountain = False#= get_mountin_info_pps2012(cloudtype_qflag)
else:
    i_flags = get_day_night_twilight_info_pps2014(cloudtype_conditions)
    (no_qflag, isPPSNight, isPPSTwilight, isPPSDay, all_dnt_flag) = i_flags
    lt_flag =get_land_coast_sea_info_pps2014(cloudtype_conditions)
    (no_qflag, isPPSLand, isPPSSea, isPPSCoast, all_lsc_flag) = lt_flag 
    isPPSIce = get_ice_info_pps2014(cloudtype_status)
    isPPSInversion = get_inversion_info_pps2014(cloudtype_status)
    isPPSSunglint = get_sunglint_info_pps2014(cloudtype_conditions)
    isPPSMountain = get_mountin_info_pps2014(cloudtype_conditions)


isClear = np.logical_and(isClear,isPPSCloudyOrClear)
isCloudy = np.logical_and(isCloudy,isPPSCloudyOrClear)
isThinCirrus = np.logical_and(isCloudy, np.logical_and(
        np.logical_and(SeesThrough,isHigh),
        np.logical_and(isOpticallyThin, isSingleLayerCloud)))
isThinCirrus = np.logical_and(isThinCirrus, notVeryThinTop)
#isCloudy = isThinCirrus
# nsidc 
# 255 = Coast
# 1-100 = Sea ice concentration %
# 101 = Permamnent ice
# 0 = Ice free
isCalipsoIce = np.logical_and(caObj.calipso.all_arrays['nsidc'] >= 1,
                       caObj.calipso.all_arrays['nsidc'] <= 101)

isPPSSeaNotIce = np.logical_and(isPPSSea,np.equal(isPPSIce,False))
isSea = np.equal(caObj.calipso.all_arrays['igbp'], 17)
isSea = isSea#np.logical_and(isSea,isPPSSea)
isWater = np.logical_and(isPPSSea, isPPSSeaNotIce)
#isIce = np.logical_and(isCalipsoIce, isPPSIce)
isIce = np.logical_and(isPPSSea, isPPSIce)
isSunglintDay = np.logical_and(np.logical_and(isPPSDay, isPPSSunglint), isPPSSea)
isSunglintTwilight = np.logical_and(np.logical_and(isPPSTwilight, isPPSSunglint), isPPSSea)
isWaterNight = np.logical_and(isWater,isPPSNight)
isWaterDay = np.logical_and(np.logical_and(isWater,isPPSDay),np.equal(isSunglintDay,False))
isWaterTwilight = np.logical_and(np.logical_and(isWater,isPPSTwilight),np.equal(isSunglintTwilight,False))
isIceNight = np.logical_and(isIce,isPPSNight)
isIceDay = np.logical_and(isIce,isPPSDay)
isIceTwilight = np.logical_and(isIce,isPPSTwilight)
isLandNight = np.logical_and(np.logical_and(isPPSLand,isPPSNight),
                             np.logical_and(np.equal(isPPSInversion,False),
                                            np.equal(isPPSMountain,False)))
isLandNightMount = np.logical_and(np.logical_or(isPPSLand,isPPSCoast),
                                  np.logical_and(isPPSNight, isPPSMountain))
isLandNightInv = np.logical_and(np.logical_and(isPPSLand,isPPSNight),
                                  np.logical_and(isPPSInversion, np.equal(isPPSMountain,False)))
isCoastNightInv = np.logical_and(np.logical_and(isPPSCoast,isPPSNight),
                                  np.logical_and(isPPSInversion, np.equal(isPPSMountain,False)))
isCoastNight = np.logical_and(np.logical_and(isPPSCoast,isPPSNight),
                                  np.logical_and(np.equal(isPPSInversion,False), 
                                                 np.equal(isPPSMountain,False)))

isLandDayMount = np.logical_and(np.logical_and(isPPSLand,isPPSDay),
                                isPPSMountain)
isLandDay = np.logical_and(np.logical_and(isPPSLand,isPPSDay),
                                np.equal(isPPSMountain,False))
isCoastDayMount = np.logical_and(np.logical_and(isPPSCoast,isPPSDay),
                                isPPSMountain)
isCoastDay = np.logical_and(np.logical_and(isPPSCoast,isPPSDay),
                                np.equal(isPPSMountain,False))

isLandTwilight = np.logical_and(np.logical_and(isPPSLand,isPPSTwilight),
                             np.logical_and(np.equal(isPPSInversion,False),
                                            np.equal(isPPSMountain,False)))
isLandTwilightMount = np.logical_and(np.logical_or(isPPSLand,isPPSCoast),
                                  np.logical_and(isPPSTwilight, isPPSMountain))
isLandTwilightInv = np.logical_and(np.logical_and(isPPSLand,isPPSTwilight),
                                  np.logical_and(isPPSInversion, np.equal(isPPSMountain,False)))
isCoastTwilightInv = np.logical_and(np.logical_and(isPPSCoast,isPPSTwilight),
                                  np.logical_and(isPPSInversion, np.equal(isPPSMountain,False)))
isCoastTwilight = np.logical_and(np.logical_and(isPPSCoast,isPPSTwilight),
                                  np.logical_and(np.equal(isPPSInversion,False), 
                                                 np.equal(isPPSMountain,False)))
#isWater =  isPPSSeaNotIce
#isIce =  np.logical_and(isPPSSea,isPPSIce)
isClearWater = np.logical_and(isWater, isClear)
isCloudyWater = np.logical_and(isWater, isCloudy)
isClearIce = np.logical_and(isIce, isClear)
isCloudyIce = np.logical_and(isIce, isCloudy)
isClearLand =   np.logical_and(isClear, isPPSLand)
isCloudyLand =   np.logical_and(isCloudy, isPPSLand)
isClearCoast =   np.logical_and(isClear, isPPSCoast)
isCloudyCoast =   np.logical_and(isCloudy, isPPSCoast)
isClearIceNight = np.logical_and(isClear, isIceNight)
isCloudyIceNight = np.logical_and(isCloudy, isIceNight)
isCloudyCirrusIceNight = np.logical_and(isThinCirrus, isIceNight)
isClearWaterNight = np.logical_and(isClear, isWaterNight)
isCloudyWaterNight = np.logical_and(isCloudy, isWaterNight)
isCloudyCirrusWaterNight = np.logical_and(isThinCirrus, isWaterNight)
isClearLandNight =   np.logical_and(isClear, isLandNight)
isCloudyLandNight =  np.logical_and(isCloudy, isLandNight)
isClearLandNightMount =  np.logical_and(isClear, isLandNightMount)
isCloudyLandNightMount =  np.logical_and(isCloudy, isLandNightMount)
isClearLandNightInv =  np.logical_and(isClear,isLandNightInv) 
isCloudyLandNightInv = np.logical_and(isCloudy, isLandNightInv)
isClearCoastNightInv =  np.logical_and(isClear,isCoastNightInv) 
isCloudyCoastNightInv = np.logical_and(isCloudy, isCoastNightInv)
isClearCoastNight =  np.logical_and(isClear,isCoastNight) 
isCloudyCoastNight = np.logical_and(isCloudy, isCoastNight)
isClearIceDay = np.logical_and(isClear, isIceDay)
isCloudyIceDay = np.logical_and(isCloudy, isIceDay)
isClearSunglintDay = np.logical_and(isClear, isSunglintDay)
isCloudySunglintDay = np.logical_and(isCloudy, isSunglintDay)
isClearWaterDay = np.logical_and(isClear, isWaterDay)
isCloudyWaterDay = np.logical_and(isCloudy, isWaterDay)
isCloudyCirrusWaterDay = np.logical_and(isThinCirrus, isWaterDay)
isClearLandDay = np.logical_and(isClear, isLandDay)
isCloudyLandDay = np.logical_and(isCloudy, isLandDay)
isClearLandDayMount = np.logical_and(isClear, isLandDayMount)
isCloudyLandDayMount = np.logical_and(isCloudy, isLandDayMount)
isClearCoastDay = np.logical_and(isClear, isCoastDay)
isCloudyCoastDay = np.logical_and(isCloudy, isCoastDay)
isClearCoastDayMount = np.logical_and(isClear, isCoastDayMount)
isCloudyCoastDayMount = np.logical_and(isCloudy, isCoastDayMount)

isClearSunglintTwilight = np.logical_and(isClear, isSunglintTwilight)
isCloudySunglintTwilight = np.logical_and(isCloudy, isSunglintTwilight)
isClearIceTwilight = np.logical_and(isClear, isIceTwilight)
isCloudyIceTwilight = np.logical_and(isCloudy, isIceTwilight)
isClearWaterTwilight = np.logical_and(isClear, isWaterTwilight)
isCloudyWaterTwilight = np.logical_and(isCloudy, isWaterTwilight)
isClearLandTwilight =   np.logical_and(isClear, isLandTwilight)
isCloudyLandTwilight =  np.logical_and(isCloudy, isLandTwilight)
isClearLandTwilightMount =  np.logical_and(isClear, isLandTwilightMount)
isCloudyLandTwilightMount =  np.logical_and(isCloudy, isLandTwilightMount)
isClearLandTwilightInv =  np.logical_and(isClear,isLandTwilightInv) 
isCloudyLandTwilightInv = np.logical_and(isCloudy, isLandTwilightInv)
isClearCoastTwilightInv =  np.logical_and(isClear,isCoastTwilightInv) 
isCloudyCoastTwilightInv = np.logical_and(isCloudy, isCoastTwilightInv)
isClearCoastTwilight =  np.logical_and(isClear,isCoastTwilight) 
isCloudyCoastTwilight = np.logical_and(isCloudy, isCoastTwilight)
isCloudtypeSnowOrIce = np.logical_or(np.equal(cloudtype,4),
                                     np.equal(cloudtype,3))
nodata = np.logical_or(t12<= -9, t37<= -9)
t37t12 = np.ma.array(t37-t12, mask = nodata)
nodata = np.logical_or(t11<= -9, t37<= -9)
t11t37 = np.ma.array(t11-t37, mask = nodata)
try:
    nodata = np.logical_or(t11<= -9, t85<= -9)
    t85t11 = np.ma.array(t85-t11, mask = nodata)
except:
    t85t11=None
nodata = np.logical_or(t11<= -9, surftemp<= -9)
t11ts = np.ma.array(t11-surftemp, mask = nodata)
nodata = t11<= -9
t11 = np.ma.array(t11, mask = nodata)
nodata = np.logical_or(t11<= -9, t12<= -9)
t11t12 = np.ma.array(t11-t12, mask = nodata)
nodata = t11text<= -9
t11text = np.ma.array(t11text, mask = nodata)
nodata = r06text<= -9
r06text = np.ma.array(r06text, mask = nodata)
nodata = t37t12text<= -9
t37t12text = np.ma.array(t37t12text, mask = nodata)
nodata = t37text<= -9
t37text = np.ma.array(t37text, mask = nodata)
nodata = r06<= -9
r06 = np.ma.array(r06, mask = nodata)
nodata = r09<= -9
r09 = np.ma.array(r09, mask = nodata)
nodata = r16<= -9
r16 = np.ma.array(r16, mask = nodata)
pseudor06_not_scaled=r06
try:
    r16=r16/np.cos(np.radians(sunz))
except:
    r16=None
r06=r06/np.cos(np.radians(sunz))
r09=r09/np.cos(np.radians(sunz))
nodata = np.logical_or(r09<=-9,np.logical_or(r06<=-0.9, r06==0))
try:
    qr09r06 = np.ma.array(r09/r06,  mask = nodata)
except:
    qr09r06 = None
nodata = np.logical_or(r16<=-9,np.logical_or(r06<=-0.9, r06==0))
try:
    qr16r06 = np.ma.array(r16/r06,  mask = nodata)
except:
    qr16r06 = None
nodata = thr_r06<= -9
thr_r06 = np.ma.array(thr_r06, mask = nodata)
cu=3.5/4.0
cl=0.5/4.0
cu=3/4.0
cl=1/4.0
#cu=0.5
#cl=0.5
t37_t12_minus_threshold = t37t12 - cu*thr_t37t12 - cl*thr_t37t12_inv
t37_t12_minus_threshold_inv = t37t12 - cu*thr_t37t12_inv - cl*thr_t37t12
t11_t37_minus_threshold = t11t37 - cu*thr_t11t37 - cl*thr_t11t37_inv
t11_t37_minus_threshold_inv = t11t37 - cu*thr_t11t37_inv - cl*thr_t11t37
t11_ts_minus_threshold = t11ts - cu*thr_t11ts - cl*thr_t11ts_inv
t11_ts_minus_threshold_inv = t11ts - cu*thr_t11ts_inv - cl*thr_t11ts
t11_t12_minus_threshold = t11t12 - cu*thr_t11t12 - cl*thr_t11t12_inv
t11_t12_minus_threshold_inv = t11t12 - cu*thr_t11t12_inv - cl*thr_t11t12
try:
    t85_t11_minus_threshold = t85t11 - thr_t85t11 
except:
    t85_t11_minus_threshold = None
t37_t12_minus_threshold = t37t12 - thr_t37t12 
t37_t12_minus_threshold_inv = t37t12 - thr_t37t12_inv 
t11_t37_minus_threshold = t11t37 - thr_t11t37 
t11_t37_minus_threshold_inv = t11t37 - thr_t11t37_inv 
t11_ts_minus_threshold = t11ts - thr_t11ts 
t11_ts_minus_threshold_inv = t11ts - thr_t11ts_inv 
t11_t12_minus_threshold = t11t12 - thr_t11t12 
t11_t12_minus_threshold_inv = t11t12 - thr_t11t12_inv 
#----------------------------------------------
#Plotting
#----------------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def keep_combined_ok(TestOk,TestOkAll):
    if hasattr(TestOk,'mask') and hasattr(TestOkAll,'mask'):               
        return np.logical_or(TestOk.data, TestOkAll.data)
    elif hasattr(TestOkAll,'mask'):               
        return np.logical_or(TestOk, TestOkAll.data)
    elif hasattr(TestOk,'mask'):
        return np.logical_or(TestOk.data, TestOkAll)
    else:
        return np.logical_or(TestOk, TestOkAll)

def plot_inner(args,isCloudy, isClear, TestOk, THRESHOLD):
    isClearPPS = args['isClearPPS']  
    isCloudyPPS = args['isCloudyPPS']    
    xvector = args['xvector']
    yvector = args['yvector']
    fig = plt.figure(figsize = (9,8))
    ax = fig.add_subplot(111)
    isClearPPSCloudy = np.logical_and(isClear,isCloudyPPS)
    isCloudyPPSClear = np.logical_and(isCloudy,isClearPPS)
    isDetectedByThisTest = np.logical_and(np.logical_and(isCloudy,isCloudyPPS),
                                          TestOk)
    isMissclassifiedByThisTest = np.logical_and(isClear,TestOk)
    isMissclassifiedByThisTestWhereOk = np.logical_and(isClearPPS,isMissclassifiedByThisTest)
    isNotDetectedByThisTest = np.logical_and(np.logical_and(isCloudy, 
                                                             TestOk),
                                              isClearPPS)
    #pixels inlcuded in plot and ok test and pps cloudy.
    isPPSCloudyForTestOk = np.logical_and(TestOk,
                                          np.logical_or(isCloudy,isClear))
    POD_cloudy=np.divide(len(xvector[isDetectedByThisTest==True])*1.0,
                         len(xvector[isCloudy==True]))
    FAR_cloudy=np.divide(len(xvector[isMissclassifiedByThisTest==True])*1.0,
                         len(xvector[isPPSCloudyForTestOk==True]))
    POD_cloudy_missed = np.divide(len(xvector[isNotDetectedByThisTest==True]),
                                  1.0*len(xvector[isCloudy==True]))
    if len(xvector[isNotDetectedByThisTest == True]) !=0:
        print "warning missed %d"%(len(xvector[isNotDetectedByThisTest == True]))
    POD_FAR_INFO =  (
        "Number of cloudy pixels considered %d \n"%(
            len(xvector[isCloudy==True]))+
        "Number of clear pixels considered %d \n"%(
            len(xvector[isClear==True]))+
        "Number of clouds detected by this test %d \n"%(
            len(xvector[isDetectedByThisTest==True]))+
        "Number of clouds misclassified by this test %d \n"%(
            len(xvector[isMissclassifiedByThisTest==True]))+
        "Number of clouds misclassified by this and clear in pps %d \n"%(
            len(xvector[isMissclassifiedByThisTestWhereOk==True]))+
        "Number of cloudy pixels could be detected by this test %d \n"%(
            len(xvector[isNotDetectedByThisTest==True]))+
        "POD cloudy by this test %2.1f\n"%(POD_cloudy*100)+
        "FAR cloudy %2.1f\n"%(FAR_cloudy*100)+
        "POD cloudy missed %2.1f"%(POD_cloudy_missed*100))
    print "-Test stat---------"
    print args['title']
    print POD_FAR_INFO
    print "-------------------"
    clfree_percentiles = np.percentile(
        yvector[ np.logical_and(isClear,yvector.mask==False)==True],
        [01, 10, 25, 50, 75, 90, 99])
    #cloudy_percentiles = np.percentile(
    #    yvector[np.logical_and(isCloudy,yvector.mask==False)==True],
    #    [01, 10, 25, 50, 75, 90, 99])isCloudyPPSCorrectAll
    plt.plot(xvector[isCloudy == True], 
             yvector[isCloudy == True], 'r.')
    plt.plot(xvector[isCloudyPPSClear == True], 
             yvector[isCloudyPPSClear == True], 'y.')
    plt.plot(xvector[isClear == True], 
             yvector[isClear == True], 'bx')
    plt.plot(xvector[isDetectedByThisTest == True], 
             yvector[isDetectedByThisTest == True], 'm.')
    plt.plot(xvector[isMissclassifiedByThisTest == True], 
             yvector[isMissclassifiedByThisTest == True], 'cx')
    plt.plot(xvector[isNotDetectedByThisTest == True], 
             yvector[isNotDetectedByThisTest == True], 'g.')
    plt.axhline(y = THRESHOLD, color = 'k',linewidth = 2)
    plt.axhline(y = clfree_percentiles[0], color = 'b',ls = '--')
    plt.axhline(y = clfree_percentiles[1], color = 'b',ls = '--')
    plt.axhline(y = clfree_percentiles[3], color = 'b',linewidth = 2)
    plt.axhline(y = clfree_percentiles[5], color = 'b',ls = '--')
    plt.axhline(y = clfree_percentiles[6], color = 'b',ls = '--')
    ax.set_ylabel(args['ylable'])
    ax.set_xlabel(args['xlable'])
    limits = ax.axis()
    ax.set_ylim(1.3*limits[2]-0.3*limits[3],limits[3])
    plt.text(0.9*limits[0]+0.1*limits[1],1.25*limits[2]-0.25*limits[3], 
             POD_FAR_INFO, backgroundcolor='w')
    fig_title= (" %s \n"%(args['title'])+
                "blue-cyan: clear (calipso), red-manga: cloudy (clipso)\n"+
                "manga: detected cloudy by this test, " +
                "cyan: missclassifed by this test)")
    ax.set_title(fig_title)
    #filename=args['title']+"_"+SATELLITE+'.png'
    #plt.savefig('/local_disk/nina_pps/'+ filename)
    #plt.show()

def plot_test(args,isCloudy, isClear, TestOk, THRESHOLD, show=False):
    plot_inner(args,isCloudy, isClear, TestOk, THRESHOLD)
    filename=args['title']+"_"+SATELLITE+'.png'
    plt.savefig(PLOT_DIR + filename)
    if show:
        plt.show()
    plt.close()

def plot_test_2_lim(args,isCloudy, isClear, TestOk, 
                    THRESHOLD1, THRESHOLD2, show=False):
    plot_inner(args,isCloudy, isClear, TestOk, THRESHOLD1)
    #Add some vertical lines with threshold and clear percentiles.
    xvector = args['xvector']
    clfree_percentiles_vert = np.percentile(
        np.array(xvector[ np.logical_and(isClear, xvector.mask==False)==True]),
        [01, 10, 25, 50, 75, 90, 99])
    plt.axvline(x = clfree_percentiles_vert[0], color = 'b',ls = '--')
    plt.axvline(x = clfree_percentiles_vert[1], color = 'b',ls = '--')
    plt.axvline(x = clfree_percentiles_vert[3], color = 'b',linewidth = 2)
    plt.axvline(x = clfree_percentiles_vert[5], color = 'b',ls = '--')
    plt.axvline(x = clfree_percentiles_vert[6], color = 'b',ls = '--')
    plt.axvline(x = THRESHOLD2, color = 'k',linewidth = 2)
    filename=args['title']+"_"+SATELLITE+'.png'
    plt.savefig(PLOT_DIR + filename) 
    if show:
        plt.show()
    plt.close()

def print_stats_snow(cloudtype, isCloudyStats, isClearStats, 
                     TestOk):
    isClearPPSCloudy = np.logical_and(isClearStats,isCloudyPPS)
    isCloudyPPSClear = np.logical_and(isCloudyStats,isClearPPS)
    isDetectedByThisTest = np.logical_and(np.logical_and(isClearStats,isClearPPS),
                                          TestOk)
    isMissclassifiedByThisTest = np.logical_and(isCloudyStats,TestOk)
    isNotDetectedByThisTest = np.logical_and(np.logical_and(isClearStats, 
                                                            TestOk),
                                             isCloudyPPS)
    #pixels inlcuded in plot and ok test and pps cloudy.
    isPPSClearForTestOk = np.logical_and(TestOk,
                                          np.logical_or(isCloudyStats,isClearStats))
    POD_clear=np.divide(len(cloudtype[isDetectedByThisTest==True])*1.0,
                         len(cloudtype[isClearStats==True]))
    FAR_clear=np.divide(len(cloudtype[isMissclassifiedByThisTest==True])*1.0,
                         len(cloudtype[isPPSClearForTestOk==True]))
    POD_clear_missed = np.divide(len(cloudtype[isNotDetectedByThisTest==True]),
                                  1.0*len(cloudtype[isClearStats==True]))
    PART_cloudy_missclassed = np.divide(len(cloudtype[isMissclassifiedByThisTest==True]),
                                        1.0*len(cloudtype[isCloudyStats==True]))
    POD_FAR_INFO =  (
        "Number of clear pixels considered %d \n"%(
            len(cloudtype[isClearStats==True]))+
        "Number of clear detected by this test %d \n"%(
            len(cloudtype[isDetectedByThisTest==True]))+
        "Number of cloudy misclassified by this test %d \n"%(
            len(cloudtype[isMissclassifiedByThisTest==True]))+
        "Number of cloudy pixels could be detected by this test %d \n"%(
            len(cloudtype[isNotDetectedByThisTest==True]))+
        "POD clear by this test %2.1f\n"%(POD_clear*100)+
        "FAR clear %2.1f\n"%(FAR_clear*100)+
        "PART cloudy mis classed %2.1f"%(PART_cloudy_missclassed*100))
    print POD_FAR_INFO 
    print "*********************************"

def print_stats(cloudtype, isCloudyStats, isClearStats, 
                isCloudyPPS, TestOkAll=None):
    print "*********************************"
    print "*** stats so far ****"
    num_cloudy =len(cloudtype[isCloudyStats==True])
    num_clear = len(cloudtype[isClearStats==True])
    num_clear_miscl = len(cloudtype[
            np.logical_and(isClearStats,isCloudyPPS)==True])
    num_cloudy_miscl = len(cloudtype[
            np.logical_and(isCloudyStats,isClearPPS)==True])
    num_cloudy_ok = len(cloudtype[
            np.logical_and(isCloudyStats,isCloudyPPS)==True])
    print "Number of cloudy pixels %d"%(num_cloudy)
    print "Number of clear pixels  %d"%(num_clear)
    print "Number of clear pixels miss classed  %d"%(num_clear_miscl)
    print "Number of cloudy pixels miss classed  %d"%(num_cloudy_miscl)
    print "POD cloudy %f"%( np.divide(num_cloudy_ok*1.0,num_cloudy))
    print "FAR cloudy %f"%(np.divide(num_clear_miscl*1.0,
                                      num_cloudy_ok+num_clear_miscl))
    print "Part clear missclassed %f"%(np.divide(num_clear_miscl*1.0,num_clear))
    if TestOkAll is not None:
        detected_so_far=np.logical_and(TestOkAll, np.logical_and(
                isCloudyStats,isCloudyPPS))
        would_have_been_detected_so_far = np.logical_and(
            TestOkAll, 
            np.logical_and(isCloudyStats,isPPSCloudyOrClear))
        would_have_been_missclassed_so_far = np.logical_and(
            TestOkAll, 
            np.logical_and(isClearStats,isPPSCloudyOrClear))
        false_so_far=np.logical_and(
            TestOkAll, 
            np.logical_and(isClearStats,isCloudyPPS))
        print "Number cloudy  detected so far %f"%(
            len(cloudtype[detected_so_far==True]))
        print "POD cloudy  detected so far %f"%(
            len(cloudtype[detected_so_far==True])*1.0/
            len(cloudtype[isCloudyStats==True]))
        print "FAR cloudy so far %f"%(np.divide(
                len(cloudtype[false_so_far==True])*1.0,
                (len(cloudtype[false_so_far==True])+
                 len(cloudtype[detected_so_far==True]))))
        print "POD cloudy  detected so far including changes%f"%(
            np.divide(len(cloudtype[would_have_been_detected_so_far==True])*1.0,
            len(cloudtype[isCloudyStats==True])))
        print "FAR cloudy so far including changes %f"%(np.divide(
                len(cloudtype[would_have_been_missclassed_so_far==True])*1.0,
                (len(cloudtype[would_have_been_detected_so_far==True])+
                 len(cloudtype[would_have_been_missclassed_so_far==True]))))
    print "*********************************"



isBadShemesCloudy=np.logical_or(isCloudyLandNightInv,np.logical_or(isCloudyLandNightMount,isCloudyLandTwilight))
isBadShemesClear=np.logical_or(isClearLandNightInv,np.logical_or(isClearLandNightMount,isClearLandTwilight))
isGoodShemesCloudy=np.logical_and(np.logical_and(isPPSCloudyOrClear,isCloudy),
                                  np.equal(isBadShemesCloudy, False))
isGoodShemesClear=np.logical_and(np.logical_and(isPPSCloudyOrClear,isClear),
                                 np.equal(isBadShemesClear, False))
clear_dict={}
cloudy_dict={}
spec_clear_dict={}
spec_cloudy_dict={}
spec_clear_dict['AllShemesCirrus'] = isClear
spec_cloudy_dict['AllShemesCirrus'] = isThinCirrus
spec_clear_dict['BadShemes'] = isBadShemesClear
spec_cloudy_dict['BadShemes'] = isBadShemesCloudy
spec_clear_dict['GoodShemes'] = isGoodShemesClear
spec_cloudy_dict['GoodShemes'] = isGoodShemesCloudy

spec_clear_dict['Ice'] = isClearIce
spec_cloudy_dict['Ice'] = isCloudyIce
spec_clear_dict['Land'] =  isClearLand
spec_cloudy_dict['Land'] = isCloudyLand
spec_cloudy_dict['Coast'] =isCloudyCoast
spec_clear_dict['Coast'] = isClearCoast
spec_clear_dict['Water'] = isClearWater 
spec_cloudy_dict['Water'] =isCloudyWater

clear_dict['IceNight'] = isClearIceNight
cloudy_dict['IceNight'] = isCloudyIceNight
clear_dict['WaterNight'] =  isClearWaterNight
cloudy_dict['WaterNight'] = isCloudyWaterNight
clear_dict['LandNight'] = isClearLandNight 
cloudy_dict['LandNight'] =isCloudyLandNight 
clear_dict['LandNightMount'] = isClearLandNightMount 
cloudy_dict['LandNightMount'] = isCloudyLandNightMount
clear_dict['LandNightInv'] =  isClearLandNightInv
cloudy_dict['LandNightInv'] = isCloudyLandNightInv
clear_dict['CoastNightInv'] = isClearCoastNightInv 
cloudy_dict['CoastNightInv'] =isCloudyCoastNightInv 
clear_dict['CoastNight'] = isClearCoastNight 
cloudy_dict['CoastNight'] =isCloudyCoastNight 
clear_dict['IceDay'] =  isClearIceDay
cloudy_dict['IceDay'] = isCloudyIceDay
clear_dict['SunglintDay'] = isClearSunglintDay 
cloudy_dict['SunglintDay'] =isCloudySunglintDay 
clear_dict['WaterDay'] = isClearWaterDay 
cloudy_dict['WaterDay'] =isCloudyWaterDay 
clear_dict['LandDay'] = isClearLandDay 
cloudy_dict['LandDay'] =isCloudyLandDay 
clear_dict['LandDayMount'] = isClearLandDayMount 
cloudy_dict['LandDayMount'] = isCloudyLandDayMount
clear_dict['CoastDay'] = isClearCoastDay 
cloudy_dict['CoastDay'] =isCloudyCoastDay 
clear_dict['CoastDayMount'] =  isClearCoastDayMount
cloudy_dict['CoastDayMount'] =isCloudyCoastDayMount
clear_dict['SunglintTwilight'] = isClearSunglintTwilight 
cloudy_dict['SunglintTwilight'] =isCloudySunglintTwilight 
clear_dict['IceTwilight'] =  isClearIceTwilight
cloudy_dict['IceTwilight'] = isCloudyIceTwilight
clear_dict['WaterTwilight'] =  isClearWaterTwilight
cloudy_dict['WaterTwilight'] = isCloudyWaterTwilight
clear_dict['LandTwilight'] =  isClearLandTwilight
cloudy_dict['LandTwilight'] = isCloudyLandTwilight
clear_dict['LandTwilightMount'] =  isClearLandTwilightMount
cloudy_dict['LandTwilightMount'] = isCloudyLandTwilightMount
clear_dict['LandTwilightInv'] = isClearLandTwilightInv 
cloudy_dict['LandTwilightInv'] = isCloudyLandTwilightInv
clear_dict['CoastTwilightInv'] =  isClearCoastTwilightInv
cloudy_dict['CoastTwilightInv'] = isCloudyCoastTwilightInv
clear_dict['CoastTwilight'] =  isClearCoastTwilight
cloudy_dict['CoastTwilight'] = isCloudyCoastTwilight

TestOkAll = None
for key in cloudy_dict.keys():
    print key
    print_stats(cloudtype, cloudy_dict[key], clear_dict[key], isCloudyPPS)
for key in spec_cloudy_dict.keys():
    print key
    print_stats(cloudtype, spec_cloudy_dict[key], spec_clear_dict[key], isCloudyPPS)

####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
dt11tsur_offset_extra = -8.0
t11t37_clfree = 0.5*(thr_t11t37 + thr_t11t37_inv)
t37t12_clfree = 0.5*(thr_t37t12 + thr_t37t12_inv)
t11t37_departure = np.abs(t11t37 - t11t37_clfree)
t37t12_departure = np.abs(t37t12 - t37t12_clfree)
t11t37_maxdep = np.abs(0.5*(thr_t11t37 - thr_t11t37_inv))
t37t12_maxdep = np.abs(0.5*(thr_t37t12 - thr_t37t12_inv))
t11t37_div = t11t37_departure/t11t37_maxdep
t37t12_div = t37t12_departure/t37t12_maxdep
T11_OFFSET=-4
offset = (T11_OFFSET + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TSUR_THRESHOLD = 0
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
       surftemp>TSUR_THRESHOLD)
TestOkAll=TestOk
plot_test(args,isCloudy,isClear,TestOk,np.min(THRESHOLD), show=True)
print_stats(cloudtype,isCloudy,isClear,isCloudyPPS,TestOkAll)

#----------------------------------
#ColdWaterCloudTest all
#----------------------------------
args = {'title': "coldwaterCloudTest_v2014_All",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']
THRESHOLD1 = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT']+ -1*OFFSETS['QUALITY_MARGIN_T11TSUR']
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = - 7.0
    THRESHOLD1 = -0
TestOk = np.logical_and(np.logical_and(surftemp>TSUR_THRESHOLD,
                        (((t11_t37_minus_threshold-THRESHOLD2)/11)**2 +((t11_ts_minus_threshold-THRESHOLD1)/20)**2)>1.0               ),
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                                       t11_ts_minus_threshold<THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudy,isClear, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudy,isClear,isCloudyPPS,TestOkAll)
#----------------------------------
#thincoldCirrusTest All
#----------------------------------
args = {'title': "thincoldCirrusTest_v2014_All",
        'xlable': 'T11-TS minus dynamic threshold',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t11_ts_minus_threshold,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}

if RunWithOutMargins:
    THRESHOLD1 = 0#OFFSETS['T11T12_OFFSET_LAND_DAY']
    THRESHOLD2 = -2.0#OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE']
TSUR_THRESHOLD = 0.0
TestOk = np.logical_and(
    np.logical_and(surftemp>TSUR_THRESHOLD,
                   (((t11_t12_minus_threshold-THRESHOLD2)/4)**2 +((t11_ts_minus_threshold-THRESHOLD1)/20)**2)>1.0               ),
    np.logical_and(t11_t12_minus_threshold>THRESHOLD1 ,
                   t11_ts_minus_threshold<THRESHOLD2))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudy, isClear, TestOk, THRESHOLD1,THRESHOLD2, show=False)
print_stats(cloudtype,isCloudy,isClear,isCloudyPPS, TestOkAll)

####################################
# NIGHT LAND 
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)
#----------------------------------
#coldcloudtest land night  
#----------------------------------
args = {'title': "coldCloudTest_All_LandNight",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_LAND_NIGHT_OPAQUE'] - OFFSETS['QUALITY_MARGIN_T11TSUR'] 
TSUR_THRESHOLD = 0
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_LAND_NIGHT_OPAQUE'] 
TestOk = t11_ts_minus_threshold<THRESHOLD
TestOkAll = TestOk
plot_test(args,isCloudyLandNight,isClearLandNight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest land night 
#limit for desert is highter
#----------------------------------
print "limit ts for desert is 7deg lower"
args = {'title': "coldwaterCloudTest_All_LandNight",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']
THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_NIGHT']+ -1*OFFSETS['QUALITY_MARGIN_T11TSUR']
TSUR_THRESHOLD = 230
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_NIGHT'] 
TestOk = np.logical_and(surftemp>TSUR_THRESHOLD,
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                                       t11_ts_minus_threshold<THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandNight,isClearLandNight, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)
#----------------------------------
#watercloudtest land night  
#----------------------------------
args = {'title': "less_generous_waterCloudTest_All_LandNight",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] + OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']
TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNight,isClearLandNight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)
#----------------------------------
args['title']= "less_generous_waterCloudTest_All_LandNight_orginal"
THRESHOLD = OFFSETS['QUALITY_MARGIN_T11T37'] 
if RunWithOutMargins:
    THRESHOLD = 0
TestOk = t11_t37_minus_threshold>THRESHOLD
plot_test(args,isCloudyLandNight,isClearLandNight,TestOk,THRESHOLD, show=False)
#----------------------------------
#ThinCirrusSecondaryTest land night 
#----------------------------------
print "only for desert can not sort out those pixels!"
args = {'title': "desert_thinCirrusSecondaryTest_All_LandNight_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T12_OFFSET_LAND_NIGHT'] + OFFSETS['QUALITY_MARGIN_T11T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_LAND_NIGHT'] 
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNight, isClearLandNight, TestOk, THRESHOLD,show=False)
print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)
#----------------------------------
#ThinCirrusPrimaryTest land night 
#----------------------------------
args = {'title': "more_generous_thinCirrusPrimaryTest_All_LandNight",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -2.0   +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -2.0 
TestOk = t37_t12_minus_threshold>THRESHOLD                      
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNight,isClearLandNight,TestOk,THRESHOLD,show=False)
print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)
args['title']= "more_generous_thinCirrusPrimaryTest_All_LandNight_original"
THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT']   +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] 
TestOk = t37_t12_minus_threshold>THRESHOLD
plot_test(args,isCloudyLandNight,isClearLandNight,TestOk,THRESHOLD,show=False)
#----------------------------------
# VERSION 2014 COLD_CLOUD_TEST
#----------------------------------
args = {'title': "coldCloudTest_v2014_All_LandNight",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
dt11tsur_offset_extra = -8
offset = (OFFSETS['T11_OFFSET_LAND_NIGHT'] + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset - OFFSETS['QUALITY_MARGIN_T11TSUR']
TSUR_THRESHOLD = 0
if RunWithOutMargins:
    THRESHOLD = offset 
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNight,isClearLandNight,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)
#----------------------------------
#HighcloudTestt85t11land land night 
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11land_All_LandNight",
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD = OFFSETS['HI_T85T11_CM_LAND']
    TestOk = t85_t11_minus_threshold>THRESHOLD
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyLandNight,isClearLandNight,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyLandNight,isClearLandNight,isCloudyPPS,TestOkAll)


####################################
# NIGHT LAND MOUNTAIN
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyLandNightMount,isClearLandNightMount,isCloudyPPS,TestOkAll)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_LandNightMount",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
TSUR_THRESHOLD = 0
offset = (-8 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSUR_THRESHOLD)
plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
#watercloudtest land night  mountain
#----------------------------------
args = {'title': "waterCloudTest_All_LandNightMount",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] + OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
TSUR_THRESHOLD=230
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD,
                        surftemp>TSUR_THRESHOLD)
TestOkAll=TestOk
plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyLandNightMount,isClearLandNightMount,isCloudyPPS,TestOkAll)
#----------------------------------
THRESHOLD = 0.0
TestOk = t11_t37_minus_threshold>THRESHOLD
args['title']= "removed_waterCloudTest_All_LandNightMount"
plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,THRESHOLD)
#----------------------------------
#ColdWaterCloudTest land night mountain
#limit for desert is highter
#----------------------------------
print "limit ts for desert is 3deg lower"
args = {'title': "coldwaterCloudTest_All_LandNightMount",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']
THRESHOLD1 = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT']+ -1*OFFSETS['QUALITY_MARGIN_T11TSUR']
TSUR_THRESHOLD = 230
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT'] 
TestOk = np.logical_and(surftemp>TSUR_THRESHOLD,
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                                       t11_ts_minus_threshold<THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandNightMount,isClearLandNightMount, TestOk, 
                THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyLandNightMount,isClearLandNightMount,isCloudyPPS,TestOkAll)
#----------------------------------
#coldcloudtest land night  mountain
#----------------------------------
args = {'title': "more_generous_not_opaque_coldCloudTest_All_LandNightMount",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT'] - OFFSETS['QUALITY_MARGIN_T11TSUR'] 
TSUR_THRESHOLD = 230
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT'] 
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,
                        surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyLandNightMount,isClearLandNightMount,isCloudyPPS,TestOkAll)
THRESHOLD = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT'] - OFFSETS['QUALITY_MARGIN_T11TSUR'] 
#----------------------------------
args['title'] = "more_generous_not_opaque_coldCloudTest_All_LandNightMount_orginal"
THRESHOLD = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT_OPAQUE'] - OFFSETS['QUALITY_MARGIN_T11TSUR'] 
TSUR_THRESHOLD = 230
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT_OPAQUE'] 
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,
                        surftemp>TSUR_THRESHOLD)
plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,THRESHOLD)
#----------------------------------
#ThinCirrusPrimaryTest land night mountain
#----------------------------------
args = {'title': "more_generous_thinCirrusPrimaryTest_All_LandNightMount",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT']  -1.0 +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -1.0
TSUR_THRESHOLD = 230
TestOk = np.logical_and(t37_t12_minus_threshold>THRESHOLD,
                        surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyLandNightMount,isClearLandNightMount,isCloudyPPS,TestOkAll)
#----------------------------------
args['title'] = "more_generous_thinCirrusPrimaryTest_All_LandNightMount_orginal_thr"
THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] 
TSUR_THRESHOLD = 230
TestOk = np.logical_and(t37_t12_minus_threshold>THRESHOLD,
                        surftemp>TSUR_THRESHOLD)
plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,THRESHOLD)
#----------------------------------
#ThinCirrusSecondaryTest land night mountain
#----------------------------------
print "only for desert can not sort out those pixels!"
args = {'title': "thinCirrusSecondaryTest_All_LandNightMount_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T12_OFFSET_LAND_NIGHT'] + OFFSETS['QUALITY_MARGIN_T11T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_LAND_NIGHT'] 
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNightMount, isClearLandNightMount, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyLandNightMount,isClearLandNightMount,isCloudyPPS,TestOkAll)
#----------------------------------
#HighcloudTestt85t11land land night mount
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11land_All_LandNightMount",
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    TSUR_THRESHOLD = 230
    THRESHOLD = OFFSETS['HI_T85T11_CM_LAND']
    TestOk = np.logical_and(t85_t11_minus_threshold>THRESHOLD,
                            surftemp>TSUR_THRESHOLD)
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyLandNightMount,isClearLandNightMount,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyLandNightMount,isClearLandNightMount,isCloudyPPS,TestOkAll)

####################################
# NIGHT COAST INVERSION
####################################
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_CoastNightInv",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-11 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TSURF_THRSHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSURF_THRSHOLD) 
plot_test(args,isCloudyCoastNightInv,isClearCoastNightInv,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
#watercloudtest coast night  inversion
#----------------------------------
args = {'title': "waterCloudTest_All_CoastNightInv",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD =OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] + OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll=TestOk
plot_test(args,isCloudyCoastNightInv,isClearCoastNightInv,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyCoastNightInv,isClearCoastNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest  coast night  inversion
#----------------------------------
print "not done if strong inversion, results might differ a bit"
args = {'title': "coldwaterCloudTest_All_CoastNightInv",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
tmp_off=np.min([OFFSETS['T11_OFFSET_SEA_NIGHT'],OFFSETS['T11_OFFSET_LAND_NIGHT']])
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.5
THRESHOLD1 = tmp_off -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = tmp_off 
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                        t11_ts_minus_threshold<THRESHOLD1)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyCoastNightInv,isClearCoastNightInv, TestOk, 
                THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyCoastNightInv,isClearCoastNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#watercloudtest coast night  inversion
#----------------------------------
args = {'title': "removed_waterCloudTest_All_CoastNightInv_without_security_offset",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD =OFFSETS['QUALITY_MARGIN_T11T37']
if RunWithOutMargins:
    THRESHOLD = 0.0 
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD,surftemp>TSUR_THRESHOLD)
plot_test(args,isCloudyCoastNightInv,isClearCoastNightInv,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyCoastNightInv,isClearCoastNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#ThinCirrusPrimaryTest coast night  inversion
#----------------------------------
args = {'title': "more_generous_thinCirrusPrimaryTest_All_CoastNightInv",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
tmp_off=np.min([OFFSETS['T37T12_OFFSET_SEA_NIGHT'],OFFSETS['T37T12_OFFSET_LAND_NIGHT']])
THRESHOLD = tmp_off -1.0 +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = tmp_off -1.0
TestOk = t37_t12_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyCoastNightInv,isClearCoastNightInv,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyCoastNightInv,isClearCoastNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
args = {'title': "more_generous_thinCirrusPrimaryTest_All_CoastNightInv_orginal",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
tmp_off=np.min([OFFSETS['T37T12_OFFSET_SEA_NIGHT'],OFFSETS['T37T12_OFFSET_LAND_NIGHT']])
THRESHOLD = tmp_off +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = tmp_off
TestOk = t37_t12_minus_threshold>THRESHOLD
plot_test(args,isCloudyCoastNightInv,isClearCoastNightInv,TestOk,THRESHOLD)
#----------------------------------
#coldcloudtest coast night  inversion
#----------------------------------
print "not done if strong inversion, results might differ a bit"
args = {'title': "coldCloudTest_All_CoastNightInv_without_tsur_threshold",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_INVERSION_WEAK'] + OFFSETS['QUALITY_MARGIN_T11TSUR'] 
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_INVERSION_WEAK'] 
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>0)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyCoastNightInv,isClearCoastNightInv,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyCoastNightInv,isClearCoastNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#HighcloudTestt85t11 coast night inv
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11land_All_CoastNightInv",
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    TSUR_THRESHOLD = 0
    THRESHOLD = OFFSETS['HI_T85T11_CM_LAND']
    TestOk = np.logical_and(t85_t11_minus_threshold>THRESHOLD,
                            surftemp>TSUR_THRESHOLD)
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyCoastNightInv,isClearCoastNightInv,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyCoastNightInv,isClearCoastNightInv,isCloudyPPS,TestOkAll)


####################################
# NIGHT LAND INVERSION
####################################
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_LandNightInv",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-11 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
#watercloudtest land night  inversion
#----------------------------------
args = {'title': "removed_waterCloudTest_All_LandNightInv",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] + OFFSETS['QUALITY_MARGIN_T11T37']#
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] 
TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyLandNightInv,isClearLandNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest land night inverison 
#----------------------------------
print "not done if strong inversion"
args = {'title': "coldwaterCloudTest_All_LandNightInv",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.5
THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_NIGHT']+ -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_NIGHT'] 
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                        t11_ts_minus_threshold<THRESHOLD1)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandNightInv,isClearLandNightInv, TestOk, 
                THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyLandNightInv,isClearLandNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#watercloudtest land night  inversion
#----------------------------------
print "this test is not done over desert"
args = {'title': "less_generous_waterCloudTest_All_LandNightInv",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] +3.0 + OFFSETS['QUALITY_MARGIN_T11T37']#
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']+3.0 
TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyLandNightInv,isClearLandNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
print "this test is not done over desert"
args['title'] = "less_generous_waterCloudTest_All_LandNightInv_orginal"
THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] + OFFSETS['QUALITY_MARGIN_T11T37']#
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']
TestOk = t11_t37_minus_threshold>THRESHOLD
plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,THRESHOLD)
#----------------------------------
#thinCirrusSecondaryTest night land inversion
#----------------------------------
print "only for desert can not sort out those pixels!"
args = {'title': "thinCirrusSecondaryTest_All_LandNightInv_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T12_OFFSET_LAND_NIGHT'] + OFFSETS['QUALITY_MARGIN_T11T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_LAND_NIGHT'] 
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNightInv, isClearLandNightInv, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyLandNightInv,isClearLandNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#ThinCirrusPrimaryTest land night inv
#----------------------------------
args = {'title': "more_generous_thinCirrusPrimaryTest_All_LandNightInv",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -1.0 +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -1.0
TestOk = t37_t12_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyLandNightInv,isClearLandNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
args['title'] = "more_generous_thinCirrusPrimaryTest_All_LandNightInv_orginal"
THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT'] +OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_NIGHT']
TestOk = t37_t12_minus_threshold>THRESHOLD
plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,THRESHOLD)
#----------------------------------
#coldcloudtest land night  inv
#----------------------------------
print "this is not done over desert and not if strong inversion"
args = {'title': "latitude_problem_coldCloudTest_All_LandNightInv_without_tsur_threshold",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_INVERSION_WEAK'] + OFFSETS['QUALITY_MARGIN_T11TSUR']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_INVERSION_WEAK'] 
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>0)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyLandNightInv,isClearLandNightInv,isCloudyPPS,TestOkAll)
#----------------------------------
#HighcloudTestt85t11land land night inv
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11land_All_LandNightInv",
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    TSUR_THRESHOLD = 0
    THRESHOLD = OFFSETS['HI_T85T11_CM_LAND']
    TestOk = np.logical_and(t85_t11_minus_threshold>THRESHOLD,
                            surftemp>TSUR_THRESHOLD)
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyLandNightInv,isClearLandNightInv,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyLandNightInv,isClearLandNightInv,isCloudyPPS,TestOkAll)


####################################
# SEA NIGHT
####################################
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "add_coldCloudTest_v2014_All_SeaNightNoIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-7 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
#t11_t37 should be inv
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
TSURF_THRSHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSURF_THRSHOLD) 
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS)
#----------------------------------
#watercloudtest sea night 
#----------------------------------
args = {'title': "waterCloudTest_All_SeaNightNoIce",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
if RunWithOutMargins:
    THRESHOLD = 0.0
TestOk = t11_t37_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD,show=False)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
#----------------------------------
#ThinCirrusPrimaryTest sea night 
#----------------------------------
args = {'title': "latitude_problem_thinCirrusPrimaryTest_All_SeaNightNoIce",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_SEA_NIGHT']+OFFSETS['QUALITY_MARGIN_T37T12']#2.0+0.3
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_SEA_NIGHT']
TestOk = t37_t12_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
#----------------------------------
args['title']= "thinCirrusPrimaryTest_ThinCirrus_SeaNightNoIce"
plot_test(args,isCloudyCirrusWaterNight,isClearWaterNight,TestOk,THRESHOLD)
#----------------------------------
args = {'title': "testThinCirrusPrimaryTest_All_SeaNightNoIce_abslat_linearthreshold",
        'xlable': 'abs Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': np.abs(latitude),
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_SEA_NIGHT']+OFFSETS['QUALITY_MARGIN_T37T12']#2.0+0.3
THRESHOLD_V=abs(latitude)*(THRESHOLD+1.0)/40-1.0
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_SEA_NIGHT']
    THRESHOLD_V=abs(latitude)*(THRESHOLD+1.0)/40-1.0
TestOk = args['yvector']>THRESHOLD_V
plot_test(args,isCloudyCirrusWaterNight,isClearWaterNight,TestOk,THRESHOLD)
#----------------------------------
args = {'title':  "testThinCirrusPrimaryTest_ThinCirrus_SeaNightNoIce_lat_modified_thr",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus 0.5*dynamic threshold',
        'xvector': latitude,
        'yvector': t37t12 -0.5*thr_t37t12,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_SEA_NIGHT']+3**OFFSETS['QUALITY_MARGIN_T37T12']#2.0+0.3
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_SEA_NIGHT']
TestOk = args['yvector']>THRESHOLD 
plot_test(args,isCloudyCirrusWaterNight, isClearWaterNight, TestOk, THRESHOLD)
#----------------------------------
#textureNightTest sea night 
#----------------------------------
args = {'title': "texturNightTest_All_SeaNightNoIce",
        'ylable': 'T37-T12text minus dynamic threshold',
        'xlable': 'T11text minus dynamic threshold',
        'yvector': t37t12text,
        'xvector': t11text,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}

THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'] + OFFSETS['QUALITY_MARGIN_T11TEXT']
THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT'] +  OFFSETS['QUALITY_MARGIN_T37T12TEXT']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'] 
    THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT'] 
if not isGAC:
    THRESHOLD1 = THRESHOLD1+0.3
TestOk = np.logical_and(t11text>THRESHOLD2, t37t12text>THRESHOLD1)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterNight, isClearWaterNight, TestOk, THRESHOLD1, THRESHOLD2, show=True)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
#----------------------------------
# Test is not perfect but better than nothing
#args['title']= "less_generous_texturNightTest_All_SeaNightNoIce_orginal"
#if isGAC==True:
#    T37T12TEXT_OFFSET_SEA_NIGHT_ORGINAL=1.7
#else:
#    T37T12TEXT_OFFSET_SEA_NIGHT_ORGINAL=0.85
#THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'] + OFFSETS['QUALITY_MARGIN_T11TEXT']
#THRESHOLD1 = T37T12TEXT_OFFSET_SEA_NIGHT_ORGINAL +  OFFSETS['QUALITY_MARGIN_T37T12TEXT']
#if RunWithOutMargins:
#    THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'] 
#    THRESHOLD1 = T37T12TEXT_OFFSET_SEA_NIGHT_ORGINAL
#TestOk = np.logical_and(t11text>THRESHOLD2, t37t12text>THRESHOLD1)
#plot_test_2_lim(args,isCloudyWaterNight, isClearWaterNight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
#----------------------------------
#coldcloudtest sea night without tsur limit
#----------------------------------
args = {'title': "without_tsur_limit_coldCloudTest_All_SeaNightNoIce_without_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] 
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE']
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest sea night 
#----------------------------------
args = {'title': "coldwaterCloudTest_All_SeaNightNoIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']
THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_NIGHT']+ -1*OFFSETS['QUALITY_MARGIN_T11TSUR']
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_NIGHT']
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                        t11_ts_minus_threshold<THRESHOLD1)
TestOkAll=keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterNight,isClearWaterNight, TestOk, 
                THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,
            isCloudyPPS,TestOkAll)
#----------------------------------
#coldcloudtest sea night 
#----------------------------------
args = {'title': "latitude_problem_coldCloudTest_All_SeaNightNoIce_with_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT'] 
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT']
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] 
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
#----------------------------------
#HighcloudTestt85t11sea sea night 
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11sea_All_SeaNightNoIce",
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD = OFFSETS['HI_T85T11_CM_SEA']
    TestOk = t85_t11_minus_threshold>THRESHOLD
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
#----------------------------------
#notexistingCirrustest sea night 
#----------------------------------
print "notexistingCirrustest sea night"
args = {'title': "non_existingT11T12test_All_SeaNightNoIce",
      'xlable': 'Latitude',
      'ylable': 'T11-T12 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t12_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.3
TestOk = t11_t12_minus_threshold>THRESHOLD
#TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD, show=True)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
print "watercloudOverWaterTest"
print "*** need warmest neighbour to plot this test***"
####################################
# DAY LAND MOUNTAIN 
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyLandDayMount,
                                                            isClearLandDayMount))
print_stats_snow(cloudtype, isCloudyLandDayMount, isClearLandDayMount, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyLandDayMount = np.logical_and(isCloudyLandDayMount, np.not_equal(isCloudtypeSnowOrIce,True))
isClearLandDayMount = np.logical_and(isClearLandDayMount, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_LandDayMount",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-9 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyLandDayMount,isClearLandDayMount,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
# newWaterBodiesTest
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND'] +OFFSETS['R06_OFFSET_LAND']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
if qr09r06 is not None:
    args = {'title': "CLEAR_newWaterBodiesTest_All_LandDayMount",
            'xlable': 'T37-T12',
            'ylable': 'qr09r06',
            'xvector': t37t12,
            'yvector': qr09r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD1 = 0.8
    THRESHOLD2 = r06_threshold
    THRESHOLD3 = 3.0
    THRESHOLD4 = 3.0
    if RunWithOutMargins:
        pass
    TestOk = np.logical_and(np.logical_and(qr09r06<THRESHOLD1,
                                           r06<THRESHOLD2),
                            np.logical_and(np.abs(t37t12)<THRESHOLD3,
                                           np.abs(t11t37)<THRESHOLD4))
    #TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudy,isClear, TestOk, 
                np.min(THRESHOLD1), np.min(THRESHOLD2), show=False)
    plot_test_2_lim(args,isCloudyLandDayMount,isClearLandDayMount, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    isCloudyLandDayMount = np.logical_and(isCloudyLandDayMount, np.equal(TestOk,False))
    isClearLandDayMount = np.logical_and(isClearLandDayMount, np.equal(TestOk,False))
    print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS)
#----------------------------------
#brightCloudTest land day  mountain
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND'] +OFFSETS['R06_OFFSET_LAND']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
args = {'title': "brightCloudTest_All_LandDayMount",
      'xlable': 'T37-T12 minus dynamic threshold',
      'ylable': 'r06 ',
      'xvector': t37_t12_minus_threshold,
      'yvector': r06,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}

THRESHOLD2 = OFFSETS['T37T12_OFFSET_LAND_DAY'] +OFFSETS['QUALITY_MARGIN_T37T12']#
THRESHOLD1 = r06_threshold  + OFFSETS['QUALITY_MARGIN_R06']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37T12_OFFSET_LAND_DAY']
    THRESHOLD1 = r06_threshold 
TestOk =   np.logical_and(t37_t12_minus_threshold>THRESHOLD2,
                          r06>THRESHOLD1)
TestOkAll = TestOk
plot_test_2_lim(args,isCloudyLandDayMount,isClearLandDayMount, TestOk, 
                np.min(THRESHOLD1), THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS,TestOkAll)
#----------------------------------
#brightCloudTestNoSunglint3A Land day mountain
#----------------------------------
if qr16r06 is not None:
    args = {'title': "notusedif37_brightCloudTestNoSunglint3A_All_LandDayMount",
            'xlable': 'r06 ',
            'ylable': 'qr16r06',
            'xvector': r06,
            'yvector': qr16r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD2 = r06_threshold  +OFFSETS['QUALITY_MARGIN_R06']
    THRESHOLD1 = (OFFSETS['R06_R16_THR_SNOW'] + 
                  1*OFFSETS['QUALITY_MARGIN_QR16R06'])
    THRESHOLD3 = 1.35 - OFFSETS['QUALITY_MARGIN_QR16R06']
    TSUR_THRESHOLD = 0.0
    if RunWithOutMargins:
        THRESHOLD2 = r06_threshold 
        THRESHOLD1 = OFFSETS['R06_R16_THR_SNOW'] 
        THRESHOLD3 = 1.35
        TestOk = np.logical_and(np.logical_and(r06>THRESHOLD2,
                                               np.logical_and(qr16r06>THRESHOLD1,
                                                              qr16r06<THRESHOLD3)),
                                t37<=-9)
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudyLandDayMount,isClearLandDayMount, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,
                isCloudyPPS,TestOkAll)
    args['title']= "dont_usealsoift37_brightCloudTestNoSunglint3A_All_LandDayMount"
    TestOk = np.logical_and(r06>THRESHOLD2,
                            np.logical_and(qr16r06>THRESHOLD1,
                                           qr16r06<THRESHOLD3))
    plot_test_2_lim(args,isCloudyLandDayMount,isClearLandDayMount, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)


#----------------------------------

#coldBrightCloudTest37 land day mount
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND'] +OFFSETS['R06_OFFSET_LAND']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))

args = {'title': "coldBrightCloudTest37_All_LandDayMount",
        'xlable': 'R06',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': r06,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = OFFSETS['TSUR_THR_SNOW_MAX'] + OFFSETS['QUALITY_MARGIN_TSUR']#
THRESHOLD2 = r06_threshold + OFFSETS['QUALITY_MARGIN_R06']
THRESHOLD3 = 8.0 + OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'] 
    TSUR_THRESHOLD = OFFSETS['TSUR_THR_SNOW_MAX']
    THRESHOLD2 = r06_threshold
    THRESHOLD3 = 8.0
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold<THRESHOLD1,
                                       surftemp>TSUR_THRESHOLD),
                        np.logical_and(r06>THRESHOLD2,
                                       t37_t12_minus_threshold>THRESHOLD3))                        
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDayMount, isClearLandDayMount, TestOk, THRESHOLD1, np.min(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS,TestOkAll)


#----------------------------------
#coldBrightCloudTest land day mountain
#----------------------------------
r06_threshold_opaque =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND_OPAQUE'] +OFFSETS['R06_OFFSET_LAND_OPAQUE']
r06_threshold_opaque = np.ma.array(
    r06_threshold_opaque, 
    mask=np.logical_or(r06.mask,thr_r06.mask))

args = {'title': "coldBrightCloudTest_All_LandDayMount",
        'xlable': 'R06',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': r06,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = 0
THRESHOLD2 = r06_threshold_opaque + OFFSETS['QUALITY_MARGIN_R06']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY'] 
    TSUR_THRESHOLD = 0
    THRESHOLD2 = r06_threshold_opaque
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold<THRESHOLD1,
                                       sunz<OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT']),
                        r06>THRESHOLD2)                        
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDayMount, isClearLandDayMount, TestOk, THRESHOLD1, np.min(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS,TestOkAll)
#----------------------------------
#thincoldCirrusTest land day mountain
#----------------------------------
args = {'title': "twice_thincoldCirrusTest_All_LandDayMount",
        'xlable': 'T11-TS minus dynamic threshold',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t11_ts_minus_threshold,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}

THRESHOLD1 = OFFSETS['T11T12_OFFSET_LAND_DAY'] + OFFSETS['QUALITY_MARGIN_T11T12']#-7.0-1.0
THRESHOLD2 = OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'] - OFFSETS['QUALITY_MARGIN_T11TSUR']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11T12_OFFSET_LAND_DAY']
    THRESHOLD2 = OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE']
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1 ,
                        t11_ts_minus_threshold<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDayMount, isClearLandDayMount, TestOk, THRESHOLD1,THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS, TestOkAll)
#----------------------------------
#ColdWaterCloudTestDay land day mountain 
#----------------------------------
args = {'title': "coldwaterCloudTestDay_All_LandDayMount",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_MOUNTAIN_DAY'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#
THRESHOLD3 = 298.15 - 1.0
THRESHOLD4 = OFFSETS['R09_R06_THR_CLOUDY_MAX']- OFFSETS['QUALITY_MARGIN_QR09R06']#
THRESHOLD5 = OFFSETS['R09_R06_THR_CLOUDY_MIN']+ OFFSETS['QUALITY_MARGIN_QR09R06']#
TSUR_THRESHOLD = 230
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_MOUNTAIN_DAY']
    THRESHOLD3 = 298.15
    THRESHOLD4 = OFFSETS['R09_R06_THR_CLOUDY_MAX']
    THRESHOLD5 = OFFSETS['R09_R06_THR_CLOUDY_MIN']

t11surftemp_cond=np.logical_and( surftemp>TSUR_THRESHOLD,
                                 t11<THRESHOLD3)
qr0906cond=np.logical_and( qr09r06<THRESHOLD4,
                           qr09r06>THRESHOLD5)
TestOk = np.logical_and(
               np.logical_and(t11surftemp_cond,qr0906cond),
               np.logical_and( t11_t37_minus_threshold>THRESHOLD2,
                               t11_ts_minus_threshold<THRESHOLD1))

TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDayMount,isClearLandDayMount, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS,TestOkAll)
#----------------------------------
#thincoldCirrusTest land day mountain
#----------------------------------
args = {'title': "thincoldCirrusTest_All_LandDayMount",
        'xlable': 'T11-TS minus dynamic threshold',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t11_ts_minus_threshold,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}

THRESHOLD1 = OFFSETS['T11T12_OFFSET_LAND_DAY'] + OFFSETS['QUALITY_MARGIN_T11T12']#-7.0-1.0
THRESHOLD2 = OFFSETS['T11_OFFSET_LAND_DAY'] - OFFSETS['QUALITY_MARGIN_T11TSUR']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11T12_OFFSET_LAND_DAY']
    THRESHOLD2 = OFFSETS['T11_OFFSET_LAND_DAY']
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1 ,
                        t11_ts_minus_threshold<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDayMount, isClearLandDayMount, TestOk, THRESHOLD1,THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS, TestOkAll)
#----------------------------------
#HighcloudTestt85t11land land day mountain
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11land_All_LandDayMount",
            'xlable': 'Latidue',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD = OFFSETS['HI_T85T11_CM_LAND']
    TestOk = t85_t11_minus_threshold>THRESHOLD
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyLandDayMount,isClearLandDayMount,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyLandDayMount,isClearLandDayMount,isCloudyPPS,TestOkAll)

####################################
# NIGHT COAST
####################################
#----------------------------------
# VERSION 2014 COLD_CLOUD_TEST
#----------------------------------
args = {'title': "coldCloudTest_v2014_All_CoastNight",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
thr_temp=np.min([OFFSETS['T11_OFFSET_LAND_NIGHT'], OFFSETS['T11_OFFSET_SEA_NIGHT']])
offset = (thr_temp + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset + OFFSETS['QUALITY_MARGIN_T11TSUR']
TSURF_THRSHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']
if RunWithOutMargins:
    THRESHOLD = offset
    TSURF_THRSHOLD = 0.0#OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSURF_THRSHOLD) 
                        
TestOkAll=TestOk
plot_test(args,isCloudyCoastNight,isClearCoastNight,TestOk,np.min(THRESHOLD), show=False)
####################################
# DAY COAST MOUNTAIN
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyCoastDayMount,isClearCoastDayMount,isCloudyPPS)
TestOk = np.logical_and(np.equal(cloudtype,4),np.logical_or(isCloudyCoastDayMount,
                                                            isClearCoastDayMount))
print_stats_snow(cloudtype, isCloudyCoastDayMount, isClearCoastDayMount, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyCoastDayMount = np.logical_and(isCloudyCoastDayMount, np.not_equal(cloudtype,4))
isClearCoastDayMount = np.logical_and(isClearCoastDayMount, np.not_equal(cloudtype,4))
print_stats(cloudtype,isCloudyCoastDayMount,isClearCoastDayMount,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_CoastDayMount",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-11 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyCoastDayMount,isClearCoastDayMount,TestOk,np.min(THRESHOLD), show=False)
###################################
# DAY LAND 
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyLandDay,
                                                            isClearLandDay))
print_stats_snow(cloudtype, isCloudyLandDay, isClearLandDay, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyLandDay = np.logical_and(isCloudyLandDay, np.not_equal(isCloudtypeSnowOrIce,True))
isClearLandDay = np.logical_and(isClearLandDay, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS)
#----------------------------------
# newWaterBodiesTest
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND'] +OFFSETS['R06_OFFSET_LAND']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
if qr09r06 is not None:
    args = {'title': "CLEAR_newWaterBodiesTest_All_LandDay",
            'xlable': 'T37-T12',
            'ylable': 'qr09r06',
            'xvector': t37t12,
            'yvector': qr09r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD1 = 0.8
    THRESHOLD2 = r06_threshold
    THRESHOLD3 = 3.0
    THRESHOLD4 = 3.0
    if RunWithOutMargins:
        pass
    TestOk = np.logical_and(np.logical_and(qr09r06<THRESHOLD1,
                                           r06<THRESHOLD2),
                            np.logical_and(np.abs(t37t12)<THRESHOLD3,
                                           np.abs(t11t37)<THRESHOLD4))
    #TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudy,isClear, TestOk, 
                np.min(THRESHOLD1), np.min(THRESHOLD2), show=False)
    plot_test_2_lim(args,isCloudyLandDay,isClearLandDay, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    isCloudyLandDay = np.logical_and(isCloudyLandDay, np.equal(TestOk,False))
    isClearLandDay = np.logical_and(isClearLandDay, np.equal(TestOk,False))
    print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS)
#----------------------------------
#brightCloudTest land day 
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND'] +OFFSETS['R06_OFFSET_LAND']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
args = {'title': "brightCloudTest_All_LandDay",
      'xlable': 'T37-T12 minus dynamic threshold',
      'ylable': 'r06 ',
      'xvector': t37_t12_minus_threshold,
      'yvector': r06,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}

THRESHOLD2 = OFFSETS['T37T12_OFFSET_LAND_DAY'] +OFFSETS['QUALITY_MARGIN_T37T12']#
THRESHOLD1 = r06_threshold  + OFFSETS['QUALITY_MARGIN_R06']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37T12_OFFSET_LAND_DAY']
    THRESHOLD1 = r06_threshold 
    TestOk =  np.logical_and(t37_t12_minus_threshold>THRESHOLD2,
                             r06>THRESHOLD1)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDay,isClearLandDay, TestOk, 
                np.min(THRESHOLD1), THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS,TestOkAll)
#----------------------------------
#brightCloudTestNoSunglint3A Land day 
#----------------------------------
if qr16r06 is not None:
    args = {'title': "notusedif37_brightCloudTestNoSunglint3A_All_LandDay",
            'xlable': 'r06 ',
            'ylable': 'qr16r06',
            'xvector': r06,
            'yvector': qr16r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD2 = r06_threshold  +OFFSETS['QUALITY_MARGIN_R06']
    THRESHOLD1 = (OFFSETS['R06_R16_THR_SNOW'] + 
                  1*OFFSETS['QUALITY_MARGIN_QR16R06'])
    THRESHOLD3 = 1.35 - OFFSETS['QUALITY_MARGIN_QR16R06']
    TSUR_THRESHOLD = 0.0
    if RunWithOutMargins:
        THRESHOLD2 = r06_threshold 
        THRESHOLD1 = OFFSETS['R06_R16_THR_SNOW'] 
        THRESHOLD3 = 1.35
        TestOk = np.logical_and(np.logical_and(r06>THRESHOLD2,
                                               np.logical_and(qr16r06>THRESHOLD1,
                                                              qr16r06<THRESHOLD3)),
                                t37<=-9)
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudyLandDay,isClearLandDay, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    print_stats(cloudtype,isCloudyLandDay,isClearLandDay,
                isCloudyPPS,TestOkAll)
    args['title'] = "usealsoift37_but_less_generous_brightCloudTestNoSunglint3A_All_LandDay"
    THRESHOLD2 = THRESHOLD2 + 10
    TestOk = np.logical_and(r06>THRESHOLD2,
                            np.logical_and(qr16r06>THRESHOLD1,
                                           qr16r06<THRESHOLD3))
    plot_test_2_lim(args,isCloudyLandDay,isClearLandDay, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
#----------------------------------
#coldBrightCloudTest37 land day 
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND'] +OFFSETS['R06_OFFSET_LAND']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))

args = {'title': "coldBrightCloudTest37_All_LandDay",
        'xlable': 'R06',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': r06,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = OFFSETS['TSUR_THR_SNOW_MAX'] + OFFSETS['QUALITY_MARGIN_TSUR']#
THRESHOLD2 = r06_threshold + OFFSETS['QUALITY_MARGIN_R06']
THRESHOLD3 = 6.0 + OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'] 
    TSUR_THRESHOLD = OFFSETS['TSUR_THR_SNOW_MAX']
    THRESHOLD2 = r06_threshold
    THRESHOLD3 = 6.0
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold<THRESHOLD1,
                                       surftemp>TSUR_THRESHOLD),
                        np.logical_and(r06>THRESHOLD2,
                                       t37_t12_minus_threshold>THRESHOLD3))                        
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDay, isClearLandDay, TestOk, THRESHOLD1, np.min(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTestDay land day
#----------------------------------
args = {'title': "coldwaterCloudTestDay_All_LandDay",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#
THRESHOLD3 = 298.15 - 1.0
THRESHOLD4 = OFFSETS['R09_R06_THR_CLOUDY_MAX']- OFFSETS['QUALITY_MARGIN_QR09R06']#
THRESHOLD5 = OFFSETS['R09_R06_THR_CLOUDY_MIN']+ OFFSETS['QUALITY_MARGIN_QR09R06']#
TSUR_THRESHOLD = 230
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY']
    THRESHOLD3 = 298.15
    THRESHOLD4 = OFFSETS['R09_R06_THR_CLOUDY_MAX']
    THRESHOLD5 = OFFSETS['R09_R06_THR_CLOUDY_MIN']

t11surftemp_cond=np.logical_and( surftemp>TSUR_THRESHOLD,
                                 t11<THRESHOLD3)
qr0906cond=np.logical_and( qr09r06<THRESHOLD4,
                           qr09r06>THRESHOLD5)
TestOk = np.logical_and(
               np.logical_and(t11surftemp_cond,qr0906cond),
               np.logical_and( t11_t37_minus_threshold>THRESHOLD2,
                               t11_ts_minus_threshold<THRESHOLD1))

TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDay,isClearLandDay, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS,TestOkAll)
#----------------------------------
#coldBrightCloudTest land day 
#----------------------------------
r06_threshold_opaque =   100.0*thr_r06*OFFSETS['R06_GAIN_LAND_OPAQUE'] +OFFSETS['R06_OFFSET_LAND_OPAQUE']
r06_threshold_opaque = np.ma.array(
    r06_threshold_opaque, 
    mask=np.logical_or(r06.mask,thr_r06.mask))

args = {'title': "coldBrightCloudTest_All_LandDay",
        'xlable': 'R06',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': r06,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = 0
THRESHOLD2 = r06_threshold_opaque + OFFSETS['QUALITY_MARGIN_R06']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'] 
    TSUR_THRESHOLD = 0
    THRESHOLD2 = r06_threshold_opaque
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold<THRESHOLD1,
                                       surftemp>TSUR_THRESHOLD),
                        r06>THRESHOLD2)                        
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDay, isClearLandDay, TestOk, THRESHOLD1, np.min(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS,TestOkAll)

#----------------------------------
#thincoldCirrusTest land day mountain
#----------------------------------
args = {'title': "thincoldCirrusTest_All_LandDay",
        'xlable': 'T11-TS minus dynamic threshold',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t11_ts_minus_threshold,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}

THRESHOLD1 = OFFSETS['T11T12_OFFSET_LAND_DAY'] + OFFSETS['QUALITY_MARGIN_T11T12']#-7.0-1.0
THRESHOLD2 = OFFSETS['T11_OFFSET_LAND_DAY'] - OFFSETS['QUALITY_MARGIN_T11TSUR']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11T12_OFFSET_LAND_DAY']
    THRESHOLD2 = OFFSETS['T11_OFFSET_LAND_DAY']
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1 ,
                        t11_ts_minus_threshold<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandDay, isClearLandDay, TestOk, THRESHOLD1,THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS, TestOkAll)
#----------------------------------
#HighcloudTestt85t11land land day 
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11land_All_LandDay",
            'xlable': 'Latidue',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD = OFFSETS['HI_T85T11_CM_LAND']
    TestOk = t85_t11_minus_threshold>THRESHOLD
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyLandDay,isClearLandDay,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS,TestOkAll)
#----------------------------------
#ThinCirrusPrimaryTest land day
#----------------------------------
args = {'title': "non_existing_thinCirrusPrimaryTest_All_LandDay",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_DAY']+OFFSETS['QUALITY_MARGIN_T37T12']#2.0+0.3
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T37T12_OFFSET_LAND_DAY']
TestOk = t37_t12_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandDay,isClearLandDay,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyLandDay,isClearLandDay,isCloudyPPS,TestOkAll)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_LandDay",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-8 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyLandDay,isClearLandDay,TestOk,np.min(THRESHOLD), show=False)
###################################
# DAY COAST 
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyCoastDay,isClearCoastDay,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyCoastDay,
                                                            isClearCoastDay))
print_stats_snow(cloudtype, isCloudyCoastDay, isClearCoastDay, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyCoastDay = np.logical_and(isCloudyCoastDay, np.not_equal(isCloudtypeSnowOrIce,True))
isClearCoastDay = np.logical_and(isClearCoastDay, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyCoastDay,isClearCoastDay,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "add_coldCloudTest_v2014_All_CoastDay",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-8 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
#t11_t37 should be inv
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyCoastDay,isClearCoastDay,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyCoastDay,isClearCoastDay,isCloudyPPS)
###################################
# TWILIGHT COAST 
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyCoastTwilight,isClearCoastTwilight,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyCoastTwilight,
                                                            isClearCoastTwilight))
print_stats_snow(cloudtype, isCloudyCoastTwilight, isClearCoastTwilight, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyCoastTwilight = np.logical_and(isCloudyCoastTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
isClearCoastTwilight = np.logical_and(isClearCoastTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyCoastTwilight,isClearCoastTwilight,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "added_coldCloudTest_v2014_All_CoastTwilight",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-5.0 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyCoastTwilight,isClearCoastTwilight,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
#thincoldCirrusTest coast twilight
#----------------------------------
args = {'title': "added_thincoldCirrusTest_All_CoastTwilight",
        'xlable': 'T11-TS minus dynamic threshold',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t11_ts_minus_threshold,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
thr1=0.5*(OFFSETS['T11T12_OFFSET_SEA_NIGHT']+OFFSETS['T11T12_OFFSET_SEA_DAY'])
thr2=0.5*(OFFSETS['T11T12_OFFSET_LAND_NIGHT']+OFFSETS['T11T12_OFFSET_LAND_DAY'])
thr3=0.5*(OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE']+OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'])
thr4=0.5*(OFFSETS['T11_OFFSET_LAND_NIGHT_OPAQUE']+OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'])
offset_t11t12 =np.max([thr1,thr2])
offset_t11 =np.max([thr3,thr4])
THRESHOLD1 = offset_t11t12 + OFFSETS['QUALITY_MARGIN_T11T12']
THRESHOLD2 = offset_t11 - OFFSETS['QUALITY_MARGIN_T11TSUR']
if RunWithOutMargins:
    THRESHOLD1 = offset_t11t12
    THRESHOLD2 = offset_t11 

TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1 ,
                        t11_ts_minus_threshold<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyCoastTwilight, isClearCoastTwilight, TestOk, THRESHOLD1,THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyCoastTwilight,isClearCoastTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#ThinCirrusPrimaryTest coast twilight
#----------------------------------
args = {'title': "offset_day_night_slope_thinCirrusPrimaryTest_All_CoastTwilight",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
tmp_night=np.max([OFFSETS['T37T12_OFFSET_LAND_NIGHT'],OFFSETS['T37T12_OFFSET_SEA_NIGHT']])
tmp_day=np.max([ OFFSETS['T37T12_OFFSET_LAND_DAY'], OFFSETS['T37T12_OFFSET_SEA_DAY']])
tmp_thr =np.where(
    sunz>=90.0,
    0.5*(tmp_night+tmp_day),
    tmp_day)
#added threshold on slope from day to night
c = 1.0/(90.0-OFFSETS['MAX_SUNZEN_DAY'])
intercept =((1-c*90.0)*tmp_night +
	       c*90.0*tmp_day)
slope = c*(tmp_night-tmp_day)
tmp_thr =np.where(
    sunz>=90.0,
    tmp_night,
    slope * sunz + intercept)
#tmp_thr=np.max([OFFSETS['T37T12_OFFSET_LAND_NIGHT'], OFFSETS['T37T12_OFFSET_LAND_DAY'],OFFSETS['T37T12_OFFSET_SEA_NIGHT'], OFFSETS['T37T12_OFFSET_SEA_DAY']])
THRESHOLD = tmp_thr+OFFSETS['QUALITY_MARGIN_T37T12']#2.0+0.3
if RunWithOutMargins:
    THRESHOLD = tmp_thr
TestOk = t37_t12_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyCoastTwilight,isClearCoastTwilight,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyCoastTwilight,isClearCoastTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#thinCirrusSecondaryTest coast twilight 
#----------------------------------
args = {'title': "offset_day_night_slope_thinCirrusSecondaryTest_All_CoastTwilight_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
#tmp_thr =np.where(
#    sunz>=90.0,
#    0.5*(OFFSETS['T11T12_OFFSET_LAND_NIGHT'] + OFFSETS['T11T12_OFFSET_LAND_DAY']),
#    OFFSETS['T11T12_OFFSET_LANd_DAY'])
tmp_night = 0.5*(OFFSETS['T11T12_OFFSET_LAND_NIGHT'] + OFFSETS['T11T12_OFFSET_SEA_NIGHT'])
tmp_day = 0.5*(OFFSETS['T11T12_OFFSET_SEA_DAY'] + OFFSETS['T11T12_OFFSET_LAND_DAY'])
tmp_thr=np.max([tmp_day,tmp_night])
#added threshold on slope from day to night
c = 1.0/(90.0-OFFSETS['MAX_SUNZEN_DAY'])
intercept =((1-c*90.0)*tmp_night +
	       c*90.0*tmp_day)
slope = c*(tmp_night-tmp_day)
tmp_thr =np.where(
    sunz>=90.0,
    tmp_night,
    slope * sunz + intercept)
THRESHOLD = tmp_thr  + OFFSETS['QUALITY_MARGIN_T11T12']#-7.0-1.0
if RunWithOutMargins:
    THRESHOLD = tmp_thr 
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyCoastTwilight, isClearCoastTwilight, TestOk, np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyCoastTwilight,isClearCoastTwilight,isCloudyPPS, TestOkAll)

###################################
# TWILIGHT COAST INV
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyCoastTwilightInv,isClearCoastTwilightInv,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyCoastTwilightInv,
                                                            isClearCoastTwilightInv))
print_stats_snow(cloudtype, isCloudyCoastTwilightInv, isClearCoastTwilightInv, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyCoastTwilightInv = np.logical_and(isCloudyCoastTwilightInv, np.not_equal(isCloudtypeSnowOrIce,True))
isClearCoastTwilightInv = np.logical_and(isClearCoastTwilightInv, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyCoastTwilightInv,isClearCoastTwilightInv,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_CoastTwilightInv",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-7 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
plot_test(args,isCloudyCoastTwilightInv,isClearCoastTwilightInv,TestOk,np.min(THRESHOLD), show=False)

####################################
# TWILIGHT LAND MOUNTAIN
####################################
print "TWILIGHT LAND MOUNTAIN"
TestOkAll=None
print_stats(cloudtype,isCloudyLandTwilightMount,isClearLandTwilightMount,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyLandTwilightMount,
                                                            isClearLandTwilightMount))
print_stats_snow(cloudtype, isCloudyLandTwilightMount, isClearLandTwilightMount, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyLandTwilightMount = np.logical_and(isCloudyLandTwilightMount, np.not_equal(isCloudtypeSnowOrIce,True))
isClearLandTwilightMount = np.logical_and(isClearLandTwilightMount, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyLandTwilightMount,isClearLandTwilightMount,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "add_coldCloudTest_v2014_All_LandTwilightMount",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-12 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
#t11_t37 should be inv
TSUR_THRESHOLD=230.0
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp > TSUR_THRESHOLD)
TestOkAll=TestOk
plot_test(args,isCloudyLandTwilightMount,isClearLandTwilightMount,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyLandTwilightMount,isClearLandTwilightMount,isCloudyPPS)
####################################
# TWILIGHT LAND INV 
####################################
print "TWILIGHT LAND INV "
TestOkAll=None
print_stats(cloudtype,isCloudyLandTwilightInv,isClearLandTwilightInv,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyLandTwilightInv,
                                                            isClearLandTwilightInv))
print_stats_snow(cloudtype, isCloudyLandTwilightInv, isClearLandTwilightInv, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyLandTwilightInv = np.logical_and(isCloudyLandTwilightInv, np.not_equal(isCloudtypeSnowOrIce,True))
isClearLandTwilightInv = np.logical_and(isClearLandTwilightInv, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyLandTwilightInv,isClearLandTwilightInv,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "added_coldCloudTest_v2014_All_LandTwilightInv",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-8.0 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyLandTwilightInv,isClearLandTwilightInv,TestOk,np.min(THRESHOLD), show=False)

####################################
# TWILIGHT LAND
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyLandTwilight,
                                                            isClearLandTwilight))
print_stats_snow(cloudtype, isCloudyLandTwilight, isClearLandTwilight, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyLandTwilight = np.logical_and(isCloudyLandTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
isClearLandTwilight = np.logical_and(isClearLandTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "added_coldCloudTest_v2014_All_LandTwilight",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-6.5 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyLandTwilight,isClearLandTwilight,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
#coldcloudtest land twilight
#----------------------------------
args = {'title': "coldCloudTest_All_LandTwilight_without_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
tmp_thr = 0.5*(OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE']+ OFFSETS['T11_OFFSET_LAND_NIGHT_OPAQUE'])
THRESHOLD = tmp_thr -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #
if RunWithOutMargins:
    THRESHOLD = tmp_thr
TestOk = t11_ts_minus_threshold<THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandTwilight,isClearLandTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#reflectingCloudTest land twilight
#----------------------------------
args = {'title': "reflectingCloudTest_All_LandTwilight",
      'xlable': 'r06 not scaled',
      'ylable': 'T37-T2 minus dynamic threshold',
      'xvector': pseudor06_not_scaled,
      'yvector': t37_t12_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
vis_static = np.where(
      sunz>=90.0,
      OFFSETS['VIS_STATIC_LAND_OFFSET'],
      OFFSETS['VIS_STATIC_LAND_GAIN']*(90.0-sunz) + OFFSETS['VIS_STATIC_LAND_OFFSET'])
c = 1.0/(90.0-OFFSETS['MAX_SUNZEN_DAY'])
intercept =((1-c*90.0)*OFFSETS['T37T12_OFFSET_LAND_NIGHT'] +
	       c*90.0*OFFSETS['T37T12_OFFSET_LAND_DAY'])
slope = c*(OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -
		       OFFSETS['T37T12_OFFSET_LAND_DAY'])
tmp_thr =np.where(
    sunz>=90.0,
    OFFSETS['T37T12_OFFSET_LAND_NIGHT'],
    slope * sunz + intercept)

THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
THRESHOLD1 = tmp_thr +1*OFFSETS['QUALITY_MARGIN_T37T12']
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = vis_static
    THRESHOLD1 = tmp_thr 
TestOk = np.logical_and(2>0,
                        np.logical_and(t37_t12_minus_threshold>THRESHOLD1,
                                       pseudor06_not_scaled>THRESHOLD2))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandTwilight,isClearLandTwilight, TestOk, 
                np.min(THRESHOLD1), np.min(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#pseudo06CloudTest3A land twilight 
#----------------------------------
if qr16r06 is not None:
    args = {'title': "notusedif37_pseudo06CloudTest3A_All_LandTwilight",
            'xlable': 'r06 not scaled',
            'ylable': 'qr16r06',
            'xvector': pseudor06_not_scaled,
            'yvector': qr16r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    vis_static = np.where(
        sunz>=90.0,
        OFFSETS['VIS_STATIC_LAND_OFFSET'],
        OFFSETS['VIS_STATIC_LAND_GAIN']*(90.0-sunz) + 
        OFFSETS['VIS_STATIC_LAND_OFFSET'])
    THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
    THRESHOLD1 = (OFFSETS['R06_R16_THR_LOWSUN'] + 
                  1*OFFSETS['QUALITY_MARGIN_QR16R06'])
    TSUR_THRESHOLD = 0.0
    if RunWithOutMargins:
        THRESHOLD2 = vis_static
        THRESHOLD1 = OFFSETS['R06_R16_THR_LOWSUN']
        TestOk = np.logical_and(2>0,
                                np.logical_and(qr16r06>THRESHOLD1,
                                               pseudor06_not_scaled>THRESHOLD2))
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudy,isClear, TestOk, 
                np.min(THRESHOLD1), np.min(THRESHOLD2), show=False)
    plot_test_2_lim(args,isCloudyLandTwilight,isClearLandTwilight, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,
                isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest land twilight 
#----------------------------------
args = {'title': "coldwaterCloudTest_All_LandTwilight",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
tmp_thr = np.min([OFFSETS['T11_OFFSET_LAND_NIGHT'], OFFSETS['T11_OFFSET_LAND_DAY']])
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#
THRESHOLD1 = tmp_thr -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = tmp_thr 
TestOk = np.logical_and(surftemp>TSUR_THRESHOLD,
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                                       t11_ts_minus_threshold<THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyLandTwilight,isClearLandTwilight, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#watercloudtest land twilight 
#----------------------------------
args = {'title': "waterCloudTest_All_LandTwilight",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37']
if RunWithOutMargins:
    THRESHOLD = 0.0
TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandTwilight,isClearLandTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#coldcloudtest land twilight
#----------------------------------
args = {'title': "less_generous_coldCloudTest_All_LandTwilight_without_Tsur_limit_should_perhaps_be_with",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
tmp_thr = np.min([OFFSETS['T11_OFFSET_LAND_DAY'], OFFSETS['T11_OFFSET_LAND_NIGHT']])
THRESHOLD = tmp_thr -1.0 -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #-15.0-1.0
if RunWithOutMargins:
    THRESHOLD = tmp_thr -1.0 #-0.5
TestOk = t11_ts_minus_threshold<THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandTwilight,isClearLandTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
args['title'] = "less_generous_coldCloudTest_All_LandTwilight_without_Tsur_limit_should_perhaps_be_with_orginal"
tmp_thr = np.min([OFFSETS['T11_OFFSET_LAND_DAY'], OFFSETS['T11_OFFSET_LAND_NIGHT']])
THRESHOLD = tmp_thr -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] 
if RunWithOutMargins:
    THRESHOLD = tmp_thr 
TestOk = t11_ts_minus_threshold<THRESHOLD
plot_test(args,isCloudyLandTwilight,isClearLandTwilight,TestOk,THRESHOLD, show=False)

#----------------------------------
#thinCirrusSecondaryTest land twilight 
#----------------------------------
args = {'title': "offset_day_night_slope_thinCirrusSecondaryTest_All_LandTwilight_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}

#day night slope?
tmp_thr = 0.5*(OFFSETS['T11T12_OFFSET_LAND_NIGHT'] + OFFSETS['T11T12_OFFSET_LAND_DAY'])
tmp_night=OFFSETS['T11T12_OFFSET_LAND_NIGHT']
tmp_day=OFFSETS['T11T12_OFFSET_LAND_DAY']
c = 1.0/(90.0-OFFSETS['MAX_SUNZEN_DAY'])
intercept =((1-c*90.0)*tmp_night +
	       c*90.0*tmp_day)
slope = c*(tmp_night-tmp_day)
tmp_thr =np.where(
    sunz>=90.0,
    tmp_night,
    slope * sunz + intercept)
THRESHOLD = tmp_thr  + OFFSETS['QUALITY_MARGIN_T11T12']#-7.0-1.0
if RunWithOutMargins:
    THRESHOLD = tmp_thr 
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandTwilight, isClearLandTwilight, TestOk, np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#ThinCirrusPrimaryTest land twilight
#----------------------------------
args = {'title': "offset_slope_night_day_thinCirrusPrimaryTest_All_LandTwilight",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
#tmp_thr =np.where(
#    sunz>=90.0,
#    0.5*(OFFSETS['T37T12_OFFSET_LAND_NIGHT']+OFFSETS['T37T12_OFFSET_LAND_DAY']),
#    OFFSETS['T37T12_OFFSET_LAND_DAY'])
tmp_night=OFFSETS['T37T12_OFFSET_LAND_NIGHT']
tmp_day=OFFSETS['T37T12_OFFSET_LAND_DAY']
c = 1.0/(90.0-OFFSETS['MAX_SUNZEN_DAY'])
intercept =((1-c*90.0)*tmp_night +
	       c*90.0*tmp_day)
slope = c*(tmp_night-tmp_day)
tmp_thr =np.where(
    sunz>=90.0,
    tmp_night,
    slope * sunz + intercept)
THRESHOLD = tmp_thr+OFFSETS['QUALITY_MARGIN_T37T12']
if RunWithOutMargins:
    THRESHOLD = tmp_thr
TestOk = t37_t12_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyLandTwilight,isClearLandTwilight,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#HighcloudTestt85t11land land twilight 
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11land_All_LandTwilight",
            'xlable': 'Latidue',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD = OFFSETS['HI_T85T11_CM_LAND']
    TestOk = t85_t11_minus_threshold>THRESHOLD
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyLandTwilight,isClearLandTwilight,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyLandTwilight,isClearLandTwilight,isCloudyPPS,TestOkAll)

####################################
# DAY ICE
####################################
TestOkAll = None
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS, TestOkAll)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyIceDay,
                                                            isClearIceDay))
print_stats_snow(cloudtype, isCloudyIceDay, isClearIceDay, TestOk)
isCloudyIceDay = np.logical_and(isCloudyIceDay, np.not_equal(isCloudtypeSnowOrIce,True))
isClearIceDay = np.logical_and(isClearIceDay, np.not_equal(isCloudtypeSnowOrIce,True))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS, TestOkAll)

#----------------------#----------------------------------
# "T12T37_T11TS_All sea day ice
#----------------------------------
args = {'title': "T12T37_T11TS_All_SeaDayIce",
      'ylable': 'T37T12text',
      'xlable': 'T37 T12 minus dynamic threshold',
      'xvector': t37_t12_minus_threshold,
      'yvector': t37t12text,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
TestOk=np.where(TestOk,False,False)
plot_test_2_lim(args,isCloudyIceDay, isClearIceDay, TestOk, 0, 0, show=False)



#----------------------------------
#coldcloudtest sea day ice
#----------------------------------
args = {'title': "coldCloudTest_All_SeaDayIce_without_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #-15.0-1.0
if RunWithOutMargins:
    THRESHOLD =  OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE']
TSUR_THRESHOLD = 0.0
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = TestOk
plot_test(args,isCloudyIceDay,isClearIceDay,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#coldBrightCloudTest sea day ice
#----------------------------------
r06_threshold_opaque =   100.0*thr_r06*OFFSETS['R06_GAIN_SEA_OPAQUE'] +OFFSETS['R06_OFFSET_SEA_OPAQUE']
r06_threshold_opaque = np.ma.array(
    r06_threshold_opaque, 
    mask=np.logical_or(r06.mask,thr_r06.mask))

args = {'title': "coldBrightCloudTest_All_SeaDayIce",
        'xlable': 'R06',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': r06,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']#
#100.0*threshold->r06*throff->sm_acmg_r06_gain_sea_opaque +
#      throff->sm_acmg_r06_offset_sea_opaque;
#100.0*threshold->r06*throff->sm_acmg_r06_gain_sea_opaque +
#      throff->sm_acmg_r06_offset_sea_opaque;
#/* The satellite reflectances are given in percent. Thus
#   multiply tabulated reflectances by 100.0!	*/
THRESHOLD2 = r06_threshold_opaque + OFFSETS['QUALITY_MARGIN_R06']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] 
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
    THRESHOLD2 = r06_threshold_opaque
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold<THRESHOLD1,
                                       surftemp>TSUR_THRESHOLD),
                        r06>THRESHOLD2)                        
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceDay, isClearIceDay, TestOk, THRESHOLD1, np.min(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#brightCloudTest sea day ice
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_SEA'] +OFFSETS['R06_OFFSET_SEA']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
args = {'title': "maybe_add_t37t12thr_brightCloudTest_All_SeaDayIce",
      'xlable': 'T37-T12 minus dynamic threshold',
      'ylable': 'r06 ',
      'xvector': t37_t12_minus_threshold,
      'yvector': r06,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}

THRESHOLD2 = 5.0 +OFFSETS['QUALITY_MARGIN_T37T12']#
THRESHOLD1 = r06_threshold  + OFFSETS['QUALITY_MARGIN_R06']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = 5.0
    THRESHOLD1 = r06_threshold 
TestOk =   np.logical_and(t37_t12_minus_threshold<18,
    np.logical_and(t37_t12_minus_threshold>THRESHOLD2,
                   r06>THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIce,isClearIce, TestOk, 
                np.min(THRESHOLD1), THRESHOLD2, show=False)
plot_test_2_lim(args,isCloudyIceDay,isClearIceDay, TestOk, 
                np.min(THRESHOLD1), THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest sea day ice
#----------------------------------
args = {'title': "coldwaterCloudTest_All_SeaDayIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#
THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] 
TestOk = np.logical_and(surftemp>TSUR_THRESHOLD,
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                                       t11_ts_minus_threshold<THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceDay,isClearIceDay, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#brightCloudTestNoSunglint3A sea day ice 
#----------------------------------
if qr16r06 is not None:
    args = {'title': "brightCloudTestNoSunglint3A_All_SeaDayIce",
            'xlable': 'r06 ',
            'ylable': 'qr16r06',
            'xvector': r06,
            'yvector': qr16r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD2 = r06_threshold  +OFFSETS['QUALITY_MARGIN_R06']
    THRESHOLD1 = (OFFSETS['R06_R16_THR_SNOW'] + 
                  1*OFFSETS['QUALITY_MARGIN_QR16R06'])
    THRESHOLD3 = 1.35 - OFFSETS['QUALITY_MARGIN_QR16R06']
    TSUR_THRESHOLD = 0.0
    if RunWithOutMargins:
        THRESHOLD2 = r06_threshold 
        THRESHOLD1 = OFFSETS['R06_R16_THR_SNOW'] 
        THRESHOLD3 = 1.35
        TestOk = np.logical_and(r06>THRESHOLD2,
                                np.logical_and(qr16r06>THRESHOLD1,
                                               qr16r06<THRESHOLD3))
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudyIceDay,isClearIceDay, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    print_stats(cloudtype,isCloudyIceDay,isClearIceDay,
                isCloudyPPS,TestOkAll)
#----------------------------------
#coldcloudtest sea day ice
#----------------------------------
args = {'title': "coldCloudTest_All_SeaDayIce_with_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #-15.0-1.0
TSUR_THRESHOLD=OFFSETS['COLDEST_SEASURFACE_TEMP']+OFFSETS['QUALITY_MARGIN_TSUR'] #-15.0-1.0
if RunWithOutMargins:
    THRESHOLD =  OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE']
    TSUR_THRESHOLD=OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceDay,isClearIceDay,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#thinCirrusSecondaryTest sea day ice
#----------------------------------
args = {'title': "thinCirrusSecondaryTest_All_SeaDayIce_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}

THRESHOLD =OFFSETS['T11T12_OFFSET_SEA_DAY']  + OFFSETS['QUALITY_MARGIN_T11T12']
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY'] 
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceDay, isClearIceDay, TestOk, THRESHOLD, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS, TestOkAll)
#----------------------------------
#brightCloudTestSea sea day ice
#----------------------------------
r06_threshold =   100.0*thr_r06*OFFSETS['R06_GAIN_SEA'] +OFFSETS['R06_OFFSET_SEA']
r06_threshold = np.ma.array(
    r06_threshold, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
args = {'title': "maybe_remove_brightCloudTestSea_All_SeaDayIce",
      'xlable': 'T11',
      'ylable': 'r06 ',
      'xvector': t11,
      'yvector': r06,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T11_SEA_MIN'] -5.0 - OFFSETS['QUALITY_MARGIN_T11']#
THRESHOLD1 = r06_threshold  + OFFSETS['QUALITY_MARGIN_R06']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T11_SEA_MIN'] -5.0
    THRESHOLD1 = r06_threshold 
TestOk =   np.logical_and(t11<THRESHOLD2,
                          r06>THRESHOLD1)
#TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIce,isClearIce, TestOk, 
                np.min(THRESHOLD1), THRESHOLD2, show=False)
plot_test_2_lim(args,isCloudyWater,isClearWater, TestOk, 
                np.min(THRESHOLD1), THRESHOLD2, show=False)
plot_test_2_lim(args,isCloudyIceDay,isClearIceDay, TestOk, 
                np.min(THRESHOLD1), THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)

####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "add_coldCloudTest_v2014_All_SeaDayIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-8 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TSURF_THRSHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSURF_THRSHOLD) 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceDay,isClearIceDay,TestOk,np.min(THRESHOLD), show=False)

##############################3333
#Add arctic tests???
#----------------------------------
#arcticwaterCloudTest sea day
#----------------------------------
args = {'title': "arcticWaterCloudTest_All_SeaDayIce",
      'xlable': 't3712text',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': t37t12text,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC']
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T37T12TEXT']
if RunWithOutMargins:
    THRESHOLD1 = 0.0
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD1,t37t12text<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceDay, isClearIceDay, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticThinCirrusPrimaryTest sea day ice
#----------------------------------
args = {'title': "arcticThinCirrusPrimaryTest_All_SeaDayIce",
        'xlable': 'T37text',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
T37TEXT_OFFSET_SEA_DAY_ICE=0.7
THRESHOLD2 = T37TEXT_OFFSET_SEA_DAY_ICE+OFFSETS['QUALITY_MARGIN_T37TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = T37TEXT_OFFSET_SEA_DAY_ICE
    THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC']
TestOk = np.logical_and(t37_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceDay,isClearIceDay,TestOk,THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticThinWaterCloudTest sea day ice
#----------------------------------
args = {'title': "arcticThinWaterCloudTest_All_SeaDayIce",
        'xlable': 'T37t12text',
        'ylable': 'T37-T12 minus dynamic threshold inv',
        'xvector': t37t12text,
        'yvector': t37_t12_minus_threshold_inv,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37T12TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC_INV'] + OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC'] 
    THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC_INV']

TestOk = np.logical_and(t37_t12_minus_threshold_inv<THRESHOLD1, t37t12text<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceDay,isClearIceDay,TestOk, THRESHOLD1, THRESHOLD2,show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticWarmcloudtest Salomons sea day ice
#----------------------------------
args = {'title': "SalomonsArcticWarmCloudTest_All_SeaDayIce",
      'ylable': 'T11-Ts minus dynamic threshold',
      'xlable': 'T11-T37 minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['OFFSET_T11TSUR_ARCTIC_INV'] + OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC'] #3.0
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T37T12TEXT'] #-0.4
THRESHOLD3 = OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC'] #0.3
THRESHOLD4 = OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC'] #0.6
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold_inv>THRESHOLD1,
                                       t37_t12_minus_threshold_inv<THRESHOLD2),
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD3,
                                       t37t12text<THRESHOLD4))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceDay, isClearIceDay, TestOk, THRESHOLD1, THRESHOLD3,show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)

#----------------------------------
#arcticWarmCirrusSecondaryTest sea day
#----------------------------------
args = {'title': "arcticWarmCirrusSecondaryTest_All_SeaDayIce_tsur",
      'xlable': 'Latitude',
      'ylable': 'T11-T12 minus dynamic threshold inv',
      'xvector': latitude,
      'yvector': t11_t12_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}

THRESHOLD = OFFSETS['OFFSET_T11T12_ARCTIC_INV'] -1* OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
if RunWithOutMargins:
    pass
    #Not run without thresholds!
    #THRESHOLD = OFFSETS['OFFSET_T11T12_ARCTIC_INV'] #MARGIN are offset
TestOk = t11_t12_minus_threshold_inv<THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceDay, isClearIceDay, TestOk, THRESHOLD,show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticThinCirrusSecondaryTest sea day ice
#----------------------------------
args = {'title': "arcticthinCirrusSecondaryTest_All_SeaDayIce",
        'xlable': 'T37text',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']
THRESHOLD1 = OFFSETS['OFFSET_T11T12_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
print OFFSETS['OFFSET_T11T12_ARCTIC']
plt.show()
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = OFFSETS['OFFSET_T11T12_ARCTIC']
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceDay, isClearIceDay, TestOk, THRESHOLD1, THRESHOLD2,show=False)
print_stats(cloudtype,isCloudyIceDay,isClearIceDay,isCloudyPPS,TestOkAll)
####################################
# TWILIGHT ICE
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS, TestOkAll)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyIceTwilight,
                                                            isClearIceTwilight))
print_stats_snow(cloudtype, isCloudyIceTwilight, isClearIceTwilight, TestOk)
isCloudyIceTwilight = np.logical_and(isCloudyIceTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
isClearIceTwilight = np.logical_and(isClearIceTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#arcticwaterCloudTest sea twilight
#----------------------------------
args = {'title': "arcticWaterCloudTest_All_SeaTwilightIce",
      'xlable': 't3712text',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': t37t12text,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC']
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T37T12TEXT']
if RunWithOutMargins:
    THRESHOLD1 = 0.0-0.5
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD1,t37t12text<THRESHOLD2)
#TestOkAll=TestOk
plot_test_2_lim(args,isCloudyIceTwilight, isClearIceTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticThinCirrusSecondaryTest sea twilight ice
#----------------------------------
args = {'title': "maybe_add_arcticthinCirrusSecondaryTest_All_SeaTwilightIce",
        'xlable': 'T37text',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']
THRESHOLD1 = OFFSETS['OFFSET_T11T12_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = 0.0#OFFSETS['OFFSET_T11T12_ARCTIC']
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
TestOkAll=TestOk
plot_test_2_lim(args,isCloudyIceTwilight, isClearIceTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS, TestOkAll)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "add_coldCloudTest_v2014_All_SeaTwilightIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
#-4.0, -4.5
offset = (-13.0 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TSURF_THRSHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSURF_THRSHOLD) 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceTwilight,isClearIceTwilight,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#thinCirrusSecondaryTest sea twilight ice
#----------------------------------
args = {'title': "thinCirrusSecondaryTest_All_SeaTwilightIce_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
#slope offset?
tmp_thr = 0.5*(OFFSETS['T11T12_OFFSET_SEA_NIGHT'] + OFFSETS['T11T12_OFFSET_SEA_DAY'])
THRESHOLD =tmp_thr  + OFFSETS['QUALITY_MARGIN_T11T12']
if RunWithOutMargins:
    THRESHOLD = tmp_thr 
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceTwilight, isClearIceTwilight, TestOk, THRESHOLD, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#textureIrvisTest sea twilight  ice
#----------------------------------
args = {'title': "textureIrvisTest_All_SeaTwilightIce",
        'xlable': 'r06text minus dynamic threshold',
        'ylable': 'T11text minus dynamic threshold',
        'xvector': t11text,
        'yvector': r06text,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['R06TEXT_OFFSET_SUNGLINT'] + OFFSETS['QUALITY_MARGIN_R06TEXT']
THRESHOLD1 = OFFSETS['T11TEXT_OFFSET_SUNGLINT'] + OFFSETS['QUALITY_MARGIN_T11TEXT']
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
SUNZ_THRESHOLD = OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['R06TEXT_OFFSET_SUNGLINT']
    THRESHOLD1 = OFFSETS['T11TEXT_OFFSET_SUNGLINT']
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
    SUNZ_THRESHOLD = OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT']
TestOk = np.logical_and(np.logical_and(t11text>THRESHOLD2, 
                                       r06text>THRESHOLD1),
                        np.logical_and(sunz<SUNZ_THRESHOLD,
                                       surftemp> TSUR_THRESHOLD))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceTwilight, isClearIceTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#textureNightTest sea twilight  ice
#----------------------------------
args = {'title': "sunz_thr_added_texturNightTest_All_SeaTwilightIce",
        'xlable': 'T11text minus dynamic threshold',
        'ylable': 'T37-T12text minus dynamic threshold',
        'xvector': t11text,
        'yvector': t37t12text,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_DAY'] + OFFSETS['QUALITY_MARGIN_T11TEXT']#0.85+0.15
THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_DAY'] + OFFSETS['QUALITY_MARGIN_T37T12TEXT']#0.75+0.15
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
SUNZ_THRESHOLD = OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_DAY']
    THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_DAY']
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(np.logical_and(sunz>SUNZ_THRESHOLD,
                                       t11text>THRESHOLD2, ),
                        np.logical_and(t37t12text>THRESHOLD1, 
                                       surftemp> TSUR_THRESHOLD))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceTwilight, isClearIceTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
args['title']= "sunz_thr_added_texturNightTest_All_SeaTwilightIce_orginal"
TestOk = np.logical_and(t11text>THRESHOLD2,
                        np.logical_and(t37t12text>THRESHOLD1, 
                                       surftemp> TSUR_THRESHOLD))
plot_test_2_lim(args,isCloudyIceTwilight, isClearIceTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
#----------------------------------
#coldcloudtest sea twilight ice
#----------------------------------
args = {'title': "latitude_problem_coldCloudTest_All_SeaTwilightIce_without_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #-15.0-1.0
if RunWithOutMargins:
    THRESHOLD =  OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE']
TSUR_THRESHOLD=0.0
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceTwilight,isClearIceTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#reflectingCloudTest sea twilight ice
#----------------------------------
print "test only done of low prob of sunglint"
args = {'title': "reflectingCloudTest_All_SeaTwilightIce",
      'xlable': 'r06 not scaled',
      'ylable': 'T37-T2 minus dynamic threshold',
      'xvector': pseudor06_not_scaled,
      'yvector': t37_t12_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
vis_static = np.where(
    sunz>=90.0,
    OFFSETS['VIS_STATIC_SEA_OFFSET'],
    OFFSETS['VIS_STATIC_SEA_GAIN']*(90.0-sunz) + OFFSETS['VIS_STATIC_SEA_OFFSET'])
tmp_thr = 0.5*(OFFSETS['T37T12_OFFSET_SEA_NIGHT'] + OFFSETS['T37T12_OFFSET_SEA_DAY'])
THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
THRESHOLD1 = tmp_thr +1*OFFSETS['QUALITY_MARGIN_T37T12']
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = vis_static
    THRESHOLD1 = tmp_thr 
TestOk = np.logical_and(2>0,
                        np.logical_and(t37_t12_minus_threshold>THRESHOLD1,
                                       pseudor06_not_scaled>THRESHOLD2))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceTwilight,isClearIceTwilight, TestOk, 
                THRESHOLD1, np.mean(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#pseudo06CloudTest3A sea twilight ice
#----------------------------------
if qr16r06 is not None:
    args = {'title': "notusedif37_pseudo06CloudTest3A_All_SeaTwilightIce",
            'xlable': 'r06 not scaled',
            'ylable': 'qr16r06',
            'xvector': pseudor06_not_scaled,
            'yvector': qr16r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    vis_static = np.where(
        sunz>=90.0,
        OFFSETS['VIS_STATIC_SEA_OFFSET'],
        (OFFSETS['VIS_STATIC_SEA_GAIN']*(90.0-sunz) + 
         OFFSETS['VIS_STATIC_SEA_OFFSET']))
    THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
    THRESHOLD1 = (OFFSETS['R06_R16_THR_LOWSUN'] + 
                  1*OFFSETS['QUALITY_MARGIN_QR16R06'])
    TSUR_THRESHOLD = 0.0
    if RunWithOutMargins:
        THRESHOLD2 = vis_static
        THRESHOLD1 = OFFSETS['R06_R16_THR_LOWSUN']
    TestOk = np.logical_and(2>0,
                            np.logical_and(qr16r06>THRESHOLD1,
                                           pseudor06_not_scaled>THRESHOLD2))
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudy,isClear, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    plot_test_2_lim(args,isCloudyIceTwilight,isClearIceTwilight, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,
                isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest sea twilight ice
#----------------------------------
args = {'title': "coldwaterCloudTest_All_SeaTwilightIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
#thr slope?
tmp_thr = np.min([OFFSETS['T11_OFFSET_SEA_NIGHT'], OFFSETS['T11_OFFSET_SEA_DAY']])
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#
THRESHOLD1 = tmp_thr -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = tmp_thr 
TestOk = np.logical_and(surftemp>TSUR_THRESHOLD,
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                                       t11_ts_minus_threshold<THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceTwilight,isClearIceTwilight, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#watercloudtest sea twilight ice
#----------------------------------
args = {'title': "waterCloudTest_All_SeaTwilightIce",
      'xlable': 'Latitdue',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
if RunWithOutMargins:
    THRESHOLD = 0.0 
TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceTwilight,isClearIceTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#coldcloudtest sea twilight ice
#----------------------------------
args = {'title': "coldCloudTest_All_SeaTwilightIce_with_tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
#tmp slope?
tmp_thr=np.min([OFFSETS['T11_OFFSET_SEA_DAY'],OFFSETS['T11_OFFSET_SEA_NIGHT']])
THRESHOLD = tmp_thr -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #-15.0-1.0
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']
if RunWithOutMargins:
    THRESHOLD = tmp_thr
    TSUR_THRESHOLD= OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceTwilight,isClearIceTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#reflectingCloudTest sea twilight ice
#----------------------------------
print "this test is for sea! not ice!"
args = {'title': "removed_reflectingCloudTestSea_All_SeaTwilightIce",
      'xlable': 'r06 not scaled',
      'ylable': 'T11',
      'xvector': pseudor06_not_scaled,
      'yvector': t11,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
vis_static = np.where(
    sunz>=90.0,
    OFFSETS['VIS_STATIC_SEA_OFFSET'],
    OFFSETS['VIS_STATIC_SEA_GAIN']*(90.0-sunz) + OFFSETS['VIS_STATIC_SEA_OFFSET'])
THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
THRESHOLD1 = OFFSETS['T11_SEA_MIN'] -1*OFFSETS['QUALITY_MARGIN_T11']
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = vis_static
    THRESHOLD1 = OFFSETS['T11_SEA_MIN'] #+1*OFFSETS['QUALITY_MARGIN_T11']
TestOk = np.logical_and(2>0,
                        np.logical_and(t11<THRESHOLD1,
                                       pseudor06_not_scaled>THRESHOLD2))
#TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudy,isClear, TestOk, 
                THRESHOLD1, np.mean(THRESHOLD2), show=False)
plot_test_2_lim(args,isCloudyIceTwilight,isClearIceTwilight, TestOk, 
                THRESHOLD1, np.mean(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyIceTwilight,isClearIceTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
# "T12T37_T11TS_All sea twilight ice
#----------------------------------
args = {'title': "T12T37_T11TS_All_SeaTwilightIce",
      'ylable': 'T12-Ts minus dynamic threshold',
      'xlable': 'T11-T37 minus dynamic threshold',
      'xvector': t37_t12_minus_threshold,
      'yvector': t11_ts_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
TestOk=np.where(TestOk,False,False)
plot_test_2_lim(args,isCloudyIceTwilight, isClearIceTwilight, TestOk, 0, 0, show=False)

####################################
# TWILIGHT SEA NO ICE 
####################################
TestOkAll = None
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS, TestOkAll)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyWaterTwilight,
                                                            isClearWaterTwilight))
print_stats_snow(cloudtype, isCloudyWaterTwilight, isClearWaterTwilight, TestOk)
isCloudyWaterTwilight = np.logical_and(isCloudyWaterTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
isClearWaterTwilight = np.logical_and(isClearWaterTwilight, np.not_equal(isCloudtypeSnowOrIce,True))
print "***snow test, not implemented***"
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "added_coldCloudTest_v2014_All_SeaTwilightNoIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
#-4.0, -4.5
offset = (-4.0 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = t11_ts_minus_threshold < THRESHOLD
TestOkAll=TestOk
plot_test(args,isCloudyWaterTwilight,isClearWaterTwilight,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
#thinCirrusSecondaryTest sea twilightno ice
#----------------------------------
args = {'title': "strange_thr_thinCirrusSecondaryTest_All_SeaTwilightNoIce_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
#slope thr?
tmp_thr = 0.5*(OFFSETS['T11T12_OFFSET_SEA_NIGHT'] + OFFSETS['T11T12_OFFSET_SEA_DAY'])
THRESHOLD =tmp_thr  + OFFSETS['QUALITY_MARGIN_T11T12']
if RunWithOutMargins:
    THRESHOLD = tmp_thr
TestOk = t11_t12_minus_threshold>THRESHOLD 
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterTwilight, isClearWaterTwilight, TestOk, THRESHOLD, show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#textureIrvisTest sea twilight no ice
#----------------------------------
args = {'title': "maybe_less_generous_textureIrvisTest_All_SeaTwilightNoIce",
        'xlable': 'r06text minus dynamic threshold',
        'ylable': 'T11text minus dynamic threshold',
        'xvector': t11text,
        'yvector': r06text,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
#Not withour margins!!
THRESHOLD2 = OFFSETS['R06TEXT_OFFSET_SUNGLINT'] + OFFSETS['QUALITY_MARGIN_R06TEXT']
THRESHOLD1 = OFFSETS['T11TEXT_OFFSET_SUNGLINT'] + OFFSETS['QUALITY_MARGIN_T11TEXT']
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
SUNZ_THRESHOLD = OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['R06TEXT_OFFSET_SUNGLINT']
    THRESHOLD1 = OFFSETS['T11TEXT_OFFSET_SUNGLINT']
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
    SUNZ_THRESHOLD = OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT']
TestOk = np.logical_and(np.logical_and(t11text>THRESHOLD2, 
                                       r06text>THRESHOLD1),
                        np.logical_and(sunz<SUNZ_THRESHOLD,
                                       surftemp> TSUR_THRESHOLD))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterTwilight, isClearWaterTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS, TestOkAll)
#----------------------------------
#textureNightTest sea twilight no ice
#----------------------------------
args = {'title': "sunz_thr_added_texturNightTest_All_SeaTwilightNoIce",
        'xlable': 'T11text minus dynamic threshold',
        'ylable': 'T37-T12text minus dynamic threshold',
        'xvector': t11text,
        'yvector': t37t12text,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
#only for almost night
THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_DAY'] + OFFSETS['QUALITY_MARGIN_T11TEXT']
THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_DAY'] + OFFSETS['QUALITY_MARGIN_T37T12TEXT']
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
SUNZ_THRESHOLD = OFFSETS['MAX_SUNZEN_TWILIGHT_LIGHT']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_DAY']
    THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_DAY']
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(np.logical_and(t11text>THRESHOLD2, 
                                       sunz>SUNZ_THRESHOLD),
                        np.logical_and(t37t12text>THRESHOLD1, 
                                       surftemp> TSUR_THRESHOLD))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterTwilight, isClearWaterTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#Not in code yet!
args['title'] = "sunz_thr_added_texturNightTest_All_SeaTwilightNoIce_orginal"
TestOk = np.logical_and(t11text>THRESHOLD2,                                      
                        np.logical_and(t37t12text>THRESHOLD1, 
                                       surftemp> TSUR_THRESHOLD))
plot_test_2_lim(args,isCloudyWaterTwilight, isClearWaterTwilight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
#----------------------------------
#coldcloudtest sea twilightno ice
#----------------------------------
args = {'title': "coldCloudTest_All_SeaTwilightNoIce_without_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] 
if RunWithOutMargins:
    THRESHOLD =  OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE']
TSUR_THRESHOLD=0.0
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterTwilight,isClearWaterTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#reflectingCloudTest sea twilightno ice
#----------------------------------
args = {'title': "reflectingCloudTest_All_SeaTwilightNoIce",
      'xlable': 'r06 not scaled',
      'ylable': 'T37-T2 minus dynamic threshold',
      'xvector': pseudor06_not_scaled,
      'yvector': t37_t12_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
#thr slope?
vis_static = np.where(
    sunz>=90.0,
    OFFSETS['VIS_STATIC_SEA_OFFSET'],
    OFFSETS['VIS_STATIC_SEA_GAIN']*(90.0-sunz) + OFFSETS['VIS_STATIC_SEA_OFFSET'])
tmp_thr = 0.5*(OFFSETS['T37T12_OFFSET_SEA_NIGHT'] + OFFSETS['T37T12_OFFSET_SEA_DAY'])
THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
THRESHOLD1 = tmp_thr +1*OFFSETS['QUALITY_MARGIN_T37T12']
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = vis_static
    THRESHOLD1 = tmp_thr +1*OFFSETS['QUALITY_MARGIN_T37T12']
#PSUNGLINT_THR= OFFSETS['SUNGLINT_PROBABILITY']
TestOk = np.logical_and(2>0,
                        np.logical_and(t37_t12_minus_threshold>THRESHOLD1,
                                       pseudor06_not_scaled>THRESHOLD2))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterTwilight,isClearWaterTwilight, TestOk, 
                THRESHOLD1, np.mean(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS, TestOkAll)

#----------------------------------
#pseudo06CloudTest3A sea twilightno ice
#----------------------------------
if qr16r06 is not None:
    args = {'title': "notusedif37_pseudo06CloudTest3A_All_SeaTwilightNoIce",
            'xlable': 'r06 not scaled',
            'ylable': 'qr16r06',
            'xvector': pseudor06_not_scaled,
            'yvector': qr16r06,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    vis_static = np.where(
        sunz>=90.0,
        OFFSETS['VIS_STATIC_SEA_OFFSET'],
        (OFFSETS['VIS_STATIC_SEA_GAIN']*(90.0-sunz) + 
         OFFSETS['VIS_STATIC_SEA_OFFSET']))
    THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
    THRESHOLD1 = OFFSETS['R06_R16_THR_LOWSUN'] + 1*OFFSETS['QUALITY_MARGIN_QR16R06']
    TSUR_THRESHOLD = 0.0
    if RunWithOutMargins:
        THRESHOLD2 = vis_static
        THRESHOLD1 = OFFSETS['R06_R16_THR_LOWSUN']
    TestOk = np.logical_and(2>0,
                            np.logical_and(qr16r06>THRESHOLD1,
                                           pseudor06_not_scaled>THRESHOLD2))
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test_2_lim(args,isCloudy,isClear, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    plot_test_2_lim(args,isCloudyWaterTwilight,isClearWaterTwilight, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), show=False)
    print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,
                isCloudyPPS,TestOkAll)
#----------------------------------
#ColdWaterCloudTest sea twilight no ice
#----------------------------------
args = {'title': "maybe_more_generous_coldwaterCloudTest_All_SeaTwilightNoIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}#lower thr ?
tmp_thr = np.min([OFFSETS['T11_OFFSET_SEA_NIGHT'], OFFSETS['T11_OFFSET_SEA_DAY']])
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#
THRESHOLD1 = tmp_thr -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = 0.0#-3.0
    THRESHOLD1 = tmp_thr 
TestOk = np.logical_and(surftemp>TSUR_THRESHOLD,
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD2,
                                       t11_ts_minus_threshold<THRESHOLD1))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterTwilight,isClearWaterTwilight, TestOk, 
                THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#watercloudtest sea twilight no ice
#----------------------------------
args = {'title': "less_generous_waterCloudTest_All_SeaTwilightNoIce",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']+ OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
#higher thr
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']
TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterTwilight,isClearWaterTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
args['title'] = "less_generous_waterCloudTest_All_SeaTwilightNoIce_orginal"
THRESHOLD = OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
if RunWithOutMargins:
    THRESHOLD = 0
TestOk = t11_t37_minus_threshold>THRESHOLD
plot_test(args,isCloudyWaterTwilight,isClearWaterTwilight,TestOk,THRESHOLD, show=False)
#----------------------------------
#coldcloudtest sea twilight no ice
#----------------------------------
args = {'title': "coldCloudTest_All_SeaTwilightNoIce_with_tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
tmp_thr=np.min([OFFSETS['T11_OFFSET_SEA_DAY'],OFFSETS['T11_OFFSET_SEA_NIGHT']])
THRESHOLD = tmp_thr -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #-15.0-1.0
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']
if RunWithOutMargins:
    THRESHOLD = tmp_thr
    TSUR_THRESHOLD= OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterTwilight,isClearWaterTwilight,TestOk,THRESHOLD, show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
#reflectingCloudTestSea sea twilight no ice
#----------------------------------
args = {'title': "less_generous_reflectingCloudTestSea_All_SeaTwilightNoIce",
      'xlable': 'r06 not scaled',
      'ylable': 'T11',
      'xvector': pseudor06_not_scaled,
      'yvector': t11,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
#lowe thr visual clear limit above
vis_static = np.where(
    sunz>=90.0,
    OFFSETS['VIS_STATIC_SEA_OFFSET'],
    OFFSETS['VIS_STATIC_SEA_GAIN']*(90.0-sunz) + OFFSETS['VIS_STATIC_SEA_OFFSET'])
THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
THRESHOLD1 = OFFSETS['T11_SEA_MIN'] -3.0 -1*OFFSETS['QUALITY_MARGIN_T11']
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = vis_static
    THRESHOLD1 = OFFSETS['T11_SEA_MIN'] -3.0 #+1*OFFSETS['QUALITY_MARGIN_T11']
TestOk = np.logical_and(2>0,
                        np.logical_and(t11<THRESHOLD1,
                                       pseudor06_not_scaled>THRESHOLD2))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterTwilight,isClearWaterTwilight, TestOk, 
                THRESHOLD1, np.mean(THRESHOLD2), show=False)
print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS,TestOkAll)
#----------------------------------
args['title'] = "less_generous_reflectingCloudTestSea_All_SeaTwilightNoIce_orginal"
vis_static = np.where(
    sunz>=90.0,
    OFFSETS['VIS_STATIC_SEA_OFFSET'],
    OFFSETS['VIS_STATIC_SEA_GAIN']*(90.0-sunz) + OFFSETS['VIS_STATIC_SEA_OFFSET'])
THRESHOLD2 = vis_static +OFFSETS['QUALITY_MARGIN_VIS_STATIC']
THRESHOLD1 = OFFSETS['T11_SEA_MIN'] -1*OFFSETS['QUALITY_MARGIN_T11']
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD2 = vis_static
    THRESHOLD1 = OFFSETS['T11_SEA_MIN'] 
TestOk = np.logical_and(2>0,
                        np.logical_and(t11<THRESHOLD1,
                                       pseudor06_not_scaled>THRESHOLD2))
plot_test_2_lim(args,isCloudyWaterTwilight,isClearWaterTwilight, TestOk, 
                THRESHOLD1, np.mean(THRESHOLD2), show=False)
#----------------------------------
# sea twilight no ice
#----------------------------------
args = {'title': "T12T37_T11TS_All_SeaTwilightNoIce",
      'ylable': 'T12-Ts minus dynamic threshold',
      'xlable': 'T11-T37 minus dynamic threshold',
      'xvector': t37_t12_minus_threshold,
      'yvector': t11_ts_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
TestOk=np.where(TestOk,False,False)
plot_test_2_lim(args,isCloudyWaterTwilight, isClearWaterTwilight, TestOk, 0, 0, show=False)

#----------------------------------
#HighcloudTestt85t11sea sea twilight 
#----------------------------------
if t85_t11_minus_threshold is not None:
    args = {'title': "HighcloudTestt85t11sea_All_SeaTwilightNoIce",
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': latitude,
            'yvector': t85_t11_minus_threshold,
            'isCloudyPPS': isCloudyPPS,
            'isClearPPS': isClearPPS,   
            }
    THRESHOLD = OFFSETS['HI_T85T11_CM_SEA']
    TestOk = t85_t11_minus_threshold>THRESHOLD
    TestOkAll = keep_combined_ok(TestOk,TestOkAll)
    plot_test(args,isCloudyWaterTwilight,isClearWaterTwilight,TestOk,THRESHOLD, show=False)
    print_stats(cloudtype,isCloudyWaterTwilight,isClearWaterTwilight,isCloudyPPS,TestOkAll)


####################################
# NIGHT ICE
####################################
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS)
TestOkAll=None
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "added_coldCloudTest_v2014_All_SeaNightIce",
      'xlable': 'surftemp',#'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
#
offset = (-6.0 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                       t11_t37_minus_threshold_inv<0)
TestOkAll=TestOk
plot_test(args,isCloudyIceNight,isClearIceNight,TestOk,np.min(THRESHOLD), show=False)
#----------------------------------
#arcticwaterCloudTest sea night
#----------------------------------
print "arcticwatercloudTest sea night ice"
args = {'title': "maybe_less_generous_arcticWaterCloudTest_All_SeaNightIce",
      'xlable': 't3712text',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': t37t12text,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC']
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T37T12TEXT']
if RunWithOutMargins:
    THRESHOLD1 = 0.0
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD1,t37t12text<THRESHOLD2)
TestOkAll=TestOk
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)
#----------------------------------
#coldcloudtest sea night ice
#----------------------------------
args = {'title': "coldCloudTest_All_SeaNightIce_without_Tsur_limit",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['OFFSET_T11TSUR_ARCTIC'] -1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC'] #-15.0-1.0
if RunWithOutMargins:
    THRESHOLD = OFFSETS['OFFSET_T11TSUR_ARCTIC'] 
TSUR_THRESHOLD=0.0
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceNight,isClearIceNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

#----------------------------------
#arcticThinCirrusPrimaryTest sea night ice
#----------------------------------
args = {'title': "maybe_more_generous_arcticThinCirrusPrimaryTest_All_SeaNightIce",
        'xlable': 'T37text',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC']
TestOk = np.logical_and(t37_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight,isClearIceNight,TestOk,THRESHOLD1, THRESHOLD2, show=True)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)
args['title']= "arcticthinCirrusPrimaryTest_all_SeaNightIce_original"
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC']
TestOk = np.logical_and(t37_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2, show=True)
#----------------------------------
#arcticThinWaterCloudTest sea night ice
#----------------------------------
args = {'title': "arcticThinWaterCloudTest_All_SeaNightIce",
        'xlable': 'T37t12text',
        'ylable': 'T37-T12 minus dynamic threshold inv',
        'xvector': t37t12text,
        'yvector': t37_t12_minus_threshold_inv,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37T12TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC_INV'] + OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC'] 
    THRESHOLD1 = OFFSETS['OFFSET_T37T12_ARCTIC_INV']

TestOk = np.logical_and(t37_t12_minus_threshold_inv<THRESHOLD1, t37t12text<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight,isClearIceNight,TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticWarmcloudtest sea night ice
#----------------------------------
args = {'title': "maybe_reintroduce_arcticWarmCloudTest_All_SeaNightIce",
      'xlable': 't37t12text',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t37t12text,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 2.0#-1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC_INV']
THRESHOLD2 = 1.5 #OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T37TEXT']
TestOk = np.logical_and(t11_ts_minus_threshold>THRESHOLD1,t37t12text<THRESHOLD2)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2, show=False)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)
#----------------------------------
args = {'title': "removed_arcticWarmCloudTest_All_SeaNightIce",
      'ylable': 'T11-Ts minus dynamic threshold',
      'xlable': 'T37text',
      'xvector': t37text,
      'yvector': t11_ts_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
sm_acmg_quality_margin_t11tsur_arctic_inv =-5.0
THRESHOLD1 = -1* sm_acmg_quality_margin_t11tsur_arctic_inv
THRESHOLD4 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC']
#-   if ((f->t11 - f->tsur > dt11tsur_thr_inv + nmargin*dt11tsur_margin_inv)
#-       && (f->t37_text < text_t37_thr - nmargin*text_t37_margin)) {
TestOk = np.logical_and(t11_ts_minus_threshold_inv>THRESHOLD1,
                        t37text<THRESHOLD4)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD4, show=False)
#----------------------------------
#arcticWarmcloudtest Salomons sea night ice
#----------------------------------
args = {'title': "SalomonsArcticWarmCloudTest_All_SeaNightIce",
      'ylable': 'T11-Ts minus dynamic threshold',
      'xlable': 'T11-T37 minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['OFFSET_T11TSUR_ARCTIC_INV'] + OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC'] #3.0
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T37T12TEXT'] #-0.4
THRESHOLD3 = OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC'] #0.3
THRESHOLD4 = OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC'] #0.6
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['OFFSET_T11TSUR_ARCTIC_INV'] + OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC'] #3.0
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T37T12TEXT'] #-0.4
    THRESHOLD3 = OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC'] #0.3
    THRESHOLD4 = OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC'] #0.6

TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold>THRESHOLD1,
                                       t37_t12_minus_threshold_inv<THRESHOLD4),
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD3,
                                       t37t12text<THRESHOLD2))
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD3, show=False)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticWarmCirrusSecondaryTest sea night
#----------------------------------
args = {'title': "arcticWarmCirrusSecondaryTest_All_SeaNightIce_tsur",
      'xlable': 'Latitude',
      'ylable': 'T11-T12 minus dynamic threshold inv',
      'xvector': latitude,
      'yvector': t11_t12_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['OFFSET_T11T12_ARCTIC_INV'] -1* OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
if RunWithOutMargins:
    pass
    #Not run without thresholds!
    #THRESHOLD = OFFSETS['OFFSET_T11T12_ARCTIC_INV'] #MARGIN are offset
TestOk = t11_t12_minus_threshold_inv<THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)
#----------------------------------
#arcticThinCirrusSecondaryTest sea night ice
#----------------------------------
args = {'title': "maybe_less_generous_arcticthinCirrusSecondaryTest_All_SeaNightIce",
        'xlable': 'T37text',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']
THRESHOLD1 = OFFSETS['OFFSET_T11T12_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = OFFSETS['OFFSET_T11T12_ARCTIC']
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2, show=True)
args['title']= "arcticthinCirrusSecondaryTest_All_SeaNightIce_original"
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']
THRESHOLD1 = OFFSETS['OFFSET_T11T12_ARCTIC'] + OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = OFFSETS['OFFSET_T11T12_ARCTIC']
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2, show=True)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)
#----------------------------------
#extra waterCloudTest sea night
#----------------------------------
args = {'title': "removed_waterCloudTest_All_SeaNightIce_tsur",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC_EXTRA']
if RunWithOutMargins:
    THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC_EXTRA']
TestOk = t11_t37_minus_threshold>THRESHOLD
plot_test(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)




####################################
# SEA DAY
####################################
TestOkAll=None
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS)
TestOk = np.logical_and(isCloudtypeSnowOrIce,np.logical_or(isCloudyWaterDay,
                                                            isClearWaterDay))
print_stats_snow(cloudtype, isCloudyWaterDay, isClearWaterDay, TestOk)
print "***snow test, not implemented as figure, remove snow/ice pixels***"
isCloudyWaterDay = np.logical_and(isCloudyWaterDay, np.not_equal(isCloudtypeSnowOrIce,True))
isClearWaterDay = np.logical_and(isClearWaterDay, np.not_equal(isCloudtypeSnowOrIce,True))
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS)
####################################
# VERSION 2014 COLD_CLOUD_TEST
####################################
args = {'title': "coldCloudTest_v2014_All_SeaDayNoIce",
      'xlable': 'surftemp',#'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': surftemp,#t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,
}
offset = (-5 + dt11tsur_offset_extra * 
    (1 - np.where(t11t37_div>1.0,1.0,t11t37_div)) * 
    (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
THRESHOLD = offset
#t11_t37 should be inv
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        t11_t37_minus_threshold_inv<-1)
TSURF_THRSHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
TestOk = np.logical_and(t11_ts_minus_threshold < THRESHOLD,
                        surftemp>TSURF_THRSHOLD) 
TestOkAll=TestOk
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,np.min(THRESHOLD), show=False)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS)
#----------------------------------
#coldCloudTest sea day
#----------------------------------
args = {'title': "coldCloudTest_All_SeaDayNoIce",
        'xlable': 'Latitude',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']#
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY']
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] 

TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll=TestOk
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,THRESHOLD)
args['title'] = "coldCloudTest_All_SeaDayNoIce_without_Tsur_limit"
THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'] + OFFSETS['QUALITY_MARGIN_T11TSUR']#-13.0-1.0
TSUR_THRESHOLD=0.0
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'] 
    TSUR_THRESHOLD=0.0
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)
#----------------------------------
#coldBrightCloudTest sea day
#----------------------------------
r06_minus_threshold_opaque =  r06 - 100.0*thr_r06*OFFSETS['R06_GAIN_SEA_OPAQUE']
r06_minus_threshold_opaque = np.ma.array(
    r06_minus_threshold_opaque, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
print len(r06_minus_threshold_opaque[r06_minus_threshold_opaque.mask==False])
print len(t11_ts_minus_threshold)
args = {'title': "coldBrightCloudTest_All_SeaDayNoIce",
        'xlable': 'R06 minus dynamic threshold',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': r06_minus_threshold_opaque,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']#
#100.0*threshold->r06*throff->sm_acmg_r06_gain_sea_opaque +
#      throff->sm_acmg_r06_offset_sea_opaque;
#100.0*threshold->r06*throff->sm_acmg_r06_gain_sea_opaque +
#      throff->sm_acmg_r06_offset_sea_opaque;
#/* The satellite reflectances are given in percent. Thus
#   multiply tabulated reflectances by 100.0!	*/
THRESHOLD2 = OFFSETS['R06_OFFSET_SEA_OPAQUE'] + OFFSETS['QUALITY_MARGIN_R06']
if RunWithOutMargins:
    HRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] 
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP']
    THRESHOLD2 = OFFSETS['R06_OFFSET_SEA_OPAQUE']
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold<THRESHOLD1,
                                       surftemp>TSUR_THRESHOLD),
                        r06_minus_threshold_opaque>THRESHOLD2)                        
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterDay, isClearWaterDay, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)


r06_minus_threshold_opaque =  r06 - 100.0*thr_r06*OFFSETS['R06_GAIN_SEA_OPAQUE']
r06_minus_threshold_opaque = np.ma.array(
    r06_minus_threshold_opaque, 
    mask=np.logical_or(r06.mask,thr_r06.mask))
print len(r06_minus_threshold_opaque[r06_minus_threshold_opaque.mask==False])
print len(t11_ts_minus_threshold)
args = {'title': "testcoldBrightCloudTest_All_SeaDayNoIce",
        'xlable': 'R06 minus dynamic threshold',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': r06_minus_threshold_opaque,
        'yvector': t11_ts_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']#

#100.0*threshold->r06*throff->sm_acmg_r06_gain_sea_opaque +
#      throff->sm_acmg_r06_offset_sea_opaque;
#100.0*threshold->r06*throff->sm_acmg_r06_gain_sea_opaque +
#      throff->sm_acmg_r06_offset_sea_opaque;
#/* The satellite reflectances are given in percent. Thus
#   multiply tabulated reflectances by 100.0!	*/
THRESHOLD2 = OFFSETS['R06_OFFSET_SEA_OPAQUE'] + OFFSETS['QUALITY_MARGIN_R06']
if RunWithOutMargins:
    THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_DAY'] 
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] 
    THRESHOLD2 = OFFSETS['R06_OFFSET_SEA_OPAQUE']
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD1,
                        r06_minus_threshold_opaque>THRESHOLD2)                        
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterDay, isClearWaterDay, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)

#----------------------------------
#thinCirrusSecondaryTest sea day
#----------------------------------
args = {'title': "latitude_problem_thinCirrusSecondaryTest_All_SeaDayNoIce_lat",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY'] + OFFSETS['QUALITY_MARGIN_T11T12']#-7.0-1.0
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY'] 
TestOk = t11_t12_minus_threshold>THRESHOLD                      
plot_test(args,isCloudyWaterDay, isClearWaterDay, TestOk, THRESHOLD, show=False)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)
#also make plot with only cirrus
args['title'] = "thinCirrusSecondaryTest_ThinCirrus_SeaDayNoIce_lat"
plot_test(args,isCloudyCirrusWaterDay, isClearWaterDay, TestOk, THRESHOLD)
args['xlable'] ='Tsur'
args['xvector'] = surftemp
args['title'] = "thinCirrusSecondaryTest_ThinCirrus_SeaDayNoIce_tsur"
plot_test(args,isCloudyCirrusWaterDay, isClearWaterDay, TestOk, THRESHOLD)
#linear tsur dep threshold
args = {'title': "testThinCirrusSecondaryTest_All_SeaDayNoIce_lat_linearthreshold",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY']+OFFSETS['QUALITY_MARGIN_T11T12']#2.0+0.3
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY']
THRESHOLD_V=THRESHOLD + (surftemp<305)*(surftemp>280)*(surftemp-280)*(THRESHOLD-1.5)/(305-280)+(surftemp>305)*(-1.5)
TestOk = args['yvector']>THRESHOLD_V
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyCirrusWaterDay,isClearWaterDay,TestOk,THRESHOLD)
args['title'] = "testThinCirrusSecondaryTest_All_SeaDayNoIce_lat_linearthreshold_1p5"
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,THRESHOLD)
args = {'title': "testThinCirrusSecondaryTest_All_SeaDayNoIce_lat_linearthreshold_1.0",
        'xlable': 'Latitude',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY']+OFFSETS['QUALITY_MARGIN_T11T12']#2.0+0.3
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY']
THRESHOLD_V=THRESHOLD + (surftemp<305)*(surftemp>280)*(surftemp-280)*(THRESHOLD-1.0)/(305-280)+(surftemp>305)*(-1.0)
TestOk = args['yvector']>THRESHOLD_V
plot_test(args,isCloudyCirrusWaterDay,isClearWaterDay,TestOk,THRESHOLD)
args['title'] = "testThinCirrusSecondaryTest_All_SeaDayNoIce_lat_linearthreshold_10_margin02"
THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY']+0.2
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11T12_OFFSET_SEA_DAY']
THRESHOLD_V=THRESHOLD + (surftemp<305)*(surftemp>280)*(surftemp-280)*(THRESHOLD-1.0)/(305-280)+(surftemp>305)*(-1.0)
TestOk = args['yvector']>THRESHOLD_V
plot_test(args,isCloudyCirrusWaterDay,isClearWaterDay,TestOk,THRESHOLD)

args['title'] = "testThinCirrusSecondaryTest_ThinCirrus_SeaDayNoIce_lat_modified_thr"
args['yvector'] =  t11t12 -0.75*thr_t11t12
TestOk = args['yvector']>THRESHOLD 
args['ylable'] = 'T11-T12 minus 0.75*dynamic threshold'
plot_test(args,isCloudyCirrusWaterDay, isClearWaterDay, TestOk, THRESHOLD)
plot_test(args,isCloudyWaterDay, isClearWaterDay, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)
#----------------------------------
#notexistingCirrustest sea day
#----------------------------------
print "notexistingCirrustest sea day no ice"
args = {'title': "non_existingT11T37test_All_SeaDayNoIce_lat",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0
TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)


'3D'
fig = plt.figure()
ax = fig.gca(projection = '3d')
plt.plot(surftemp[isClearWaterDay == True],
         ciwv[isClearWaterDay == True], 
         t37_t12_minus_threshold[isClearWaterDay == True], 'bx')
plt.plot(surftemp[isCloudyWaterDay == True],
         ciwv[isCloudyWaterDay == True],
         t37_t12_minus_threshold[isCloudyWaterDay == True], 'r.')

#ax.legend()
#ax.set_xlim3d(0, 1)
#ax.set_ylim3d(0, 1)
#ax.set_zlim3d(-3, 10)
#ax.set_zlable('T37-T12')
#ax.set_ylable('Ciwv')
#ax.set_xlable('Tsur')
ax.set_title('T37-T12 minus t37-T12_threshold: \nblue clear, \nred-manga thin cirrus, \nmanga(detected by pps)')
#plt.show()
