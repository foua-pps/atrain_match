#!/usr/bin/env python
# -*- coding: utf-8 -*-


from pps_threshold_functions import (print_stats,
                                     print_stats_snow,
                                     get_clear_and_cloudy_vectors,
                                     get_feature_values_and_thresholds,
                                     keep_combined_ok)
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from pps_cloudtests import(coldCloudTest_v2014,
                           coldWatercloudTest_v2014,
                           thincoldCirrusTest_v2014,
                           coldCloudTest_no_tsurf_lim,
                           coldCloudTest_tsurf_lim,
                           coldWatercloudTest,
                           watercloudTest,
                           thinCirrusSecondaryTest,
                           thinCirrusPrimaryTest,
                           HighcloudTestt85t11land,
                           HighcloudTestt85t11Sea,
                           textureNightTest,
                           newClearWaterBodiesTest,
                           brightCloudTest,
                           brightCloudTest3A,
                           brightCloudTestSea,
                           brightCloudTestNoSunglint3A,
                           sunglintTest3A,
                           coldBrightCloudTest3A,
                           coldBrightCloudTest37,
                           coldBrightCloudTest,
                           thincoldCirrusTest,
                           coldWatercloudTestDay,
                           reflectingCloudTest,
                           reflectingCloudTestSea,
                           pseudo06CloudTest3A,
                           textureIrVisTest,
                           arcticWarmCloudTestSalomon,
                           arcticWaterCloudTest,
                           arcticThinWaterCloudTest,
                           arcticWarmCirrusSecondaryTest,
                           arcticThinCirrusSecondaryTest,
                           arcticThinCirrusPrimaryTest,
                           DriedOutLakesAndRiversTest)

"""Read all matched data and make some plotting of thresholds
"""

from glob import glob
import os.path
import os
import numpy as np
from scipy import histogram
import re
OFFSET_DIR=os.environ.get('SM_CST_DIR')
OFFSET_FILE=os.environ.get('SM_THROFFSETS_NAME')   
#OFFSET_FILE="threshold_offsets_gac.cfg"
OFFSET_FILE = os.path.join(OFFSET_DIR,OFFSET_FILE)

isNPP = True
isGAC_v2012 = False
RunWithOutMargins=True
RunEdited= True
PLOT_DIR = '/local_disk/nina_pps/Pictures/threshold_plots_3/'
#PLOT_DIR = '/home/rustan/Nina/threshold_plots/'
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
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_12_and_half_shemes_18v2014coldtestsadded_sunz_thr_on_snow_tests_texture_and_ice_thr_all/Reshaped_Files/npp/1km/"
        
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit_tuning4shcemes_4/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit_tuning4shcemes_4_newtunings_2mountandseadayandnightinv/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit_tuning4shcemes_4_newtunings_2mountandseadayandnightinv_tunedarctic2_coastwb/Reshaped_Files/npp/1km/"
    	#ROOT_DIR ="thresholds/Reshaped_Files_2014/Reshaped_Files/npp6sh85/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit_tuning4shcemes_4_newtunings_2mountandseadayandnightinv_tunedarctic2_coastwb_newarc_coast3A_4/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit_tuning4shcemes_4_newtunings_2mountandseadayandnightinv_tunedarctic2_coastwb_newarc_coast3A_4_next_arc_1p6/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit_tuning4shcemes_4_newtunings_2mountandseadayandnightinv_tunedarctic2_coastwb_newarc_coast3A_4_next_arc_1p6_newlog/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_coastalzone0percent_4_37cond_test_newlogic_day_Sea_tohre_rogh_limit_tuning4shcemes_4_newtunings_2mountandseadayandnightinv_tunedarctic2_coastwb_newarc_coast3A_4_next_arc_1p6_newlog_errcorr_2/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_including_new_sst_newlogic_t11ts_errcorr_5/Reshaped_Files/npp/1km/"
        #ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_npp_as_metop/Reshaped_Files/npp/1km/"
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_ice_coast_int_ice_scheme_extratests_2/Reshaped_Files/npp/1km/"
        ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20131014_Test_CM_ice_coast_int_ice_scheme_extratests_2_new_emiss_4/Reshaped_Files/npp/1km/"
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




offset_file = open(OFFSET_FILE, 'r')  
re_offset  = re.compile(r"^SM_ACMG_(\w+)[:=]\s*(-*\d+\.*\d*)\D")
OFFSETS={}
OFFSETS['MAX_SUNZEN_DAY'] = 80.0 #from pps_common.h

for line in offset_file:
    match_offset = re_offset.match(line)
    if match_offset:
        print "read_offset %s=%s"%(match_offset.group(1),match_offset.group(2))
        OFFSETS[match_offset.group(1)] = 1.0*float(match_offset.group(2))

caobj = CalipsoAvhrrTrackObject()
for filename in files:
    print os.path.basename(filename)
    caobj = caobj + readCaliopAvhrrMatchObj(filename)
caobj.avhrr.mask_nodata(nodata=-9)
fthr = get_feature_values_and_thresholds(caobj)

cloudObj = get_clear_and_cloudy_vectors(caobj, isACPGv2012, isGAC)

import matplotlib.pyplot as plt
fig = plt.figure(figsize = (9,8))
ax = fig.add_subplot(111)
plt.plot(caobj.avhrr.thr_t11ts, 
         caobj.avhrr.thr_t11ts_inv, 'b.')
t11ts_per = np.percentile(
     caobj.avhrr.thr_t11ts,
     [01, 10, 25, 50, 75, 90, 99])
t11ts_inv_per = np.percentile(
     caobj.avhrr.thr_t11ts_inv,
     [01, 10, 25, 50, 75, 90, 99])
plt.axhline(y = t11ts_inv_per[0], color = 'g',ls = '--')
plt.axhline(y = t11ts_inv_per[2], color = 'g',ls = '--')
plt.axhline(y = t11ts_inv_per[3], color = 'g',linewidth = 2)
plt.axhline(y = t11ts_inv_per[4], color = 'g',ls = '--')
plt.axhline(y = t11ts_inv_per[6], color = 'g',ls = '--')
ax.set_ylabel('t11ts_inv')
ax.set_xlabel('t11ts')
plt.axvline(x = t11ts_per[0], color = 'r',ls = '--')
plt.axvline(x = t11ts_per[2], color = 'r',ls = '--')
plt.axvline(x = t11ts_per[3], color = 'r',linewidth = 2)
plt.axvline(x = t11ts_per[4], color = 'r',ls = '--')
plt.axvline(x = t11ts_per[6], color = 'r',ls = '--')
plt.show()

TestOkAll = None
args={'PLOT_DIR': PLOT_DIR,
      'SATELLITE': SATELLITE,
      'USE_MARGINS':False
      }
LittleLight = fthr.sunz>OFFSETS['MAX_SUNZEN_TWILIGHT_VIS']
EnoughLight = fthr.sunz<OFFSETS['MAX_SUNZEN_TWILIGHT_VIS']
TempAboveIce = fthr.surftemp>OFFSETS['COLDEST_SEASURFACE_TEMP']
TempAbove230 = fthr.surftemp>230
TempBelowMinus5 = fthr.surftemp<268.15
TempBelow10C = fthr.surftemp<283.15
TempBelow20C = fthr.surftemp<293.15
TempBelow30C = fthr.surftemp<303.15
NoChannel3_7 = caobj.avhrr.all_arrays['bt37micron']<=-9
TempAboveColdestSea = fthr.surftemp>OFFSETS['COLDEST_SEASURFACE_TEMP']
DiffT1137Above1K = np.logical_or(fthr.t11_t37_minus_threshold_inv<-1.0,
                                 fthr.t11_t37_minus_threshold>1.0)
t37Above240 = caobj.avhrr.all_arrays['bt37micron']>240


###################################
#ALL CLOUDS
###################################
SchemeName = 'All'
TestOkAll = coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   show=False, TestOkAll=TestOkAll)
print_stats(SchemeName, caobj, cloudObj,args,TestOkAll)
TestOk = coldWatercloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
print_stats(SchemeName, caobj, cloudObj, args, TestOkAll)
TestOk = thincoldCirrusTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
print_stats(SchemeName, caobj, cloudObj, args, TestOkAll)
TestOkClear = DriedOutLakesAndRiversTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",show=True)
if fthr.qr16r06 is not None:
     TestOkClear = newClearWaterBodiesTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True)
###################################################
#Testing stuff
###################################################
if fthr.qr16r06 is not None:
     coldBrightCloudTest3A('LandDay', caobj, cloudObj, fthr, OFFSETS, args,   
                           info="extra_",   show=False, TestOkAll=TestOkAll)
     brightCloudTestNoSunglint3A('LandDay', caobj, cloudObj, fthr, OFFSETS, args,   
                                 info="extra_",  show=True, TestOkAll=TestOkAll)
     brightCloudTest3A('LandDay', caobj, cloudObj, fthr, OFFSETS, args,   
                       info="extra_",   ExtraCond=TempBelow30C, show=True, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = sunglintTest3A('CoastDay', caobj, cloudObj, fthr, 
                                OFFSETS, args,  info="",   show=False, TestOkAll=TestOkAll)
     TestOkAll = sunglintTest3A('CoastTwilight', caobj, cloudObj, fthr, 
                                OFFSETS, args,  info="",   show=False, TestOkAll=TestOkAll)
     TestOkAll = sunglintTest3A('CoastTwilightInv', caobj, cloudObj, fthr, 
                                OFFSETS, args,  info="",   show=False, TestOkAll=TestOkAll)
brightCloudTest('CoastDay', caobj, cloudObj, fthr, OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
brightCloudTest('WaterDay', caobj, cloudObj, fthr, OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
brightCloudTest('IceDay', caobj, cloudObj, fthr, OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
brightCloudTest('WaterTwilight', caobj, cloudObj, fthr, OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
watercloudTest('LandNightInv', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="test_",show=False, TestOkAll=TestOkAll)
watercloudTest('LandNightInv', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_",NEW_THRESHOLD=0, ExtraCond=TempBelowMinus5, show=False, TestOkAll=TestOkAll)
watercloudTest('All', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="all_",  show=False, TestOkAll=TestOkAll)

NEW_THRESHOLD = np.where(caobj.avhrr.all_arrays['bt37micron'] <263.15, 0.0, (caobj.avhrr.all_arrays['bt37micron']-263.15)/5.0*1.0)
NEW_THRESHOLD = np.where(caobj.avhrr.all_arrays['bt37micron'] <268.15, NEW_THRESHOLD, 1+(caobj.avhrr.all_arrays['bt37micron']-268.15)/15.0*1.0)
NEW_THRESHOLD = np.where(caobj.avhrr.all_arrays['bt37micron'] <283.15, NEW_THRESHOLD, 2+(caobj.avhrr.all_arrays['bt37micron']-268.15)/5.0*1.0)
NEW_THRESHOLD = np.where(caobj.avhrr.all_arrays['bt37micron'] >293.15, 3.0,  NEW_THRESHOLD)
watercloudTest('LandNightInv', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_",NEW_THRESHOLD=NEW_THRESHOLD,ExtraCond=t37Above240,  show=False, TestOkAll=TestOkAll)
NEW_THRESHOLD = np.where(fthr.surftemp <263.15, 0.0, (fthr.surftemp-263.15)/30.0*3.0)
#NEW_THRESHOLD = np.where(fthr.surftemp <268.15, NEW_THRESHOLD, 1+(fthr.surftemp-268.15)/15.0*1.0)
#NEW_THRESHOLD = np.where(fthr.surftemp <283.15, NEW_THRESHOLD, 2+(fthr.surftemp-268.15)/5.0*1.0)
NEW_THRESHOLD = np.where(fthr.surftemp >293.15, 3.0,  NEW_THRESHOLD)
watercloudTest('LandNightInv', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_",NEW_THRESHOLD=NEW_THRESHOLD,  show=False, TestOkAll=TestOkAll)

NEW_THRESHOLD = np.where(fthr.surftemp <268.15, 0.0, (fthr.surftemp-268.15)/15.0*1.0)
NEW_THRESHOLD = np.where(fthr.surftemp >283.15, 1.0,  NEW_THRESHOLD)
watercloudTest('LandNight', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_",NEW_THRESHOLD=NEW_THRESHOLD,  show=False, TestOkAll=TestOkAll)
watercloudTest('LandNightMount', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_", ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)

NEW_THRESHOLD = np.where(fthr.surftemp <263.15, 0.0, (fthr.surftemp-263.15)/5.0*1.0)
NEW_THRESHOLD = np.where(fthr.surftemp <268.15, NEW_THRESHOLD, 1+(fthr.surftemp-268.15)/15.0*1.0)
NEW_THRESHOLD = np.where(fthr.surftemp <283.15, NEW_THRESHOLD, 2+(fthr.surftemp-268.15)/5.0*1.0)
NEW_THRESHOLD = np.where(fthr.surftemp >293.15, 3.0,  NEW_THRESHOLD)
watercloudTest('LandNightMount', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_", ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)
watercloudTest('LandNightMount', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_",NEW_THRESHOLD=NEW_THRESHOLD, ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)
watercloudTest('CoastNight', caobj, cloudObj, fthr, OFFSETS, args, 
                        info="extra_overcoldsurfaces_",NEW_THRESHOLD=NEW_THRESHOLD, ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)




if fthr.qr16r06 is not None:
     TestOk = pseudo06CloudTest3A('LandDay', caobj, cloudObj, fthr, OFFSETS, args,   
                                  info="non_existing_", show=False, TestOkAll=TestOkAll)

TestOkClear = DriedOutLakesAndRiversTest('LandDay', caobj, cloudObj, fthr, OFFSETS, args, info="",show=True)

####################################
# DAY COAST MOUNTAIN
#################################
SchemeName = 'CoastDayMount'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.CoastDayMount = np.logical_and(cloudObj.isClear.CoastDayMount,
                                                np.equal(TestOkSnow,False))
cloudObj.isCloudy.CoastDayMount = np.logical_and(cloudObj.isCloudy.CoastDayMount,
                                                 np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
TestOkClear = DriedOutLakesAndRiversTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",show=True)
cloudObj.isClear.CoastDayMount = np.logical_and(cloudObj.isClear.CoastDayMount,
                                                np.equal(TestOkClear,False))
cloudObj.isCloudy.CoastDayMount = np.logical_and(cloudObj.isCloudy.CoastDayMount,
                                                 np.equal(TestOkClear,False))
if fthr.qr16r06 is not None:
     TestOkClear = newClearWaterBodiesTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True)
     cloudObj.isClear.CoastDayMount = np.logical_and(cloudObj.isClear.CoastDayMount,
                                                     np.equal(TestOkClear,False))
     cloudObj.isCloudy.CoastDayMount = np.logical_and(cloudObj.isCloudy.CoastDayMount,
                                                      np.equal(TestOkClear,False))

print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll)
TestOkAll = coldBrightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",  show=True, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTestDay(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   show=True, TestOkAll=TestOkAll)
TestOkAll = brightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",  show=True, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     DontUse = brightCloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",  show=True, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     DontUse = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, ADD_TO_ROG_THR=10.0, info="non_existing_", show=True, TestOkAll=TestOkAll)
TestOkAll = coldBrightCloudTest37(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  info="",   show=True, TestOkAll=TestOkAll)
TestOkAll = thincoldCirrusTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll)
print "not implemented coldclouds in sunglint yet  cloudsInSunglint"
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  info="", NEW_THRESHOLD=None, onlyCirrus=False, show=True, TestOkAll=TestOkAll)


####################################
# DAY SEA NO ICE
####################################
SchemeName = 'WaterDay'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.WaterDay = np.logical_and(cloudObj.isClear.WaterDay,
                                               np.equal(TestOkSnow, False))
cloudObj.isCloudy.WaterDay = np.logical_and(cloudObj.isCloudy.WaterDay,
                                                np.equal(TestOkSnow, False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
TestOkClear = DriedOutLakesAndRiversTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",show=True, ExtraCond=None)#cloudObj.isPpsClear)

cloudObj.isClear.WaterDay = np.logical_and(cloudObj.isClear.WaterDay,
                                               np.equal(TestOkClear, False))
cloudObj.isCloudy.WaterDay = np.logical_and(cloudObj.isCloudy.WaterDay,
                                                np.equal(TestOkClear, False))
if fthr.qr16r06 is not None:
     TestOkClear = newClearWaterBodiesTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True)
     cloudObj.isClear.WaterDay = np.logical_and(cloudObj.isClear.WaterDay,
                                               np.equal(TestOkClear, False))
     cloudObj.isCloudy.WaterDay = np.logical_and(cloudObj.isCloudy.WaterDay,
                                                np.equal(TestOkClear, False))
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="non_exisitng_", show=False, TestOkAll=TestOkAll)#, ExtraCond=DiffT1137Above1K)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                                    args, show=False, TestOkAll=TestOkAll)                                       
TestOkAll = coldBrightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                ExtraCond=EnoughLight, show=False, TestOkAll=TestOkAll)
print "pixels with sunglint risk have higher thresholds"
TestOkAll = brightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="use_opaquer06_thr_",  show=False, TestOkAll=TestOkAll)
print "should theis be coldWatercloudTestDay??"
DontUseIt = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                               info="removed_", show=False, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTestDay(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                               info="add_non_existing", show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = brightCloudTest3A('WaterDay', caobj, cloudObj, fthr, OFFSETS, args,   
                                   info="extra_",  show=False, ExtraCond=NoChannel3_7, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args, info="notusedif37_",  
                                 ExtraCond=NoChannel3_7, show=False, TestOkAll=TestOkAll)
     DontUseIt = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args,   
                                 info="maybe_use_notusedif37_",  show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = textureIrVisTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                          info="less_generous_and_maybe_seasurftemp_thr_", show=False, TestOkAll=TestOkAll)
if not isGAC:
     textureIrVisTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                      info="original_less_generous_and_", 
                      REMOVE_2times_THR=True, show=False, TestOkAll=TestOkAll)
print "not implemented texturevistest"
TestOkAll = textureNightTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                          info="less_generous_and_seasurftemp_thr_", show=False, TestOkAll=TestOkAll)
if not isGAC:
     TestOkAll=textureNightTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                      info="original_less_generous_and_seasurftemp_thr_", 
                      REMOVE_2times_THR=True, show=False, TestOkAll=TestOkAll)
TestOkAll = brightCloudTestSea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                    info="less_generous", show=False, TestOkAll=TestOkAll)
Orgiginal = brightCloudTestSea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                            info="original_less_generous", 
                            OTHER_T11_THR=OFFSETS['T11_SEA_MIN'], show=False, TestOkAll=TestOkAll)

print "**not implemented brightCloudOverWaterTest"
print "**not implemented sstDaytimeTest"
print "*** not implemented cirrusCloudOverWaterTest**"
coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                                    args, info="non_existing_",  show=False, TestOkAll=TestOkAll)                 
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11Sea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                     show=False, TestOkAll=TestOkAll)
###################################
# TWILIGHT SEA NO ICE 
####################################
SchemeName = 'WaterTwilight'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.WaterTwilight = np.logical_and(cloudObj.isClear.WaterTwilight,
                                          np.equal(TestOkSnow,False))
cloudObj.isCloudy.WaterTwilight = np.logical_and(cloudObj.isCloudy.WaterTwilight,
                                           np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                              args, show=True, TestOkAll=TestOkAll)#, ExtraCond=DiffT1137Above1K)
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                                 info="strange_dynamic_thr_add_slope_nina_", show=False, TestOkAll=TestOkAll)
TestOkAll = textureIrVisTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                          info="add_margins_ordianrynotsunglintthr_nina", 
                          show=False, TestOkAll=TestOkAll)
TestOkAll = textureNightTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                          info="add2margins", show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                                    args, show=False, TestOkAll=TestOkAll)                                       
TestOkAll = reflectingCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                             info="offset_slope_", show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     DontUseIt = pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  info="notudedif37_", ExtraCond=NoChannel3_7, show=False, TestOkAll=TestOkAll)
     TestOkAll = pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                         info="maybe_use_alsoif37_", show=False, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                       show=False, TestOkAll=TestOkAll)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                        info="even_less_generous_added_security_offset_",show=False, TestOkAll=TestOkAll)
watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                        info="original_added_security_offset_",NEW_THRESHOLD=0,show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,
                                 show=False, TestOkAll=TestOkAll)  
TestOkAll = reflectingCloudTestSea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                info="less_generous_", show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                               info="non_existing_",  show=False, TestOkAll=TestOkAll)
#nodifference: brightCloudTest('WaterTwilight', caobj, cloudObj, fthr, OFFSETS, args, info="non_existing_",  show=False, TestOkAll=TestOkAll)
#add original reflectingCloudTestSea
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11Sea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                     show=False, TestOkAll=TestOkAll)

###################################
# DAY COAST 
#################################
SchemeName = 'CoastDay'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.CoastDay = np.logical_and(cloudObj.isClear.CoastDay,
                                           np.equal(TestOkSnow,False))
cloudObj.isCloudy.CoastDay = np.logical_and(cloudObj.isCloudy.CoastDay,
                                            np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
if fthr.qr16r06 is not None:
     TestOkClear = newClearWaterBodiesTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True)
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTestDay(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                               args,   show=False, TestOkAll=TestOkAll)
TestOkAll = brightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     DontUse = brightCloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
     TestOkAll = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args, ADD_TO_ROG_THR=10.0, info="non_existing_",  
                                 show=False, TestOkAll=TestOkAll)
TestOkAll = coldBrightCloudTest37(SchemeName, caobj, cloudObj, fthr, 
                               OFFSETS, args,  info="",   show=False, TestOkAll=TestOkAll)
TestOkAll = thincoldCirrusTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                   show=False, TestOkAll=TestOkAll)
print "***coldCloudsInSunglint, cloudsInSunglint3A***"
TestOkAll = coldBrightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                    info="nonexisting", show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = sunglintTest3A(SchemeName, caobj, cloudObj, fthr, 
                                OFFSETS, args,  info="",   show=False, TestOkAll=TestOkAll)
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, 
                                      OFFSETS, args, show=False, TestOkAll=TestOkAll)
####################################
# TWILIGHT LAND
####################################
SchemeName = 'LandTwilight'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.LandTwilight = np.logical_and(cloudObj.isClear.LandTwilight,
                                          np.equal(TestOkSnow,False))
cloudObj.isCloudy.LandTwilight = np.logical_and(cloudObj.isCloudy.LandTwilight,
                                           np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                              args, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                                    args, show=False, TestOkAll=TestOkAll)                                       
TestOkAll = reflectingCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  info="notudedif37_", ExtraCond=NoChannel3_7, show=False, TestOkAll=TestOkAll)
     pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                         info="maybe_use_alsoif37_",  show=False, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTestDay(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                       show=False, TestOkAll=TestOkAll)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                                 info="maybe_add_tsur_limit_less_generous_also", 
                                  show=False, TestOkAll=TestOkAll)        
tmpthr=np.min([OFFSETS['T11_OFFSET_LAND_DAY'], 
               OFFSETS['T11_OFFSET_LAND_NIGHT']])
coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,
                                 info="original_less_generous_also", 
                                  NEW_THRESHOLD=tmpthr, show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                                 info="maybe_offset_slope_",  show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                               info="maybe_offset_slope_",  show=False, TestOkAll=TestOkAll)
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, 
                                      OFFSETS, args, show=False, TestOkAll=TestOkAll)
###################################
# DAY LAND 
####################################
SchemeName = 'LandDay'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.LandDay = np.logical_and(cloudObj.isClear.LandDay,
                                          np.equal(TestOkSnow,False))
cloudObj.isCloudy.LandDay = np.logical_and(cloudObj.isCloudy.LandDay,
                                           np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
if fthr.qr16r06 is not None:
     TestOkClear = newClearWaterBodiesTest(SchemeName, caobj, cloudObj, fthr, 
                                           OFFSETS, args, show=False)
     cloudObj.isClear.LandDay = np.logical_and(cloudObj.isClear.LandDay,
                                               np.equal(TestOkClear,False))
     cloudObj.isCloudy.LandDay = np.logical_and(cloudObj.isCloudy.LandDay,
                                                np.equal(TestOkClear,False))
     print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
TestOkAll = brightCloudTest(SchemeName, caobj, cloudObj, fthr, 
                        OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args, info="notusedif37_",  
                                 ExtraCond=NoChannel3_7, show=False, TestOkAll=TestOkAll)
     brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args,# ADD_TO_ROG_THR=10,   
                                 info="usealsoif37_notusedif37_",  show=False, TestOkAll=TestOkAll)
TestOkAll = coldBrightCloudTest37(SchemeName, caobj, cloudObj, fthr, 
                               OFFSETS, args,  info="",   show=False, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTestDay(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                                  show=False, TestOkAll=TestOkAll)
TestOkAll = coldBrightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,
                    ExtraCond=EnoughLight, show=False, TestOkAll=TestOkAll)
TestOkAll = thincoldCirrusTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                   show=False, TestOkAll=TestOkAll)
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, 
                                      OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                              args, show=False, TestOkAll=TestOkAll)


###################################
# TWILIGHT COAST INVERSION
#################################
SchemeName = 'CoastTwilightInv'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.CoastTwilightInv = np.logical_and(
     cloudObj.isClear.CoastTwilightInv,np.equal(TestOkSnow,False))
cloudObj.isCloudy.CoastTwilightInv = np.logical_and(
     cloudObj.isCloudy.CoastTwilightInv,np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                                 args, ExtraCond=TempAbove230,show=False, TestOkAll=TestOkAll)
###################################
# TWILIGHT LAND INVERSION
#################################
SchemeName = 'LandTwilightInv'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.LandTwilightInv = np.logical_and(
     cloudObj.isClear.LandTwilightInv,np.equal(TestOkSnow,False))
cloudObj.isCloudy.LandTwilightInv = np.logical_and(
     cloudObj.isCloudy.LandTwilightInv,np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                                 args, ExtraCond=TempAbove230,show=False, TestOkAll=TestOkAll)
###################################
# TWILIGHT COAST 
#################################
SchemeName = 'CoastTwilight'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.CoastTwilight = np.logical_and(cloudObj.isClear.CoastTwilight,
                                                np.equal(TestOkSnow,False))
cloudObj.isCloudy.CoastTwilight = np.logical_and(cloudObj.isCloudy.CoastTwilight,
                                                 np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = thincoldCirrusTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                   info="added_", show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="offset_slope_",  show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="offset_slope_",  show=False, TestOkAll=TestOkAll)

####################################
# TWILIGHT LAND MOUNTAIN 
####################################
SchemeName = 'LandTwilightMount'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.LandTwilightMount = np.logical_and(cloudObj.isClear.LandTwilightMount,
                                                    np.equal(TestOkSnow,False))
cloudObj.isCloudy.LandTwilightMount = np.logical_and(cloudObj.isCloudy.LandTwilightMount,
                                                     np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
TestOkAll = coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = thincoldCirrusTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                   info="non_exisiting_", show=False, TestOkAll=TestOkAll)
####################################
# NIGHT COAST
####################################
SchemeName = 'CoastNight'
TestOkAll=None
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
####################################
# NIGHT LAND 
####################################
SchemeName = 'LandNight'
TestOkAll=None
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                                    args, show=False, TestOkAll=TestOkAll)                                       
#TestOkAll = coldCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=False, TestOkAll=TestOkAll)  
print "coldWatercloudTest: limit ts for desert is 7deg lower"
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                        info="added_security_offset_",show=False, TestOkAll=TestOkAll)
watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                        info="original_added_security_offset_",NEW_THRESHOLD=0,show=False, TestOkAll=TestOkAll)
print "only for desert can not sort out those pixels!"
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                               info="desert_only_",  show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                               info="more_generous_",  show=False, TestOkAll=TestOkAll)
thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="original_more_generous_", 
                      NEW_THRESHOLD=OFFSETS['T37T12_OFFSET_LAND_NIGHT'],  show=False, TestOkAll=TestOkAll)
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  info="", NEW_THRESHOLD=None, onlyCirrus=False, show=False, TestOkAll=TestOkAll)

####################################
# NIGHT LAND MOUNTAIN
####################################
SchemeName = 'LandNightMount'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                        info="", ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)
watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                        info="removed_",NEW_THRESHOLD=0, ExtraCond=TempAbove230,show=False, TestOkAll=TestOkAll)
print "limit ts for desert is 3deg lower"
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                            ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS,
                                    args, info="more_generous_opaque_", 
                                    ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)  
coldCloudTest_no_tsurf_lim(
     SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
     info="orginal_more_generous_opaque_", 
     NEW_THRESHOLD=OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT_OPAQUE'], 
     ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll) 
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, 
                               OFFSETS, args,   
                               info="more_generous_",  show=False, TestOkAll=TestOkAll)
thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                      info="original_more_generous_", 
                      NEW_THRESHOLD=OFFSETS['T37T12_OFFSET_LAND_NIGHT'],  
                      ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)
print "only for desert can not sort out those pixels!"
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                               info="desert_only_",  show=False, TestOkAll=TestOkAll)
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, 
                              args, show=False, TestOkAll=TestOkAll)
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)
####################################
# NIGHT COAST INVERSION
####################################
SchemeName = 'CoastNightInv'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
Original = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                         ExtraCond=fthr.latitude>55,info="latitude_limit_", NEW_THRESHOLD=0, show=False, TestOkAll=TestOkAll)
print "not done if strong inversion, results might differ a bit"
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, NEW_THRESHOLD=0.0, info="removed_without_security_thr_", show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                               info="more_generous_",  show=False, TestOkAll=TestOkAll)
thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="original_more_generous_", 
                      NEW_THRESHOLD=np.min([OFFSETS['T37T12_OFFSET_SEA_NIGHT'],OFFSETS['T37T12_OFFSET_LAND_NIGHT']]),  show=False, TestOkAll=TestOkAll)
print "not done if strong inversion, results might differ a bit"
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)  
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  show=False, TestOkAll=TestOkAll)
####################################
# NIGHT LAND INVERSION
####################################
SchemeName = 'LandNightInv'
TestOkAll=None
TestOkAll =  coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="removed_",
                        NEW_THRESHOLD=OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'], show=False, TestOkAll=TestOkAll)
TestOkAll = keep_combined_ok(TestOk,TestOkAll)
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
print "not done if strong inversion"
print "coldWatercloudTest: limit ts for desert is lower ?????CHECK???"
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
print "this test is not done over desert, do over desert??"
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="less_generous_",
                         show=False, TestOkAll=TestOkAll)
watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, 
                         ExtraCond=fthr.latitude>55,info="latitude_limit_", NEW_THRESHOLD=0, show=False, TestOkAll=TestOkAll)
print "only for desert can not sort out those pixels!"
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                               info="desert_only_",  show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                               info="more_generous_",  show=False, TestOkAll=TestOkAll)
thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="original_more_generous_", 
                      NEW_THRESHOLD=OFFSETS['T37T12_OFFSET_LAND_NIGHT'],  show=False, TestOkAll=TestOkAll)
print "this is not done over desert and not if strong inversion"
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)  
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  info="", NEW_THRESHOLD=None, onlyCirrus=False, show=False, TestOkAll=TestOkAll)
####################################
# SEA NIGHT
####################################
SchemeName = 'WaterNight'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkAll = coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)#, ExtraCond=DiffT1137Above1K)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="latitude_problem",  show=False, TestOkAll=TestOkAll)
thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="latitude_problem", onlyCirrus=True,  show=False, TestOkAll=TestOkAll)
TestOkAll = textureNightTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                          info="less_generous_t37t12text_", show=False, TestOkAll=TestOkAll)
textureNightTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                          info="original_less_generous_t37t12text_", NEW_THRESHOLD_T37T12=OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT']-0.15, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)  
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="latitude_problem", show=False, TestOkAll=TestOkAll)  
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11Sea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                  info="", NEW_THRESHOLD=None, onlyCirrus=False, show=False, TestOkAll=TestOkAll)
print "watercloudOverWaterTest"
print "*** need warmest neighbour to plot this test***"

####################################
# DAY LAND MOUNTAIN 
####################################
SchemeName = 'LandDayMount'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.LandDayMount = np.logical_and(cloudObj.isClear.LandDayMount,
                                               np.equal(TestOkSnow,False))
cloudObj.isCloudy.LandDayMount = np.logical_and(cloudObj.isCloudy.LandDayMount,
                                                np.equal(TestOkSnow,False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll) 
if fthr.qr16r06 is not None:
     TestOkClear = newClearWaterBodiesTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False)
     cloudObj.isClear.LandDayMount = np.logical_and(cloudObj.isClear.LandDayMount,
                                                    np.equal(TestOkClear,False))
     cloudObj.isCloudy.LandDayMount = np.logical_and(cloudObj.isCloudy.LandDayMount,
                                                     np.equal(TestOkClear,False))
     print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)  
TestOkAll = coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = brightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="",  show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args, info="notusedif37_",  
                                 ExtraCond=NoChannel3_7, show=False, TestOkAll=TestOkAll)
     brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args,   
                                 info="dont_notusedif37_",  show=False, TestOkAll=TestOkAll)
TestOkAll = coldBrightCloudTest37(SchemeName, caobj, cloudObj, fthr, 
                               OFFSETS, args,  info="",   show=False, TestOkAll=TestOkAll)
TestOkAll = coldBrightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                    ExtraCond=EnoughLight, show=False, TestOkAll=TestOkAll)
TestOkAll = thincoldCirrusTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                   show=False, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTestDay(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  
                      ExtraCond=TempAbove230, show=False, TestOkAll=TestOkAll)
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11land(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                      show=False, TestOkAll=TestOkAll)

####################################
# ICE DAY
####################################
SchemeName = 'IceDay'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.IceDay = np.logical_and(cloudObj.isClear.IceDay,
                                               np.equal(TestOkSnow, False))
cloudObj.isCloudy.IceDay = np.logical_and(cloudObj.isCloudy.IceDay,
                                                np.equal(TestOkSnow, False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll) 
TestOkAll = coldBrightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll)
TestOkAll = brightCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="use_opaquer06_thr_",  show=True, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     TestOkAll = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args, info="notusedif37_",  
                                 ExtraCond=NoChannel3_7, show=True, TestOkAll=TestOkAll)
     DontUseIt = brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, fthr, 
                                 OFFSETS, args,   
                                 info="maybe_use_notusedif37_",  show=True, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll) 
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=True, TestOkAll=TestOkAll)
TestOkAll = brightCloudTestSea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=True, TestOkAll=TestOkAll)

TestOkAll = textureIrVisTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll, ExtraCond=np.logical_and(EnoughLight,TempAboveIce))
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll)
TestOkAll =coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=True, TestOkAll=TestOkAll)

 
TestOkAll = reflectingCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=True, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     DontUseIt = pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, ExtraCond=NoChannel3_7, show=True, TestOkAll=TestOkAll)
     pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="also_for37_", show=True, TestOkAll=TestOkAll)

TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=True, TestOkAll=TestOkAll)

TestOkAll = arcticThinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="extra_", show=True, TestOkAll=TestOkAll)


####################################
# ICE TWILIGHT
####################################
SchemeName = 'IceTwilight'
TestOkAll=None
print_stats(SchemeName,caobj, cloudObj, args, TestOkAll)
TestOkSnow = print_stats_snow(SchemeName,caobj, cloudObj, args)
cloudObj.isClear.IceTwilight = np.logical_and(cloudObj.isClear.IceTwilight,
                                               np.equal(TestOkSnow, False))
cloudObj.isCloudy.IceTwilight = np.logical_and(cloudObj.isCloudy.IceTwilight,
                                                np.equal(TestOkSnow, False))
print "***snow test, not implemented as figure, remove snow/ice pixels***"
TestOkAll = thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=False, TestOkAll=TestOkAll)
TestOkAll = textureIrVisTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll, ExtraCond=np.logical_and(EnoughLight,TempAboveIce))
TestOkAll = thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)  
TestOkAll = reflectingCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=False, TestOkAll=TestOkAll)
if fthr.qr16r06 is not None:
     DontUseIt = pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, ExtraCond=NoChannel3_7, show=False, TestOkAll=TestOkAll)
     pseudo06CloudTest3A(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="also_for37_", show=False, TestOkAll=TestOkAll)
TestOkAll = coldWatercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = watercloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll) 
TestOkAll = coldCloudTest_v2014(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,  show=False, TestOkAll=TestOkAll)#, ExtraCond=DiffT1137Above1K)
TestOkAll = arcticThinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, info="extra_", show=False, TestOkAll=TestOkAll)
TestOkAll = arcticWarmCloudTestSalomon(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)

####################################
# ICE NIGHT
####################################
SchemeName = 'IceNight'
TestOkAll=None
TestOkAll = arcticWaterCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = arcticWaterCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll, LessSecureText=True)
TestOkAll = coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)  
TestOkAll = arcticThinCirrusPrimaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = arcticThinWaterCloudTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = arcticWarmCloudTestSalomon(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = arcticWarmCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
TestOkAll = arcticThinCirrusSecondaryTest(SchemeName, caobj, cloudObj, fthr, OFFSETS, args, show=False, TestOkAll=TestOkAll)
if fthr.t85_t11_minus_threshold is not None:
     TestOkAll = HighcloudTestt85t11Sea(SchemeName, caobj, cloudObj, fthr, OFFSETS, args,   
                                     show=False, TestOkAll=TestOkAll)
