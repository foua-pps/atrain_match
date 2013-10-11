#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
OFFESET_FILE = os.path.join(OFFSET_DIR,OFFSET_FILE)
isNPP = False
isGAC_v2012 = False
RunWithOutMargins=True
if isNPP:
    isACPGv2012=False
    ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20130908_TEST1B_2/Reshaped_Files/npp/1km/"
    files = glob(ROOT_DIR + "/????/??/*/*h5")
    SATELLITE='npp'
elif isGAC_v2012:
    isACPGv2012=True
    ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_GAC_2012_el2010/Reshaped_Files/noaa18/5km/"
    files = glob(ROOT_DIR + "/????/??/cea*/*noaa18*h5")
    SATELLITE='avhrr18'
else:
    isACPGv2012=False
    ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_GAC_2014/Reshaped_Files/noaa18/5km/"
    files = glob(ROOT_DIR + "/????/??/cea*/*noaa18*h5")
    SATELLITE='avhrr18_v2014'

print OFFESET_FILE 


  #print "hi", filename
offset_file = open(OFFESET_FILE, 'r')  
re_offset  = re.compile(r"^SM_ACMG_(\w+)[:=]\s*(-*\d+\.*\d*)\D")
OFFSETS={}
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
    get_ice_info_pps2012)
caObj = CalipsoAvhrrTrackObject()
for filename in files:
    print os.path.basename(filename)
    caObj = caObj + readCaliopAvhrrMatchObj(filename)

isClear = caObj.calipso.all_arrays['number_of_layers_found'] == 0
isCloudy = caObj.calipso.all_arrays['number_of_layers_found'] >0
isSingleLayerCloud = caObj.calipso.all_arrays['number_of_layers_found'] == 1
SeesThrough = caObj.calipso.all_arrays['lidar_surface_elevation'][0] >-999
isHigh = caObj.calipso.all_arrays['cloud_top_profile'][0] >5.0
isOpticallyThin = caObj.calipso.all_arrays['optical_depth_top_layer5km']<1.0
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
thr_t11ts_inv = caObj.avhrr.all_arrays['thr_t11ts_inv']
thr_t11t12_inv = caObj.avhrr.all_arrays['thr_t11t12_inv']
thr_t37t12_inv = caObj.avhrr.all_arrays['thr_t37t12_inv']
thr_r06 = caObj.avhrr.all_arrays['thr_r06']
surftemp = caObj.avhrr.all_arrays['surftemp']
latitude = caObj.avhrr.all_arrays['latitude']
ciwv = caObj.avhrr.all_arrays['ciwv']
t11 = caObj.avhrr.all_arrays['bt11micron']
t12 = caObj.avhrr.all_arrays['bt12micron']
t37 = caObj.avhrr.all_arrays['bt37micron']
r06 = caObj.avhrr.all_arrays['r06micron']
t11text = caObj.avhrr.all_arrays['text_t11']
t37t12text = caObj.avhrr.all_arrays['text_t37t12']
t37text = caObj.avhrr.all_arrays['text_t37']

cloudtype_conditions = caObj.avhrr.all_arrays['cloudtype_conditions']
cloudtype_status = caObj.avhrr.all_arrays['cloudtype_status']
cloudtype = caObj.avhrr.all_arrays['cloudtype']
cloudtype_qflag = caObj.avhrr.all_arrays['cloudtype_qflag']
if isACPGv2012:
    i_flags = get_day_night_twilight_info_pps2012(cloudtype_qflag)
    (no_qflag, isPPSNight, twilight_flag, isPPSDay, all_dnt_flag) =  i_flags
    lt_flag = get_land_coast_sea_info_pps2012(cloudtype_qflag)
    (no_qflag, land_flag, isPPSSea, coast_flag, all_lsc_flag) = lt_flag
    isPPSIce = get_ice_info_pps2012(cloudtype_qflag)
else:
    i_flags = get_day_night_twilight_info_pps2014(cloudtype_conditions)
    (no_qflag, isPPSNight, twilight_flag, isPPSDay, all_dnt_flag) = i_flags
    lt_flag =get_land_coast_sea_info_pps2014(cloudtype_conditions)
    (no_qflag, land_flag, isPPSSea, coast_flag, all_lsc_flag) = lt_flag 
    isPPSIce = get_ice_info_pps2014(cloudtype_status)

#isClear = np.logical_and(isClear,isPPSCloudyOrClear)
#isCloudy = np.logical_and(isCloudy,isPPSCloudyOrClear)
isThinCirrus = np.logical_and(isCloudy, np.logical_and(
        np.logical_and(SeesThrough,isHigh),
        np.logical_and(isOpticallyThin, isSingleLayerCloud)))

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
isIce =np.logical_and(isPPSSea, isPPSIce)
isWaterNight = np.logical_and(isWater,isPPSNight)
isWaterDay = np.logical_and(isWater,isPPSDay)
isIceNight = np.logical_and(isIce,isPPSNight)
isIceDay = np.logical_and(isIce,isPPSDay)
#isWater =  isPPSSeaNotIce
#isIce =  np.logical_and(isPPSSea,isPPSIce)

isClearIce = np.logical_and(isIce, isClear)
isClearWaterNight = np.logical_and(isClear, isWaterNight)
isCloudyWaterNight = np.logical_and(isCloudy, isWaterNight)
isCloudyCirrusWaterNight = np.logical_and(isThinCirrus, isWaterNight)
isClearIceNight = np.logical_and(isClear, isIceNight)
isCloudyIceNight = np.logical_and(isCloudy, isIceNight)
isCloudyCirrusIceNight = np.logical_and(isThinCirrus, isIceNight)
isClearWaterDay = np.logical_and(isClear, isWaterDay)
isCloudyWaterDay = np.logical_and(isCloudy, isWaterDay)
isCloudyCirrusWaterDay = np.logical_and(isThinCirrus, isWaterDay)
nodata = np.logical_or(t12<= -9, t37<= -9)
t37t12 = np.ma.array(t37-t12, mask = nodata)
nodata = np.logical_or(t11<= -9, t37<= -9)
t11t37 = np.ma.array(t11-t37, mask = nodata)
nodata = np.logical_or(t11<= -9, surftemp<= -9)
t11ts = np.ma.array(t11-surftemp, mask = nodata)
nodata = np.logical_or(t11<= -9, t12<= -9)
t11t12 = np.ma.array(t11-t12, mask = nodata)
nodata = t11text<= -9
t11text = np.ma.array(t11text, mask = nodata)
nodata = t37t12text<= -9
t37t12text = np.ma.array(t37t12text, mask = nodata)
nodata = t37text<= -9
t37text = np.ma.array(t37text, mask = nodata)
nodata = r06<= -9
r06 = np.ma.array(r06, mask = nodata)
nodata = thr_r06<= -9
thr_r06 = np.ma.array(thr_r06, mask = nodata)
t37_t12_minus_threshold = t37t12 - thr_t37t12
t37_t12_minus_threshold_inv = t37t12 - thr_t37t12_inv
t11_t37_minus_threshold = t11t37 - thr_t11t37
t11_ts_minus_threshold = t11ts - thr_t11ts
t11_ts_minus_threshold_inv = t11ts - thr_t11ts_inv
t11_t12_minus_threshold = t11t12 - thr_t11t12
t11_t12_minus_threshold_inv = t11t12 - thr_t11t12_inv
#----------------------------------------------
#Plotting
#----------------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
    isNotDetectedByThisTest = np.logical_and(np.logical_and(isCloudy, 
                                                             TestOk),
                                              isClearPPS)
    #pixels inlcuded in plot and ok test and pps cloudy.
    isPPSCloudyForTestOk = np.logical_and(TestOk,
                                          np.logical_or(isCloudy,isClear))
    POD_cloudy=len(xvector[isDetectedByThisTest==True])*1.0/len(xvector[isCloudy==True])
    FAR_cloudy=np.divide(len(xvector[isMissclassifiedByThisTest==True])*1.0,len(xvector[isPPSCloudyForTestOk==True]))
    POD_cloudy_missed = len(isNotDetectedByThisTest[isNotDetectedByThisTest == True])*1.0/len(xvector[isCloudy==True])
    if len(isNotDetectedByThisTest[isNotDetectedByThisTest == True]) !=0:
        print "warning missed %d"%(len(isNotDetectedByThisTest[isNotDetectedByThisTest == True]))
    POD_FAR_INFO =  ("Number of cloudy pixels considered %d \n"%(len(xvector[isCloudy==True]))+
                     "Number of clouds detected by this test %d \n"%(len(xvector[isDetectedByThisTest==True]))+
                     "Number of clouds misclassified by this test %d \n"%(len(xvector[isMissclassifiedByThisTest==True]))+
                     "Number of cloudy pixels could be detected by this test %d \n"%(len(xvector[isNotDetectedByThisTest==True]))+
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
    plt.plot(xvector[isCloudy == True], yvector[isCloudy == True], 'r.')
    plt.plot(xvector[isCloudyPPSClear == True], yvector[isCloudyPPSClear == True], 'y.')
    plt.plot(xvector[isClear == True], yvector[isClear == True], 'bx')
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
    ax.set_ylim(1.2*limits[2]-0.2*limits[3],limits[3])
    plt.text(0.9*limits[0]+0.1*limits[1],1.15*limits[2]-0.15*limits[3], POD_FAR_INFO, backgroundcolor='w')
    fig_title= (" %s \n"%(args['title'])+
              "blue-cyan: clear (calipso), red-manga: cloudy (clipso)\n"+
              "manga: detected cloudy by this test, cyan: missclassifed by this test)")
    ax.set_title(fig_title)
    #filename=args['title']+"_"+SATELLITE+'.png'
    #plt.savefig('/local_disk/nina_pps/'+ filename)
    #plt.show()

def plot_test(args,isCloudy, isClear, TestOk, THRESHOLD, show=False):
    plot_inner(args,isCloudy, isClear, TestOk, THRESHOLD)
    filename=args['title']+"_"+SATELLITE+'.png'
    plt.savefig('/local_disk/nina_pps/threshold_plots/' + filename)
    if show:
        plt.show()
    plt.close()

def plot_test_2_lim(args,isCloudy, isClear, TestOk, THRESHOLD1, THRESHOLD2, show=False):
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
    plt.savefig('/local_disk/nina_pps/threshold_plots/' + filename)
 
    if show:
        plt.show()
    plt.close()





def print_stats(cloudtype,isCloudyStats,isClearStats,isCloudyPPS,TestOkAll=None):
    print "*********************************"
    print "*** stats so far ****"
    print "Number of cloudy pixels %d"%(
        len(cloudtype[isCloudyStats==True]))
    print "Number of clear pixels  %d"%(
        len(cloudtype[isClearStats==True]))
    print "Part of clouds sea night detected %f"%(sum(np.logical_and(isCloudyStats,isCloudyPPS))*1.0/len(cloudtype[isCloudyStats==True]))
    print "Part of clear sea night missclassified %f"%(sum(np.logical_and(isClearStats,isCloudyPPS))*1.0/len(cloudtype[isCloudyStats==True]))
    if TestOkAll is not None:
        detected_so_far=np.logical_and(TestOkAll, np.logical_and(
                isCloudyStats,isCloudyPPS))
        would_have_been_detected_so_far = np.logical_and(TestOkAll, np.logical_and(
                isCloudyStats,isPPSCloudyOrClear))
        would_have_been_missclassed_so_far = np.logical_and(TestOkAll, np.logical_and(
                isClearStats,isPPSCloudyOrClear))
        false_so_far=np.logical_and(TestOkAll, np.logical_and(isClearStats,isCloudyPPS))
        print "POD cloudy  detected so far %f"%(
            len(cloudtype[detected_so_far==True])*1.0/
            len(cloudtype[isCloudyStats==True]))
        print "FAR cloudy so far %f"%(
            len(cloudtype[false_so_far==True])*1.0/(len(cloudtype[false_so_far==True])+len(cloudtype[detected_so_far==True])))
        print "POD cloudy  detected so far including extra test and new thresholds%f"%(
            len(cloudtype[would_have_been_detected_so_far==True])*1.0/
        len(cloudtype[isCloudyStats==True]))
        print "FAR cloudy so far including extra test and new thresholds %f"%(
            len(cloudtype[would_have_been_missclassed_so_far==True])*1.0/(len(cloudtype[would_have_been_detected_so_far==True])+len(cloudtype[would_have_been_missclassed_so_far==True])))
    print "*********************************"

####################################
# SEA NIGHT
####################################
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS)
#----------------------------------
#ColdWaterCloudTest sea night 
#----------------------------------
print "coldwatercloudtest sea night"
args = {'title': "coldwaterCloudTest_All_SeaNightNoIce",
      'xlable': 'T11-T37 minus dynamic threshold',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0.0 +OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.5
THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_NIGHT']+ -1*OFFSETS['QUALITY_MARGIN_T11TSUR']#-7.0-1.0
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T11_OFFSET_SEA_NIGHT']
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD2,t11_ts_minus_threshold<THRESHOLD1)
TestOkAll=TestOk
plot_test_2_lim(args,isCloudyWaterNight,isClearWaterNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)


#----------------------------------
#ThinCirrusPrimaryTest sea night 
#----------------------------------
print "thinCirrusPrimarytest all"
args = {'title': "thinCirrusPrimaryTest_All_SeaNightNoIce",
        'xlable': 'Tsur',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': surftemp,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T37T12_OFFSET_SEA_NIGHT']+OFFSETS['QUALITY_MARGIN_T37T12']#2.0+0.3
TestOk = t37_t12_minus_threshold>THRESHOLD
#TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
print "thinCirrusPrimarytest cirrus night sea"
args['title']= "thinCirrusPrimaryTest_ThinCirrus_SeaNightNoIce"
plot_test(args,isCloudyCirrusWaterNight,isClearWaterNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
print ciwv
args = {'title': "thinCirrusPrimaryTest_All_SeaNightNoIce_ciwv",
        'xlable': 'Ciwv',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': ciwv,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
plot_test(args,isCloudyCirrusWaterNight,isClearWaterNight,TestOk,THRESHOLD)
args = {'title': "thinCirrusPrimaryTest_All_SeaNightNoIce_lat",
        'xlable': 'Latitude',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': latitude,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
plot_test(args,isCloudyCirrusWaterNight,isClearWaterNight,TestOk,THRESHOLD)
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
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyCirrusWaterNight,isClearWaterNight,TestOk,THRESHOLD)

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
    OFFSETS['T37T12_OFFSET_SEA_NIGHT']
TestOk = args['yvector']>THRESHOLD 
#TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyCirrusWaterNight, isClearWaterNight, TestOk, THRESHOLD)

#----------------------------------
#textureNightTest sea night 
#----------------------------------
print "textureNightTest sea night"
args = {'title': "texturNightTest_All_SeaNightNoIce",
        'xlable': 'T37-T12text minus dynamic threshold',
        'ylable': 'T11text minus dynamic threshold',
        'xvector': t11text,
        'yvector': t37t12text,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'] + OFFSETS['QUALITY_MARGIN_T11TEXT']#0.85+0.15
THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT'] + OFFSETS['QUALITY_MARGIN_T37T12TEXT']#0.75+0.15
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T11TEXT_OFFSET_SEA_NIGHT']
    THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT']
TestOk = np.logical_and(t11text>THRESHOLD2, t37t12text>THRESHOLD1)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterNight, isClearWaterNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)

print "textureNightTest sea night"
args = {'title': "testtexturNightTest_All_SeaNightNoIce",
        'xlable': 'T37-T12text minus dynamic threshold',
        'ylable': 'T11text minus dynamic threshold',
        'xvector': t11text,
        'yvector': t37t12text,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = 0#OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'] + OFFSETS['QUALITY_MARGIN_T11TEXT']#0.85+0.15
THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT'] + OFFSETS['QUALITY_MARGIN_T37T12TEXT']#0.75+0.15
if RunWithOutMargins:
    THRESHOLD2 = 0.0
    THRESHOLD1 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT']

TestOk = np.logical_and(t11text>THRESHOLD2, t37t12text>THRESHOLD1)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterNight, isClearWaterNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)

#----------------------------------
#notexistingCirrustest sea night 
#----------------------------------
print "notexistingCirrustest sea night"
args = {'title': "non_existingT11T12test_All_SeaNightNoIce",
      'xlable': 'Tsur',
      'ylable': 'T11-T12 minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_t12_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.3
TestOk = t11_t12_minus_threshold>THRESHOLD
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)

#----------------------------------
#coldcloudtest sea night 
#----------------------------------
print "coldcloudtest sea night"
args = {'title': "coldCloudTest_All_SeaNightNoIce_with_Tsur_limit",
      'xlable': 'Tsur',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT'] #+ OFFSETS[] # -7.0-1.0
TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] + OFFSETS['QUALITY_MARGIN_TSUR']#
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT'] 
    TSUR_THRESHOLD = OFFSETS['COLDEST_SEASURFACE_TEMP'] 
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
args = {'title': "coldCloudTest_All_SeaNightNoIce_with_Tsur_limit_lat",
      'xlable': 'Latitude',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
args = {'title': "coldCloudTest_All_SeaNightNoIce_with_Tsur_limit_ciwv",
      'xlable': 'Ciwv',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': ciwv,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)

#----------------------------------
#watercloudtest sea night 
#----------------------------------
print "watercloudtest sea night "
args = {'title': "waterCloudTest_All_SeaNightNoIce",
      'xlable': 'Tsur',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37']#0.0+0.3
if RunWithOutMargins:
    THRESHOLD = 0.0

TestOk = t11_t37_minus_threshold>THRESHOLD
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)
args = {'title': "waterCloudTest_All_SeaNightNoIce_lat",
      'xlable': 'Latitude',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': latitude,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
args = {'title': "waterCloudTest_All_SeaNightNoIce_ciwv",
      'xlable': 'Ciwv',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': ciwv,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)

#----------------------------------
#coldcloudtest sea night without tsur limit
#----------------------------------
print "coldcloudtest sea night"
args = {'title': "coldCloudTest_All_SeaNightNoIce_without_Tsur_limit",
      'xlable': 'Tsur',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE'] + -1*OFFSETS['QUALITY_MARGIN_T11TSUR'] #-15.0-1.0
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD = OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE']
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyWaterNight,isClearWaterNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterNight,isClearWaterNight,isCloudyPPS,TestOkAll)


print "WatercloudTestt85t11sea"
print      "**not done yet**"
print "watercloudOverWaterTest"
print "*** need warmest neighbour to plot this test***"
####################################
# NIGHT ICE
####################################
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS)
TestOkAll=None
#----------------------------------
#arcticwaterCloudTest sea night
#----------------------------------
print "arcticwatercloudTest sea night ice"
args = {'title': "arcticWaterCloudTest_All_SeaNightIce",
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
#TestOkAll=TestOk
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

print "testarcticwatercloudTest sea night ice"
args = {'title': "testarcticWaterCloudTest_All_SeaNightIce",
      'xlable': 't3712text',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': t37t12text,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC']+0.3
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T37T12TEXT']+0.3
if RunWithOutMargins:
    THRESHOLD1 = 0.3
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']#MARGIN are offset 
TestOk = np.logical_and(t11_t37_minus_threshold>THRESHOLD1,t37t12text<THRESHOLD2)
TestOkAll=TestOk
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

#----------------------------------
#coldcloudtest sea night ice
#----------------------------------
print "coldcloudtest sea night"
args = {'title': "coldCloudTest_All_SeaNightIce_without_Tsur_limit",
      'xlable': 'Tsur',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 -1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC'] #-15.0-1.0
if RunWithOutMargins:
    THRESHOLD = -1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC']# MARGIN are offset 
TSUR_THRESHOLD = 0.0
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
#TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyIceNight,isClearIceNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

print "coldcloudtest sea night"
args = {'title': "testcoldCloudTest_All_SeaNightIce_without_Tsur_limit",
      'xlable': 'Tsur',
      'ylable': 'T11-Ts minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 -1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC'] +2.0 
TSUR_THRESHOLD = 0.0
if RunWithOutMargins:
    THRESHOLD = -1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC'] + 4.0# MARGIN are offset 
TestOk = np.logical_and(t11_ts_minus_threshold<THRESHOLD,surftemp>TSUR_THRESHOLD)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyIceNight,isClearIceNight,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)


#----------------------------------
#arcticThinCirrusPrimaryTest sea night ice
#----------------------------------
print "arcticthinCirrusPrimarytest all night sea ice "
args = {'title': "arcticThinCirrusPrimaryTest_All_SeaNightIce",
        'xlable': 'T37text',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = 0.0#OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC'] #MARGIN are offset??
TestOk = np.logical_and(t37_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
#TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight,isClearIceNight,TestOk,THRESHOLD1, THRESHOLD2)
print "arcticthinCirrusPrimarytest cirrus night sea ice"
args['title']= "arcticthinCirrusPrimaryTest_ThinCirrus_SeaNightIce"
plot_test_2_lim(args,isCloudyCirrusIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

print "testarcticthinCirrusPrimarytest all night sea ice "
args = {'title': "testarcticThinCirrusPrimaryTest_All_SeaNightIce",
        'xlable': 'T37text',
        'ylable': 'T37-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t37_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 0.3
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']
if RunWithOutMargins:
    THRESHOLD1 = 0.3
    THRESHOLD2 = 16#OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']

TestOk = np.logical_and(t37_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight,isClearIceNight,TestOk,THRESHOLD1, THRESHOLD2)
print "arcticthinCirrusPrimarytest cirrus night sea ice"
args['title']= "testarcticthinCirrusPrimaryTest_ThinCirrus_SeaNightIce"
plot_test_2_lim(args,isCloudyCirrusIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

#----------------------------------
#arcticThinWaterCloudTest sea night ice
#----------------------------------
print "arcticThinWaterCloudTest all"
args = {'title': "arcticThinWaterCloudTest_All_SeaNightIce",
        'xlable': 'T37t12text',
        'ylable': 'T37-T12 minus dynamic threshold inv',
        'xvector': t37t12text,
        'yvector': t37_t12_minus_threshold_inv,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37T12TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC_INV']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT_ARCTIC'] 
    THRESHOLD1 = OFFSETS['QUALITY_MARGIN_T37T12_ARCTIC_INV'] #MARGIN are offset

TestOk = np.logical_and(t37_t12_minus_threshold_inv<THRESHOLD1, t37t12text<THRESHOLD2)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight,isClearIceNight,TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

#----------------------------------
#arcticWarmcloudtest sea night ice
#----------------------------------
print "arcticWarmcloudtest sea night"
args = {'title': "arcticWarmCloudTest_All_SeaNightIce",
      'xlable': 't37text',
      'ylable': 'T11-Ts minus dynamic threshold inv',
      'xvector': t37text,
      'yvector': t11_ts_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 0.0 -1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC_INV']
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']  + OFFSETS['QUALITY_MARGIN_T37TEXT']
if RunWithOutMargins:
    THRESHOLD1 = 0.0 -1*OFFSETS['QUALITY_MARGIN_T11TSUR_ARCTIC_INV'] #MARGIN are offset
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC'] 
TSUR_THRESHOLD = 0.0
TestOk = np.logical_and(t11_ts_minus_threshold_inv>THRESHOLD1,t37text<THRESHOLD2)
#TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

print "testSalomonsarcticWarmcloudtest sea night"
args = {'title': "testSalomonsArcticWarmCloudTest_All_SeaNightIce",
      'ylable': 'T11-Ts minus dynamic threshold',
      'xlable': 'T11-T37 minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 3.0
THRESHOLD2 = -0.3
THRESHOLD3 = 0.3
THRESHOLD4 = 0.6
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold>THRESHOLD1,
                                       t37_t12_minus_threshold<THRESHOLD2),
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD3,
                                       t37t12text<THRESHOLD4))
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

print "testSalomonsarcticWarmcloudtest sea night"
args = {'title': "testSalomonsArcticWarmCloudTest_All_SeaNightIce_changed_off",
      'ylable': 'T11-Ts minus dynamic threshold',
      'xlable': 'T11-T37 minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = 1.5
THRESHOLD2 = -0.4
THRESHOLD3 = -0.3
THRESHOLD4 = 0.6
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold>THRESHOLD1,
                                       t37_t12_minus_threshold<THRESHOLD2),
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD3,
                                       t37t12text<THRESHOLD4))
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

print "testSalomonsarcticWarmcloudtest sea night"
args = {'title': "testSalomonsArcticWarmCloudTest_All_SeaNightIce_changed_off",
      'ylable': 'T11-Ts minus dynamic threshold',
      'xlable': 'T37-T12 minus dynamic threshold',
      'xvector': t11_t37_minus_threshold,
      'yvector': t11_ts_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD1 = -1.0
THRESHOLD2 = -0.4
THRESHOLD3 = -0.3
THRESHOLD4 = 0.6
TestOk = np.logical_and(np.logical_and(t11_ts_minus_threshold>THRESHOLD1,
                                       t37_t12_minus_threshold<THRESHOLD2),
                        np.logical_and(t11_t37_minus_threshold>THRESHOLD3,
                                       t37t12text<THRESHOLD4))
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)




#----------------------------------
#arcticWarmCirrusSecondaryTest sea night
#----------------------------------
print "arcticWarmcloudtest sea night"
args = {'title': "arcticWarmCirrusSecondaryTest_All_SeaNightIce_tsur",
      'xlable': 'Tsur',
      'ylable': 'T11-T12 minus dynamic threshold inv',
      'xvector': surftemp,
      'yvector': t11_t12_minus_threshold_inv,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}

THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC_INV']
if RunWithOutMargins:
    THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC_INV'] #MARGIN are offset
TestOk = t11_t12_minus_threshold_inv<THRESHOLD
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

#----------------------------------
#arcticThinCirrusSecondaryTest sea night ice
#----------------------------------
print "thinCirrusSecondarytest all"
args = {'title': "arcticthinCirrusSecondaryTest_All_SeaNightIce",
        'xlable': 'T37text',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']#MARGIN are offset
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
#TestOkAll=np.logical_or(TestOk,TestOkAll)
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print "thinCirrusSecondarytest cirrus night sea ice"
args['title']= "arcticthinCirrusSecondaryTest_ThinCirrus_SeaNightIce"
plot_test_2_lim(args,isCloudyCirrusIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

print "testthinCirrusSecondarytest all"
args = {'title': "testarcticthinCirrusSecondaryTest_All_SeaNightIce",
        'xlable': 'T37text',
        'ylable': 'T11-T12 minus dynamic threshold',
        'xvector': t37text,
        'yvector': t11_t12_minus_threshold,
        'isCloudyPPS': isCloudyPPS,
        'isClearPPS': isClearPPS,   
}
THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']+OFFSETS['QUALITY_MARGIN_T37TEXT']#2.0+0.3
THRESHOLD1 = OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']
if RunWithOutMargins:
    THRESHOLD2 = OFFSETS['T37TEXT_OFFSET_SEA_NIGHT_ARCTIC']
    THRESHOLD1 = 0.4# OFFSETS['QUALITY_MARGIN_T11T12_ARCTIC']-0.3#MARGIN are offset
TestOk = np.logical_and(t11_t12_minus_threshold>THRESHOLD1, t37text<THRESHOLD2)
plot_test_2_lim(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print "thinCirrusSecondarytest cirrus night sea ice"
args['title']= "testarcticthinCirrusSecondaryTest_ThinCirrus_SeaNightIce"
plot_test_2_lim(args,isCloudyCirrusIceNight, isClearIceNight, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)

#----------------------------------
#extra arcticwaterCloudTest sea night
#----------------------------------
print "watercloudTest sea night ice"
args = {'title': "waterCloudTest_All_SeaNightIce_tsur",
      'xlable': 'Tsur',
      'ylable': 'T11-T37 minus dynamic threshold',
      'xvector': surftemp,
      'yvector': t11_t37_minus_threshold,
      'isCloudyPPS': isCloudyPPS,
      'isClearPPS': isClearPPS,   
}
THRESHOLD = 0.0 + OFFSETS['QUALITY_MARGIN_T11T37_ARCTIC_EXTRA']
if RunWithOutMargins:
    THRESHOLD = 0.0
TestOk = t11_t37_minus_threshold>THRESHOLD
#TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyIceNight, isClearIceNight, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyIceNight,isClearIceNight,isCloudyPPS,TestOkAll)


####################################
# SEA DAY
####################################

print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS)
print "snowSeaTest sea day"
print "not done yet"

#----------------------------------
#coldCloudTest sea day
#----------------------------------
print "coldcloudtest sea day"
args = {'title': "coldCloudTest_All_SeaDayNoIce",
        'xlable': 'Tsur',
        'ylable': 'T11-Ts minus dynamic threshold',
        'xvector': surftemp,
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
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,THRESHOLD)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)
#----------------------------------
#coldBrightCloudTest sea day
#----------------------------------
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)
print "coldBrightCloudTest sea day no ice"
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
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterDay, isClearWaterDay, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)

print "testcoldBrightCloudTest sea day no ice"
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
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test_2_lim(args,isCloudyWaterDay, isClearWaterDay, TestOk, THRESHOLD1, THRESHOLD2)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)

#----------------------------------
#thinCirrusSecondaryTest sea day
#----------------------------------
print "thinCirrusSecondaryTest sea day no ice"
args = {'title': "thinCirrusSecondaryTest_All_SeaDayNoIce_lat",
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
plot_test(args,isCloudyWaterDay, isClearWaterDay, TestOk, THRESHOLD)
print_stats(cloudtype,isCloudyWaterDay,isClearWaterDay,isCloudyPPS,TestOkAll)
#also make plot with only cirrus
args['title'] = "thinCirrusSecondaryTest_ThinCirrus_SeaDayNoIce_lat"
plot_test(args,isCloudyCirrusWaterDay, isClearWaterDay, TestOk, THRESHOLD)
args['xlable'] ='Tsur'
args['xvector'] = surftemp
args['title'] = "thinCirrusSecondaryTest_ThinCirrus_SeaDayNoIce_tsur"
plot_test(args,isCloudyCirrusWaterDay, isClearWaterDay, TestOk, THRESHOLD)
#linear tsur dep threshold
args = {'title': "testThinCirrusPrimaryTest_All_SeaDayNoIce_lat_linearthreshold",
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
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyCirrusWaterDay,isClearWaterDay,TestOk,THRESHOLD)
args['title'] = "testThinCirrusSecondaryTest_All_SeaDayNoIce_lat_linearthreshold_1p5"
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,THRESHOLD)
args = {'title': "testThinCirrusPrimaryTest_All_SeaDayNoIce_lat_linearthreshold_1.0",
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
#TestOkAll=np.logical_or(TestOk,TestOkAll)
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
TestOkAll=np.logical_or(TestOk,TestOkAll)
plot_test(args,isCloudyWaterDay,isClearWaterDay,TestOk,THRESHOLD, show=True)
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

ax.legend()
#ax.set_xlim3d(0, 1)
#ax.set_ylim3d(0, 1)
#ax.set_zlim3d(-3, 10)
#ax.set_zlable('T37-T12')
#ax.set_ylable('Ciwv')
#ax.set_xlable('Tsur')
ax.set_title('T37-T12 minus t37-T12_threshold: \nblue clear, \nred-manga thin cirrus, \nmanga(detected by pps)')
#plt.show()
