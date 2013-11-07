import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matchobject_io import (DataObject)

from cloudsat_calipso_avhrr_statistics import (
    get_day_night_twilight_info_pps2014,
    get_land_coast_sea_info_pps2014,
    get_ice_info_pps2014,
    get_day_night_twilight_info_pps2012,    
    get_land_coast_sea_info_pps2012,
    get_ice_info_pps2012,
    get_inversion_info_pps2014,
    get_sunglint_info_pps2014,
    get_mountin_info_pps2014,
    get_high_terrain_info_pps2014)

def keep_combined_ok(TestOk,TestOkAll):
    if hasattr(TestOk,'mask') and hasattr(TestOkAll,'mask'):               
        return np.logical_or(TestOk.data, TestOkAll.data)
    elif hasattr(TestOkAll,'mask'):               
        return np.logical_or(TestOk, TestOkAll.data)
    elif hasattr(TestOk,'mask'):
        return np.logical_or(TestOk.data, TestOkAll)
    else:
        return np.logical_or(TestOk, TestOkAll)


def plot_inner_clear(SchemeName, args_test, cloudobj, TestOk, THRESHOLD):

    isClearPPS = cloudobj.isPpsClear
    isCloudyPPS = cloudobj.isPpsCloudy
    isClear = cloudobj.isClear.all_arrays[SchemeName]
    isCloudy = cloudobj.isCloudy.all_arrays[SchemeName]
    xvector = args_test['xvector']
    yvector = args_test['yvector']
    fig = plt.figure(figsize = (9,8))
    ax = fig.add_subplot(111)
    isCloudyPPSClear = np.logical_and(isCloudy,isClearPPS)
    isDetectedByThisTest = np.logical_and(np.logical_and(isClear,isClearPPS),
                                          TestOk)
    isMissclassifiedByThisTest = np.logical_and(isCloudy,TestOk)
    isMissclassifiedByThisTestWhereOk = np.logical_and(isCloudyPPS,isMissclassifiedByThisTest)
    isNotDetectedByThisTest = np.logical_and(np.logical_and(isClear, 
                                                             TestOk),
                                             isCloudyPPS)
    #pixels inlcuded in plot and ok test and pps cloudy.
    isPPSClearForTestOk = np.logical_and(TestOk,
                                          np.logical_or(isCloudy,isClear))
    POD_clear=np.divide(len(xvector[isDetectedByThisTest==True])*1.0,
                         len(xvector[isCloudy==True]))
    FAR_clear=np.divide(len(xvector[isMissclassifiedByThisTest==True])*1.0,
                         len(xvector[isPPSClearForTestOk==True]))
    POD_clear_missed = np.divide(len(xvector[isNotDetectedByThisTest==True]),
                                  1.0*len(xvector[isClear==True]))
    POD_FAR_INFO =  (
        "Number of cloudy pixels considered %d \n"%(
            len(xvector[isCloudy==True]))+
        "Number of clear pixels considered %d \n"%(
            len(xvector[isClear==True]))+
        "Number of clear detected by this test %d \n"%(
            len(xvector[isDetectedByThisTest==True]))+
        "Number of clear misclassified by this test %d \n"%(
            len(xvector[isMissclassifiedByThisTest==True]))+
        "Number of clear misclassified by this and cloudy in pps %d \n"%(
            len(xvector[isMissclassifiedByThisTestWhereOk==True]))+
        "Number of clear pixels could be detected by this test %d \n"%(
            len(xvector[isNotDetectedByThisTest==True]))+
        "POD clear by this test %2.1f\n"%(POD_clear*100)+
        "FAR clear %2.1f\n"%(FAR_clear*100)+
        "POD clear missed %2.1f"%(POD_clear_missed*100))
    print "-Test stat---------"
    print args_test['title']
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
             yvector[isDetectedByThisTest == True], 'cx')
    plt.plot(xvector[isMissclassifiedByThisTest == True], 
             yvector[isMissclassifiedByThisTest == True], 'm.')
    plt.plot(xvector[isNotDetectedByThisTest == True], 
             yvector[isNotDetectedByThisTest == True], 'kx')
    plt.axhline(y = THRESHOLD, color = 'k',linewidth = 2)
    plt.axhline(y = clfree_percentiles[0], color = 'b',ls = '--')
    plt.axhline(y = clfree_percentiles[1], color = 'b',ls = '--')
    plt.axhline(y = clfree_percentiles[3], color = 'b',linewidth = 2)
    plt.axhline(y = clfree_percentiles[5], color = 'b',ls = '--')
    plt.axhline(y = clfree_percentiles[6], color = 'b',ls = '--')
    ax.set_ylabel(args_test['ylable'])
    ax.set_xlabel(args_test['xlable'])
    limits = ax.axis()
    ax.set_ylim(1.3*limits[2]-0.3*limits[3],limits[3])
    plt.text(0.9*limits[0]+0.1*limits[1],1.25*limits[2]-0.25*limits[3], 
             POD_FAR_INFO, backgroundcolor='w')
    fig_title= (" %s \n"%(args_test['title'])+
                "blue-cyan: clear (calipso), red-manga: cloudy (clipso)\n"+
                "manga: detected cloudy by this test, " +
                "cyan: missclassifed by this test)")
    ax.set_title(fig_title)


def plot_inner(SchemeName, args_test, cloudobj, TestOk, THRESHOLD, onlyCirrus=None):

    isClearPPS = cloudobj.isPpsClear
    isCloudyPPS = cloudobj.isPpsCloudy
    isClear = cloudobj.isClear.all_arrays[SchemeName]
    isCloudy = cloudobj.isCloudy.all_arrays[SchemeName]
    if onlyCirrus:
        isCloudy= np.logical_and(isCloudy,cloudobj.isCirrus)
    xvector = args_test['xvector']
    yvector = args_test['yvector']
    fig = plt.figure(figsize = (9,8))
    ax = fig.add_subplot(111)
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
        "POD cloudy by this test %2.1f + %2.1f\n"%(POD_cloudy*100, POD_cloudy_missed*100)+
        "FAR cloudy %2.1f\n"%(FAR_cloudy*100))
        #"POD cloudy missed %2.1f"%(POD_cloudy_missed*100))
    print "-Test stat---------"
    print args_test['title']
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
    ax.set_ylabel(args_test['ylable'])
    ax.set_xlabel(args_test['xlable'])
    limits = ax.axis()
    ax.set_ylim(1.3*limits[2]-0.3*limits[3],limits[3])
    plt.text(0.9*limits[0]+0.1*limits[1],1.25*limits[2]-0.25*limits[3], 
             POD_FAR_INFO, backgroundcolor='w')
    fig_title= (" %s \n"%(args_test['title'])+
                "blue-cyan: clear (calipso), red-manga: cloudy (clipso)\n"+
                "manga: detected cloudy by this test, " +
                "cyan: missclassifed by this test)")
    ax.set_title(fig_title)


def plot_test(SchemeName, args_test, args, cloudobj, TestOk, THRESHOLD1, onlyCirrus=None, show=False):
    plot_inner(SchemeName, args_test, cloudobj, TestOk, THRESHOLD1, onlyCirrus=onlyCirrus)
    filename = args_test['title']+"_"+ args['SATELLITE']+'.png'
    plt.savefig(args['PLOT_DIR'] + filename)
    if show:
        plt.show()
    plt.close()

def plot_test_clear(SchemeName, args_test, args, cloudobj, TestOk, THRESHOLD1,  show=False):
    plot_inner_clear(SchemeName, args_test, cloudobj, TestOk, THRESHOLD1)
    filename = args_test['title']+"_"+ args['SATELLITE']+'.png'
    plt.savefig(args['PLOT_DIR'] + filename)
    if show:
        plt.show()
    plt.close()


def plot_test_2_lim(SchemeName, args_test, args, cloudobj, TestOk, 
                    THRESHOLD1, THRESHOLD2, onlyCirrus=None, show=False):
    plot_inner(SchemeName, args_test, cloudobj, TestOk, THRESHOLD1, onlyCirrus=onlyCirrus)
    #Add some vertical lines with threshold and clear percentiles.
    xvector = args_test['xvector']
    clfree_percentiles_vert = np.percentile(
        np.array(xvector[np.logical_and(cloudobj.isClear.all_arrays[SchemeName], 
                                         xvector.mask==False)==True]),
        [01, 10, 25, 50, 75, 90, 99])
    plt.axvline(x = clfree_percentiles_vert[0], color = 'b',ls = '--')
    plt.axvline(x = clfree_percentiles_vert[1], color = 'b',ls = '--')
    plt.axvline(x = clfree_percentiles_vert[3], color = 'b',linewidth = 2)
    plt.axvline(x = clfree_percentiles_vert[5], color = 'b',ls = '--')
    plt.axvline(x = clfree_percentiles_vert[6], color = 'b',ls = '--')
    plt.axvline(x = THRESHOLD2, color = 'k',linewidth = 2)
    filename=args_test['title']+"_" +args['SATELLITE']+'.png'
    plt.savefig(args['PLOT_DIR'] + filename) 
    if show:
        plt.show()
    plt.close()

def print_stats_snow(Scheme,caobj, cloudObj, args):
    cloudtype = caobj.avhrr.cloudtype     
    isClearPPS = cloudObj.isPpsClear
    isCloudyPPS = cloudObj.isPpsCloudy
    isClearStats = cloudObj.isClear.all_arrays[Scheme]
    isCloudyStats = cloudObj.isCloudy.all_arrays[Scheme]
    isPPSCloudyOrClear = cloudObj.isPpsProcessed

    TestOk = np.logical_or(cloudObj.isPpsCloudtypeSnow,
                           cloudObj.isPpsCloudtypeIce)

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
        "PART of all cloudy mis-classed %2.1f"%(PART_cloudy_missclassed*100))
    print POD_FAR_INFO 
    print "*********************************"
    return TestOk

def print_stats(Scheme,caobj, cloudObj, args=None, TestOkAll=None):
    cloudtype = caobj.avhrr.cloudtype   
    isClearPPS = cloudObj.isPpsClear
    isCloudyPPS = cloudObj.isPpsCloudy
    isClearStats = cloudObj.isClear.all_arrays[Scheme]
    isCloudyStats = cloudObj.isCloudy.all_arrays[Scheme]
    isPPSCloudyOrClear = cloudObj.isPpsProcessed
    
    isSnowIce = np.logical_or(cloudObj.isPpsCloudtypeSnow,
                              cloudObj.isPpsCloudtypeIce)
    isSnowIceAndCalipsoCloudy = np.logical_and(isSnowIce, isCloudyStats)
    print "*********************************"
    print "*** stats so far ****"
    num_cloudy =len(cloudtype[isCloudyStats==True])
    num_clear = len(cloudtype[isClearStats==True])
    num_clear_miscl = len(cloudtype[
            np.logical_and(isClearStats,isCloudyPPS)==True])
    num_cloudy_miscl = len(cloudtype[
            np.logical_and(isCloudyStats,isClearPPS)==True])
    num_cloudy_miscl_snowice = len(
       cloudtype[isSnowIceAndCalipsoCloudy==True])
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
    print "Part cloudy missclassed as snow/ice %f "%(
        np.divide(num_cloudy_miscl_snowice*1.0,num_cloudy))
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


class ppsSchemesObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'All': None,
            'BadShemes': None,
            'GoodShemes': None,
            'Coast': None,
            'CoastDay': None,
            'CoastDayMount': None,
            'CoastNight': None,
            'CoastNightInv': None,
            'CoastTwilight': None,
            'CoastTwilightInv': None,
            'Ice': None,
            'IceCoast': None,
            'IceCoastNight': None,
            'IceDay': None,
            'IceNight': None,
            'IceTwilight': None,
            'Land': None,
            'LandDay': None,
            'LandDayMount': None,
            'LandNight': None,
            'LandNightInv': None,
            'LandNightMount': None,
            'LandTwilight': None,
            'LandTwilightInv': None,
            'LandTwilightMount': None,
            'SunglintDay': None,
            'SunglintTwilight': None,
            'Water': None,
            'WaterDay': None,
            'WaterNight': None,
            'WaterTwilight': None,
            }

class cloudClearInfo:
    def __init__(self):
        self.isCloudy = ppsSchemesObject()
        self.isClear = ppsSchemesObject()
        self.isPpsCloudtypeSnow = None
        self.isPpsCloudtypeIce = None
        self.isPpsProcessed = None
        self.isCirrus = None
        self.isPpsClear =   None
        self.isPpsCloudy =  None

def get_clear_and_cloudy_vectors(caObj, isACPGv2012, isGAC):


    isClear = caObj.calipso.all_arrays['number_of_layers_found'] == 0
    isCloudy = caObj.calipso.all_arrays['number_of_layers_found'] >0
    isSingleLayerCloud = caObj.calipso.all_arrays['number_of_layers_found'] == 1
    SeesThrough = caObj.calipso.all_arrays['lidar_surface_elevation'][0] >-999
    isHigh = caObj.calipso.all_arrays['cloud_top_profile'][0] >5.0
    if isGAC:
        isOpticallyThinTopLay = caObj.calipso.all_arrays['optical_depth'][0]<1.0
        notVeryThinTopLay = caObj.calipso.all_arrays['optical_depth'][0]>0.2
        #isOpticallyThin = caObj.calipso.all_arrays['optical_depth'][0]<1.0
        #notVeryThin = caObj.calipso.all_arrays['optical_depth'][0]>0.2
    else:
        isOpticallyThinTopLay = caObj.calipso.all_arrays['optical_depth_top_layer5km']<5.0
        notVeryThinTopLay = caObj.calipso.all_arrays['optical_depth_top_layer5km']>0.2
#isTooThin = caObj.calipso.all_arrays['total_optical_depth5km']<0.2
    isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                                 caObj.avhrr.all_arrays['cloudtype']<21)
    isClearPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                caObj.avhrr.all_arrays['cloudtype']<5)
    isPPSCloudyOrClear = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                        caObj.avhrr.all_arrays['cloudtype']<21)

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
        if isPPSMountain.all()==False:
            isPPSMountain = get_high_terrain_info_pps2014(cloudtype_conditions)

    isClear = np.logical_and(isClear,isPPSCloudyOrClear)
    isCloudy = np.logical_and(isCloudy,isPPSCloudyOrClear)
    isThinCirrus = np.logical_and(isCloudy, np.logical_and(
            np.logical_and(SeesThrough,isHigh),
            np.logical_and(isOpticallyThinTopLay, isSingleLayerCloud)))
    isThinCirrus = np.logical_and(isThinCirrus, notVeryThinTopLay)
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
    isIceCoast = np.logical_and(isPPSIce,isPPSCoast)
    isIceCoastNight = np.logical_and(isPPSNight,isIceCoast)
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
    
    isClearWater = np.logical_and(isWater, isClear)
    isCloudyWater = np.logical_and(isWater, isCloudy)
    isClearIce = np.logical_and(isIce, isClear)
    isCloudyIce = np.logical_and(isIce, isCloudy)
    isClearLand =   np.logical_and(isClear, isPPSLand)
    isCloudyLand =   np.logical_and(isCloudy, isPPSLand)
    isClearCoast =   np.logical_and(isClear, isPPSCoast)
    isCloudyCoast =   np.logical_and(isCloudy, isPPSCoast)
    isClearIceCoast = np.logical_and(isClear, isIceCoast)
    isCloudyIceCoast = np.logical_and(isCloudy, isIceCoast)
    isClearIceCoastNight = np.logical_and(isClear, isIceCoastNight)
    isCloudyIceCoastNight = np.logical_and(isCloudy, isIceCoastNight)
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
    isPpsCloudtypeIce =np.equal(cloudtype,4)
    isPpsCloudtypeSnow = np.equal(cloudtype,3)

    isBadShemesCloudy=np.logical_or(isCloudyLandNightInv,np.logical_or(isCloudyLandNightMount,isCloudyLandTwilight))
    isBadShemesClear=np.logical_or(isClearLandNightInv,np.logical_or(isClearLandNightMount,isClearLandTwilight))
    isGoodShemesCloudy=np.logical_and(np.logical_and(isPPSCloudyOrClear,isCloudy),
                                          np.equal(isBadShemesCloudy, False))
    isGoodShemesClear=np.logical_and(np.logical_and(isPPSCloudyOrClear,isClear),
                                     np.equal(isBadShemesClear, False))
    
    cloudObj=cloudClearInfo()
    
    cloudObj.isPpsCloudtypeSnow = isPpsCloudtypeSnow 
    cloudObj.isPpsCloudtypeIce =   isPpsCloudtypeIce
    cloudObj.isPpsProcessed = isPPSCloudyOrClear
    cloudObj.isCirrus = isThinCirrus
    cloudObj.isPpsClear =  isClearPPS
    cloudObj.isPpsCloudy = isCloudyPPS

    cloudObj.isClear.BadShemes = isBadShemesClear
    cloudObj.isCloudy.BadShemes = isBadShemesCloudy
    cloudObj.isClear.GoodShemes = isGoodShemesClear
    cloudObj.isCloudy.GoodShemes = isGoodShemesCloudy
                    
    cloudObj.isClear.All = isClear 
    cloudObj.isClear.Coast = isClearCoast
    cloudObj.isClear.CoastDay = isClearCoastDay 
    cloudObj.isClear.CoastDayMount =  isClearCoastDayMount
    cloudObj.isClear.CoastNight = isClearCoastNight 
    cloudObj.isClear.CoastNightInv = isClearCoastNightInv 
    cloudObj.isClear.CoastTwilight =  isClearCoastTwilight
    cloudObj.isClear.CoastTwilightInv =  isClearCoastTwilightInv
    cloudObj.isClear.Ice = isClearIce
    cloudObj.isClear.IceCoast = isClearIceCoast
    cloudObj.isClear.IceCoastNight = isClearIceCoastNight
    cloudObj.isClear.IceDay =  isClearIceDay
    cloudObj.isClear.IceNight = isClearIceNight
    cloudObj.isClear.IceTwilight =  isClearIceTwilight
    cloudObj.isClear.Land =  isClearLand
    cloudObj.isClear.LandDay = isClearLandDay 
    cloudObj.isClear.LandDayMount = isClearLandDayMount 
    cloudObj.isClear.LandNight = isClearLandNight 
    cloudObj.isClear.LandNightInv =  isClearLandNightInv
    cloudObj.isClear.LandNightMount = isClearLandNightMount 
    cloudObj.isClear.LandTwilight =  isClearLandTwilight
    cloudObj.isClear.LandTwilightInv = isClearLandTwilightInv 
    cloudObj.isClear.LandTwilightMount =  isClearLandTwilightMount
    cloudObj.isClear.SunglintDay = isClearSunglintDay 
    cloudObj.isClear.SunglintTwilight = isClearSunglintTwilight 
    cloudObj.isClear.Water = isClearWater 
    cloudObj.isClear.WaterDay = isClearWaterDay 
    cloudObj.isClear.WaterNight =  isClearWaterNight
    cloudObj.isClear.WaterTwilight =  isClearWaterTwilight
    cloudObj.isCloudy.All = isCloudy
    cloudObj.isCloudy.Coast =isCloudyCoast
    cloudObj.isCloudy.CoastDay =isCloudyCoastDay 
    cloudObj.isCloudy.CoastDayMount =isCloudyCoastDayMount
    cloudObj.isCloudy.CoastNight =isCloudyCoastNight 
    cloudObj.isCloudy.CoastNightInv =isCloudyCoastNightInv 
    cloudObj.isCloudy.CoastTwilight = isCloudyCoastTwilight
    cloudObj.isCloudy.CoastTwilightInv = isCloudyCoastTwilightInv
    cloudObj.isCloudy.Ice = isCloudyIce
    cloudObj.isCloudy.IceCoast = isCloudyIceCoast
    cloudObj.isCloudy.IceCoastNight = isCloudyIceCoastNight
    cloudObj.isCloudy.IceDay = isCloudyIceDay
    cloudObj.isCloudy.IceNight = isCloudyIceNight
    cloudObj.isCloudy.IceTwilight = isCloudyIceTwilight
    cloudObj.isCloudy.Land = isCloudyLand
    cloudObj.isCloudy.LandDay =isCloudyLandDay 
    cloudObj.isCloudy.LandDayMount = isCloudyLandDayMount
    cloudObj.isCloudy.LandNight =isCloudyLandNight 
    cloudObj.isCloudy.LandNightInv = isCloudyLandNightInv
    cloudObj.isCloudy.LandNightMount = isCloudyLandNightMount
    cloudObj.isCloudy.LandTwilight = isCloudyLandTwilight
    cloudObj.isCloudy.LandTwilightInv = isCloudyLandTwilightInv
    cloudObj.isCloudy.LandTwilightMount = isCloudyLandTwilightMount
    cloudObj.isCloudy.SunglintDay =isCloudySunglintDay 
    cloudObj.isCloudy.SunglintTwilight =isCloudySunglintTwilight 
    cloudObj.isCloudy.Water =isCloudyWater
    cloudObj.isCloudy.WaterDay =isCloudyWaterDay 
    cloudObj.isCloudy.WaterNight = isCloudyWaterNight
    cloudObj.isCloudy.WaterTwilight = isCloudyWaterTwilight
 
    return  cloudObj



class ppsTestFeature(DataObject):
    def __init__(self):
        DataObject.__init__(self)  
        self.all_arrays = {              
            'pseudor06' : None,
            'r06': None,
            'r09': None,
            'r16': None,
            't85_t11_minus_threshold': None,
            't85_t11_minus_threshold_inv': None,
            't37_t12_minus_threshold': None,
            't37_t12_minus_threshold_inv': None,
            't11_t37_minus_threshold': None,
            't11_t37_minus_threshold_inv': None,
            't11_ts_minus_threshold': None,
            't11_ts_minus_threshold_inv': None,
            't11_t12_minus_threshold': None,
            't11_t12_minus_threshold_inv': None,
            'thr_r06': None,
            'text_r06': None,
            'text_t11': None,
            'text_t37t12': None,
            'text_t37': None,
            'sunz': None,
            'surftemp': None,
            'ciwv': None,
            'latitude': None,
            't11t37': None,
            't37t12': None
          }

def get_feature_values_and_thresholds(caObj):
    feature_minus_thr = ppsTestFeature()

    sunz = caObj.avhrr.all_arrays['sunz']
    surftemp = caObj.avhrr.all_arrays['surftemp']
    feature_minus_thr.sunz = caObj.avhrr.all_arrays['sunz']
    feature_minus_thr.surftemp = caObj.avhrr.all_arrays['surftemp']
    feature_minus_thr.latitude = caObj.avhrr.all_arrays['latitude']
    feature_minus_thr.ciwv = caObj.avhrr.all_arrays['ciwv']

    feature_minus_thr.thr_r06 = caObj.avhrr.all_arrays['thr_r06']

    #thresholds
    thr_t11t12 = caObj.avhrr.thr_t11t12
    thr_t11t12 = caObj.avhrr.thr_t11t12
    thr_t11t37 = caObj.avhrr.thr_t11t37
    thr_t37t12 = caObj.avhrr.thr_t37t12
    thr_t11ts = caObj.avhrr.thr_t11ts
    thr_t85t11 = caObj.avhrr.thr_t85t11
    thr_t11ts_inv = caObj.avhrr.thr_t11ts_inv
    thr_t11t12_inv = caObj.avhrr.thr_t11t12_inv
    thr_t37t12_inv = caObj.avhrr.thr_t37t12_inv
    thr_t11t37_inv = caObj.avhrr.thr_t11t37_inv
    thr_t85t11_inv = caObj.avhrr.thr_t85t11_inv
    #brightness temperature
    t11 = caObj.avhrr.all_arrays['bt11micron']
    t85 = caObj.avhrr.all_arrays['bt86micron']
    t12 = caObj.avhrr.all_arrays['bt12micron']
    t37 = caObj.avhrr.all_arrays['bt37micron']
    #radiances
    r06 = caObj.avhrr.all_arrays['r06micron']
    r09 = caObj.avhrr.all_arrays['r09micron']
    r16 = caObj.avhrr.all_arrays['r16micron']
    #texture
    feature_minus_thr.t11text = caObj.avhrr.all_arrays['text_t11']
    feature_minus_thr.r06text = caObj.avhrr.all_arrays['text_r06']
    feature_minus_thr.t37t12text = caObj.avhrr.all_arrays['text_t37t12']
    feature_minus_thr.t37text = caObj.avhrr.all_arrays['text_t37']

    #feature_diff
    nodata = np.logical_or(t12<= -9, t37<= -9)
    t37t12 = np.ma.array(t37-t12, mask = nodata)
    feature_minus_thr.t37t12 = t37t12
    nodata = np.logical_or(t11<= -9, t37<= -9)
    t11t37 = np.ma.array(t11-t37, mask = nodata)
    feature_minus_thr.t11t37 = t11t37
    try:
        nodata = np.logical_or(t11<= -9, t85<= -9)
        t85t11 = np.ma.array(t85-t11, mask = nodata)
    except:
        t85t11=None
    nodata = np.logical_or(t11<= -9, surftemp<= -9)
    t11ts = np.ma.array(t11-surftemp, mask = nodata)
    nodata = np.logical_or(t11<= -9, t12<= -9)
    t11t12 = np.ma.array(t11-t12, mask = nodata)

    #scale radiances
    feature_minus_thr.pseudor06 = r06
    feature_minus_thr.r06=r06/np.cos(np.radians(sunz))
    feature_minus_thr.r09=r09/np.cos(np.radians(sunz))
    try:
        feature_minus_thr.r16=r16/np.cos(np.radians(sunz))
    except:
        feature_minus_thr.r16=None
    nodata = np.logical_or(r09<=-9,np.logical_or(r06<=-0.9, r06==0))
    try:
        feature_minus_thr.qr09r06 = np.ma.array(r09/r06,  mask = nodata)
    except:
        feature_minus_thr.qr09r06 = None
    nodata = np.logical_or(r16<=-9,np.logical_or(r06<=-0.9, r06==0))
    try:
        feature_minus_thr.qr16r06 = np.ma.array(r16/r06,  mask = nodata)
    except:
        feature_minus_thr.qr16r06 = None


    #test with 1 and 1.5 std thresholds:
    cu=3.5/4.0
    cl=0.5/4.0
    cu=3/4.0
    cl=1/4.0
    t37_t12_minus_threshold = t37t12 - cu*thr_t37t12 - cl*thr_t37t12_inv
    t37_t12_minus_threshold_inv = t37t12 - cu*thr_t37t12_inv - cl*thr_t37t12
    t11_t37_minus_threshold = t11t37 - cu*thr_t11t37 - cl*thr_t11t37_inv
    t11_t37_minus_threshold_inv = t11t37 - cu*thr_t11t37_inv - cl*thr_t11t37
    t11_ts_minus_threshold = t11ts - cu*thr_t11ts - cl*thr_t11ts_inv
    t11_ts_minus_threshold_inv = t11ts - cu*thr_t11ts_inv - cl*thr_t11ts
    t11_t12_minus_threshold = t11t12 - cu*thr_t11t12 - cl*thr_t11t12_inv
    t11_t12_minus_threshold_inv = t11t12 - cu*thr_t11t12_inv - cl*thr_t11t12

    #feature minus thresholds
    try:
        feature_minus_thr.t85_t11_minus_threshold = t85t11 - thr_t85t11 
        feature_minus_thr.t85_t11_minus_threshold_inv = t85t11 - thr_t85t11_inv
    except:
        pass
    feature_minus_thr.t37_t12_minus_threshold = t37t12 - thr_t37t12 
    feature_minus_thr.t37_t12_minus_threshold_inv = t37t12 - thr_t37t12_inv 
    feature_minus_thr.t11_t37_minus_threshold = t11t37 - thr_t11t37 
    feature_minus_thr.t11_t37_minus_threshold_inv = t11t37 - thr_t11t37_inv 
    feature_minus_thr.t11_ts_minus_threshold = t11ts - thr_t11ts 
    feature_minus_thr.t11_ts_minus_threshold_inv = t11ts - thr_t11ts_inv 
    feature_minus_thr.t11_t12_minus_threshold = t11t12 - thr_t11t12 
    feature_minus_thr.t11_t12_minus_threshold_inv = t11t12 - thr_t11t12_inv 

    feature_minus_thr.t85t11 = t85t11
    feature_minus_thr.t37t12 = t37t12
    feature_minus_thr.t11t37 = t11t37
    feature_minus_thr.t11ts = t11ts
    feature_minus_thr.t11t12 = t11t12

    return feature_minus_thr


