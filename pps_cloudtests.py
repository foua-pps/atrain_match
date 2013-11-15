import numpy as np
from pps_threshold_functions import (plot_test,
                                     plot_test_2_lim,
                                     plot_test_clear)

def thr_slope_between_90_and_day(thr_night, thr_day, OFFSETS, thr):
    c = 1.0/(90.0-OFFSETS['MAX_SUNZEN_DAY'])
    intercept =((1-c*90.0)*thr_night +
                c*90.0*thr_day)
    slope = c*(thr_night-thr_day)
    thr_thr =np.where(
        thr.sunz>=90.0,
        thr_night,
        slope * thr.sunz + intercept)
    return thr_thr

def coldCloudTest_v2014(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  info="", onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_coldCloudTest_v2014 = {}
    OFFSETS_coldCloudTest_v2014['t11ts_OFFSETS'] = { 
        'All': -4,
        'CoastDay': -6,  #-8 check with avhrr
        'CoastDayMount': -11, 
        'CoastNight': np.min([OFFSETS['T11_OFFSET_LAND_NIGHT'], 
                              OFFSETS['T11_OFFSET_SEA_NIGHT']]), 
        'CoastNightInv': -8, #-11 check with avhrr
        'CoastTwilight': -5, 
        'CoastTwilightInv': -7, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': -8, 
        'LandDayMount': -9, 
        'LandNight': OFFSETS['T11_OFFSET_LAND_NIGHT'], 
        'LandNightInv': -8,#Same as land night and land night mountain
        'LandNightMount': -8, 
        'LandTwilight': -6.5, 
        'LandTwilightInv': -6,#-8, check with avhrr 
        'LandTwilightMount': -12, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': -8,#-5 
        'WaterNight': -7, 
        'WaterTwilight': -4#-4 
        }
    OFFSETS_coldCloudTest_v2014['tsur_threshold'] = {  
        'All': 0,
        'CoastDay': OFFSETS['COLDEST_SEASURFACE_TEMP'], #add??
        'CoastDayMount': 0, 
        'CoastNight': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'CoastNightInv': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'CoastTwilight': 0, 
        'CoastTwilightInv': 0, 
        'IceDay': 0, 
        'IceNight': 0, 
        'IceTwilight': 0, 
        'LandDay': 0, 
        'LandDayMount': 0, 
        'LandNight': 0, 
        'LandNightInv': 0, 
        'LandNightMount': 0, 
        'LandTwilight': 0, 
        'LandTwilightInv': 0, 
        'LandTwilightMount': 230.0, 
        'SunglintDay': 0, 
        'SunglintTwilight': 0, 
        'WaterDay': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'WaterNight': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'WaterTwilight': 0#OFFSETS['COLDEST_SEASURFACE_TEMP']
        }
    if onlyCirrus:
        title =info+"coldCloudTest_v2014_Cirrus_"+ SchemeName
    else:
        title =info+"coldCloudTest_v2014_All_"+ SchemeName

    args_test = {'title':  title,
                 'xlable': 'T11-T37 minus dynamic threshold_inv',
                 'ylable': 'T11-Ts minus dynamic threshold',
                 'xvector': thr.t11_t37_minus_threshold_inv,
                 'yvector': thr.t11_ts_minus_threshold,
            }
    T11_OFFSET = OFFSETS_coldCloudTest_v2014['t11ts_OFFSETS'][SchemeName] 
    TSUR_THRESHOLD = OFFSETS_coldCloudTest_v2014['tsur_threshold'][SchemeName]

    dt11tsur_offset_extra = -8.0
    t11t37_clfree = 0.5*(caobj.avhrr.thr_t11t37 + caobj.avhrr.thr_t11t37_inv)
    t37t12_clfree = 0.5*(caobj.avhrr.thr_t37t12 + caobj.avhrr.thr_t37t12_inv)
    t11t37_departure = np.abs(thr.t11t37 - t11t37_clfree)
    t37t12_departure = np.abs(thr.t37t12 - t37t12_clfree)
    t11t37_maxdep = np.abs(0.5*(caobj.avhrr.thr_t11t37 - caobj.avhrr.thr_t11t37_inv))
    t37t12_maxdep = np.abs(0.5*(caobj.avhrr.thr_t37t12 - caobj.avhrr.thr_t37t12_inv))
    t11t37_div = t11t37_departure/t11t37_maxdep
    t37t12_div = t37t12_departure/t37t12_maxdep
    offset = (T11_OFFSET + dt11tsur_offset_extra * 
              (1 - np.where(t11t37_div>1.0,1.0, t11t37_div)) * 
              (1 - np.where(t37t12_div>1.0,1.0, t37t12_div)))
    THRESHOLD = offset
    if args['USE_MARGINS']:
        THRESHOLD = offset - OFFSETS['QUALITY_MARGIN_T11TSUR']
    TestOk = np.logical_and(thr.t11_ts_minus_threshold < THRESHOLD,
                            thr.surftemp>TSUR_THRESHOLD)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, np.min(THRESHOLD), 
              onlyCirrus=onlyCirrus, show=show)

    return TestOk


def coldWatercloudTest_v2014(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  
                               info="", onlyCirrus=False, ExtraCond=None, show=False):
    if onlyCirrus:
        title =info+"coldWatercloudTest_v2014_Cirrus_"+ SchemeName
    else:
        title =info+"coldWatercloudTest_v2014_All_"+ SchemeName

    args_test = {'title': "coldwaterCloudTest_v2014"+ SchemeName,
            'xlable': 'T11-T37 minus dynamic threshold',
            'ylable': 'T11-Ts minus dynamic threshold',
            'xvector': thr.t11_t37_minus_threshold,
            'yvector': thr.t11_ts_minus_threshold,
            }

    TSUR_THRESHOLD = 0#OFFSETS_coldWatercloudTest_v2014['tsur_threshold'][SchemName] 
    THRESHOLD2 = -7.00#OFFSETS_coldWatercloudTest_v2014['off_t11t37'][SchemName]
    THRESHOLD1 = 0#OFFSETS_coldWatercloudTest_v2014['off_t11tsur'][SchemeName] 
    if  args['USE_MARGINS']:
        THRESHOLD2 += thresholds_dict['MARGIN_T37T12'] 
        THRESHOLD1 -= thresholds_dict['MARGIN_T11TSUR']
    TestOk = np.logical_and(
        np.logical_and(thr.surftemp>TSUR_THRESHOLD,
                       (((thr.t11_t37_minus_threshold-THRESHOLD2)/11)**2 +
                        ((thr.t11_ts_minus_threshold-THRESHOLD1)/20)**2)>1.0),
        np.logical_and(thr.t11_t37_minus_threshold>THRESHOLD2,
                       thr.t11_ts_minus_threshold<THRESHOLD1))
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def thincoldCirrusTest_v2014(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  
                             info="", onlyCirrus=False, ExtraCond=None, show=False):
    if onlyCirrus:
        title =info+"thincoldCirrusTest_v2014_Cirrus_"+ SchemeName
    else:
        title =info+"thincoldCirrusTest_v2014_All_ "+ SchemeName

    args_test = {'title': title ,
            'xlable': 'T11-TS minus dynamic threshold',
            'ylable': 'T11-T12 minus dynamic threshold',
            'xvector': thr.t11_ts_minus_threshold,
            'yvector': thr.t11_t12_minus_threshold, 
            }
    THRESHOLD1 = 0#OFFSETS['T11T12_OFFSET_LAND_DAY']
    THRESHOLD2 = -2.0#OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE']
    TSUR_THRESHOLD = 0.0
    TestOk = np.logical_and(
        np.logical_and(thr.surftemp>TSUR_THRESHOLD,
                       (((thr.t11_t12_minus_threshold-THRESHOLD2)/4)**2 +
                        ((thr.t11_ts_minus_threshold-THRESHOLD1)/20)**2)>1.0  ),
        np.logical_and(thr.t11_t12_minus_threshold>THRESHOLD1 ,
                       thr.t11_ts_minus_threshold<THRESHOLD2))
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def coldWatercloudTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  
                               info="", onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_coldWatercloudTest = {}
    OFFSETS_coldWatercloudTest['t11ts_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': OFFSETS['T11_OFFSET_SEA_NIGHT'], 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': OFFSETS['T11_OFFSET_LAND_NIGHT'], 
        'LandNightInv': OFFSETS['T11_OFFSET_LAND_NIGHT']+2.0, #try 
        'LandNightMount': OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT'], 
        'LandTwilight': np.min([OFFSETS['T11_OFFSET_LAND_NIGHT'], 
                                OFFSETS['T11_OFFSET_LAND_DAY']]), 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11_OFFSET_SEA_DAY'], 
        'WaterNight': OFFSETS['T11_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': np.min([OFFSETS['T11_OFFSET_SEA_NIGHT'], 
                                 OFFSETS['T11_OFFSET_SEA_DAY']])
        }    
    OFFSETS_coldWatercloudTest['t11t37_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 0.0, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': 0.0, 
        'LandNightInv': 0.0, 
        'LandNightMount': 0.0, 
        'LandTwilight': 0.0, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 0.0, 
        'WaterNight': 0.0, 
        'WaterTwilight': 0.0 
        }
    if onlyCirrus:
        title =info+"coldWatercloudTest_Cirrus_"+ SchemeName
    else:
        title =info+"coldWatercloudTest_All_ "+ SchemeName

    args_test = {'title':title,
                 'xlable': 'T11-T37 minus dynamic threshold',
                 'ylable': 'T11-Ts minus dynamic threshold',
                 'xvector': thr.t11_t37_minus_threshold,
                 'yvector': thr.t11_ts_minus_threshold,  
                 }
    THRESHOLD2 = OFFSETS_coldWatercloudTest['t11t37_OFFSETS'][SchemeName]    
    THRESHOLD1 = OFFSETS_coldWatercloudTest['t11ts_OFFSETS'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_T11T37']
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR']
    TestOk =  np.logical_and(thr.t11_t37_minus_threshold>THRESHOLD2,
                             thr.t11_ts_minus_threshold<THRESHOLD1)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk


def watercloudTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  
                               info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_watercloudTest = {}   
    OFFSETS_watercloudTest['t11t37_OFFSETS'] = { 
        'All': 3.0,
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'], 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'], 
        'LandNightInv':  OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] + 2.0, 
        'LandNightMount': OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'], 
        'LandTwilight': 0.0, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 0+OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']+1.0, 
        'WaterTwilight': OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET']#+1.0 
        }

    if onlyCirrus:
        title =info+"watercloudTest_Cirrus_"+ SchemeName
    else:
        title =info+"watercloudTest_All_ "+ SchemeName

    args_test = {'title':title,
                 'xlable': 'Latitude',
                 'ylable': 'T11-T37 minus dynamic threshold',
                 'xvector': thr.latitude,
                 'yvector': thr.t11_t37_minus_threshold
                 }
    THRESHOLD = OFFSETS_watercloudTest['t11t37_OFFSETS'][SchemeName]
    if NEW_THRESHOLD is not None:
        THRESHOLD = NEW_THRESHOLD 
    if  args['USE_MARGINS']:
        THRESHOLD +=  OFFSETS['QUALITY_MARGIN_T11T37']
    TestOk =  thr.t11_t37_minus_threshold>THRESHOLD
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, 
                    np.mean(THRESHOLD), onlyCirrus=onlyCirrus, show=show)
    return TestOk

def coldCloudTest_tsurf_lim(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                               info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):   
    OFFSETS_coldCloudTest_tsurf_lim = {}
    OFFSETS_coldCloudTest_tsurf_lim['t11ts_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': OFFSETS['T11_OFFSET_INVERSION_WEAK'], 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': 100, 
        'LandNightInv': 100, 
        'LandNightMount': 100, 
        'LandTwilight': np.min([OFFSETS['T11_OFFSET_LAND_DAY'], 
                                OFFSETS['T11_OFFSET_LAND_NIGHT']])-1.0, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11_OFFSET_SEA_DAY'], 
        'WaterNight': OFFSETS['T11_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': np.min([OFFSETS['T11_OFFSET_SEA_DAY'],
                                 OFFSETS['T11_OFFSET_SEA_NIGHT']])
        }     
    OFFSETS_coldCloudTest_tsurf_lim['tsur_threshold'] = {  
        'All': 0,
        'CoastDay': 0, 
        'CoastDayMount': 0, 
        'CoastNight': 0, 
        'CoastNightInv': 0, 
        'CoastTwilight': 0, 
        'CoastTwilightInv': 0, 
        'IceDay': 0, 
        'IceNight': 0, 
        'IceTwilight': 0, 
        'LandDay': 0, 
        'LandDayMount': 0, 
        'LandNight': 0, 
        'LandNightInv': 0, 
        'LandNightMount': 0, 
        'LandTwilight': 0, 
        'LandTwilightInv': 0, 
        'LandTwilightMount': 0, 
        'SunglintDay': 0, 
        'SunglintTwilight': 0, 
        'WaterDay': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'WaterNight': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'WaterTwilight': OFFSETS['COLDEST_SEASURFACE_TEMP']        }
    if onlyCirrus:
        title =info+"coldCloudTest_tsurf_lim_Cirrus_"+ SchemeName
    else:
        title =info+"coldCloudTest_tsurf_lim_All_ "+ SchemeName
    args_test = {'title': title,
            'xlable': 'Latitude',
            'ylable': 'T11-Ts minus dynamic threshold',
            'xvector': thr.latitude,
            'yvector': thr.t11_ts_minus_threshold
            }
    THRESHOLD1 = OFFSETS_coldCloudTest_tsurf_lim['t11ts_OFFSETS'][SchemeName]
    TSUR_THRESHOLD = OFFSETS_coldCloudTest_tsurf_lim['tsur_threshold'][SchemeName]
    if NEW_THRESHOLD is not None:
        THRESHOLD = NEW_THRESHOLD
    if  args['USE_MARGINS']:
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR'] 
        TSUR_THRESHOLD += OFFSETS['QUALITY_MARGIN_TSUR']
    TestOk = np.logical_and(thr.t11_ts_minus_threshold<THRESHOLD1,
                            thr.surftemp>TSUR_THRESHOLD)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD1, 
              onlyCirrus=onlyCirrus, show=show)
    return TestOk


def coldCloudTest_no_tsurf_lim(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                               info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):   
    OFFSETS_coldCloudTest = {}
    OFFSETS_coldCloudTest['t11ts_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': OFFSETS['T11_OFFSET_INVERSION_WEAK'], 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': OFFSETS['T11_OFFSET_LAND_NIGHT_OPAQUE'], 
        'LandNightInv': OFFSETS['T11_OFFSET_INVERSION_WEAK'], 
        'LandNightMount': OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT'], 
        'LandTwilight': 0.5*(OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'] + 
                             OFFSETS['T11_OFFSET_LAND_NIGHT_OPAQUE']), 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'], 
        'WaterNight': OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE'], 
        'WaterTwilight': OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE'] 
        }      
    if onlyCirrus:
        title =info+"coldCloudTest_no_tsurf_lim_Cirrus_"+ SchemeName
    else:
        title =info+"coldCloudTest_no_tsurf_lim_All_ "+ SchemeName
    args_test = {'title': title,
            'xlable': 'Latitude',
            'ylable': 'T11-Ts minus dynamic threshold',
            'xvector': thr.latitude,
            'yvector': thr.t11_ts_minus_threshold,
            }
    THRESHOLD1 = OFFSETS_coldCloudTest['t11ts_OFFSETS'][SchemeName]
    if NEW_THRESHOLD is not None:
        THRESHOLD1 = NEW_THRESHOLD
    if  args['USE_MARGINS']:
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR'] 
    TestOk = thr.t11_ts_minus_threshold<THRESHOLD1
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD1, 
              onlyCirrus=onlyCirrus, show=show)
    return TestOk



def thinCirrusSecondaryTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                               info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):
    tmp_night=0.5*(OFFSETS['T11T12_OFFSET_LAND_NIGHT']+
                      OFFSETS['T11T12_OFFSET_SEA_NIGHT'])
    tmp_day=0.5*( OFFSETS['T11T12_OFFSET_LAND_DAY']+
                     OFFSETS['T11T12_OFFSET_SEA_DAY'])
    OFFSETS_thinCirrusSecondaryTest = {}
    OFFSETS_thinCirrusSecondaryTest['t11t12_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': thr_slope_between_90_and_day(tmp_night, tmp_day, OFFSETS, thr), 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': OFFSETS['T11T12_OFFSET_LAND_NIGHT'], 
        'LandNightInv': OFFSETS['T11T12_OFFSET_LAND_NIGHT'], 
        'LandNightMount':  OFFSETS['T11T12_OFFSET_LAND_NIGHT'], 
        'LandTwilight': thr_slope_between_90_and_day(OFFSETS['T11T12_OFFSET_LAND_NIGHT'],
                                                     OFFSETS['T11T12_OFFSET_LAND_DAY'], 
                                                     OFFSETS, thr), 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11T12_OFFSET_SEA_DAY'] , 
        'WaterNight': 100, 
        'WaterTwilight': thr_slope_between_90_and_day(OFFSETS['T11T12_OFFSET_SEA_NIGHT'], 
                                                      OFFSETS['T11T12_OFFSET_SEA_DAY'],
                                                      OFFSETS, thr),
        #'WaterTwilight': 0.5*(OFFSETS['T11T12_OFFSET_SEA_NIGHT'] + 
        #                      OFFSETS['T11T12_OFFSET_SEA_DAY'])
        } 
    if onlyCirrus:
        title =info+"thinCirrusSecondaryTest_Cirrus_"+ SchemeName
    else:
        title =info+"thinCirrusSecondaryTest_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'Latitude',
                 'ylable': 'T11-T12 minus dynamic threshold',
                 'xvector': thr.latitude,
                 'yvector': thr.t11_t12_minus_threshold,
                 }
    THRESHOLD = OFFSETS_thinCirrusSecondaryTest['t11t12_OFFSETS'][SchemeName]
    if NEW_THRESHOLD is not None:
        THRESHOLD = NEW_THRESHOLD 
    if  args['USE_MARGINS']:
        THRESHOLD += OFFSETS['T11T12_OFFSET_LAND_NIGHT'] 
    TestOk = thr.t11_t12_minus_threshold>THRESHOLD 
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, np.mean(THRESHOLD), 
              onlyCirrus=onlyCirrus, show=show)
    return TestOk



def thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                          info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):
    tmp_night=np.max([OFFSETS['T37T12_OFFSET_LAND_NIGHT'],
                      OFFSETS['T37T12_OFFSET_SEA_NIGHT']])
    tmp_day=np.max([ OFFSETS['T37T12_OFFSET_LAND_DAY'], 
                     OFFSETS['T37T12_OFFSET_SEA_DAY']])
    OFFSETS_thinCirrusPrimaryTest = {}
    OFFSETS_thinCirrusPrimaryTest['t37t12_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': np.max([OFFSETS['T37T12_OFFSET_SEA_NIGHT'],
                                 OFFSETS['T37T12_OFFSET_LAND_NIGHT']])-1.0, 
        'CoastTwilight': thr_slope_between_90_and_day(tmp_night, tmp_day, OFFSETS, thr), 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -2.0,
        'LandNightInv': OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -1.0, 
        'LandNightMount': OFFSETS['T37T12_OFFSET_LAND_NIGHT']-1.0, 
        'LandTwilight': thr_slope_between_90_and_day(OFFSETS['T37T12_OFFSET_LAND_NIGHT'],
                                                     OFFSETS['T37T12_OFFSET_LAND_DAY'], 
                                                     OFFSETS, thr), 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['T37T12_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': thr_slope_between_90_and_day(OFFSETS['T37T12_OFFSET_SEA_NIGHT'],
                                                     OFFSETS['T37T12_OFFSET_SEA_DAY'], 
                                                     OFFSETS, thr)+2.0 
        } 
    if onlyCirrus:
        title =info+"thinCirrusPrimaryTest_Cirrus_"+ SchemeName
    else:
        title =info+"thinCirrusPrimaryTest_All_ "+ SchemeName

    args_test = {'title':  title,
            'xlable': 'Latitude',
            'ylable': 'T37-T12 minus dynamic threshold',
            'xvector': thr.latitude,
            'yvector': thr.t37_t12_minus_threshold, 
            }
    THRESHOLD = OFFSETS_thinCirrusPrimaryTest['t37t12_OFFSETS'][SchemeName]
    if NEW_THRESHOLD is  not None:
        THRESHOLD = NEW_THRESHOLD
    if  args['USE_MARGINS']:
        THRESHOLD +=   OFFSETS['QUALITY_MARGIN_T37T12']
    TestOk = thr.t37_t12_minus_threshold>THRESHOLD                      
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, np.mean(THRESHOLD), 
              onlyCirrus=onlyCirrus, show=show)
    return TestOk

def HighcloudTestt85t11Sea(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                          info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_HighcloudTestt85t11Sea = {}
    OFFSETS_HighcloudTestt85t11Sea['t85t11_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['HI_T85T11_CM_SEA'], 
        'WaterNight': OFFSETS['HI_T85T11_CM_SEA'], 
        'WaterTwilight': OFFSETS['HI_T85T11_CM_SEA'] 
        } 
    if onlyCirrus:
        title =info+"HighcloudTestt85t11Sea_Cirrus_"+ SchemeName
    else:
        title =info+"HighcloudTestt85t11Sea_All_ "+ SchemeName

    args_test = {'title': title,
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': thr.latitude,
            'yvector': thr.t85_t11_minus_threshold,
            }
    THRESHOLD = OFFSETS_HighcloudTestt85t11Sea['t85t11_OFFSETS'][SchemeName]
    TestOk = thr.t85_t11_minus_threshold>THRESHOLD                      
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD, 
              onlyCirrus=onlyCirrus, show=show)
    return TestOk

def HighcloudTestt85t11land(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                          info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_HighcloudTestt85t11land = {}
    OFFSETS_HighcloudTestt85t11land['t85t11_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': OFFSETS['HI_T85T11_CM_LAND'], 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['HI_T85T11_CM_LAND'], 
        'LandDayMount': OFFSETS['HI_T85T11_CM_LAND'], 
        'LandNight':      OFFSETS['HI_T85T11_CM_LAND'],
        'LandNightInv':   OFFSETS['HI_T85T11_CM_LAND'], 
        'LandNightMount': OFFSETS['HI_T85T11_CM_LAND'], 
        'LandTwilight': OFFSETS['HI_T85T11_CM_LAND'], 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 100, 
        'WaterTwilight': 100 
        } 
    if onlyCirrus:
        title =info+"HighcloudTestt85t11land_Cirrus_"+ SchemeName
    else:
        title =info+"HighcloudTestt85t11land_All_ "+ SchemeName

    args_test = {'title': title,
            'xlable': 'Latitude',
            'ylable': 'T85-T11 minus dynamic threshold',
            'xvector': thr.latitude,
            'yvector': thr.t85_t11_minus_threshold,
            }
    THRESHOLD = OFFSETS_HighcloudTestt85t11land['t85t11_OFFSETS'][SchemeName]
    TestOk = thr.t85_t11_minus_threshold>THRESHOLD                      
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond,TestOk)
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD, 
              onlyCirrus=onlyCirrus, show=show)
    return TestOk

def textureNightTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                          info="", REMOVE_2times_THR=False, NEW_THRESHOLD_T37T12=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_textureNight = {}
    OFFSETS_textureNight['t11text_OFFSETS'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11TEXT_OFFSET_SEA_DAY'], 
        'WaterNight': OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': OFFSETS['T11TEXT_OFFSET_SEA_DAY'] 
        }
    OFFSETS_textureNight['t37t12text_OFFSETS'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T37T12TEXT_OFFSET_SEA_DAY'], 
        'WaterNight': OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': OFFSETS['T37T12TEXT_OFFSET_SEA_DAY']
        }
    OFFSETS_textureNight['tsur_threshold'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['COLDEST_SEASURFACE_TEMP'],#OFFSETS['TSUR_THR_SNOW_MAX'], 
        'WaterNight': 0,
        'WaterTwilight':OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        }
    if onlyCirrus:
        title =info+"textureNight_Cirrus_"+ SchemeName
    else:
        title =info+"textureNight_All_ "+ SchemeName
    args_test = {'title': title,
                 'ylable': 'T37T12text minus dynamic threshold',
                 'xlable': 'T11text minus dynamic threshold',
                 'yvector': thr.t37t12text,
                 'xvector': thr.t11text,
                 }
    THRESHOLD2 = OFFSETS_textureNight['t11text_OFFSETS'][SchemeName]
    THRESHOLD1 = OFFSETS_textureNight['t37t12text_OFFSETS'][SchemeName]
    TSUR_THRESHOLD = OFFSETS_textureNight['tsur_threshold'][SchemeName]
    if NEW_THRESHOLD_T37T12 is  not None:
        THRESHOLD1 = NEW_THRESHOLD_T37T12
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_T11TEXT']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_T37T12TEXT']
    if  REMOVE_2times_THR:
        THRESHOLD2 -= 2*OFFSETS['QUALITY_MARGIN_T11TEXT']
        THRESHOLD1 -= 2*OFFSETS['QUALITY_MARGIN_T37T12TEXT']
    TestOk = np.logical_and(np.logical_and(thr.t11text>THRESHOLD2, 
                                           thr.t37t12text>THRESHOLD1),
                            thr.surftemp>TSUR_THRESHOLD)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def textureIrVisTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                          info="", REMOVE_2times_THR=False, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_textureIrVisTest = {}
    OFFSETS_textureIrVisTest['t11text_OFFSETS'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11TEXT_OFFSET_SEA_DAY'], 
        'WaterNight': 100, 
        'WaterTwilight': OFFSETS['T11TEXT_OFFSET_SEA_DAY'],#OFFSETS['T11TEXT_OFFSET_SUNGLINT'], 
        }
    OFFSETS_textureIrVisTest['r06text_OFFSETS'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06TEXT_OFFSET_SEA'], 
        'WaterNight': 100, 
        'WaterTwilight': OFFSETS['R06TEXT_OFFSET_SEA'],#OFFSETS['R06TEXT_OFFSET_SUNGLINT'], 
        }
    OFFSETS_textureIrVisTest['tsur_threshold'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['COLDEST_SEASURFACE_TEMP'],#OFFSETS['TSUR_THR_SNOW_MAX'], 
        'WaterNight': 100, 
        'WaterTwilight': OFFSETS['COLDEST_SEASURFACE_TEMP']
        }
    OFFSETS_textureIrVisTest['sunz_threshold'] = {         
        'IceDay': 800, 
        'IceNight': 800, 
        'IceTwilight': 800,  
        'SunglintDay': 800, 
        'SunglintTwilight': 800, 
        'WaterDay': 800, 
        'WaterNight': 800, 
        'WaterTwilight': OFFSETS['MAX_SUNZEN_TWILIGHT_VIS']
        }
    if onlyCirrus:
        title =info+"textureIrVisTest_Cirrus_"+ SchemeName
    else:
        title =info+"textureIrVisTest_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'r06text minus dynamic threshold',
                 'ylable': 'T11text minus dynamic threshold',
                 'yvector': thr.t11text,
                 'xvector': thr.r06text,
                 }
    THRESHOLD1 = OFFSETS_textureIrVisTest['t11text_OFFSETS'][SchemeName]
    THRESHOLD2 = OFFSETS_textureIrVisTest['r06text_OFFSETS'][SchemeName]
    TSUR_THRESHOLD = OFFSETS_textureIrVisTest['tsur_threshold'][SchemeName]
    SUNZ_THRESHOLD = OFFSETS_textureIrVisTest['sunz_threshold'][SchemeName]

    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_R06TEXT']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_T11TEXT']
    if   REMOVE_2times_THR:
        THRESHOLD1 -= 2*OFFSETS['QUALITY_MARGIN_T11TEXT']
    TestOk = np.logical_and(np.logical_and(thr.t11text>THRESHOLD1, 
                                           thr.r06text>THRESHOLD2),
                            np.logical_and(thr.sunz<SUNZ_THRESHOLD,
                                           thr.surftemp> TSUR_THRESHOLD))
    #TestOk = np.logical_and(TestOk, 
    #                        ((thr.t11text-THRESHOLD1)**2+
    #                        (thr.r06text-THRESHOLD2)**2)>4)#(np.min(THRESHOLD1,THRESHOLD2)+1.0))
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def newClearWaterBodiesTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args, info="",show=False):
    r06_threshold =   100.0*thr.thr_r06*OFFSETS['R06_GAIN_LAND'] + OFFSETS['R06_OFFSET_LAND']
    r06_threshold = np.ma.array(
        r06_threshold, 
        mask=np.logical_or(thr.r06.mask,thr.thr_r06.mask))
    title =info+"newClearWaterBodiesTest_All"+ SchemeName
    args_test = {'title': title,
                 'xlable': 'T37-T12',
                 'ylable': 'qr09r06',
                 'xvector': thr.t37t12,
                 'yvector': thr.qr09r06,
                 }
    THRESHOLD1 = 0.8
    THRESHOLD2 = r06_threshold
    THRESHOLD3 = 3.0
    THRESHOLD4 = 3.0
    TestOk = np.logical_and(np.logical_and(thr.qr09r06<THRESHOLD1,
                                           thr.r06<THRESHOLD2),
                            np.logical_and(np.abs(thr.t37t12)<THRESHOLD3,
                                           np.abs(thr.t11t37)<THRESHOLD4))
    TestOk = np.logical_and(TestOk, caobj.avhrr.all_arrays['bt37micron']>-9)
    plot_test_clear(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD1,  show=show)
    return TestOk

def brightCloudTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_brightCloudTest = {}
    OFFSETS_brightCloudTest['t37t12_OFFSETS'] = { 
        'CoastDay': 12.0, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 5.0, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['T37T12_OFFSET_LAND_DAY'], 
        'LandDayMount': OFFSETS['T37T12_OFFSET_LAND_DAY']+3.0, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T37T12_OFFSET_SEA_DAY'], 
        'WaterTwilight': OFFSETS['T37T12_OFFSET_SEA_DAY'] #non existing
        }     
    OFFSETS_brightCloudTest['r06gain_OFFSETS'] = { 
        'CoastDay': OFFSETS['R06_GAIN_SUNGLINT'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': OFFSETS['R06_GAIN_SEA_OPAQUE'], 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_GAIN_LAND'], 
        'LandDayMount': OFFSETS['R06_GAIN_LAND_OPAQUE'], #try opaque limit
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_GAIN_SEA_OPAQUE'],#OFFSETS['R06_GAIN_SEA'], 
        'WaterTwilight': OFFSETS['R06_GAIN_SEA_OPAQUE']
        } 
    OFFSETS_brightCloudTest['r06offset_OFFSETS'] = { 
        'CoastDay': OFFSETS['R06_OFFSET_SUNGLINT'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': OFFSETS['R06_OFFSET_SEA_OPAQUE'],
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_OFFSET_LAND'], 
        'LandDayMount': OFFSETS['R06_OFFSET_LAND_OPAQUE'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_OFFSET_SEA'],#OFFSETS['R06_OFFSET_SEA'],  
        'WaterTwilight': OFFSETS['R06_OFFSET_SEA_OPAQUE']
        } 

    R06_OFFSET = OFFSETS_brightCloudTest['r06offset_OFFSETS'][SchemeName]

    R06_GAIN = OFFSETS_brightCloudTest['r06gain_OFFSETS'][SchemeName]

    r06_threshold =   100.0*thr.thr_r06*R06_GAIN + R06_OFFSET
    r06_threshold = np.ma.array(
        r06_threshold, 
        mask=np.logical_or(thr.r06.mask,thr.thr_r06.mask))
    if onlyCirrus:
        title =info+"brightCloudTest_Cirrus_"+ SchemeName
    else:
        title =info+"brightCloudTest_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'T37-T12 minus dynamic threshold',
                 'ylable': 'r06 ',
                 'xvector': thr.t37_t12_minus_threshold,
                 'yvector': thr.r06,
                 }
    THRESHOLD2 = OFFSETS_brightCloudTest['t37t12_OFFSETS'][SchemeName]
    THRESHOLD1 = r06_threshold  
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_T37T12']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_R06']
    #if scheme in ['LandDayMount']:
    #    THRESHOLD1=THRESHOLD1+10
    TestOk =   np.logical_and(thr.t37_t12_minus_threshold>THRESHOLD2,
                              thr.r06>THRESHOLD1)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                np.mean(THRESHOLD1), THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def brightCloudTest3A(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    satsec = 1.0/np.cos(np.radians(caobj.avhrr.satz))
    threshold_r16= np.where(
        caobj.avhrr.azidiff > 90.0,
        OFFSETS['R16_THR_SNOW_OFFSET'] + OFFSETS['R16_THR_SNOW_SATZ_FORWARD_GAIN'] * satsec,
        OFFSETS['R16_THR_SNOW_OFFSET'] + OFFSETS['R16_THR_SNOW_SATZ_BACKWARD_GAIN'] * satsec)
   

    OFFSETS_brightCloudTest3A = {}   
    OFFSETS_brightCloudTest3A['r06gain_OFFSETS'] = { 
        'All': OFFSETS['R06_GAIN_LAND_OPAQUE'],
        'CoastDay': OFFSETS['R06_GAIN_LAND_OPAQUE'],
        'CoastDayMount': OFFSETS['R06_GAIN_LAND'], 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_GAIN_LAND_OPAQUE'],
        'LandDayMount':100, #try opaque limit
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': OFFSETS['R06_GAIN_LAND'], 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_GAIN_SEA_OPAQUE'],
        'WaterTwilight': 100
        } 
    OFFSETS_brightCloudTest3A['r06offset_OFFSETS'] = {
        'All': OFFSETS['R06_OFFSET_LAND_OPAQUE'],
        'CoastDay': OFFSETS['R06_OFFSET_LAND_OPAQUE'],
        'CoastDayMount': OFFSETS['R06_OFFSET_LAND'], 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_OFFSET_LAND_OPAQUE'],
        'LandDayMount': 100,
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': OFFSETS['R06_OFFSET_LAND'], 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_OFFSET_SEA_OPAQUE'],
        'WaterTwilight': 100,
        }

    R06_OFFSET = OFFSETS_brightCloudTest3A['r06offset_OFFSETS'][SchemeName]

    R06_GAIN = OFFSETS_brightCloudTest3A['r06gain_OFFSETS'][SchemeName]

    r06_threshold =   100.0*thr.thr_r06*R06_GAIN + R06_OFFSET
    r06_threshold = np.ma.array(
        r06_threshold, 
        mask=np.logical_or(thr.r06.mask,thr.thr_r06.mask))
    if onlyCirrus:
        title =info+"brightCloudTest3A_Cirrus_"+ SchemeName
    else:
        title =info+"brightCloudTest3A_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'r16',
                 'ylable': 'r06 ',
                 'xvector': thr.r16,
                 'yvector': thr.r06,
                 }
    THRESHOLD2 = threshold_r16
    THRESHOLD1 = r06_threshold*(1+thr.r16/50.0)  
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_R16']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_R06']
    TestOk =   np.logical_and(thr.r16>THRESHOLD2,
                              thr.r06>THRESHOLD1)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                np.mean(THRESHOLD1), np.mean(THRESHOLD2), onlyCirrus=onlyCirrus, show=show)
    return TestOk


def brightCloudTestSea(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  OTHER_T11_THR=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_brightCloudTestSea = {}
    OFFSETS_brightCloudTestSea['t11_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': OFFSETS['T11_SEA_MIN'] -5.0, 
        'IceTwilight': 100, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11_SEA_MIN'] -5.0,
        'WaterTwilight': 100 
        }     
    OFFSETS_brightCloudTestSea['r06gain_OFFSETS'] = { 
        'CoastDay':100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,
        'IceTwilight': 100, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_GAIN_SEA'], 
        'WaterTwilight': 100 
        } 
    OFFSETS_brightCloudTestSea['r06offset_OFFSETS'] = { 
        'CoastDay': OFFSETS['R06_OFFSET_SUNGLINT'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': OFFSETS['R06_OFFSET_SEA'],
        'IceTwilight': 100,  
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_OFFSET_SEA'],  
        'WaterTwilight': 100 
        } 

    R06_OFFSET = OFFSETS_brightCloudTestSea['r06offset_OFFSETS'][SchemeName]

    R06_GAIN = OFFSETS_brightCloudTestSea['r06gain_OFFSETS'][SchemeName]

    r06_threshold =   100.0*thr.thr_r06*R06_GAIN + R06_OFFSET
    r06_threshold = np.ma.array(
        r06_threshold, 
        mask=np.logical_or(thr.r06.mask,thr.thr_r06.mask))
    if onlyCirrus:
        title =info+"brightCloudTestSea_Cirrus_"+ SchemeName
    else:
        title =info+"brightCloudTestSea_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'T11',
                 'ylable': 'r06 ',
                 'xvector': caobj.avhrr.all_arrays['bt11micron'],
                 'yvector': thr.r06,
                 }

    THRESHOLD2 = OFFSETS_brightCloudTestSea['t11_OFFSETS'][SchemeName]
    THRESHOLD1 = r06_threshold  
    if OTHER_T11_THR:
        THRESHOLD2 = OTHER_T11_THR 
    if  args['USE_MARGINS']:
        THRESHOLD2 -= OFFSETS['QUALITY_MARGIN_T11']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_R06']
    TestOk =   np.logical_and(caobj.avhrr.all_arrays['bt11micron']<THRESHOLD2,
                              thr.r06>THRESHOLD1)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                np.mean(THRESHOLD1), THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk



def brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="", ADD_TO_ROG_THR=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_brightCloudTestNoSunglint3A = {}  
     
    OFFSETS_brightCloudTestNoSunglint3A['r06gain_OFFSETS'] = { 
        'CoastDay': 0.5*(OFFSETS['R06_GAIN_LAND']+OFFSETS['R06_GAIN_SEA']), 
        'CoastDayMount': 0.5*(OFFSETS['R06_GAIN_LAND']+OFFSETS['R06_GAIN_SEA']), 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_GAIN_LAND'], 
        'LandDayMount': OFFSETS['R06_GAIN_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_GAIN_SEA'],  
        'WaterTwilight': 100 
        } 
    OFFSETS_brightCloudTestNoSunglint3A['r06offset_OFFSETS'] = { 
        'CoastDay': 0.5*(OFFSETS['R06_OFFSET_LAND']+OFFSETS['R06_OFFSET_SEA']), 
        'CoastDayMount': 0.5*(OFFSETS['R06_OFFSET_LAND']+OFFSETS['R06_OFFSET_SEA']),  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_OFFSET_LAND'], 
        'LandDayMount': OFFSETS['R06_OFFSET_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_OFFSET_SEA'], 
        'WaterTwilight': 100 
        } 

    R06_OFFSET = OFFSETS_brightCloudTestNoSunglint3A['r06offset_OFFSETS'][SchemeName]

    R06_GAIN = OFFSETS_brightCloudTestNoSunglint3A['r06gain_OFFSETS'][SchemeName]

    r06_threshold =   100.0*thr.thr_r06*R06_GAIN + R06_OFFSET
    r06_threshold = np.ma.array(
        r06_threshold, 
        mask=np.logical_or(thr.r06.mask,thr.thr_r06.mask))
    if onlyCirrus:
        title =info+"brightCloudTestNoSunglint3A_Cirrus_"+ SchemeName
    else:
        title =info+"brightCloudTestNoSunglint3A_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'r06 ',
                 'ylable': 'qr16r06',
                 'xvector': thr.r06,
                 'yvector': thr.qr16r06,
            }
    THRESHOLD2 = r06_threshold  
    THRESHOLD1 = OFFSETS['R06_R16_THR_SNOW'],
    THRESHOLD3 = 1.35 
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_R06']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_QR16R06']
        THRESHOLD3 -= OFFSETS['QUALITY_MARGIN_QR16R06']
    if ADD_TO_ROG_THR is not None:
        THRESHOLD2 += ADD_TO_ROG_THR
    TestOk = np.logical_and(thr.r06>THRESHOLD2,
                            np.logical_and(thr.qr16r06>THRESHOLD1,
                                           thr.qr16r06<THRESHOLD3))
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                THRESHOLD1, np.mean(THRESHOLD2), onlyCirrus=onlyCirrus, show=show)
    return TestOk


def coldBrightCloudTest37(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_coldBrightCloudTest37 = {}
    OFFSETS_coldBrightCloudTest37['t37t12_OFFSETS'] = { 
        'CoastDay': 7.0, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 6.0, 
        'LandDayMount': 8.0 ,  
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        }     
    OFFSETS_coldBrightCloudTest37['t11ts_OFFSETS'] = { 
        'CoastDay': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'], 
        'LandDayMount': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        }          
    OFFSETS_coldBrightCloudTest37['r06gain_OFFSETS'] = { 
        'CoastDay': OFFSETS['R06_GAIN_COAST'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_GAIN_LAND'], 
        'LandDayMount': OFFSETS['R06_GAIN_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldBrightCloudTest37['r06offset_OFFSETS'] = { 
        'CoastDay': OFFSETS['R06_OFFSET_COAST'], 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_OFFSET_LAND'], 
        'LandDayMount': OFFSETS['R06_OFFSET_LAND'], 
        'LandNight': 100, 
        'LandNightInv': 100, 
        'LandNightMount': 100, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldBrightCloudTest37['tsur_threshold'] = { 
        'CoastDay': OFFSETS['TSUR_THR_SNOW_MAX'], 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['TSUR_THR_SNOW_MAX'], 
        'LandDayMount': OFFSETS['TSUR_THR_SNOW_MAX'], 
        'LandNight': 100, 
        'LandNightInv': 100, 
        'LandNightMount': 100, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 100, 
        'WaterTwilight': 100 
        }

    R06_OFFSET = OFFSETS_coldBrightCloudTest37['r06offset_OFFSETS'][SchemeName]
    R06_GAIN = OFFSETS_coldBrightCloudTest37['r06gain_OFFSETS'][SchemeName]
    r06_threshold =   100.0*thr.thr_r06*R06_GAIN + R06_OFFSET
    r06_threshold = np.ma.array(
        r06_threshold, 
        mask=np.logical_or(thr.r06.mask,thr.thr_r06.mask))

    if onlyCirrus:
        title =info+"coldBrightCloudTest37_Cirrus_"+ SchemeName
    else:
        title =info+"coldBrightCloudTest37_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'R06',
                 'ylable': 'T11-Ts minus dynamic threshold',
                 'xvector': thr.r06,
                 'yvector': thr.t11_ts_minus_threshold,
            }
    THRESHOLD2 = r06_threshold  
    THRESHOLD1 = OFFSETS_coldBrightCloudTest37['t11ts_OFFSETS'][SchemeName]
    THRESHOLD3 = OFFSETS_coldBrightCloudTest37['t37t12_OFFSETS'][SchemeName]
    TSUR_THRESHOLD = OFFSETS_coldBrightCloudTest37['tsur_threshold'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_R06']
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR']
        THRESHOLD3 += OFFSETS['QUALITY_MARGIN_T37T12']
    TestOk = np.logical_and(
        np.logical_and(thr.t11_ts_minus_threshold<THRESHOLD1,
                       thr.surftemp>TSUR_THRESHOLD),
        np.logical_and(thr.r06>THRESHOLD2,
                       thr.t37_t12_minus_threshold>THRESHOLD3))
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), onlyCirrus=onlyCirrus, 
                    show=show)
    return TestOk

def coldBrightCloudTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_coldBrightCloudTest = {}  
    OFFSETS_coldBrightCloudTest['t11ts_OFFSETS'] = { 
        'CoastDay': OFFSETS['T11_OFFSET_LAND_DAY'], #test currentlynot used for coast_day 
        'CoastDayMount': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'],
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,  
        'IceTwilight': 100, 
        'LandDay': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'], 
        'LandDayMount': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'], #use opaque isntead 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11_OFFSET_SEA_DAY'], 
        'WaterTwilight': 100 
        }          
    OFFSETS_coldBrightCloudTest['r06gain_OFFSETS'] = { 
        'CoastDay': OFFSETS['R06_GAIN_LAND_OPAQUE'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_GAIN_LAND_OPAQUE'], 
        'LandDayMount': OFFSETS['R06_GAIN_LAND_OPAQUE'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_GAIN_SEA_OPAQUE'], 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldBrightCloudTest['r06offset_OFFSETS'] = { 
        'CoastDay': OFFSETS['R06_OFFSET_LAND_OPAQUE'], 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R06_OFFSET_LAND_OPAQUE'], 
        'LandDayMount': OFFSETS['R06_OFFSET_LAND_OPAQUE'],  
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R06_OFFSET_SEA_OPAQUE'], 
        'WaterTwilight': 100 
        }
    OFFSETS_coldBrightCloudTest['tsur_threshold'] = { 
        'CoastDay': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 0, 
        'LandDayMount': 0,
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'WaterTwilight': 100 
        } 

    R06_OFFSET = OFFSETS_coldBrightCloudTest['r06offset_OFFSETS'][SchemeName]
    R06_GAIN = OFFSETS_coldBrightCloudTest['r06gain_OFFSETS'][SchemeName]
    
    r06_threshold =   100.0*thr.thr_r06*R06_GAIN + R06_OFFSET
    r06_threshold = np.ma.array(
        r06_threshold, 
        mask=np.logical_or(thr.r06.mask,thr.thr_r06.mask))

    if onlyCirrus:
        title =info+"coldBrightCloudTest_Cirrus_"+ SchemeName
    else:
        title =info+"coldBrightCloudTest_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'R06',
                 'ylable': 'T11-Ts minus dynamic threshold',
                 'xvector': thr.r06,
                 'yvector': thr.t11_ts_minus_threshold,
            }
    THRESHOLD2 = r06_threshold  
    THRESHOLD1 = OFFSETS_coldBrightCloudTest['t11ts_OFFSETS'][SchemeName]
    TSUR_THRESHOLD = OFFSETS_coldBrightCloudTest['tsur_threshold'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_R06']
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR']
        TSUR_THRESHOLD += OFFSETS['QUALITY_MARGIN_TSUR']
    TestOk = np.logical_and(thr.t11_ts_minus_threshold<THRESHOLD1,
                            thr.r06>THRESHOLD2)                    
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    THRESHOLD1, np.mean(THRESHOLD2), onlyCirrus=onlyCirrus, 
                    show=show)
    return TestOk


def thincoldCirrusTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  
                       info="", onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_thincoldCirrusTest = {}
    tmp_day=0.5*(OFFSETS['T11T12_OFFSET_SEA_NIGHT']+ OFFSETS['T11T12_OFFSET_SEA_DAY'])
    tmp_night=0.5*(OFFSETS['T11T12_OFFSET_LAND_NIGHT']+ OFFSETS['T11T12_OFFSET_LAND_DAY'])
    OFFSETS_thincoldCirrusTest['t11t12_OFFSETS'] = { 
        'CoastDay': np.max([ OFFSETS['T11T12_OFFSET_SEA_DAY'],
                             OFFSETS['T11T12_OFFSET_LAND_DAY']]), 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': thr_slope_between_90_and_day(tmp_night, tmp_day, OFFSETS, thr), 
        #np.max([0.5*(OFFSETS['T11T12_OFFSET_SEA_NIGHT']+
        #                              OFFSETS['T11T12_OFFSET_SEA_DAY']),
        #                         0.5*(OFFSETS['T11T12_OFFSET_LAND_NIGHT']+
        #                              OFFSETS['T11T12_OFFSET_LAND_DAY'])]), 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['T11T12_OFFSET_LAND_DAY'], 
        'LandDayMount': OFFSETS['T11T12_OFFSET_LAND_DAY'], 
        'LandNight': 100,
        'LandNightInv': 100,
        'LandNightMount':  100,
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': thr_slope_between_90_and_day(
            OFFSETS['T11T12_OFFSET_LAND_NIGHT'],
            OFFSETS['T11T12_OFFSET_LAND_DAY'],
            OFFSETS, thr), 
        #np.max([OFFSETS['T11T12_OFFSET_LAND_DAY'], 
         #       OFFSETS['T11T12_OFFSET_LAND_NIGHT']]),
        #test currentlynot used for land_twilight_inv
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_thincoldCirrusTest['t11ts_OFFSETS'] = { 
        'CoastDay': np.max([ OFFSETS['T11_OFFSET_SEA_DAY'],
                             OFFSETS['T11_OFFSET_LAND_DAY']]),
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': np.min([0.5*(OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE']+
                                      OFFSETS['T11_OFFSET_SEA_DAY_OPAQUE']),
                                 0.5*(OFFSETS['T11_OFFSET_LAND_NIGHT_OPAQUE']+
                                      OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'])]), 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['T11_OFFSET_LAND_DAY'], 
        'LandDayMount': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'], 
        'LandNight': 100,
        'LandNightInv': 100,
        'LandNightMount':  100,
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': np.max([OFFSETS['T11_OFFSET_LAND_DAY'], 
                                     OFFSETS['T11_OFFSET_LAND_NIGHT']]),#test currentlynot used for land_twilight_inv 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 100, 
        'WaterTwilight': 100 
        } 
    if onlyCirrus:
        title =info+"thincoldCirrusTest_Cirrus_"+ SchemeName
    else:
        title =info+"thincoldCirrusTest_All_ "+ SchemeName
    args_test = {'title': title ,
            'xlable': 'T11-TS minus dynamic threshold',
            'ylable': 'T11-T12 minus dynamic threshold',
            'xvector': thr.t11_ts_minus_threshold,
            'yvector': thr.t11_t12_minus_threshold, 
            }
    THRESHOLD1 = OFFSETS_thincoldCirrusTest['t11t12_OFFSETS'][SchemeName]
    THRESHOLD2 = OFFSETS_thincoldCirrusTest['t11ts_OFFSETS'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD2 -= OFFSETS['QUALITY_MARGIN_T11TSUR']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_T11T12']
    TestOk = np.logical_and(thr.t11_t12_minus_threshold>THRESHOLD1,
                            thr.t11_ts_minus_threshold<THRESHOLD2)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    np.mean(THRESHOLD1), THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def coldWatercloudTestDay(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  
                               info="", onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_coldWatercloudTestDay = {} 
    OFFSETS_coldWatercloudTestDay['t11ts_OFFSETS'] = { 
        'CoastDay': np.min([OFFSETS['T11_OFFSET_LAND_DAY'],
                            OFFSETS['T11_OFFSET_SEA_DAY']]), 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['T11_OFFSET_LAND_DAY'], 
        'LandDayMount': OFFSETS['T11_OFFSET_MOUNTAIN_DAY'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['T11_OFFSET_SEA_DAY'], 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldWatercloudTestDay['t11t37_OFFSETS'] = { 
        'CoastDay': 0.0, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 0.0, 
        'LandDayMount': 0.0, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 0.0, 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldWatercloudTestDay['r0609_clo_max_OFFSETS'] = { 
        'CoastDay': OFFSETS['R09_R06_THR_CLOUDY_MAX'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R09_R06_THR_CLOUDY_MAX'], 
        'LandDayMount': OFFSETS['R09_R06_THR_CLOUDY_MAX'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R09_R06_THR_CLOUDY_MAX'], 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldWatercloudTestDay['r0609_clo_min_OFFSETS'] = { 
        'CoastDay': OFFSETS['R09_R06_THR_CLOUDY_MIN'], 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': OFFSETS['R09_R06_THR_CLOUDY_MIN'], 
        'LandDayMount': OFFSETS['R09_R06_THR_CLOUDY_MIN'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': OFFSETS['R09_R06_THR_CLOUDY_MIN'],
        'WaterTwilight': 100 
        } 

    if onlyCirrus:
        title =info+"coldWatercloudTestDay_Cirrus_"+ SchemeName
    else:
        title =info+"coldWatercloudTestDay_All_ "+ SchemeName

    args_test = {'title':title,
                 'xlable': 'T11-T37 minus dynamic threshold',
                 'ylable': 'T11-Ts minus dynamic threshold',
                 'xvector': thr.t11_t37_minus_threshold,
                 'yvector': thr.t11_ts_minus_threshold,  
                 }
    THRESHOLD1 = OFFSETS_coldWatercloudTestDay['t11ts_OFFSETS'][SchemeName]
    THRESHOLD2 = OFFSETS_coldWatercloudTestDay['t11t37_OFFSETS'][SchemeName]    
    THRESHOLD3 = 298.5
    THRESHOLD4 = OFFSETS_coldWatercloudTestDay['r0609_clo_max_OFFSETS'][SchemeName]
    THRESHOLD5 = OFFSETS_coldWatercloudTestDay['r0609_clo_min_OFFSETS'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR']
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_T11T37']
        THRESHOLD3 -= 1.0
        THRESHOLD4 += OFFSETS['QUALITY_MARGIN_QR09R06']
        THRESHOLD5 -= OFFSETS['QUALITY_MARGIN_QR09R06']
    qr0906cond=np.logical_and( thr.qr09r06<THRESHOLD4,
                               thr.qr09r06>THRESHOLD5)
    TestOk = np.logical_and(
        np.logical_and(caobj.avhrr.all_arrays['bt11micron']<THRESHOLD3,
                       qr0906cond),
        np.logical_and( thr.t11_t37_minus_threshold>THRESHOLD2,
                        thr.t11_ts_minus_threshold<THRESHOLD1))
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def reflectingCloudTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_reflectingCloudTest = {}  
    OFFSETS_reflectingCloudTest['t37t12_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,  
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': thr_slope_between_90_and_day(OFFSETS['T37T12_OFFSET_LAND_NIGHT'],
                                                     OFFSETS['T37T12_OFFSET_LAND_DAY'],
                                                     OFFSETS, thr), 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': thr_slope_between_90_and_day(OFFSETS['T37T12_OFFSET_SEA_NIGHT'],
                                                      OFFSETS['T37T12_OFFSET_SEA_DAY'],
                                                     OFFSETS, thr),

            #0.5*(OFFSETS['T37T12_OFFSET_SEA_NIGHT'] + 
            #                  OFFSETS['T37T12_OFFSET_SEA_DAY'])
        }          
    OFFSETS_reflectingCloudTest['vis_gain_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': OFFSETS['VIS_STATIC_LAND_GAIN'], 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['VIS_STATIC_SEA_GAIN'] 
        } 
    OFFSETS_reflectingCloudTest['vis_offset_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': OFFSETS['VIS_STATIC_LAND_OFFSET'], 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['VIS_STATIC_SEA_OFFSET'] 
        } 

    VIS_OFFSET = OFFSETS_reflectingCloudTest['vis_offset_OFFSETS'][SchemeName]
    VIS_GAIN = OFFSETS_reflectingCloudTest['vis_gain_OFFSETS'][SchemeName]
    vis_static = np.where(
        thr.sunz>=90.0,
        VIS_OFFSET,
        VIS_GAIN*(90.0-thr.sunz) + VIS_OFFSET)

    if onlyCirrus:
        title =info+"reflectingCloudTest_Cirrus_"+ SchemeName
    else:
        title =info+"reflectingCloudTest_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'r06 not scaled',
                 'ylable': 'T37-T2 minus dynamic threshold',
                 'xvector': thr.pseudor06,
                 'yvector': thr.t37_t12_minus_threshold,
            }
    THRESHOLD2 = vis_static
    THRESHOLD1 = OFFSETS_reflectingCloudTest['t37t12_OFFSETS'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_VIS_STATIC']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_T37T12']

    TestOk = np.logical_and(thr.t37_t12_minus_threshold>THRESHOLD1,
                            thr.pseudor06>THRESHOLD2)

    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    np.mean(THRESHOLD1), np.mean(THRESHOLD2), onlyCirrus=onlyCirrus, 
                    show=show)
    return TestOk
def reflectingCloudTestSea(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  NEW_T11_THR=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_reflectingCloudTestSea = {}  
    OFFSETS_reflectingCloudTestSea['t11_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,  
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': 100,
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['T11_SEA_MIN'] -3.0
        }          
    OFFSETS_reflectingCloudTestSea['vis_gain_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': 100,
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['VIS_STATIC_SEA_GAIN'] 
        } 
    OFFSETS_reflectingCloudTestSea['vis_offset_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': 100,
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['VIS_STATIC_SEA_OFFSET'] 
        } 

    VIS_OFFSET = OFFSETS_reflectingCloudTestSea['vis_offset_OFFSETS'][SchemeName]
    VIS_GAIN = OFFSETS_reflectingCloudTestSea['vis_gain_OFFSETS'][SchemeName]
    vis_static = np.where(
        thr.sunz>=90.0,
        VIS_OFFSET,
        VIS_GAIN*(90.0-thr.sunz) + VIS_OFFSET)

    if onlyCirrus:
        title =info+"reflectingCloudTestSea_Cirrus_"+ SchemeName
    else:
        title =info+"reflectingCloudTestSea_All_ "+ SchemeName
    args_test = {'title': title,
                 'ylable': 't11',
                 'xlable': 'T37-T2 minus dynamic threshold',
                 'xvector': thr.pseudor06,
                 'yvector': caobj.avhrr.all_arrays['bt11micron'],
            }
    THRESHOLD2 = vis_static
    THRESHOLD1 = OFFSETS_reflectingCloudTestSea['t11_OFFSETS'][SchemeName]
    if NEW_T11_THR is not None:
        THRESHOLD1 = NEW_T11_THR
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_VIS_STATIC']
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11']

    TestOk = np.logical_and(caobj.avhrr.all_arrays['bt11micron']<THRESHOLD1,
                            thr.pseudor06>THRESHOLD2)

    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    np.mean(THRESHOLD1), np.mean(THRESHOLD2), onlyCirrus=onlyCirrus, 
                    show=show)
    return TestOk


def pseudo06CloudTest3A(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_pseudo06CloudTest3A = {}  
    OFFSETS_pseudo06CloudTest3A['qr06r16_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,  
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': OFFSETS['R06_R16_THR_LOWSUN'],
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['R06_R16_THR_LOWSUN']
        }          
    OFFSETS_pseudo06CloudTest3A['vis_gain_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': OFFSETS['VIS_STATIC_LAND_GAIN'], 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['VIS_STATIC_SEA_OFFSET']
        } 
    OFFSETS_pseudo06CloudTest3A['vis_offset_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100,
        'LandDayMount': 100,
        'LandTwilight': OFFSETS['VIS_STATIC_LAND_OFFSET'], 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': OFFSETS['VIS_STATIC_SEA_GAIN'] 
        } 

    VIS_OFFSET = OFFSETS_pseudo06CloudTest3A['vis_offset_OFFSETS'][SchemeName]
    VIS_GAIN = OFFSETS_pseudo06CloudTest3A['vis_gain_OFFSETS'][SchemeName]
    vis_static = np.where(
        thr.sunz>=90.0,
        VIS_OFFSET,
        VIS_GAIN*(90.0-thr.sunz) + VIS_OFFSET)

    if onlyCirrus:
        title =info+"pseudo06CloudTest3A_Cirrus_"+ SchemeName
    else:
        title =info+"pseudo06CloudTest3A_All_ "+ SchemeName
    args_test = {'title': title,
                 'xlable': 'r06 not scaled',
                 'ylable': 'qr16r06',
                 'xvector': thr.pseudor06,
                 'yvector': thr.qr16r06,
                 }
    THRESHOLD2 = vis_static
    THRESHOLD1 = OFFSETS_pseudo06CloudTest3A['qr06r16_OFFSETS'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_VIS_STATIC']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_QR16R06']

    TestOk = np.logical_and(thr.qr16r06>THRESHOLD1,
                            thr.pseudor06>THRESHOLD2)

    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    np.mean(THRESHOLD1), np.mean(THRESHOLD2), onlyCirrus=onlyCirrus, 
                    show=show)
    return TestOk
