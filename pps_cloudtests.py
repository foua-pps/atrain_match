import numpy as np
from pps_threshold_functions import (plot_test,
                                     plot_test_2_lim,
                                     plot_test_clear)


def coldCloudTest_v2014(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  info="", onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_coldCloudTest_v2014 = {}
    OFFSETS_coldCloudTest_v2014['t11ts_OFFSETS'] = { 
        'All': -4,
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': -11, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': -9, 
        'LandNight': OFFSETS['T11_OFFSET_LAND_NIGHT'], 
        'LandNightInv': -11, 
        'LandNightMount': -8, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': -7, 
        'WaterTwilight': 100 
        }
    OFFSETS_coldCloudTest_v2014['tsur_threshold'] = {  
        'All': 0,
        'CoastDay': 0, 
        'CoastDayMount': 0, 
        'CoastNight': 0, 
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
        'LandTwilightMount': 0, 
        'SunglintDay': 0, 
        'SunglintTwilight': 0, 
        'WaterDay': 0, 
        'WaterNight': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'WaterTwilight': 0        }
    if onlyCirrus:
        title =info+"coldCloudTest_v2014_Cirrus_"+ SchemeName
    else:
        title =info+"coldCloudTest_v2014_All_"+ SchemeName

    args_test = {'title':  title,
                 'xlable': 'T11-T37 minus dynamic threshold',
                 'ylable': 'T11-Ts minus dynamic threshold',
                 'xvector': thr.t11_t37_minus_threshold,
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
        'LandNightInv': OFFSETS['T11_OFFSET_LAND_NIGHT'], 
        'LandNightMount': OFFSETS['T11_OFFSET_MOUNTAIN_NIGHT'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['T11_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': 100 
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
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 0.0, 
        'WaterTwilight': 100 
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
        'LandNightInv':  OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'] + 3.0, 
        'LandNightMount': OFFSETS['T11T37_WATERCLOUD_SECURITY_OFFSET'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 0, 
        'WaterTwilight': 100 
        }

    if onlyCirrus:
        title =info+"watercloudTest_Cirrus_"+ SchemeName
    else:
        title =info+"watercloudTest_All_ "+ SchemeName

    args_test = {'title':title,
                 'xlable': 'Latitude',
                 'ylable': 'T11-T37 minus dynamic threshold',
                 'xvector': thr.latitude,
                 'yvector': thr.t11_t37_minus_threshold,   
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
                    THRESHOLD, onlyCirrus=onlyCirrus, show=show)
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
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['T11_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': 100 
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
        'WaterDay': 0, 
        'WaterNight': OFFSETS['COLDEST_SEASURFACE_TEMP'], 
        'WaterTwilight': 0        }
    if onlyCirrus:
        title =info+"coldCloudTest_tsurf_lim_Cirrus_"+ SchemeName
    else:
        title =info+"coldCloudTest_tsurf_lim_All_ "+ SchemeName
    args_test = {'title': title,
            'xlable': 'Latitude',
            'ylable': 'T11-Ts minus dynamic threshold',
            'xvector': thr.latitude,
            'yvector': thr.t11_ts_minus_threshold,
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
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['T11_OFFSET_SEA_NIGHT_OPAQUE'], 
        'WaterTwilight': 100 
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
    OFFSETS_thinCirrusSecondaryTest = {}
    OFFSETS_thinCirrusSecondaryTest['t11t12_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': OFFSETS['T11T12_OFFSET_LAND_NIGHT'], 
        'LandNightInv': OFFSETS['T11T12_OFFSET_LAND_NIGHT'], 
        'LandNightMount':  OFFSETS['T11T12_OFFSET_LAND_NIGHT'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 100, 
        'WaterTwilight': 100 
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
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD, 
              onlyCirrus=onlyCirrus, show=show)
    return TestOk

def thinCirrusPrimaryTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                          info="", NEW_THRESHOLD=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_thinCirrusPrimaryTest = {}
    OFFSETS_thinCirrusPrimaryTest['t37t12_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': np.min([OFFSETS['T37T12_OFFSET_SEA_NIGHT'],
                                 OFFSETS['T37T12_OFFSET_LAND_NIGHT']])-1.0, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight': OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -2.0,
        'LandNightInv': OFFSETS['T37T12_OFFSET_LAND_NIGHT'] -1.0, 
        'LandNightMount': OFFSETS['T37T12_OFFSET_LAND_NIGHT']-1.0, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['T37T12_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': 100 
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
    plot_test(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD, 
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
        'LandDay': 100, 
        'LandDayMount': 100, 
        'LandNight':      100,
        'LandNightInv':  100, 
        'LandNightMount': 100, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['HI_T85T11_CM_SEA'], 
        'WaterTwilight': 100 
        } 
    if onlyCirrus:
        title =info+"HighcloudTestt85t11Sea_Cirrus_"+ SchemeName
    else:
        title =info+"HighcloudTestt85t11Sea_All_ "+ SchemeName

    args_test = {'title': title,
            'xlable': 'Latitude',
            'ylable': 'T37-T12 minus dynamic threshold',
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
        'LandDay': 100, 
        'LandDayMount': OFFSETS['HI_T85T11_CM_LAND'], 
        'LandNight':      OFFSETS['HI_T85T11_CM_LAND'],
        'LandNightInv':   OFFSETS['HI_T85T11_CM_LAND'], 
        'LandNightMount': OFFSETS['HI_T85T11_CM_LAND'], 
        'LandTwilight': 100, 
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
            'ylable': 'T37-T12 minus dynamic threshold',
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
                          info="", NEW_THRESHOLD_T37T12=None, onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_textureNight = {}
    OFFSETS_textureNight['t11text_OFFSETS'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['T11TEXT_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': 100 
        }
    OFFSETS_textureNight['t37t12text_OFFSETS'] = {         
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100,  
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': OFFSETS['T37T12TEXT_OFFSET_SEA_NIGHT'], 
        'WaterTwilight': 100 
        }
    if onlyCirrus:
        title =info+"textureNight_Cirrus_"+ SchemeName
    else:
        title =info+"textureNight_All_ "+ SchemeName
    args_test = {'title': title,
                 'ylable': 'T37-T12text minus dynamic threshold',
                 'xlable': 'T11text minus dynamic threshold',
                 'yvector': thr.t37t12text,
                 'xvector': thr.t11text,
                 }
    THRESHOLD2 = OFFSETS_textureNight['t11text_OFFSETS'][SchemeName]
    THRESHOLD1 = OFFSETS_textureNight['t37t12text_OFFSETS'][SchemeName]
    if NEW_THRESHOLD_T37T12 is  not None:
        THRESHOLD1 = NEW_THRESHOLD_T37T12
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_T11TEXT']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_T37T12TEXT']
    TestOk = np.logical_and(thr.t11text>THRESHOLD2, 
                            thr.t37t12text>THRESHOLD1)
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
    plot_test_clear(SchemeName, args_test, args, cloudObj, TestOk, THRESHOLD1,  show=False)
    return TestOk

def brightCloudTest(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_brightCloudTest = {}
    OFFSETS_brightCloudTest['t37t12_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['T37T12_OFFSET_LAND_DAY'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        }     
    OFFSETS_brightCloudTest['r06gain_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R06_GAIN_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_brightCloudTest['r06offset_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R06_OFFSET_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100,  
        'WaterTwilight': 100 
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
    THRESHOLD1 = r06_threshold  + OFFSETS['QUALITY_MARGIN_R06']
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_T37T12']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_R06']
    TestOk =   np.logical_and(thr.t37_t12_minus_threshold>THRESHOLD2,
                              thr.r06>THRESHOLD1)
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                np.mean(THRESHOLD1), THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk


def brightCloudTestNoSunglint3A(SchemeName, caobj, cloudObj, thr, OFFSETS, args,   
                    info="",  onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_brightCloudTestNoSunglint3A = {}
    OFFSETS_brightCloudTestNoSunglint3A['r06r16_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R06_R16_THR_SNOW'],  
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100,  
        'WaterTwilight': 100 
        }     
     
    OFFSETS_brightCloudTestNoSunglint3A['r06gain_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R06_GAIN_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100,  
        'WaterTwilight': 100 
        } 
    OFFSETS_brightCloudTestNoSunglint3A['r06offset_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R06_OFFSET_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
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
    THRESHOLD1 = OFFSETS_brightCloudTestNoSunglint3A['r06r16_OFFSETS'][SchemeName]
    THRESHOLD3 = 1.35 
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_R06']
        THRESHOLD1 += OFFSETS['QUALITY_MARGIN_QR16R06']
        THRESHOLD3 -= OFFSETS['QUALITY_MARGIN_QR16R06']
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
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
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
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 8.0 , 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        }          
    OFFSETS_coldBrightCloudTest37['r06gain_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
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
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
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
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
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
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100,  
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['T11_OFFSET_LAND_DAY'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        }          
    OFFSETS_coldBrightCloudTest['r06gain_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R06_GAIN_LAND'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldBrightCloudTest['r06offset_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100,  
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R06_OFFSET_LAND'],  
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
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
    if  args['USE_MARGINS']:
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_R06']
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR']

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
    OFFSETS_thincoldCirrusTest['t11t12_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['T11T12_OFFSET_LAND_DAY'], 
        'LandNight': 100,
        'LandNightInv': 100,
        'LandNightMount':  100,
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterNight': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_thincoldCirrusTest['t11ts_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastNight': 100, 
        'CoastNightInv': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceNight': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['T11_OFFSET_LAND_DAY_OPAQUE'], 
        'LandNight': 100,
        'LandNightInv': 100,
        'LandNightMount':  100,
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
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
                    THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk

def coldWatercloudTestDay(SchemeName, caobj, cloudObj, thr, OFFSETS, args,  
                               info="", onlyCirrus=False, ExtraCond=None, show=False):
    OFFSETS_coldWatercloudTestDay = {}
    OFFSETS_coldWatercloudTestDay['t11_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 298.15, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        }    
    OFFSETS_coldWatercloudTestDay['t11ts_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['T11_OFFSET_MOUNTAIN_DAY'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldWatercloudTestDay['t11t37_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': 0.0, 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldWatercloudTestDay['r0609_clo_max_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R09_R06_THR_CLOUDY_MAX'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
        } 
    OFFSETS_coldWatercloudTestDay['r0609_clo_min_OFFSETS'] = { 
        'CoastDay': 100, 
        'CoastDayMount': 100, 
        'CoastTwilight': 100, 
        'CoastTwilightInv': 100, 
        'IceDay': 100, 
        'IceTwilight': 100, 
        'LandDay': 100, 
        'LandDayMount': OFFSETS['R09_R06_THR_CLOUDY_MIN'], 
        'LandTwilight': 100, 
        'LandTwilightInv': 100, 
        'LandTwilightMount': 100, 
        'SunglintDay': 100, 
        'SunglintTwilight': 100, 
        'WaterDay': 100, 
        'WaterTwilight': 100 
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
    THRESHOLD1 = OFFSETS_coldWatercloudTestDay['t11ts_OFFSETS'][SchemeName]
    THRESHOLD2 = OFFSETS_coldWatercloudTestDay['t11t37_OFFSETS'][SchemeName]    
    THRESHOLD3 = OFFSETS_coldWatercloudTestDay['t11_OFFSETS'][SchemeName]
    THRESHOLD4 = OFFSETS_coldWatercloudTestDay['r0609_clo_max_OFFSETS'][SchemeName]
    THRESHOLD5 = OFFSETS_coldWatercloudTestDay['r0609_clo_min_OFFSETS'][SchemeName]
    if  args['USE_MARGINS']:
        THRESHOLD1 -= OFFSETS['QUALITY_MARGIN_T11TSUR']
        THRESHOLD2 += OFFSETS['QUALITY_MARGIN_T11T37']
        THRESHOLD3 -= 1.0
        THRESHOLD4 += OFFSETS['QUALITY_MARGIN_QR09R06']
        THRESHOLD5 -= OFFSETS['QUALITY_MARGIN_QR09R06']
    qr0906cond=np.logical_and( qr09r06<THRESHOLD4,
                               qr09r06>THRESHOLD5)
    TestOk = np.logical_and(
        np.logical_and(caobj.avhrr.all_arrays['bt11micron']<THRESHOLD3,
                       thr.qr0906cond),
        np.logical_and( thr.t11_t37_minus_threshold>THRESHOLD2,
                        thr.t11_ts_minus_threshold<THRESHOLD1))
    if ExtraCond is not None:
        TestOk=np.logical_and(ExtraCond, TestOk)
    plot_test_2_lim(SchemeName, args_test, args, cloudObj, TestOk, 
                    THRESHOLD1, THRESHOLD2, onlyCirrus=onlyCirrus, show=show)
    return TestOk
