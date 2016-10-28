import numpy as np
import logging
logger = logging.getLogger(__name__)
from config import PPS_VALIDATION
CHANNEL_MICRON_DESCRIPTIONS = {'11': ["avhrr channel 4 - 11um",
                                      "Avhrr channel channel4.",
                                      "AVHRR ch4",
                                      "AVHRR ch 4",
                                      "channel4",
                                      "AVHRR 4",
                                      "MODIS 31",
                                      "VIIRS M15",
                                      "Avhrr channel channel4."],
                               '12': ["avhrr channel 5 - 12um",
                                      "Avhrr channel channel5.",
                                      "AVHRR ch5",
                                      "AVHRR ch 5",
                                      "channel5",
                                      "AVHRR 5",
                                      "MODIS 32",
                                      "VIIRS M16",
                                      "Avhrr channel channel5."],
                               '06': [ "VIIRS M05",
                                       "AVHRR ch 1", 
                                       "AVHRR ch1",
                                       "AVHRR 1",
                                       "MODIS 1"],
                               '09': [ "VIIRS M07",
                                       "AVHRR ch 2",
                                       "AVHRR ch2",
                                       "AVHRR 2",
                                       "MODIS 2"],
                               '16': [ "VIIRS M10",
                                       "AVHRR ch 3a",
                                       "3a",
                                       "AVHRR ch3a",
                                       "AVHRR 3A",
                                       "MODIS 6"],
                               '37': [ "VIIRS M12",
                                       "AVHRR ch 3b",
                                       "AVHRR ch3b",
                                       "MODIS 20",
                                       "3b",
                                       "AVHRR 3B"],
                               '22': [ "VIIRS M11"],   
                               '13': [ "VIIRS M09",
                                       "MODIS 26"],
                               '86': [ "VIIRS M14",
                                       "MODIS 29"],
                               'modis_3': ['MODIS 3'],
                               'modis_4': ['MODIS 4'],
                               'modis_5': ['MODIS 5'],
                               'modis_7': ['MODIS 7'],
                               'modis_8': ['MODIS 8'],
                               'modis_9': ['MODIS 9'],
                               'modis_10': ['MODIS 10'],
                               'modis_11': ['MODIS 11'],
                               'modis_12': ['MODIS 12'],
                               'modis_13lo': ['MODIS 13lo'],
                               'modis_13hi': ['MODIS 13hi'],
                               'modis_14lo': ['MODIS 14lo'],
                               'modis_14hi': ['MODIS 14hi'],
                               'modis_15': ['MODIS 15'],
                               'modis_16': ['MODIS 16'],
                               'modis_17': ['MODIS 17'],
                               'modis_18': ['MODIS 18'],
                               'modis_19': ['MODIS 19'],
                               'modis_21': ['MODIS 21'],
                               'modis_22': ['MODIS 22'],
                               'modis_23': ['MODIS 23'],
                               'modis_24': ['MODIS 24'],
                               'modis_25': ['MODIS 25'],
                               'modis_27': ['MODIS 27'],
                               'modis_28': ['MODIS 28'],
                               'modis_30': ['MODIS 30'],
                               'modis_33': ['MODIS 33'],
                               'modis_34': ['MODIS 34'],
                               'modis_35': ['MODIS 35'],
                               'modis_36': ['MODIS 36']                     
                               }
CHANNEL_MICRON_AVHRR_PPS = {'11': 3,
                            '12': 4,
                            '06': 0,
                            '09': 1,
                            '37': 2,
                            '86': -1,               
                            '16': 5,
                            '22': -1,
                            '13': -1}  
CURRENTLY_UNUSED_MODIS_CHANNELS = ['modis_3',
                                   'modis_4',
                                   'modis_5',
                                   'modis_7',
                                   'modis_8',
                                   'modis_9',
                                   'modis_10',
                                   'modis_11',
                                   'modis_12',
                                   'modis_13lo',
                                   'modis_13hi',
                                   'modis_14lo',
                                   'modis_14hi',
                                   'modis_15',
                                   'modis_16',
                                   'modis_17',
                                   'modis_18',
                                   'modis_19',
                                   'modis_21',
                                   'modis_22',
                                   'modis_23',
                                   'modis_24',
                                   'modis_25',
                                   'modis_27',
                                   'modis_28',
                                   'modis_30',
                                   'modis_33',
                                   'modis_34',
                                   'modis_35',
                                   'modis_36']

def get_channel_data_from_object(dataObj, chn_des, matched, nodata=-9):
    """Get the AVHRR/VIIRS channel data on the track

    matched: dict of matched indices (row, col)

    """
    try:
        channels = dataObj.channels
    except:
        channels = dataObj.channel
    
    numOfChannels = len(channels)
    chnum=-1
    for ich in range(numOfChannels):
        if channels[ich].des in CHANNEL_MICRON_DESCRIPTIONS[chn_des]:
            chnum = ich
    if chnum ==-1:
        #chnum = CHANNEL_MICRON_AVHRR_PPS[chn_des]
        if chnum ==-1:
            return None
        logger.warning(  "Using pps channel numbers to find "
              "corresponding avhrr channel")
    temp = [channels[chnum].data[matched['row'][idx], 
                                 matched['col'][idx]]
            for idx in range(matched['row'].shape[0])] 
    chdata_on_track = [channels[chnum].data[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])]
    return np.array(chdata_on_track)



def insert_nwp_segments_data(nwp_segments, row_matched, col_matched, obt):
        npix = row_matched.shape[0]
        """
        #obt.avhrr.segment_nwgeoheight
        obt.avhrr.segment_nwp_moist
        obt.avhrr.segment_nwp_pressure
        obt.avhrr.segment_nwp_temp
        obt.avhrr.segment_surfaceLandTemp
        obt.avhrr.segment_surfaceSeaTemp
        obt.avhrr.segment_surfaceGeoHeight
        obt.avhrr.segment_surfaceMoist
        obt.avhrr.segment_surfacePressure
        obt.avhrr.segment_fractionOfLand
        obt.avhrr.segment_meanElevation
        obt.avhrr.segment_ptro
        obt.avhrr.segment_ttro
        #obt.avhrr.segment_t850
        obt.avhrr.segment_tb11clfree_sea
        obt.avhrr.segment_tb12clfree_sea
        obt.avhrr.segment_tb11clfree_land
        obt.avhrr.segment_tb12clfree_land
        obt.avhrr.segment_tb11cloudy_surface
        obt.avhrr.segment_tb12cloudy_surface  
        """
        def get_segment_row_col_idx(nwp_segments, row_matched, col_matched):
            segment_colidx = nwp_segments['colidx']
            segment_rowidx = nwp_segments['rowidx']
            seg_row = np.zeros(np.size(row_matched)) -9
            seg_col = np.zeros(np.size(col_matched)) -9
            for s_col in xrange(nwp_segments['norows']):
                for s_row in xrange(nwp_segments['nocols']):
                    within_segment = np.logical_and(
                        np.logical_and(row_matched>=segment_rowidx[s_row,s_col]-nwp_segments['segSizeX']/2,
                                       row_matched<segment_rowidx[s_row,s_col]+nwp_segments['segSizeX']/2),
                        np.logical_and(col_matched>=segment_colidx[s_row,s_col]-nwp_segments['segSizeY']/2,
                                       col_matched<segment_colidx[s_row,s_col]+nwp_segments['segSizeY']/2))
                    seg_row[within_segment] = s_row
                    seg_col[within_segment] = s_col
            return  seg_row, seg_col   
        seg_row, seg_col = get_segment_row_col_idx(nwp_segments, row_matched, col_matched)
        for data_set in ['surfaceLandTemp',
                         'surfaceSeaTemp',
                         'surfaceGeoHeight',
                         'surfaceMoist',
                         'surfacePressure',
                         'fractionOfLand',
                         'meanElevation',
                         'ptro',
                         'ttro',
                         't850',
                         'tb11clfree_sea',
                         'tb12clfree_sea',
                         'tb11clfree_land',
                         'tb12clfree_land']:
                         #'tb11cloudy_surface',
                         #'tb12cloudy_surface ',
            setattr(obt.avhrr,'segment_nwp_' + data_set, 
                    np.array([nwp_segments[data_set][seg_row[idx], seg_col[idx]]
                              for idx in range(npix)]))
        #obt.avhrr.segment_t850 = np.array([nwp_segments['t850'][seg_row[idx], seg_col[idx]]
        #                                   for idx in range(npix)])

        for data_set in ['moist', 'pressure', 'geoheight', 'temp']:
            setattr(obt.avhrr,'segment_nwp_' + data_set, 
                              np.array([nwp_segments[data_set][seg_row[idx], seg_col[idx]]
                                        for idx in range(npix)]))
        #obt.avhrr.segment_nwp_geoheight = np.array([nwp_segments['geoheight'][seg_row[idx], seg_col[idx]]
        #                                            for idx in range(npix)])

        return obt

#---------------------------------------------------------------------------
def avhrr_track_from_matched(obt, GeoObj, dataObj, AngObj, 
                             nwp_obj, ctth, ctype, cma, 
                             row_matched, col_matched, 
                             avhrrLwp=None, avhrrCph=None,
                             nwp_segments=None):
    value_track = []
    row_col = {'row': row_matched, 'col': col_matched} 
    npix = row_matched.shape[0]

    #Find lat/lon and cloudtype
    value_track = [GeoObj.latitude[row_matched[idx], col_matched[idx]] 
                     for idx in range(npix)]
    obt.avhrr.latitude = np.array(value_track)
    value_track = [GeoObj.longitude[row_matched[idx], col_matched[idx]]
                     for idx in range(npix)]
    obt.avhrr.longitude = np.array(value_track)
    value_track = [ctype.cloudtype[row_matched[idx], col_matched[idx]]
                 for idx in range(npix)]
    obt.avhrr.cloudtype = np.array(value_track)
    value_track = [cma.cma_ext[row_matched[idx], col_matched[idx]]
                   for idx in range(npix)]
    obt.avhrr.cloudmask = np.array(value_track)
    #cloud-type and ctth flags
    if hasattr(ctype, 'ct_quality') and PPS_VALIDATION:
        value_track = [ctype.ct_quality[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_quality = np.array(value_track)
    if hasattr(ctype, 'ct_conditions') and ctype.ct_conditions is not None:
        value_track = [ctype.ct_conditions[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_conditions = np.array(value_track)
    if hasattr(ctype, 'ct_statusflag') and PPS_VALIDATION:
        value_track = [ctype.ct_statusflag[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_status = np.array(value_track)
    if hasattr(ctth, 'ctth_statusflag') and PPS_VALIDATION:
        value_track = [ctth.ctth_statusflag[row_matched[idx], col_matched[idx]]
                                       for idx in range(npix)]
        obt.avhrr.ctth_status = np.array(value_track)
    if hasattr(ctype, 'qualityflag') and PPS_VALIDATION:
        value_track = [ctype.qualityflag[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_qflag = np.array(value_track)
    if  ctype.phaseflag != None and PPS_VALIDATION:
        value_track = [ctype.phaseflag[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_pflag = np.array(value_track)
    for nwp_info in ["surftemp", "t500", "t700", "t850", "t950", "ttro", "ciwv",
                     "t900", "t1000", "t800", "t250", "t2m", "ptro", "psur"]:
        if hasattr(nwp_obj, nwp_info):
            data = getattr(nwp_obj, nwp_info)
            if np.size(data)>1:
                value_track = [data[row_matched[idx], col_matched[idx]]
                               for idx in range(npix)]
                setattr(obt.avhrr, nwp_info, np.array(value_track))
        else:
            print "missing", nwp_info
    from pps_prototyping_util import (get_t11t12_texture_data_from_object,
                                      get_coldest_values,get_darkest_values,
                                      get_warmest_values)
    if dataObj != None:
        nwp_obj = get_t11t12_texture_data_from_object(dataObj, nwp_obj, '11','12', 
                                                      'text_t11t12_square') 
        nwp_obj = get_t11t12_texture_data_from_object(dataObj, nwp_obj, '37','12',
                                                      'text_t37t12_square')
    for texture in ["text_r06", "text_t11", "text_t37", "text_t37t12", 
                    "text_t37t12_square", "text_t11t12_square", "text_t11t12"]:
        if hasattr(nwp_obj, texture):
            data = getattr(nwp_obj, texture)
            if np.size(data)>1:
                value_track = [data[row_matched[idx], col_matched[idx]]
                               for idx in range(npix)]
                setattr(obt.avhrr, texture, np.array(value_track))
    if dataObj != None:
        neighbour_obj =get_warmest_values(dataObj, row_col)
        for key in ["warmest_t11", "warmest_t12", "warmest_t37",
                    "warmest_r06", "warmest_r09", "warmest_r16"]:
            setattr(obt.avhrr, key, np.array(getattr(neighbour_obj,key)))
        neighbour_obj =get_darkest_values(dataObj, row_col)
        for key in ["darkest_t11", "darkest_t12", "darkest_t37",
                    "darkest_r06", "darkest_r09", "darkest_r16"]:
            setattr(obt.avhrr, key, np.array(getattr(neighbour_obj,key)))
        neighbour_obj =get_coldest_values(dataObj, row_col)
        for key in ["coldest_t11", "coldest_t12", "coldest_t37",
                    "coldest_r06", "coldest_r09", "coldest_r16"]:
            setattr(obt.avhrr, key, np.array(getattr(neighbour_obj,key)))
    #Thresholds:    
    for thr in ["thr_t11ts_inv",
                "thr_t85t11_inv",
                "thr_t11t37_inv",
                "thr_t37t12_inv",
                "thr_t11t12_inv",
                "thr_t85t11",
                "thr_t11ts",
                "thr_t11t37",
                "thr_t37t12",
                "thr_t11t12",
                "thr_r06",
                "thr_r09"]:
        if hasattr(nwp_obj, thr):
            data = getattr(nwp_obj, thr)
            if np.size(data)>1:
                value_track = [data[row_matched[idx], col_matched[idx]]
                               for idx in range(npix)]
                setattr(obt.avhrr, thr, np.array(value_track))
    for emis in ["emis1", "emis6", "emis8", "emis9"]:
        if hasattr(nwp_obj, emis):
            data = getattr(nwp_obj, emis)
            if np.size(data)>1:
                value_track = [data[row_matched[idx], col_matched[idx]]
                               for idx in range(npix)]
                setattr(obt.avhrr, emis, np.array(value_track))
    if dataObj != None:
        obt.avhrr.r06micron = get_channel_data_from_object(dataObj, '06', row_col)
        # r09   
        obt.avhrr.r09micron = get_channel_data_from_object(dataObj, '09', row_col)
        # bt37   
        obt.avhrr.bt37micron = get_channel_data_from_object(dataObj, '37', row_col)
        # b11
        obt.avhrr.bt11micron = get_channel_data_from_object(dataObj, '11', row_col)
        # b12
        obt.avhrr.bt12micron = get_channel_data_from_object(dataObj, '12', row_col)
        # b86
        obt.avhrr.bt86micron = get_channel_data_from_object(dataObj, '86', row_col)
        # b16
        obt.avhrr.r16micron = get_channel_data_from_object(dataObj, '16', row_col)
        # b22
        obt.avhrr.r22micron = get_channel_data_from_object(dataObj, '22', row_col)
        #b13
        obt.avhrr.r13micron = get_channel_data_from_object(dataObj, '13', row_col)
        for modis_channel in CURRENTLY_UNUSED_MODIS_CHANNELS:
            modis_track = get_channel_data_from_object(dataObj, 
                                                       modis_channel, row_col)
            setattr(obt.avhrr, modis_channel, modis_track)
    #Angles, scale with gain and intercept when reading
    obt.avhrr.satz = [AngObj.satz.data[row_matched[idx], col_matched[idx]] 
                      for idx in range(npix)]
    obt.avhrr.sunz = [AngObj.sunz.data[row_matched[idx], col_matched[idx]] 
                      for idx in range(npix)]
    if AngObj.azidiff is not None:
        obt.avhrr.azidiff = [AngObj.azidiff.data[row_matched[idx], col_matched[idx]] 
                             for idx in range(npix)]
    if ctth == None:
        logger.info("Not extracting ctth")
    else:
        logger.info("Extracting ctth along track ")
        #scale with gain and intercept when reading!
        obt.avhrr.ctth_height = [ctth.height[row_matched[idx], col_matched[idx]]
                                 for idx in range(npix)]
        obt.avhrr.ctth_temperature = [ctth.temperature[row_matched[idx], col_matched[idx]]
                                      for idx in range(npix)]
        obt.avhrr.ctth_pressure = [ctth.pressure[row_matched[idx], col_matched[idx]]
                                   for idx in range(npix)]

        if (PPS_VALIDATION and hasattr(ctth, 'processingflag')):
            is_opaque = np.bitwise_and(np.right_shift(ctth.processingflag, 2), 1)
            value_track = [is_opaque[row_matched[idx], col_matched[idx]]
                           for idx in range(npix)]
            obt.avhrr.ctth_opaque = np.array(value_track)
    #NWP on ctth resolution        
    if nwp_segments != None:
        obt = insert_nwp_segments_data(nwp_segments, row_matched, col_matched, obt)
  
    #: TODO Do not use fix nodata-values but instead something.no_data
    #CPP-products
    if avhrrLwp != None:
        if PPS_FORMAT_2012_OR_EARLIER:
            nodata_temp = -1
        else:
            nodata_temp = 65535
        lwp_temp = [avhrrLwp[row_matched[idx], col_matched[idx]]
                    for idx in range(npix)]
        obt.avhrr.lwp = np.where(np.equal(lwp_temp, nodata_temp), -9, lwp_temp)
    if avhrrCph != None:
        if PPS_FORMAT_2012_OR_EARLIER:
            nodata_temp = -1
        else:
            nodata_temp = 255
        cph_temp = [avhrrCph[row_matched[idx], col_matched[idx]]
                    for idx in range(npix)]
        obt.avhrr.cph = np.where(np.equal(cph_temp, nodata_temp), -9, cph_temp)

    return obt


    
