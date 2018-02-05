import numpy as np
import logging
logger = logging.getLogger(__name__)
from config import PPS_VALIDATION, SAVE_NEIGHBOUR_INFO, CTTH_TYPES
CHANNEL_MICRON_DESCRIPTIONS = {'11': ["avhrr channel 4 - 11um",
                                      "Avhrr channel channel4.",
                                      "AVHRR ch4",
                                      "AVHRR ch 4",
                                      "AVHRR ch_4",
                                      "SEVIRI IR_108",
                                      "channel4",
                                      "AVHRR 4",
                                      "MODIS 31",
                                      "VIIRS M15",
                                      "Avhrr channel channel4."],
                               '12': ["avhrr channel 5 - 12um",
                                      "Avhrr channel channel5.",
                                      "AVHRR ch5",
                                      "AVHRR ch 5",
                                      "AVHRR ch_5",
                                      "SEVIRI IR_120",
                                      "channel5",
                                      "AVHRR 5",
                                      "MODIS 32",
                                      "VIIRS M16",
                                      "Avhrr channel channel5."],
                               '06': [ "VIIRS M05",
                                       "AVHRR ch 1", 
                                       "SEVIRI VIS006",
                                       "AVHRR ch1",
                                       "AVHRR ch_1",
                                       "AVHRR 1",
                                       "MODIS 1"],
                               '09': [ "VIIRS M07",
                                       "AVHRR ch 2",
                                       "AVHRR ch2",
                                       "AVHRR ch_2",
                                       "AVHRR 2",
                                       "SEVIRI VIS008",
                                       "MODIS 2"],
                               '16': [ "VIIRS M10",
                                       "AVHRR ch 3a",
                                       "3a",
                                       "SEVIRI IR_016",
                                       "AVHRR ch3a",
                                       "AVHRR ch_3a",
                                       "AVHRR 3A",
                                       "MODIS 6"],
                               '37': [ "VIIRS M12",
                                       "AVHRR ch 3b",
                                       "SEVIRI IR_039",
                                       "AVHRR ch_3b",
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
    chnum = -1
    for ich in range(numOfChannels):
        if channels[ich].des in CHANNEL_MICRON_DESCRIPTIONS[chn_des]:
            chnum = ich
    if chnum == -1:
        #chnum = CHANNEL_MICRON_AVHRR_PPS[chn_des]
        if chnum == -1:
            logger.debug("Did not find pps channel number for channel "
                         "{:s}".format(chn_des))
            return None, "" 
    if matched is None:
        return channels[chnum].data, ""
    
    chdata_on_track = np.array([channels[chnum].data[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])])
    extra_info = ""
    if channels[chnum].SZA_corr_done:
        extra_info = "_sza_correction_done"
    return np.array(chdata_on_track), extra_info


def _interpolate_height_and_temperature_from_pressure(imagerObj,
                                                      level):
    """ Function to find height att pressure level (level)
    from segment_nwp, pressure and height vectors.
    High means high in pressure. The level closest to ground i hi, and lo is at lower 
    pressure further up in atmosphere.
    """
    values_h =  imagerObj.segment_nwp_geoheight
    pressure_v=  imagerObj.segment_nwp_pressure
    surface_h = imagerObj.segment_nwp_surfaceGeoHeight
    psur = imagerObj.segment_nwp_surfacePressure
    nlev = pressure_v.shape[1]
    npix = pressure_v.shape[0]
    k = np.arange(npix)
    higher_index = np.array([nlev -1 - np.searchsorted(pressure_v[ind,:], level, side='right',
                                                       sorter=xrange(nlev -1, -1, -1)) 
                             for ind in xrange(npix)])
    higher_index[higher_index >= (nlev - 1)] = nlev - 2
    lower_index = higher_index + 1
    # update "lo" where level is between surface and first level in array
    below_level_1 = level > pressure_v[:,0]
    lower_index[below_level_1] = 0
    # get pressure and height for layer below and above level
    hi = pressure_v[k,higher_index]
    lo = pressure_v[k,lower_index]
    height_hi_ = values_h[k,higher_index]*1.0
    height_lo_ = values_h[k,lower_index]*1.0
    # update "hi" where level is between surface and first level in array
    hi[below_level_1] = psur[below_level_1]
    height_hi_[below_level_1] = surface_h[below_level_1]
    # log pressures
    hi = np.log(hi)
    lo = np.log(lo)
    level = np.log(level)
    # interpolate
    out_h = height_hi_ - (hi - level) * (height_hi_ - height_lo_) / (hi - lo)
    return out_h 


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
                        np.logical_and(
                            row_matched >= (segment_rowidx[s_row,s_col] 
                                            - nwp_segments['segSizeX']/2),
                            row_matched < (segment_rowidx[s_row,s_col] 
                                           + nwp_segments['segSizeX']/2)),
                        np.logical_and(
                            col_matched >= (segment_colidx[s_row,s_col]
                                            - nwp_segments['segSizeY']/2),
                            col_matched < (segment_colidx[s_row,s_col]
                                           + nwp_segments['segSizeY']/2)))
                    seg_row[within_segment] = s_row
                    seg_col[within_segment] = s_col
            return  seg_row.astype(np.int16), seg_col.astype(np.int16)   
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
                         'tb12clfree_land',
                         'tb11lowcloud_sea',
                         'tb12lowcloud_sea',
                         'tb11lowcloud_land',
                         'tb12lowcloud_land']:
            if data_set in nwp_segments.keys():
                #'tb11cloudy_surface',
                #'tb12cloudy_surface ',
                setattr(obt.avhrr,'segment_nwp_' + data_set, 
                        np.array([nwp_segments[data_set][seg_row[idx], seg_col[idx]]
                                  for idx in range(npix)]))
            elif 'clfree' in data_set or 'lowcloud' in data_set:
                #these are nor always present
                pass

                
        for data_set in ['moist', 'pressure', 'geoheight', 'temp']:
            setattr(obt.avhrr,'segment_nwp_' + data_set, 
                              np.array([nwp_segments[data_set][seg_row[idx], seg_col[idx]]
                                        for idx in range(npix)]))
        #Remove nodata and not used upper part of atmosphere   
        N = obt.avhrr.segment_nwp_pressure.shape[1]
        pressure_n_to_keep = np.sum(np.max(obt.avhrr.segment_nwp_pressure,axis=0)>50)
        logger.debug("Not saving upper %d levels of 3-D nwp from segment file"%(N-pressure_n_to_keep))
        logger.debug("Keeping %d lower levels of 3-D nwp from segment file"%(pressure_n_to_keep))
        for data_set in ['segment_nwp_moist', 'segment_nwp_pressure', 
                          'segment_nwp_geoheight', 'segment_nwp_temp']:
            data = getattr(obt.avhrr, data_set) 
            setattr(obt.avhrr, data_set, data[:,0:pressure_n_to_keep])

        

        #obt.avhrr.segment_nwp_geoheight = np.array([nwp_segments['geoheight'][seg_row[idx], seg_col[idx]]
        #                                            for idx in range(npix)])
        #Extract h440: hight at 440 hPa, and h680
        data = _interpolate_height_and_temperature_from_pressure(obt.avhrr, 440)
        setattr(obt.avhrr, 'segment_nwp_h440', data)
        data = _interpolate_height_and_temperature_from_pressure(obt.avhrr, 680)
        setattr(obt.avhrr, 'segment_nwp_h680', data)
        
        return obt

#---------------------------------------------------------------------------
def avhrr_track_from_matched(obt, GeoObj, dataObj, AngObj, 
                             nwp_obj, ctth, ctype, cma, 
                             row_matched, col_matched, 
                             cpp=None,
                             nwp_segments=None):
    value_track = []
    row_col = {'row': row_matched, 'col': col_matched} 
    npix = row_matched.shape[0]


    value_track = [GeoObj.latitude[row_matched[idx], col_matched[idx]] 
                   for idx in range(npix)]
    obt.avhrr.latitude = np.array(value_track)
    value_track = [GeoObj.longitude[row_matched[idx], col_matched[idx]]
                     for idx in range(npix)]
    obt.avhrr.longitude = np.array(value_track)
    if ctype is not None:
        value_track = [ctype.cloudtype[row_matched[idx], col_matched[idx]]
                       for idx in range(npix)]
        obt.avhrr.cloudtype = np.array(value_track)
    if cma is not None:
        value_track = [cma.cma_ext[row_matched[idx], col_matched[idx]]
                       for idx in range(npix)]
        obt.avhrr.cloudmask = np.array(value_track)
    #cma flags or cma prob
    if hasattr(cma, 'cma_prob') and cma.cma_prob is not None:
        obt.avhrr.cma_prob = np.array(
            [cma.cma_prob[row_matched[idx], col_matched[idx]]
             for idx in range(npix)])
    if hasattr(cma, 'cma_aerosolflag') and cma.cma_aerosolflag is not None:
        obt.avhrr.cma_aerosolflag = np.array(
            [cma.cma_aerosolflag[row_matched[idx], col_matched[idx]]
             for idx in range(npix)])
    if hasattr(cma, 'cma_testlist0') and cma.cma_testlist0 is not None:
        obt.avhrr.cma_testlist0 = np.array(
            [cma.cma_testlist0[row_matched[idx], col_matched[idx]]
            for idx in range(npix)])
        obt.avhrr.cma_testlist1 = np.array(
            [cma.cma_testlist1[row_matched[idx], col_matched[idx]]
             for idx in range(npix)])
        obt.avhrr.cma_testlist2 = np.array(
            [cma.cma_testlist2[row_matched[idx], col_matched[idx]]
             for idx in range(npix)])
        obt.avhrr.cma_testlist3 = np.array(
            [cma.cma_testlist3[row_matched[idx], col_matched[idx]]
             for idx in range(npix)])
        obt.avhrr.cma_testlist4 = np.array(
            [cma.cma_testlist4[row_matched[idx], col_matched[idx]]
            for idx in range(npix)])
        obt.avhrr.cma_testlist5 = np.array(
            [cma.cma_testlist5[row_matched[idx], col_matched[idx]]
             for idx in range(npix)])
    #cloud-type and ctth flags
    if hasattr(ctype, 'ct_quality') and PPS_VALIDATION and ctype.ct_quality is not None:
        value_track = [ctype.ct_quality[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_quality = np.array(value_track)
    if hasattr(ctype, 'ct_conditions') and ctype.ct_conditions is not None:
        value_track = [ctype.ct_conditions[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_conditions = np.array(value_track)
    if hasattr(ctype, 'ct_statusflag') and PPS_VALIDATION and ctype.ct_statusflag is not None:
        value_track = [ctype.ct_statusflag[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_status = np.array(value_track)
    if hasattr(ctth, 'ctth_statusflag') and PPS_VALIDATION and ctth.ctth_statusflag is not None:
        value_track = [ctth.ctth_statusflag[row_matched[idx], col_matched[idx]]
                                       for idx in range(npix)]
        obt.avhrr.ctth_status = np.array(value_track)
    if hasattr(ctype, 'qualityflag') and PPS_VALIDATION and ctype.ct_qualityflag is not None:
        value_track = [ctype.qualityflag[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_qflag = np.array(value_track)
    if  hasattr(ctype, 'phaseflag') and PPS_VALIDATION and ctype.phaseflag is not None:
        value_track = [ctype.phaseflag[row_matched[idx], col_matched[idx]]
                             for idx in range(npix)]
        obt.avhrr.cloudtype_pflag = np.array(value_track)
    for nwp_info in ["surftemp", "t500", "t700", "t850", "t950", "ttro", "ciwv",
                     "t900", "t1000", "t800", "t250", "t2m", "ptro", "psur", 
                     "snowa", "snowd", "seaice", "landuse", "fractionofland", "elevation",
                     "r37_sza_correction_done"]:
        if hasattr(nwp_obj, nwp_info):
            data = getattr(nwp_obj, nwp_info)
            if np.size(data)>1:
                value_track = [data[row_matched[idx], col_matched[idx]]
                               for idx in range(npix)]
                setattr(obt.avhrr, nwp_info, np.array(value_track))
        else:
            logger.debug("missing {:s}".format(nwp_info))
    if len(CTTH_TYPES)>1:        
        for ctth_type in CTTH_TYPES[1:]:
            ctth_obj = getattr(nwp_obj,ctth_type)
            for data_set in ["pressure", "temperature", "height"]:
                data = getattr(ctth_obj, data_set)
                value_track = [data[row_matched[idx], col_matched[idx]]
                               for idx in range(npix)]
                name = "%s_%s"%(ctth_type.lower(),data_set)
                setattr(obt.avhrr, name, np.array(value_track))
            
        
    from pps_prototyping_util import (get_t11t12_texture_data_from_object,
                                      get_coldest_values,get_darkest_values,
                                      get_warmest_values)
    if dataObj is not None:
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
    if dataObj is not None and SAVE_NEIGHBOUR_INFO:
        neighbour_obj = get_warmest_values(dataObj, row_col)
        for key in ["warmest_r06", "warmest_r09", "warmest_r16"]:
            setattr(obt.avhrr, key + neighbour_obj.extra_info_sza_corr, 
                    getattr(neighbour_obj,key))
        for key in ["warmest_t11", "warmest_t12", "warmest_t37"]:
            setattr(obt.avhrr, key, 
                    getattr(neighbour_obj,key))
        neighbour_obj = get_darkest_values(dataObj, row_col)
        for key in ["darkest_r06", "darkest_r09", "darkest_r16"]:
            setattr(obt.avhrr, key + neighbour_obj.extra_info_sza_corr, 
                    getattr(neighbour_obj,key))
        for key in ["darkest_t11", "darkest_t12", "darkest_t37"]:
            setattr(obt.avhrr, key, 
                    getattr(neighbour_obj,key))
        neighbour_obj = get_coldest_values(dataObj, row_col)
        for key in ["coldest_r06", "coldest_r09", "coldest_r16"]:
            setattr(obt.avhrr, key + neighbour_obj.extra_info_sza_corr, 
                    getattr(neighbour_obj,key))
        for key in ["coldest_t11", "coldest_t12", "coldest_t37"]:
            setattr(obt.avhrr, key, 
                    getattr(neighbour_obj,key))
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
    if dataObj is not None:
        temp_data, info = get_channel_data_from_object(dataObj, '06', row_col)
        setattr(obt.avhrr, "r06micron" + info, temp_data) 
        # r09   
        temp_data, info =  get_channel_data_from_object(dataObj, '09', row_col)
        setattr(obt.avhrr, "r09micron" + info, temp_data)
        # bt37   
        temp_data, info =  get_channel_data_from_object(dataObj, '37', row_col)
        setattr(obt.avhrr, "bt37micron" + info, temp_data)
        # b11
        temp_data, info =  get_channel_data_from_object(dataObj, '11', row_col)
        setattr(obt.avhrr, "bt11micron" + info, temp_data)
        # b12
        temp_data, info =  get_channel_data_from_object(dataObj, '12', row_col)
        setattr(obt.avhrr, "bt12micron" + info, temp_data)
        # b86
        temp_data, info =  get_channel_data_from_object(dataObj, '86', row_col)
        setattr(obt.avhrr, "bt86micron" + info, temp_data)
        # b16
        temp_data, info =  get_channel_data_from_object(dataObj, '16', row_col)
        setattr(obt.avhrr, "r16micron" + info, temp_data)
        # b22
        temp_data, info =  get_channel_data_from_object(dataObj, '22', row_col)
        setattr(obt.avhrr, "r22micron" + info, temp_data)
        #b13
        temp_data, info =  get_channel_data_from_object(dataObj, '13', row_col)
        setattr(obt.avhrr, "r13micron" + info, temp_data)
        for modis_channel in CURRENTLY_UNUSED_MODIS_CHANNELS:
            modis_track, info = get_channel_data_from_object(dataObj, 
                                                       modis_channel, row_col)
            setattr(obt.avhrr, modis_channel + info, modis_track)
    #Angles, scale with gain and intercept when reading
    obt.avhrr.satz = [AngObj.satz.data[row_matched[idx], col_matched[idx]] 
                      for idx in range(npix)]
    obt.avhrr.sunz = [AngObj.sunz.data[row_matched[idx], col_matched[idx]] 
                      for idx in range(npix)]
    if AngObj.azidiff is not None:
        obt.avhrr.azidiff = [AngObj.azidiff.data[row_matched[idx], col_matched[idx]] 
                             for idx in range(npix)]
    if ctth is None:
        logger.info("Not extracting ctth")
    else:
        logger.debug("Extracting ctth along track ")
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
    if nwp_segments is not None:
        obt = insert_nwp_segments_data(nwp_segments, row_matched, col_matched, obt)
    if cpp is None:    
        logger.debug("Not extracting cpp")
    else:
        logger.debug("Extracting cpp along track ")
        for data_set_name in cpp.__dict__.keys():
            data = getattr(cpp, data_set_name)
            if data is not None:
                extracted_data = [data[row_matched[idx], col_matched[idx]] 
                                  for idx in range(npix)]
                setattr(obt.avhrr, data_set_name, extracted_data)

    return obt


    
