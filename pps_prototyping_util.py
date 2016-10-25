import numpy as np
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
                               'modis 12': ['MODIS 12'],
                               'modis_13lo': ['MODIS 13lo'],
                               'modis_13hi': ['MODIS 13hi'],
                               'modis_14lo': ['MODIS 14lo'],
                               'modis_14hi': ['MODIS 14hi'],
                               'modis_15': ['MODIS 15'],
                               'modis 16': ['MODIS 16'],
                               'modis_17': ['MODIS 17'],
                               'modis_18': ['MODIS 18'],
                               'modis_19': ['MODIS 19'],
                               'modis 21': ['MODIS 21'],
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

def get_t11t12_texture_data_from_object(dataObj, nwp_obj, ch11, ch12, text_name):
    #https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    t11 = get_channel_data_from_object_all(dataObj, ch11, nodata=-9)
    t12 = get_channel_data_from_object_all(dataObj, ch12, nodata=-9)    
    t11t12 = 1.0*np.array(t11-t12)
    K=np.median(t11t12) #K is trick to get better accurracy maybe not needed as differences are often small
    t11t12 = (t11t12-K)
    from scipy.ndimage import uniform_filter
    mean = uniform_filter(t11t12, size=(5,5), mode='mirror')
    mean_of_squared = uniform_filter(t11t12**2, size=(5,5), mode='mirror')    
    t11t12_texture = mean_of_squared - mean**2 
    setattr(nwp_obj, text_name, t11t12_texture)
    return nwp_obj


class NeighbourObj(object):
    def __init__(self):
        self.warmest_t11 = None
        self.warmest_t12 = None
        self.warmest_t37 = None
        self.warmest_r16 = None
        self.warmest_r09 = None
        self.warmest_r06 = None
        self.darkest_t11 = None
        self.darkest_t12 = None
        self.darkest_t37 = None
        self.darkest_r16 = None
        self.darkest_r09 = None
        self.darkest_r06 = None
        self.coldest_t11 = None
        self.coldest_t12 = None
        self.coldest_t37 = None
        self.coldest_r16 = None
        self.coldest_r09 = None
        self.coldest_r06 = None
def get_warmest_values(dataObj, matched):
    nobj = NeighbourObj()
    t11 = get_channel_data_from_object_all(dataObj, '11', nodata=-9)
    #t11 = np.random.rand(10,10)*10+270
    row, col = np.indices(t11.shape)
    from scipy.ndimage.filters import generic_filter
    flat_index = generic_filter(t11,
                                function=np.argmax,
                                size=5,
                                mode='constant',
                                cval=-9999999999999)
    flat_index = np.array(flat_index, dtype=np.int)
    delta_row, delta_col = np.unravel_index(flat_index, (5,5))
    delta_row = delta_row - 2
    delta_col = delta_col - 2
    
    new_row = row+delta_row
    new_col = col+delta_col
    new_row_matched = np.array([new_row[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])]) 
    new_col_matched = np.array([new_col[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])]) 
    new_row_col = {'row': new_row_matched, 'col': new_col_matched}
    from calipso import get_channel_data_from_object
    nobj.warmest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)
    nobj.warmest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)
    nobj.warmest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)
    nobj.warmest_r06=get_channel_data_from_object(dataObj, '06', new_row_col)
    nobj.warmest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)
    nobj.warmest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)
    return nobj

def get_coldest_values(dataObj, matched):
    nobj = NeighbourObj()
    t11 = get_channel_data_from_object_all(dataObj, '11', nodata=-9)
    #t11 = np.random.rand(10,10)*10+270
    row, col = np.indices(t11.shape)
    from scipy.ndimage.filters import generic_filter
    flat_index = generic_filter(t11,
                                function=np.argmin,
                                size=5,
                                mode='constant',
                                cval=9999999999999)
    flat_index = np.array(flat_index, dtype=np.int)
    delta_row, delta_col = np.unravel_index(flat_index, (5,5))
    delta_row = delta_row - 2
    delta_col = delta_col - 2
    
    new_row = row+delta_row
    new_col = col+delta_col
    new_row_matched = np.array([new_row[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])]) 
    new_col_matched = np.array([new_col[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])]) 
    new_row_col = {'row': new_row_matched, 'col': new_col_matched}
    from calipso import get_channel_data_from_object
    nobj.coldest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)
    nobj.coldest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)
    nobj.coldest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)
    nobj.coldest_r06=get_channel_data_from_object(dataObj, '06', new_row_col)
    nobj.coldest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)
    nobj.coldest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)
    return nobj

def get_darkest_values(dataObj, matched):
    nobj = NeighbourObj()
    r09 = get_channel_data_from_object_all(dataObj, '09', nodata=-9)
    r06 = get_channel_data_from_object_all(dataObj, '06', nodata=-9)
    darkest = np.where(r09>r06,r06,r09)
    row, col = np.indices(r06.shape)
    from scipy.ndimage.filters import generic_filter
    flat_index = generic_filter(darkest,
                                function=np.argmin,
                                size=5,
                                mode='constant',
                                cval=9999999999999)
    flat_index = np.array(flat_index, dtype=np.int)
    delta_row, delta_col = np.unravel_index(flat_index, (5,5))
    delta_row = delta_row - 2
    delta_col = delta_col - 2
    
    new_row = row+delta_row
    new_col = col+delta_col
    new_row_matched = np.array([new_row[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])]) 
    new_col_matched = np.array([new_col[matched['row'][idx], matched['col'][idx]]
                       for idx in range(matched['row'].shape[0])]) 
    new_row_col = {'row': new_row_matched, 'col': new_col_matched}
    from calipso import get_channel_data_from_object
    nobj.darkest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)
    nobj.darkest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)
    nobj.darkest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)
    nobj.darkest_r06=get_channel_data_from_object(dataObj, '06', new_row_col)
    nobj.darkest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)
    nobj.darkest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)
    return nobj
    



def get_channel_data_from_object_all(dataObj, chn_des, nodata=-9):
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
        raise ValueError

    temp = channels[chnum].data
    chdata = channels[chnum].data* channels[chnum].gain + channels[chnum].intercept      

    chdata= np.where(
        np.logical_or(
            np.equal(temp, dataObj.nodata),
            np.equal(temp, dataObj.missing_data)),
        nodata, chdata)
    return chdata
