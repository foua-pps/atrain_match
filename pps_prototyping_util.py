import numpy as np
import logging
logger = logging.getLogger(__name__)

from extract_imager_along_track import CHANNEL_MICRON_AVHRR_PPS, CHANNEL_MICRON_DESCRIPTIONS
from extract_imager_along_track import get_channel_data_from_object

def get_t11t12_texture_data_from_object(dataObj, nwp_obj, ch11, ch12, text_name):
    #https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    t11 = get_channel_data_from_objectfull_resolution(dataObj, ch11, nodata=-9)
    t12 = get_channel_data_from_objectfull_resolution(dataObj, ch12, nodata=-9)    
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
    t11 = get_channel_data_from_objectfull_resolution(dataObj, '11', nodata=-9)
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
    nobj.warmest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)[0]
    nobj.warmest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)[0]
    nobj.warmest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)[0]
    nobj.warmest_r06=get_channel_data_from_object(dataObj, '06', new_row_col)[0]
    nobj.warmest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)[0]
    nobj.warmest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)[0]
    return nobj

def get_coldest_values(dataObj, matched):
    nobj = NeighbourObj()
    t11 = get_channel_data_from_objectfull_resolution(dataObj, '11', nodata=-9)
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
    nobj.coldest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)[0]
    nobj.coldest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)[0]
    nobj.coldest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)[0]
    nobj.coldest_r06=get_channel_data_from_object(dataObj, '06', new_row_col)[0]
    nobj.coldest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)[0]
    nobj.coldest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)[0]
    return nobj

def get_darkest_values(dataObj, matched):
    nobj = NeighbourObj()
    r09 = get_channel_data_from_objectfull_resolution(dataObj, '09', nodata=-9)
    r06 = get_channel_data_from_objectfull_resolution(dataObj, '06', nodata=-9)
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
    nobj.darkest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)[0]
    nobj.darkest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)[0]
    nobj.darkest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)[0]
    nobj.darkest_r06=get_channel_data_from_object(dataObj, '06', new_row_col)[0]
    nobj.darkest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)[0]
    nobj.darkest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)[0]
    return nobj
    
def get_channel_data_from_objectfull_resolution(dataObj, chn_des, nodata=-9):
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
    chdata[np.logical_or(np.equal(temp, dataObj.nodata),
                         np.equal(temp, dataObj.missing_data))]= nodata
    return chdata

