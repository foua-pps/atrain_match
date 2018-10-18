import numpy as np
import logging
logger = logging.getLogger(__name__)

from libs.extract_imager_along_track import CHANNEL_MICRON_IMAGER_PPS, CHANNEL_MICRON_DESCRIPTIONS
from libs.extract_imager_along_track import get_channel_data_from_object, get_data_from_array

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
        self.extra_info_sza_corr = ""

def get_data_from_array_fill_outside(array, matched, Fill=0):
    row_index = matched['row'].copy()
    col_index = matched['col'].copy()
    row_lim, col_lim = array.shape
    outside = np.logical_or(np.logical_or(row_index<0,
                                          row_index>=row_lim),
                            np.logical_or(col_index<0,
                                          col_index>=col_lim))
    row_index[outside] = 0
    col_index[outside] = 0
    temp = np.array([array[row_index[idx], col_index[idx]]
                     for idx in range(matched['row'].shape[0])]) 
    return np.where(outside, Fill, temp)


def get_warmest_or_coldest_index(t11, matched, warmest=True):
    FILL = 999999.9  #coldest
    if warmest:
        FILL = -99
        
    steps = [(i ,j) for i in [-2,-1,0,1,2] for j in  [-2,-1,0,1,2] ]
    t11_neighbour_i = np.zeros((25, matched['row'].shape[0]) )
    for i, (step_r, step_c) in enumerate(steps):
        new_row_col = {'row': matched['row'] + step_r, 
                       'col': matched['col'] + step_c}
        t11_neighbour_i[i,:] = get_data_from_array_fill_outside(t11, 
                                                                new_row_col, 
                                                                Fill=FILL)
    if warmest:    
        neigbour_index = np.argmax(t11_neighbour_i, axis=0)
    else: #coldest
        neigbour_index = np.argmin(t11_neighbour_i, axis=0)
    new_row_matched = np.array(
        [matched['row'][idx] + steps[neigbour_index[idx]][0] 
         for idx in  range(matched['row'].shape[0])] )
    new_col_matched = np.array(
        [matched['col'][idx] + steps[neigbour_index[idx]][1] 
         for idx in  range(matched['row'].shape[0])] )
    new_row_col = {'row': new_row_matched, 'col': new_col_matched}
    return new_row_col


def get_warmest_values(dataObj, matched):
    nobj = NeighbourObj()
    t11 = get_channel_data_from_objectfull_resolution(dataObj, '11', nodata=-9)
    new_row_col = get_warmest_or_coldest_index(t11,  matched)
    nobj.warmest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)[0]
    nobj.warmest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)[0]
    nobj.warmest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)[0]
    nobj.warmest_r06, nobj.extra_info_sza_corr = get_channel_data_from_object(
        dataObj, '06', new_row_col)
    nobj.warmest_r16 = get_channel_data_from_object(dataObj, '16', new_row_col)[0]
    nobj.warmest_r09 = get_channel_data_from_object(dataObj, '09', new_row_col)[0]
    return nobj

def get_coldest_values(dataObj, matched):
    nobj = NeighbourObj()
    t11 = get_channel_data_from_objectfull_resolution(dataObj, '11', nodata=-9)
    new_row_col = get_warmest_or_coldest_index(t11, matched, warmest=False)
    nobj.coldest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)[0]
    nobj.coldest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)[0]
    nobj.coldest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)[0]
    nobj.coldest_r06, nobj.extra_info_sza_corr = get_channel_data_from_object(
        dataObj, '06', new_row_col)
    nobj.coldest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)[0]
    nobj.coldest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)[0]
    return nobj

def get_darkest_values(dataObj, matched):
    nobj = NeighbourObj()
    r09 = get_channel_data_from_objectfull_resolution(dataObj, '09', nodata=-9)
    r06 = get_channel_data_from_objectfull_resolution(dataObj, '06', nodata=-9)
    darkest = np.where(r09>r06,r06,r09)
    new_row_col = get_warmest_or_coldest_index(darkest, matched, warmest=False)
    nobj.darkest_t11=get_channel_data_from_object(dataObj, '11', new_row_col)[0]
    nobj.darkest_t12=get_channel_data_from_object(dataObj, '12', new_row_col)[0]
    nobj.darkest_t37=get_channel_data_from_object(dataObj, '37', new_row_col)[0]
    nobj.darkest_r06, nobj.extra_info_sza_corr=get_channel_data_from_object(dataObj, '06', new_row_col)
    nobj.darkest_r16=get_channel_data_from_object(dataObj, '16', new_row_col)[0]
    nobj.darkest_r09=get_channel_data_from_object(dataObj, '09', new_row_col)[0]
    return nobj
    
def get_channel_data_from_objectfull_resolution(dataObj, chn_des, nodata=-9):
    """Get the IMAGER/VIIRS channel data on the track
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

