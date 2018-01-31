#Change log is found in git
import time
import numpy as np
import logging
logger = logging.getLogger(__name__)
from matchobject_io import (DataObject,
                            ppsAvhrrObject,
                            CloudsatObject,
                            CloudsatAvhrrTrackObject)                            
from config import (AREA, sec_timeThr, RESOLUTION,
                    NODATA, CLOUDSAT_CLOUDY_THR, 
                    _validation_results_dir)
from common import (MatchupError, 
                    elements_within_range)
from extract_imager_along_track import avhrr_track_from_matched

from calipso import (find_break_points, calipso_track_from_matched,
                     time_reshape_calipso, do_some_logging)

def add_validation_ctth_cloudsat(cloudsat):
    import config
    #CLOUDSAT VALIDATION HEIGHT!
    #The restriction to use only pixels where imager and cloudsat
    #both is cloudy is done later. This makes it possbile to find also
    #POD-cloudy for different cloud heights
    LARGE_POSITIVE = 99999.0
    validation_height = -9 + 0*np.zeros(cloudsat.latitude.shape)
    validation_height_base =  LARGE_POSITIVE + 0*np.zeros(cloudsat.latitude.shape)
    for i in range(125):
        height = cloudsat.Height[:,i]
        cmask_ok = cloudsat.CPR_Cloud_mask[:,i]
        top_height = height+120
        #top_height[height<240*4] = -9999 #Do not use not sure why these were not used Nina 20170317
        is_cloudy = cmask_ok > CLOUDSAT_CLOUDY_THR
        top_height[~is_cloudy] = -9999
        validation_height[validation_height<top_height] =  top_height[
            validation_height<top_height]
        cloudsat.validation_height= validation_height
        update_base = np.logical_and(top_height>0, validation_height_base>top_height)
        validation_height_base[update_base] =  top_height[update_base]
        cloudsat.validation_height_base = validation_height_base
        cloudsat.validation_height_base[cloudsat.validation_height_base>=LARGE_POSITIVE] =-9 
    return cloudsat

def add_cloudsat_cloud_fraction(cloudsat):
    cloudsat_cloud_mask = cloudsat.CPR_Cloud_mask
    cloudsat_cloud_mask = np.greater_equal(cloudsat_cloud_mask, 
                                           CLOUDSAT_CLOUDY_THR)
    cloudsat_cloud_fraction = np.zeros(cloudsat.latitude.shape[0])
    sum_cloudsat_cloud_mask = np.sum(cloudsat_cloud_mask, axis=1)
    if len(sum_cloudsat_cloud_mask) != (len(cloudsat_cloud_fraction)):
        raise ValueError('Boolean index-array should have same lenght as array!')
    cloudsat_cloud_fraction[sum_cloudsat_cloud_mask > 2] = 1 # requires at least two cloudy bins
    cloudsat.cloud_fraction = cloudsat_cloud_fraction
    return cloudsat

def get_cloudsat(filename):
    # Read CLOUDSAT Radar data for calipso something is done in this function:
    if '.h5' in filename: 
        cloudsat = read_cloudsat(filename)
    else:
        cloudsat = read_cloudsat_hdf4(filename)
    cloudsat = add_validation_ctth_cloudsat(cloudsat)
    cloudsat = add_cloudsat_cloud_fraction(cloudsat)      
    return cloudsat

def clsat_name_conversion(dataset_name_in_cloudsat_file):
    am_name = dataset_name_in_cloudsat_file
    if dataset_name_in_cloudsat_file == 'DEM_elevation':
        am_name = 'elevation'
    if dataset_name_in_cloudsat_file == 'Sigma-Zero':
        am_name = 'SigmaZero'
    if dataset_name_in_cloudsat_file == 'Longitude':
        am_name = 'longitude'
    if dataset_name_in_cloudsat_file == 'Latitude':
        am_name = 'latitude'
    return am_name    
def read_cloudsat_hdf4(filename):
    from pyhdf.SD import SD, SDC
    from pyhdf.HDF import HDF, HC
    import pyhdf.VS 
    def convert_data(data):
        if len(data.shape) == 2:
            if data.shape[1] == 1:
                return data[:, 0]
            elif data.shape[0] == 1:
                return data[0, :]
        return data
    retv = CloudsatObject()
    h4file = SD(filename, SDC.READ)
    datasets = h4file.datasets()
    attributes = h4file.attributes()
    #for idx,attr in enumerate(attributes.keys()):
    #    print idx, attr
    for idx,sds in enumerate(datasets.keys()):
        #print idx, sds
        data = h4file.select(sds).get()
        #print h4file.select(sds).attributes().keys()
        retv.all_arrays[sds] = convert_data(data)
        #print h4file.select(sds).info()
        

    h4file = HDF(filename, SDC.READ)
    vs = h4file.vstart()
    data_info_list = vs.vdatainfo()
    for item in data_info_list:
        name = item[0]
        data_handle = vs.attach(name)
        data = np.array(data_handle[:])
        attrinfo_dic = data_handle.attrinfo()
        factor = data_handle.findattr('factor')
        offset = data_handle.findattr('offset')
        #print data_handle.factor
        am_name = clsat_name_conversion(name)
        if factor is None and offset is None:
            retv.all_arrays[am_name] = convert_data(data)
        else:
            if factor is None:
                factor = 1.0
            if offset is None:
                offset = 0.0
            raise MatchupError("Not default offset and factor. Fix code")
            #The code below is probably ok, but make sure:
            the_data_scaled = convert_data(data)*factor + offset
            retv.all_arrays[am_name] = the_data_scaled 
        
        
        data_handle.detach()
    #print data_handle.attrinfo()
    h4file.close()
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    retv.sec_1970 = retv.Profile_time.ravel() + retv.TAI_start + dsec
    return retv

def read_cloudsat(filename):
    import h5py 
    from config import CLOUDSAT_TYPE
    def get_data(dataset):
        type_name = dataset.value.dtype.names
        try:
            data = dataset.value[type_name[0]]
        except TypeError:
            data = dataset.value
        # Convert 1-dimensional matrices to 1-d arrays
        if len(data.shape) == 2:
            if data.shape[1] == 1:
                return data[:, 0]
            elif data.shape[0] == 1:
                return data[0, :]
        return data
    retv = CloudsatObject()
    h5file = h5py.File(filename, 'r')
    root="2B-" + CLOUDSAT_TYPE
    for group in ['Geolocation Fields', 'Data Fields']:
        tempG = h5file["%s/%s" % (root, group)]
        for dataset in tempG.keys():
            if dataset in retv.all_arrays.keys():
                retv.all_arrays[dataset] = get_data(tempG[dataset])
            elif dataset.lower() in retv.all_arrays.keys():
                retv.all_arrays[dataset.lower()] = get_data(tempG[dataset])
            elif dataset == 'DEM_elevation':
                retv.all_arrays['elevation'] = get_data(tempG[dataset])
            elif dataset == 'Sigma-Zero':           
                retv.all_arrays['SigmaZero'] = get_data(tempG[dataset])
    h5file.close()    
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    retv.sec_1970 = retv.Profile_time.ravel() + retv.TAI_start + dsec
    return retv


def match_cloudsat_avhrr(cloudsatObj,imagerGeoObj,imagerObj,ctype,cma,ctth,nwp,imagerAngObj, 
                         cpp, nwp_segments):
    import string
    import os
    retv = CloudsatAvhrrTrackObject()
    lonCloudsat = cloudsatObj.longitude.ravel()
    latCloudsat = cloudsatObj.latitude.ravel()
    timeCloudsat = cloudsatObj.sec_1970.ravel()

    #Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    from common import map_avhrr
    cal, cap = map_avhrr(imagerGeoObj, lonCloudsat.ravel(), latCloudsat.ravel(),
                         radius_of_influence=RESOLUTION*0.7*1000.0) # somewhat larger than radius...
    calnan = np.where(cal == NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal,cap)]
        imager_lines_sec_1970 = np.where(cal != NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != NODATA, imagerGeoObj.time[cal], np.nan)
    # Find all matching Cloudsat pixels within +/- sec_timeThr from the AVHRR data
    idx_match = elements_within_range(cloudsatObj.sec_1970, imager_lines_sec_1970, sec_timeThr)

    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)  

    retv.cloudsat = calipso_track_from_matched(retv.cloudsat, cloudsatObj, idx_match)
 
    # Cloudsat line,pixel inside AVHRR swath:
    cal_on_avhrr = np.repeat(cal, idx_match)
    cap_on_avhrr = np.repeat(cap, idx_match)
    retv.cloudsat.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.cloudsat.avhrr_pixnum = cap_on_avhrr.astype('i')

   # Imager time
    if len(imagerGeoObj.time.shape)>1:
        retv.avhrr.sec_1970= [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal_on_avhrr,cap_on_avhrr)]
    else:
        retv.avhrr.sec_1970 = imagerGeoObj.time[cal_on_avhrr]
    retv.diff_sec_1970 = retv.cloudsat.sec_1970 - retv.avhrr.sec_1970
    do_some_logging(retv, cloudsatObj)
    logger.debug("Generate the latitude,cloudtype tracks!")
    retv = avhrr_track_from_matched(retv, imagerGeoObj, imagerObj, imagerAngObj, 
                                    nwp, ctth, ctype, cma, cal_on_avhrr, cap_on_avhrr, 
                                    cpp=cpp, nwp_segments=nwp_segments)
    return retv

def reshapeCloudsat(cloudsatfiles, avhrr):
    import sys
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    clsat = get_cloudsat(cloudsatfiles[0])
    for i in range(len(cloudsatfiles)-1):
        newCloudsat = get_cloudsat(cloudsatfiles[i+1])
        clsat_start_all = clsat.sec_1970.ravel()
        clsat_new_all = newCloudsat.sec_1970.ravel()
        if not clsat_start_all[0]<clsat_new_all[0]:
            print "cloudsat files are in the wrong order"
            print("Program cloudsat.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9)
        clsat_break = np.argmin(np.abs(clsat_start_all - clsat_new_all[0]))+1
        # Concatenate the feature values
        #arname = array name from cloudsatObj
        for arname, value in clsat.all_arrays.items(): 
            if value is not None:
                if value.size != 1:
                    clsat.all_arrays[arname] = np.concatenate((value[0:clsat_break,...],newCloudsat.all_arrays[arname]))
    # Finds Break point
    startBreak, endBreak = find_break_points(clsat, avhrr)
    clsat = time_reshape_calipso(clsat, startBreak, endBreak)
    return clsat
    
    

if __name__ == "__main__":
    # Testing:
#    import string
#    import epshdf
#    import pps_io
    
    from config import CLOUDSAT_DIR
    cloudsatfile = "%s/2007151082929_05796_CS_2B-GEOPROF_GRANULE_P_R04_E02.h5"%(CLOUDSAT_DIR)

    # --------------------------------------------------------------------
    logger.info("Read CLOUDSAT data")
    # Read CLOUDSAT Radar data:
    cloudsat = get_cloudsat(cloudsatfile)
    #cloudsat = read_cloudsat(cloudsatfile)

    lonCloudsat = cloudsat.longitude.ravel()
    latCloudsat = cloudsat.latitude.ravel()

    # Test:
    ndim = lonCloudsat.ravel().shape[0]
    idx_match = np.zeros((ndim,),'b')
    idx_match[0:10] = 1

    x = np.repeat(cloudsat.Height[::,0],idx_match)
    for i in range(1,125):
        x = np.concatenate((x,np.repeat(cloudsat.Height[::,i],idx_match)))
    N = x.shape[0]/125
    cloudsat.Height = np.reshape(x,(125,N))
