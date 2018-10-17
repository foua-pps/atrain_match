#Change log is found in git
import time
import numpy as np
import os
import logging
logger = logging.getLogger(__name__)
from matchobject_io import (CloudsatObject,
                            CloudsatImagerTrackObject)                            
from config import (sec_timeThr, RESOLUTION,
                    NODATA, CLOUDSAT_CLOUDY_THR, CLOUDSAT_REQUIRED)
from common import (MatchupError, ProcessingError,
                    elements_within_range)
from extract_imager_along_track import imager_track_from_matched

from calipso import (find_break_points, calipso_track_from_matched,
                     time_reshape_calipso)
from runutils import do_some_logging

def add_validation_ctth_cloudsat(cloudsat):
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
    # Read CLOUDSAT Radar data and add some variables
    if '.hdf' in filename:
        cloudsat = read_cloudsat_hdf4(filename)
    elif '.h5' in filename: 
        cloudsat = read_cloudsat(filename)
    else:
        raise MatchupError(
            "Missing reader for CloudSat file type {:s}".format(
                os.path.basename(fielname)))    
    if 'GEOPROF' in os.path.basename(filename):     
        cloudsat = add_validation_ctth_cloudsat(cloudsat)
        cloudsat = add_cloudsat_cloud_fraction(cloudsat)      
    return cloudsat

def clsat_name_conversion(dataset_name_in_cloudsat_file, retv):
    am_name = dataset_name_in_cloudsat_file
    if dataset_name_in_cloudsat_file == 'DEM_elevation':
        am_name = 'elevation'
    elif dataset_name_in_cloudsat_file == 'Sigma-Zero':
        am_name = 'SigmaZero'
    elif dataset_name_in_cloudsat_file == 'Longitude':
        am_name = 'longitude'
    elif dataset_name_in_cloudsat_file == 'Latitude':
        am_name = 'latitude'
    elif am_name.lower() in retv.all_arrays.keys():    
        am_name = am_name.lower()    
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
        #2D data, print idx, sds 
        data = h4file.select(sds).get()
        #print h4file.select(sds).attributes().keys()
        am_name = clsat_name_conversion(sds, retv)
        if am_name in retv.all_arrays.keys():
            retv.all_arrays[am_name] = convert_data(data)
        #print h4file.select(sds).info()
    h4file = HDF(filename, SDC.READ)
    vs = h4file.vstart()
    data_info_list = vs.vdatainfo()
    for item in data_info_list:
        #1D data compound/Vdata
        name = item[0]
        data_handle = vs.attach(name)
        data = np.array(data_handle[:])
        attrinfo_dic = data_handle.attrinfo()
        factor = data_handle.findattr('factor')
        offset = data_handle.findattr('offset')
        #print data_handle.factor
        am_name = clsat_name_conversion(name, retv)
        if am_name in retv.all_arrays.keys():
            #To save RAM and disk only read what we use!
            if factor is None and offset is None:
                retv.all_arrays[am_name] = convert_data(data)
            elif np.float(factor.get()) == 1.0 and np.float(offset.get()) == 0.0:
                retv.all_arrays[am_name] = convert_data(data)
            else:
                if factor is None:
                    factor = 1.0
                if offset is None:
                    offset = 0.0

                raise MatchupError("Not default offset and factor. Fix code")
                #The code below is probably ok, but make sure:
                #the_data_scaled = convert_data(data)*factor + offset
                #retv.all_arrays[am_name] = the_data_scaled 
        data_handle.detach()
    #print data_handle.attrinfo()
    h4file.close()
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    retv.sec_1970 = retv.Profile_time.ravel() + retv.TAI_start + dsec
    return retv

def read_cloudsat(filename):
    import h5py 
    CLOUDSAT_TYPE = "GEOPROF"
    if 'CWC-RVOD' in os.path.basename(filename):
        CLOUDSAT_TYPE = 'CWC-RVOD'
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


def match_cloudsat_imager(cloudsatObj,imagerGeoObj,imagerObj,ctype,cma,ctth,nwp,imagerAngObj, 
                         cpp, nwp_segments):
    retv = CloudsatImagerTrackObject()

    #Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    from common import map_imager
    cal, cap = map_imager(imagerGeoObj, 
                         cloudsatObj.longitude.ravel(),
                         cloudsatObj.latitude.ravel(),
                         radius_of_influence=RESOLUTION*0.7*1000.0) # somewhat larger than radius...
    calnan = np.where(cal == NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        if CLOUDSAT_REQUIRED:
            raise MatchupError("No matches within region.")
        else:
            logger.warning("No matches within region.")
            return None
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal,cap)]
        imager_lines_sec_1970 = np.where(cal != NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != NODATA, imagerGeoObj.time[cal], np.nan)
    # Find all matching Cloudsat pixels within +/- sec_timeThr from the IMAGER data
    idx_match = elements_within_range(cloudsatObj.sec_1970, imager_lines_sec_1970, sec_timeThr)

    if idx_match.sum() == 0:
        if CLOUDSAT_REQUIRED:            
            raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)  
        else:
            logger.warning("No matches in region within time threshold %d s.", sec_timeThr)
            return None
    retv.cloudsat = calipso_track_from_matched(retv.cloudsat, cloudsatObj, idx_match)
 
    # Cloudsat line,pixel inside IMAGER swath:
    retv.cloudsat.imager_linnum = np.repeat(cal, idx_match).astype('i')
    retv.cloudsat.imager_pixnum = np.repeat(cap, idx_match).astype('i')

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.cloudsat.sec_1970 - retv.imager.sec_1970
    do_some_logging(retv, cloudsatObj)
    logger.debug("Generate the latitude,cloudtype tracks!")
    retv = imager_track_from_matched(retv, imagerGeoObj, imagerObj, imagerAngObj, 
                                    nwp, ctth, ctype, cma, 
                                    cpp=cpp, nwp_segments=nwp_segments)
    return retv

def mergeCloudsat(cloudsat, cloudsatlwp):
    #map cloudsat_lwp to cloudsat
    from amsr_imager.match import match_lonlat
    source = (cloudsatlwp.longitude.astype(np.float64).reshape(-1,1), 
              cloudsatlwp.latitude.astype(np.float64).reshape(-1,1))
    target = (cloudsat.longitude.astype(np.float64).reshape(-1,1), 
              cloudsat.latitude.astype(np.float64).reshape(-1,1))
    mapper, dummy = match_lonlat(source, target, radius_of_influence=10, n_neighbours=1)
    cloudsat_lwp_index = mapper.rows.filled(NODATA).ravel()
    # Transfer CloudSat LWP to ordinary cloudsat obj
    index = cloudsat_lwp_index.copy()
    index[index<0] = 0
    cloudsat.RVOD_liq_water_path = np.where(
        cloudsat_lwp_index>=0, 
        cloudsatlwp.RVOD_liq_water_path[index],-9)
    cloudsat.RVOD_ice_water_path = np.where(
        cloudsat_lwp_index>=0, 
        cloudsatlwp.RVOD_ice_water_path[index],-9)
    cloudsat.LO_RVOD_liquid_water_path = np.where(
        cloudsat_lwp_index>=0, 
        cloudsatlwp.LO_RVOD_liquid_water_path[index],-9)
    cloudsat.IO_RVOD_ice_water_path = np.where(
        cloudsat_lwp_index>=0, 
        cloudsatlwp.IO_RVOD_ice_water_path[index],-9)
    cloudsat.RVOD_CWC_status = np.where(
        cloudsat_lwp_index>=0, 
        cloudsatlwp.RVOD_CWC_status[index],-9)
    return cloudsat


def reshapeCloudsat(cloudsatfiles, imager):
    #imager_end = imager.sec1970_end
    #imager_start = imager.sec1970_start
    clsat = get_cloudsat(cloudsatfiles[0])
    for i in range(len(cloudsatfiles)-1):
        newCloudsat = get_cloudsat(cloudsatfiles[i+1])
        clsat_start_all = clsat.sec_1970.ravel()
        clsat_new_all = newCloudsat.sec_1970.ravel()
        if not clsat_start_all[0]<clsat_new_all[0]:
            raise ProcessingError("CloudSat files are in the wrong order!")
        clsat_break = np.argmin(np.abs(clsat_start_all - clsat_new_all[0]))+1
        # Concatenate the feature values
        #arname = array name from cloudsatObj
        for arname, value in clsat.all_arrays.items(): 
            if value is not None:
                if value.size != 1:
                    clsat.all_arrays[arname] = np.concatenate((value[0:clsat_break,...],
                                                               newCloudsat.all_arrays[arname]))
    # Finds Break point
    startBreak, endBreak = find_break_points(clsat, imager)
    clsat = time_reshape_calipso(clsat, startBreak, endBreak)
    return clsat
    
    

if __name__ == "__main__":
    # Testing not tested for many years! 
    pass 

