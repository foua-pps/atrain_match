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
                    NODATA, NLINES, SWATHWD, 
                    _validation_results_dir)
from common import (MatchupError, 
                    elements_within_range)
from extract_imager_along_track import avhrr_track_from_matched

from calipso import (find_break_points, calipso_track_from_matched,
                     time_reshape_calipso, do_some_logging)


def get_cloudsat(filename):
    # Read CLOUDSAT Radar data for calipso something is done in this function:
    cloudsat = read_cloudsat(filename)
    return cloudsat

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
    logger.info("Generate the latitude,cloudtype tracks!")
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
