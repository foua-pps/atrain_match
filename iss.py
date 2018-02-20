#Change log is found in git

import logging
import numpy as np
logger = logging.getLogger(__name__)
from matchobject_io import (DataObject,
                            ppsAvhrrObject,
                            IssObject,
                            IssAvhrrTrackObject)                            
from config import (AREA, sec_timeThr, RESOLUTION,
                    NODATA,
                    _validation_results_dir)
from common import (MatchupError, ProcessingError,
                    elements_within_range)
from extract_imager_along_track import avhrr_track_from_matched

from calipso import (find_break_points, calipso_track_from_matched,
                     time_reshape_calipso)
from runutils import do_some_logging
import time
import datetime
import calendar

def get_iss(filename):
    # Read ISS Radar data for calipso something is done in this function:
    limit = 0.02
    iss = read_iss(filename)
    #0 clear or aerosol layer.
    #1 water cloud
    #2 unkknown phase
    #3 ice cloud
    iss.cloud_fraction = 0*iss.cloud_phase_fore_fov[:,0].copy()
    for layer in range(9,-1,-1): #layer index 9.0
        is_cloudy = np.greater(iss.cloud_phase_fore_fov[:,layer],0)
        is_not_very_thin = np.greater(iss.feature_optical_depth_1064_fore_fov[:,layer], limit)
        is_cloudy = np.logical_and(is_cloudy, is_not_very_thin) 
        iss.cloud_fraction[is_cloudy] = 1.0
    #0 clear
    #1 aersol only
    #2 cloud only
    #3 cloud and aerosol    
    iss.cloud_fraction = np.where(iss.sky_condition_fore_fov>1.5, 1.0,0.0)
    #print "iss cf", iss.cloud_fraction
    logger.warning("Currently using sky_condition_fore_fov to set cloudfraction,"
                   "\n \t use cloud_phase_fore_fov instead?")

    #used for cloud height validation, at cirtain modis it might be updated.
    # Start from bottom and find hight of highest cloud. 
    # Remember that layer_top_altitude contain also aerosols
    # This was not the case for CALIPSO-data
    # There can be a thin aerosol layer above the highest cloud!
    iss.validation_height = -9 + 0 * iss.cloud_fraction.copy()
    for layer in range(9,-1,-1): #layer index 9..0
        is_cloudy = np.greater(iss.cloud_phase_fore_fov[:,layer],0)
        is_not_very_thin = np.greater(iss.feature_optical_depth_1064_fore_fov[:,layer], limit)
        is_cloudy = np.logical_and(is_cloudy, is_not_very_thin)
        height_layer = iss.layer_top_altitude_fore_fov[:,layer]
        iss.validation_height[is_cloudy] = height_layer[is_cloudy]
    logger.warning("Currently not considering cloudheight for cloud layers "
                   "thinner than %3.2f (feature_optical_depth_1064_fore_fov)", limit)
    iss.validation_height[iss.validation_height>=0] = iss.validation_height[iss.validation_height>=0]*1000
    iss.validation_height[iss.validation_height<0] = -9
    iss.validation_height[iss.cloud_fraction<0.5] = -9 #should not be needed
    return iss


def read_iss(filename):
    import h5py
    import pdb, sys
    logger.info("Reading file %s", filename)
    scip_these_larger_variables_until_needed = {
        # if any of these are needed remove them from the dictionary!  
        # might need more work as they are 3D or 4D variables.
        "Constrained_Lidar_Ratio_Flag": True,
        "Attenuated_Backscatter_Statistics_1064_Fore_FOV": True,
        "Attenuated_Backscatter_Statistics_532_Fore_FOV":  True,
        "Attenuated_Total_Color_Ratio_Statistics_Fore_FOV": True, 
        "Volume_Depolarization_Ratio_Statistics_1064_Fore_FOV": True,
        }
    retv = IssObject()
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        for group in h5file.keys():
            if group in "metadata_parameters":
                continue
            my_group = h5file[group]

            for dataset in my_group.keys():
                if dataset in scip_these_larger_variables_until_needed.keys():
                    continue
                name = dataset.lower()
                data = my_group[dataset].value
                data = np.array(data)
                setattr(retv, name, data) 
        h5file.close()
        retv.latitude = retv.cats_fore_fov_latitude[:,1]
        retv.longitude = retv.cats_fore_fov_longitude[:,1] 
        # Elevation is given in km's. Convert to meters:
        retv.elevation = retv.dem_surface_altitude_fore_fov*1000.0 
        seconds_per_day = datetime.timedelta(days=1).total_seconds()
        midnight_sec_1970 = [ calendar.timegm(time.strptime(str(date_as_integer), "%Y%m%d")) for date_as_integer in retv.profile_utc_date]
        retv.sec_1970 = seconds_per_day*retv.profile_utc_time[:,1] + np.array(midnight_sec_1970)
    return retv  


def match_iss_avhrr(issObj,imagerGeoObj,imagerObj,ctype,cma,ctth,nwp,imagerAngObj, 
                         cpp, nwp_segments):
    import string
    import os
    retv = IssAvhrrTrackObject()
    lonIss = issObj.longitude.ravel()
    latIss = issObj.latitude.ravel()
    timeIss = issObj.sec_1970.ravel()

    from common import map_avhrr
    cal, cap = map_avhrr(imagerGeoObj, lonIss.ravel(), latIss.ravel(),
                         radius_of_influence=RESOLUTION*0.7*1000.0) # larger than radius...
    calnan = np.where(cal == NODATA, np.nan, cal)

    import matplotlib.pyplot as plt
    print lonIss.ravel(), latIss.ravel()
    plt.plot(lonIss.ravel(), latIss.ravel(),'b.')
    plt.plot(imagerGeoObj.longitude.ravel(), imagerGeoObj.latitude.ravel(),'r.')

    plt.savefig('iss.png')

    if (~np.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel 
                              in zip(cal,cap)]
        imager_lines_sec_1970 = np.where(cal != NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != NODATA, imagerGeoObj.time[cal], np.nan)
    # Find all matching Iss pixels within +/- sec_timeThr from the AVHRR data
    idx_match = elements_within_range(issObj.sec_1970, imager_lines_sec_1970, sec_timeThr)

    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." 
                           % sec_timeThr)  

    retv.iss = calipso_track_from_matched(retv.iss, issObj, idx_match)
 
    # Iss line,pixel inside AVHRR swath:
    cal_on_avhrr = np.repeat(cal, idx_match)
    cap_on_avhrr = np.repeat(cap, idx_match)
    retv.iss.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.iss.avhrr_pixnum = cap_on_avhrr.astype('i')

   # Imager time
    if len(imagerGeoObj.time.shape)>1:
        retv.avhrr.sec_1970= [imagerGeoObj.time[line,pixel] for line, pixel 
                              in zip(cal_on_avhrr,cap_on_avhrr)]
    else:
        retv.avhrr.sec_1970 = imagerGeoObj.time[cal_on_avhrr]
    retv.diff_sec_1970 = retv.iss.sec_1970 - retv.avhrr.sec_1970

    do_some_logging(retv, issObj)
    logger.info("Generate the latitude,cloudtype tracks!")
    retv = avhrr_track_from_matched(retv, imagerGeoObj, imagerObj, imagerAngObj, 
                                    nwp, ctth, ctype, cma, cal_on_avhrr, cap_on_avhrr, 
                                    cpp=cpp, nwp_segments=nwp_segments)
    return retv



def reshapeIss(issfiles, avhrr):
    import sys
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    iss = get_iss(issfiles[0])
    for i in range(len(issfiles)-1):
        newIss = get_iss(issfiles[i+1])
        iss_start_all = iss.sec_1970.ravel()
        iss_new_all = newIss.sec_1970.ravel()
        if not iss_start_all[0]<iss_new_all[0]:
            raise ProcessingError("Iss files are in the wrong order")
        iss_break = np.argmin(np.abs(iss_start_all - iss_new_all[0]))+1
        # Concatenate the feature values
        #arname = array name from issObj
        for arname, value in iss.all_arrays.items(): 
            if value is not None:
                if value.size != 1:
                    iss.all_arrays[arname] = np.concatenate((iss.all_arrays[arname],
                                                             newIss.all_arrays[arname]),axis=0)          
    # Finds Break point
    startBreak, endBreak = find_break_points(iss, avhrr)
    iss = time_reshape_calipso(iss, startBreak, endBreak)
    return iss
    

if __name__ == "__main__":
    pass
