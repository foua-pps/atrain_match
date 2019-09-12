# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
#Change log is found in git

import logging
import numpy as np
logger = logging.getLogger(__name__)
from matchobject_io import (DataObject,
                            ppsImagerObject,
                            IssObject,
                            TruthImagerTrackObject)                            
import config
from utils.common import (MatchupError, ProcessingError,
                    elements_within_range)
from libs.extract_imager_along_track import imager_track_from_matched

from truths.calipso import (find_break_points)
from utils.runutils import do_some_logging
import time
import datetime
import calendar

def get_iss(filename):
    # Read ISS Radar data for calipso something is done in this function:
    limit = 0.02
    iss = read_iss(filename)
    iss.cloud_fraction = 0.0 + 0*iss.cloud_phase_fore_fov[:,0].copy()
    iss.total_optical_depth_5km = 0.0 + 0*iss.cloud_phase_fore_fov[:,0].copy()
    iss.validation_height = -9 + 0 * iss.cloud_fraction.copy()
    #  feature_optical_depth_1064_fore_fov
    #1.       Handle layers with values -1 och -999.9 as optical thick (= 5)
    #2.       Check that each layer is a cloud layer (variable Feature_Type)
    #3.       For safty, check also Feature_Type_Score (>5.0) and Extinction_QC_Flag_1064==0

    #used for cloud height validation, at certain modes it might be updated.
    # Start from bottom and find hight of highest cloud. 
    # Remember that layer_top_altitude contain also aerosols
    # This was not the case for CALIPSO-data
    # There can be a thin aerosol layer above the highest cloud!
    for layer in range(9,-1,-1): #layer index 9.0
        #is_cloudy = np.greater(iss.cloud_phase_fore_fov[:,layer],0)
        is_cloudy = np.equal(iss.feature_type_fore_fov[:,layer],1)
        is_cloudy = np.logical_and(
            is_cloudy,
            np.equal(iss.extinction_qc_flag_1064_fore_fov[:,layer],0))
        is_cloudy = np.logical_and(
            is_cloudy,
            np.greater(iss.feature_type_score_fore_fov[:,layer],5))
        od_layer = iss.feature_optical_depth_1064_fore_fov[:,layer].copy()
        od_layer[od_layer==-1] = 5.0
        od_layer[od_layer==-999.9] = 5.0
        od_layer[od_layer<0] = 0.0
        iss.total_optical_depth_5km += od_layer
        is_not_very_thin = np.greater(od_layer, limit)
        #is_cloudy = np.logical_and(is_cloudy, is_not_very_thin) 
        height_layer = iss.layer_top_altitude_fore_fov[:,layer]
        iss.validation_height[is_cloudy] = height_layer[is_cloudy]

    logger.warning("Currently not considering cloudheight for cloud layers "
                   "thinner than %3.2f (feature_optical_depth_1064_fore_fov)", limit)
    iss.validation_height[iss.validation_height>=0] = iss.validation_height[iss.validation_height>=0]*1000
    iss.validation_height[iss.validation_height<0] = -9
    #  --- cloud_phase_fore_fov ---
    #0 clear or aerosol layer.
    #1 water cloud
    #2 unkknown phase
    #3 ice cloud
    # --- feature_type_score ---
    #| 10 | = high confidence
    #| 1 | = low confidence
    #0 = zero confidence
    # --- feature_type ---
    #0 = invalid
    #1 = cloud    
    #2 = undetermined    
    #3 = aerosol
    # --- sky_condition_fore_fov ---
    #0 = clean skies (no clouds/aerosols)
    #1 = clear skies (no clouds)       
    #2 = cloudy skies (no aerosols)       
    #3 = hazy/cloudy (both clouds/aerosols)
    iss.cloud_fraction = np.where(iss.sky_condition_fore_fov>1.5, 1.0,0.0)
    #print "iss cf", iss.cloud_fraction
    logger.warning("Currently using sky_condition_fore_fov to set cloudfraction,"
                   "\n \t use cloud_phase_fore_fov instead?")

    iss.validation_height[iss.cloud_fraction<0.5] = -9 #should not be needed
    return iss



scip_these_larger_variables_until_needed_iss = {
        # if any of these are needed remove them from the dictionary!  
        # might need more work as they are 3D or 4D variables.
        "Constrained_Lidar_Ratio_Flag": True,
        "Attenuated_Backscatter_Statistics_1064_Fore_FOV": True,
        "Attenuated_Backscatter_Statistics_532_Fore_FOV":  True,
        "Attenuated_Total_Color_Ratio_Statistics_Fore_FOV": True, 
        "Volume_Depolarization_Ratio_Statistics_1064_Fore_FOV": True,
        }
def read_iss(filename):
    import h5py
    import pdb, sys
    logger.info("Reading file %s", filename)
    retv = IssObject()
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        for group in h5file.keys():
            if group in "metadata_parameters":
                continue
            my_group = h5file[group]

            for dataset in my_group.keys():
                name = dataset.lower()
                if dataset in scip_these_larger_variables_until_needed_iss.keys():
                    continue
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


def match_iss_imager(issObj,imagerGeoObj,imagerObj,ctype,cma,ctth,nwp,imagerAngObj, 
                         cpp, nwp_segments, SETTINGS):
    retv = TruthImagerTrackObject(truth='iss')
    retv.imager_instrument = imagerGeoObj.instrument.lower()
    retv.iss = issObj

    from utils.common import map_imager
    cal, cap = map_imager(imagerGeoObj, issObj.longitude.ravel(),
                         issObj.latitude.ravel(),
                         radius_of_influence=config.RESOLUTION*0.7*1000.0) # larger than radius...
    calnan = np.where(cal == config.NODATA, np.nan, cal)

    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel 
                              in zip(cal,cap)]
        imager_lines_sec_1970 = np.where(cal != config.NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != config.NODATA, imagerGeoObj.time[cal], np.nan)
    # Find all matching Iss pixels within +/- sec_timeThr from the IMAGER data
    idx_match = elements_within_range(issObj.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr"])

    if idx_match.sum() == 0:
        logger.warning("No matches in region within time threshold %d s.", SETTINGS["sec_timeThr"])
        return None

    retv.iss = retv.iss.extract_elements(idx=idx_match)

     # Calipso line,pixel inside IMAGER swath:
    retv.iss.imager_linnum = np.repeat(cal, idx_match).astype('i')
    retv.iss.imager_pixnum = np.repeat(cap, idx_match).astype('i')
    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.iss.sec_1970 - retv.imager.sec_1970                        

    do_some_logging(retv, issObj)
    logger.info("Generate the latitude,cloudtype tracks!")
    retv = imager_track_from_matched(retv, SETTINGS, imagerGeoObj, imagerObj, imagerAngObj, 
                                    nwp, ctth, ctype, cma,  
                                    cpp=cpp, nwp_segments=nwp_segments)
    return retv



def reshapeIss(issfiles, imager, SETTINGS):
    import sys
    imager_end = imager.sec1970_end
    imager_start = imager.sec1970_start
    iss = get_iss(issfiles[0])
    for i in range(len(issfiles)-1):
        newIss = get_iss(issfiles[i+1])
        iss_start_all = iss.sec_1970.ravel()
        iss_new_all = newIss.sec_1970.ravel()
        if not iss_start_all[0]<iss_new_all[0]:
            raise ProcessingError("Iss files are in the wrong order")
        iss_break = np.argmin(np.abs(iss_start_all - iss_new_all[0]))+1
        # Concatenate the feature values
        iss = iss + newIss
 
    # Finds Break point
    startBreak, endBreak = find_break_points(iss, imager, SETTINGS)
    iss = iss.extract_elements(starti=startBreak, 
                               endi=endBreak) 
    return iss
    

if __name__ == "__main__":
    pass
