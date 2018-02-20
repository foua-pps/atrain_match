import h5py
import numpy as np
from datetime import datetime
from calendar import timegm
TAI93 = datetime(1993, 1, 1)
from config import (RESOLUTION, sec_timeThr, NODATA)
from read_cloudproducts_and_nwp_pps import NWPObj
from matchobject_io import (AmsrAvhrrTrackObject, 
                            AmsrObject)
from calipso import (find_break_points, calipso_track_from_matched,
                     time_reshape_calipso, do_some_logging)

from common import (ProcessingError, MatchupError, elements_within_range)
from extract_imager_along_track import avhrr_track_from_matched
from amsr_avhrr.validate_lwp_util import LWP_THRESHOLD
import logging
logger = logging.getLogger(__name__)
AMSR_RADIUS = 10e3
def get_amsr(filename):
    retv = AmsrObject()
    with h5py.File(filename, 'r') as f:
        #ravel AMSR-E data to 1 dimension

        lon = f['Swath1/Geolocation Fields/Longitude'][:] 
        lat = f['Swath1/Geolocation Fields/Latitude'][:]
        n_lat_scans = lat.shape[1] #= 242!
        sec1993 = f['Swath1/Geolocation Fields/Time']['Time'][:]
        #if lon.shape[1]>0:
        #    print "Warning, expected 1D vectors"
        density = 1e3 # Density of water [kg m**-3]
        lwp_mm = f['Swath1/Data Fields/High_res_cloud'][:]
                            #description='lwp (mm)',
        lwp_gain = f['Swath1/Data Fields/High_res_cloud'].attrs['Scale']#.ravel()   
    epoch_diff = timegm(TAI93.utctimetuple())    
    nadir_sec_1970 = sec1993 + epoch_diff
    retv.sec_1970 = np.repeat(nadir_sec_1970.ravel(), n_lat_scans)
    retv.longitude = lon.ravel()
    retv.latitude = lat.ravel()
    retv.lwp = lwp_mm.ravel() *  lwp_gain * density # [mm * kg m**-3 = g m**-2]
    if f:
        f.close()   
    logger.info("Extract AMSR-E lwp between 0 and %d g/m-2", LWP_THRESHOLD)
    use_amsr = np.logical_and(retv.lwp >0 ,
                              retv.lwp < LWP_THRESHOLD)
    retv = calipso_track_from_matched(retv, retv, use_amsr)

    return retv 


def reshapeAmsr(amsrfiles, avhrr):
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    amsr = get_amsr(amsrfiles[0])
    for i in range(len(amsrfiles)-1):
        newAmsr = get_amsr(amsrfiles[i+1])
        amsr_start_all = amsr.sec_1970.ravel()
        amsr_new_all = newAmsr.sec_1970.ravel()
        if not amsr_start_all[0]<amsr_new_all[0]:
            raise ProcessingError("AMSR files are in the wrong order")
        amsr_break = np.argmin(np.abs(amsr_start_all - amsr_new_all[0]))+1
        # Concatenate the feature values
        #arname = array name from amsrObj
        for arname, value in amsr.all_arrays.items(): 
            if value is not None:
                if value.size != 1:
                    amsr.all_arrays[arname] = np.concatenate(
                        (amsr.all_arrays[arname],
                         newAmsr.all_arrays[arname]),axis=0)          
    # Finds Break point
    #import pdb; pdb.set_trace()
    startBreak, endBreak = find_break_points(amsr, avhrr)
    amsr = time_reshape_calipso(amsr, startBreak, endBreak)
    return amsr



def match_amsr_avhrr(amsrObj, imagerGeoObj, imagerObj, ctype, cma, ctth, nwp,
                     imagerAngObj, cpp, nwp_segments):
    retv = AmsrAvhrrTrackObject()
    from common import map_avhrr
    cal, cap = map_avhrr(imagerGeoObj, 
                         amsrObj.longitude.ravel(), 
                         amsrObj.latitude.ravel(),
                         radius_of_influence=AMSR_RADIUS,
                         n_neighbours=8)
    cal_1 = cal[:,0]
    cap_1 = cap[:,0]

    calnan = np.where(cal_1 == NODATA, np.nan, cal_1)
    if (~np.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal_1,cap_1)]
        imager_lines_sec_1970 = np.where(cal_1 != NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal_1 != NODATA, imagerGeoObj.time[cal_1], np.nan)
    # Find all matching Amsr pixels within +/- sec_timeThr from the AVHRR data
    idx_match = elements_within_range(amsrObj.sec_1970, imager_lines_sec_1970, sec_timeThr)

    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)  

    retv.amsr = calipso_track_from_matched(retv.amsr, amsrObj, idx_match)
 
    # Amsr line,pixel inside AVHRR swath (one neighbour):
    retv.amsr.imager_linnum = np.repeat(cal_1, idx_match).astype('i')
    retv.amsr.imager_pixnum = np.repeat(cap_1, idx_match).astype('i')
    retv.amsr.imager_linnum_nneigh = np.repeat(cal, idx_match, axis=0)
    retv.amsr.imager_pixnum_nneigh = np.repeat(cap, idx_match, axis=0)

    # Imager time
    retv.avhrr.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.amsr.sec_1970 - retv.avhrr.sec_1970

    do_some_logging(retv, amsrObj)
    logger.debug("Extract imager lwp along track!")
    
    nwp_small = NWPObj({'fractionofland': getattr(nwp,'fractionofland') ,
                        'landuse': getattr(nwp, 'landuse')})
                    
    retv = avhrr_track_from_matched(retv, imagerGeoObj, None, imagerAngObj, 
                                    #nwp, ctth, ctype, cma,  
                                    nwp_small, None, None, cma,
                                    cpp=cpp, nwp_segments=None,
                                    extract_some_data_for_x_neighbours=True)
    return retv
