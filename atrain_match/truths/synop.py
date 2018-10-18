import numpy as np
import pandas as pd
import time
import calendar
from datetime import datetime, timedelta
from calendar import timegm
TAI93 = datetime(1993, 1, 1)
import  config
from matchobject_io import (SynopImagerTrackObject, 
                            SynopObject)
from truths.calipso import (find_break_points, calipso_track_from_matched,
                     do_some_logging)

from utils.common import (ProcessingError, MatchupError, elements_within_range)
from libs.extract_imager_along_track import imager_track_from_matched
import logging
logger = logging.getLogger(__name__)

def reshapeSynop(synopfiles, imager,  SETTINGS):
    start_t = datetime.utcfromtimestamp(imager.sec1970_start)
    end_t = datetime.utcfromtimestamp(imager.sec1970_end)
    #datetime.datetime.fromtimestamp(
    from truths.read_synop_dwd import get_synop_data
    items = [get_synop_data(filename) for filename in synopfiles]
    panda_synops = pd.concat(items, ignore_index=True)
    #import pdb
    #pdb.set_trace()
    dt_ = timedelta(seconds=SETTINGS["sec_timeThr_synop"])
    newsynops = panda_synops[panda_synops['date'] < end_t + dt_]
    panda_synops = newsynops[newsynops['date']  >  start_t - dt_ ]
    retv = SynopObject()
    retv.longitude = np.array(panda_synops['lon'])
    retv.latitude = np.array(panda_synops['lat'])
    retv.cloud_fraction = np.array(panda_synops['total_cloud_cover'])/8.0

#    clmask_obs = np.logical_and(
#        matchups['obs'] > 2. / 8., matchups['obs'] < 6. / 8.)
#    clmask_sat = np.logical_and(
#        matchups['sat'] > 5. / 16., matchups['sat'] < 10. / 16.)

    retv.sec_1970 = np.array([calendar.timegm(tobj.timetuple()) for tobj in panda_synops['date']])
    return retv

def match_synop_imager(synopObj, imagerGeoObj, imagerObj, ctype, cma, ctth, nwp,
                       imagerAngObj, cpp, nwp_segments, SETTINGS):
    retv = SynopImagerTrackObject()
    retv.imager_instrument = imagerGeoObj.instrument.lower()
    from utils.common import map_imager_distances
    n_neighbours = 250
    if config.RESOLUTION == 5:
        n_neighbours = 16
    mapper_and_dist = map_imager_distances(imagerGeoObj, 
                                          synopObj.longitude.ravel(), 
                                          synopObj.latitude.ravel(), 
                                          radius_of_influence=SETTINGS["SYNOP_RADIUS"], 
                                          n_neighbours=n_neighbours)
    #pdb.set_trace()
    cal, cap = mapper_and_dist["mapper"]
    distances = mapper_and_dist["distances"]
    cal_1 = cal[:,0]
    cap_1 = cap[:,0]

    calnan = np.where(cal_1 == config.NODATA, np.nan, cal_1)
    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None   
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal_1,cap_1)]
        imager_lines_sec_1970 = np.where(cal_1 != config.NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal_1 != config.NODATA, imagerGeoObj.time[cal_1], np.nan)
    idx_match = elements_within_range(synopObj.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr_synop"])
    if idx_match.sum() == 0:
        logger.warning("No  matches in region within time threshold %d s.", SETTINGS["sec_timeThr_synop"])
        return None
    retv.synop = calipso_track_from_matched(retv.synop, synopObj, idx_match)
 
    # Synop line,pixel inside IMAGER swath (one nearest neighbour):
    retv.synop.imager_linnum = np.repeat(cal_1, idx_match).astype('i')
    retv.synop.imager_pixnum = np.repeat(cap_1, idx_match).astype('i')
    # Synop line,pixel inside IMAGER swath (several neighbours):
    retv.synop.imager_linnum_nneigh = np.repeat(cal, idx_match, axis=0)
    retv.synop.imager_pixnum_nneigh = np.repeat(cap, idx_match, axis=0)
    retv.synop.imager_synop_dist = np.repeat(distances, idx_match, axis=0)

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.synop.sec_1970 - retv.imager.sec_1970

    do_some_logging(retv, synopObj)
    logger.debug("Extract imager along track!")
    
    retv = imager_track_from_matched(retv, SETTINGS,
                                     imagerGeoObj, imagerObj, imagerAngObj, 
                                     nwp, ctth, ctype, cma,  
                                     #nwp_small, ctth, ctype, cma,
                                     cpp=cpp, nwp_segments=None,
                                     find_mean_data_for_x_neighbours=True)
    return retv
