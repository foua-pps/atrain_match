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
import numpy as np
import pandas as pd
import time
import calendar
from datetime import datetime, timedelta
from calendar import timegm
TAI93 = datetime(1993, 1, 1)
from matchobject_io import (TruthImagerTrackObject, 
                            MoraObject)
from utils.runutils import do_some_logging
import config
from utils.common import (ProcessingError, MatchupError, elements_within_range)
from libs.extract_imager_along_track import imager_track_from_matched
import logging
logger = logging.getLogger(__name__)

TEST_FILE = "/DATA_MISC/atrain_match_testcases/mora/cb_2010.dat"


def get_mora_data(filename):

    convert_datefunc = lambda x: datetime.strptime(x.decode("utf-8"), '%Y%m%dT%H%M%S')
    #convert_datefunc = lambda x: x.decode("utf-8")
    
    dtype = [('station', '|S5'),
             ('lat', 'f8'), 
             ('lon', 'f8'),
             ('x', 'f8'), 
             ('date', object), 
             ('cloud_base_height', 'i4')]

    data = np.genfromtxt(filename,
                         skip_header=0,
                         skip_footer=0,
                         usecols=(
                             0, 1, 2, 3, 4, 5),
                         dtype=dtype,
                         unpack=True,
                         converters={
                                     #2: lambda x: float(x) / 100.,
                                     #3: lambda x: float(x) / 100.,
                             4: convert_datefunc,
                                     #6: lambda x: float(x) / 10.,
                                     #7: lambda x: float(x) / 10.,
                                     #8: lambda x: float(x) / 10., 
                         })
    
    return pd.DataFrame(data)
    
def reshapeMora(morafiles, imager,  SETTINGS):
    start_t = datetime.utcfromtimestamp(imager.sec1970_start)
    end_t = datetime.utcfromtimestamp(imager.sec1970_end)
    #datetime.datetime.fromtimestamp(
    items = [get_mora_data(filename) for filename in morafiles]
    panda_moras = pd.concat(items, ignore_index=True)
    #import pdb
    #pdb.set_trace()
    dt_ = timedelta(seconds=SETTINGS["sec_timeThr_synop"])
    newmoras = panda_moras[panda_moras['date'] < end_t + dt_]
    panda_moras = newmoras[newmoras['date']  >  start_t - dt_ ]
    retv = MoraObject()
    retv.longitude = np.array(panda_moras['lon'])
    retv.latitude = np.array(panda_moras['lat'])
    retv.cloud_base_height = np.array(panda_moras['cloud_base_height'])
    retv.sec_1970 = np.array([calendar.timegm(tobj.timetuple()) for tobj in panda_moras['date']])
    return retv

def match_mora_imager(moraObj, cloudproducts, SETTINGS):
    retv = TruthImagerTrackObject('mora')
    retv.imager_instrument = cloudproducts.instrument.lower()
    retv.mora = moraObj

    from utils.common import map_imager

    """
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    from plotting.histogram_plotting import distribution_map
    plt.plot(cloudproducts.longitude.ravel(),
              cloudproducts.latitude.ravel(), 'r.')
    plt.plot(moraObj.longitude.ravel(), moraObj.latitude.ravel(),'b.')
    plt.show()
    """
    cal, cap = map_imager(cloudproducts, 
                         moraObj.longitude.ravel(),
                         moraObj.latitude.ravel(),
                          radius_of_influence=config.RESOLUTION*0.7*1000.0)
    calnan = np.where(cal == config.NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None   
    #check if it is within time limits:
    if len(cloudproducts.time.shape)>1:
        imager_time_vector = [cloudproducts.time[line,pixel] for line, pixel in zip(cal,cap)]
        imager_lines_sec_1970 = np.where(cal != config.NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != config.NODATA, cloudproducts.time[cal], np.nan)
    idx_match = elements_within_range(moraObj.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr_synop"])
    if idx_match.sum() == 0:
        logger.warning("No matches in region within time threshold %d s.", SETTINGS["sec_timeThr_synop"])
        return None
    retv.mora = retv.mora.extract_elements(idx=idx_match)
    # Mora line,pixel inside IMAGER swath (one nearest neighbour):
    retv.mora.imager_linnum = np.repeat(cal, idx_match).astype('i')
    retv.mora.imager_pixnum = np.repeat(cap, idx_match).astype('i')
    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.mora.sec_1970 - retv.imager.sec_1970

    do_some_logging(retv, moraObj)
    logger.debug("Extract imager along track!")
    
    retv = imager_track_from_matched(retv, 
                                     SETTINGS,
                                     cloudproducts,
                                     extract_nwp_segments = False)
    return retv

if __name__ == "__main__":
    get_mora_data(TEST_FILE)
