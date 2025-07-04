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
# Change log is found in git
from atrain_match.utils.runutils import do_some_logging
from atrain_match.truths.calipso import (find_break_points)
from atrain_match.libs.extract_imager_along_track import imager_track_from_matched
from atrain_match.utils.common import (MatchupError, ProcessingError,
                                       elements_within_range)
import atrain_match.config as config
from atrain_match.matchobject_io import (EarthCareObject,
                                         TruthImagerTrackObject)
import time
import numpy as np
import os
import logging
import netCDF4
logger = logging.getLogger(__name__)

def read_earthcare(filename):
    retv = EarthCareObject()
    ncf = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    retv.latitude = ncf["ScienceData"]["latitude"][:]
    retv.longitude = ncf["ScienceData"]["longitude"][:]
    validation_height = ncf["ScienceData"]["ATLID_cloud_top_height"][:]
    validation_height[validation_height > 99999] = 1 
    retv.cloud_fraction = np.where(validation_height > 0, 1, 0)
    retv.validation_height = validation_height
    import datetime
    ep = datetime.datetime(1970,1,1,0,0,0)
    ep2 = datetime.datetime(2000,1,1,0,0,0)
    delta_seconds = ep2 - ep
    retv.sec_1970 = ncf["ScienceData"]["time"][:] + delta_seconds.total_seconds()
    return retv
    
def get_earthcare(filename):
    # Read EARTHCARE Radar data and add some variables
    earthcare = read_earthcare(filename)
    return earthcare

def match_earthcare_imager(earthcare, cloudproducts, SETTINGS):
    retv = TruthImagerTrackObject(truth='earthcare')
    retv.imager_instrument = cloudproducts.instrument.lower()
    retv.earthcare = earthcare
    # Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    from atrain_match.utils.common import map_imager
    cal, cap = map_imager(cloudproducts,
                          earthcare.longitude.ravel(),
                          earthcare.latitude.ravel(),
                          radius_of_influence=config.RESOLUTION * 0.7 * 1000.0)  # somewhat larger than radius...
    calnan = np.where(cal == config.NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None
    # check if it is within time limits:
    if len(cloudproducts.time.shape) > 1:
        imager_time_vector = [cloudproducts.time[line, pixel] for line, pixel in zip(cal, cap)]
        imager_lines_sec_1970 = np.where(cal != config.NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != config.NODATA, cloudproducts.time[cal], np.nan)
    # Find all matching Earthcare pixels within +/- sec_timeThr from the IMAGER data
    idx_match = elements_within_range(earthcare.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr"])

    if idx_match.sum() == 0:
        logger.warning("No matches in region within time threshold %d s.", SETTINGS["sec_timeThr"])
        return None
    retv.earthcare = retv.earthcare.extract_elements(idx=idx_match)

    # Earthcare line, pixel inside IMAGER swath:
    retv.earthcare.imager_linnum = np.repeat(cal, idx_match).astype('i')
    retv.earthcare.imager_pixnum = np.repeat(cap, idx_match).astype('i')

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.earthcare.sec_1970 - retv.imager.sec_1970
    do_some_logging(retv, earthcare)
    logger.debug("Generate the latitude, cloudtype tracks!")
    retv = imager_track_from_matched(retv, SETTINGS,
                                     cloudproducts)
    return retv


def reshape_earthcare(earthcarefiles, imager, SETTINGS):
    # imager_end = imager.sec1970_end
    # imager_start = imager.sec1970_start
    clsat = get_earthcare(earthcarefiles[0])
    for i in range(len(earthcarefiles) - 1):
        newEarthcare = get_earthcare(earthcarefiles[i + 1])
        clsat_start_all = clsat.sec_1970.ravel()
        clsat_new_all = newEarthcare.sec_1970.ravel()
        if not clsat_start_all[0] < clsat_new_all[0]:
            raise ProcessingError("Earthcare files are in the wrong order!")
        # clsat_break = np.argmin(np.abs(clsat_start_all - clsat_new_all[0]))+1
        # Concatenate the feature values
        clsat = clsat + newEarthcare
        # print("taistart", clsat.TAI_start)

    # Finds Break point
    start_break, end_break = find_break_points(clsat, imager, SETTINGS)
    clsat = clsat.extract_elements(starti=start_break,
                                   endi=end_break)
    return clsat
