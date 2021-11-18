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
from atrain_match.matchobject_io import (DardarObject,
                                         TruthImagerTrackObject)
import time
import datetime as dt
import numpy as np
import os
import logging
from netCDF4 import Dataset
logger = logging.getLogger(__name__)


def get_dardar(filename):
    if filename.endswith('.nc'):
        dardar = read_dardar_nc(filename)
    return dardar


def read_dardar_nc(filename):

    retv = DardarObject()
    ds = Dataset(filename, 'r')
    for v in ds.variables:
        if v in retv.all_arrays.keys():
            #retv.all_arrays[v] = ds[v][:]
            setattr(retv, v, ds[v][:])

    time_unit = ds['time'].units.split(' ')  # e.g. seconds since 2019 02 03 00:00:00 UTC
    year = int(time_unit[2])
    month = int(time_unit[3])
    day = int(time_unit[4])
    hour = int(time_unit[5].split(':')[0])
    min = int(time_unit[5].split(':')[1])
    sec = int(time_unit[5].split(':')[2])

    # set_auto_mask set to False because time valid limit in DARDAR files is too low (8000)
    ds['time'].set_auto_mask(False)
    time = ds['time'][:]

    # number of seconds between start of day of aquesition time and 1970/01/01 00:00:00
    dsec = dt.datetime(year, month, day, hour, min, sec) - dt.datetime(1970, 1, 1, 0, 0, 0)
    dsec = dsec.total_seconds()
    # seconds since 1970 are seconds between start of day of aquesition and 1970
    # + seconds between start of day of aquesition and measurement
    retv.sec_1970 = time + dsec

    ds.close()
    return retv


def reshape_dardar(dardarfiles, imager, SETTINGS):
    dardar = get_dardar(dardarfiles[0])
    for i in range(len(dardarfiles) - 1):
        print(dardarfiles[i+1])
        new_dardar = get_dardar(dardarfiles[i + 1])
        dardar_start_all = dardar.sec_1970.ravel()
        dardar_new_all = new_dardar.sec_1970.ravel()
        if not dardar_start_all[0] < dardar_new_all[0]:
            raise ProcessingError("DARDAR files are in the wrong order!")
        dardar = dardar + new_dardar

    # Finds Break point
    start_break, end_break = find_break_points(dardar, imager, SETTINGS)
    dardar = dardar.extract_elements(starti=start_break,
                                    endi=end_break)
    return dardar


def match_dardar_imager(dardar, cloudproducts, SETTINGS):
    retv = TruthImagerTrackObject(truth='dardar')
    retv.imager_instrument = cloudproducts.instrument.lower()
    retv.dardar = dardar
    # Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    from atrain_match.utils.common import map_imager
    cal, cap = map_imager(cloudproducts,
                          dardar.longitude.ravel(),
                          dardar.latitude.ravel(),
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
    # Find all matching Cloudsat pixels within +/- sec_timeThr from the IMAGER data
    idx_match = elements_within_range(dardar.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr"])

    if idx_match.sum() == 0:
        logger.warning("No matches in region within time threshold %d s.", SETTINGS["sec_timeThr"])
        return None

    print(idx_match.shape)
    retv.dardar = retv.dardar.extract_elements(idx=idx_match)

    # Cloudsat line, pixel inside IMAGER swath:
    retv.dardar.imager_linnum = np.repeat(cal, idx_match).astype('i')
    retv.dardar.imager_pixnum = np.repeat(cap, idx_match).astype('i')

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.dardar.sec_1970 - retv.imager.sec_1970
    do_some_logging(retv, dardar)
    logger.debug("Generate the latitude, cloudtype tracks!")
    retv = imager_track_from_matched(retv, SETTINGS,
                                     cloudproducts)
    return retv

if __name__ == '__main__':
    get_dardar('/cmsaf/cmsaf-cld8/EXTERNAL_DATA/DARDAR/DARDAR-CLOUD.v3.00/2019/02/03/DARDAR-CLOUD_2019034002856_68008_V3-00.nc')
