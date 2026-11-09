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
from enum import Enum
from pathlib import Path
from atrain_match.utils.runutils import do_some_logging
from atrain_match.truths.calipso import find_break_points
from atrain_match.libs.extract_imager_along_track import imager_track_from_matched
from atrain_match.utils.common import (
    MatchupError,
    ProcessingError,
    elements_within_range,
)
import atrain_match.config as config
from atrain_match.matchobject_io import EarthCareObject, TruthImagerTrackObject
import time
import numpy as np

import datetime

import os
import logging
import netCDF4

logger = logging.getLogger(__name__)


class EarthCareProduct(Enum):
    TYPE_1 = "ATL_CTH"
    TYPE_2 = "CPR_TC"
    TYPE_3 = "AC_CLP"

    @classmethod
    def from_path(cls, filename) -> "EarthCareProduct":
        name = Path(filename).name

        if "AC__CLP_2B" in name:
            return cls.TYPE_3

        if "CPR_TC__2A" in name:
            return cls.TYPE_2

        return cls.TYPE_1


def read_earthcare(filename):
    """Read earthcare data based on the input file."""

    product_type = EarthCareProduct.from_path(filename)

    retv = EarthCareObject()
    ncf = netCDF4.Dataset(filename, "r", format="NETCDF4")

    retv.latitude = ncf["ScienceData"]["latitude"][:]
    retv.longitude = ncf["ScienceData"]["longitude"][:]

    retv.sec_1970 = get_timestamps(ncf["ScienceData"]["time"][:])

    if product_type == "ATL_CTH":
        # Cloud height only from ATLID

        atlid_data = get_validation_height_atlid(ncf)
        retv.validation_height = atlid_data["validation_height"]
        retv.thick_height = atlid_data["thick_height"]
        retv.quality_status = atlid_data["quality_status"]
        retv.cloud_fraction = atlid_data["cloud_fraction"]

        retv.cloud_fraction = get_earthcare_classification(
            ncf["ScienceData"]["simplified_uppermost_cloud_classification"][:]
        )

    if product_type == "CPR_TC":
        # cloud top and base height from CPR

        retv.validation_height_base = get_cloud_height(
            ncf["ScienceData"]["cloud_base"][:], mode="base"
        )
        retv.validation_height = get_cloud_height(
            ncf["ScienceData"]["cloud_top"][:], mode="top"
        )

    ncf.close()

    return retv


def get_timestamps(raw_time):
    """Convert satellite timestamps from the 2000-01-01 epoch to the 1970-01-01 Unix epoch."""

    ep = datetime.datetime(1970, 1, 1, 0, 0, 0)
    ep2 = datetime.datetime(2000, 1, 1, 0, 0, 0)
    delta_seconds = ep2 - ep

    return raw_time.data + delta_seconds.total_seconds()


def get_cloud_height(layer_data, mode="base"):
    """Extract profile-wise cloud base or top heights."""

    if mode == "base":
        result = layer_data.min(axis=1)
    elif mode == "top":
        result = layer_data.max(axis=1)
    else:
        raise ValueError("mode must be either 'base' or 'top'")

    result.set_fill_value(-9)

    return result


def get_validation_height_atlid(ncf):
    """Get validation height from ATLID."""

    validation_height = ncf["ScienceData"]["ATLID_cloud_top_height"][:]
    thick_height = ncf["ScienceData"]["ATLID_thick_cloud_top_height"][:]
    thick_height[thick_height > 99999] = -9

    validation_height[validation_height > 99999] = -9

    quality_status = ncf["ScienceData"]["quality_status"][:]

    validation_height[quality_status == 3] = -9
    validation_height[quality_status == 2] = -9
    validation_height[quality_status == 4] = -9
    update = np.logical_and(validation_height == -9, thick_height > 0)
    # validation_height[update] = thick_height[update]
    # thick_validation_height = ncf["ScienceData"]["ATLID_thick_cloud_top_height"][:]
    # thick_validation_height[thick_validation_height >9999 ] = -9

    return {
        "validation_height": validation_height,
        "thick_height": thick_height,
        "quality_status": quality_status,
    }


def get_earthcare_classification(simplified_uppermost_cloud_classification):
    """Earthcare cloud type."""

    # 0 – no cloud,
    # 1 – thick cloud,
    # 2 – thin cloud,
    # 3 – thin over thick cloud,
    # 4 – thick over thick cloud,
    # 5 – thin over thin cloud
    # 6 – no cloud, but probably cloud influenced
    # import pdb;pdb.set_trace()
    cloud_fraction = np.where(
        simplified_uppermost_cloud_classification.data == 0, 0.0, 1.0
    )
    cloud_fraction[simplified_uppermost_cloud_classification.data == 6] = 0.7
    cloud_fraction[simplified_uppermost_cloud_classification.data == 2] = 1.0
    cloud_fraction[simplified_uppermost_cloud_classification.data == 5] = 1.0
    cloud_fraction[simplified_uppermost_cloud_classification.data == 3] = 1.0
    cloud_fraction[simplified_uppermost_cloud_classification.data == 4] = 1.0
    cloud_fraction[simplified_uppermost_cloud_classification.data == 1] = 1.0
    cloud_fraction[simplified_uppermost_cloud_classification.mask] = -9
    return cloud_fraction


def get_earthcare(filename):
    # Read EARTHCARE Radar data and add some variables

    try:
        earthcare = read_earthcare(filename)
    except:
        return None

    return earthcare


def match_earthcare_imager(earthcare, cloudproducts, SETTINGS):
    retv = TruthImagerTrackObject(truth="earthcare")
    retv.imager_instrument = cloudproducts.instrument.lower()
    retv.earthcare = earthcare
    # Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    from atrain_match.utils.common import map_imager

    cal, cap = map_imager(
        cloudproducts,
        earthcare.longitude.ravel(),
        earthcare.latitude.ravel(),
        radius_of_influence=config.RESOLUTION * 0.7 * 1000.0,
    )  # somewhat larger than radius...
    calnan = np.where(cal == config.NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None
    # check if it is within time limits:
    if len(cloudproducts.time.shape) > 1:
        imager_time_vector = [
            cloudproducts.time[line, pixel] for line, pixel in zip(cal, cap)
        ]
        imager_lines_sec_1970 = np.where(
            cal != config.NODATA, imager_time_vector, np.nan
        )
    else:
        imager_lines_sec_1970 = np.where(
            cal != config.NODATA, cloudproducts.time[cal], np.nan
        )
    # Find all matching Earthcare pixels within +/- sec_timeThr from the IMAGER data
    idx_match = elements_within_range(
        earthcare.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr"]
    )

    if idx_match.sum() == 0:
        logger.warning(
            "No matches in region within time threshold %d s.", SETTINGS["sec_timeThr"]
        )
        return None
    retv.earthcare = retv.earthcare.extract_elements(idx=idx_match)

    # Earthcare line, pixel inside IMAGER swath:
    retv.earthcare.imager_linnum = np.repeat(cal, idx_match).astype("i")
    retv.earthcare.imager_pixnum = np.repeat(cap, idx_match).astype("i")

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.earthcare.sec_1970 - retv.imager.sec_1970
    do_some_logging(retv, earthcare)
    logger.debug("Generate the latitude, cloudtype tracks!")
    retv = imager_track_from_matched(retv, SETTINGS, cloudproducts)
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
    clsat = clsat.extract_elements(starti=start_break, endi=end_break)
    return clsat
