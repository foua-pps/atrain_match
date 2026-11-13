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
from dataclasses import dataclass
from scipy.interpolate import interp1d
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


@dataclass
class AcClpReader:
    """Reader for EarthCARE CPR + ATLID based cloud microphysics product."""

    dataset: netCDF4.Dataset

    @property
    def optical_thickness_1km(self):
        """Get 1km optical thickness."""

        return self.dataset["ScienceData/Data/optical_thickness_1km"][:]

    @property
    def optical_thickness_10km(self):
        """Get 10km optical thickness."""

        return self.dataset["ScienceData/Data/optical_thickness_10km"][:]

    @property
    def ice_water_path_1km(self):
        """Get 1km ice water path along track."""

        return self.dataset["ScienceData/Data/ice_water_path_1km"][:]

    @property
    def ice_water_path_10km(self):
        """Get 10km ice water path along track."""

        return self.dataset["ScienceData/Data/ice_water_path_10km"][:]

    @property
    def liquid_water_path_1km(self):
        """Get 1km liquid water path along track."""

        return self.dataset["ScienceData/Data/liquid_water_path_1km"][:]

    @property
    def liquid_water_path_10km(self):
        """Get 10km liquid water path along track."""

        return self.dataset["ScienceData/Data/liquid_water_path_10km"][:]


@dataclass
class CprTcReader:
    """Reader for EarthCARE Cloud Profiling Radar Target Classification (CPR-TC)"""

    dataset: netCDF4.Dataset

    @property
    def cloud_top(self):
        """Returns the maximum cloud top height profile-wise along the track."""

        return self.get_cloud_height("top")

    @property
    def cloud_base(self):
        """Returns the minimum cloud base height profile-wise along the track."""

        return self.get_cloud_height("base")

    @property
    def cloud_fraction(self):
        """
        Extract a simple binary cloud fraction (0/1) from target classifictaion.
        Hydrometeor Classification, where:
        -1 = no data,
        0 = ground and underground levels,
        1 = clear,
        2 = no drizzling liquid cloud,
        3 = drizzling liquid clouds,
        4 = drizzle or warm rain,
        5 = cold rain,
        6 = melting snow,
        7 = rimed snow,
        8 = snow,
        9 = ice,
        10 = stratospheric cloud (ice),
        11 = insects (or artifacts),
        12 = heavy rain likely present,
        13 = heavy mixed-phase precipitation likely present,
        14 = heavy rain,
        15 = heavy mixed-phase precipitation,
        16 = removed ground clutter, possible rain in clutter,
        17 = removed ground clutter, possible snow in clutter,
        18 = removed ground clutter, possible cloud in clutter,
        19 = removed ground clutter, cloud/precip unlikely present in this layer,
        20 = uncertain
        """

        hyd_class = self.dataset["ScienceData"]["hydrometeor_classification"][:]
        cloud_codes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15]

        is_clear = hyd_class == 1
        is_cloud = np.isin(hyd_class, cloud_codes)

        has_cloud = is_cloud.any(axis=1)
        has_clear = is_clear.all(axis=1)

        along_track_size = hyd_class.shape[0]
        cloud_fraction = np.full(along_track_size, -9, dtype=int)

        cloud_fraction[has_clear] = 0
        cloud_fraction[has_cloud] = 1

        return cloud_fraction

    def get_cloud_height(self, mode="base"):
        """Extract profile-wise cloud base or top heights."""

        if mode == "base":
            height_wgs84 = (self.dataset["ScienceData"]["cloud_base"][:]).min(axis=1)
        elif mode == "top":
            height_wgs84 = (self.dataset["ScienceData"]["cloud_top"][:]).max(axis=1)
        else:
            raise ValueError("mode must be either 'base' or 'top'")

        # above sea level = geodetic_altitude - geoid_offset
        geoid_offset = self.dataset["ScienceData"]["geoid_offset"][:]
        height_asl = height_wgs84 - geoid_offset
        height_asl.set_fill_value(-9)

        return height_asl

    def get_cloud_pressure_or_temperature(self, var_name, mode="base"):
        """Interpolates height to pressure or temperature, from profiles."""

        track_heights = self.dataset["ScienceData"]["height"][:]
        track_values = self.dataset["ScienceData"][var_name][:]
        n_track = track_heights.shape[0]
        if mode == "base":
            target_heights = self.cloud_base
        elif mode == "top":
            target_heights = self.cloud_top
        else:
            raise ValueError("mode must be either 'base' or 'top'")

        var_interpolated = np.full(n_track, -9, dtype=float)

        for i in range(n_track):
            if target_heights[i] > 0:
                var_interpolated[i] = _interpolate_single_profile(
                    track_heights[i],
                    track_values[i],
                    target_heights[i],
                )

        return var_interpolated


def _interpolate_single_profile(heights, values, target_height):
    """Helper function to interpolate from profile."""

    sort_idx = np.argsort(heights)

    return np.interp(
        target_height, heights[sort_idx], values[sort_idx], left=-9, right=-9
    )


@dataclass
class AtlCthReader:
    dataset: netCDF4.Dataset

    @property
    def cloud_top(self):
        """Returns quality controlled cloud top height from ATLID."""

        return self.get_cloud_top_height()

    @property
    def thick_height(self):
        """Returns height of thick clouds from ATLID."""

        thick_height = self.dataset["ScienceData"]["ATLID_thick_cloud_top_height"][:]
        thick_height[thick_height > 99999] = -9
        return thick_height

    @property
    def quality_status(self):
        """Returns quality flag."""

        return self.dataset["ScienceData"]["quality_status"][:]

    @property
    def cloud_fraction(self):
        """Returns simplified cloud fraction from ATLID.
        # 0 – no cloud,
        # 1 – thick cloud,
        # 2 – thin cloud,
        # 3 – thin over thick cloud,
        # 4 – thick over thick cloud,
        # 5 – thin over thin cloud
        # 6 – no cloud, but probably cloud influenced
        # import pdb;pdb.set_trace()
        """

        simplified_uppermost_cloud_classification = self.dataset["ScienceData"][
            "simplified_uppermost_cloud_classification"
        ][:]

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

    def get_cloud_top_height(self):
        """Extract the cloud top height with quality control."""

        height_wgs84 = self.dataset["ScienceData"]["ATLID_cloud_top_height"][:]
        geoid_offset = self.dataset["ScienceData"]["geoid_offset"][:]

        validation_height = height_wgs84 - geoid_offset
        validation_height[validation_height > 99999] = -9

        validation_height[self.quality_status == 3] = -9
        validation_height[self.quality_status == 2] = -9
        validation_height[self.quality_status == 4] = -9

        # update = np.logical_and(validation_height == -9, self.thick_height > 0)
        # validation_height[update] = thick_height[update]
        # thick_validation_height = self.dataset["ScienceData"]["ATLID_thick_cloud_top_height"][:]
        # thick_validation_height[thick_validation_height >9999 ] = -9

        return validation_height


def read_earthcare(filename):
    """
    Read earthcare data based on the input file.
    For each product type supported, exists a dataclass to read in the data.
    To add new products, simply make a new dataclass.
    """

    product_type = EarthCareProduct.from_path(filename)

    retv = EarthCareObject()
    ncf = netCDF4.Dataset(filename, "r", format="NETCDF4")

    # Helper function to handle variable location differences across file versions
    def _get_var(name):
        try:
            return ncf["ScienceData"]["Geo"][name][:]
        except Exception as e:
            return ncf["ScienceData"][name][:]

    # Fetch navigation and time variables dynamically
    retv.latitude = _get_var("latitude")
    retv.longitude = _get_var("longitude")
    retv.sec_1970 = get_timestamps(_get_var("time"))

    if product_type.value == "ATL_CTH":
        # Cloud height only from ATLID
        logger.debug("Reading Earthcare product ATL_CTH.")
        atl_cth = AtlCthReader(ncf)

        retv.validation_height = atl_cth.cloud_top
        retv.thick_height = atl_cth.thick_height
        retv.quality_status = atl_cth.quality_status
        retv.cloud_fraction = atl_cth.cloud_fraction

    elif product_type.value == "CPR_TC":
        # cloud top and base height from CPR
        logger.debug("Reading Earthcare product CPR_TC.")
        cpr_tc = CprTcReader(ncf)

        retv.validation_height_base = cpr_tc.cloud_base
        retv.validation_height = cpr_tc.cloud_top
        retv.cloud_fraction = cpr_tc.cloud_fraction
        retv.layer_base_pressure = cpr_tc.get_cloud_pressure_or_temperature(
            "pressure", "base"
        )
        retv.layer_base_temperature = cpr_tc.get_cloud_pressure_or_temperature(
            "temperature", "base"
        )
        retv.layer_top_pressure = cpr_tc.get_cloud_pressure_or_temperature(
            "pressure", "top"
        )
        retv.layer_top_temperature = cpr_tc.get_cloud_pressure_or_temperature(
            "temperature", "top"
        )

    elif product_type.value == "AC_CLP":
        # LWP, IWP, Cloud mask from CPR + ATLID
        logger.debug("Reading Earthcare product AC_CLP.")
        ac_clp = AcClpReader(ncf)

        retv.ice_water_path_1km = ac_clp.ice_water_path_1km
        retv.liquid_water_path_1km = ac_clp.liquid_water_path_1km
        retv.optical_thickness_1km = ac_clp.optical_thickness_1km

        # similar can be added for 10km, if required!

    else:
        logger.warning(
            "No Earthcare product identified! The process will proceed with None."
        )

    ncf.close()

    return retv


def get_timestamps(raw_time):
    """Convert satellite timestamps from the 2000-01-01 epoch to the 1970-01-01 Unix epoch."""

    ep = datetime.datetime(1970, 1, 1, 0, 0, 0)
    ep2 = datetime.datetime(2000, 1, 1, 0, 0, 0)
    delta_seconds = ep2 - ep

    return raw_time.data + delta_seconds.total_seconds()


def get_earthcare(filename):
    # Read EARTHCARE Radar data and add some variables

    try:
        earthcare = read_earthcare(filename)
    except Exception as e:
        logger.error(f"Something went wrong, check {e}")
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

    # Finds Break point
    start_break, end_break = find_break_points(clsat, imager, SETTINGS)
    clsat = clsat.extract_elements(starti=start_break, endi=end_break)
    return clsat
