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
"""
Match TRUTH and IMAGER data

"""

# from __future__ import with_statement
import numpy as np
import logging
import h5py
from atrain_match.config import RESOLUTION, NODATA
logger = logging.getLogger(__name__)


class MatchMapper(object):
    """
    Map arrays from one swath to another.

    Note that MatchMapper always works with an extra dimension: neighbour

    """

    def __init__(self, rows, cols, pixel_mask, time_diff=None,
                 time_threshold=None):
        self._rows = np.array(rows).astype(np.int)
        self._cols = np.array(cols).astype(np.int)
        self._pixel_mask = pixel_mask
        self._time_diff = time_diff
        self.time_threshold = time_threshold

    def __call__(self, array):
        """
        Maps *array* to target swath.

        """
        return np.ma.array(array[self.rows, self.cols], mask=self.mask)

    @property
    def rows(self):
        # self._rows = np.array(self._rows, dtype=np.int64)
        return np.ma.array(self._rows, mask=self.mask, fill_value=NODATA,
                           hard_mask=True)

    @property
    def cols(self):

        # self._cols = np.array(self._cols, dtype=np.int64)
        return np.ma.array(self._cols, mask=self.mask, fill_value=-NODATA,
                           hard_mask=True)

    @property
    def time_diff(self):
        """Time difference in seconds"""
        if self._time_diff is None:
            return None
        # Only use pixel mask
        return np.ma.array(self._time_diff, mask=self._pixel_mask,
                           fill_value=np.inf, hard_mask=True)

    @time_diff.setter
    def time_diff(self, value):
        self._time_diff = value

    @property
    def mask(self):
        if not None in (self.time_diff, self.time_threshold):
            return (self._pixel_mask +
                    (abs(self.time_diff) > self.time_threshold))
        return self._pixel_mask


def match_lonlat(source, target,
                 radius_of_influence=0.7*RESOLUTION*1000.0,
                 n_neighbours=1):
    """
    Produce a masked array of the same shape as the arrays in *target*, with
    indices of nearest neighbours in *source*. *source* and *target* should be
    tuples (lon, lat) of the source and target swaths, respectively.

    Note::

        * Fastest matching is obtained when *target* has lower resolution than
        *source*.

        * *source* should have 2-dimensional lon and lat arrays.

    """
    from pyresample.geometry import SwathDefinition
    from pyresample.kd_tree import get_neighbour_info
    from pyresample.kd_tree import get_sample_from_neighbour_info

    lon, lat = source
    mask_out_lat = np.logical_or(lat < -90, lat > 90)
    mask_out_lon = np.logical_or(lon > 180, lon < -180)
    mask_out = np.logical_or(mask_out_lat, mask_out_lon)
    lat = np.ma.masked_array(lat, mask=mask_out)
    lon = np.ma.masked_array(lon, mask=mask_out)
    # lat = np.around(lat, decimals=4)
    # lon = np.around(lon, decimals=4)
    source_def = SwathDefinition(*(lon, lat))
    target_def = SwathDefinition(*target)
    logger.debug("Matching %d nearest neighbours", n_neighbours)
    valid_in, valid_out, indices, distances = get_neighbour_info(
        source_def, target_def, radius_of_influence, neighbours=n_neighbours)
    # import pdb; pdb.set_trace()
    # Use pyresampe code to find colmun and row numbers for each pixel
    # This is works also with no-data in imager lat/lon.
    cols_matrix, rows_matrix = np.meshgrid(np.array(range(0, lat.shape[1])),
                                           np.array(range(0, lat.shape[0])))
    if n_neighbours == 1:
        first_indices = indices
    else:
        first_indices = indices[:, 0]

    cols = get_sample_from_neighbour_info('nn', target_def.shape,
                                          cols_matrix,
                                          valid_in,
                                          valid_out,
                                          first_indices)
    rows = get_sample_from_neighbour_info('nn', target_def.shape,
                                          rows_matrix,
                                          valid_in, valid_out,
                                          first_indices)
    if n_neighbours > 1:
        rows_0 = rows.copy()
        cols_0 = cols.copy()
        rows = NODATA + np.zeros((len(rows_0), n_neighbours))
        cols = NODATA + np.zeros((len(cols_0), n_neighbours))
        rows[:, 0] = rows_0
        cols[:, 0] = cols_0
        for i in range(1, n_neighbours):
            cols[:, i] = get_sample_from_neighbour_info('nn', target_def.shape,
                                                        cols_matrix,
                                                        valid_in, valid_out,
                                                        indices[:, i])
            rows[:, i] = get_sample_from_neighbour_info('nn', target_def.shape,
                                                        rows_matrix,
                                                        valid_in, valid_out,
                                                        indices[:, i])
            test = (distances[:, 0] - distances[:, i])
            if sum(~np.isnan(test)) > 0 and np.max(test[~np.isnan(test)]) > 0:
                raise ValueError(
                    'We count on the first neighbour beeing the closest')

    rows = np.array(rows)
    cols = np.array(cols)
    # import pdb; pdb.set_trace()
    """ Code used during debugging, leaving it here for now
    if indices.dtype in ['uint32']:
        # With pykdtree installed get_neighbour_info returns indices
        # as type uint32
        # This does not combine well with a nodata value of -9.
        indices = np.array(indices, dtype=np.int64)
    """
    # Make sure all indices are valid
    # import pdb; pdb.set_trace()
    rows[rows >= source_def.shape[0]] = NODATA
    cols[cols >= source_def.shape[1]] = NODATA
    mask = np.logical_or(distances > radius_of_influence,
                         indices >= len(valid_in))
    distances[distances > radius_of_influence] = -9
    # import pdb; ipdb.set_trace()
    rows[mask] = NODATA
    cols[mask] = NODATA
    # import pdb; pdb.set_trace()
    return MatchMapper(rows, cols, mask), distances
