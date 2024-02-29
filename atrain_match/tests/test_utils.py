#!/usr/bin/env python
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
#

"""Test utils functions."""

import numpy as np
import unittest
from atrain_match.utils.pps_prototyping_util import (get_warmest_or_coldest_index)
from atrain_match.utils import match


def get_warmest_index_old(t11, matched):
    from scipy.ndimage.filters import generic_filter
    row, col = np.indices(t11.shape)
    flat_index = generic_filter(t11,
                                function=np.argmax,
                                size=5,
                                mode='constant',
                                cval=-9999999999999)
    flat_index = np.array(flat_index, dtype=np.int64)
    delta_row, delta_col = np.unravel_index(flat_index, (5, 5))
    delta_row = delta_row - 2
    delta_col = delta_col - 2
    new_row = row + delta_row
    new_col = col + delta_col
    new_row_matched = np.array([new_row[matched['row'][idx], matched['col'][idx]]
                                for idx in range(matched['row'].shape[0])])
    new_col_matched = np.array([new_col[matched['row'][idx], matched['col'][idx]]
                                for idx in range(matched['row'].shape[0])])
    return new_row_matched, new_col_matched


class test_prototyping_utils(unittest.TestCase):

    def setUp(self):

        self.t11 = np.array([[1, 43, 3, 4, 5, 6, 7, 8, 9, 35],
                             [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                             [1, 2, 3, 4, 25, 6, 7, 8, 9, 10],
                             [1, 2, 3, 4, 5, 6, 7, 8, 9, 101],
                             [11, 2, 3, 4, 5, 6, 7, 8, 9, 101]])

        self.matched = {'row': np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4]),
                        'col': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])}

    def test_warmest(self):
        retv = get_warmest_or_coldest_index(self.t11, self.matched)
        out_r = retv['row']
        out_c = retv['col']
        out_r2, out_c2 = get_warmest_index_old(self.t11, self.matched)

        self.assertTrue((out_r[0:4] == [0, 0, 0, 0]).all())
        self.assertTrue((out_c[0:4] == [1, 1, 1, 1]).all())
        self.assertTrue((out_r[4:7] == [2, 2, 2]).all())
        self.assertTrue((out_c[4:7] == [4, 4, 4]).all())
        self.assertTrue((out_c[7:10] == [9, 9, 9]).all())
        self.assertTrue((out_c == out_c2).all())
        self.assertTrue((out_r == out_r2).all())

class test_match_lon_lat(unittest.TestCase):

    def setUp(self):
        lon = np.array([[999, 10, 15],
                        [-999, 20, 25]]).astype(np.float64)
        lat = np.array([[999, 10, 15],
                        [-999, 20, 25]]).astype(np.float64)
        lon_t = np.array([31, 10, 15]).astype(np.float64)
        lat_t = np.array([89, 10, 15]).astype(np.float64)
        lon3 = np.array([[25, 10, 15], 
                         [10, 10, 15]]).astype(np.float64)
        lat3 = np.array([[10, 11, 16], 
                         [12, 10, 15]]).astype(np.float64)
        self.source = (lon, lat)
        self.source3 = (lon3, lat3)
        self.target = (lon_t, lat_t)
        self.RESOLUTION = 5

    def test_match_integers(self):
        from pyresample.geometry import SwathDefinition
        from pyresample.kd_tree import get_neighbour_info
        lon = np.array([0, 10, 25])
        source_def  = SwathDefinition(*(lon, lon))
        target_def = SwathDefinition(*(lon, lon))
        valid_in, valid_out, indices_int, distances = get_neighbour_info(source_def, target_def, 1000, neighbours=1)
        lon = np.array([0, 10, 25]).astype(np.float64)
        source_def  = SwathDefinition(*(lon, lon))
        target_def = SwathDefinition(*(lon, lon))
        valid_in, valid_out, indices_float, distances = get_neighbour_info(source_def, target_def, 1000, neighbours=1)


    def test_match(self):
        mapper, _ = match.match_lonlat(self.source3, self.target,
                                       radius_of_influence=0.7 * 5 * 1000.0,
                                       n_neighbours=1)
  
        self.assertEqual(mapper.cols.data[0], -9)
        self.assertEqual(mapper.cols.data[2], 2)

    def test_match_nodata(self):
        mapper, _ = match.match_lonlat(self.source, self.target,
                                       radius_of_influence=0.7 * 5 * 1000.0,
                                       n_neighbours=1)
        self.assertEqual(mapper.cols.data[1], 1)

def suite():
    """Create the suite for test_utils."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(test_prototyping_utils))
    mysuite.addTest(loader.loadTestsFromTestCase(test_match_lon_lat))
    return mysuite


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
