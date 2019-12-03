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

"""Test function for the atrain_match.__init__ file."""

from atrain_match import get_fibonacci_grid, centered_modulus
import unittest
import numpy as np
from atrain_match.reshaped_files_scr.plot_kuipers_on_area_util import get_fibonacci_spread_points_on_earth


class TestInit(unittest.TestCase):
    """Test function for the angle calculation."""

    def test_get_fibonacci_grid_compare_to_old_method(self):
        """Test get_fibonacci_grid compared to old code."""
        lons, lats = get_fibonacci_grid(200, half_num_points=2010)
        lonsold, latsold = get_fibonacci_spread_points_on_earth(200, num_points=4021)
        np.testing.assert_allclose(lats, latsold, atol=3.0)
        np.testing.assert_allclose(lons[:-1], lonsold[1:], atol=3.0)

    def test_get_fibonacci_grid(self):
        """Test get_fibonacci_grid."""
        expected_lons = np.array([-32.46117975, -169.9689438, 52.52329215, -84.9844719,
                                  137.50776405, 0, -137.50776405, 84.9844719,
                                  -52.52329215, 169.9689438, 32.46117975])
        expected_lats = np.array([-90, -53.13010235, -36.86989765, -23.57817848,
                                  -11.53695903, 0, 11.53695903, 23.57817848,
                                  36.86989765, 53.13010235, 90.])
        lons, lats = get_fibonacci_grid(200, half_num_points=5)
        np.testing.assert_almost_equal(lats, expected_lats, decimal=3)
        np.testing.assert_almost_equal(lons, expected_lons, decimal=3)
        self.assertTrue(get_fibonacci_grid(200)[0].shape[0] - 4000 < 200)
        self.assertTrue(get_fibonacci_grid(100)[0].shape[0] - 16000 < 500)
        self.assertEqual(get_fibonacci_grid(200, half_num_points=1000)[0].shape[0], 2001)

    def test_centered_modulus(self):
        """Test centered_modulus."""
        angles = np.ma.array(
            [[180.0,  -180.0, -179.9,  201.0],
             [80.0,  360.0, -360.0, 604.0],
             [-80.0,  -88.0, -796.0, -104.0],
             [-3.0,  -188.0, -196.0, -204.0]], mask=False)
        expected = np.ma.array(
            [[180.0,  180.0, -179.9,  -159.0],
             [80.0,  0.0, 0.0, -116.0],
             [-80.0,  -88.0, -76.0, -104.0],
             [-3.0,  172.0, 164.0, 156.0]], mask=False)
        transformed = centered_modulus(angles, 360.0)
        np.testing.assert_allclose(transformed, expected)


def suite():
    """Test non angle functions in pygac init file."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestInit))

    return mysuite


if __name__ == '__main__':
    unittest.main()
