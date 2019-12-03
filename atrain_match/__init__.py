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
"""Package Initializer for level1c4pps."""

import numpy as np

from atrain_match.version import __version__  # noqa F401


def centered_modulus(array, divisor):
    """Transform array to half open range ]-divisor/2, divisor/2]."""
    arr = array % divisor
    arr[arr > divisor / 2] -= divisor
    return arr


def get_fibonacci_grid(radius_km, half_num_points=None):
    """Return lon lat for fibonacci grid with requested radius, or half_num_points.

    Method as in:
    # http://arxiv.org/pdf/0912.4540.pdf
    # Alvaro Gonzalez: Measurement of areas on sphere usig Fibonacci and latitude-longitude grid.

    Relation number of points and radius of the points:
    number of points    radius
        4000             200km
       16000             100km
       65000              50km
     1625000              10km

    """
    if half_num_points is None:
        EARTH_AREA = 510072000  # Earth area (km2)
        POINT_AREA = radius_km * radius_km * np.pi
        half_num_points = int(0.5 * EARTH_AREA / POINT_AREA)
    num_points = 2 * half_num_points + 1
    latitude = np.rad2deg(np.arcsin(np.linspace(-1.0, 1.0, num_points)))
    lin_space = np.linspace(-half_num_points, half_num_points, num_points)
    theta = (1 + np.sqrt(5)) * 0.5
    longitude = 360 * (lin_space % theta) / theta
    longitude = centered_modulus(longitude, 360)
    return longitude, latitude
