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
# -*- coding: utf-8 -*-
#
# test_match.py
"""for the match module match.py """
#
import numpy as np
import unittest
import match 

class test_match_lon_lat(unittest.TestCase): 

    def setUp(self):
        lat = np.array([[999, 10, 15],
                        [-999, 20, 25]])
        lon = np.array([[999, 10, 15],
                        [-999, 20, 25]])
        lat2 = np.array([10, 10, 15])
        lon2 = np.array([10, 10, 15])
        lat3 = np.array([[10, 10, 15],[10, 10, 15]])
        lon3 = np.array([[10, 11, 16],[12, 10, 15]])
        self.source = (lat, lon)
        self.source3 = (lat3, lon3)
        self.target = (lat2, lon2)
        self.RESOLUTION=5

    def test_match(self):
        mapper = match.match_lonlat(self.source3, self.target,
                                      radius_of_influence=0.7*5*1000.0, 
                                      n_neighbours=1)
        print mapper.rows
        print mapper.cols.data[1]
        self.assertEqual(mapper.cols.data[0], 0)
        self.assertEqual(mapper.cols.data[2], 2)

    def test_match_nodata(self):
         mapper = match.match_lonlat(self.source, self.target,
                                      radius_of_influence=0.7*5*1000.0, 
                                      n_neighbours=1)
         print mapper.rows
         print mapper.cols.data[1]
         self.assertEqual(mapper.cols.data[1], 1)


def suite():
    """The suite for test_utils.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(test_match_lon_lat))
    return mysuite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
