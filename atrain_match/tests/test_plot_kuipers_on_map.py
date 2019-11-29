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

"""Tests for the calculation measurements on fibonacci grid."""


import numpy as np
import unittest
from atrain_match.reshaped_files_scr.plot_kuipers_on_area_util import StatsOnFibonacciLattice


class test_kuipers_plot_on_map(unittest.TestCase):
    """Test the calcualtions made on fibbonaci grid."""

    def setUp(self):
        """Set up a finbonacci gird with some data."""
        self.lattice = StatsOnFibonacciLattice()
        self.lattice.num_detected_clouds = np.array([10, 1000, 1000, 1000, 1000, 1000])
        self.lattice.num_undetected_clouds = np.array([1000, 100, 100, 100, 100, 100])
        self.lattice.num_detected_clear = np.array([20, 80, 100, 0, 200, 100])
        self.lattice.num_false_clouds = np.array([80, 20, 1, 1, 200, 300])
        self.lattice.lats = np.array([80, 20, 1, 1, 200, 300])

        dummy = 0 * np.array([80, 20, 1, 1, 200, 300])
        self.lattice.sum_ctth_mae_low = dummy
        self.lattice.sum_ctth_mae_high = dummy
        self.lattice.sum_ctth_mae = dummy
        self.lattice.sum_ctth_mae_diff = dummy

        self.lattice.sum_ctth_bias_low = dummy
        self.lattice.sum_lapse_bias_low = dummy
        self.lattice.sum_ctth_bias_temperature_low = dummy
        self.lattice.sum_ctth_bias_temperature_low_t11 = dummy
        self.lattice.Min_lapse_rate = dummy
        self.lattice.num_new_detected_clouds = dummy
        self.lattice.num_new_false_clouds = dummy
        self.lattice.num_detected_height_low = dummy

    def test_calculate_kuipers(self):
        """Test kuipers calculation."""
        self.lattice.calculate_kuipers()
        self.assertTrue(np.abs(self.lattice.Kuipers[0] + 0.79) < 0.01)
        self.assertTrue(np.abs(self.lattice.Kuipers[1] - 0.709) < 0.01)
        self.assertTrue(np.abs(self.lattice.Kuipers[2] - 0.899) < 0.01)
        self.assertTrue(self.lattice.Kuipers.mask[3])
        self.assertTrue(np.abs(self.lattice.Kuipers[4] - 0.40) < 0.01)
        self.assertTrue(np.abs(self.lattice.Kuipers[5] - 0.15) < 0.01)


def suite():
    """Test suite for test remap measurements on fibonacci grid."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(test_kuipers_plot_on_map))
    return mysuite


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
