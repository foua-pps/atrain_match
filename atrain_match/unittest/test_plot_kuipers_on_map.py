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
# test_match.py
"""for the match module match.py """
#
import numpy as np
import unittest
from reshaped_files_plotting.plot_kuipers_on_area_util import ppsStatsOnFibLatticeObject

class test_kuipers_plot_on_map(unittest.TestCase): 

    def setUp(self):
        self.area = ppsStatsOnFibLatticeObject()
        self.area.N_detected_clouds =   np.array([10,  1000, 1000, 1000, 1000, 1000])
        self.area.N_undetected_clouds = np.array([1000, 100,  100,  100,  100,  100]) 
        self.area.N_detected_clear =    np.array([20,    80,  100,    0,  200,  100])
        self.area.N_false_clouds =      np.array([80,    20,    1,    1,  200,  300])

        dummy = 0*np.array([80,    20,    1,    1,  200,  300])
        self.area.Sum_ctth_mae_low = dummy
        self.area.Sum_ctth_mae_high = dummy
        self.area.Sum_ctth_mae = dummy
        self.area.Sum_ctth_mae_diff = dummy

        self.area.Sum_ctth_bias_low = dummy
        self.area.Sum_lapse_bias_low = dummy
        self.area.Sum_ctth_bias_temperature_low = dummy
        self.area.Sum_ctth_bias_temperature_low_t11 = dummy
        self.area.Min_lapse_rate = dummy
        self.area.N_new_detected_clouds = dummy
        self.area.N_new_false_clouds = dummy
        self.area.N_detected_height_low = dummy


    def test_calculate_kuipers(self):
        self.area.calculate_kuipers()
        #print 
        #print self.area.Kuipers                                                  
        self.assertTrue(np.abs(self.area.Kuipers[0]+0.79) <0.01)
        self.assertTrue(np.abs(self.area.Kuipers[1]-0.709)<0.01)
        self.assertTrue(np.abs(self.area.Kuipers[2]-0.899)<0.01)
        self.assertTrue(self.area.Kuipers.mask[3])
        self.assertTrue(np.abs(self.area.Kuipers[4]-0.40) <0.01)
        self.assertTrue(np.abs(self.area.Kuipers[5]-0.15) <0.01)


def suite():
    """The suite for test_utils.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(test_kuipers_plot_on_map))
    return mysuite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
