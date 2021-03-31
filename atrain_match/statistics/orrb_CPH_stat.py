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
# Program orrb_CPH_stat

# This program calculates basic statistics for the cloud amount (CPH) product for
# each month

import numpy as np

from atrain_match.statistics.orrb_stat_class import OrrbStats
from atrain_match.statistics.scores import ScoreUtils
# -----------------------------------------------------


class CloudPhaseStats(OrrbStats):
    def do_stats(self):
        OrrbStats.do_stats(self)

        n_ice_ice_cal = self.ac_data["n_ice_ice_cal"]  # a 
        n_ice_water_cal = self.ac_data["n_ice_water_cal"]  #c 
        n_water_ice_cal =  self.ac_data["n_water_ice_cal"]  # b
        n_water_water_cal = self.ac_data["n_water_water_cal"]  # d
        Num = self.ac_data["Num"]

        scu = ScoreUtils(n_ice_ice_cal,
                         n_water_ice_cal,
                         n_ice_water_cal,
                         n_water_water_cal)
        
        # Store values of interest as attributes
        self.Num = Num
        self.pod_water_cal = 100. * scu.pod_0()
        self.bias_cal = scu.bias()
        self.pod_ice_cal = 100. * scu.pod_1()
        self.far_water_cal = 100. * scu.far_0()
        self.far_ice_cal = 100 * scu.far_1()
        self.cph_kuipers = scu.kuiper()
        self.cph_heidke = scu.heidke()
        self.cph_hitrate = scu.hitrate()
        self.n_water = n_water_water_cal + n_water_ice_cal
        self.n_ice = n_ice_ice_cal + n_ice_water_cal
        self.n_water_water = n_water_water_cal
        self.n_ice_ice = n_ice_ice_cal

    def printout(self):
        lines = []
        if self.Num == 0:
            return lines
        lines.append("Total number of matched scenes is: %s"
                     % self.ac_data["scenes"])
        lines.append("Total number of %s matched FOVs: %d"
                     % (self.truth_sat.upper(), self.Num))
        lines.append("Mean error: %.2f" % self.bias_cal)
        lines.append("POD water: %.2f" % self.pod_water_cal)
        lines.append("FAR water: %.2f" % self.far_water_cal)
        lines.append("POD ice: %.2f" % self.pod_ice_cal)
        lines.append("FAR ice: %.2f" % self.far_ice_cal)
        lines.append("CPH Hitrate: %.2f" % self.cph_hitrate)
        lines.append("CPH Kuipers: %.2f" % self.cph_kuipers)
        lines.append("CPH Heidke: %.2f" % self.cph_heidke)
        lines.append("CPH N Water: %.2f" % self.n_water)
        lines.append("CPH N Ice: %.2f" % self.n_ice)
        lines.append("CPH N Water ok: %.2f" % self.n_water_water)
        lines.append("CPH N Ice ok: %.2f" % self.n_ice_ice)
        lines.append("")
        return lines


if __name__ == "__main__":
    stats = CloudPhaseStats()
