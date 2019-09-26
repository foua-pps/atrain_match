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

from statistics.orrb_stat_class import OrrbStats

# -----------------------------------------------------


class CloudPhaseStats(OrrbStats):
    def do_stats(self):
        OrrbStats.do_stats(self)

        n_ice_ice_cal = self.ac_data["n_ice_ice_cal"]
        n_ice_water_cal = self.ac_data["n_ice_water_cal"]
        n_water_ice_cal = self.ac_data["n_water_ice_cal"]
        n_water_water_cal = self.ac_data["n_water_water_cal"]
        Num = self.ac_data["Num"]

        pod_water_cal = np.divide(100.0 * n_water_water_cal,
                                  (n_water_water_cal + n_water_ice_cal))
        pod_ice_cal = np.divide(100.0 * n_ice_ice_cal,
                                n_ice_ice_cal + n_ice_water_cal)
        far_water_cal = np.divide(100.0 * n_ice_water_cal,
                                  n_water_water_cal + n_ice_water_cal)
        far_ice_cal = np.divide(100.0 * n_water_ice_cal,
                                n_ice_ice_cal + n_water_ice_cal)

        kuipers = np.divide(
            1.0 * (n_ice_ice_cal * n_water_water_cal -
                   n_water_ice_cal * n_ice_water_cal),
            ((n_ice_ice_cal + n_ice_water_cal) *
             (n_water_ice_cal + n_water_water_cal)))

        hitrate = np.divide(
            1.0 * (n_ice_ice_cal + n_water_water_cal),
            (n_ice_ice_cal + n_ice_water_cal +
             n_water_ice_cal + n_water_water_cal))

        # Store values of interest as attributes
        self.Num = Num
        self.pod_water_cal = pod_water_cal
        self.pod_ice_cal = pod_ice_cal
        self.far_water_cal = far_water_cal
        self.far_ice_cal = far_ice_cal
        self.cph_kuipers = kuipers
        self.cph_hitrate = hitrate
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
        lines.append("POD water: %.2f" % self.pod_water_cal)
        lines.append("FAR water: %.2f" % self.far_water_cal)
        lines.append("POD ice: %.2f" % self.pod_ice_cal)
        lines.append("FAR ice: %.2f" % self.far_ice_cal)
        lines.append("CPH Hitrate: %.2f" % self.cph_hitrate)
        lines.append("CPH Kuipers: %.2f" % self.cph_kuipers)
        lines.append("CPH N Water: %.2f" % self.n_water)
        lines.append("CPH N Ice: %.2f" % self.n_ice)
        lines.append("CPH N Water ok: %.2f" % self.n_water_water)
        lines.append("CPH N Ice ok: %.2f" % self.n_ice_ice)
        lines.append("")
        return lines


if __name__ == "__main__":
    stats = CloudPhaseStats()
