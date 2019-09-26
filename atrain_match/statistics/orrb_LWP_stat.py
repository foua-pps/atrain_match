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
# Program orrb_CTH_stat

# This program calculates basic statistics for the cloud top (CTH) product for
# each month

import math
from statistics.orrb_stat_class import OrrbStats
import numpy as np
# -----------------------------------------------------


def bias_corrected_rms(rms, bias, N):
    # ie formula should be bcRMS= sqrt(RMS^2-c*bias^2),
    # where c=N/(N-1). However our Ns are usually large
    if N < 2:
        print "Warning too few elements to calculate bc-RMS"
        return -9
    cnn1 = N / (N - 1)
    return np.sqrt(rms * rms - cnn1 * bias * bias)


class CloudLwpStats(OrrbStats):

    def do_stats(self):
        OrrbStats.do_stats(self)
        # for key in self.ac_data.keys():
        #    print key

        amsr_all_samples = self.ac_data["amsr_all_samples"]
        mean_error_amsr_all_sum = self.ac_data["mean_error_amsr_all_sum"]
        rms_error_amsr_all_sum = self.ac_data["rms_error_amsr_all_sum"]

        self.amsr_all_samples = amsr_all_samples
        self.bias_amsr_all = np.divide(
                mean_error_amsr_all_sum, amsr_all_samples )
        self.rms_amsr_all = np.sqrt(
            np.divide(rms_error_amsr_all_sum, amsr_all_samples))

        self.bcrms_amsr_all = bias_corrected_rms(
                self.rms_amsr_all,
                self.bias_amsr_all,
                amsr_all_samples)

    def printout(self):
        lines = []
        lines.append("========== Cloud lwp ===========")
        lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
        lines.append("")
        lines.append("Total number of %s matched cloudtops: %d" %( self.truth_sat.upper(), self.amsr_all_samples))
        lines.append("Mean error total cases: %.1f" % self.bias_amsr_all)
        lines.append("RMS error total cases: %.1f" % self.rms_amsr_all)
        lines.append("bc-RMS error total cases: %.1f" % self.bcrms_amsr_all)

        return lines

if __name__ == "__main__":
    stats = CloudTopStats()

