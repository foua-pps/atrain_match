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
# Program orrb_CTY_stat

# This program calculates basic statistics for the cloud type (CTY) product for
# each month

import numpy as np
from statistics.orrb_stat_class import OrrbStats

# -----------------------------------------------------


class CloudTypeStats(OrrbStats):

    def do_stats(self):
        OrrbStats.do_stats(self)

        Num_cma = self.ac_data["Num"]
        CFC_TRUTH = np.divide(100.0 * (self.ac_data["n_cloudy_cloudy_cal"] +
                                       self.ac_data["n_cloudy_clear_cal"]),
                              Num_cma)
        CFC_PPS = np.divide(100.0 * (self.ac_data["n_cloudy_cloudy_cal"] +
                                      self.ac_data["n_clear_cloudy_cal"]),
                            Num_cma)

        n_low_low = self.ac_data["n_low_low"]
        n_low_medium = self.ac_data["n_low_medium"]
        n_low_high = self.ac_data["n_low_high"]
        n_medium_low = self.ac_data["n_medium_low"]
        n_medium_medium = self.ac_data["n_medium_medium"]
        n_medium_high = self.ac_data["n_medium_high"]
        n_high_low = self.ac_data["n_high_low"]
        n_high_medium = self.ac_data["n_high_medium"]
        n_high_high = self.ac_data["n_high_high"]
        n_cirrus_low = self.ac_data["n_cirrus_low"]
        n_cirrus_medium_tp = self.ac_data["n_cirrus_medium_tp"]
        n_cirrus_high_tp = self.ac_data["n_cirrus_high_tp"]
        n_cirrus_medium_op = self.ac_data["n_cirrus_medium_op"]
        n_cirrus_high_op = self.ac_data["n_cirrus_high_op"]
        # This is really CMA buisness!
        pps_undetected_low = self.ac_data["n_clear_low"]
        pps_undetected_medium = self.ac_data["n_clear_medium"]
        pps_undetected_high = self.ac_data["n_clear_high"]
        pps_false_low = self.ac_data["n_low_clear"]  # This is really CMA buisness!
        pps_false_medium = self.ac_data["n_medium_clear"]
        pps_false_high = self.ac_data["n_high_clear"]
        pps_false_cirrus = self.ac_data["n_cirrus_clear"]
        pps_false = (pps_false_low + pps_false_medium +
                     pps_false_high + pps_false_cirrus )
        pps_undetected = (pps_undetected_low + pps_undetected_medium +
                          pps_undetected_high)
        common_cloud_free = self.ac_data["n_clear_clear_cal"]

        pps_cirrus = (n_cirrus_low +
                      n_cirrus_medium_op + n_cirrus_medium_tp +
                      n_cirrus_high_op + n_cirrus_high_tp)

        Num_ct = (n_low_low + n_low_medium + n_low_high +
                  n_medium_low + n_medium_medium + n_medium_high +
                  n_high_low + n_high_medium + n_high_high +
                  n_cirrus_low +
                  n_cirrus_medium_op + n_cirrus_high_op +
                  n_cirrus_medium_tp + n_cirrus_high_tp )

        N_low_cal = (n_low_low + n_medium_low + n_high_low
                     + n_cirrus_low)
        N_medium_cal = (n_low_medium + n_medium_medium + n_high_medium +
                        + n_cirrus_medium_tp + + n_cirrus_medium_op)
        N_high_cal = (n_low_high + n_medium_high + n_high_high +
                       + n_cirrus_high_op + n_cirrus_high_tp)

        N_low_pps = (n_low_low + n_low_medium + n_low_high )
        N_medium_pps = (n_medium_low + n_medium_medium +
                        n_medium_high)
        N_high_pps = (n_high_low + n_high_medium + n_high_high)

        N_cirrus_pps = pps_cirrus

        # numpy.divide handles potential division by zero
        pod_low = 100.0 * np.divide(float(n_low_low), N_low_cal)
        far_low = 100.0 * np.divide(float(n_low_medium +
                                          n_low_high), N_low_pps)
        pod_medium = 100.0 * np.divide(
            float(n_medium_medium + n_cirrus_medium_tp), N_medium_cal)
        far_medium = 100.0 * np.divide(
            float(n_medium_low + n_medium_high), N_medium_pps)
        pod_high = 100.0 * np.divide(float(n_high_high + n_cirrus_high_tp), N_high_cal)
        far_high = 100.0 * np.divide(
            float(n_high_low + n_high_medium), N_high_pps)
        far_cirrus = 100.0 * np.divide(float(n_cirrus_low + n_cirrus_high_op + n_cirrus_medium_op ), N_cirrus_pps)

        low_fraction_pps_rel = np.divide(float(N_low_pps), Num_ct) * 100.0
        low_fraction_cal_rel = np.divide(float(N_low_cal), Num_ct) * 100.0
        medium_fraction_pps_rel = np.divide(float(N_medium_pps), Num_ct) * 100.0
        medium_fraction_cal_rel = np.divide(float(N_medium_cal), Num_ct) * 100.0
        high_fraction_pps_rel = np.divide(float(N_high_pps), Num_ct) * 100.0
        high_fraction_cal_rel = np.divide(float(N_high_cal), Num_ct) * 100.0
        cirrus_fraction_pps_rel = np.divide(float(pps_cirrus), Num_ct) * 100.0

        low_fraction_pps_abs = low_fraction_pps_rel * 0.01 * CFC_PPS
        low_fraction_cal_abs = low_fraction_cal_rel * 0.01 * CFC_TRUTH
        medium_fraction_pps_abs = medium_fraction_pps_rel * 0.01 * CFC_PPS
        medium_fraction_cal_abs = medium_fraction_cal_rel * 0.01 * CFC_TRUTH
        high_fraction_pps_abs = high_fraction_pps_rel * 0.01 * CFC_PPS
        high_fraction_cal_abs = high_fraction_cal_rel * 0.01 * CFC_TRUTH
        cirrus_fraction_pps_abs = cirrus_fraction_pps_rel * 0.01 * CFC_PPS

        # Removed bc-RMS Not sure how ot interpret it

        # POD, FAR, HR and KSS calculations =============================================================

        hitrate = np.divide(
            100.0 * (
                n_low_low + n_medium_medium + n_high_high +
                n_cirrus_medium_tp + n_cirrus_high_tp),
            Num_ct)

        self.Num_cma = Num_cma
        self.common_cloud_free = common_cloud_free
        self.Num_ct = Num_ct
        self.pod_low = pod_low
        self.pod_medium = pod_medium
        self.pod_high = pod_high
        self.far_low = far_low
        self.far_medium = far_medium
        self.far_high = far_high
        self.low_fraction_pps_rel = low_fraction_pps_rel
        self.low_fraction_cal_rel = low_fraction_cal_rel
        self.medium_fraction_pps_rel = medium_fraction_pps_rel
        self.medium_fraction_cal_rel = medium_fraction_cal_rel
        self.high_fraction_pps_rel = high_fraction_pps_rel
        self.high_fraction_cal_rel = high_fraction_cal_rel
        self.low_fraction_pps_abs = low_fraction_pps_abs
        self.low_fraction_cal_abs = low_fraction_cal_abs
        self.medium_fraction_pps_abs = medium_fraction_pps_abs
        self.medium_fraction_cal_abs = medium_fraction_cal_abs
        self.high_fraction_pps_abs = high_fraction_pps_abs
        self.high_fraction_cal_abs = high_fraction_cal_abs
        self.cirrus_fraction_pps_rel = cirrus_fraction_pps_rel
        self.cirrus_fraction_pps_abs = cirrus_fraction_pps_abs
        self.pod_low = pod_low
        self.pod_medium = pod_medium
        self.pod_high = pod_high
        self.far_low = far_low
        self.far_medium = far_medium
        self.far_high = far_high
        self.far_cirrus = far_cirrus
        self.hitrate = hitrate
        self.pps_undetected = pps_undetected
        self.pps_undetected_low = pps_undetected_low
        self.pps_undetected_medium = pps_undetected_medium
        self.pps_undetected_high = pps_undetected_high
        self.pps_false = pps_false
        self.pps_false_low = pps_false_low
        self.pps_false_medium = pps_false_medium
        self.pps_false_high = pps_false_high
        self.pps_false_cirrus = pps_false_cirrus
    def printout(self):
        lines = []
        if self.Num_cma == 0 or self.Num_ct == 0:
            return lines
        lines.append( "Total pixels: %d" % self.Num_cma)
        lines.append( "Common cloud-free: %d" % self.common_cloud_free)
        lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
        lines.append("Total number of matched cloud types: %d" % self.Num_ct)
        lines.append("")
        lines.append("Note: There is a separate script to get CT-statitcs for validaion report.")
        lines.append("Probability of detecting LOW, MEDIUM and HIGH: %.2f %.2f %.2f" % \
                     (self.pod_low, self.pod_medium, self.pod_high))
        lines.append("False alarm rate LOW, MEDIUM and HIGH: %.2f %.2f %.2f %.2f" % \
                     (self.far_low, self.far_medium, self.far_high, self.far_cirrus))

        lines.append("Rel. Fraction LOW for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(), self.low_fraction_pps_rel, self.low_fraction_cal_rel))
        lines.append("Rel. Fraction MEDIUM for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(), self.medium_fraction_pps_rel, self.medium_fraction_cal_rel))
        lines.append("Rel. Fraction HIGH for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(), self.high_fraction_pps_rel, self.high_fraction_cal_rel))
        lines.append("Rel. Fraction CIRRUS for PPS: %.2f" % self.cirrus_fraction_pps_rel)

        lines.append("Abs. Fraction LOW for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(), self.low_fraction_pps_abs, self.low_fraction_cal_abs))
        lines.append("Abs. Fraction MEDIUM for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(), self.medium_fraction_pps_abs, self.medium_fraction_cal_abs))
        lines.append("Abs. Fraction HIGH for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(), self.high_fraction_pps_abs, self.high_fraction_cal_abs))
        lines.append("Abs. Fraction CIRRUS for PPS: %.2f" % self.cirrus_fraction_pps_abs)
        lines.append("POD (Low, Medium, High): %.2f %.2f %.2f" % \
                     (self.pod_low, self.pod_medium, self.pod_high))
        lines.append("FAR (Low, Medium, High, Cirrus): %.2f %.2f %.2f %.2f" % \
                     (self.far_low, self.far_medium, self.far_high, self.far_cirrus))
        lines.append("HR: %.3f " % (self.hitrate))
        lines.append("This is really CMA buisness here!")
        lines.append("Missclassified cloudy: %d (low:%d, medium:%d, high:%d, cirrus:%d)" % \
                     (self.pps_false, self.pps_false_low,
                      self.pps_false_medium, self.pps_false_high,
                      self.pps_false_cirrus))
        lines.append("Missclassified clear: %d (low:%d, medium:%d, high:%d)" % \
                     (self.pps_undetected, self.pps_undetected_low,
                      self.pps_undetected_medium, self.pps_undetected_high))
        lines.append("")

        return lines


if __name__ == "__main__":
    stats = CloudTypeStats()

