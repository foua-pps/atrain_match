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
from atrain_match.statistics.orrb_stat_class import OrrbStats
import numpy as np
# -----------------------------------------------------


def bias_corrected_rms(rms, bias, N):
    # ie formula should be bcRMS= sqrt(RMS^2-c*bias^2),
    # where c=N/(N-1). However our Ns are usually large
    if N < 2:
        print("Warning too few elements to calculate bc-RMS")
        return -9
    cnn1 = N / (N - 1)
    return np.sqrt(rms * rms - cnn1 * bias * bias)


class CloudTopStats(OrrbStats):

    def do_stats(self):
        OrrbStats.do_stats(self)
        # for key in self.ac_data.keys():
        #    print key

        cal_all_samples = self.ac_data["cal_all_samples"]
        cal_low_samples = self.ac_data["cal_low_samples"]
        cal_medium_samples = self.ac_data["cal_medium_samples"]
        cal_high_samples = self.ac_data["cal_high_samples"]
        mean_error_cal_all_sum = self.ac_data["mean_error_cal_all_sum"]
        mean_error_cal_low_sum = self.ac_data["mean_error_cal_low_sum"]
        mean_error_cal_medium_sum = self.ac_data["mean_error_cal_medium_sum"]
        mean_error_cal_high_sum = self.ac_data["mean_error_cal_high_sum"]
        rms_error_cal_all_sum = self.ac_data["rms_error_cal_all_sum"]
        rms_error_cal_low_sum = self.ac_data["rms_error_cal_low_sum"]
        rms_error_cal_medium_sum = self.ac_data["rms_error_cal_medium_sum"]
        rms_error_cal_high_sum = self.ac_data["rms_error_cal_high_sum"]
        mae_error_cal_all_sum = self.ac_data["mae_error_cal_all_sum"]
        mae_error_cal_low_sum = self.ac_data["mae_error_cal_low_sum"]
        mae_error_cal_medium_sum = self.ac_data["mae_error_cal_medium_sum"]
        mae_error_cal_high_sum = self.ac_data["mae_error_cal_high_sum"]

        n_missed_ctth_all = self.ac_data["n_missed_ctth_all"]
        n_missed_cma_all = self.ac_data["n_missed_cma_all"]
        n_missed_ctth_low = self.ac_data["n_missed_ctth_low"]
        n_missed_cma_low = self.ac_data["n_missed_cma_low"]
        n_missed_ctth_medium = self.ac_data["n_missed_ctth_medium"]
        n_missed_cma_medium = self.ac_data["n_missed_cma_medium"]
        n_missed_ctth_high = self.ac_data["n_missed_ctth_high"]
        n_missed_cma_high = self.ac_data["n_missed_cma_high"]

        self.retrieval_rate_all = {}
        self.retrieval_rate_low = {}
        self.retrieval_rate_medium = {}
        self.retrieval_rate_high = {}
        self.estimate_pod_cloudy_all = {}
        self.estimate_pod_cloudy_low = {}
        self.estimate_pod_cloudy_medium = {}
        self.estimate_pod_cloudy_high = {}
        self.cal_all_samples = cal_all_samples
        self.cal_low_samples = cal_low_samples
        self.cal_medium_samples = cal_medium_samples
        self.cal_high_samples = cal_high_samples
        self.bias_cal_all = {}
        self.bias_cal_low = {}
        self.bias_cal_medium = {}
        self.bias_cal_high = {}
        self.rms_cal_all = {}
        self.rms_cal_low = {}
        self.rms_cal_medium = {}
        self.rms_cal_high = {}
        self.bcrms_cal_all = {}
        self.bcrms_cal_low = {}
        self.bcrms_cal_medium = {}
        self.bcrms_cal_high = {}
        self.mae_cal_all = {}
        self.mae_cal_low = {}
        self.mae_cal_medium = {}
        self.mae_cal_high = {}
        self.pe250_cal_all = {}
        self.pe250_cal_low = {}
        self.pe250_cal_medium = {}
        self.pe250_cal_high = {}
        self.pe2500_cal_all = {}
        self.pe2500_cal_low = {}
        self.pe2500_cal_medium = {}
        self.pe2500_cal_high = {}
        self.pe500_cal_all = {}
        self.pe500_cal_low = {}
        self.pe500_cal_medium = {}
        self.pe500_cal_high = {}
        self.pe1000_cal_all = {}
        self.pe1000_cal_low = {}
        self.pe1000_cal_medium = {}
        self.pe1000_cal_high = {}

        for tc in cal_all_samples.keys():
            for clev in ["all", "low", "medium", "high"]:
            # numpy.divide handles potential division by zero
                N_CTTH = self.ac_data["cal_{clev}_samples".format(clev=clev)][tc]
                stats = getattr(self, "retrieval_rate_{clev}".format(clev=clev))
                stats[tc] = (
                    N_CTTH /
                    (N_CTTH + 
                     self.ac_data["n_missed_ctth_{clev}".format(clev=clev)][tc]))
                stats = getattr(self, "estimate_pod_cloudy_{clev}".format(clev=clev))
                stats[tc] = (
                    N_CTTH /
                    (N_CTTH + 
                     self.ac_data["n_missed_cma_{clev}".format(clev=clev)][tc]))
                stats = getattr(self, "bias_cal_{clev}".format(clev=clev))
                stats[tc] = np.divide(
                    self.ac_data["mean_error_cal_{clev}_sum".format(clev=clev)][tc], 
                    N_CTTH)
                stats = getattr(self, "mae_cal_{clev}".format(clev=clev))
                stats[tc] = np.divide(
                    self.ac_data["mae_error_cal_{clev}_sum".format(clev=clev)][tc], 
                    N_CTTH)
                stats = getattr(self, "rms_cal_{clev}".format(clev=clev))
                stats[tc] = math.sqrt(
                    np.divide(self.ac_data["rms_error_cal_{clev}_sum".format(clev=clev)][tc], 
                              N_CTTH))
                stats = getattr(self, "bcrms_cal_{clev}".format(clev=clev))
                stats[tc] = bias_corrected_rms(
                    getattr(self, "rms_cal_{clev}".format(clev=clev))[tc],
                    getattr(self, "bias_cal_{clev}".format(clev=clev))[tc],
                    N_CTTH)
                stats = getattr(self, "pe250_cal_{clev}".format(clev=clev))
                stats[tc] = np.divide(
                     100.0*self.ac_data["n_over_250_cal_{clev}".format(clev=clev)][tc], 
                    N_CTTH)
                stats = getattr(self, "pe500_cal_{clev}".format(clev=clev))
                stats[tc] = np.divide(
                     100.0*self.ac_data["n_over_500_cal_{clev}".format(clev=clev)][tc], 
                    N_CTTH)
                stats = getattr(self, "pe2500_cal_{clev}".format(clev=clev))
                stats[tc] = np.divide(
                     100.0*self.ac_data["n_over_2500_cal_{clev}".format(clev=clev)][tc], 
                    N_CTTH)    
                stats = getattr(self, "pe1000_cal_{clev}".format(clev=clev))
                stats[tc] = np.divide(
                    100.0*self.ac_data["n_over_1000_cal_{clev}".format(clev=clev)][tc], 
                    N_CTTH)
                

    def printout(self):
        lines = []
        for tc in sorted(self.cal_low_samples.keys()):
            lines.append("========== Cloud top height ===========")
            lines.append("====== %s ======" % (tc))
            lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
            lines.append("")
            lines.append("Imager retrieval rate: %3.2f" % self.retrieval_rate_all[tc])
            lines.append("Total number of %s matched cloudtops: %d" %
                         (self.truth_sat.upper(), self.cal_all_samples[tc]))
            lines.append("Number of %s matched low cloudtops: %d" % (self.truth_sat.upper(), self.cal_low_samples[tc]))
            lines.append("Number of %s matched medium cloudtops: %d" %
                         (self.truth_sat.upper(), self.cal_medium_samples[tc]))
            lines.append("Number of %s matched high cloudtops: %d" %
                         (self.truth_sat.upper(), self.cal_high_samples[tc]))
            lines.append("Mean absolute error total cases: %.0f" % self.mae_cal_all[tc])
            lines.append("Mean absolute error low-level cases: %.0f" % self.mae_cal_low[tc])
            lines.append("Mean absolute error medium-level cases: %.0f" % self.mae_cal_medium[tc])
            lines.append("Mean absolute error high-level cases: %.0f" % self.mae_cal_high[tc])
            lines.append("Part of error above 500m total cases: %.0f" % self.pe500_cal_all[tc])
            lines.append("Part of error above 500m low-level cases: %.0f" % self.pe500_cal_low[tc])
            lines.append("Part of error above 500m medium-level cases: %.0f" % self.pe500_cal_medium[tc])
            lines.append("Part of error above 500m high-level cases: %.0f" % self.pe500_cal_high[tc])
            lines.append("Mean error total cases: %.0f" % self.bias_cal_all[tc])
            lines.append("Mean error low-level cases: %.0f" % self.bias_cal_low[tc])
            lines.append("Mean error medium-level cases: %.0f" % self.bias_cal_medium[tc])
            lines.append("Mean error high-level cases: %.0f" % self.bias_cal_high[tc])
            lines.append("RMS error total cases: %.0f" % self.rms_cal_all[tc])
            lines.append("RMS error low-level cases: %.0f" % self.rms_cal_low[tc])
            lines.append("RMS error medium-level cases: %.0f" % self.rms_cal_medium[tc])
            lines.append("RMS error high-level cases: %.0f" % self.rms_cal_high[tc])
            lines.append("bc-RMS error total cases: %.0f" % self.bcrms_cal_all[tc])
            lines.append("bc-RMS error low-level cases: %.0f" % self.bcrms_cal_low[tc])
            lines.append("bc-RMS error medium-level cases: %.0f" % self.bcrms_cal_medium[tc])
            lines.append("bc-RMS error high-level cases: %.0f" % self.bcrms_cal_high[tc])
            lines.append("Imager estimated POD-cloudy total-level cases: %3.2f" % self.estimate_pod_cloudy_all[tc])
            lines.append("Imager estimated POD-cloudy low-level cases: %3.2f" % self.estimate_pod_cloudy_low[tc])
            lines.append("Imager estimated POD-cloudy medium-level cases: %3.2f" % self.estimate_pod_cloudy_medium[tc])
            lines.append("Imager estimated POD-cloudy high-level cases: %3.2f" % self.estimate_pod_cloudy_high[tc])
            lines.append("Part of error above 2500m total cases: %.0f" % self.pe2500_cal_all[tc])
            lines.append("Part of error above 2500m low-level cases: %.0f" % self.pe2500_cal_low[tc])
            lines.append("Part of error above 2500m medium-level cases: %.0f" % self.pe2500_cal_medium[tc])
            lines.append("Part of error above 2500m high-level cases: %.0f" % self.pe2500_cal_high[tc])
            lines.append("Part of error above 250m total cases: %.0f" % self.pe250_cal_all[tc])
            lines.append("Part of error above 250m low-level cases: %.0f" % self.pe250_cal_low[tc])
            lines.append("Part of error above 250m medium-level cases: %.0f" % self.pe250_cal_medium[tc])
            lines.append("Part of error above 250m high-level cases: %.0f" % self.pe250_cal_high[tc])
            lines.append("Part of error above 1000m total cases: %.0f" % self.pe1000_cal_all[tc])
            lines.append("Part of error above 1000m low-level cases: %.0f" % self.pe1000_cal_low[tc])
            lines.append("Part of error above 1000m medium-level cases: %.0f" % self.pe1000_cal_medium[tc])
            lines.append("Part of error above 1000m high-level cases: %.0f" % self.pe1000_cal_high[tc])
            lines.append("")
        for tc in self.cal_all_samples.keys():
            if tc in self.cal_low_samples.keys():
                # already go these results printed
                continue
            else:
                lines.append("========== Cloud top height ===========")
                lines.append("====== %s ======" % (tc))
                lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
                lines.append("")
                lines.append("Imager retrieval rate: %3.2f" % self.retrieval_rate_all[tc])
                lines.append("Total number of %s matched cloudtops: %d" %
                             (self.truth_sat.upper(), self.cal_all_samples[tc]))
                lines.append("Mean absolute error total cases: %.0f" % self.mae_cal_all[tc])
                lines.append("Mean error total cases: %.0f" % self.bias_cal_all[tc])
                lines.append("RMS error total cases: %.0f" % self.rms_cal_all[tc])
                lines.append("bc-RMS error total cases: %.0f" % self.bcrms_cal_all[tc])
                lines.append("Imager estimated POD-cloudy total-level cases: %3.2f" % self.estimate_pod_cloudy_all[tc])
                lines.append("")
        return lines


if __name__ == "__main__":
    stats = CloudTopStats()
