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
# Program orrb_CFC_stat

# This program calculates basic statistics for the cloud amount (CFC) product for
# each month

import math
import numpy as np

from statistics.orrb_stat_class import OrrbStats

# -----------------------------------------------------
class CloudFractionStats(OrrbStats):
    def do_stats(self):
        OrrbStats.do_stats(self)

        n_clear_clear_cal = self.ac_data["n_clear_clear_cal"]
        n_clear_cloudy_cal = self.ac_data["n_clear_cloudy_cal"]
        n_cloudy_clear_cal = self.ac_data["n_cloudy_clear_cal"]
        n_cloudy_cloudy_cal = self.ac_data["n_cloudy_cloudy_cal"]
        n_clear_clear_cal_MODIS = self.ac_data["n_clear_clear_cal_MODIS"]
        n_clear_cloudy_cal_MODIS = self.ac_data["n_clear_cloudy_cal_MODIS"]
        n_cloudy_clear_cal_MODIS = self.ac_data["n_cloudy_clear_cal_MODIS"]
        n_cloudy_cloudy_cal_MODIS = self.ac_data["n_cloudy_cloudy_cal_MODIS"]

        if "step_cmaprob" in self.ac_data.keys ():
            step_cmaprob = int(self.ac_data["step_cmaprob"])
            n_clear_cmaprob = np.array(self.ac_data["n_clear_cmaprob"])
            n_cloudy_cmaprob = np.array(self.ac_data["n_cloudy_cmaprob"])
            min_prob = np.array([percent*1.0 for percent in range(0, 100, step_cmaprob)])
            max_prob = np.array([percent*1.0 for percent in range(step_cmaprob, 100 + step_cmaprob, step_cmaprob)])
            limit_v =  np.array([percent*1.0 for percent in range(0, 101, step_cmaprob)])
            max_prob[-1] = 100
            Num_cloudy_tot = np.sum(n_cloudy_cmaprob)
            Num_clear_tot = np.sum(n_clear_cmaprob)
            percent_cloudy_prob = np.array([100.0/Num_cloudy_tot*np.int(nc) for nc in n_cloudy_cmaprob])
            print percent_cloudy_prob, Num_cloudy_tot, n_cloudy_cmaprob
            percent_clear_prob =  np.array([100.0/Num_clear_tot*np.int(nc)  for nc in n_clear_cmaprob])
            print percent_clear_prob


            detected_clouds = np.array([np.sum(n_cloudy_cmaprob[min_prob>=limit]) for limit in limit_v])
            undetected_clouds = np.array([np.sum(n_cloudy_cmaprob[min_prob<limit]) for limit in limit_v])
            detected_clear = np.array([np.sum(n_clear_cmaprob[max_prob<=limit]) for limit in limit_v])
            false_clouds = np.array([np.sum(n_clear_cmaprob[max_prob>limit]) for limit in limit_v])

            print max_prob
            print limit_v

            print false_clouds
            print detected_clear
            print detected_clouds
            print undetected_clouds

            pod_cloudy_prob = 100.0/Num_cloudy_tot*detected_clouds
            pod_clear_prob = 100.0/Num_clear_tot*detected_clear
            far_cloudy_prob = 100.0*false_clouds/(false_clouds + detected_clouds)
            far_clear_prob = 100.0*undetected_clouds/(undetected_clouds + detected_clear)

            hitrate_prob = (1.0/(Num_cloudy_tot + Num_clear_tot) *
                            (detected_clouds + detected_clear))
            if Num_cloudy_tot*Num_clear_tot ==0:
                kuipers_prob = [-9  for limit in limit_v]
            else:
                kuipers_prob = (
                    1.0 *(detected_clouds * detected_clear - undetected_clouds * false_clouds)/
                    (Num_cloudy_tot*Num_clear_tot))



        Num = self.ac_data["Num"]

        mean_CFC_cal = np.divide(100.0*(n_cloudy_cloudy_cal+n_cloudy_clear_cal), Num)
        bias_cal = np.divide(1.*(n_clear_cloudy_cal - n_cloudy_clear_cal), Num)
        bias_cal_perc = np.divide(100.0*(n_clear_cloudy_cal - n_cloudy_clear_cal), Num)
        square_sum_cal =  (n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_cal**2 + \
                            n_cloudy_clear_cal*(-1.0-bias_cal)**2 + \
                            n_clear_cloudy_cal*(1.0-bias_cal)**2
        rms_cal = 100.0*math.sqrt(np.divide(square_sum_cal, Num-1.))
        pod_cloudy_cal = np.divide(100.0*n_cloudy_cloudy_cal, (n_cloudy_cloudy_cal+n_cloudy_clear_cal))
        pod_clear_cal = np.divide(100.0*n_clear_clear_cal, n_clear_clear_cal+n_clear_cloudy_cal)
        far_cloudy_cal = np.divide(100.0*n_clear_cloudy_cal, n_cloudy_cloudy_cal+n_clear_cloudy_cal)
        far_clear_cal = np.divide(100.0*n_cloudy_clear_cal, n_clear_clear_cal+n_cloudy_clear_cal)

        kuipers = np.divide(1.0*(n_clear_clear_cal*n_cloudy_cloudy_cal-n_cloudy_clear_cal*n_clear_cloudy_cal),
                         ((n_clear_clear_cal+n_clear_cloudy_cal)*(n_cloudy_clear_cal+n_cloudy_cloudy_cal)))

        hitrate = np.divide(1.0*(n_clear_clear_cal+n_cloudy_cloudy_cal),
                         (n_clear_clear_cal+n_clear_cloudy_cal+n_cloudy_clear_cal+n_cloudy_cloudy_cal))

        # MODIS
        if self.ac_data["got_cloudsat_modis_flag"]:
            bias_modis = np.divide(1.*(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS), Num-1.)
            bias_modis_perc = np.divide(100.0*(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS), Num-1.)
            square_sum_modis =  (n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_modis**2 + \
                                n_cloudy_clear_cal*(-1.0-bias_modis)**2 + \
                                n_clear_cloudy_cal*(1.0-bias_modis)**2
            rms_modis = 100.0*math.sqrt(np.divide(square_sum_modis, Num-1.))
            pod_cloudy_cal_MODIS = np.divide(
                100.0 * n_cloudy_cloudy_cal_MODIS,
                n_cloudy_cloudy_cal_MODIS + n_cloudy_clear_cal_MODIS)
            pod_clear_cal_MODIS = np.divide(
                100.0 * n_clear_clear_cal_MODIS,
                n_clear_clear_cal_MODIS + n_clear_cloudy_cal_MODIS)
            far_cloudy_cal_MODIS = np.divide(
                100.0 * n_clear_cloudy_cal_MODIS,
                n_cloudy_cloudy_cal_MODIS+n_clear_cloudy_cal_MODIS)
            far_clear_cal_MODIS = np.divide(
                100.0 * n_cloudy_clear_cal_MODIS,
                n_clear_clear_cal_MODIS+n_cloudy_clear_cal_MODIS)
            kuipers_MODIS = np.divide(
                1.0 * (n_clear_clear_cal_MODIS * n_cloudy_cloudy_cal_MODIS -
                     n_cloudy_clear_cal_MODIS * n_clear_cloudy_cal_MODIS),
                ((n_clear_clear_cal_MODIS + n_clear_cloudy_cal_MODIS) *
                 (n_cloudy_clear_cal_MODIS + n_cloudy_cloudy_cal_MODIS)))
            hitrate_MODIS = np.divide(
                1.0 * (n_clear_clear_cal_MODIS + n_cloudy_cloudy_cal_MODIS),
                (n_clear_clear_cal_MODIS + n_clear_cloudy_cal_MODIS +
                 n_cloudy_clear_cal_MODIS + n_cloudy_cloudy_cal_MODIS))

        # Store values of interest as attributes
        if self.ac_data["got_cloudsat_modis_flag"]:
            self.bias_modis_perc = bias_modis_perc
            self.rms_modis = rms_modis
            self.pod_cloudy_cal_MODIS = pod_cloudy_cal_MODIS
            self.pod_clear_cal_MODIS = pod_clear_cal_MODIS
            self.far_cloudy_cal_MODIS = far_cloudy_cal_MODIS
            self.far_clear_cal_MODIS = far_clear_cal_MODIS
            self.kuipers_MODIS = kuipers_MODIS
            self.hitrate_MODIS = hitrate_MODIS
        self.Num = Num
        self.mean_CFC_cal = mean_CFC_cal
        self.bias_cal_perc = bias_cal_perc
        self.rms_cal = rms_cal
        self.pod_cloudy_cal = pod_cloudy_cal
        self.pod_clear_cal = pod_clear_cal
        self.far_cloudy_cal = far_cloudy_cal
        self.far_clear_cal = far_clear_cal
        self.kuipers = kuipers
        self.hitrate = hitrate
        if "step_cmaprob" in self.ac_data.keys():
            self.pod_clear_prob = pod_clear_prob
            self.pod_cloudy_prob = pod_cloudy_prob
            self.far_clear_prob = far_clear_prob
            self.far_cloudy_prob = far_cloudy_prob
            self.min_prob = min_prob
            self.max_prob = max_prob
            self.hitrate_prob = hitrate_prob
            self.kuipers_prob = kuipers_prob
            self.percent_clear_prob = percent_clear_prob
            self.percent_cloudy_prob = percent_cloudy_prob
    def printout(self):

        lines = []
        if self.Num == 0:
            return lines
        lines.append("Total number of matched scenes is: {:d}".format( self.ac_data["scenes"]))
        lines.append("")
        lines.append("Total number of {:s} matched FOVs: {:.0f}".format( self.truth_sat.upper(), self.Num))
        lines.append("Mean CFC {:s}: {:6.2f} ".format( self.truth_sat.upper(), self.mean_CFC_cal))
        lines.append("Mean error: {:6.2f}".format( self.bias_cal_perc))
        lines.append("RMS error:  {:6.2f}".format( self.rms_cal))
        lines.append("POD cloudy: {:6.2f}".format( self.pod_cloudy_cal))
        lines.append("POD clear:  {:6.2f}".format( self.pod_clear_cal))
        lines.append("FAR cloudy: {:6.2f}".format( self.far_cloudy_cal))
        lines.append("FAR clear:  {:6.2f}".format( self.far_clear_cal))
        lines.append("Kuipers: {:5.3f}".format( self.kuipers))
        lines.append("Hitrate: {:5.3f}".format( self.hitrate))
        lines.append("")
        if "step_cmaprob" in self.ac_data.keys():
            lines.append("Results for CMAPROB")
            lines.append('Limit  POD-cloudy  POD-clear FAR-cloudy  FAR-clear Hitrate(%)  Kuipers   ')
            for ind, limit in enumerate(np.append(self.min_prob, [100])):
                lines.append("L{:3.0f}: {:9.2f} {:9.2f} {:9.2f} {:9.2f} {:9.3f} {:9.3f}".format(
                    limit,
                    self.pod_cloudy_prob[ind],
                    self.pod_clear_prob[ind],
                    self.far_cloudy_prob[ind],
                    self.far_clear_prob[ind],
                    self.hitrate_prob[ind],
                    self.kuipers_prob[ind]))
            lines.append('CMAPROB-interval   Truth Clouds(%)  Truth Clears(%)')
            for ind, limit in enumerate(self.min_prob):
                upper_limit = "<=100"
                if (ind + 1) < len(self.min_prob):
                    upper_limit = "<{:4.0f}".format(self.min_prob[ind+1])
                lines.append("{:3.0f}<= cmaprob {:s}: {:9.2f}    {:9.3f}".format(
                    limit,
                    upper_limit,
                    self.percent_cloudy_prob[ind],
                    self.percent_clear_prob[ind]))
            lines.append("")

        if self.ac_data["got_cloudsat_modis_flag"]:
            lines.append("Results for MODIS flag from CPR (CloudSat)")
            lines.append("Mean error MODIS: {:6.2f}".format( self.bias_modis_perc))
            lines.append("RMS error MODIS: {:6.2f}".format( self.rms_modis))
            lines.append("POD cloudy MODIS: {:6.2f}".format( self.pod_cloudy_cal_MODIS))
            lines.append("POD clear MODIS: {:6.2f}".format( self.pod_clear_cal_MODIS))
            lines.append("FAR cloudy MODIS: {:6.2f}".format( self.far_cloudy_cal_MODIS))
            lines.append("FAR clear MODIS: {:6.2f}".format( self.far_clear_cal_MODIS))
            lines.append("Kuipers MODIS: {:5.3f}".format( self.kuipers_MODIS))
            lines.append("Hitrate MODIS: {:5.3f}".format( self.hitrate_MODIS))
        return lines


if __name__ == "__main__":
    stats = CloudFractionStats()

