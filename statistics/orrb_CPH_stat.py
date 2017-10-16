# Program orrb_CPH_stat

# This program calculates basic statistics for the cloud amount (CPH) product for
# each month

import string
import math
import numpy as np

from orrb_stat_class import OrrbStats

# -----------------------------------------------------
class CloudPhaseStats(OrrbStats):
    def do_stats(self):
        OrrbStats.do_stats(self)
        
        n_ice_ice_cal = self.ac_data["n_ice_ice_cal"]
        n_ice_water_cal = self.ac_data["n_ice_water_cal"]
        n_water_ice_cal = self.ac_data["n_water_ice_cal"]
        n_water_water_cal = self.ac_data["n_water_water_cal"]
        Num = self.ac_data["Num"]

        pod_water_cal = np.divide(100.0*n_water_water_cal, (n_water_water_cal+n_water_ice_cal))
        pod_ice_cal = np.divide(100.0*n_ice_ice_cal, n_ice_ice_cal+n_ice_water_cal)
        far_water_cal = np.divide(100.0*n_ice_water_cal, n_water_water_cal+n_ice_water_cal)
        far_ice_cal = np.divide(100.0*n_water_ice_cal, n_ice_ice_cal+n_water_ice_cal)
    
        kuipers = np.divide(1.0*(n_ice_ice_cal*n_water_water_cal-n_water_ice_cal*n_ice_water_cal),
                         ((n_ice_ice_cal+n_ice_water_cal)*(n_water_ice_cal+n_water_water_cal)))
    
        hitrate = np.divide(1.0*(n_ice_ice_cal+n_water_water_cal),
                         (n_ice_ice_cal+n_ice_water_cal+n_water_ice_cal+n_water_water_cal))
    
        # Store values of interest as attributes
        self.Num = Num
        self.pod_water_cal = pod_water_cal
        self.pod_ice_cal = pod_ice_cal
        self.far_water_cal = far_water_cal
        self.far_ice_cal = far_ice_cal
        self.cph_kuipers = kuipers
        self.cph_hitrate = hitrate
        

    def printout(self):
        lines = []
        if self.Num == 0:
            return lines
        lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
        lines.append("Total number of %s matched FOVs: %d" % (self.truth_sat.upper(), self.Num))
        lines.append("POD water: %.2f" % self.pod_water_cal)
        lines.append("FAR water: %.2f" % self.far_water_cal)
        lines.append("POD ice: %.2f" % self.pod_ice_cal)
        lines.append("FAR ice: %.2f" % self.far_ice_cal)
        lines.append("CPH Hitrate: %.2f" % self.cph_hitrate)
        lines.append("")
        return lines


if __name__ == "__main__":
    stats = CloudFractionStats()

