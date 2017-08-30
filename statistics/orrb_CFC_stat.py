# Program orrb_CFC_stat

# This program calculates basic statistics for the cloud amount (CFC) product for
# each month

import string
import math
import numpy as np

from orrb_stat_class import OrrbStats

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
    
        #MODIS    
        if self.ac_data["got_cloudsat_modis_flag"]:
            bias_modis = np.divide(1.*(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS), Num-1.)
            bias_modis_perc = np.divide(100.0*(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS), Num-1.)
            square_sum_modis =  (n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_modis**2 + \
                                n_cloudy_clear_cal*(-1.0-bias_modis)**2 + \
                                n_clear_cloudy_cal*(1.0-bias_modis)**2
            rms_modis = 100.0*math.sqrt(np.divide(square_sum_modis, Num-1.))
            pod_cloudy_cal_MODIS = np.divide(100.0*n_cloudy_cloudy_cal_MODIS, n_cloudy_cloudy_cal_MODIS+n_cloudy_clear_cal_MODIS)
            pod_clear_cal_MODIS = np.divide(100.0*n_clear_clear_cal_MODIS, n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS)
            far_cloudy_cal_MODIS = np.divide(100.0*n_clear_cloudy_cal_MODIS, n_cloudy_cloudy_cal_MODIS+n_clear_cloudy_cal_MODIS)
            far_clear_cal_MODIS = np.divide(100.0*n_cloudy_clear_cal_MODIS, n_clear_clear_cal_MODIS+n_cloudy_clear_cal_MODIS)
            kuipers_MODIS = np.divide(1.0*(n_clear_clear_cal_MODIS*n_cloudy_cloudy_cal_MODIS-n_cloudy_clear_cal_MODIS*n_clear_cloudy_cal_MODIS),
                                      ((n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS)*(n_cloudy_clear_cal_MODIS+n_cloudy_cloudy_cal_MODIS)))        
            hitrate_MODIS = np.divide(1.0*(n_clear_clear_cal_MODIS+n_cloudy_cloudy_cal_MODIS),
                                      (n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS+n_cloudy_clear_cal_MODIS+ n_cloudy_cloudy_cal_MODIS))
        
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
        

    def printout(self):
        if self.ac_data["got_cloudsat_modis_flag"]:
            lines = []
            lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
            lines.append("Total number of %s matched FOVs: %d" % (self.truth_sat.upper(), self.Num))
            lines.append("Mean CFC %s: %.2f " %( self.truth_sat.upper(), self.mean_CFC_cal))
            lines.append("Mean error: %.2f" % self.bias_cal_perc)
            lines.append("RMS error: %.2f" % self.rms_cal)
            lines.append("Mean error MODIS: %.2f" % self.bias_modis_perc)
            lines.append("RMS error MODIS: %.2f" % self.rms_modis)
            lines.append("POD cloudy: %.2f" % self.pod_cloudy_cal)
            lines.append("POD cloudy MODIS: %.2f" % self.pod_cloudy_cal_MODIS)
            lines.append("POD clear: %.2f" % self.pod_clear_cal)
            lines.append("POD clear MODIS: %.2f" % self.pod_clear_cal_MODIS)
            lines.append("FAR cloudy: %.2f" % self.far_cloudy_cal)
            lines.append("FAR cloudy MODIS: %.2f" % self.far_cloudy_cal_MODIS)
            lines.append("FAR clear: %.2f" % self.far_clear_cal)
            lines.append("FAR clear MODIS: %.2f" % self.far_clear_cal_MODIS)
            lines.append("Kuipers: %.2f" % self.kuipers)
            lines.append("Kuipers MODIS: %.2f" % self.kuipers_MODIS)
            lines.append("Hitrate: %.3f" % self.hitrate)
            lines.append("Hitrate MODIS: %.3f" % self.hitrate_MODIS)
            lines.append("")
        else:
            lines = []
            lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
            lines.append("")
            lines.append("Total number of %s matched FOVs: %d" % (self.truth_sat.upper(), self.Num))
            lines.append("Mean CFC %s: %.2f " % (self.truth_sat.upper(), self.mean_CFC_cal))
            lines.append("Mean error: %.2f" % self.bias_cal_perc)
            lines.append("RMS error: %.2f" % self.rms_cal)
            lines.append("POD cloudy: %.2f" % self.pod_cloudy_cal)
            lines.append("POD clear: %.2f" % self.pod_clear_cal)
            lines.append("FAR cloudy: %.2f" % self.far_cloudy_cal)
            lines.append("FAR clear: %.2f" % self.far_clear_cal)
            lines.append("Kuipers: %.3f" % self.kuipers)
            lines.append("Hitrate: %.3f" % self.hitrate)
            lines.append("")
        return lines


if __name__ == "__main__":
    stats = CloudFractionStats()

