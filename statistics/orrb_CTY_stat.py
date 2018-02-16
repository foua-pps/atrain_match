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
        n_frac_low = self.ac_data["n_frac_low"] 
        n_frac_medium = self.ac_data["n_frac_medium"] 
        n_frac_high = self.ac_data["n_frac_high"] 
        n_cirrus_low = self.ac_data["n_cirrus_low"] 
        n_cirrus_medium = self.ac_data["n_cirrus_medium"] 
        n_cirrus_high = self.ac_data["n_cirrus_high"] 
        # This is really CMA buisness!
        pps_undetected_low = self.ac_data["n_clear_low"] 
        pps_undetected_medium = self.ac_data["n_clear_medium"] 
        pps_undetected_high = self.ac_data["n_clear_high"] 
        pps_false_low = self.ac_data["n_low_clear"] # This is really CMA buisness!
        pps_false_medium = self.ac_data["n_medium_clear"] 
        pps_false_high = self.ac_data["n_high_clear"] 
        pps_false_frac = self.ac_data["n_frac_clear"] 
        pps_false_cirrus = self.ac_data["n_cirrus_clear"] 
        pps_false = (pps_false_low + pps_false_medium + 
                     pps_false_high + pps_false_frac + pps_false_cirrus )
        pps_undetected = (pps_undetected_low + pps_undetected_medium + 
                          pps_undetected_high)
        number_excluded_cirrus_medium = n_cirrus_medium 
        #n_cirrus_medium = 0 # Lets not include these!
        common_cloud_free = self.ac_data["n_clear_clear_cal"]
   
        pps_frac = n_frac_low + n_frac_medium + n_frac_high
        pps_cirrus = n_cirrus_low + n_cirrus_medium + n_cirrus_high
      

        Num_ct = (n_low_low + n_low_medium + n_low_high + 
                  n_medium_low + n_medium_medium + n_medium_high + 
                  n_high_low + n_high_medium + n_high_high + 
                  n_frac_low + n_frac_medium +n_frac_high +
                  n_cirrus_low + n_cirrus_medium +n_cirrus_high)

        N_low_cal = (n_low_low + n_medium_low + n_high_low 
                     + n_frac_low + n_cirrus_low)
        N_medium_cal = (n_low_medium + n_medium_medium + n_high_medium + 
                        n_frac_medium + n_cirrus_medium)
        N_high_cal = (n_low_high + n_medium_high + n_high_high + 
                      n_frac_high + n_cirrus_high)
     
        N_low_pps = (n_low_low + n_low_medium + n_low_high +
                     n_frac_low + n_frac_medium + n_frac_high)
        N_medium_pps = (n_medium_low + n_medium_medium + 
                        n_medium_high + n_cirrus_medium)
        N_high_pps = (n_high_low + n_high_medium + n_high_high + 
                      n_cirrus_low + n_cirrus_high) #n_cirrus_medium belong to medium class!
    
     
        # numpy.divide handles potential division by zero
        pod_low = 100.0 * np.divide(float(n_low_low + n_frac_low), N_low_cal)
        far_low = 100.0 * np.divide(float(n_low_medium + n_frac_medium + 
                                          n_low_high + n_frac_high), N_low_pps)
        pod_medium = 100.0 * np.divide(
            float(n_medium_medium + n_cirrus_medium), N_medium_cal)
        far_medium = 100.0 * np.divide(
            float(n_medium_low + n_medium_high), N_medium_pps)
        pod_high = 100.0 * np.divide(float(n_high_high + n_cirrus_high), N_high_cal)
        far_high = 100.0 * np.divide(
            float(n_high_low + n_high_medium + n_cirrus_low), N_high_pps)
    
        low_fraction_pps_rel = np.divide(float(N_low_pps), Num_ct) * 100.0
        low_fraction_cal_rel = np.divide(float(N_low_cal), Num_ct) * 100.0
        medium_fraction_pps_rel = np.divide(float(N_medium_pps), Num_ct) * 100.0
        medium_fraction_cal_rel = np.divide(float(N_medium_cal), Num_ct) * 100.0
        high_fraction_pps_rel = np.divide(float(N_high_pps), Num_ct) * 100.0
        high_fraction_cal_rel = np.divide(float(N_high_cal), Num_ct) * 100.0
        frac_fraction_pps_rel = np.divide(pps_frac, Num_ct) * 100.0

        low_fraction_pps_abs = low_fraction_pps_rel * 0.01 *CFC_PPS
        low_fraction_cal_abs = low_fraction_cal_rel * 0.01 *CFC_TRUTH
        medium_fraction_pps_abs = medium_fraction_pps_rel * 0.01 *CFC_PPS
        medium_fraction_cal_abs = medium_fraction_cal_rel * 0.01 *CFC_TRUTH
        high_fraction_pps_abs = high_fraction_pps_rel * 0.01 *CFC_PPS
        high_fraction_cal_abs = high_fraction_cal_rel * 0.01 *CFC_TRUTH
        frac_fraction_pps_abs =  frac_fraction_pps_rel * 0.01 *CFC_PPS
        
        bias_low = low_fraction_pps_abs - low_fraction_cal_abs
        bias_low_noperc = bias_low / 100.0
        bias_medium = medium_fraction_pps_abs - medium_fraction_cal_abs
        bias_medium_noperc = bias_medium / 100.0
        bias_high = high_fraction_pps_abs - high_fraction_cal_abs
        bias_high_noperc = bias_high / 100.0

        #Removed bc-RMS Not sure how ot interpret it                     
        
        #POD,FAR,HR and KSS calculations =============================================================

        hitrate = np.divide(
            100.0*(n_low_low + n_frac_low + n_medium_medium + n_high_high),
            Num_ct)
        
        #Finally, calculate distribution of CALIOP categories among the PPS Fractional category
  
                     
        cal_low_fractional = np.divide(100.0 * n_frac_low, pps_frac )
        cal_medium_fractional = np.divide(100.0 * n_frac_medium, pps_frac)
        cal_high_fractional = np.divide(100.0 * n_frac_high, pps_frac)
        
        #n_medium_yes_no = n_medium_low+n_medium_high
        #far_medium_low = np.divide(100.0 * n_medium_low, n_medium_yes_no )
        #far_medium_high = np.divide(100.0 * n_medium_high, n_medium_yes_no )
       
        #n_high_yes_no = n_high_low+n_high_medium
        #far_high_low = np.divide(100.0 * n_high_low, n_high_yes_no )
        #far_high_medium = np.divide(100.0 * n_high_medium, n_high_yes_no )
   
    
        self.Num_cma = Num_cma
        self.common_cloud_free = common_cloud_free
        self.number_excluded_cirrus_medium = number_excluded_cirrus_medium
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
        self.frac_fraction_pps_rel = frac_fraction_pps_rel
        self.low_fraction_pps_abs = low_fraction_pps_abs
        self.low_fraction_cal_abs = low_fraction_cal_abs
        self.medium_fraction_pps_abs = medium_fraction_pps_abs
        self.medium_fraction_cal_abs = medium_fraction_cal_abs
        self.high_fraction_pps_abs = high_fraction_pps_abs
        self.high_fraction_cal_abs = high_fraction_cal_abs
        self.frac_fraction_pps_abs = frac_fraction_pps_abs
        self.bias_low = bias_low
        self.bias_medium = bias_medium
        self.bias_high = bias_high
        self.pod_low = pod_low
        self.pod_medium = pod_medium
        self.pod_high = pod_high
        self.far_low = far_low
        self.far_medium = far_medium
        self.far_high = far_high
        self.hitrate = hitrate
        self.cal_low_fractional = cal_low_fractional
        self.cal_medium_fractional = cal_medium_fractional
        self.cal_high_fractional = cal_high_fractional 
        self.pps_undetected = pps_undetected
        self.pps_undetected_low = pps_undetected_low
        self.pps_undetected_medium = pps_undetected_medium
        self.pps_undetected_high = pps_undetected_high
        self.pps_false = pps_false
        self.pps_false_low = pps_false_low
        self.pps_false_medium = pps_false_medium
        self.pps_false_high = pps_false_high
        self.pps_false_frac = pps_false_frac
    def printout(self):
        lines = []
        if self.Num_cma == 0 or self.Num_ct==0:
            return lines
        lines.append( "Total pixels: %d" % self.Num_cma)
        lines.append( "Common cloud-free: %d" % self.common_cloud_free)
        lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
        lines.append("Total number of matched cloud types: %d" % self.Num_ct)
        lines.append("Some %d pixels in the category PPS-cirrus CALIPSO-clear are not excluded\n" 
                     %self.number_excluded_cirrus_medium)
        lines.append("")
        lines.append("Note: There is a separate script to get CT-statitcs for validaion report.")
        lines.append("Probability of detecting LOW, MEDIUM and HIGH: %.2f %.2f %.2f" % \
                     (self.pod_low, self.pod_medium, self.pod_high))
        lines.append("False alarm rate LOW, MEDIUM and HIGH: %.2f %.2f %.2f" % \
                     (self.far_low, self.far_medium, self.far_high))

        lines.append("Rel. Fraction LOW for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(), self.low_fraction_pps_rel, self.low_fraction_cal_rel))
        lines.append("Rel. Fraction MEDIUM for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(),self.medium_fraction_pps_rel, self.medium_fraction_cal_rel))
        lines.append("Rel. Fraction HIGH for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(),self.high_fraction_pps_rel, self.high_fraction_cal_rel))
        lines.append("Rel. Fraction FRACTIONAL for PPS: %.2f" % self.frac_fraction_pps_rel)
        lines.append("Abs. Fraction LOW for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(),self.low_fraction_pps_abs, self.low_fraction_cal_abs))
        lines.append("Abs. Fraction MEDIUM for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(),self.medium_fraction_pps_abs, self.medium_fraction_cal_abs))
        lines.append("Abs. Fraction HIGH for PPS and for %s: %.2f %.2f" % \
                     (self.truth_sat.upper(),self.high_fraction_pps_abs, self.high_fraction_cal_abs))
        lines.append("Abs. Fraction FRACTIONAL for PPS: %.2f" % self.frac_fraction_pps_abs)
        lines.append("Bias Low: %.2f" % self.bias_low)
        lines.append("Bias Medium: %.2f" % self.bias_medium)
        lines.append("Bias High: %.2f" % self.bias_high)
        lines.append("POD (Low, Medium, High): %.2f %.2f %.2f" % \
                     (self.pod_low, self.pod_medium, self.pod_high))
        lines.append("FAR (Low, Medium, High): %.2f %.2f %.2f" % \
                     (self.far_low, self.far_medium, self.far_high))
        lines.append("HR: %.3f " % (self.hitrate))
        lines.append("%s parts of Fractional (Low, Medium, High): %.2f %.2f %.2f" % \
                     (self.truth_sat.upper(), self.cal_low_fractional, self.cal_medium_fractional, self.cal_high_fractional))
        lines.append("This is really CMA buisness here!")
        lines.append("Missclassified cloudy: %d (low:%d, medium:%d, high:%d, frac:%d)" % \
                     (self.pps_false, self.pps_false_low,
                      self.pps_false_medium, self.pps_false_high,
                      self.pps_false_frac))
        lines.append("Missclassified clear: %d (low:%d, medium:%d, high:%d)" % \
                     (self.pps_undetected, self.pps_undetected_low,
                      self.pps_undetected_medium, self.pps_undetected_high))
        lines.append("")
        
        return lines


if __name__ == "__main__":
    stats = CloudTypeStats()

