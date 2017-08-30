# Program orrb_CTY_stat

# This program calculates basic statistics for the cloud type (CTY) product for
# each month

import string
import math
import numpy as np
from orrb_stat_class import OrrbStats

# -----------------------------------------------------
class CloudTypeStats(OrrbStats):
    
    def do_stats(self):
        OrrbStats.do_stats(self)

        Num = self.ac_data["Num"]
        CFC_TRUTH = np.divide(100.0*(self.ac_data["n_cloudy_cloudy_cal"]+
                                        self.ac_data["n_cloudy_clear_cal"]), 
                              Num)
        CFC_PPS = np.divide(100.0*(self.ac_data["n_cloudy_cloudy_cal"]+
                                   self.ac_data["n_clear_cloudy_cal"]), 
                            Num)

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
        n_clear_low = self.ac_data["n_clear_low"] 
        n_clear_medium = self.ac_data["n_clear_medium"] 
        n_clear_high = self.ac_data["n_clear_high"] 
        n_low_clear = self.ac_data["n_low_clear"] 
        n_medium_clear = self.ac_data["n_medium_clear"] 
        n_high_clear = self.ac_data["n_high_clear"] 
        n_frac_clear = self.ac_data["n_frac_clear"] 
            
        samples_pps_missed = n_clear_low + n_clear_medium + n_clear_high
        samples_pps_missed_low = n_clear_low
        samples_pps_missed_medium = n_clear_medium
        samples_pps_missed_high = n_clear_high
        samples_pps_misclass = n_low_clear + n_medium_clear + n_high_clear + n_frac_clear
        samples_pps_misclass_low = n_low_clear
        samples_pps_misclass_medium = n_medium_clear
        samples_pps_misclass_high = n_high_clear
        samples_pps_misclass_frac = n_frac_clear
    
        samples_tot = n_low_low + n_low_medium + n_low_high + n_medium_low + n_medium_medium + \
                        n_medium_high + n_high_low + n_high_medium + n_high_high + n_frac_low + \
                        n_frac_medium +n_frac_high
        samples_low_cal = n_low_low + n_medium_low + n_high_low + n_frac_low
        samples_medium_cal = n_low_medium + n_medium_medium + n_high_medium + n_frac_medium
        samples_high_cal = n_low_high + n_medium_high + n_high_high + n_frac_high
    
        common_cloud_free = self.ac_data["n_clear_clear_cal"]
        
        # numpy.divide handles potential division by zero
        pod_low = np.divide(float(n_low_low + n_frac_low),float(samples_low_cal))
        far_low = 1-pod_low
        pod_medium = np.divide(float(n_medium_medium),float(samples_medium_cal))
        far_medium = 1-pod_medium
        pod_high = np.divide(float(n_high_high),float(samples_high_cal))
        far_high = 1-pod_high
    
        low_fraction_pps_rel = np.divide(float(n_low_low+n_low_medium+n_low_high),float(samples_tot))*100.0
        low_fraction_cal_rel = np.divide(float(samples_low_cal),float(samples_tot))*100.0
        medium_fraction_pps_rel = np.divide(float(n_medium_medium+n_medium_low+n_medium_high),float(samples_tot))*100.0
        medium_fraction_cal_rel = np.divide(float(samples_medium_cal),float(samples_tot))*100.0
        high_fraction_pps_rel = np.divide(float(n_high_high+n_high_low+n_high_medium),float(samples_tot))*100.0
        high_fraction_cal_rel = np.divide(float(samples_high_cal),float(samples_tot))*100.0
        frac_fraction_pps_rel = np.divide(float(n_frac_low+n_frac_medium+n_frac_high),float(samples_tot))*100.0
        low_fraction_pps_abs = np.divide(float(n_low_low+n_low_medium+n_low_high),float(samples_tot))*CFC_PPS
        low_fraction_cal_abs = np.divide(float(samples_low_cal),float(samples_tot))*CFC_TRUTH
        medium_fraction_pps_abs = np.divide(float(n_medium_medium+n_medium_low+n_medium_high),float(samples_tot))*CFC_PPS
        medium_fraction_cal_abs = np.divide(float(samples_medium_cal),float(samples_tot))*CFC_TRUTH
        high_fraction_pps_abs = np.divide(float(n_high_high+n_high_low+n_high_medium),float(samples_tot))*CFC_PPS
        high_fraction_cal_abs = np.divide(float(samples_high_cal),float(samples_tot))*CFC_TRUTH
        frac_fraction_pps_abs = np.divide(float(n_frac_low+n_frac_medium+n_frac_high),float(samples_tot))*CFC_PPS
        
        bias_low=low_fraction_pps_abs-low_fraction_cal_abs
        bias_low_noperc = bias_low/100.0
        bias_medium=medium_fraction_pps_abs-medium_fraction_cal_abs
        bias_medium_noperc = bias_medium/100.0
        bias_high=high_fraction_pps_abs-high_fraction_cal_abs
        bias_high_noperc = bias_high/100.0
        
        #bc-RMS calculations ===============================================================================
        
        n_low_no_no = n_medium_medium+n_medium_high+n_high_medium+n_high_high+n_frac_medium+n_frac_high+\
                        (Num-samples_low_cal-samples_medium_cal-samples_high_cal)
        n_low_yes_yes = n_low_low + n_frac_low
        n_low_yes_no = n_low_medium+n_low_high+n_low_clear
        n_low_no_yes = n_medium_low+n_high_low+n_clear_low
        square_sum_cal_low = float(n_low_yes_yes+n_low_no_no)*bias_low_noperc*bias_low_noperc + \
                                n_low_no_yes*(-1.0-bias_low_noperc)*(-1.0-bias_low_noperc) + \
                                n_low_yes_no*(1.0-bias_low_noperc)*(1.0-bias_low_noperc)
        rms_low = 100.0*math.sqrt(square_sum_cal_low/(Num-1))
    
        n_medium_no_no = n_low_low+n_low_high+n_high_low+n_high_high+n_frac_low+n_frac_high+\
                        (Num-samples_low_cal-samples_medium_cal-samples_high_cal)
        n_medium_yes_yes = n_medium_medium
        n_medium_yes_no = n_medium_low+n_medium_high+n_medium_clear
        n_medium_no_yes = n_low_medium+n_high_medium+n_clear_medium+n_frac_medium
        square_sum_cal_medium = float(n_medium_yes_yes+n_medium_no_no)*bias_medium_noperc*bias_medium_noperc+ \
                                n_medium_no_yes*(-1.0-bias_medium_noperc)*(-1.0-bias_medium_noperc)+ \
                                n_medium_yes_no*(1.0-bias_medium_noperc)*(1.0-bias_medium_noperc)
        rms_medium = 100.0*math.sqrt(square_sum_cal_medium/(Num-1))
    
        n_high_no_no = n_medium_medium+n_medium_low+n_low_low+n_low_medium+n_frac_medium+n_frac_low+\
                        (Num-samples_low_cal-samples_medium_cal-samples_high_cal)
        n_high_yes_yes = n_high_high
        n_high_yes_no = n_high_medium+n_high_low+n_high_clear
        n_high_no_yes = n_medium_high+n_low_high+n_clear_high+n_frac_high
        square_sum_cal_high = float(n_high_yes_yes+n_high_no_no)*bias_high_noperc*bias_high_noperc + \
                                n_high_no_yes*(-1.0-bias_high_noperc)*(-1.0-bias_high_noperc) + \
                                n_high_yes_no*(1.0-bias_high_noperc)*(1.0-bias_high_noperc)
        rms_high = 100.0*math.sqrt(square_sum_cal_high/(Num-1))
    
        #=============================================================================================
        
        #POD,FAR,HR and KSS calculations =============================================================
        
        pod_low = np.divide(100.0*n_low_yes_yes,samples_low_cal )
        pod_medium = np.divide(100.0*n_medium_yes_yes,samples_medium_cal )
        pod_high = np.divide(100.0*n_high_yes_yes,samples_high_cal )
    
        far_low = np.divide(100.0*(n_low_yes_no),(n_low_yes_yes+n_low_yes_no) )
        far_medium = np.divide(100.0*(n_medium_yes_no),(n_medium_yes_yes+n_medium_yes_no) )
        far_high = np.divide(100.0*(n_high_yes_no),(n_high_yes_yes+n_high_yes_no) )
    
        hitrate_low = np.divide(1.0*(n_low_yes_yes+n_low_no_no),Num )
        hitrate_medium = np.divide(1.0*(n_medium_yes_yes+n_low_no_no),Num )
        hitrate_high = np.divide(1.0*(n_high_yes_yes+n_high_no_no),Num )
    
        kuipers_low = np.divide(1.0*(n_low_no_no*n_low_yes_yes-n_low_yes_no*n_low_no_yes), \
                        ((n_low_no_no+n_low_no_yes)*(n_low_yes_no+n_low_yes_yes)) )
        kuipers_medium = np.divide(1.0*(n_medium_no_no*n_medium_yes_yes-n_medium_yes_no*n_medium_no_yes), \
                        ((n_medium_no_no+n_medium_no_yes)*(n_medium_yes_no+n_medium_yes_yes)) )
        kuipers_high = np.divide(1.0*(n_high_no_no*n_high_yes_yes-n_high_yes_no*n_high_no_yes), \
                        ((n_high_no_no+n_high_no_yes)*(n_high_yes_no+n_high_yes_yes)) )
    
        #=============================================================================================
        
        #Finally, calculate distribution of CALIOP categories among the PPS Fractional category
    
        cal_low_fractional = np.divide(100.0*n_frac_low,(n_frac_low+n_frac_medium+n_frac_high) )
        cal_medium_fractional = np.divide(100.0*n_frac_medium,(n_frac_low+n_frac_medium+n_frac_high) )
        cal_high_fractional = np.divide(100.0*n_frac_high,(n_frac_low+n_frac_medium+n_frac_high) )
        
        #=============================================================================================
        
        #Special analysis: Check composition of false alarm categories for medium and high clouds
    
        n_medium_yes_no = n_medium_low+n_medium_high+n_medium_clear
        far_medium_low = np.divide(100.0*n_medium_low,n_medium_yes_no )
        far_medium_high = np.divide(100.0*n_medium_high,n_medium_yes_no )
        far_medium_clear = np.divide(100.0*n_medium_clear,n_medium_yes_no )
    
        n_high_yes_no = n_high_low+n_high_medium+n_high_clear
        far_high_low = np.divide(100.0*n_high_low,n_high_yes_no )
        far_high_medium = np.divide(100.0*n_high_medium,n_high_yes_no )
        far_high_clear = np.divide(100.0*n_high_clear,n_high_yes_no  )       
    
        self.Num = Num
        self.common_cloud_free = common_cloud_free
        self.samples_tot = samples_tot
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
        self.rms_low = rms_low
        self.rms_medium = rms_medium
        self.rms_high = rms_high
        self.pod_low = pod_low
        self.pod_medium = pod_medium
        self.pod_high = pod_high
        self.far_low = far_low
        self.far_medium = far_medium
        self.far_high = far_high
        self.far_medium_low = far_medium_low
        self.far_medium_high = far_medium_high
        self.far_medium_clear = far_medium_clear
        self.far_high_low = far_high_low
        self.far_high_medium = far_high_medium
        self.far_high_clear = far_high_clear 
        self.hitrate_low = hitrate_low
        self.hitrate_medium = hitrate_medium
        self.hitrate_high = hitrate_high
        self.kuipers_low = kuipers_low
        self.kuipers_medium = kuipers_medium
        self.kuipers_high = kuipers_high
        self.cal_low_fractional = cal_low_fractional
        self.cal_medium_fractional = cal_medium_fractional
        self.cal_high_fractional = cal_high_fractional 
        self.samples_pps_missed = samples_pps_missed
        self.samples_pps_missed_low = samples_pps_missed_low
        self.samples_pps_missed_medium = samples_pps_missed_medium
        self.samples_pps_missed_high = samples_pps_missed_high
        self.samples_pps_misclass = samples_pps_misclass
        self.samples_pps_misclass_low = samples_pps_misclass_low
        self.samples_pps_misclass_medium = samples_pps_misclass_medium
        self.samples_pps_misclass_high = samples_pps_misclass_high
        self.samples_pps_misclass_frac = samples_pps_misclass_frac
    
    def printout(self):
        lines = []
        if self.Num == 0:
            return lines
        lines.append( "Total pixels: %s" % self.Num)
        lines.append( "Common cloud-free: %s" % self.common_cloud_free)
        lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
        lines.append("Total number of matched cloud types: %s" % self.samples_tot)
        lines.append("")
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
        lines.append("Bc-RMS Low: %.2f" % self.rms_low)
        lines.append("Bc-RMS Medium: %.2f" % self.rms_medium)
        lines.append("Bc_RMS High: %.2f" % self.rms_high)
        lines.append("POD (Low,Medium,High): %.2f %.2f %.2f" % \
                     (self.pod_low, self.pod_medium, self.pod_high))
        lines.append("FAR (Low,Medium,High): %.2f %.2f %.2f" % \
                     (self.far_low, self.far_medium, self.far_high))
        lines.append("FAR Medium fraction Low,High,Clear: %.2f %.2f %.2f" % \
                     (self.far_medium_low, self.far_medium_high, self.far_medium_clear))
        lines.append("FAR High fraction Low,Medium,Clear: %.2f %.2f %.2f" % \
                     (self.far_high_low, self.far_high_medium, self.far_high_clear) )
        lines.append("HR (Low,Medium,High): %.3f %.3f %.3f" % \
                     (self.hitrate_low, self.hitrate_medium, self.hitrate_high))
        lines.append("KSS (Low,Medium,High): %.3f %.3f %.3f" % \
                     (self.kuipers_low, self.kuipers_medium, self.kuipers_high))
        lines.append("%s parts of Fractional (Low, Medium, High): %.2f %.2f %.2f" % \
                     (self.truth_sat.upper(), self.cal_low_fractional, self.cal_medium_fractional, self.cal_high_fractional))
        lines.append("Missclassified clear: %d (low:%d, medium:%d, high:%d)" % \
                     (self.samples_pps_missed, self.samples_pps_missed_low,
                      self.samples_pps_missed_medium, self.samples_pps_missed_high))
        lines.append("Missclassified cloudy: %d (low:%d, medium:%d, high:%d, frac:%d)" % \
                     (self.samples_pps_misclass, self.samples_pps_misclass_low,
                      self.samples_pps_misclass_medium, self.samples_pps_misclass_high,
                      self.samples_pps_misclass_frac))
        lines.append("")
        
        return lines


if __name__ == "__main__":
    stats = CloudTypeStats()

