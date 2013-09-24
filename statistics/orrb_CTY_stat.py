# Program orrb_CTY_stat

# This program calculates basic statistics for the cloud type (CTY) product for
# each month

import string
import math
from numpy import divide
from orrb_stat_class import OrrbStats

# -----------------------------------------------------
class CloudTypeStats(OrrbStats):
    
    def __init__(self, results_files=None, cfc_stats=None):
        self.cfc_stats = cfc_stats
        OrrbStats.__init__(self, results_files)
    
    def do_stats(self):
        OrrbStats.do_stats(self)
        
        try:
            CFC_CALIOP = self.cfc_stats.mean_CFC_cal
            CFC_MEAN_ERROR = self.cfc_stats.bias_cal_perc
            Totpix = self.cfc_stats.samples_cal
        except (AttributeError, TypeError):
            from orrb_CFC_stat import CloudFractionStats
            cfc_stats = CloudFractionStats(self.results_files)
            CFC_CALIOP = cfc_stats.mean_CFC_cal
            CFC_MEAN_ERROR = cfc_stats.bias_cal_perc
            Totpix = cfc_stats.samples_cal
        CFC_PPS = CFC_CALIOP + CFC_MEAN_ERROR
        
        n_low_low = 0
        n_low_medium = 0
        n_low_high = 0
        n_medium_low = 0
        n_medium_medium = 0
        n_medium_high = 0
        n_high_low = 0
        n_high_medium = 0
        n_high_high = 0
        n_frac_low = 0
        n_frac_medium = 0
        n_frac_high = 0
    
        n_clear_low = 0
        n_clear_medium = 0
        n_clear_high = 0
        n_low_clear = 0
        n_medium_clear = 0
        n_high_clear = 0
        n_frac_clear = 0
    
        scenes = len(self.results_files)
    
        for datafile in self.results_files:
            current_datafile = open(datafile, "r")
            datalist = current_datafile.readlines()
            #print "Datafile: ", datafile
            cal_data = string.split(datalist[12])
            cal_data_missed = string.split(datalist[14])
    
            # Accumulate CALIOP statistics
    
            n_low_low = n_low_low + int(cal_data[4])
            n_low_medium = n_low_medium + int(cal_data[5])
            n_low_high = n_low_high + int(cal_data[6])
            n_medium_low = n_medium_low + int(cal_data[7])
            n_medium_medium = n_medium_medium + int(cal_data[8])
            n_medium_high = n_medium_high + int(cal_data[9])
            n_high_low = n_high_low + int(cal_data[10])
            n_high_medium = n_high_medium + int(cal_data[11])
            n_high_high = n_high_high + int(cal_data[12])
            n_frac_low = n_frac_low + int(cal_data[13])
            n_frac_medium = n_frac_medium + int(cal_data[14])
            n_frac_high = n_frac_high + int(cal_data[15])
    
            n_clear_low = n_clear_low + int(cal_data_missed[5]) 
            n_clear_medium = n_clear_medium + int(cal_data_missed[6]) 
            n_clear_high = n_clear_high + int(cal_data_missed[7])
            n_low_clear = n_low_clear + int(cal_data_missed[8])
            n_medium_clear = n_medium_clear + int(cal_data_missed[9])
            n_high_clear = n_high_clear + int(cal_data_missed[10])
            n_frac_clear = n_frac_clear + int(cal_data_missed[11])
    
            current_datafile.close()
            
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
    
        common_cloud_free = Totpix-samples_low_cal-samples_medium_cal-samples_high_cal
        
        # numpy.divide handles potential division by zero
        pod_low = divide(float(n_low_low),float(samples_low_cal))
        far_low = 1-pod_low
        pod_medium = divide(float(n_medium_medium),float(samples_medium_cal))
        far_medium = 1-pod_medium
        pod_high = divide(float(n_high_high),float(samples_high_cal))
        far_high = 1-pod_high
    
        low_fraction_pps_rel = divide(float(n_low_low+n_low_medium+n_low_high),float(samples_tot))*100.0
        low_fraction_cal_rel = divide(float(samples_low_cal),float(samples_tot))*100.0
        medium_fraction_pps_rel = divide(float(n_medium_medium+n_medium_low+n_medium_high),float(samples_tot))*100.0
        medium_fraction_cal_rel = divide(float(samples_medium_cal),float(samples_tot))*100.0
        high_fraction_pps_rel = divide(float(n_high_high+n_high_low+n_high_medium),float(samples_tot))*100.0
        high_fraction_cal_rel = divide(float(samples_high_cal),float(samples_tot))*100.0
        frac_fraction_pps_rel = divide(float(n_frac_low+n_frac_medium+n_frac_high),float(samples_tot))*100.0
        low_fraction_pps_abs = divide(float(n_low_low+n_low_medium+n_low_high),float(samples_tot))*CFC_PPS
        low_fraction_cal_abs = divide(float(samples_low_cal),float(samples_tot))*CFC_CALIOP
        medium_fraction_pps_abs = divide(float(n_medium_medium+n_medium_low+n_medium_high),float(samples_tot))*CFC_PPS
        medium_fraction_cal_abs = divide(float(samples_medium_cal),float(samples_tot))*CFC_CALIOP
        high_fraction_pps_abs = divide(float(n_high_high+n_high_low+n_high_medium),float(samples_tot))*CFC_PPS
        high_fraction_cal_abs = divide(float(samples_high_cal),float(samples_tot))*CFC_CALIOP
        frac_fraction_pps_abs = divide(float(n_frac_low+n_frac_medium+n_frac_high),float(samples_tot))*CFC_PPS
        
        bias_low=low_fraction_pps_abs-low_fraction_cal_abs
        bias_low_noperc = bias_low/100.0
        bias_medium=medium_fraction_pps_abs-medium_fraction_cal_abs
        bias_medium_noperc = bias_medium/100.0
        bias_high=high_fraction_pps_abs-high_fraction_cal_abs
        bias_high_noperc = bias_high/100.0
        
        #bc-RMS calculations ===============================================================================
        
        n_low_no_no = n_medium_medium+n_medium_high+n_high_medium+n_high_high+n_frac_medium+n_frac_high+\
                        (Totpix-samples_low_cal-samples_medium_cal-samples_high_cal)
        n_low_yes_yes = n_low_low
        n_low_yes_no = n_low_medium+n_low_high+n_low_clear
        n_low_no_yes = n_medium_low+n_high_low+n_clear_low+n_frac_low
        square_sum_cal_low = float(n_low_yes_yes+n_low_no_no)*bias_low_noperc*bias_low_noperc + \
                                n_low_no_yes*(-1.0-bias_low_noperc)*(-1.0-bias_low_noperc) + \
                                n_low_yes_no*(1.0-bias_low_noperc)*(1.0-bias_low_noperc)
        rms_low = 100.0*math.sqrt(square_sum_cal_low/(Totpix-1))
    
        n_medium_no_no = n_low_low+n_low_high+n_high_low+n_high_high+n_frac_low+n_frac_high+\
                        (Totpix-samples_low_cal-samples_medium_cal-samples_high_cal)
        n_medium_yes_yes = n_medium_medium
        n_medium_yes_no = n_medium_low+n_medium_high+n_medium_clear
        n_medium_no_yes = n_low_medium+n_high_medium+n_clear_medium+n_frac_medium
        square_sum_cal_medium = float(n_medium_yes_yes+n_medium_no_no)*bias_medium_noperc*bias_medium_noperc+ \
                                n_medium_no_yes*(-1.0-bias_medium_noperc)*(-1.0-bias_medium_noperc)+ \
                                n_medium_yes_no*(1.0-bias_medium_noperc)*(1.0-bias_medium_noperc)
        rms_medium = 100.0*math.sqrt(square_sum_cal_medium/(Totpix-1))
    
        n_high_no_no = n_medium_medium+n_medium_low+n_low_low+n_low_medium+n_frac_medium+n_frac_low+\
                        (Totpix-samples_low_cal-samples_medium_cal-samples_high_cal)
        n_high_yes_yes = n_high_high
        n_high_yes_no = n_high_medium+n_high_low+n_high_clear
        n_high_no_yes = n_medium_high+n_low_high+n_clear_high+n_frac_high
        square_sum_cal_high = float(n_high_yes_yes+n_high_no_no)*bias_high_noperc*bias_high_noperc + \
                                n_high_no_yes*(-1.0-bias_high_noperc)*(-1.0-bias_high_noperc) + \
                                n_high_yes_no*(1.0-bias_high_noperc)*(1.0-bias_high_noperc)
        rms_high = 100.0*math.sqrt(square_sum_cal_high/(Totpix-1))
    
        #=============================================================================================
        
        #POD,FAR,HR and KSS calculations =============================================================
        
        pod_low = divide(100.0*n_low_yes_yes,samples_low_cal )
        pod_medium = divide(100.0*n_medium_yes_yes,samples_medium_cal )
        pod_high = divide(100.0*n_high_yes_yes,samples_high_cal )
    
        far_low = divide(100.0*(n_low_yes_no),(n_low_yes_yes+n_low_yes_no) )
        far_medium = divide(100.0*(n_medium_yes_no),(n_medium_yes_yes+n_medium_yes_no) )
        far_high = divide(100.0*(n_high_yes_no),(n_high_yes_yes+n_high_yes_no) )
    
        hitrate_low = divide(1.0*(n_low_yes_yes+n_low_no_no),Totpix )
        hitrate_medium = divide(1.0*(n_medium_yes_yes+n_low_no_no),Totpix )
        hitrate_high = divide(1.0*(n_high_yes_yes+n_high_no_no),Totpix )
    
        kuipers_low = divide(1.0*(n_low_no_no*n_low_yes_yes-n_low_yes_no*n_low_no_yes), \
                        ((n_low_no_no+n_low_no_yes)*(n_low_yes_no+n_low_yes_yes)) )
        kuipers_medium = divide(1.0*(n_medium_no_no*n_medium_yes_yes-n_medium_yes_no*n_medium_no_yes), \
                        ((n_medium_no_no+n_medium_no_yes)*(n_medium_yes_no+n_medium_yes_yes)) )
        kuipers_high = divide(1.0*(n_high_no_no*n_high_yes_yes-n_high_yes_no*n_high_no_yes), \
                        ((n_high_no_no+n_high_no_yes)*(n_high_yes_no+n_high_yes_yes)) )
    
        #=============================================================================================
        
        #Finally, calculate distribution of CALIOP categories among the PPS Fractional category
    
        cal_low_fractional = divide(100.0*n_frac_low,(n_frac_low+n_frac_medium+n_frac_high) )
        cal_medium_fractional = divide(100.0*n_frac_medium,(n_frac_low+n_frac_medium+n_frac_high) )
        cal_high_fractional = divide(100.0*n_frac_high,(n_frac_low+n_frac_medium+n_frac_high) )
        
        #=============================================================================================
        
        #Special analysis: Check composition of false alarm categories for medium and high clouds
    
        n_medium_yes_no = n_medium_low+n_medium_high+n_medium_clear
        far_medium_low = divide(100.0*n_medium_low,n_medium_yes_no )
        far_medium_high = divide(100.0*n_medium_high,n_medium_yes_no )
        far_medium_clear = divide(100.0*n_medium_clear,n_medium_yes_no )
    
        n_high_yes_no = n_high_low+n_high_medium+n_high_clear
        far_high_low = divide(100.0*n_high_low,n_high_yes_no )
        far_high_medium = divide(100.0*n_high_medium,n_high_yes_no )
        far_high_clear = divide(100.0*n_high_clear,n_high_yes_no  )       
    
        self.Totpix = Totpix
        self.common_cloud_free = common_cloud_free
        self.scenes = scenes
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
        lines.append( "Total pixels: %s" % self.Totpix)
        lines.append( "Common cloud-free: %s" % self.common_cloud_free)
        lines.append("Total number of matched scenes is: %s" % self.scenes)
        lines.append("Total number of matched cloud types: %s" % self.samples_tot)
        lines.append("")
        lines.append("Probability of detecting LOW, MEDIUM and HIGH: %f %f %f" % \
                     (self.pod_low, self.pod_medium, self.pod_high))
        lines.append("False alarm rate LOW, MEDIUM and HIGH: %f %f %f" % \
                     (self.far_low, self.far_medium, self.far_high))
        lines.append("Rel. Fraction LOW for PPS and for CALIOP: %f %f" % \
                     (self.low_fraction_pps_rel, self.low_fraction_cal_rel))
        lines.append("Rel. Fraction MEDIUM for PPS and for CALIOP: %f %f" % \
                     (self.medium_fraction_pps_rel, self.medium_fraction_cal_rel))
        lines.append("Rel. Fraction HIGH for PPS and for CALIOP: %f %f" % \
                     (self.high_fraction_pps_rel, self.high_fraction_cal_rel))
        lines.append("Rel. Fraction FRACTIONAL for PPS: %f" % self.frac_fraction_pps_rel)
        lines.append("Abs. Fraction LOW for PPS and for CALIOP: %f %f" % \
                     (self.low_fraction_pps_abs, self.low_fraction_cal_abs))
        lines.append("Abs. Fraction MEDIUM for PPS and for CALIOP: %f %f" % \
                     (self.medium_fraction_pps_abs, self.medium_fraction_cal_abs))
        lines.append("Abs. Fraction HIGH for PPS and for CALIOP: %f %f" % \
                     (self.high_fraction_pps_abs, self.high_fraction_cal_abs))
        lines.append("Abs. Fraction FRACTIONAL for PPS: %f" % self.frac_fraction_pps_abs)
        lines.append("Bias Low: %f" % self.bias_low)
        lines.append("Bias Medium: %f" % self.bias_medium)
        lines.append("Bias High: %f" % self.bias_high)
        lines.append("Bc-RMS Low: %f" % self.rms_low)
        lines.append("Bc-RMS Medium: %f" % self.rms_medium)
        lines.append("Bc_RMS High: %f" % self.rms_high)
        lines.append("POD (Low,Medium,High): %f %f %f" % \
                     (self.pod_low, self.pod_medium, self.pod_high))
        lines.append("FAR (Low,Medium,High): %f %f %f" % \
                     (self.far_low, self.far_medium, self.far_high))
        lines.append("FAR Medium fraction Low,High,Clear: %f %f %f" % \
                     (self.far_medium_low, self.far_medium_high, self.far_medium_clear))
        lines.append("FAR High fraction Low,Medium,Clear: %f %f %f" % \
                     (self.far_high_low, self.far_high_medium, self.far_high_clear) )
        lines.append("HR (Low,Medium,High): %f %f %f" % \
                     (self.hitrate_low, self.hitrate_medium, self.hitrate_high))
        lines.append("KSS (Low,Medium,High): %f %f %f" % \
                     (self.kuipers_low, self.kuipers_medium, self.kuipers_high))
        lines.append("CALIOP parts of Fractional (Low, Medium, High): %f %f %f" % \
                     (self.cal_low_fractional, self.cal_medium_fractional, self.cal_high_fractional))
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
    stats.old_interface(modes=['BASIC'], output_file_desc="cty_results_summary")
