# Program orrb_CTH_stat

# This program calculates basic statistics for the cloud top (CTH) product for
# each month

import string
import math
import pdb
import config
from orrb_stat_class import OrrbStats
import numpy as np
# -----------------------------------------------------

def bias_corrected_rms(rms, bias, N):
    # ie formula should be bcRMS= sqrt(RMS^2-c*bias^2), 
    # where c=N/(N-1). However our Ns are usually large
    if N < 2:
        print "Warning too few elements to calculate bc-RMS"
        return -9
    cnn1 = N/(N-1)
    return np.sqrt(rms*rms-cnn1*bias*bias)
    

class CloudTopStats(OrrbStats):

    def do_stats(self):
        OrrbStats.do_stats(self)
        #for key in self.ac_data.keys():
        #    print key
        
        cal_all_samples = self.ac_data["cal_all_samples"] 
        cal_low_samples = self.ac_data["cal_low_samples"] 
        cal_medium_samples = self.ac_data["cal_medium_samples"] 
        cal_high_samples = self.ac_data["cal_high_samples"] 
        mean_error_csa_sum = self.ac_data["mean_error_csa_sum"] 
        mean_error_cal_all_sum = self.ac_data["mean_error_cal_all_sum"] 
        mean_error_cal_low_sum = self.ac_data["mean_error_cal_low_sum"] 
        mean_error_cal_medium_sum = self.ac_data["mean_error_cal_medium_sum"] 
        mean_error_cal_high_sum = self.ac_data["mean_error_cal_high_sum"] 
        rms_error_csa_sum = self.ac_data["rms_error_csa_sum"] 
        rms_error_cal_all_sum = self.ac_data["rms_error_cal_all_sum"] 
        rms_error_cal_low_sum = self.ac_data["rms_error_cal_low_sum"] 
        rms_error_cal_medium_sum = self.ac_data["rms_error_cal_medium_sum"] 
        rms_error_cal_high_sum = self.ac_data["rms_error_cal_high_sum"]

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

        for tc in cal_all_samples.keys():
        # numpy.divide handles potential division by zero
            self.retrieval_rate_all[tc] = cal_all_samples[tc]/(cal_all_samples[tc] +
                                                           n_missed_ctth_all[tc])
            self.estimate_pod_cloudy_all[tc] = cal_all_samples[tc]/(cal_all_samples[tc] +
                                                           n_missed_cma_all[tc])
            self.bias_cal_all[tc] = np.divide(mean_error_cal_all_sum[tc], 
                                              cal_all_samples[tc])
            self.rms_cal_all[tc] = math.sqrt(np.divide(rms_error_cal_all_sum[tc], 
                                                       cal_all_samples[tc]))
            self.bcrms_cal_all[tc] = bias_corrected_rms(self.rms_cal_all[tc], 
                                                        self.bias_cal_all[tc], 
                                                        cal_all_samples[tc])
        for tc in cal_low_samples.keys():
            self.retrieval_rate_low[tc] = cal_low_samples[tc]/(cal_low_samples[tc] +
                                                           n_missed_ctth_low[tc])
            self.retrieval_rate_medium[tc] = cal_medium_samples[tc]/(cal_medium_samples[tc] +
                                                           n_missed_ctth_medium[tc])
            self.retrieval_rate_high[tc] = cal_high_samples[tc]/(cal_high_samples[tc] +
                                                           n_missed_ctth_high[tc])

            self.estimate_pod_cloudy_low[tc] = cal_low_samples[tc]/(cal_low_samples[tc] +
                                                           n_missed_cma_low[tc])
            self.estimate_pod_cloudy_medium[tc] = cal_medium_samples[tc]/(cal_medium_samples[tc] +
                                                           n_missed_cma_medium[tc])
            self.estimate_pod_cloudy_high[tc] = cal_high_samples[tc]/(cal_high_samples[tc] +
                                                           n_missed_cma_high[tc])

            self.bias_cal_low[tc] = np.divide(mean_error_cal_low_sum[tc], 
                                         cal_low_samples[tc])
            self.bias_cal_medium[tc] = np.divide(mean_error_cal_medium_sum[tc], 
                                            cal_medium_samples[tc])
            self.bias_cal_high[tc] = np.divide(mean_error_cal_high_sum[tc], 
                                          cal_high_samples[tc])
    

            self.rms_cal_low[tc] = math.sqrt(np.divide(rms_error_cal_low_sum[tc], 
                                                  cal_low_samples[tc]))
            self.rms_cal_medium[tc]= math.sqrt(np.divide(rms_error_cal_medium_sum[tc], 
                                                     cal_medium_samples[tc]))
            self.rms_cal_high[tc] = math.sqrt(np.divide(rms_error_cal_high_sum[tc], 
                                                   cal_high_samples[tc]))


            self.bcrms_cal_low[tc] = bias_corrected_rms( self.rms_cal_low[tc], 
                                                     self.bias_cal_low[tc], 
                                                    cal_low_samples[tc])
            self.bcrms_cal_medium[tc] = bias_corrected_rms( self.rms_cal_medium[tc], 
                                                        self.bias_cal_medium[tc], 
                                                       cal_medium_samples[tc])
            self.bcrms_cal_high[tc] = bias_corrected_rms( self.rms_cal_high[tc], 
                                                      self.bias_cal_high[tc], 
                                                     cal_high_samples[tc])


    
    def printout(self):            
        lines = []
        for tc in sorted(self.cal_low_samples.keys()):
            lines.append("========== Cloud top height ===========")
            lines.append("====== %s ======"%(tc))
            lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
            lines.append("")
            lines.append("Imager retrieval rate: %3.2f" % self.retrieval_rate_all[tc])
            lines.append("Total number of %s matched cloudtops: %d" %( self.truth_sat.upper(), self.cal_all_samples[tc]))
            lines.append("Number of %s matched low cloudtops: %d" % (self.truth_sat.upper(), self.cal_low_samples[tc]))
            lines.append("Number of %s matched medium cloudtops: %d" %( self.truth_sat.upper(), self.cal_medium_samples[tc]))
            lines.append("Number of %s matched high cloudtops: %d" %( self.truth_sat.upper(), self.cal_high_samples[tc]))
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
            lines.append("")
        for tc in  self.cal_all_samples.keys():   
            if tc in self.cal_low_samples.keys():
                #already go these results printed
                continue
            else:
                lines.append("========== Cloud top height ===========")
                lines.append("====== %s ======"%(tc))
                lines.append("Total number of matched scenes is: %s" % self.ac_data["scenes"])
                lines.append("")
                lines.append("Imager retrieval rate: %3.2f" % self.retrieval_rate_all[tc])
                lines.append("Total number of %s matched cloudtops: %d" %( self.truth_sat.upper(), self.cal_all_samples[tc]))
                lines.append("Mean error total cases: %.0f" % self.bias_cal_all[tc])
                lines.append("RMS error total cases: %.0f" % self.rms_cal_all[tc])
                lines.append("bc-RMS error total cases: %.0f" % self.bcrms_cal_all[tc])
                lines.append("Imager estimated POD-cloudy total-level cases: %3.2f" % self.estimate_pod_cloudy_all[tc]) 
                lines.append("")  
        return lines


if __name__ == "__main__":
    stats = CloudTopStats()

