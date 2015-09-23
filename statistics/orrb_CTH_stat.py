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

    def __init__(self, results_files=None, cal_ctth_rows=None, note=None):
        self.cal_ctth_rows = cal_ctth_rows
        self.note = note
        OrrbStats.__init__(self, results_files)

    def do_stats(self):
        OrrbStats.do_stats(self)
        
        from numpy import divide
        
        scenes = len(self.results_files)
        
        # First calculate total number of cases
        csa_samples = 0
        cal_all_samples = 0
        cal_low_samples = 0
        cal_medium_samples = 0
        cal_high_samples = 0
        mean_error_csa_sum = 0.0
        mean_error_cal_all_sum = 0.0
        mean_error_cal_low_sum = 0.0
        mean_error_cal_medium_sum = 0.0
        mean_error_cal_high_sum = 0.0
        rms_error_csa_sum = 0.0
        rms_error_cal_all_sum = 0.0
        rms_error_cal_low_sum = 0.0
        rms_error_cal_medium_sum = 0.0
        rms_error_cal_high_sum = 0.0
        for datafile in self.results_files:
            current_datafile = open(datafile, "r")
            datalist = current_datafile.readlines()
            current_datafile.close()

            csa_data = string.split(datalist[15])
            cal_all_data = string.split(datalist[self.cal_ctth_rows[0]])#16])
            cal_low_data = string.split(datalist[self.cal_ctth_rows[1]])#17])
            cal_medium_data = string.split(datalist[self.cal_ctth_rows[2]])#18])
            cal_high_data = string.split(datalist[self.cal_ctth_rows[3]])#19])
            if ' '.join(csa_data) != "No CloudSat":
                # Accumulate CloudSat statistics
                csa_samples += int(csa_data[6])
                mean_error_csa_sum += int(csa_data[6])*float(csa_data[4])
                rms_error_csa_sum += int(csa_data[6])*float(csa_data[5])
            
            # Accumulate CALIOP statistics
            cal_all_samples += int(cal_all_data[7])
            cal_low_samples += int(cal_low_data[7])
            cal_medium_samples += int(cal_medium_data[7])
            cal_high_samples += int(cal_high_data[7])
            mean_error_cal_all_sum += int(cal_all_data[7])*float(cal_all_data[5])
            mean_error_cal_low_sum += int(cal_low_data[7])*float(cal_low_data[5])
            mean_error_cal_medium_sum += int(cal_medium_data[7])*float(cal_medium_data[5])
            mean_error_cal_high_sum += int(cal_high_data[7])*float(cal_high_data[5])
            
            rms_error_cal_all_sum += int(cal_all_data[7])*float(cal_all_data[6])**2
            rms_error_cal_low_sum += int(cal_low_data[7])*float(cal_low_data[6])**2
            rms_error_cal_medium_sum += int(cal_medium_data[7])*float(cal_medium_data[6])**2
            rms_error_cal_high_sum += int(cal_high_data[7])*float(cal_high_data[6])**2
        
        # numpy.divide handles potential division by zero
        if ' '.join(csa_data) != "No CloudSat":
            bias_csa = divide(1.*mean_error_csa_sum, csa_samples)
            rms_csa = divide(1.*rms_error_csa_sum, csa_samples)
        bias_cal_all = divide(1.*mean_error_cal_all_sum, cal_all_samples)
        bias_cal_low = divide(1.*mean_error_cal_low_sum, cal_low_samples)
        bias_cal_medium = divide(1.*mean_error_cal_medium_sum, cal_medium_samples)
        bias_cal_high = divide(1.*mean_error_cal_high_sum, cal_high_samples)
    
        rms_cal_all = math.sqrt(divide(1.*rms_error_cal_all_sum, cal_all_samples))
        rms_cal_low = math.sqrt(divide(1.*rms_error_cal_low_sum, cal_low_samples))
        rms_cal_medium = math.sqrt(divide(1.*rms_error_cal_medium_sum, cal_medium_samples))
        rms_cal_high = math.sqrt(divide(1.*rms_error_cal_high_sum, cal_high_samples))

        bc_rms_cal_all = bias_corrected_rms(rms_cal_all, bias_cal_all, cal_all_samples)
        bc_rms_cal_low = bias_corrected_rms(rms_cal_low, bias_cal_low, cal_low_samples)
        bc_rms_cal_medium = bias_corrected_rms(rms_cal_medium, bias_cal_medium, cal_medium_samples)
        bc_rms_cal_high = bias_corrected_rms(rms_cal_high, bias_cal_high, cal_high_samples)
        self.scenes = scenes
        if ' '.join(csa_data) != "No CloudSat":
            self.csa_samples = csa_samples
            self.bias_csa = bias_csa
            self.rms_csa = rms_csa
        self.cal_all_samples = cal_all_samples
        self.cal_low_samples = cal_low_samples
        self.cal_medium_samples = cal_medium_samples
        self.cal_high_samples = cal_high_samples
        self.bias_cal_all = bias_cal_all
        self.bias_cal_low = bias_cal_low
        self.bias_cal_medium = bias_cal_medium
        self.bias_cal_high = bias_cal_high
        self.rms_cal_all = rms_cal_all
        self.rms_cal_low = rms_cal_low
        self.rms_cal_medium = rms_cal_medium
        self.rms_cal_high = rms_cal_high
        self.bcrms_cal_all = bc_rms_cal_all
        self.bcrms_cal_low = bc_rms_cal_low
        self.bcrms_cal_medium = bc_rms_cal_medium
        self.bcrms_cal_high = bc_rms_cal_high
    
    def printout(self):
        try:
            lines = []
            if self.note is not None:
                lines.append(self.note)
            lines.append("Total number of matched scenes is: %s" % self.scenes)
            lines.append("Total number of Cloudsat matched cloudtops: %d " % self.csa_samples)
            lines.append("Mean error: %f" % self.bias_csa)
            lines.append("RMS error: %f" % self.rms_csa)
            #lines.append("RMS error: %f" % self.rms_csa)
            lines.append("")
            lines.append("Total number of CALIOP matched cloudtops: %d" % self.cal_all_samples)
            lines.append("Number of CALIOP matched low cloudtops: %d" % self.cal_low_samples)
            lines.append("Number of CALIOP matched medium cloudtops: %d" % self.cal_medium_samples)
            lines.append("Number of CALIOP matched high cloudtops: %d" % self.cal_high_samples)
            lines.append("Mean error total cases: %f" % self.bias_cal_all)
            lines.append("Mean error low-level cases: %f" % self.bias_cal_low)
            lines.append("Mean error medium-level cases: %f" % self.bias_cal_medium)
            lines.append("Mean error high-level cases: %f" % self.bias_cal_high)
            lines.append("RMS error total cases: %f" % self.rms_cal_all)
            lines.append("RMS error low-level cases: %f" % self.rms_cal_low)
            lines.append("RMS error medium-level cases: %f" % self.rms_cal_medium)
            lines.append("RMS error high-level cases: %f" % self.rms_cal_high)
            lines.append("bc-RMS error total cases: %f" % self.bcrms_cal_all)
            lines.append("bc-RMS error low-level cases: %f" % self.bcrms_cal_low)
            lines.append("bc-RMS error medium-level cases: %f" % self.bcrms_cal_medium)
            lines.append("bc-RMS error high-level cases: %f" % self.bcrms_cal_high)
            lines.append("")
        except AttributeError:
            lines = []
            if self.note is not None:
                lines.append(self.note)
            lines.append("Total number of matched scenes is: %s" % self.scenes)
            #lines.append("RMS error: %f" % self.rms_csa)
            lines.append("")
            lines.append("Total number of CALIOP matched cloudtops: %d" % self.cal_all_samples)
            lines.append("Number of CALIOP matched low cloudtops: %d" % self.cal_low_samples)
            lines.append("Number of CALIOP matched medium cloudtops: %d" % self.cal_medium_samples)
            lines.append("Number of CALIOP matched high cloudtops: %d" % self.cal_high_samples)
            lines.append("Mean error total cases: %f" % self.bias_cal_all)
            lines.append("Mean error low-level cases: %f" % self.bias_cal_low)
            lines.append("Mean error medium-level cases: %f" % self.bias_cal_medium)
            lines.append("Mean error high-level cases: %f" % self.bias_cal_high)
            lines.append("RMS error total cases: %f" % self.rms_cal_all)
            lines.append("RMS error low-level cases: %f" % self.rms_cal_low)
            lines.append("RMS error medium-level cases: %f" % self.rms_cal_medium)
            lines.append("RMS error high-level cases: %f" % self.rms_cal_high)
            lines.append("bc-RMS error total cases: %f" % self.bcrms_cal_all)
            lines.append("bc-RMS error low-level cases: %f" % self.bcrms_cal_low)
            lines.append("bc-RMS error medium-level cases: %f" % self.bcrms_cal_medium)
            lines.append("bc-RMS error high-level cases: %f" % self.bcrms_cal_high)
            lines.append("")
        return lines


if __name__ == "__main__":
    stats = CloudTopStats()

