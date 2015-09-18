# Program orrb_CFC_stat

# This program calculates basic statistics for the cloud amount (CFC) product for
# each month

import string
import math

from orrb_stat_class import OrrbStats

# -----------------------------------------------------
class CloudFractionStats(OrrbStats):
    def do_stats(self):
        OrrbStats.do_stats(self)
        
        from numpy import divide
    
        scenes = len(self.results_files)
        
        n_clear_clear_csa = 0
        n_clear_cloudy_csa = 0
        n_cloudy_clear_csa = 0
        n_cloudy_cloudy_csa = 0
        n_clear_clear_cal = 0
        n_clear_cloudy_cal = 0
        n_cloudy_clear_cal = 0
        n_cloudy_cloudy_cal = 0
        n_clear_clear_cal_MODIS = 0
        n_clear_cloudy_cal_MODIS = 0
        n_cloudy_clear_cal_MODIS = 0
        n_cloudy_cloudy_cal_MODIS = 0
        for datafile in self.results_files:
            current_datafile = open(datafile, "r")
            datalist = current_datafile.readlines()
            current_datafile.close()
            #print "Datafile: ", datafile
            csa_data = string.split(datalist[4])
            cal_data = string.split(datalist[8])
            modis_data = string.split(datalist[10])

            # Accumulate CloudSat statistics
            if ' '.join(csa_data) != "No CloudSat":
                n_clear_clear_csa += int(csa_data[4])
                n_clear_cloudy_csa += int(csa_data[5])
                n_cloudy_clear_csa += int(csa_data[6])
                n_cloudy_cloudy_csa += int(csa_data[7])
            
            # Accumulate CALIOP statistics
            
            n_clear_clear_cal += int(cal_data[4])
            n_clear_cloudy_cal += int(cal_data[5])
            n_cloudy_clear_cal += int(cal_data[6])
            n_cloudy_cloudy_cal += int(cal_data[7])
    
            # Accumulate CALIOP-MODIS statistics
            if ' '.join(csa_data) != "No CloudSat":
                n_clear_clear_cal_MODIS += int(modis_data[4])
                n_clear_cloudy_cal_MODIS += int(modis_data[5])
                n_cloudy_clear_cal_MODIS += int(modis_data[6])
                n_cloudy_cloudy_cal_MODIS += int(modis_data[7])
        if ' '.join(csa_data) != "No CloudSat":
            samples_csa = n_clear_clear_csa + n_clear_cloudy_csa + n_cloudy_clear_csa +\
                            n_cloudy_cloudy_csa
        samples_cal = n_clear_clear_cal + n_clear_cloudy_cal + n_cloudy_clear_cal +\
                        n_cloudy_cloudy_cal
        
        # numpy.divide handles potential division by zero
        if ' '.join(csa_data) != "No CloudSat":
            bias_csa = divide(1.*(n_clear_cloudy_csa - n_cloudy_clear_csa), samples_csa)
            bias_csa_perc = divide(100.0*(n_clear_cloudy_csa - n_cloudy_clear_csa), samples_csa)
            mean_CFC_csa = divide(100.0*(n_cloudy_cloudy_csa+n_cloudy_clear_csa), samples_csa)
            bias_modis = divide(1.*(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS), samples_cal-1.)
            bias_modis_perc = divide(100.0*(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS), samples_cal-1.)
        mean_CFC_cal = divide(100.0*(n_cloudy_cloudy_cal+n_cloudy_clear_cal), samples_cal)
        bias_cal = divide(1.*(n_clear_cloudy_cal - n_cloudy_clear_cal), samples_cal)        
        bias_cal_perc = divide(100.0*(n_clear_cloudy_cal - n_cloudy_clear_cal), samples_cal)
        
        if ' '.join(csa_data) != "No CloudSat":
            square_sum_csa =  (n_clear_clear_csa+n_cloudy_cloudy_csa)*bias_csa**2 + \
                                n_cloudy_clear_csa*(-1.0-bias_csa)**2 + \
                                n_clear_cloudy_csa*(1.0-bias_csa)**2
            rms_csa = 100.0*math.sqrt(divide(square_sum_csa, samples_csa-1.))
        
        square_sum_cal =  (n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_cal**2 + \
                            n_cloudy_clear_cal*(-1.0-bias_cal)**2 + \
                            n_clear_cloudy_cal*(1.0-bias_cal)**2
        rms_cal = 100.0*math.sqrt(divide(square_sum_cal, samples_cal-1.))
        
        if ' '.join(csa_data) != "No CloudSat":
            square_sum_modis =  (n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_modis**2 + \
                                n_cloudy_clear_cal*(-1.0-bias_modis)**2 + \
                                n_clear_cloudy_cal*(1.0-bias_modis)**2
            rms_modis = 100.0*math.sqrt(divide(square_sum_modis, samples_cal-1.))
        if ' '.join(csa_data) != "No CloudSat":
            pod_cloudy_cal_MODIS = divide(100.0*n_cloudy_cloudy_cal_MODIS, n_cloudy_cloudy_cal_MODIS+n_cloudy_clear_cal_MODIS)
            pod_clear_cal_MODIS = divide(100.0*n_clear_clear_cal_MODIS, n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS)
            far_cloudy_cal_MODIS = divide(100.0*n_clear_cloudy_cal_MODIS, n_cloudy_cloudy_cal_MODIS+n_clear_cloudy_cal_MODIS)
            far_clear_cal_MODIS = divide(100.0*n_cloudy_clear_cal_MODIS, n_clear_clear_cal_MODIS+n_cloudy_clear_cal_MODIS)
        pod_cloudy_cal = divide(100.0*n_cloudy_cloudy_cal, (n_cloudy_cloudy_cal+n_cloudy_clear_cal))
        pod_clear_cal = divide(100.0*n_clear_clear_cal, n_clear_clear_cal+n_clear_cloudy_cal)
        far_cloudy_cal = divide(100.0*n_clear_cloudy_cal, n_cloudy_cloudy_cal+n_clear_cloudy_cal)
        far_clear_cal = divide(100.0*n_cloudy_clear_cal, n_clear_clear_cal+n_cloudy_clear_cal)
    
        kuipers = divide(1.0*(n_clear_clear_cal*n_cloudy_cloudy_cal-n_cloudy_clear_cal*n_clear_cloudy_cal),
                         ((n_clear_clear_cal+n_clear_cloudy_cal)*(n_cloudy_clear_cal+n_cloudy_cloudy_cal)))
    
        hitrate = divide(1.0*(n_clear_clear_cal+n_cloudy_cloudy_cal),
                         (n_clear_clear_cal+n_clear_cloudy_cal+n_cloudy_clear_cal+n_cloudy_cloudy_cal))
        if ' '.join(csa_data) != "No CloudSat":
            kuipers_MODIS = divide(1.0*(n_clear_clear_cal_MODIS*n_cloudy_cloudy_cal_MODIS-n_cloudy_clear_cal_MODIS*n_clear_cloudy_cal_MODIS),
                                   ((n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS)*(n_cloudy_clear_cal_MODIS+n_cloudy_cloudy_cal_MODIS)))
        
            hitrate_MODIS = divide(1.0*(n_clear_clear_cal_MODIS+n_cloudy_cloudy_cal_MODIS),
                                   (n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS+n_cloudy_clear_cal_MODIS+ n_cloudy_cloudy_cal_MODIS))
        
        # Store values of interest as attributes
        self.scenes = scenes
        if ' '.join(csa_data) != "No CloudSat":
            self.samples_csa = samples_csa
            self.mean_CFC_csa = mean_CFC_csa
            self.bias_csa_perc = bias_csa_perc
            self.rms_csa = rms_csa
            self.bias_modis_perc = bias_modis_perc
            self.rms_modis = rms_modis
            self.pod_cloudy_cal_MODIS = pod_cloudy_cal_MODIS
            self.pod_clear_cal_MODIS = pod_clear_cal_MODIS
            self.far_cloudy_cal_MODIS = far_cloudy_cal_MODIS
            self.far_clear_cal_MODIS = far_clear_cal_MODIS
            self.kuipers_MODIS = kuipers_MODIS
            self.hitrate_MODIS = hitrate_MODIS
        self.samples_cal = samples_cal
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
        try:
            lines = []
            lines.append("Total number of matched scenes is: %s" % self.scenes)
            lines.append("Total number of Cloudsat matched FOVs: %d " % self.samples_csa)
            lines.append("Mean CFC Cloudsat: %f " % self.mean_CFC_csa)
            lines.append("Mean error: %f" % self.bias_csa_perc)
            lines.append("RMS error: %f" % self.rms_csa)
            lines.append("")
            lines.append("Total number of CALIOP matched FOVs: %d" % self.samples_cal)
            lines.append("Mean CFC CALIOP: %f " % self.mean_CFC_cal)
            lines.append("Mean error: %f" % self.bias_cal_perc)
            lines.append("RMS error: %f" % self.rms_cal)
            lines.append("Mean error MODIS: %f" % self.bias_modis_perc)
            lines.append("RMS error MODIS: %f" % self.rms_modis)
            lines.append("POD cloudy: %f" % self.pod_cloudy_cal)
            lines.append("POD cloudy MODIS: %f" % self.pod_cloudy_cal_MODIS)
            lines.append("POD clear: %f" % self.pod_clear_cal)
            lines.append("POD clear MODIS: %f" % self.pod_clear_cal_MODIS)
            lines.append("FAR cloudy: %f" % self.far_cloudy_cal)
            lines.append("FAR cloudy MODIS: %f" % self.far_cloudy_cal_MODIS)
            lines.append("FAR clear: %f" % self.far_clear_cal)
            lines.append("FAR clear MODIS: %f" % self.far_clear_cal_MODIS)
            lines.append("Kuipers: %f" % self.kuipers)
            lines.append("Kuipers MODIS: %f" % self.kuipers_MODIS)
            lines.append("Hitrate: %f" % self.hitrate)
            lines.append("Hitrate MODIS: %f" % self.hitrate_MODIS)
            lines.append("")
        except AttributeError:
            lines = []
            lines.append("Total number of matched scenes is: %s" % self.scenes)
            lines.append("")
            lines.append("Total number of CALIOP matched FOVs: %d" % self.samples_cal)
            lines.append("Mean CFC CALIOP: %f " % self.mean_CFC_cal)
            lines.append("Mean error: %f" % self.bias_cal_perc)
            lines.append("RMS error: %f" % self.rms_cal)
            lines.append("POD cloudy: %f" % self.pod_cloudy_cal)
            lines.append("POD clear: %f" % self.pod_clear_cal)
            lines.append("FAR cloudy: %f" % self.far_cloudy_cal)
            lines.append("FAR clear: %f" % self.far_clear_cal)
            lines.append("Kuipers: %f" % self.kuipers)
            lines.append("Hitrate: %f" % self.hitrate)
            lines.append("")
        return lines


if __name__ == "__main__":
    stats = CloudFractionStats()

