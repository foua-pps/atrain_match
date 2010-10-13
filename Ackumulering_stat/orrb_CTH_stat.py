# Program orrb_CTH_stat

# This program calculates basic statistics for the cloud top (CTH) product for
# each month

import string, glob
import math

from orrb_conf import SATELLITE, RESOLUTION, STUDIED_YEAR, STUDIED_MONTHS, MAP, MAIN_DATADIR, OUTPUT_DIR

# -----------------------------------------------------
if __name__ == "__main__":
    result=[]
    for i in range(len(STUDIED_MONTHS)):
        month=("%s%s" %(STUDIED_YEAR[0],STUDIED_MONTHS[i]))
        basic_indata_dir = "%s/Results/%s/%s/%s/%s/%s/BASIC/" % (MAIN_DATADIR, SATELLITE[0], RESOLUTION[0] , STUDIED_YEAR[0], STUDIED_MONTHS[i], MAP[0])
        #datafiles = glob.glob("%s*%s*1.dat" % (BASIC_INDATA_DIR,month))
        #datafiles2 = glob.glob("%s*%s*2.dat" % (BASIC_INDATA_DIR,month))
        datafiles = glob.glob("%s*.dat" % basic_indata_dir)
        scenes = len(datafiles)
        
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
        for datafile in datafiles:
            current_datafile = open(datafile, "r")
            datalist = current_datafile.readlines()
            csa_data = string.split(datalist[15])
            cal_all_data = string.split(datalist[16])
            cal_low_data = string.split(datalist[17])
            cal_medium_data = string.split(datalist[18])
            cal_high_data = string.split(datalist[19])

            # Accumulate CloudSat statistics
                                
            csa_samples = csa_samples + int(csa_data[6])
            mean_error_csa_sum = mean_error_csa_sum + int(csa_data[6])*float(csa_data[4])  
            rms_error_csa_sum = rms_error_csa_sum + int(csa_data[6])*float(csa_data[5])  

            

            # Accumulate CALIOP statistics

            cal_all_samples = cal_all_samples + int(cal_all_data[7])
            cal_low_samples = cal_low_samples + int(cal_low_data[7])
            cal_medium_samples = cal_medium_samples + int(cal_medium_data[7])
            cal_high_samples = cal_high_samples + int(cal_high_data[7])
            mean_error_cal_all_sum = mean_error_cal_all_sum + \
                                        int(cal_all_data[7])*float(cal_all_data[5])  
            mean_error_cal_low_sum = mean_error_cal_low_sum + \
                                        int(cal_low_data[7])*float(cal_low_data[5])  
            mean_error_cal_medium_sum = mean_error_cal_medium_sum + \
                                        int(cal_medium_data[7])*float(cal_medium_data[5])  
            mean_error_cal_high_sum = mean_error_cal_high_sum + \
                                        int(cal_high_data[7])*float(cal_high_data[5])  

# Notice that the original linear averaging of RMS values is wrong! We have now changed to a correct averaging!/KG 20091126
##             rms_error_cal_all_sum = rms_error_cal_all_sum + \
##                                      int(cal_all_data[7])*float(cal_all_data[6])  
##             rms_error_cal_low_sum = rms_error_cal_low_sum + \
##                                      int(cal_low_data[7])*float(cal_low_data[6])  
##             rms_error_cal_medium_sum = rms_error_cal_medium_sum + \
##                                         int(cal_medium_data[7])*float(cal_medium_data[6])  
##             rms_error_cal_high_sum = rms_error_cal_high_sum + \
##                                       int(cal_high_data[7])*float(cal_high_data[6])  
            rms_error_cal_all_sum = rms_error_cal_all_sum + \
                                    int(cal_all_data[7])*float(cal_all_data[6])*float(cal_all_data[6])  
            rms_error_cal_low_sum = rms_error_cal_low_sum + \
                                    int(cal_low_data[7])*float(cal_low_data[6])*float(cal_low_data[6])  
            rms_error_cal_medium_sum = rms_error_cal_medium_sum + \
                                        int(cal_medium_data[7])*float(cal_medium_data[6])*float(cal_medium_data[6])  
            rms_error_cal_high_sum = rms_error_cal_high_sum + \
                                        int(cal_high_data[7])*float(cal_high_data[6])*float(cal_high_data[6])  
            
            current_datafile.close()
            
        bias_csa = mean_error_csa_sum/csa_samples
        bias_cal_all = mean_error_cal_all_sum/cal_all_samples
        bias_cal_low = mean_error_cal_low_sum/cal_low_samples
        bias_cal_medium = mean_error_cal_medium_sum/cal_medium_samples
        bias_cal_high = mean_error_cal_high_sum/cal_high_samples
        rms_csa = rms_error_csa_sum/csa_samples


# Notice that the original linear averaging of RMS values is wrong! We have now changed to a correct averaging!/KG 20091126
##         rms_cal_all = rms_error_cal_all_sum/cal_all_samples
##         rms_cal_low = rms_error_cal_low_sum/cal_low_samples
##         rms_cal_medium = rms_error_cal_medium_sum/cal_medium_samples
##         rms_cal_high = rms_error_cal_high_sum/cal_high_samples

        rms_cal_all = math.sqrt(rms_error_cal_all_sum/cal_all_samples)
        rms_cal_low = math.sqrt(rms_error_cal_low_sum/cal_low_samples)
        rms_cal_medium = math.sqrt(rms_error_cal_medium_sum/cal_medium_samples)
        rms_cal_high = math.sqrt(rms_error_cal_high_sum/cal_high_samples)

##         square_sum_csa =  float(n_clear_clear_csa+n_cloudy_cloudy_csa)*bias_csa*bias_csa + \
##                          n_cloudy_clear_csa*(-1.0-bias_csa)*(-1.0-bias_csa) + \
##                          n_clear_cloudy_csa*(1.0-bias_csa)*(1.0-bias_csa)
##         rms_csa = 100.0*math.sqrt(square_sum_csa/(samples_csa-1))
##         square_sum_cal =  float(n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_cal*bias_cal + \
##                          n_cloudy_clear_cal*(-1.0-bias_cal)*(-1.0-bias_cal) + \
##                          n_clear_cloudy_cal*(1.0-bias_cal)*(1.0-bias_cal)
##         rms_cal = 100.0*math.sqrt(square_sum_cal/(samples_cal-1))
        
        print "Month is:  %s" % month
        print "Total number of matched scenes is: %s" % scenes
        print "Total number of Cloudsat matched cloudtops: %d " % csa_samples
        print "Mean error: %f" % bias_csa
        print "Weighted RMS error: %f" % rms_csa
        #print "RMS error: %f" % rms_csa
        print
        print "Total number of CALIOP matched cloudtops: %d" % cal_all_samples
        print "Number of CALIOP matched low cloudtops: %d" % cal_low_samples
        print "Number of CALIOP matched medium cloudtops: %d" % cal_medium_samples
        print "Number of CALIOP matched high cloudtops: %d" % cal_high_samples
        print "Mean error total cases: %f" % bias_cal_all
        print "Mean error low-level cases: %f" % bias_cal_low
        print "Mean error medium-level cases: %f" % bias_cal_medium
        print "Mean error high-level cases: %f" % bias_cal_high
        print "Weighted RMS error total cases: %f" % rms_cal_all
        print "Weighted RMS error low-level cases: %f" % rms_cal_low
        print "Weighted RMS error medium-level cases: %f" % rms_cal_medium
        print "Weighted RMS error high-level cases: %f" % rms_cal_high
        #print "RMS error: %f" % rms_cal
        print
        
        result.append("Month is:  %s\n" % month)
        result.append("Total number of matched scenes is: %s\n" % scenes)
        result.append("Total number of Cloudsat matched cloudtops: %d \n" % csa_samples)
        result.append("Mean error: %f\n" % bias_csa)
        result.append("Weighted RMS error: %f\n" % rms_csa)
        #result.append("RMS error: %f\n" % rms_csa)
        result.append("\n")
        result.append("Total number of CALIOP matched cloudtops: %d\n" % cal_all_samples)
        result.append("Number of CALIOP matched low cloudtops: %d\n" % cal_low_samples)
        result.append("Number of CALIOP matched medium cloudtops: %d\n" % cal_medium_samples)
        result.append("Number of CALIOP matched high cloudtops: %d\n" % cal_high_samples)
        result.append("Mean error total cases: %f\n" % bias_cal_all)
        result.append("Mean error low-level cases: %f\n" % bias_cal_low)
        result.append("Mean error medium-level cases: %f\n" % bias_cal_medium)
        result.append("Mean error high-level cases: %f\n" % bias_cal_high)
        result.append("Weighted RMS error total cases: %f\n" % rms_cal_all)
        result.append("Weighted RMS error low-level cases: %f\n" % rms_cal_low)
        result.append("Weighted RMS error medium-level cases: %f\n" % rms_cal_medium)
        result.append("Weighted RMS error high-level cases: %f\n" % rms_cal_high)
        result.append("\n")
    fd=open("%s/cth_results_summary_%s%s-%s%s.dat" %(OUTPUT_DIR, STUDIED_YEAR[0],STUDIED_MONTHS[0],STUDIED_YEAR[-1],STUDIED_MONTHS[-1]),'w')
    fd.writelines(result)
    fd.close()
