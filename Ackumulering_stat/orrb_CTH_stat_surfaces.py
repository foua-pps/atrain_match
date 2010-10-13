# Program orrb_CTH_stat

# This program calculates basic statistics for the cloud top (CTH) product for
# each month

import string, glob

from orrb_conf import SATELLITE, RESOLUTION, STUDIED_YEAR, STUDIED_MONTHS, MAP, MAIN_DATADIR, OUTPUT_DIR, SURFACES

# -----------------------------------------------------
if __name__ == "__main__":
    result=[]
    for surface in SURFACES:
        print ""
        print ""
        print "Surface is: ", surface
        print
        #result.append("%s\n" %("-"*50))
        result.append("Surface is: %s\n" % surface)
        result.append("\n")
        no_value=0
        for i in range(len(STUDIED_MONTHS)):
            month="%s%s" %(STUDIED_YEAR[0],STUDIED_MONTHS[i])
            surface_indata_dir = "%s/Results/%s/%s/%s/%s/%s/%s/" % (MAIN_DATADIR, SATELLITE[0], RESOLUTION[0] , STUDIED_YEAR[0], STUDIED_MONTHS[i], MAP[0],surface)
            datafiles = glob.glob("%s*.dat" % surface_indata_dir)
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
                items = len(datalist)
                csa_data = string.split(datalist[15])
                if items > 16:
                    cal_all_data = string.split(datalist[16])
                if items > 17:
                    cal_low_data = string.split(datalist[17])
                if items > 18:
                    cal_medium_data = string.split(datalist[18])
                if items > 19:
                    cal_high_data = string.split(datalist[19])

                # Accumulate CALIOP statistics

                if items > 16:
                    cal_all_samples = cal_all_samples + int(cal_all_data[7])
                    mean_error_cal_all_sum = mean_error_cal_all_sum + \
                                             int(cal_all_data[7])*float(cal_all_data[5])  
                    rms_error_cal_all_sum = rms_error_cal_all_sum + \
                                            int(cal_all_data[7])*float(cal_all_data[6])  
                if items > 17:
                    cal_low_samples = cal_low_samples + int(cal_low_data[7])
                    mean_error_cal_low_sum = mean_error_cal_low_sum + \
                                             int(cal_low_data[7])*float(cal_low_data[5])  
                    rms_error_cal_low_sum = rms_error_cal_low_sum + \
                                            int(cal_low_data[7])*float(cal_low_data[6])  
                if items > 18:
                    cal_medium_samples = cal_medium_samples + int(cal_medium_data[7])
                    mean_error_cal_medium_sum = mean_error_cal_medium_sum + \
                                                int(cal_medium_data[7])*float(cal_medium_data[5])  
                    rms_error_cal_medium_sum = rms_error_cal_medium_sum + \
                                               int(cal_medium_data[7])*float(cal_medium_data[6])  
                if items > 19:
                    cal_high_samples = cal_high_samples + int(cal_high_data[7])
                    mean_error_cal_high_sum = mean_error_cal_high_sum + \
                                              int(cal_high_data[7])*float(cal_high_data[5])  
                    rms_error_cal_high_sum = rms_error_cal_high_sum + \
                                             int(cal_high_data[7])*float(cal_high_data[6])  
            
                current_datafile.close()
            
            if cal_all_samples > 0:
                bias_cal_all = mean_error_cal_all_sum/cal_all_samples
                rms_cal_all = rms_error_cal_all_sum/cal_all_samples
            else:
                bias_cal_all = -9.0
                rms_cal_all = -9.0
            if cal_low_samples > 0: 
                bias_cal_low = mean_error_cal_low_sum/cal_low_samples
                rms_cal_low = rms_error_cal_low_sum/cal_low_samples
            else:
                bias_cal_low = -9.0
                rms_cal_low = -9.0
            if cal_medium_samples > 0:  
                bias_cal_medium = mean_error_cal_medium_sum/cal_medium_samples
                rms_cal_medium = rms_error_cal_medium_sum/cal_medium_samples
            else: 
                bias_cal_medium = -9.0
                rms_cal_medium = -9.0
            if cal_high_samples > 0:                
                bias_cal_high = mean_error_cal_high_sum/cal_high_samples
                rms_cal_high = rms_error_cal_high_sum/cal_high_samples
            else:
                bias_cal_high = -9.0
                rms_cal_high = -9.0
            

            print "Month is:  %s" % month
            print "Total number of matched scenes is: %s" % scenes
            #print "Total number of Cloudsat matched cloudtops: %d " % csa_samples
            #print "Mean error: %f" % bias_csa
            #print "Weighted RMS error: %f" % rms_csa
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
        result.append("\n")
        result.append("\n")
    result.append("\n")            
    fd=open("%s/cth_results_surfaces_summary_%s%s-%s%s.dat" %(OUTPUT_DIR, STUDIED_YEAR[0],STUDIED_MONTHS[0],STUDIED_YEAR[-1],STUDIED_MONTHS[-1]),'w')
    fd.writelines(result)
    fd.close()        
            
            
            
            
            
            
