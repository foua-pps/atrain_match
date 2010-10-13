# Program orrb_CFC_stat

# This program calculates basic statistics for the cloud amount (CFC) product for
# each month

import string, glob
import math

from orrb_conf import SATELLITE, RESOLUTION, STUDIED_YEAR, STUDIED_MONTHS, MAP, MAIN_DATADIR, OUTPUT_DIR

# -----------------------------------------------------
if __name__ == "__main__":
    result=[]
    for i in range(len(STUDIED_MONTHS)):
        month=(STUDIED_YEAR[0],STUDIED_MONTHS[i])
        basic_indata_dir = "%s/Results/%s/%s/%s/%s/%s/BASIC/" % (MAIN_DATADIR, SATELLITE[0], RESOLUTION[0] , STUDIED_YEAR[0], STUDIED_MONTHS[i], MAP[0])
        
        datafiles = glob.glob("%s*.dat" % basic_indata_dir)
        #pdb.set_trace()
        #datafiles = glob.glob("%s*%s*1.dat" % (BASIC_INDATA_DIR,month))
        #datafiles2 = glob.glob("%s*%s*2.dat" % (BASIC_INDATA_DIR,month))
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
        samples_csa = 0
        samples_cal = 0

        scenes = len(datafiles)
        cfc_sum_csa = 0
        cfc_sum_cal = 0
        
        for datafile in datafiles:

            current_datafile = open(datafile, "r")
            datalist = current_datafile.readlines()
            current_datafile.close()
            #print "Datafile: ", datafile
            csa_data = string.split(datalist[4])
            cal_data = string.split(datalist[8])
            modis_data = string.split(datalist[10])

            # Accumulate CloudSat statistics
            
            n_clear_clear_csa = n_clear_clear_csa + int(csa_data[4])
            n_clear_cloudy_csa = n_clear_cloudy_csa + int(csa_data[5])
            n_cloudy_clear_csa = n_cloudy_clear_csa + int(csa_data[6])
            n_cloudy_cloudy_csa = n_cloudy_cloudy_csa + int(csa_data[7])
            
            # Accumulate CALIOP statistics
            
            n_clear_clear_cal = n_clear_clear_cal + int(cal_data[4])
            n_clear_cloudy_cal = n_clear_cloudy_cal + int(cal_data[5])
            n_cloudy_clear_cal = n_cloudy_clear_cal + int(cal_data[6])
            n_cloudy_cloudy_cal = n_cloudy_cloudy_cal + int(cal_data[7])

            # Accumulate CALIOP-MODIS statistics
            
            n_clear_clear_cal_MODIS = n_clear_clear_cal_MODIS + int(modis_data[4])
            n_clear_cloudy_cal_MODIS = n_clear_cloudy_cal_MODIS + int(modis_data[5])
            n_cloudy_clear_cal_MODIS = n_cloudy_clear_cal_MODIS + int(modis_data[6])
            n_cloudy_cloudy_cal_MODIS = n_cloudy_cloudy_cal_MODIS + int(modis_data[7])

            current_datafile.close()
            
        samples_csa = n_clear_clear_csa + n_clear_cloudy_csa + n_cloudy_clear_csa +\
                        n_cloudy_cloudy_csa
        samples_cal = n_clear_clear_cal + n_clear_cloudy_cal + n_cloudy_clear_cal +\
                        n_cloudy_cloudy_cal
        mean_CFC_csa = 100.0*(n_cloudy_cloudy_csa+n_cloudy_clear_csa)/samples_csa
        mean_CFC_cal = 100.0*(n_cloudy_cloudy_cal+n_cloudy_clear_cal)/samples_cal
        bias_csa = float(n_clear_cloudy_csa - n_cloudy_clear_csa)/float(samples_csa)
        bias_cal = float(n_clear_cloudy_cal - n_cloudy_clear_cal)/float(samples_cal)
        bias_modis = float(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS)/float(samples_cal-1)
        bias_csa_perc = 100.0*float(n_clear_cloudy_csa - n_cloudy_clear_csa)/float(samples_csa)
        bias_cal_perc = 100.0*float(n_clear_cloudy_cal - n_cloudy_clear_cal)/float(samples_cal)
        bias_modis_perc = 100.0*float(n_clear_cloudy_cal_MODIS - n_cloudy_clear_cal_MODIS)/float(samples_cal-1)

        square_sum_csa =  float(n_clear_clear_csa+n_cloudy_cloudy_csa)*bias_csa*bias_csa + \
                            n_cloudy_clear_csa*(-1.0-bias_csa)*(-1.0-bias_csa) + \
                            n_clear_cloudy_csa*(1.0-bias_csa)*(1.0-bias_csa)
        rms_csa = 100.0*math.sqrt(square_sum_csa/(samples_csa-1))
        square_sum_cal =  float(n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_cal*bias_cal + \
                            n_cloudy_clear_cal*(-1.0-bias_cal)*(-1.0-bias_cal) + \
                            n_clear_cloudy_cal*(1.0-bias_cal)*(1.0-bias_cal)
        rms_cal = 100.0*math.sqrt(square_sum_cal/(samples_cal-1))
        square_sum_modis =  float(n_clear_clear_cal+n_cloudy_cloudy_cal)*bias_modis*bias_modis + \
                            n_cloudy_clear_cal*(-1.0-bias_modis)*(-1.0-bias_modis) + \
                            n_clear_cloudy_cal*(1.0-bias_modis)*(1.0-bias_modis)
        rms_modis = 100.0*math.sqrt(square_sum_modis/(samples_cal-1))

        pod_cloudy_cal = 100.0*n_cloudy_cloudy_cal/(n_cloudy_cloudy_cal+n_cloudy_clear_cal)
        pod_cloudy_cal_MODIS = 100.0*n_cloudy_cloudy_cal_MODIS/(n_cloudy_cloudy_cal_MODIS+n_cloudy_clear_cal_MODIS)
        pod_clear_cal = 100.0*n_clear_clear_cal/(n_clear_clear_cal+n_clear_cloudy_cal)
        pod_clear_cal_MODIS = 100.0*n_clear_clear_cal_MODIS/(n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS)
##         far_cloudy_cal = 100.0*n_cloudy_clear_cal/(n_cloudy_cloudy_cal+n_cloudy_clear_cal) #Not correct!
##         far_clear_cal = 100.0*n_clear_cloudy_cal/(n_clear_clear_cal+n_clear_cloudy_cal)    #Not correct!
        far_cloudy_cal = 100.0*n_clear_cloudy_cal/(n_cloudy_cloudy_cal+n_clear_cloudy_cal)
        far_clear_cal = 100.0*n_cloudy_clear_cal/(n_clear_clear_cal+n_cloudy_clear_cal)
        far_cloudy_cal_MODIS = 100.0*n_clear_cloudy_cal_MODIS/(n_cloudy_cloudy_cal_MODIS+n_clear_cloudy_cal_MODIS)
        far_clear_cal_MODIS = 100.0*n_cloudy_clear_cal_MODIS/(n_clear_clear_cal_MODIS+n_cloudy_clear_cal_MODIS)

        kuipers = 1.0*(n_clear_clear_cal*n_cloudy_cloudy_cal-n_cloudy_clear_cal*n_clear_cloudy_cal)/\
                    ((n_clear_clear_cal+n_clear_cloudy_cal)*(n_cloudy_clear_cal+n_cloudy_cloudy_cal))

        hitrate = 1.0*(n_clear_clear_cal+n_cloudy_cloudy_cal)/(n_clear_clear_cal+n_clear_cloudy_cal+\
                                                            n_cloudy_clear_cal+n_cloudy_cloudy_cal)

        kuipers_MODIS = 1.0*(n_clear_clear_cal_MODIS*n_cloudy_cloudy_cal_MODIS-n_cloudy_clear_cal_MODIS*n_clear_cloudy_cal_MODIS)/((n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS)*(n_cloudy_clear_cal_MODIS+n_cloudy_cloudy_cal_MODIS))

        hitrate_MODIS = 1.0*(n_clear_clear_cal_MODIS+n_cloudy_cloudy_cal_MODIS)/\
                        (n_clear_clear_cal_MODIS+n_clear_cloudy_cal_MODIS+n_cloudy_clear_cal_MODIS+\
                            n_cloudy_cloudy_cal_MODIS)
        
        print "Month is:  %s%s" % month
        print "Total number of matched scenes is: %s" % scenes
        print "Total number of Cloudsat matched FOVs: %d " % samples_csa
        print "Mean CFC Cloudsat: %f " % mean_CFC_csa
        print "Mean error: %f" % bias_csa_perc
        print "RMS error: %f" % rms_csa
        print
        print "Total number of CALIOP matched FOVs: %d" % samples_cal
        print "Mean CFC CALIOP: %f " % mean_CFC_cal
        print "Mean error: %f" % bias_cal_perc
        print "RMS error: %f" % rms_cal
        print "Mean error MODIS: %f" % bias_modis_perc
        print "RMS error MODIS: %f" % rms_modis
        print "POD cloudy: %f" % pod_cloudy_cal
        print "POD cloudy MODIS: %f" % pod_cloudy_cal_MODIS
        print "POD clear: %f" % pod_clear_cal
        print "POD clear MODIS: %f" % pod_clear_cal_MODIS
        print "FAR cloudy: %f" % far_cloudy_cal
        print "FAR cloudy MODIS: %f" % far_cloudy_cal_MODIS
        print "FAR clear: %f" % far_clear_cal
        print "FAR clear MODIS: %f" % far_clear_cal_MODIS
        print "Kuipers: %f" % kuipers
        print "Kuipers MODIS: %f" % kuipers_MODIS
        print "Hitrate: %f" % hitrate
        print "Hitrate MODIS: %f" % hitrate_MODIS
        print
        result.append("Month is:  %s%s\n" % month)
        result.append("Total number of matched scenes is: %s\n" % scenes)
        result.append("Total number of Cloudsat matched FOVs: %d \n" % samples_csa)
        result.append("Mean CFC Cloudsat: %f \n" % mean_CFC_csa)
        result.append("Mean error: %f\n" % bias_csa_perc)
        result.append("RMS error: %f\n" % rms_csa)
        result.append("\n")
        result.append("Total number of CALIOP matched FOVs: %d\n" % samples_cal)
        result.append("Mean CFC CALIOP: %f \n" % mean_CFC_cal)
        result.append("Mean error: %f\n" % bias_cal_perc)
        result.append("RMS error: %f\n" % rms_cal)
        result.append("Mean error MODIS: %f\n" % bias_modis_perc)
        result.append("RMS error MODIS: %f\n" % rms_modis)
        result.append("POD cloudy: %f\n" % pod_cloudy_cal)
        result.append("POD cloudy MODIS: %f\n" % pod_cloudy_cal_MODIS)
        result.append("POD clear: %f\n" % pod_clear_cal)
        result.append("POD clear MODIS: %f\n" % pod_clear_cal_MODIS)
        result.append("FAR cloudy: %f\n" % far_cloudy_cal)
        result.append("FAR cloudy MODIS: %f\n" % far_cloudy_cal_MODIS)
        result.append("FAR clear: %f\n" % far_clear_cal)
        result.append("FAR clear MODIS: %f\n" % far_clear_cal_MODIS)
        result.append("Kuipers: %f\n" % kuipers)
        result.append("Kuipers MODIS: %f\n" % kuipers_MODIS)
        result.append("Hitrate: %f\n" % hitrate)
        result.append("Hitrate MODIS: %f\n" % hitrate_MODIS)
        result.append("\n")
    fd=open("%s/cfc_results_summary_%s%s-%s%s.dat" %(OUTPUT_DIR, STUDIED_YEAR[0],STUDIED_MONTHS[0],STUDIED_YEAR[-1],STUDIED_MONTHS[-1]),'w')
    fd.writelines(result)
    fd.close()

#pdb.set_trace()
