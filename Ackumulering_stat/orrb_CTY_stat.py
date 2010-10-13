# Program orrb_CTY_stat

# This program calculates basic statistics for the cloud type (CTY) product for
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
        datafiles = glob.glob("%s*.dat" % basic_indata_dir)
        if month == "200706":
            CFC_CALIOP = 67.54
            CFC_PPS = 66.24
            Totpix = 92116
        elif month == "200707":
            CFC_CALIOP = 75.26
            CFC_PPS = 73.63
            Totpix = 115606
        elif month == "200708":
            CFC_CALIOP = 79.59
            CFC_PPS = 72.59
            Totpix = 116332
        elif month == "200712":
            CFC_CALIOP = 62.44
            CFC_PPS = 32.11
            Totpix = 73785
        elif month == "200808":
            CFC_CALIOP = 81.45
            CFC_MEAN_ERROR = 0.024868
            CFC_PPS = CFC_CALIOP + CFC_MEAN_ERROR
            Totpix = 40212
        elif month == "200809": 
            CFC_CALIOP = 75.57
            CFC_MEAN_ERROR = -0.92
            CFC_PPS = CFC_CALIOP + CFC_MEAN_ERROR
            Totpix = 3574
        #datafiles = glob.glob("%s*%s*1.dat" % (BASIC_INDATA_DIR,month))
        #datafiles2 = glob.glob("%s*%s*2.dat" % (BASIC_INDATA_DIR,month))
        
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

        samples_cal = 0
        scenes = len(datafiles)

        for datafile in datafiles:
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
    
        pod_low = float(n_low_low)/float(samples_low_cal)
        far_low = 1-pod_low
        pod_medium = float(n_medium_medium)/float(samples_medium_cal)
        far_medium = 1-pod_medium
        pod_high = float(n_high_high)/float(samples_high_cal)
        far_high = 1-pod_high
    
        low_fraction_pps_rel = float(n_low_low+n_low_medium+n_low_high)/float(samples_tot)*100.0
        low_fraction_cal_rel = float(samples_low_cal)/float(samples_tot)*100.0
        medium_fraction_pps_rel = float(n_medium_medium+n_medium_low+n_medium_high)/float(samples_tot)*100.0
        medium_fraction_cal_rel = float(samples_medium_cal)/float(samples_tot)*100.0
        high_fraction_pps_rel = float(n_high_high+n_high_low+n_high_medium)/float(samples_tot)*100.0
        high_fraction_cal_rel = float(samples_high_cal)/float(samples_tot)*100.0
        frac_fraction_pps_rel = float(n_frac_low+n_frac_medium+n_frac_high)/float(samples_tot)*100.0
        low_fraction_pps_abs = float(n_low_low+n_low_medium+n_low_high)/float(samples_tot)*CFC_PPS
        low_fraction_cal_abs = float(samples_low_cal)/float(samples_tot)*CFC_CALIOP
        medium_fraction_pps_abs = float(n_medium_medium+n_medium_low+n_medium_high)/float(samples_tot)*CFC_PPS
        medium_fraction_cal_abs = float(samples_medium_cal)/float(samples_tot)*CFC_CALIOP
        high_fraction_pps_abs = float(n_high_high+n_high_low+n_high_medium)/float(samples_tot)*CFC_PPS
        high_fraction_cal_abs = float(samples_high_cal)/float(samples_tot)*CFC_CALIOP
        frac_fraction_pps_abs = float(n_frac_low+n_frac_medium+n_frac_high)/float(samples_tot)*CFC_PPS
        
        bias_low=low_fraction_pps_abs-low_fraction_cal_abs
        bias_low_noperc = bias_low/100.0
        bias_medium=medium_fraction_pps_abs-medium_fraction_cal_abs
        bias_medium_noperc = bias_medium/100.0
        bias_high=high_fraction_pps_abs-high_fraction_cal_abs
        bias_high_noperc = bias_high/100.0
        
        #bc-RMS calculations ===============================================================================
        
        n_low_no_no = n_medium_medium+n_medium_high+n_high_medium+n_high_high+n_frac_medium+n_frac_high+\
                        (Totpix-samples_low_cal-samples_medium_cal-samples_high_cal)
        print "Total pixels: ",Totpix
        print "Common cloud-free: ",Totpix-samples_low_cal-samples_medium_cal-samples_high_cal
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
        
        pod_low = 100.0*n_low_yes_yes/samples_low_cal
        pod_medium = 100.0*n_medium_yes_yes/samples_medium_cal
        pod_high = 100.0*n_high_yes_yes/samples_high_cal

        far_low = 100.0*(n_low_yes_no)/(n_low_yes_yes+n_low_yes_no)
        far_medium = 100.0*(n_medium_yes_no)/(n_medium_yes_yes+n_medium_yes_no)
        far_high = 100.0*(n_high_yes_no)/(n_high_yes_yes+n_high_yes_no)

        hitrate_low = 1.0*(n_low_yes_yes+n_low_no_no)/Totpix
        hitrate_medium = 1.0*(n_medium_yes_yes+n_low_no_no)/Totpix
        hitrate_high = 1.0*(n_high_yes_yes+n_high_no_no)/Totpix

        kuipers_low = 1.0*(n_low_no_no*n_low_yes_yes-n_low_yes_no*n_low_no_yes)/ \
                        ((n_low_no_no+n_low_no_yes)*(n_low_yes_no+n_low_yes_yes))
        kuipers_medium = 1.0*(n_medium_no_no*n_medium_yes_yes-n_medium_yes_no*n_medium_no_yes)/ \
                        ((n_medium_no_no+n_medium_no_yes)*(n_medium_yes_no+n_medium_yes_yes))
        kuipers_high = 1.0*(n_high_no_no*n_high_yes_yes-n_high_yes_no*n_high_no_yes)/ \
                        ((n_high_no_no+n_high_no_yes)*(n_high_yes_no+n_high_yes_yes))

        #=============================================================================================
        
        #Finally, calculate distribution of CALIOP categories among the PPS Fractional category

        cal_low_fractional = 100.0*n_frac_low/(n_frac_low+n_frac_medium+n_frac_high)
        cal_medium_fractional = 100.0*n_frac_medium/(n_frac_low+n_frac_medium+n_frac_high)
        cal_high_fractional = 100.0*n_frac_high/(n_frac_low+n_frac_medium+n_frac_high)
        
        #=============================================================================================
        
        #Special analysis: Check composition of false alarm categories for medium and high clouds

        n_medium_yes_no = n_medium_low+n_medium_high+n_medium_clear
        far_medium_low = 100.0*n_medium_low/n_medium_yes_no
        far_medium_high = 100.0*n_medium_high/n_medium_yes_no
        far_medium_clear = 100.0*n_medium_clear/n_medium_yes_no

        n_high_yes_no = n_high_low+n_high_medium+n_high_clear
        far_high_low = 100.0*n_high_low/n_high_yes_no
        far_high_medium = 100.0*n_high_medium/n_high_yes_no
        far_high_clear = 100.0*n_high_clear/n_high_yes_no        

        print "Month is:  %s" % month
        print "Total number of matched scenes is: %s" % scenes
        print "Total number of matched cloud types: %s" % samples_tot
        
        print "Probability of detecting LOW, MEDIUM and HIGH: %f %f %f" % (pod_low,pod_medium,pod_high)
        print "False alarm rate LOW, MEDIUM and HIGH: %f %f %f" % (far_low,far_medium,far_high)
        print "Rel. Fraction LOW for PPS and for CALIOP: %f %f" % (low_fraction_pps_rel,low_fraction_cal_rel)
        print "Rel. Fraction MEDIUM for PPS and for CALIOP: %f %f" % (medium_fraction_pps_rel,medium_fraction_cal_rel)
        print "Rel. Fraction HIGH for PPS and for CALIOP: %f %f" % (high_fraction_pps_rel,high_fraction_cal_rel)
        print "Rel. Fraction FRACTIONAL for PPS: %f" % (frac_fraction_pps_rel)
        print "Abs. Fraction LOW for PPS and for CALIOP: %f %f" % (low_fraction_pps_abs,low_fraction_cal_abs)
        print "Abs. Fraction MEDIUM for PPS and for CALIOP: %f %f" % (medium_fraction_pps_abs,medium_fraction_cal_abs)
        print "Abs. Fraction HIGH for PPS and for CALIOP: %f %f" % (high_fraction_pps_abs,high_fraction_cal_abs)
        print "Abs. Fraction FRACTIONAL for PPS: %f" % (frac_fraction_pps_abs)
        print "Bias Low: %f" % (bias_low)
        print "Bias Medium: %f" % (bias_medium)
        print "Bias High: %f" % (bias_high)
        print "Bc-RMS Low: %f" % (rms_low)
        print "Bc-RMS Medium: %f" % (rms_medium)
        print "Bc_RMS High: %f" % (rms_high)
        print "POD (Low,Medium,High): %f %f %f" %(pod_low,pod_medium,pod_high)
        print "FAR (Low,Medium,High): %f %f %f" %(far_low,far_medium,far_high)
        print "FAR Medium fraction Low,High,Clear: %f %f %f" %(far_medium_low,far_medium_high,far_medium_clear)
        print "FAR High fraction Low,Medium,Clear: %f %f %f" %(far_high_low,far_high_medium,far_high_clear) 
        print "HR (Low,Medium,High): %f %f %f" %(hitrate_low,hitrate_medium,hitrate_high)
        print "KSS (Low,Medium,High): %f %f %f" %(kuipers_low,kuipers_medium,kuipers_high)
        print "CALIOP parts of Fractional (Low, Medium, High): %f %f %f" %(cal_low_fractional,cal_medium_fractional,cal_high_fractional) 
        print "Missclassified clear: %d (low:%d, medium:%d, high:%d)" % (samples_pps_missed,samples_pps_missed_low,samples_pps_missed_medium,samples_pps_missed_high)
        print "Missclassified cloudy: %d (low:%d, medium:%d, high:%d, frac:%d)" % (samples_pps_misclass,samples_pps_misclass_low,samples_pps_misclass_medium,samples_pps_misclass_high,samples_pps_misclass_frac)
        print
        
        result.append("Month is:  %s\n" % month)
        result.append( "Total pixels: %s\n" %Totpix)
        result.append( "Common cloud-free: %s\n" %(Totpix-samples_low_cal-samples_medium_cal-samples_high_cal))
        result.append("Total number of matched scenes is: %s\n" % scenes)
        result.append("Total number of matched cloud types: %s\n" % samples_tot)
        result.append("\n")
        result.append("Probability of detecting LOW, MEDIUM and HIGH: %f %f %f\n" % (pod_low,pod_medium,pod_high))
        result.append("False alarm rate LOW, MEDIUM and HIGH: %f %f %f\n" % (far_low,far_medium,far_high))
        result.append("Rel. Fraction LOW for PPS and for CALIOP: %f %f\n" % (low_fraction_pps_rel,low_fraction_cal_rel))
        result.append("Rel. Fraction MEDIUM for PPS and for CALIOP: %f %f\n" % (medium_fraction_pps_rel,medium_fraction_cal_rel))
        result.append("Rel. Fraction HIGH for PPS and for CALIOP: %f %f\n" % (high_fraction_pps_rel,high_fraction_cal_rel))
        result.append("Rel. Fraction FRACTIONAL for PPS: %f\n" % (frac_fraction_pps_rel))
        result.append("Abs. Fraction LOW for PPS and for CALIOP: %f %f\n" % (low_fraction_pps_abs,low_fraction_cal_abs))
        result.append("Abs. Fraction MEDIUM for PPS and for CALIOP: %f %f\n" % (medium_fraction_pps_abs,medium_fraction_cal_abs))
        result.append("Abs. Fraction HIGH for PPS and for CALIOP: %f %f\n" % (high_fraction_pps_abs,high_fraction_cal_abs))
        result.append("Abs. Fraction FRACTIONAL for PPS: %f\n" % (frac_fraction_pps_abs))
        result.append("Bias Low: %f\n" % (bias_low))
        result.append("Bias Medium: %f\n" % (bias_medium))
        result.append("Bias High: %f\n" % (bias_high))
        result.append("Bc-RMS Low: %f\n" % (rms_low))
        result.append("Bc-RMS Medium: %f\n" % (rms_medium))
        result.append("Bc_RMS High: %f\n" % (rms_high))
        result.append("POD (Low,Medium,High): %f %f %f\n" %(pod_low,pod_medium,pod_high))
        result.append("FAR (Low,Medium,High): %f %f %f\n" %(far_low,far_medium,far_high))
        result.append("FAR Medium fraction Low,High,Clear: %f %f %f\n" %(far_medium_low,far_medium_high,far_medium_clear))
        result.append("FAR High fraction Low,Medium,Clear: %f %f %f\n" %(far_high_low,far_high_medium,far_high_clear) )
        result.append("HR (Low,Medium,High): %f %f %f\n" %(hitrate_low,hitrate_medium,hitrate_high))
        result.append("KSS (Low,Medium,High): %f %f %f\n" %(kuipers_low,kuipers_medium,kuipers_high))
        result.append("CALIOP parts of Fractional (Low, Medium, High): %f %f %f\n" %(cal_low_fractional,cal_medium_fractional,cal_high_fractional))
        result.append("Missclassified clear: %d (low:%d, medium:%d, high:%d)\n" % (samples_pps_missed,samples_pps_missed_low,samples_pps_missed_medium,samples_pps_missed_high))
        result.append("Missclassified cloudy: %d (low:%d, medium:%d, high:%d, frac:%d)\n" % (samples_pps_misclass,samples_pps_misclass_low,samples_pps_misclass_medium,samples_pps_misclass_high,samples_pps_misclass_frac))
        result.append("\n")
    fd=open("%s/cty_results_summary_%s%s-%s%s.dat" %(OUTPUT_DIR, STUDIED_YEAR[0],STUDIED_MONTHS[0],STUDIED_YEAR[-1],STUDIED_MONTHS[-1]),'w')
    fd.writelines(result)
    fd.close()
        
        
