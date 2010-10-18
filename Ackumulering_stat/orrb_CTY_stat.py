# Program orrb_CTY_stat

# This program calculates basic statistics for the cloud type (CTY) product for
# each month

import string
import math

# -----------------------------------------------------
def do_stats(results_files, cfc_results=None):
    if cfc_results is not None:
        CFC_CALIOP = cfc_results['mean_CFC_cal']
        CFC_MEAN_ERROR = cfc_results['bias_cal_perc']
        Totpix = cfc_results['samples_cal']
    else:
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
            Totpix = 40212
        elif month == "200809": 
            CFC_CALIOP = 75.57
            CFC_MEAN_ERROR = -0.92
            Totpix = 3574
        elif month == "200907": 
            CFC_CALIOP = 67.121972
            CFC_MEAN_ERROR = -9.289297
            Totpix = 56226
    # The above is crazy!!! I will keep it for now, but extend the functionality
    # to make it possible to pass the values to the function.
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

    scenes = len(results_files)

    for datafile in results_files:
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

    results = {'Totpix': Totpix,
               'common_cloud_free': common_cloud_free,
               'scenes': scenes,
               'samples_tot': samples_tot,
               'pod_low': pod_low,
               'pod_medium': pod_medium,
               'pod_high': pod_high,
               'far_low': far_low,
               'far_medium': far_medium,
               'far_high': far_high,
               'low_fraction_pps_rel': low_fraction_pps_rel,
               'low_fraction_cal_rel': low_fraction_cal_rel,
               'medium_fraction_pps_rel': medium_fraction_pps_rel,
               'medium_fraction_cal_rel': medium_fraction_cal_rel,
               'high_fraction_pps_rel': high_fraction_pps_rel,
               'high_fraction_cal_rel': high_fraction_cal_rel,
               'frac_fraction_pps_rel': frac_fraction_pps_rel,
               'low_fraction_pps_abs': low_fraction_pps_abs,
               'low_fraction_cal_abs': low_fraction_cal_abs,
               'medium_fraction_pps_abs': medium_fraction_pps_abs,
               'medium_fraction_cal_abs': medium_fraction_cal_abs,
               'high_fraction_pps_abs': high_fraction_pps_abs,
               'high_fraction_cal_abs': high_fraction_cal_abs,
               'frac_fraction_pps_abs': frac_fraction_pps_abs,
               'bias_low': bias_low,
               'bias_medium': bias_medium,
               'bias_high': bias_high,
               'rms_low': rms_low,
               'rms_medium': rms_medium,
               'rms_high': rms_high,
               'pod_low': pod_low,
               'pod_medium': pod_medium,
               'pod_high': pod_high,
               'far_low': far_low,
               'far_medium': far_medium,
               'far_high': far_high,
               'far_medium_low': far_medium_low,
               'far_medium_high': far_medium_high,
               'far_medium_clear': far_medium_clear,
               'far_high_low': far_high_low,
               'far_high_medium': far_high_medium,
               'far_high_clear': far_high_clear, 
               'hitrate_low': hitrate_low,
               'hitrate_medium': hitrate_medium,
               'hitrate_high': hitrate_high,
               'kuipers_low': kuipers_low,
               'kuipers_medium': kuipers_medium,
               'kuipers_high': kuipers_high,
               'cal_low_fractional': cal_low_fractional,
               'cal_medium_fractional': cal_medium_fractional,
               'cal_high_fractional': cal_high_fractional, 
               'samples_pps_missed': samples_pps_missed,
               'samples_pps_missed_low': samples_pps_missed_low,
               'samples_pps_missed_medium': samples_pps_missed_medium,
               'samples_pps_missed_high': samples_pps_missed_high,
               'samples_pps_misclass': samples_pps_misclass,
               'samples_pps_misclass_low': samples_pps_misclass_low,
               'samples_pps_misclass_medium': samples_pps_misclass_medium,
               'samples_pps_misclass_high': samples_pps_misclass_high,
               'samples_pps_misclass_frac': samples_pps_misclass_frac}
    
    return results

def printout(results):
    lines = []
    try:
        lines.append("Month is:  %s" % results['month'])
    except KeyError:
        pass
    lines.append( "Total pixels: %s" % results['Totpix'])
    lines.append( "Common cloud-free: %s" % results['common_cloud_free'])
    lines.append("Total number of matched scenes is: %s" % results['scenes'])
    lines.append("Total number of matched cloud types: %s" % results['samples_tot'])
    lines.append("")
    lines.append("Probability of detecting LOW, MEDIUM and HIGH: %f %f %f" % (results['pod_low'], results['pod_medium'], results['pod_high']))
    lines.append("False alarm rate LOW, MEDIUM and HIGH: %f %f %f" % (results['far_low'], results['far_medium'], results['far_high']))
    lines.append("Rel. Fraction LOW for PPS and for CALIOP: %f %f" % (results['low_fraction_pps_rel'], results['low_fraction_cal_rel']))
    lines.append("Rel. Fraction MEDIUM for PPS and for CALIOP: %f %f" % (results['medium_fraction_pps_rel'], results['medium_fraction_cal_rel']))
    lines.append("Rel. Fraction HIGH for PPS and for CALIOP: %f %f" % (results['high_fraction_pps_rel'], results['high_fraction_cal_rel']))
    lines.append("Rel. Fraction FRACTIONAL for PPS: %f" % results['frac_fraction_pps_rel'])
    lines.append("Abs. Fraction LOW for PPS and for CALIOP: %f %f" % (results['low_fraction_pps_abs'], results['low_fraction_cal_abs']))
    lines.append("Abs. Fraction MEDIUM for PPS and for CALIOP: %f %f" % (results['medium_fraction_pps_abs'], results['medium_fraction_cal_abs']))
    lines.append("Abs. Fraction HIGH for PPS and for CALIOP: %f %f" % (results['high_fraction_pps_abs'], results['high_fraction_cal_abs']))
    lines.append("Abs. Fraction FRACTIONAL for PPS: %f" % results['frac_fraction_pps_abs'])
    lines.append("Bias Low: %f" % results['bias_low'])
    lines.append("Bias Medium: %f" % results['bias_medium'])
    lines.append("Bias High: %f" % results['bias_high'])
    lines.append("Bc-RMS Low: %f" % results['rms_low'])
    lines.append("Bc-RMS Medium: %f" % results['rms_medium'])
    lines.append("Bc_RMS High: %f" % results['rms_high'])
    lines.append("POD (Low,Medium,High): %f %f %f" % (results['pod_low'], results['pod_medium'], results['pod_high']))
    lines.append("FAR (Low,Medium,High): %f %f %f" % (results['far_low'], results['far_medium'], results['far_high']))
    lines.append("FAR Medium fraction Low,High,Clear: %f %f %f" % (results['far_medium_low'], results['far_medium_high'], results['far_medium_clear']))
    lines.append("FAR High fraction Low,Medium,Clear: %f %f %f" % (results['far_high_low'], results['far_high_medium'], results['far_high_clear']) )
    lines.append("HR (Low,Medium,High): %f %f %f" % (results['hitrate_low'], results['hitrate_medium'], results['hitrate_high']))
    lines.append("KSS (Low,Medium,High): %f %f %f" % (results['kuipers_low'], results['kuipers_medium'], results['kuipers_high']))
    lines.append("CALIOP parts of Fractional (Low, Medium, High): %f %f %f" % (results['cal_low_fractional'], results['cal_medium_fractional'], results['cal_high_fractional']))
    lines.append("Missclassified clear: %d (low:%d, medium:%d, high:%d)" % (results['samples_pps_missed'], results['samples_pps_missed_low'], results['samples_pps_missed_medium'], results['samples_pps_missed_high']))
    lines.append("Missclassified cloudy: %d (low:%d, medium:%d, high:%d, frac:%d)" % (results['samples_pps_misclass'], results['samples_pps_misclass_low'], results['samples_pps_misclass_medium'], results['samples_pps_misclass_high'], results['samples_pps_misclass_frac']))
    lines.append("")
    
    for l in lines:
        print(l)
    
    return lines


if __name__ == "__main__":
    from glob import glob
    from orrb_conf import SATELLITE, RESOLUTION, STUDIED_YEAR, STUDIED_MONTHS, MAP, MAIN_DATADIR, OUTPUT_DIR
    lines = []
    
    for i in range(len(STUDIED_MONTHS)):
        month=("%s%s" %(STUDIED_YEAR[0],STUDIED_MONTHS[i]))
        results_dir = "%s/Results/%s/%s/%s/%s/%s/BASIC" % (MAIN_DATADIR, SATELLITE[0], RESOLUTION[0] , STUDIED_YEAR[0], STUDIED_MONTHS[i], MAP[0])
        
        results_files = glob("%s/*.dat" % results_dir)
        results = do_stats(results_files)
        results['month'] = month
        
        lines.extend(printout(results))
    fd=open("%s/cty_results_summary_%s%s-%s%s.dat" %(OUTPUT_DIR, STUDIED_YEAR[0],STUDIED_MONTHS[0],STUDIED_YEAR[-1],STUDIED_MONTHS[-1]),'w')
    for l in lines:
        fd.writelines(l + '\n')
    fd.close()