# Program orrb_CTH_stat

# This program calculates basic statistics for the cloud top (CTH) product for
# each month

import string

# -----------------------------------------------------
def do_stats(results_files):
    scenes = len(results_files)

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
    for datafile in results_files:
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
    
    results = {}
    results['month'] = month
    results['scenes'] = scenes
    #results['csa_samples'] = csa_samples
    #results['bias_csa'] = bias_csa
    #results['rms_csa'] = rms_csa
    #results['rms_csa'] = rms_csa
    results['cal_all_samples'] = cal_all_samples
    results['cal_low_samples'] = cal_low_samples
    results['cal_medium_samples'] = cal_medium_samples
    results['cal_high_samples'] = cal_high_samples
    results['bias_cal_all'] = bias_cal_all
    results['bias_cal_low'] = bias_cal_low
    results['bias_cal_medium'] = bias_cal_medium
    results['bias_cal_high'] = bias_cal_high
    results['rms_cal_all'] = rms_cal_all
    results['rms_cal_low'] = rms_cal_low
    results['rms_cal_medium'] = rms_cal_medium
    results['rms_cal_high'] = rms_cal_high
    #results['rms_cal'] = rms_cal
    
    return results


def printout(results):
    lines = []
    try:
        lines.append("Month is:  %s" % results['month'])
    except KeyError:
        pass
    lines.append("Month is:  %s" % results['month'])
    lines.append("Total number of matched scenes is: %s" % results['scenes'])
    lines.append("")
    lines.append("Total number of CALIOP matched cloudtops: %d" % results['cal_all_samples'])
    lines.append("Number of CALIOP matched low cloudtops: %d" % results['cal_low_samples'])
    lines.append("Number of CALIOP matched medium cloudtops: %d" % results['cal_medium_samples'])
    lines.append("Number of CALIOP matched high cloudtops: %d" % results['cal_high_samples'])
    lines.append("Mean error total cases: %f" % results['bias_cal_all'])
    lines.append("Mean error low-level cases: %f" % results['bias_cal_low'])
    lines.append("Mean error medium-level cases: %f" % results['bias_cal_medium'])
    lines.append("Mean error high-level cases: %f" % results['bias_cal_high'])
    lines.append("Weighted RMS error total cases: %f" % results['rms_cal_all'])
    lines.append("Weighted RMS error low-level cases: %f" % results['rms_cal_low'])
    lines.append("Weighted RMS error medium-level cases: %f" % results['rms_cal_medium'])
    lines.append("Weighted RMS error high-level cases: %f" % results['rms_cal_high'])
    lines.append("")
    
    for l in lines:
        print(l)
    
    return lines


if __name__ == "__main__":
    from glob import glob
    from orrb_conf import SATELLITE, RESOLUTION, STUDIED_YEAR, STUDIED_MONTHS, MAP, MAIN_DATADIR, OUTPUT_DIR, SURFACES
    lines = []
    
    for surface in SURFACES:
        lines.append("Surface is: %s" % surface)
        lines.append("")
        for i in range(len(STUDIED_MONTHS)):
            month="%s%s" %(STUDIED_YEAR[0],STUDIED_MONTHS[i])
            results_dir = "%s/Results/%s/%s/%s/%s/%s/%s" % (MAIN_DATADIR, SATELLITE[0], RESOLUTION[0] , STUDIED_YEAR[0], STUDIED_MONTHS[i], MAP[0],surface)
            
            results_files = glob("%s/*.dat" % results_dir)
            results = do_stats(results_files)
            results['month'] = month
            
            lines.extend(printout(results))
        lines.append("")
        lines.append("")
    lines.append("")
    fd=open("%s/cth_results_surfaces_summary_%s%s-%s%s.dat" %(OUTPUT_DIR, STUDIED_YEAR[0],STUDIED_MONTHS[0],STUDIED_YEAR[-1],STUDIED_MONTHS[-1]),'w')
    for l in lines:
        fd.writelines(l + '\n')
    fd.close()