# Program orrb_CFC_stat_emissfilt

# This program calculates basic statistics for the cloud amount (CFC) product for
# each month after filtering of t

import string
import math

# -----------------------------------------------------
def do_stats(results_files):
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

    scenes = len(results_files)
    cfc_sum_csa = 0
    cfc_sum_cal = 0
    
    for datafile in results_files:
        current_datafile = open(datafile, "r")
        datalist = current_datafile.readlines()
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
    
    # Check to see that we have some samples...
    if 0 in [samples_csa, samples_cal]:
        raise ValueError("samples_cal: %d, samples_csa: %d" % (samples_cal, samples_csa))
    
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

    results = {}
    results['scenes'] = scenes
    results['samples_csa'] = samples_csa
    results['mean_CFC_csa'] = mean_CFC_csa
    results['bias_csa_perc'] = bias_csa_perc
    results['rms_csa'] = rms_csa
    results['samples_cal'] = samples_cal
    results['mean_CFC_cal'] = mean_CFC_cal
    results['bias_cal_perc'] = bias_cal_perc
    results['rms_cal'] = rms_cal
    results['bias_modis_perc'] = bias_modis_perc
    results['rms_modis'] = rms_modis
    
    return results


def printout(results):
    lines = []
    try:
        lines.append("Month is:  %s" % results['month'])
    except KeyError:
        pass
    lines.append("Total number of matched scenes is: %s" % results['scenes'])
    lines.append("Total number of Cloudsat matched FOVs: %d" % results['samples_csa'])
    lines.append("Mean CFC Cloudsat: %f" % results['mean_CFC_csa'])
    lines.append("Mean error: %f" % results['bias_csa_perc'])
    lines.append("RMS error: %f" % results['rms_csa'])
    lines.append("")
    lines.append("Total number of CALIOP matched FOVs: %d" % results['samples_cal'])
    lines.append("Mean CFC CALIOP: %f" % results['mean_CFC_cal'])
    lines.append("Mean error: %f" % results['bias_cal_perc'])
    lines.append("RMS error: %f" % results['rms_cal'])
    lines.append("Mean error MODIS: %f" % results['bias_modis_perc'])
    lines.append("RMS error MODIS: %f" % results['rms_modis'])
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
        results_dir = "%s/Results/%s/%s/%s/%s/%s/EMISSFILT" % (MAIN_DATADIR, SATELLITE[0], RESOLUTION[0] , STUDIED_YEAR[0], STUDIED_MONTHS[i], MAP[0])
        
        results_files = glob("%s/*.dat" % results_dir)
        results = do_stats(results_files)
        results['month'] = month
        
        lines.extend(printout(results))
    fd=open("%s/cfc_results_emissfilt_summary_%s%s-%s%s.dat" %(OUTPUT_DIR, STUDIED_YEAR[0],STUDIED_MONTHS[0],STUDIED_YEAR[-1],STUDIED_MONTHS[-1]),'w')
    for l in lines:
        fd.writelines(l + '\n')
    fd.close()