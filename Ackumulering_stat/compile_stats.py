'''
Created on Oct 13, 2010

@author: a001696
'''

from orrb_conf import CASES, MAIN_DATADIR, MAP, RESOLUTION


def compile_stats(results_files):
    import orrb_CFC_stat, orrb_CFC_stat_emissfilt, orrb_CFC_stat_surfaces
    import orrb_CTH_stat, orrb_CTH_stat_emissfilt, orrb_CTH_stat_surfaces
    import orrb_CTY_stat
    
    """Run through all stats."""
    cfc_results = orrb_CFC_stat.do_stats(results_files)
    print("=========== Cloud fraction ============")
    orrb_CFC_stat.printout(cfc_results)
    #orrb_CFC_stat_emissfilt.do_stats()
    #orrb_CFC_stat_surfaces.do_stats()
    
    cth_results = orrb_CTH_stat.do_stats(results_files)
    print("========== Cloud top height ===========")
    orrb_CTH_stat.printout(cth_results)
    #orrb_CTH_stat_emissfilt.do_stats()
    #orrb_CTH_stat_surfaces.do_stats()
    
    cty_results = orrb_CTY_stat.do_stats(results_files, cfc_results)
    print("============= Cloud type ==============")
    orrb_CTY_stat.printout(cty_results)
    
    return cfc_results, cth_results, cty_results


def gather_values_all_months(cfc_results, cth_results):
    all_values = {}
    for results_type, name in [(cfc_results, 'CFC'), (cth_results, 'CTH')]:
        for results in results_type.values():
            for k, v in results.items():
                if type(v) is float:
                    try:
                        all_values['%s.%s' % (name, k)].append(v)
                    except KeyError:
                        all_values['%s.%s' % (name, k)] = [v]
    return all_values, tot_samples


def weighted_averages(vals, num_key):
    """Calculate weighted averages of all values in dict *vals*, based on *num_key*."""
    raise NotImplementedError("This function needs correction and is currently disabled.")
    import numpy
    tot_samples = vals[num_key]
    grand_tot_samples = sum(tot_samples)
    averages = {}
    for k, v_list in vals.items():
        averages[k] = sum(numpy.multiply(v_list, tot_samples)) / grand_tot_samples
    
    return averages

if __name__ == '__main__':
    from glob import glob
    results_files = []
    print("Gathering statistics from all validation results files in the "
          "following directories:")
    for case in CASES:
        basic_indata_dir = "%s/Results/%s/%s/%s/%02d/%s/BASIC" % \
            (MAIN_DATADIR, case['satname'], RESOLUTION[0] , case['year'], case['month'], MAP[0])
        print("-> " + basic_indata_dir)
        results_files.extend(glob("%s/*.dat" % basic_indata_dir))
    print('')
    
    cfc_results, cth_results, cty_results = compile_stats(results_files)
    #all_values, tot_samples = gather_values_all_months(cfc_results, cth_results)
    #averages = weighted_averages(all_values, tot_samples)
    #for k in sorted(averages.keys()):
    #    print("%s: %.3f" % (k, averages[k]))

