'''
Created on Oct 13, 2010

@author: a001696
'''

from setup import CASES, MAIN_DATADIR, MAP, RESOLUTION


def compile_basic_stats(results_files):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    import orrb_CFC_stat
    cfc_results = orrb_CFC_stat.CloudFractionStats(results_files)
    for l in cfc_results.printout():
        print(l)
    
    print("========== Cloud top height ===========")
    import orrb_CTH_stat
    cth_results = orrb_CTH_stat.do_stats(results_files)
    orrb_CTH_stat.printout(cth_results)
    
    print("============= Cloud type ==============")
    import orrb_CTY_stat
    cty_results = orrb_CTY_stat.do_stats(results_files, cfc_results)
    orrb_CTY_stat.printout(cty_results)
    
    return cfc_results, cth_results, cty_results


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
    
    cfc_results, cth_results, cty_results = compile_basic_stats(results_files)

