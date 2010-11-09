'''
Created on Oct 13, 2010

@author: a001696
'''

from config import CASES, MAIN_DATADIR, MAP, RESOLUTION, COMPILED_STATS_FILENAME


def compile_basic_stats(results_files, write=False):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    import orrb_CFC_stat
    cfc_stats = orrb_CFC_stat.CloudFractionStats(results_files)
    if write:
        cfs_stats.write(COMPILED_STATS_FILENAME + '_cfc')
    for l in cfc_stats.printout():
        print(l)
    
    print("========== Cloud top height ===========")
    import orrb_CTH_stat
    cth_stats = orrb_CTH_stat.CloudTopStats(results_files)
    if write:
        cfs_stats.write(COMPILED_STATS_FILENAME + '_cth')
    cth_stats.printout()
    for l in cth_stats.printout():
        print(l)
    
    print("============= Cloud type ==============")
    import orrb_CTY_stat
    cty_stats = orrb_CTY_stat.CloudTypeStats(results_files, cfc_stats)
    if write:
        cfs_stats.write(COMPILED_STATS_FILENAME + '_cty')
    cty_stats.printout()
    for l in cty_stats.printout():
        print(l)
    
    return cfc_stats, cth_stats, cty_stats


if __name__ == '__main__':
    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option('-w', '--write', action='store_true',
                      help="Write results to file")
    (options, args) = parser.parse_args()
    
    from glob import glob
    results_files = []
    print("Gathering statistics from all validation results files in the "
          "following directories:")
    for case in CASES:
        basic_indata_dir = "%s/Results/%s/%skm/%s/%02d/%s/BASIC" % \
            (MAIN_DATADIR, case['satname'], RESOLUTION[0] , case['year'], case['month'], MAP[0])
        print("-> " + basic_indata_dir)
        results_files.extend(glob("%s/*.dat" % basic_indata_dir))
    print('')
    
    cfc_stats, cth_stats, cty_stats = compile_basic_stats(results_files, write=options.write)

