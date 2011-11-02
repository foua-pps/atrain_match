'''
Created on Oct 13, 2010

@author: a001696
'''
import pdb
import sys
sys.path.append('/home/sm_erjoh/Projects/atrain_match')
from config import CASES, MAIN_DATADIR, MAP, RESOLUTION, COMPILED_STATS_FILENAME

def compile_filtered_stats(results_files, filttype, write=False):
    """Run through all summary statistics."""
    from config import MIN_OPTICAL_DEPTH
    print("=========== Cloud fraction ============")
    import orrb_CFC_stat_emissfilt
    cfc_stats = orrb_CFC_stat_emissfilt.CloudFractionFilteredStats(results_files)
    if write:
        if filttype == 'OPTICAL_DEPTH':
            cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_%s-%.1f.txt' %(filttype, MIN_OPTICAL_DEPTH))
        else:
            cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_%s.txt' %filttype)
    for l in cfc_stats.printout():
        print(l)
    
    print("========== Cloud top height ===========")
    import orrb_CTH_stat_emissfilt
    cth_stats = orrb_CTH_stat_emissfilt.CloudTopFilteredStats(results_files)
    if write:
        if filttype == 'OPTICAL_DEPTH':
            cth_stats.write(COMPILED_STATS_FILENAME + '_cth_%s-%.1f.txt' %(filttype, MIN_OPTICAL_DEPTH))
        else:
            cth_stats.write(COMPILED_STATS_FILENAME + '_cth_%s.txt' %filttype)
    cth_stats.printout()
    for l in cth_stats.printout():
        print(l)
    return cfc_stats, cth_stats

def compile_surface_stats(results_files, surfacetype, write=False):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    import orrb_CFC_stat_surfaces
    cfc_stats = orrb_CFC_stat_surfaces.CloudFractionSurfacesStats(results_files)
    if write:
        cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_%s.txt' %surfacetype)
    for l in cfc_stats.printout():
        print(l)
    
    print("========== Cloud top height ===========")
    import orrb_CTH_stat_surfaces
    cth_stats = orrb_CTH_stat_surfaces.CloudTopSurfacesStats(results_files)
    if write:
        cth_stats.write(COMPILED_STATS_FILENAME + '_cth_%s.txt' %surfacetype)
    cth_stats.printout()
    for l in cth_stats.printout():
        print(l)
    return cfc_stats, cth_stats

def compile_basic_stats(results_files, write=False):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    import orrb_CFC_stat
    cfc_stats = orrb_CFC_stat.CloudFractionStats(results_files)
    if write:
        cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_BASIC.txt')
    for l in cfc_stats.printout():
        print(l)
    
    print("========== Cloud top height ===========")
    import orrb_CTH_stat
    cth_stats = orrb_CTH_stat.CloudTopStats(results_files)
    if write:
        cth_stats.write(COMPILED_STATS_FILENAME + '_cth_BASIC.txt')
    cth_stats.printout()
    for l in cth_stats.printout():
        print(l)
    
    print("============= Cloud type ==============")
    import orrb_CTY_stat
    cty_stats = orrb_CTY_stat.CloudTypeStats(results_files, cfc_stats)
    if write:
        cty_stats.write(COMPILED_STATS_FILENAME + '_cty_BASIC.txt')
    cty_stats.printout()
    for l in cty_stats.printout():
        print(l)
    
    return cfc_stats, cth_stats, cty_stats


if __name__ == '__main__':
    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option('-w', '--write', action='store_true',
                      help="Write results to file")
    parser.add_option('-b', '--basic', action='store_true',
                      help='calculate the statistic for the mode BASIC')
    parser.add_option('-s', '--surface', action='store_true',
                      help='calculate the statistic for the different surfaces')
    parser.add_option('-f', '--filter', action='store_true',
                      help='calculate the statistic for the different filters')
    
    (options, args) = parser.parse_args()
    
    from glob import glob
    results_files = []
    print("Gathering statistics from all validation results files in the "
          "following directories:")
    if options.basic == True:
        results_files = []
        print('Calculate statistic for mode BASIC')
        for case in CASES:
            print(RESOLUTION)
            basic_indata_dir = "%s/Results/%s/%skm/%s/%02d/%s/BASIC" % \
                (MAIN_DATADIR, case['satname'], RESOLUTION , case['year'], case['month'], MAP[0])
            print("-> " + basic_indata_dir)
            results_files.extend(glob("%s/*.dat" % basic_indata_dir))
        cfc_stats, cth_stats, cty_stats = compile_basic_stats(results_files, write=options.write)
        print('')
    if options.surface == True:
        from config import SURFACES
        print('Calculate statistic for the different surfaces')
        for surface in SURFACES:
            results_files = []
            for case in CASES:
                print(RESOLUTION)
                basic_indata_dir = "%s/Results/%s/%skm/%s/%02d/%s/%s" % \
                (MAIN_DATADIR, case['satname'], RESOLUTION , case['year'], case['month'], MAP[0], surface)
                print("-> " + basic_indata_dir)
                results_files.extend(glob("%s/*.dat" % basic_indata_dir))
            cfc_stats, cth_stats = compile_surface_stats(results_files, surface, write=options.write)
        print('')
    if options.filter == True:
        from config import FILTERTYPE, MIN_OPTICAL_DEPTH
        print('Calculate statistic for the different surfaces')
        for filttype in FILTERTYPE:
            results_files = []
            for case in CASES:
                print(RESOLUTION)
                basic_indata_dir = "%s/Results/%s/%skm/%s/%02d/%s/%s" % \
                (MAIN_DATADIR, case['satname'], RESOLUTION , case['year'], case['month'], MAP[0], filttype)
                if filttype == 'OPTICAL_DEPTH':
                    basic_indata_dir = "%s-%.1f" %(basic_indata_dir, MIN_OPTICAL_DEPTH)
                print("-> " + basic_indata_dir)
                results_files.extend(glob("%s/*.dat" % basic_indata_dir))
            cfc_stats, cth_stats = compile_filtered_stats(results_files, filttype, write=options.write)
        print('')
    pdb.set_trace()

