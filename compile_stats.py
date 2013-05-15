'''
Created on Oct 13, 2010

@author: a001696
'''
from config import CASES, MAIN_DATADIR, MAP, RESOLUTION, COMPILED_STATS_FILENAME,\
    DNT_FLAG, CALIPSO_CLOUD_FRACTION,COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC,\
    ALSO_USE_5KM_FILES

def compile_filtered_stats(results_files, filttype, write=False):
    """Run through all summary statistics."""
    print("=========== Cloud fraction ============")
    from statistics import orrb_CFC_stat_emissfilt
    cfc_stats = orrb_CFC_stat_emissfilt.CloudFractionFilteredStats(results_files)
    if write:
        cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_%s.txt' %filttype)
    for l__ in cfc_stats.printout():
        print(l__)
    
    print("========== Cloud top height ===========")
    from statistics import orrb_CTH_stat_emissfilt
    cth_stats = orrb_CTH_stat_emissfilt.CloudTopFilteredStats(results_files)
    if write:
        cth_stats.write(COMPILED_STATS_FILENAME + '_cth_%s.txt' %filttype)
    cth_stats.printout()
    for l__ in cth_stats.printout():
        print(l__)
    return cfc_stats, cth_stats

def compile_surface_stats(results_files, surfacetype, write=False):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    from statistics import orrb_CFC_stat_surfaces
    cfc_stats = orrb_CFC_stat_surfaces.CloudFractionSurfacesStats(results_files)
    if write:
        cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_%s.txt' %surfacetype)
    for l in cfc_stats.printout():
        print(l)
    
    print("========== Cloud top height ===========")
    from statistics import orrb_CTH_stat_surfaces
    cth_stats = orrb_CTH_stat_surfaces.CloudTopSurfacesStats(results_files)
    if write:
        cth_stats.write(COMPILED_STATS_FILENAME + '_cth_%s.txt' %surfacetype)
    cth_stats.printout()
    for l in cth_stats.printout():
        print(l)
    return cfc_stats, cth_stats

def compile_basic_stats(results_files, result_end, write=False):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    from statistics import orrb_CFC_stat
    cfc_stats = orrb_CFC_stat.CloudFractionStats(results_files)
    if write:
        cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_BASIC%s.txt' %result_end)
    for l in cfc_stats.printout():
        print(l)
    
    print("========== Cloud top height ===========")
    from statistics import orrb_CTH_stat
    cth_stats = orrb_CTH_stat.CloudTopStats(results_files, [16,17,18,19])
    if write:
        cth_stats.write(COMPILED_STATS_FILENAME + '_cth_BASIC%s.txt' %result_end)
    cth_stats.printout()
    for l in cth_stats.printout():
        print(l)

    if COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC:
            print("========== Cloud top height ===========")
            print("========== Single Layer     ===========")
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, [21,22,23,24])
            if write:
                cth_stats.write(COMPILED_STATS_FILENAME + '_cth_BASIC%s.txt' %result_end)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)

    if (COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC and
        (ALSO_USE_5KM_FILES or RESOLUTION==5)): 
            print("========== Cloud top height ===========")
            print("=== Single Layer, Not optically thin ==")
            print("========== Expect good results here! ==")
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, [26,27,28,29])
            if write:
                cth_stats.write(COMPILED_STATS_FILENAME + '_cth_BASIC%s.txt' %result_end)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            print("========== Cloud top height ===========")
            print("===== Top Layer, Not optically thin ===")
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, [31,32,33,34])
            if write:
                cth_stats.write(COMPILED_STATS_FILENAME + '_cth_BASIC%s.txt' %result_end)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            print("========== Cloud top height ===========")
            print("==== Top Layer, Optically very thin ===")
            print("========== Expect bad results here! ===")
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, [36,37,38,39])
            if write:
                cth_stats.write(COMPILED_STATS_FILENAME + '_cth_BASIC%s.txt' %result_end)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)

    
    print("============= Cloud type ==============")
    from statistics import orrb_CTY_stat
    cty_stats = orrb_CTY_stat.CloudTypeStats(results_files, cfc_stats)
    if write:
        cty_stats.write(COMPILED_STATS_FILENAME + '_cty_BASIC%s.txt' %result_end)
    cty_stats.printout()
    for l in cty_stats.printout():
        print(l)
    
    return cfc_stats, cth_stats, cty_stats


if __name__ == '__main__':
    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option('-w', '--write', action='store_true',
                      help="Write results to file")
    parser.add_option('-b', '--no-basic', action='store_true',
                      help="don't calculate the statistic for mode BASIC")
    parser.add_option('-s', '--surface', action='store_true',
                      help='calculate the statistic for the different surfaces')
    parser.add_option('-f', '--filter', action='store_true',
                      help='calculate the statistic for the different filters')
    
    (options, args) = parser.parse_args()
    
    from glob import glob
    results_files = []
    print("Gathering statistics from all validation results files in the "
          "following directories:")
    if not options.no_basic:
        print('Calculate statistic for mode BASIC')
        print(RESOLUTION)
        for dnt in DNT_FLAG:
            results_files = []
            for case in CASES:
                basic_indata_dir = "%s/Results/%s/%skm/%s/%02d/%s/BASIC" % \
                    (MAIN_DATADIR, case['satname'], RESOLUTION , case['year'], case['month'], MAP[0])
                if dnt in [None, 'ALL', 'all']:
                    result_end = ''
                else:
                    basic_indata_dir = basic_indata_dir + '_' + dnt
                    result_end = '_%s' %dnt
                print("-> " + basic_indata_dir)
                if CALIPSO_CLOUD_FRACTION == True:
                    results_files.extend(glob("%s/CCF_%skm*.dat" %(basic_indata_dir, RESOLUTION)))
                else:
                    results_files.extend(glob("%s/%skm*.dat" %(basic_indata_dir, RESOLUTION)))
            cfc_stats, cth_stats, cty_stats = compile_basic_stats(results_files, result_end, write=options.write)
            print('')
            
    if options.surface == True:
        from config import SURFACES
        print('Calculate statistic for the different surfaces')
        for surface in SURFACES:
            for dnt in DNT_FLAG:
                results_files = []
                for case in CASES:
                    print(RESOLUTION)
                    basic_indata_dir = "%s/Results/%s/%skm/%s/%02d/%s/%s" % \
                    (MAIN_DATADIR, case['satname'], RESOLUTION , case['year'], case['month'], MAP[0], surface)
                    if dnt in [None, 'ALL', 'all']:
                        result_end = surface
                    else:
                        basic_indata_dir = basic_indata_dir + '_' + dnt
                        result_end = '%s_%s' %(surface, dnt)
                    print("-> " + basic_indata_dir)
                    if CALIPSO_CLOUD_FRACTION == True:
                        results_files.extend(glob("%s/CCF_%skm*.dat" %(basic_indata_dir, RESOLUTION)))
                    else:
                        results_files.extend(glob("%s/%skm*.dat" %(basic_indata_dir, RESOLUTION)))
                cfc_stats, cth_stats = compile_surface_stats(results_files, result_end, write=options.write)
                print('')
    
    if options.filter == True:
        from config import FILTERTYPE, MIN_OPTICAL_DEPTH
        print('Calculate statistic for the different surfaces')
        for filttype in FILTERTYPE:
            for dnt in DNT_FLAG:
                results_files = []
                for case in CASES:
                    print(RESOLUTION)
                    basic_indata_dir = "%s/Results/%s/%skm/%s/%02d/%s/%s" % \
                    (MAIN_DATADIR, case['satname'], RESOLUTION , case['year'], case['month'], MAP[0], filttype)
                    if dnt in [None, 'ALL', 'all']:
                        result_end = filttype
                    else:
                        basic_indata_dir = basic_indata_dir + '_' + dnt
                        result_end = '%s_%s' %(filttype, dnt)
                    if filttype == 'OPTICAL_DEPTH':
                        basic_indata_dir = "%s-%.1f" %(basic_indata_dir, MIN_OPTICAL_DEPTH)
                        result_end = "%s-%.1f" %(result_end, MIN_OPTICAL_DEPTH)
                    print("-> " + basic_indata_dir)
                    if CALIPSO_CLOUD_FRACTION == True:
                        results_files.extend(glob("%s/CCF_%skm*.dat" %(basic_indata_dir, RESOLUTION)))
                    else:
                        results_files.extend(glob("%s/%skm*.dat" %(basic_indata_dir, RESOLUTION)))
                cfc_stats, cth_stats = compile_filtered_stats(results_files, result_end, write=options.write)
                print('')
