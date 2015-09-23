'''
Created on Oct 13, 2010

@author: a001696
'''
from config import (CASES, MAIN_DATADIR, MAP, RESOLUTION, 
                    COMPILED_STATS_FILENAME,
                    MIN_OPTICAL_DEPTH, CCI_CLOUD_VALIDATION,
                    DNT_FLAG, CALIPSO_CLOUD_FRACTION,
                    COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC,
                    COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE,
                    ALSO_USE_5KM_FILES, RESULT_DIR)


def compile_stats(results_files, result_end, write=False, name_in_file='BASIC'):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    from statistics import orrb_CFC_stat
    cfc_stats = orrb_CFC_stat.CloudFractionStats(results_files)
    if write:
        cfc_stats.write(COMPILED_STATS_FILENAME + '_cfc_%s%s.txt' %(name_in_file,result_end))
    for l in cfc_stats.printout():
        print(l)
    
    compiled_cth_file_name = COMPILED_STATS_FILENAME + '_cth_%s%s.txt' %(name_in_file,result_end) 
    # First all clouds, then split results
    note = "========== Cloud top height ==========="
    from statistics import orrb_CTH_stat
    cth_stats = orrb_CTH_stat.CloudTopStats(results_files, [16,17,18,19], note)
    if write:
        cth_stats.write(compiled_cth_file_name)
    cth_stats.printout()
    for l in cth_stats.printout():
        print(l)

    if COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC:
            note = "========== Single Layer     ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [21,22,23,24], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)

    if (COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC and
        (ALSO_USE_5KM_FILES or RESOLUTION==5)): 
            note =("=== Single Layer, Not optically thin ==\n"
                   "========== Expect good results here! ==")
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [26,27,28,29], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "===== Top Layer, Not optically very thin ==="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [31,32,33,34], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = ("==== Top Layer, Optically very thin ===\n"
                    "========== Expect bad results here! ===")
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [36,37,38,39], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
    if COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE:
        try:
            note = "========== Opaque All     ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [41,42,43,44], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            DID_FIND_SEMI_OPAQUE_RESULTS = True     
        except IndexError:
            DID_FIND_SEMI_OPAQUE_RESULTS = False
            print "Can not find opaque/semi ctth results"

        if (COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE and
            DID_FIND_SEMI_OPAQUE_RESULTS):          
            note = "========== Semi-transparent All ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [46,47,48,49], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Opaque Not thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [51,52,53,54], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Semi-transparent Not thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [56,57,58,59], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Opaque Thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [61,62,63,64], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Semi-transparent Thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [66,67,68,69], note)
            if write:
                cth_stats.write(compiled_cth_file_name)
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)


    
    if not CCI_CLOUD_VALIDATION:
        print("============= Cloud type ==============")
        from statistics import orrb_CTY_stat
        cty_stats = orrb_CTY_stat.CloudTypeStats(results_files, cfc_stats)
        if write:
            cty_stats.write(COMPILED_STATS_FILENAME + '_cty_%s%s.txt' %(
                name_in_file,result_end))
        cty_stats.printout()
        for l in cty_stats.printout():
            print(l)
        return cfc_stats, cth_stats, cty_stats
    else:
        return cfc_stats, cth_stats



if __name__ == '__main__':
    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option('-w', '--write', action='store_true',
                      help="Write results to file")
    parser.add_option('-b', '--basic', action='store_true',  
                      help="Calculate the statistic for mode BASIC")
    parser.add_option('-t', '--standard', action='store_true',  
                      help="Calculate the statistic for mode STANDARD")
    parser.add_option('-o', '--odticfilter', action='store_true',  
                      help="Calculate the statistic for mode OPTICAL_DEPTH_THIN_IS_CLEAR")
    parser.add_option('-n', '--surface_new_way', action='store_true',
                      help='calculate the statistic for the different surfaces with basic method')
    parser.add_option('-c', '--cotfilter', action='store_true', # Added for handling a list of cot values/KG
                      help='calculate the statistic for the optical thickness filters')
    
    (options, args) = parser.parse_args()
    
    from glob import glob
    results_files = []
    print("Gathering statistics from all validation results files in the "
          "following directories:")
    #if not options.no_basic:
    if options.basic == True: # I prefer this one! /KG
        print('Calculate statistic for mode BASIC')
        print(RESOLUTION)
        for dnt in DNT_FLAG:
            results_files = []
            for case in CASES:
                basic_indata_dir = "%s/%s/%s/%skm/%s/%02d/%s/BASIC" % \
                    (MAIN_DATADIR, RESULT_DIR, case['satname'],  RESOLUTION , 
                     case['year'], case['month'], MAP[0])
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
            if not CCI_CLOUD_VALIDATION:
                cfc_stats, cth_stats, cty_stats = compile_stats(results_files, 
                                                                      result_end, 
                                                                      write=options.write)
            else:
                cfc_stats, cth_stats = compile_stats(results_files, result_end, write=options.write)
            print('')

    if options.standard == True: # I prefer this one! /KG
        print('Calculate statistic for mode STANDARD')
        print(RESOLUTION)
        for dnt in DNT_FLAG:
            results_files = []
            for case in CASES:
                basic_indata_dir = "%s/%s/%s/%skm/%s/%02d/%s/STANDARD" % \
                    (MAIN_DATADIR, RESULT_DIR, case['satname'],  RESOLUTION , case['year'], case['month'], MAP[0])
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
            if not CCI_CLOUD_VALIDATION:
                cfc_stats, cth_stats, cty_stats = compile_stats(results_files, 
                                                                      result_end, 
                                                                      write=options.write, 
                                                                      name_in_file='STANDARD')
            else:
                cfc_stats, cth_stats = compile_stats(results_files, 
                                                           result_end, 
                                                           write=options.write, 
                                                           name_in_file='STANDARD')        
            print('')

    if options.odticfilter == True: # I prefer this one! /KG
        print('Calculate statistic for mode OPTICAL_DEPTH_THIN_IS_CLEAR')
        print(RESOLUTION)
        for dnt in DNT_FLAG:
            results_files = []
            for case in CASES:
                basic_indata_dir = "%s/%s/%s/%skm/%s/%02d/%s/OPTICAL_DEPTH_THIN_IS_CLEAR" % \
                    (MAIN_DATADIR, RESULT_DIR, case['satname'],  RESOLUTION , 
                     case['year'], case['month'], MAP[0])
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
            if not CCI_CLOUD_VALIDATION:
                cfc_stats, cth_stats, cty_stats = compile_stats(
                    results_files, result_end, 
                    write=options.write, name_in_file='OPTICAL_DEPTH_THIN_IS_CLEAR')
            else:
                cfc_stats, cth_stats = compile_stats(
                    results_files, result_end, 
                    write=options.write, name_in_file='OPTICAL_DEPTH_THIN_IS_CLEAR')        
            print('')

    if options.surface_new_way == True:
        from config import SURFACES
        print('Calculate statistic for the different surfaces with basic method')
        for surface in SURFACES:
            for dnt in DNT_FLAG:
                results_files = []
                for case in CASES:
                    print(RESOLUTION)
                    basic_indata_dir = "%s/%s/%s/%skm/%s/%02d/%s/%s" % \
                    (MAIN_DATADIR, RESULT_DIR, case['satname'], RESOLUTION , 
                     case['year'], case['month'], MAP[0], surface)
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
                if not CCI_CLOUD_VALIDATION:
                    cfc_stats, cth_stats, cty_stats = compile_stats(
                        results_files, result_end, 
                        write=options.write, name_in_file=surface)
                else:
                    cfc_stats, cth_stats = compile_stats(
                        results_files, result_end, write=options.write, 
                        name_in_file=surface) 
                print('')            
    

    if options.cotfilter == True:
        print('Calculate statistic for mode COT-filter')
        print(RESOLUTION)

        New_DNT_FLAG = ['ALL', 'DAY', 'NIGHT', 'TWILIGHT']
        
        #for dnt in DNT_FLAG:
        for dnt in New_DNT_FLAG:
            for cot in MIN_OPTICAL_DEPTH:
                results_files = []
                for case in CASES:
                    basic_indata_dir = "%s/%s/%s/%skm/%s/%02d/%s/OPTICAL_DEPTH-%0.2f" % \
                    (MAIN_DATADIR, RESULT_DIR, case['satname'], RESOLUTION , 
                     case['year'], case['month'], MAP[0], cot)
                    if dnt in [None, 'ALL', 'all']:
                        result_end = '_%0.2f' % cot
                    else:
                        basic_indata_dir = "%s/%s/%s/%skm/%s/%02d/%s/OPTICAL_DEPTH_%s-%0.2f" % \
                                           (MAIN_DATADIR, RESULT_DIR, case['satname'], 
                                            RESOLUTION , case['year'], case['month'], MAP[0], dnt,cot)
                        #basic_indata_dir = basic_indata_dir + '_' + dnt
                        result_end = '_%s_%0.2f' % (dnt,cot)
                    print("-> " + basic_indata_dir)
                    if CALIPSO_CLOUD_FRACTION == True:
                        results_files.extend(glob("%s/CCF_%skm*.dat" %(basic_indata_dir, RESOLUTION)))
                        print "Length of result file list: ",len(results_files)
                    else:
                        results_files.extend(glob("%s/%skm*.dat" %(basic_indata_dir, RESOLUTION)))
                        print "Length of result file list: ",len(results_files)
                if not CCI_CLOUD_VALIDATION:
                    cfc_stats, cth_stats, cty_stats = compile_stats(results_files, 
                                                                    result_end, 
                                                                    write=options.write, 
                                                                    name_in_file='cotfilter')
                else:
                    cfc_stats, cth_stats = compile_stats(results_files, 
                                                                   result_end, 
                                                                   write=options.write, 
                                                                   name_in_file='cotfilter')        
                print('')
                        
             
            print('')
            
