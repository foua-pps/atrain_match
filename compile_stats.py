'''
Created on Oct 13, 2010

@author: a001696
'''
from config import (CASES, RESOLUTION, _validation_results_dir, AREA,
                    MIN_OPTICAL_DEPTH, CCI_CLOUD_VALIDATION,
                    DNT_FLAG, CONFIG_PATH, SURFACES,
                    COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC,
                    COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE,
                    ALSO_USE_5KM_FILES)

def compile_stats(results_files, write=False, outfile_cfc=None):
    """Run through all summary statistics."""
    
    print("=========== Cloud fraction ============")
    from statistics import orrb_CFC_stat
    cfc_stats = orrb_CFC_stat.CloudFractionStats(results_files)
    if write:
        cfc_stats.write(outfile_cfc)
    for l in cfc_stats.printout():
        print(l)
    
    compiled_cth_file_name = outfile_cfc.replace('_cfc_','_cth_')
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
                cth_stats.write(compiled_cth_file_name, mode='a')
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
                cth_stats.write(compiled_cth_file_name, mode='a')
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "===== Top Layer, Not optically very thin ==="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [31,32,33,34], note)
            if write:
                cth_stats.write(compiled_cth_file_name, mode='a')
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = ("==== Top Layer, Optically very thin ===\n"
                    "========== Expect bad results here! ===")
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [36,37,38,39], note)
            if write:
                cth_stats.write(compiled_cth_file_name, mode='a')
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
                cth_stats.write(compiled_cth_file_name, mode='a')
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
                cth_stats.write(compiled_cth_file_name, mode='a')
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Opaque Not thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [51,52,53,54], note)
            if write:
                cth_stats.write(compiled_cth_file_name, mode='a')
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Semi-transparent Not thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [56,57,58,59], note)
            if write:
                cth_stats.write(compiled_cth_file_name, mode='a')
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Opaque Thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [61,62,63,64], note)
            if write:
                cth_stats.write(compiled_cth_file_name, mode='a')
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
            note = "========== Semi-transparent Thin top layer ==========="
            from statistics import orrb_CTH_stat
            cth_stats = orrb_CTH_stat.CloudTopStats(results_files, 
                                                    [66,67,68,69], note)
            if write:
                cth_stats.write(compiled_cth_file_name, mode='a')
            cth_stats.printout()
            for l in cth_stats.printout():
                print(l)
    
    if not CCI_CLOUD_VALIDATION:
        print("============= Cloud type ==============")
        compiled_cty_file_name = outfile_cfc.replace('_cfc_','_cty_')
        from statistics import orrb_CTY_stat
        cty_stats = orrb_CTY_stat.CloudTypeStats(results_files, cfc_stats)
        if write:
            cty_stats.write(compiled_cty_file_name)
        cty_stats.printout()
        for l in cty_stats.printout():
            print(l)

if __name__ == '__main__':

    import os
    from glob import glob
    import ConfigParser


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--write', '-w', const=True, nargs='?', required=False,
                      help="Write results to file")
    parser.add_argument('--basic', '-b',  const=True, nargs='?', required=False,  
                      help="Calculate the statistic for mode BASIC")
    parser.add_argument('--standard', '-t', const=True, nargs='?', required=False,  
                      help="Calculate the statistic for mode STANDARD")
    parser.add_argument('--odticfilter', '-o',  const=True, nargs='?', required=False,  
                      help="Calculate the statistic for mode "
                      "OPTICAL_DEPTH_THIN_IS_CLEAR")
    parser.add_argument('--surface_new_way', '-n',  const=True, nargs='?', required=False,
                      help='calculate the statistic for the different '
                      'surfaces with basic method')
    parser.add_argument('--cotfilter', '-c',  const=True, nargs='?', required=False, 
                      help='calculate the statistic for the optical '
                      'thickness filters')    
    (options) = parser.parse_args()

    #Get filename-templates form atrain_match.cfg
    CONF = ConfigParser.ConfigParser()
    CONF.read(os.path.join(CONFIG_PATH, "atrain_match.cfg"))
    config_options = {}
    for option, value in CONF.items('general', raw = True):
        config_options[option] = value
        
    #Find all wanted modes (dnt)    
    modes_list =[]
    New_DNT_FLAG = ['', '_DAY', '_NIGHT', '_TWILIGHT']
    if options.basic == True:
        modes_list.append('BASIC')
    if options.standard == True:
        modes_list.append('STANDARD')
    if options.odticfilter == True: # I prefer this one! /KG
        print('Will calculate statistic for mode OPTICAL_DEPTH_THIN_IS_CLEAR')
        modes_list.append('OPTICAL_DEPTH_THIN_IS_CLEAR')
    if options.surface_new_way == True:
        print('Will calculate statistic for the different surfaces')
        for surface in SURFACES:
            modes_list.append(surface)
        #Add dnt flag to all modes so far
    print modes_list
    modes_dnt_list = []
    for mode in modes_list:
          for dnt in New_DNT_FLAG:
              modes_dnt_list.append(mode + dnt)
    #Treat cotfilter separately as those directories have not dnt-flag at end          
    if options.cotfilter == True:
        print('Will calculate statistic for mode COT-filter')
        for dnt in New_DNT_FLAG:
            for cot in MIN_OPTICAL_DEPTH:
                #modes_list.append("OPTICAL_DEPTH-%0.2_%sf"(dnt,cot))# if like this
                modes_dnt_list.append("OPTICAL_DEPTH%s-%0.2f"%(dnt,cot))

    #For each mode calcualte the statistics            
    for process_mode_dnt in modes_dnt_list:
        print("Gathering statistics from all validation results files in the "
              "following directories:")    
        print(RESOLUTION)
        results_files = []
        #get result files for all cases
        for case in CASES:
            indata_dir = config_options['result_dir'].format(
                val_dir=_validation_results_dir,
                satellite=case['satname'],
                resolution=str(RESOLUTION),
                area=AREA,
                month="%02d"%(case['month']),
                year=case['year'],
                mode=process_mode_dnt,
                min_opt_depth="")    
            print("-> " + indata_dir)
            results_files.extend(glob("%s/*%skm*.dat" %(indata_dir, RESOLUTION)))
        #compile and write results    
        compiled_dir = config_options['compiled_stats_dir'].format(
            val_dir=_validation_results_dir)
        if not os.path.exists(compiled_dir):
            os.makedirs(compiled_dir)
        compiled_file_cfc = config_options['compiled_stats_filename'].format(
            satellite=case['satname'],
            resolution=str(RESOLUTION),
            area=AREA,
            month=case['month'],
            year=case['year'], 
            mode=process_mode_dnt,
            stat_type='cfc',
            min_opt_depth="")
        compiled_file_cfc = os.path.join(compiled_dir, compiled_file_cfc)
        compile_stats(results_files,  write=options.write, outfile_cfc=compiled_file_cfc)




            
