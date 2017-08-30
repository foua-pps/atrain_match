'''
Created on Oct 13, 2010

@author: a001696
'''
from config import (CASES, COMPILE_STATISTICS_TRUTH, RESOLUTION, _validation_results_dir, AREA,
                    MIN_OPTICAL_DEPTH, CCI_CLOUD_VALIDATION,
                    DNT_FLAG, CONFIG_PATH, SURFACES,
                    COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC,
                    COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE,
                    ALSO_USE_5KM_FILES)

import logging
logging.basicConfig(level=logging.INFO) 
logger = logging.getLogger(__name__)

def compile_stats(results_files, write=True, outfile_cfc="merged_sat_file_cfc", truth_sat='calipso'):
    """Run through all summary statistics."""
    
    #print("=========== Cloud fraction ============")
    from statistics import orrb_CFC_stat
    cfc_stats = orrb_CFC_stat.CloudFractionStats(results_files=results_files, truth_sat=truth_sat)
    cfc_stats.write(outfile_cfc)
    
    compiled_cth_file_name = outfile_cfc.replace('_cfc_','_cth_')
    # First all clouds, then split results
    note = "========== Cloud top height ==========="
    from statistics import orrb_CTH_stat
    cth_stats = orrb_CTH_stat.CloudTopStats(ac_data=cfc_stats.ac_data, truth_sat=truth_sat)
    cth_stats.write(compiled_cth_file_name)
    
    if not CCI_CLOUD_VALIDATION:
        #print("============= Cloud type ==============")
        compiled_cty_file_name = outfile_cfc.replace('_cfc_','_cty_')
        from statistics import orrb_CTY_stat
        cty_stats = orrb_CTY_stat.CloudTypeStats(ac_data=cfc_stats.ac_data, truth_sat=truth_sat)
        cty_stats.write(compiled_cty_file_name)

if __name__ == '__main__':

    import os
    from glob import glob
    import ConfigParser


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--nodnt', '-d',  const=True, nargs='?', 
                        required=False,  
                        help="Don't calculate statistics for all DNTs")
    parser.add_argument('--basic', '-b',  const=True, nargs='?', 
                        required=False,  
                        help="Calculate the statistic for mode BASIC")
    parser.add_argument('--standard', '-t', const=True, nargs='?', 
                        required=False,  
                        help="Calculate the statistic for mode STANDARD")
    parser.add_argument('--odticfilter', '-o',  const=True, nargs='?', 
                        required=False,  
                        help="Calculate the statistic for mode "
                        "OPTICAL_DEPTH_THIN_IS_CLEAR")
    parser.add_argument('--surface_new_way', '-n',  const=True, 
                        nargs='?', required=False,
                        help='calculate the statistic for the different '
                        'surfaces')
    parser.add_argument('--cotfilter', '-c',  const=True, 
                        nargs='?', required=False, 
                        help='calculate the statistic for the optical '
                        'thickness filters')    
    parser.add_argument( '--write', '-w', const=True, nargs='?', 
                         required=False, 
                         help="Depricated: write results to file")
    (options) = parser.parse_args()

    #Get filename-templates form atrain_match.cfg
    CONF = ConfigParser.ConfigParser()
    CONF.read(os.path.join(CONFIG_PATH, "atrain_match.cfg"))
    config_options = {}
    for option, value in CONF.items('general', raw = True):
        config_options[option] = value
        
    #Find all wanted modes (dnt)    
    modes_list =[]
    if options.nodnt == True:
        New_DNT_FLAG = ['']
    else:    
        New_DNT_FLAG = ['', '_DAY', '_NIGHT', '_TWILIGHT']
    if options.write == True:
        logger.warning("Results always written to file, -w flag is depricated")
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

    if len(modes_dnt_list)==0 :
        logger.warning("No modes selected!")
        parser.print_help()
    #For each mode calcualte the statistics            
    for process_mode_dnt in modes_dnt_list:
        print("Gathering statistics from all validation results files in the "
              "following directories:")    
        print(RESOLUTION)
        #get result files for all cases
        for truth_sat in COMPILE_STATISTICS_TRUTH:
            results_files = []
            for case in CASES:
                indata_dir = config_options['result_dir'].format(
                    val_dir=_validation_results_dir,
                    satellite=case['satname'],
                    resolution=str(RESOLUTION),
                    area=AREA,
                    month="%02d"%(case['month']),
                    year=case['year'],
                    mode=process_mode_dnt,
                    truth_sat=truth_sat,
                    min_opt_depth="")  
                indata_dir =indata_dir.replace("%d_%H","*")
                print("-> " + indata_dir)
                results_files.extend(glob("%s/*%skm*%s*.dat" %(indata_dir, RESOLUTION, truth_sat.lower())))
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
                truth_sat=truth_sat,
                stat_type='cfc',
                min_opt_depth="")
            compiled_file_cfc = os.path.join(compiled_dir, compiled_file_cfc)
            compile_stats(results_files,  outfile_cfc=compiled_file_cfc,   truth_sat=truth_sat)
    if options.write == True:
        logger.warning("Results always written to file, -w flag is depricated")



            
