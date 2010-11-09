# program process_master.py

from pps_error_messages import write_log
import config


def FindFiles(avhrrfile):
    import file_finders
    
    pps_finder = file_finders.PpsFileFinder(config.PPS_DATA_DIR, time_window=5*60)
    try:
        pps_finder.set_subdir_method(config.subdir)
    except AttributeError:
        pass
    parsed = pps_finder.parse(avhrrfile)
    satname = parsed['satellite']
    datetime = parsed['datetime']
    cloudsat_finder = file_finders.CloudsatFileFinder(config.CLOUDSAT_DIR,
                                                      config.RESOLUTION,
                                                      config.CLOUDSAT_TYPE)
    cloudsat_finder.set_time_window(config.SAT_ORBIT_DURATION + config.sec_timeThr)
    cloudsat_files = sorted(cloudsat_finder.find(datetime))
    calipso_finder = file_finders.CalipsoFileFinder(config.CALIPSO_DIR,
                                                    config.RESOLUTION)
    calipso_finder.set_time_window(config.SAT_ORBIT_DURATION + config.sec_timeThr)
    calipso_files = sorted(calipso_finder.find(datetime))
    cloudtype_file = pps_finder.find(satname, datetime, ending='cloudtype.h5')[0]
    ctth_file = pps_finder.find(satname, datetime,
                                ending='%s.h5' % config.CTTH_FILE)[0]
    avhrr_file = avhrrfile
    nwp_tsur_file = pps_finder.find(satname, datetime, ending='nwp_tsur.h5')[0]
    sunsatangles_file = pps_finder.find(satname, datetime, ending='sunsatangles.h5')[0]
    
    return (cloudsat_files, calipso_files, cloudtype_file, ctth_file, 
            avhrr_file, nwp_tsur_file, sunsatangles_file)


def process_matchups(matchups):
    """Run the given SNO *matchups* through the validation system."""
    problem_files = set()
    no_matchup_files = []
    no_matchup_cross = []
    processed_avhrr_files = set()
    for match in matchups:
        if match.time_window is not None:
            avhrr_finder.set_time_window(-(config.SAT_ORBIT_DURATION + match.time_window), match.time_window)
        else:
            avhrr_finder.set_time_window(-(config.SAT_ORBIT_DURATION + config.sec_timeThr), config.sec_timeThr)
        
        try:
            avhrr_file = avhrr_finder.find(match.satellite1.lower(), match.time1)[0]
        except IndexError:
            write_log('WARNING', "No matching AVHRR file found for %s %s." % \
                      (match.satellite1, match.time1.strftime("%F %T")))
            no_matchup_cross.append(match)
            continue
        if avhrr_file in processed_avhrr_files:
            write_log('WARNING', "Matching AVHRR file has already been processed: %s" % avhrr_file)
            continue
        processed_avhrr_files.add(avhrr_file)
        
        # This is what I would like to do instead...
        # import file_finders
        # pps_finder = file_finders.PpsFileFinder(basedir=???, ending='avhrr.h5')
        #avhrr_file = pps_finder.find(match)
        try:
            cloudsat_file, calipso_file, cloudtype_file, ctth_file, avhrr_file, nwp_tsur_file, sunsatangles_file = FindFiles(avhrr_file)
            if len(cloudsat_file) == 0 or len(calipso_file) == 0:
                raise IndexError # For lack of a more descriptive error
        except IndexError:
            problem_files.add(avhrr_file)
            write_log('WARNING', "Couldn't find all needed files for %s." % avhrr_file)
            if options.debug is True:
                raise
            continue
        
        for mode in run_modes:
            try:
                cloudsat_calipso_avhrr_match.run(cloudsat_file, calipso_file,
                                                 cloudtype_file, ctth_file,
                                                 avhrr_file, nwp_tsur_file,
                                                 sunsatangles_file, mode, 
                                                 resolution)
            except MatchupError, err:
                write_log('WARNING', "Matchup problem for %s %s: %s" % \
                          (match.satellite1, match.time1.strftime("%F %T"),
                           err.message))
                no_matchup_files.append(avhrr_file)
                break
            except:
                problem_files.add(avhrr_file)
                write_log('WARNING', "Couldn't run cloudsat_calipso_avhrr_match.")
                if options.debug is True:
                    raise
                break
    if len(no_matchup_cross) > 0:
        write_log('WARNING', "%d of %d cases had no matching AVHRR files:\n%s" % \
                  (len(no_matchup_cross), len(matchups),
                   '\n'.join([str(c) for c in no_matchup_cross])))
    if len(no_matchup_files) > 0:
        write_log('WARNING', "%d of %d cases had no matchups in region, within the time window:\n%s" % \
                  (len(no_matchup_files), len(matchups), '\n'.join(no_matchup_files)))
    if len(problem_files) > 0:
        write_log('WARNING', "%d of %d cases had unknown problems:\n%s" % \
                  (len(problem_files), len(matchups), '\n'.join(problem_files)))


#------------------------------------------------------------------------------------------------------------------  
if __name__=='__main__':
    #pdb.set_trace()
    from optparse import OptionParser
    import find_crosses
    import file_finders
    import cloudsat_calipso_avhrr_match
    from common import MatchupError
    
    parser = OptionParser()
    parser.set_usage("usage: %prog [options] sno_output_files...\n"
                     "Some influential environment variables:\n"
                     "SAT_DIR                   Base directory where satellite data files are stored.\n"
                     "VALIDATION_RESULTS_DIR    Base directory where results will be stored.\n"
                     "                          Used indirectly by cloudsat_calipso_avhrr_match.py.\n"
                     "CTTH_FILE                 CTTH file type to use (one of 'ctth', 'ctth_opaque, \n"
                     "                          and 'ctth_semitransparent.\n")
    parser.add_option('-M', '--mode', type='string', action='append',
                      help="Run validation software in MODE (valid modes are %s)" % \
                      ', '.join(config.ALLOWED_MODES))
    parser.add_option('-d', '--debug', action='store_true')
    (options, args) = parser.parse_args()
    
    if options.mode is not None:
        run_modes = options.mode
    else:
        run_modes = config.ALLOWED_MODES
    
    if len(args) > 0:
        sno_output_files = args
    else:
        sno_output_files = config.match_files
    
    resolution = "%i.%i" %(config.RESOLUTION,config.clsat_type)
    avhrr_finder = file_finders.PpsFileFinder(config.PPS_DATA_DIR, 'avhrr.h5')
    try:
        avhrr_finder.set_subdir_method(config.subdir)
    except AttributeError:
        pass
    
    matchups = []
    for sno_output_file in sno_output_files:
        found_matchups = find_crosses.parse_crosses_file(sno_output_file)
            
        if len(found_matchups) == 0:
            write_log('WARNING', "No matchups found in SNO output file %s" % sno_output_file)
            if options.debug is True:
                raise Warning
            continue
        else:
            matchups.extend(found_matchups)
    
    process_matchups(matchups)
