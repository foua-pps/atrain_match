#!/usr/bin/python
"""
This module is the main entry point for matching noaa and/or metop data with
Cloudsat and Calipso data, and produce statistics for validation of PPS
cloud mask, cloud type and cloud top temperature and height products.

It is really a wrapper to :func:`cloudsat_calipso_avhrr_match.run`, for running
through a set of SNO matchups.

.. note::

    This module has a command line interface.

"""

from pps_error_messages import write_log
import config

def process_matchups(matchups, run_modes, reprocess=False, debug=False):
    """
    Run the given SNO *matchups* through the validation system.
    
    *matchups* should be a list of :class:`find_crosses.Cross` instances.
    
    If *reprocess* is True, disregard any previously generated Cloudsat- and
    Calipso-AVHRR matchup files.
    
    """
    import cloudsat_calipso_avhrr_match
    from common import MatchupError
    
    problematic = set()
    no_matchup_files = []
    for match in sorted(matchups):
#        match = matchups[i + 50]
        for mode in run_modes:
            try:
                cloudsat_calipso_avhrr_match.run(match, mode, reprocess)
            except MatchupError, err:
                write_log('WARNING', "Matchup problem: %s" % str(err))
                no_matchup_files.append(match)
                break
            except:
                problematic.add(match)
                write_log('WARNING', "Couldn't run cloudsat_calipso_avhrr_match.")
                if debug is True:
                    raise
                break
    if len(no_matchup_files) > 0:
        write_log('WARNING', "%d of %d cases had no matchups in region, within the time window:\n%s" % \
                  (len(no_matchup_files), len(matchups),
                   '\n'.join([str(m) for m in no_matchup_files])))
    if len(problematic) > 0:
        write_log('WARNING', "%d of %d cases had unknown problems:\n%s" % \
                  (len(problematic), len(matchups),
                   '\n'.join([str(m) for m in problematic])))
    

def main(args=None):
    """
    Process command line options and run matchup and validation.
    
    If *args* is provided, it should be a list of command line arguments (e.g.
    sys.argv[1:]).
    
    For a complete usage description, run 'python process_master -h'.
    
    """
    from optparse import OptionParser
    import find_crosses
    
    parser = OptionParser()
    parser.set_usage("usage: %prog [options] sno_output_files...\n"
                     "Some influential environment variables:\n"
                     "SAT_DIR                   Base directory where satellite data files are stored.\n"
                     "VALIDATION_RESULTS_DIR    Base directory where results will be stored.\n"
                     "                          Used indirectly by cloudsat_calipso_avhrr_match.py.\n"
                     "CTTH_FILE                 CTTH file type to use (one of 'ctth', 'ctth_opaque, \n"
                     "                          and 'ctth_semitransparent.\n")
    parser.add_option('-M', '--mode', type='choice', action='append', choices=config.ALLOWED_MODES,
                      help="Run validation software in MODE (valid modes are %s)" % \
                      ', '.join(config.ALLOWED_MODES))
    parser.add_option('-r', '--reprocess', action='store_true', default=False,
                      help="Disregard any previously generated Cloudsat- and "
                      "Calipso-AVHRR matchup files.")
    parser.add_option('-d', '--debug', action='store_true', default=False)
    (options, args) = parser.parse_args(args)
    
    if options.mode is not None:
        run_modes = options.mode
    else:
        run_modes = config.ALLOWED_MODES
    
    if len(args) > 0:
        sno_output_files = args
    else:
        parser.error("No snotimes output files provided")
    
    config.DEBUG = options.debug
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
    process_matchups(matchups, run_modes, options.reprocess, options.debug)
    
    return 0

#------------------------------------------------------------------------------------------------------------------  
if __name__=='__main__':
    main()
