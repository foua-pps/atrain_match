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
import logging
logging.basicConfig(level=logging.INFO) 
logger = logging.getLogger(__name__)

import config
from common import InputError

def process_matchups(matchups, run_modes, reprocess=False, debug=False):
    """
    Run the given SNO *matchups* through the validation system.
    
    *matchups* should be a list of :class:`find_crosses.Cross` instances.
    
    If *reprocess* is True, disregard any previously generated Cloudsat- and
    Calipso-AVHRR matchup files.
    
    """
    import cloudsat_calipso_avhrr_match
    from common import MatchupError
    import os
    import ConfigParser
    from config import CONFIG_PATH
    CONF = ConfigParser.ConfigParser()
    config_file = os.path.join(CONFIG_PATH, "atrain_match.cfg")
    if not os.path.isfile(config_file):
        raise IOError("Couldn't find config file %s."%(config_file))
    CONF.read(config_file)
    OPTIONS = {}    
    for option, value in CONF.items('general', raw = True):
        OPTIONS[option] = value

    problematic = set()
    no_matchup_files = []
    outstatus = 0
    for match in sorted(matchups):
        try:
            cloudsat_calipso_avhrr_match.run(match, run_modes,  OPTIONS, reprocess)
        except MatchupError, err:
            logger.warning("Matchup problem: %s" % str(err))
            import traceback
            traceback.print_exc()
            no_matchup_files.append(match)
            outstatus = 5
        except:
            import traceback
            outstatus = 5
            traceback.print_exc()
            problematic.add(match)
            logger.warning("Couldn't run cloudsat_calipso_avhrr_match.")
            if debug is True:
                raise
                
    if len(no_matchup_files) > 0:
        logger.warning(
                  "%d of %d cases had no matchups in region, within the time window:\n%s" % \
                  (len(no_matchup_files), len(matchups),
                   '\n'.join([str(m) for m in no_matchup_files])))
    if len(problematic) > 0:
        logger.warning("%d of %d cases had unknown problems:\n%s" % \
                  (len(problematic), len(matchups),
                   '\n'.join([str(m) for m in problematic])))    
    return outstatus

def main(args=None):
    """
    Process command line options and run matchup and validation.
    
    If *args* is provided, it should be a list of command line arguments (e.g.
    sys.argv[1:]).
    
    For a complete usage description, run 'python process_master -h'.
    
    """
    import find_crosses
    import argparse
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('--mode', '-M', type=str, required=False, choices=config.ALLOWED_MODES,
                      help=("Run validation software in MODE "))
    parser.add_argument('--reprocess', '-r', const=True, nargs='?', required=False,
                        help="Disregard any previously generated Cloudsat- and "
                        "Calipso-AVHRR matchup files.")
    parser.add_argument('-d', '--debug', const=True, nargs='?', required=False, 
                        help="Get debug logging")
    group.add_argument( '--pps_okay_scene', '-os', 
                      help="Interpret arguments as PPS okay scenes instead of "
                      "sno_output_files (e.g. noaa19_20101201_1345_27891*)")
    group.add_argument( '--pps_product_file', '-pf', 
                      help="Interpret arguments as inputfile with "  
                      "list of pps files")
    group.add_argument( '--cci_product_file', '-cf', 
                      help="Interpret arguments as inputfile with "  
                      "list of cci files")
    group.add_argument( '--maia_product_file', '-mf', 
                      help="Interpret arguments as inputfile with "  
                      "list of maia files")
    group.add_argument('--sno_file', '-sf', 
                      help="Interpret arguments as sno_output_file")
    group.add_argument('--reshaped_product_file', '-rf', 
                      help="Interpret arguments as reshaped_output_file")
    options = parser.parse_args()
    
    if options.mode is not None:
        run_modes = [options.mode]
    else:
        run_modes = config.ALLOWED_MODES

    reprocess = False    
    if options.reprocess is not None:
        reprocess = options.reprocess

    config.DEBUG = options.debug
    if options.debug:
        import logging
        logging.getLogger().setLevel(logging.DEBUG)

    matchups = []
    if options.pps_okay_scene:
        # Simulate crosses from PPS scenes
        from find_crosses import Cross
        from runutils import parse_scene
        scene = options.pps_okay_scene
        satname, time, orbit = parse_scene(scene) #@UnusedVariable
        matchups.append(Cross(satname, '', time, time, -999, -999))
    elif options.sno_file is not None:
        if config.USE_ORBITS_THAT_STARTS_EXACTLY_AT_CROSS:
            #print config.USE_ORBITS_THAT_STARTS_EXACTLY_AT_CROSS
            raise InputError(
                " Expected " +
                "USE_ORBITS_THAT_STARTS_EXACTLY_AT_CROSS=False (config.py) ")
        sno_output_file = options.sno_file
        found_matchups = find_crosses.parse_crosses_file(sno_output_file)
        if len(found_matchups) == 0:
            logger.warning("No matchups found in SNO output file %s" %
                      sno_output_file)
            if options.debug is True:
                raise Warning()
        else:
            matchups.extend(found_matchups)
    elif options.pps_product_file is not None:
        from find_crosses import Cross
        from runutils import  parse_scenesfile_v2014
        pps_output_file = options.pps_product_file
        read_from_file = open(pps_output_file,'r')
        for line in read_from_file:      
            if line.rstrip() in "":
                pass
            else:   
                satname, time = parse_scenesfile_v2014(line)
                matchups.append(Cross(satname, '', time, time, -999, -999))
    elif options.cci_product_file is not None:
        from find_crosses import Cross
        from runutils import  parse_scenesfile_cci
        cci_output_file = options.cci_product_file
        read_from_file = open(cci_output_file,'r')
        for line in read_from_file:
            if line.rstrip() in "":
                pass
            else:    
                satname, time = parse_scenesfile_cci(line)
                matchups.append(Cross(satname, '', time, time, -999, -999))
    elif options.maia_product_file is not None:
        from find_crosses import Cross
        from runutils import  parse_scenesfile_maia
        maia_output_file = options.maia_product_file
        read_from_file = open(maia_output_file,'r')
        for line in read_from_file:
            if line.rstrip() in "":
                pass
            else:    
                satname, time = parse_scenesfile_maia(line)
                matchups.append(Cross(satname, '', time, time, -999, -999))
    elif options.reshaped_product_file is not None:
        from find_crosses import Cross
        from runutils import  parse_scenesfile_reshaped
        reshaped_output_file = options.reshaped_product_file
        read_from_file = open(reshaped_output_file,'r')
        for line in read_from_file:
            if line.rstrip() in "":
                pass
            else: 
                satname, time = parse_scenesfile_reshaped(line)
                #print time
                matchups.append(Cross(satname, '', time, time, -999, -999))

    process_matchups(matchups, run_modes, reprocess, options.debug)
    
    return 0

#------------------------------------------------------------------------------
if __name__=='__main__':
    main()
