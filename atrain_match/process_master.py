#!/usr/bin/python
"""
This module is the main entry point for matching noaa and/or metop data with
Cloudsat and Calipso data, and produce statistics for validation of PPS
cloud mask, cloud type and cloud top temperature and height products.

It is really a wrapper to :func:`truth_imager_match.run`, for running
through a set of files with SNO matchups.

.. note::

    This module has a command line interface.

"""
import logging
logging.basicConfig(
    format='%(levelname)s |%(asctime)s|: %(message)s',
    level=logging.INFO,
    #datefmt='%Y-%m-%d %H:%M:%S')
    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)
from utils.common import (InputError, MatchupError)
import config
from libs import truth_imager_match
from utils.common import Cross
from utils.runutils import  parse_scenesfile_reshaped
from utils.runutils import  parse_scenesfile_maia
from utils.runutils import parse_scene
from utils.runutils import  parse_scenesfile_cci
from utils.runutils import  parse_scenesfile_v2014


def process_matchups(matchups, run_modes, reprocess=False, debug=False):
    """
    Run the given *matchups* through the validation system.
    
    *matchups* should be a list of :class:`common.Cross` instances.
    
    If *reprocess* is True, disregard any previously generated matchup files.
    
    """

    from utils.runutils import read_config_info
    AM_PATHS, SETTINGS = read_config_info()

    problematic = set()
    no_matchup_files = []
    outstatus = 0
    for match in matchups:
        try:
            truth_imager_match.run(match, run_modes,  AM_PATHS, SETTINGS, reprocess)
        except MatchupError as err:
            logger.warning("Matchup problem: %s", str(err))
            import traceback
            traceback.print_exc()
            no_matchup_files.append(match)
            outstatus = 5
        except:
            import traceback
            outstatus = 5
            traceback.print_exc()
            problematic.add(match)
            logger.warning("Couldn't run truth_imager_match.")
            if debug is True:
                raise
                
    if len(no_matchup_files) > 0:
        logger.warning(
            "%d of %d cases had no matchups in region, within the time window:\n%s",
            len(no_matchup_files), len(matchups),
            '\n'.join([str(m) for m in no_matchup_files]))
    if len(problematic) > 0:
        logger.warning("%d of %d cases had unknown problems:\n%s",
                       len(problematic), len(matchups),
                       '\n'.join([str(m) for m in problematic]))
    return outstatus

def main():
    """
    Process command line options and run matchup and validation.
    
    If *args* is provided, it should be a list of command line arguments (e.g.
    sys.argv[1:]).
    
    For a complete usage description, run 'python process_master -h'.
    
    """
    import argparse
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('--mode', '-M', type=str, required=False, choices=config.ALLOWED_MODES,
                      help=("Run validation software in MODE "))
    parser.add_argument('--reprocess', '-r', const=True, nargs='?', required=False,
                        help="Disregard any previously generated Cloudsat- and "
                        "Calipso-IMAGER matchup files.")
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
                      help="Depricated")
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
        logging.getLogger().setLevel(logging.DEBUG)

    matchups = []
    if options.pps_okay_scene:
        # Simulate crosses from PPS scenes
        scene = options.pps_okay_scene
        satname, time, orbit = parse_scene(scene) 
        matchups.append(Cross(satname, time))
    elif options.pps_product_file is not None:
        pps_output_file = options.pps_product_file
        read_from_file = open(pps_output_file,'r')
        for line in read_from_file:      
            if line.rstrip() in "":
                pass
            else:   
                satname, time = parse_scenesfile_v2014(line)
                matchups.append(Cross(satname, time))
    elif options.cci_product_file is not None:
        cci_output_file = options.cci_product_file
        read_from_file = open(cci_output_file,'r')
        for line in read_from_file:
            if line.rstrip() in "":
                pass
            else:    
                satname, time = parse_scenesfile_cci(line)
                matchups.append(Cross(satname, time))
    elif options.maia_product_file is not None:
        maia_output_file = options.maia_product_file
        read_from_file = open(maia_output_file,'r')
        for line in read_from_file:
            if line.rstrip() in "":
                pass
            else:    
                satname, time = parse_scenesfile_maia(line)
                matchups.append(Cross(satname, time))
    elif options.reshaped_product_file is not None:
        reshaped_output_file = options.reshaped_product_file
        read_from_file = open(reshaped_output_file,'r')
        for line in read_from_file:
            if line.rstrip() in "":
                pass
            else: 
                satname, time = parse_scenesfile_reshaped(line)
                #print time
                matchups.append(Cross(satname, time))

    process_matchups(matchups, run_modes, reprocess, options.debug)
    
    return 0

#------------------------------------------------------------------------------
if __name__=='__main__':
    main()
