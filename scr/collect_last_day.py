"""
Collect data from /data/24 for satellite passes crossing paths with Calipso. The
collected data can then be used for further validation processing.
"""

from urllib2 import urlopen
import re
import logging


#: URL to Cloudsat and Calipso two-line elements (TLEs)
SCIENCE_URL = 'http://celestrak.com/NORAD/elements/science.txt'

#: URL to NOAA two-line elements (TLEs)
NOAA_URL = 'http://celestrak.com/NORAD/elements/noaa.txt'

#: Base directory where data should be taken from
PPS_DATA24_DIR = '/data/24/saf/pps'

#: Base directory where data should be stored
PPS_ARKIV_DIR = '/data/arkiv/proj/safworks/data/pps'

#: PPS file types to collect, ordered by importance. If one file is missing, all following files will be skipped
PPS_FILE_TYPES = ['lvl1b', 'NWP', 'avhrr', 'nwp_tsur', 'cloudtype',
                  'ctth', 'ctth_opaque', 'ctth_semitransparent', 'sunsatangles']

#: Default time window, in minutes
TIME_WINDOW = 30

logging.basicConfig()
logger = logging.getLogger('collect')


def get_tles():
    """Get latest TLEs from NORAD, and store them in a file ``TLINSET``."""
    from find_crosses import SNO_EXECUTABLE
    import os
    
    try:
        science = urlopen(SCIENCE_URL)
        noaa = urlopen(NOAA_URL)
    except:
        logger.critical("Could not get TLEs.")
        raise
    
    try:
        out = open(os.path.join(os.path.dirname(SNO_EXECUTABLE), 'TLINSET'), 'w')
        for l in science.readlines():
            l_unix = l.replace('\r', '')
            out.write(l_unix)
        
        for l in noaa.readlines():
            l_unix = l.replace('\r', '')
            l_names_corrected = re.sub('NOAA ([0-9]+).*', 'NOAA\\1', l_unix)
            out.write(l_names_corrected)
        out.close()
    except IOError:
        logger.critical("Could not write to 'TLINSET'.")
        raise


def get_files(satellites=['noaa18', 'noaa19'], time_window=TIME_WINDOW,
              data24dir=PPS_DATA24_DIR, arkivdir=PPS_ARKIV_DIR, file_types=PPS_FILE_TYPES):
    """
    Find crosses and copy all files needed to run validation for found crosses.
    
    *satellites* should be a list of satellite names to consider. (Default: noaa18
    and noaa19. (metop???))
    *time_window* is the time in minutes to pass to SNO ``findtimes``.
    *basedir* is the base directory for PPS data.
    """
    from find_crosses import find_crosses
    from datetime import date, timedelta
    from file_finders import PpsExtendedFileFinder
    import shutil
    import os
    from glob import glob
    
    start = (date.today() - timedelta(1)).strftime('%y%m%d')
    end = (date.today() + timedelta(1)).strftime('%y%m%d')
    
    crosses = set()
    for satellite in satellites:
        crosses.update(find_crosses(satellite, start, end, 'CALIPSO', time_window))
    
    pps_finder = PpsExtendedFileFinder(data24dir)
    for cross in sorted(crosses):
        for file_type in file_types:
            try:
                found_file = pps_finder.find(cross.time1, file_type=file_type,
                                             satname=cross.satellite1)[0]
            except IndexError:
                logger.debug("No %s file found for %s at %s." % (file_type, cross, pps_finder.pattern(cross.time1, file_type=file_type, satname=cross.satellite1)))
                continue
            
            dst = os.path.join(arkivdir, pps_finder.subdir(cross.time1, file_type=file_type,
                                                           satname=cross.satellite1))
            if file_type == 'lvl1b':
                # We need to copy the lvl1b directory + '.okay' files
                dst_lvl1b = os.path.join(dst, os.path.basename(found_file))
                try:
                    shutil.copytree(found_file, dst_lvl1b, symlinks=False)
                except IOError:
                    logger.critical("Couldn't copy directory %s to %s." % (found_file, dst_lvl1b))
                    raise
                for f in glob(found_file + '*.okay'):
                    try:
                        shutil.copy(f, dst)
                    except IOError:
                        logger.critical("Couldn't copy %s to %s." % (f, dst))
                        raise
            else:
                try:
                    shutil.copy(found_file, dst)
                except IOError:
                    logger.critical("Couldn't copy %s to %s." % (found_file, dst))
                    raise


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-d', '--debug', action='store_true', default=False)
    
    (options, arguments) = parser.parse_args()
    if options.debug is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    get_tles()
    get_files()