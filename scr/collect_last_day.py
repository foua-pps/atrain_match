"""
Collect data from /data/24 for satellite passes crossing paths with Calipso. The
collected data can then be used for further validation processing.
"""

from urllib2 import urlopen
import re
import logging
import os


#: URL to Cloudsat and Calipso two-line elements (TLEs)
SCIENCE_URL = 'http://celestrak.com/NORAD/elements/science.txt'

#: URL to NOAA two-line elements (TLEs)
NOAA_URL = 'http://celestrak.com/NORAD/elements/noaa.txt'

#: Base directory where data should be taken from
PPS_DATA24_DIR = '/data/24/saf/pps'

#: Base directory where data should be stored
PPS_ARKIV_DIR = os.getenv('PPS_ARKIV_DIR', '/data/arkiv/proj/safworks/data/pps')

#: Required PPS file types to collect, ordered by importance. If one file is missing, processing of the cross is stopped
REQUIRED_FILE_TYPES = ['lvl1b', 'NWP', 'satpos']
#: Reproduceable PPS file types to either collect, or reproduce
REPRODUCIBLE_FILE_TYPES = ['avhrr', 'nwp_tsur', 'cloudtype', 'ctth',
                            'ctth_opaque', 'ctth_semitransparent', 'sunsatangles']

PPS_FILE_TYPES = []
PPS_FILE_TYPES.extend(REQUIRED_FILE_TYPES)
PPS_FILE_TYPES.extend(REPRODUCIBLE_FILE_TYPES)

#: Default time window, in minutes
TIME_WINDOW = 30

logging.basicConfig()
logger = logging.getLogger('collect')


def run_pps(cross, file_type=None):
    """
    Produce cloud mask, cloud type and CTTH products for the given *cross*.
    """
    from file_finders.pps import Lvl1bDirFinder
    lvl1b_finder = Lvl1bDirFinder(PPS_ARKIV_DIR)
    found_dir = lvl1b_finder.find(cross.time1, cross.satellite1)[0]
    parsed = lvl1b_finder.parse(found_dir)
    orbit = parsed['orbit']
    datetime = parsed['datetime']
    satday = "%d%02d%02d" % (datetime.year, datetime.month, datetime.day)
    sathour = "%02d%02d" % (datetime.hour, datetime.minute)
    
    try:
        from pps_basic_configure import NWP_TIME_RESOLUTION
        
        try:
            logger.debug("Running create_pps_argument")
            import pps_arguments
            ppsarg = pps_arguments.create_pps_argument(pps_arguments.USE_AAPP_LVL1B_GLOB,
                                                       pps_arguments.SATELLITE_PROJ,
                                                       cross.satellite1, orbit)
            import pps_ioproxy
            ioproxy = pps_ioproxy.pps_ioproxy(ppsarg)
        except:
            logger.warning("Couldn't build pps arguments")
            raise
        
        if file_type in ['avhrr', 'sunsatangles', None]:
            logger.debug("file_type = '%s': Making avhrr and sunsatangles" % file_type)
            try:
                logger.debug("Running extractAvhrrFromLvl1b")
                retv = ioproxy.extractAvhrrFromLvl1b()
                if retv is not None:
                    ioproxy.writeAvhrrData(retv)
                    ioproxy.writeSunSatAngles(retv)
            except:
                logger.warning("Couldn't make avhrr and sunsatangles")
                raise
        
        elif 'nwp' in file_type or file_type is None:
            logger.debug("file_type = '%s': Making remapped nwp" % file_type)
            try:
                # Physiography is needed for NWP remapping
                logger.debug("Running extractPhys")
                import ppsMakePhysiography
                status = ppsMakePhysiography.extractPhys(ppsarg)
                if status != 0:
                    raise RuntimeError("Status returned from extractPhys: %d" % status)
                
                logger.debug("Running extractNwp")
                import ppsMakeNwp
                status = ppsMakeNwp.extractNwp(ppsarg, NWP_TIME_RESOLUTION)
                if status != 0:
                    raise RuntimeError("Status returned from extractNwp: %d" % status)
            except:
                logger.warning("Couldn't make remapped nwp")
                raise
        
        elif file_type in ['cloudtype', None]:
            logger.debug("file_type = '%s': Making cloudtype" % file_type)
            try:
                logger.debug("Running extractEmiss")
                import ppsMakeEmissivity
                status = ppsMakeEmissivity.extractEmiss(ppsarg)
                if status != 0:
                    raise RuntimeError("Status returned from extractEmiss: %d" % status)
                
                logger.debug("Running pps_prepare_pge0102")
                import CloudMaskCloudTypePrepare
                status = CloudMaskCloudTypePrepare.pps_prepare_pge0102(ppsarg)
                if status != 0:
                    raise RuntimeError("Status returned from pps_prepare_pge0102: %d" % status)
                
                logger.debug("Running pps_pge01")
                import CloudMask
                status = CloudMask.pps_pge01(ppsarg,None,1)
                if status != 0:
                    raise RuntimeError("Status returned from pps_pge01: %d" % status)
                
                logger.debug("Running pps_pge02")
                import ppsCtype
                status = ppsCtype.pps_pge02(satday, sathour, orbit, cross.satellite1, None)
                if status != 0:
                    raise RuntimeError("Status returned from pps_pge02: %d" % status)
            except:
                logger.warning("Couldn't make cloudtype")
                raise
        
        elif 'ctth' in file_type or file_type is None:
            logger.debug("file_type = '%s': Making ctth" % file_type)
            try:
                logger.debug("Running pps_prepare_ctth")
                import ppsCtthPrepare
                status = ppsCtthPrepare.pps_prepare_ctth(ppsarg, NWP_TIME_RESOLUTION)
                if status != 0:
                    raise RuntimeError("Status returned from pps_prepare_ctth: %d" % status)
                
                logger.debug("Running pps_pge03")
                import ppsCtth
                status = ppsCtth.pps_pge03(ppsarg)
                if status != 0:
                    raise RuntimeError("Status returned from pps_pge03: %d" % status)
                
                logger.debug("Running pps_pge03_semi")
                import ppsCtthPlus
                status = ppsCtthPlus.pps_pge03_semi(ppsarg)
                if status != 0:
                    raise RuntimeError("Status returned from pps_pge03_semi: %d" % status)
            except:
                logger.warning("Couldn't make ctth")
                raise
        
        else:
            raise NotImplementedError("file_type = %s not yet supported in run_pps" % file_type)
    except:
        logger.info("Could not reproduce PPS products.")
        raise


def get_tles():
    """Get latest TLEs from NORAD, and store them in a file ``TLINSET``."""
    from find_crosses import SNO_EXECUTABLE
    
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


def NWP_subdir(self, time, *args, **kwargs):
    """NWP data is stored in a non-standard subdir."""
    return "import/NWP_data/global_out"


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
    from glob import glob
    
    start = (date.today() - timedelta(1)).strftime('%y%m%d')
    end = (date.today() + timedelta(1)).strftime('%y%m%d')
    
    crosses = set()
    for satellite in satellites:
        crosses.update(find_crosses(satellite, start, end, 'CALIPSO', time_window))
    
    pps_finder = PpsExtendedFileFinder(data24dir)
    pps_finder.get_finder('NWP').set_subdir_method(NWP_subdir)
    
    logger.info("Found %d crosses of satellites %s with Calipso" % (len(crosses), ', '.join(satellites)))
    for cross in sorted(crosses):
        logger.debug("Copying files for %s" % str(cross))
        for file_type in file_types:
            try:
                if file_type == 'NWP':
                    time_window = 3*60*60 # 3 h
                else:
                    time_window = None # Use default time_window
                found_file = pps_finder.find(cross.time1, file_type=file_type,
                                                satname=cross.satellite1,
                                                time_window=time_window)[0]
            except IndexError:
                logger.debug("No %s file found in source path" % file_type)
                logger.debug("Search pattern was '%s'" % pps_finder.pattern(cross.time1, file_type=file_type,
                                                                            satname=cross.satellite1))
                if file_type in REQUIRED_FILE_TYPES:
                    logger.debug("Stopping further processing of %s" % cross)
                    break # Don't continue with the rest of the file_types
                elif file_type in REPRODUCIBLE_FILE_TYPES:
                    if len(pps_finder.find(cross.time1, file_type=file_type,
                                              satname=cross.satellite1,
                                              basedir=arkivdir)) == 0:
                        logger.debug("No %s file found at destination" % file_type)
                        logger.debug("Search pattern was '%s'" % pps_finder.pattern(cross.time1, file_type=file_type,
                                                                                       satname=cross.satellite1,
                                                                                       basedir=arkivdir))
                        logger.info("Regenerating %s data" % file_type)
                        try:
                            run_pps(cross, file_type)
                        except:
                            continue
                    else:
                        logger.debug("%s file already exists at destination" % file_type)
            
            dst_dir = os.path.join(arkivdir, pps_finder.subdir(cross.time1, file_type=file_type,
                                                           satname=cross.satellite1))
            dst_path = os.path.join(dst_dir, os.path.basename(found_file))
            
            # Hack to change NWP subdir to standard
            if file_type == 'NWP':
                dst_path = dst_path.replace(pps_finder.get_finder('NWP').subdir(None),
                                            pps_finder.get_finder('NWP').__class__.subdir(pps_finder.get_finder('NWP'), None))
            
            if file_type == 'lvl1b':
                # We need to copy the lvl1b directory + '.okay' files
                if not os.path.exists(dst_path):
                    try:
                        shutil.copytree(found_file, dst_path, symlinks=False)
                    except IOError:
                        logger.critical("Couldn't copy directory %s to %s." % (found_file, dst_path))
                        raise
                else:
                    logger.debug("Destination already exists: %s", dst_path)
                for f in glob(found_file + '*.okay'):
                    dst_path = os.path.join(dst_dir, os.path.basename(f))
                    if not os.path.exists(dst_path):
                        try:
                            shutil.copy(f, dst_path)
                        except IOError:
                            logger.critical("Couldn't copy %s to %s." % (f, dst_path))
                            raise
                    else:
                        logger.debug("Destination already exists: %s", dst_path)
            else:
                if not os.path.exists(dst_path):
                    try:
                        shutil.copy(found_file, dst_path)
                    except IOError:
                        logger.critical("Couldn't copy %s to %s." % (found_file, dst_path))
                        raise
                else:
                    logger.debug("Destination already exists: %s", dst_path)


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
