"""
Utilities for running matching

"""

import os
import logging
logger = logging.getLogger(__name__)


def unzip_file(filename):
    """Unzip the file if file is bzipped = ending with 'bz2'"""

    import tempfile
    import bz2
    import gzip
    if filename.endswith('.bz2'):
        bz2file = bz2.BZ2File(filename)
        #tmpfilename = tempfile.mktemp()
        tmpdir = tempfile.mkdtemp()
        tmpfilename = os.path.join(tmpdir,
                                   os.path.basename(filename).strip('.bz2'))
        try:
            ofpt = open(tmpfilename, 'wb')
            ofpt.write(bz2file.read())
            ofpt.close()
        except IOError:
            import traceback
            traceback.print_exc()
            logger.info("Failed to read bzipped file %s", str(filename))
            os.remove(tmpfilename)
            return None

        return tmpfilename

    elif filename.endswith('.gz'):
        #tmpfilename = tempfile.mktemp()
        tmpdir = tempfile.mkdtemp()
        tmpfilename = os.path.join(tmpdir,
                                   os.path.basename(filename).strip('.gz'))
        with gzip.open(filename, 'rb') as f:
            try:
                ofpt = open(tmpfilename, 'wb')
                ofpt.write(f.read())
                ofpt.close()
            except IOError:
                import traceback
                traceback.print_exc()
                logger.info("Failed to read bzipped file %s", str(filename))
                os.remove(tmpfilename)
                return None

        return tmpfilename

    return None


def process_scenes(scenes, fun, options, ignore_errors=True, *args, **kwargs):
    """
    For each string in *scenes*, process corresponding noaa scene, by calling
    function *fun* with arguments satname, orbit, as parsed from *scenes*. Each
    element in *scenes* should look like
    'some/path/noaa19_20100110_1045_04767...'.

    If *ignore_errors* is True (default), exceptions will be ignored and
    problematic scenes reported after all scenes have been processed.

    Any additional arguments and keyword arguments are passed to *fun*.

    """
    errors = []
    for _file in scenes:
        print _file, fun
        filename = os.path.basename(_file)
        # print filename, fun
        #satname, _datetime, orbit = parse_scene(filename)
        logger.info("File : %s" % (filename))
        try:
            fun(_file, options, *args, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception, err:
            if ignore_errors:
                errors.append((filename, err))
            else:
                raise

    if errors:
        logger.info("The following errors were caught")
        for _file, err in errors:
            logger.error("Error processing scene %r: %r" %
                         (_file.strip('.okay'), err))


def parse_scene(filename):
    """
    Parse scene string (e.g. 'noaa19_20100110_1045_04767') and return
    (satname, `datetime.datetime`, orbit).

    Examples:

    >>> parse_scene('noaa19_20100110_1045_04767')
    ('noaa19', datetime.datetime(2010, 1, 10, 10, 45), 4767)
    >>> parse_scene('noaa19_20100110_1045_04767.okay')
    ('noaa19', datetime.datetime(2010, 1, 10, 10, 45), 4767)

    Leading directories are removed:

    >>> parse_scene('some/dir/noaa19_20100110_1045_04767')
    ('noaa19', datetime.datetime(2010, 1, 10, 10, 45), 4767)

    """
    from datetime import datetime
    import re
    filename = os.path.basename(filename)
    if not filename:
        raise ValueError("No file %r" % filename)
    match = re.match('(\w+)_([0-9]{8})_([0-9]{4})_([0-9]+)', filename)
    if not match:
        raise ValueError("Couldn't parse \"okay\" file %r" % filename)

    satname, date_s, time_s, orbit_s = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M')
    orbit = int(orbit_s)

    return satname, _datetime, orbit


def parse_scenesfile_v2014(filename):
    """
    Parse pps file =S_NWC_CT_{satellite}_{orbit}_%Y%m%dT%H%M???Z_*.h5 or .nc
    """
    from datetime import datetime
    import re
    filename = os.path.basename(filename)
    if not filename:
        raise ValueError("No file %r" % filename)
    match = re.match(r"S_NWC_.+_([^_]+)_\d+_(\d+)T(\d\d\d\d).+", filename)
    if not match:
        raise ValueError("Couldn't parse pps file %r" % filename)
    satname, date_s, time_s = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M')

    return satname, _datetime


def parse_scenesfile_cci(filename):
    """
    Parse cci file: 20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc
    """
    from datetime import datetime
    import re
    filename = os.path.basename(filename)
    if not filename:
        raise ValueError("No file %r" % filename)
    match = re.match(
        r"(\d\d\d\d\d\d\d\d)(\d\d\d\d).+AVHRRGAC-([^-]+)-", filename)
    if not match:
        raise ValueError("Couldn't parse cci file %r" % filename)
    date_s, time_s, satname = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M')

    return satname, _datetime


def parse_scenesfile_maia(filename):
    """
    Parse maia file:  #viiCT_npp_DB_20120817_S035411_E035535_DES_N_La052_Lo-027_00001.h5
    """
    from datetime import datetime
    import re
    filename = os.path.basename(filename)
    if not filename:
        raise ValueError("No file %r" % filename)
    match = re.match(r"\S\S\SCT_(\S\S\S)_\S\S_(\d+)_S(\d+)_", filename)
    if not match:
        raise ValueError("Couldn't parse maia file %r" % filename)
    satname, date_s, time_s = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M%S')

    return satname, _datetime

def parse_scenesfile_reshaped(filename):
    """
    Parse maia file:  #5km_noaa18_20090328_1855_99999_caliop_avhrr_match.h5
    """
    from datetime import datetime
    import re
    filename = os.path.basename(filename)
    if not filename:
        raise ValueError("No file %r" % filename)
    match = re.match(r"\d+km_([^_]+)_(\d+)_(\d+)_", filename)
    if not match:
        raise ValueError("Couldn't parse reshaped file %r" % filename)
    satname, date_s, time_s = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M')
    return satname, _datetime
