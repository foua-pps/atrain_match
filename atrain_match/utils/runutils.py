# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
"""
Utilities for running matching

"""

import re
import time
import numpy as np
import os
import logging
logger = logging.getLogger(__name__)


def read_config_info():
    import os
    from configparser import ConfigParser
    CONF = ConfigParser()
    from atrain_match.config import ATRAIN_MATCH_CONFIG_PATH
    from atrain_match.config import ATRAIN_MATCH_CONFIG_FILE
    config_file = os.path.join(ATRAIN_MATCH_CONFIG_PATH, ATRAIN_MATCH_CONFIG_FILE)
    if not os.path.isfile(config_file):
        raise IOError("Couldn't find config file %s." % (config_file))
    CONF.read(config_file)
    AM_PATHS = {}
    for option, value in CONF.items('files', raw=True):
        AM_PATHS[option] = value
    SETTINGS = {}
    for name, value in CONF.items('general', raw=True):
        name = name.upper()
        value = value.strip()
        while ' ' in value:
            value = value.replace(' ', '')
        values = value.split(',')
        if name in ['MIN_OPTICAL_DEPTH']:
            value_ = [np.float(val_i) for val_i in values]
        elif name in ["COMPILE_STATISTICS_TRUTH", "PLOT_MODES",
                      "PLOT_TYPES", "CTTH_TYPES",
                      'SATELLITES', 'YEARS', 'MONTHS']:
            value_ = values
        elif name in ['CNN_PCKL_PATH']:
            value_ = values[0]

        elif len(values) == 1 and 'true' in values[0].lower():
            value_ = True
        elif len(values) == 1 and 'false' in values[0].lower():
            value_ = False
        elif len(values) == 1 and re.match(r"\d+.*\d*", values[0]):
            value_ = np.float(values[0])

        SETTINGS[name.upper()] = value_

    if (SETTINGS['COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC'] or
            SETTINGS['CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA']):
        logger.info("Setting ALSO_USE_5KM_FILES = True as "
                    "COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC = True or "
                    "CALCULATE_DETECTION_HEIGHT_FROM_5KM_DATA = True")
        SETTINGS['ALSO_USE_5KM_FILES'] = True  # 5km data is required also for 1km processing
    SETTINGS['sec_timeThr'] = SETTINGS['MINUTES_TIMETHR']*60.0
    SETTINGS['sec_timeThr_synop'] = SETTINGS['MINUTES_TIMETHR_SYNOP']*60.0
    SETTINGS['SAT_ORBIT_DURATION'] = SETTINGS['SAT_ORBIT_DURATION_MINUTES']*60.0
    return AM_PATHS, SETTINGS


def unzip_file(filename):
    """Unzip the file if file is bzipped = ending with 'bz2'"""

    import tempfile
    import bz2
    import gzip
    if filename.endswith('.bz2'):
        bz2file = bz2.BZ2File(filename)
        # tmpfilename = tempfile.mktemp()
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
        # tmpfilename = tempfile.mktemp()
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
    match = re.match(r'(\w+)_([0-9]{8})_([0-9]{4})_([0-9]+)', filename)
    if not match:
        raise ValueError("Couldn't parse \"okay\" file %r" % filename)

    satname, date_s, time_s, orbit_s = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M')
    orbit = int(orbit_s)

    return satname, _datetime, orbit


def parse_scenesfile_v2014(filename):
    """
    Parse pps file =S_NWC_CT_{satellite}_{orbit}_%Y%m%dT%H%M%S?Z_*.h5 or .nc
    """
    from datetime import datetime
    import re
    filename = os.path.basename(filename)
    if not filename:
        raise ValueError("No file %r" % filename)
    match = re.match(r"S_NWC_.+_([^_]+)_\d+_(\d+)T(\d\d\d\d\d\d).+", filename)
    if not match:
        raise ValueError("Couldn't parse pps file %r" % filename)
    satname, date_s, time_s = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M%S')

    return satname, _datetime


def parse_scenesfile_cci(filename):
    """
    Parse cci file: 20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-IMAGERGAC-NOAA18-fv1.0.nc
    OR
    Parse cci file: 20190713002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-SEVIRI-MSG4-fv1.0.nc
    """
    from datetime import datetime
    import re
    filename = os.path.basename(filename)
    if not filename:
        raise ValueError("No file %r" % filename)

    # CCI data
    if "IMAGERGAC" in filename:
        match = re.match(
            r"(\d\d\d\d\d\d\d\d)(\d\d\d\d).+IMAGERGAC-([^-]+)-", filename)
    # CCI+ data
    elif "SEVIRI" in filename:
        match = re.match(
            r"(\d\d\d\d\d\d\d\d)(\d\d\d\d).+SEVIRI-([^-]+)-", filename)
    else:
        raise ValueError("atrain_match not able to handle %r files" % filename)

    if not match:
        raise ValueError("Couldn't parse cci file %r" % filename)

    date_s, time_s, satname = match.groups()
    _datetime = datetime.strptime(date_s + time_s, '%Y%m%d%H%M')

    return satname.lower(), _datetime


def parse_scenesfile_maia(filename):
    """
    Parse maia file:  # viiCT_npp_DB_20120817_S035411_E035535_DES_N_La052_Lo-027_00001.h5
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
    Parse maia file:  # 5km_noaa18_20090328_1855_99999_caliop_imager_match.h5
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


def do_some_geo_obj_logging(GeoObj):
    import time
    tim1 = time.strftime("%Y%m%d %H:%M",
                         time.gmtime(GeoObj.sec1970_start))
    tim2 = time.strftime("%Y%m%d %H:%M",
                         time.gmtime(GeoObj.sec1970_end))
    logger.debug("Starttime: %s, end time: %s", tim1, tim2)
    logger.debug("Min lon: %f, max lon: %d",
                 np.min(np.where(
                     np.equal(GeoObj.longitude, GeoObj.nodata),
                     99999,
                     GeoObj.longitude)),
                 np.max(GeoObj.longitude))
    logger.debug("Min lat: %d, max lat: %d",
                 np.min(np.where(
                     np.equal(GeoObj.latitude, GeoObj.nodata),
                     99999,
                     GeoObj.latitude)),
                 np.max(GeoObj.latitude))


def do_some_logging(retv, match_obj):
    logger.debug("Start and end times: %s %s",
                 time.gmtime(match_obj.sec_1970[0]),
                 time.gmtime(match_obj.sec_1970[-1]))
    logger.debug("Maximum and minimum time differences in sec (imager-reference): %d %d",
                 np.max(retv.diff_sec_1970), np.min(retv.diff_sec_1970))
    logger.debug("IMAGER observation time of first imager-reference match: %s",
                 time.gmtime(retv.imager.sec_1970[0]))
    logger.debug("IMAGER observation time of last imager-reference match: %s",
                 time.gmtime(retv.imager.sec_1970[-1]))


if __name__ == "__main__":
    pass
