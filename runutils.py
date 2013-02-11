"""
Utilities for running matching

"""

import os
import logging
logger = logging.getLogger(__name__)


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
        print filename, fun
        satname, _datetime, orbit = parse_scene(filename)
        try:
            fun(satname, orbit, options, *args, **kwargs)
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
