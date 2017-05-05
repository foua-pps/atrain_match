"""
Calculate how well the tracks of two satellites correlate.

"""

import ephem
from merge_tles import get_tles

#: Degrees to radians
_d2r = 3.141597 / 180


#: Text file containing two-line elements (TLEs)
_TLE_FILE = 'TLINSET'

def set_tle_file(filename):
    global _TLE_FILE
    _TLE_FILE = filename

def get_tle_file():
    return _TLE_FILE


def get_tle(filename, satellite, time=None):
    """Get the two-line element for *satellite* from file *filename*, closest
    in time to *time*. If *time* is None, return the latest TLE.
    
    Returns a 2-tuple (line1, line2).
    
    """
    from datetime import datetime
    from merge_tles import tle_time
    
    if time is None:
        time = datetime(3000, 1, 1) # Far enough in the future?
    
    tles = []
    for tle in get_tles(filename):
        tle_time, name, l1, l2 = tle
        if name.lower() == satellite.lower():
            tles.append((abs(tle_time(l1) - time), l1, l2))
    
    tles.sort()
    time_diff, l1, l2 = tles[0] #@UnusedVariable
    return l1, l2


def satpos(satellite, time):
    """Get the (sub) position of *satellite* at *time*.
    
    Returns a tuple (lon, lat), in degrees.
    
    """
    sat = ephem.readtle(satellite, *get_tle(_TLE_FILE, satellite, time))
    sat.compute(time)
    
    return sat.sublong / _d2r, sat.sublat / _d2r

def angle(p1, p2, radians=False):
    """Calculate the angle in degrees between points *p1* and *p2*.
    
    If *radians* is False, *p1* and *p2* are assumed to be in degrees, and
    the returned angle will be in degrees. Otherwise, radians are assumed, and
    returned.
    
    """
    import numpy as np
    
    if radians:
        factor = 1.
    else:
        factor = _d2r
        
    lon1, lat1 = np.array(p1) * factor
    lon2, lat2 = np.array(p2) * factor
    
    from math import acos, sin, cos
    
    # From http://en.wikipedia.org/wiki/Great-circle_distance
    return acos(sin(lon1) * sin(lon2) +
                cos(lon1) * cos(lon2) * cos(lat1 - lat2)) / factor
    
    # WGS84
    """
    Nina 20170504 commented out unused code
    a = 6378137.0
    b = 6356752.3142    
    from greatcircle import GreatCircle
    circle = GreatCircle(a, b, lon1, lat1, lon2, lat2)    
    return circle.distance
    """

def satangle(satellite1, satellite2, time):
    """Get the earth central angle (in degrees) between *satellite1* and 
    *satellite2* at *time*.
    
    """
    p1 = satpos(satellite1, time)
    p2 = satpos(satellite2, time)
    
    return angle(p1, p2)
