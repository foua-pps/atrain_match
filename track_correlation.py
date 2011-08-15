"""
Calculate how well the tracks of two satellites correlate.

"""

import ephem

#: Degrees to radians
_d2r = 3.141597 / 180


#: Text file containing two-line elements (TLEs)
_TLE_FILE = 'TLINSET'

def set_tle_file(filename):
    global _TLE_FILE
    _TLE_FILE = filename

def get_tle_file():
    return _TLE_FILE


def tle_time(l1):
    """Given first line in a two-line element, TLE, return a `datetime.datetime`
    object for when the TLE was valid.
    
    """
    from datetime import datetime, timedelta
    from math import floor
    
    epoch = l1.split()[3]
    year = int(epoch[:2])
    if year >= 57: # See https://www.space-track.org/tle_format.html#epoch2
        year += 1900
    else:
        year += 2000
    fractional_day = float(epoch[2:])
    day = floor(fractional_day)
    seconds = (fractional_day - day) * 24 * 60 * 60
    
    return datetime(year, 1, 1) + timedelta(floor(day), seconds)

def get_tle(filename, satellite, time=None):
    """Get the two-line element for *satellite* from file *filename*, closest
    in time to *time*. If *time* is None, return the latest TLE.
    
    Returns a 2-tuple (line1, line2).
    
    """
    from datetime import datetime
    
    if time is None:
        time = datetime(3000, 1, 1) # Far enough in the future?
    
    with open(filename, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    
    tles = []
    for ix in range(0, len(lines), 3):
        name, l1, l2 = lines[ix:ix + 3]
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
    
    return acos(sin(lon1) * sin(lon2) +
                cos(lon1) * cos(lon2) * cos(lat1 - lat2)) / factor
    
    # WGS84
    a = 6378137.0
    b = 6356752.3142
    
    from greatcircle import GreatCircle
    circle = GreatCircle(a, b, lon1, lat1, lon2, lat2)
    
    return circle.distance

def satangle(satellite1, satellite2, time):
    """Get the earth central angle (in degrees) between *satellite1* and 
    *satellite2* at *time*.
    
    """
    p1 = satpos(satellite1, time)
    p2 = satpos(satellite2, time)
    
    return angle(p1, p2)
