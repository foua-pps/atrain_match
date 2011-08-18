#!/usr/bin/python
"""
Find satellite passes where the given satellites cross paths within a given
time window.

.. note::

    This module has a command line interface.

"""

from config import get_environ

# SNO executable
SNO_EXECUTABLE = get_environ('SNO_EXECUTABLE', "/data/proj/saf/ejohansson/SNO_tools/snotimes")
# Time window between satellite data reception time and sat1-sat2 crossing time, in seconds
SATELLITE_CROSS_RECEPTION_TIME_WINDOW = int(get_environ('SATELLITE_CROSS_RECEPTION_TIME_WINDOW', 60*60))


class Cross:
    """A cross where two satellites (almost) meet."""
    def __init__(self, satellite1, satellite2, time1, time2, lon, lat, time_window=None):
        self.satellite1 = satellite1
        self.satellite2 = satellite2
        self.time1 = time1
        self.time2 = time2
        self.lon = lon
        self.lat = lat
        self.time_window = time_window
    
    def as_sno_line(self):
        """Return a string representation of the matchup, such as it would be
        output from SNO."""
        sno_line = ''
        for time in [self.time1, self.time2]:
            sno_line += "  %d %2d %2d %2d %2d %5.1f   " % (time.year, time.month, 
                     time.day, time.hour, time.minute, time.second)
        sno_line += "   % 7.2f  % 7.2f" % (self.lat, self.lon)
        return sno_line
    
    def __repr__(self):
        if self.time1 > self.time2:
            string_formatting_args = (self.satellite2, self.satellite1, self.time2,
                                      self.time1 - self.time2, self.lon, self.lat)
        else:
            string_formatting_args = (self.satellite1, self.satellite2, self.time1,
                                      self.time2 - self.time1, self.lon, self.lat)
        return "Cross(%s x %s: %s (+ %s), (%6.1f,%6.1f))" % string_formatting_args
    
    def __str__(self):
        return self.__repr__()
    
    def __lt__(self, other):
        return self.time1 < other.time1
    
    def __eq__(self, other):
        return self.time1 == other.time1 and self.time2 == other.time2 and \
               self.satellite1 == other.satellite1 and self.satellite2 == other.satellite2
    
    def __hash__(self):
        return hash(self.time1) ^ hash(self.time2) ^ hash(self.satellite1) ^ hash(self.satellite2)



def parse_sno_matchup_line(l):
    """Parse output line *l* from snotimes, and return a tuple with 
    (*time1*, *time2*, *lon*, *lat*)."""
    from datetime import datetime
    
    a = l.split()
    time1 = datetime(int(a[0]), int(a[1]), int(a[2]),
                     int(a[3]), int(a[4]), int(float(a[5])))
    time2 = datetime(int(a[6]), int(a[7]), int(a[8]),
                     int(a[9]), int(a[10]), int(float(a[11])))
    lat = float(a[12])
    lon = float(a[13])
    return (time1, time2, lon, lat)


def find_crosses(satellite1, start, end, satellite2='calipso', time_window=20, lon_range=None,
                 lat_range=None):
    """Use snotimes to find satellite passes where the given *satellite* crosses 
    CALIPSO's path within a time window *t_window*. Sort out any cross points
    which are outside *lon_range* and/or *lat_range*."""
    from subprocess import Popen, PIPE
    import os
    
    wd = os.getcwd()
    os.chdir(os.path.dirname(SNO_EXECUTABLE))
    cmd = [os.path.abspath(os.path.basename(SNO_EXECUTABLE)), satellite1, satellite2, start, end, str(time_window)]
    print(' '.join(cmd))
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    process.stderr.close()
    os.chdir(wd)
    return parse_crosses_file(process.stdout, lon_range, lat_range)


def parse_crosses_file(f, lon_range=None, lat_range=None):
    """Parse output from snotimes in file *f*. Filter on *lon_range* and 
    *lat_range*, if provided.
    
    Examples
    --------
    
    >>> crosses = parse_crosses_file(filename) # doctest: +SKIP
    >>> print(crosses[0]) # doctest: +SKIP
    Cross(aqua x envisat: 2009-01-28 19:35:49 (+ 0:09:03), (-117.1,  72.9))
    
    """
    satellite1 = None
    satellite2 = None
    time_window = None
    crosses = []
    if type(f) == str:
        f = open(f, 'r')
        f_needs_closing = True
    for l in f.readlines():
        try:
            # Assume normal matchup line from SNO
            (time1, time2, lon, lat) = parse_sno_matchup_line(l)
        except:
            import math
            try:
                # Catch line describing satellite names and time window
                a = l.split("Satellites and Time Window:")[1].split()
                satellite1 = a[0].replace('METOPA', 'METOP02').lower()
                satellite2 = a[1].replace('METOPA', 'METOP02').lower()
                time_window = int(math.ceil(float(a[2]))) * 60 # Time window in integer seconds
                continue
            except:
                continue
        
        if lon_range is not None:
            if lon_range[0] > lon or lon > lon_range[1]:
                continue
        if lat_range is not None:
            if lat_range[0] > lat or lat > lat_range[1]:
                continue
        crosses.append(Cross(satellite1, satellite2, time1, time2, lon, lat, time_window))
    
    try:
        if f_needs_closing:
            f.close()
    except:
        pass
    return crosses


def _daytime(cross):
    """
    Returns True if *cross* is during day, False otherwise.
    
    Examples
    --------
    
    >>> from datetime import datetime
    >>> t = datetime(2010, 6, 21, 12, 0, 0) # noon
    >>> cross = Cross('s1', 's2', t, t, 0, 0) # at equator
    >>> _daytime(cross) # => daytime
    True
    
    >>> cross.time1 = datetime(2010, 1, 1, 0, 0, 0) # midnight at equator
    >>> _daytime(cross) # => night time
    False
    
    >>> cross.lat = -89.9 # midnight at south pole, in antarctic summer
    >>> _daytime(cross) # => daytime
    True
    
    >>> cross.lat = 89.9 # midnight at north pole, in arctic winter
    >>> _daytime(cross) # => night time
    False
    
    """
    import ephem
    d2r = 3.141596 / 180 # Degrees to radians
    
    obs = ephem.Observer()
    obs.lon, obs.lat = cross.lon * d2r, cross.lat * d2r
    obs.date = cross.time1
    
    sun = ephem.Sun() #@UndefinedVariable
    try:
        day = obs.next_setting(sun) < obs.next_rising(sun)
    except ephem.NeverUpError:
        day = False
    except ephem.AlwaysUpError:
        day = True
    
    return day


def parse_range(range):
    """Parse *range*, which should look like -35.7354:25.1. Return (lower, upper)."""
    l = range.split(':')
    return (float(l[0]), float(l[1]))


if __name__ == '__main__':
    from optparse import OptionParser
    
    # Set up and handle command line arguments
    parser = OptionParser()
    parser.set_usage("usage: %prog [options] <start YYMMDD> <end YYMMDD> "
                     "satellite1 [satellite2]\n"
                     "Find times and locations where satellite1 and satellite2 "
                     "(or calipso, default) cross paths.")
    parser.add_option('-t', '--time_window', type='int', default=20,
                      help="Time window for crossing. Default: 20 min.")
    parser.add_option('-x', '--longitudes', type='string',
                      help="Range of acceptable longitudes (e.g. -35.7:25.12).")
    parser.add_option('-y', '--latitudes', type='string',
                      help="Range of acceptable latitudes (e.g. -35.7:25.12).")
    parser.add_option('-d', '--daytime', action='store_true',
                      help="Only daytime crosses")
    
    (options, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("at least one satellite must be specified")
    
    if options.longitudes is not None:
        try:
            lon_range = parse_range(options.longitudes)
        except ValueError:
            parser.error("longitude range should be specified as lower:upper, e.g. -35.7:25.12")
    else:
        lon_range = None
    
    if options.latitudes is not None:
        try:
            lat_range = parse_range(options.latitudes)
        except ValueError:
            parser.error("latitude range should be specified as lower:upper, e.g. -35.7:25.12")
    else:
        lat_range = None
    
    start = args[0]
    end = args[1]
    satellite1 = args[2]
    try:
        satellite2 = args[3]
    except IndexError:
        satellite2 = 'calipso'
    
    crosses = set()
    crosses.update(find_crosses(satellite1, start, end, satellite2=satellite2,
                                time_window=options.time_window,
                                lon_range=lon_range, lat_range=lat_range))
    
    for cross in sorted(crosses):
        if not options.daytime or _daytime(cross):
            print(cross)
