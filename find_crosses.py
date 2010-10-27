#!/usr/bin/python
"""Find satellite passes where the given satellites cross paths within a given
time window."""

import os

def get_environ(name, default=None):
    """Get the environment variable *name*. If it is not defined, return 
    *default*."""
    try:
        return os.environ[name]
    except KeyError:
        return default

# SNO executable
SNO_EXECUTABLE = get_environ('SNO_EXECUTABLE', "/data/proj/saf/ejohansson/SNO_tools/snotimes")
# Time window between satellite data reception time and sat1-sat2 crossing time, in seconds
SATELLITE_CROSS_RECEPTION_TIME_WINDOW = int(get_environ('SATELLITE_CROSS_RECEPTION_TIME_WINDOW', 60*60))
AAPP_WORKING_DIR = get_environ('AAPP_WORKING_DIR')


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


def find_crosses(satellite1, start, end, satellite2='calipso', t_window=2, lon_range=None,
                 lat_range=None):
    """Use snotimes to find satellite passes where the given *satellite* crosses 
    CALIPSO's path within a time window *t_window*. Sort out any cross points
    which are outside *lon_range* and/or *lat_range*."""
    import subprocess
    
    cmd = [SNO_EXECUTABLE, satellite1, satellite2, start, end, str(t_window)]
    print(' '.join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return parse_crosses_file(process.stdout, lon_range, lat_range)


def parse_crosses_file(f, lon_range=None, lat_range=None):
    """Parse output from snotimes in file *f*. Filter on *lon_range* and 
    *lat_range*, if provided."""
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
                satellite1 = a[0]
                satellite2 = a[1]
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


def parse_range(range):
    """Parse *range*, which should look like -35.7354:25.1. Return (lower, upper)."""
    l = range.split(':')
    return (float(l[0]), float(l[1]))


def prepare_aapp(hrpt_filename, aapp_working_dir=AAPP_WORKING_DIR):
    """Create symbolic link to *hrpt_filename* in *aapp_working_dir/incoming*,
    in preparation for running AAPP to create Level 1B AVHRR files, suitable for
    PPS processing."""
    import os
    os.symlink(hrpt_filename, os.path.join(aapp_working_dir, 'incoming', 
               os.path.basename(hrpt_filename)))


if __name__ == '__main__':
    from optparse import OptionParser
    
    # Set up and handle command line arguments
    parser = OptionParser()
    parser.set_usage("usage: %prog [options] satellite1 [satellite2 [...]]")
    parser.add_option('-t', '--time_window', type='int', default=2,
                      help="Time window for crossing.")
    parser.add_option('-x', '--longitudes', type='string',
                      help="Range of acceptable longitudes (e.g. -35.7:25.12).")
    parser.add_option('-y', '--latitudes', type='string',
                      help="Range of acceptable latitudes (e.g. -35.7:25.12).")
    
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
    
    for satellite in args:
        find_crosses(satellite, options.time_window, lon_range, lat_range)
