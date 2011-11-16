"""
Fetch NWP files from ecaccess.ecmwf.int

"""

import ftplib
import os
import sys

#: ftp server to fetch NWP from
SERVER = 'ecaccess.ecmwf.int'

#: username for ftp login, taken from $ECMWF_USER
USER = os.environ.get('ECMWF_USER', 'su3')

#: server directory from which to fetch NWP, taken from $ECMWF_DIR
DIRECTORY = os.environ.get('ECMWF_DIR', 'ECSCRATCH/nwp')


class ParseError(ValueError):
    """Raised when parsing failed."""

def parse_nwp(filename):
    """
    Parse *filename* and return a `datetime.datetime` object of the
    forecast/analysis time.
    
    @param filename: filename to parse
    @type filename: string
    
    @return: `datetime.datetime` object or None if *filename* couldn't be
        parsed
    
    """
    import re
    from datetime import timedelta
    match = re.match('(.*)_(\d{12})\+(\d{3})H(\d{2})M', filename)
    if not match:
        return None
    (name, #@UnusedVariable
     analysis_time, forecast_hours, forecast_minutes) = match.groups()
    base = parse_datetime_arg(analysis_time)
    forecast_length = timedelta(0, (int(forecast_hours) * 60 +
                                    int(forecast_minutes)) * 60)
    return base + forecast_length

def fetch_for_times(times, N=2, maxdiff_hours=12):
    """
    Fetch the *N* nearest-in-time NWP files to each entry in *times*
    
    @param times: list of times for which to get NWP files
    @type times: list of `datetime.datetime` instances
    @param maxdiff_hours: maximum difference in time
        (abs(*time* - nwptime) < maxdiff_hours)
    @type maxdiff_hours: float
    
    @return: None
    
    """
    from getpass import getpass
    
    ftp = ftplib.FTP(SERVER)
    
    passwd = getpass("Enter password for %s@%s: " % (USER, SERVER))
    print(ftp.login(USER, passwd))
    ftp.cwd(DIRECTORY)
    available = ftp.nlst()
    
    nwp_times = [(parse_nwp(f), f) for f in available]
    for time in times:
        diff_times = [(abs(time - nwp_time), f) for nwp_time, f in nwp_times]
        diff_times.sort()
        for diff_time, filename in diff_times[:N]: #@UnusedVariable
            if os.path.exists(filename):
                print("%r already exists" % filename)
            else:
                ftp_filesize = ftp.size(filename)
                with open(filename, 'wb') as f:
                    print "Fetching %r" % filename,
                    points = [x * .1 for x in range(0, 11)]
                    def write_chunk(chunk):
                        f.write(chunk)
                        for point in points:
                            if os.path.getsize(filename) > ftp_filesize * point:
                                print "%d %% " % (point * 100,),
                                points.pop(0)
                    ftp.retrbinary('RETR ' + filename, write_chunk)
                    print(' done')


def parse_datetime_arg(s):
    """
    Parse string *s* as a datetime YYYYMMDDhhmm and return a `datetime.datetime`
    object.
    
    @param s: YYYYMMDDhhmm
    @type s: string
    
    @return: `datetime.datetime`
    
    """
    from datetime import datetime
    
    try:
        year = int(s[:4])
        month = int(s[4:6])
        day = int(s[6:8])
        hour = int(s[8:10])
        minute = int(s[10:])
    except (IndexError, ValueError):
        raise ParseError("%r !~ YYYYMMDDhhmm" % s)
    
    return datetime(year, month, day, hour, minute)


def usage():
    print("Usage: %s YYYYMMDDhhmm [...]\n"
          "or     %s scene [...]" % sys.argv[0])
    sys.exit(-1)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        usage()
    
    times = []
    for dt_string in sys.argv[1:]:
        try:
            times.append(parse_datetime_arg(dt_string))
        except ParseError, err:
            try:
                from runutils import parse_scene
                satname, _datetime, orbit = parse_scene(
                                                    os.path.basename(dt_string))
                times.append(_datetime)
            except ValueError:
                print(err)
                usage()
    fetch_for_times(times)
