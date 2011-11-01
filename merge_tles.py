"""
Merge TLE files

"""
import sys


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


def _parse_celestrak(lines):
    # Get the satellite name from the first header line (not counting '===='...)
    h1 = lines[1].split()
    satname = h1[2:]
    name = ''.join(satname)
    
    tles = []
    for ix in range(5, len(lines), 2):
        if lines[ix] != '<End of file>':
            l1, l2 = lines[ix:ix + 2]
            tles.append((tle_time(l1), name, l1, l2))
    
    return tles

def _parse_tlinset(lines):
    tles = []
    for ix in range(0, len(lines), 3):
        name, l1, l2 = lines[ix:ix + 3]
        tles.append((tle_time(l1), name, l1, l2))
    
    return tles


def get_tles(filename):
    """
    Read *filename* and return a list of all TLEs as tuples (datetime,
    satellite name, line 1, line2).
    
    """
    import re
    
    with open(filename, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    
    # Check format of file
    if lines[0][:20] == '=' * 20:
        # TLE request file from celestrack
        return _parse_celestrak(lines)
    elif (re.match('[\w]+$', lines[0]) and lines[1][:2] == '1 ' and
          lines[2][:2] == '2 '):
        # TLINSET file suitable for snotimes
        return _parse_tlinset(lines)
    else:
        raise ValueError("Format of %r was not recognized" % filename)


def merge_files(filenames, outfile=None):
    """
    Merge all TLEs in all files in *filenames*. If *outfile* is None, print to
    stdout.
    
    """
    tles = set()
    for filename in filenames:
        tles.update(get_tles(filename))
    
    if not outfile:
        outfile = sys.stdout
    else:
        outfile = open(outfile, 'w')
    
    for tle in sorted(tles):
        tle_time, name, l1, l2 = tle
        outfile.write(name.upper() + '\n')
        outfile.write(l1 + '\n')
        outfile.write(l2 + '\n')
        outfile.flush()
    
    if outfile is not sys.stdout:
        outfile.close()


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser("Usage: %prog [-h] [options] FILE [FILE ...]")
    parser.add_option('-o', '--outfile',
                      help="write merged TLEs to OUTFILE instead of stdout")
    
    opts, args = parser.parse_args()
    
    filenames = args
    if len(filenames) == 0:
        parser.print_usage()
        sys.exit(1)
    
    merge_files(filenames, opts.outfile)
