"""
Run AAPP on one or more EARS AVHRR segments

"""

import os

_months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct',
           'nov', 'dec']
_MONTHS = [month.upper() for month in _months]

#: Directory to put symlink to hrpt file in
AAPP_incoming = "/data/proj/safworks/jakob/install/AAPP/working_dir/incoming"


#: lower case mappings of short platform names
_platform_shorts = {'n': 'noaa',
                    'm': 'metop'}

def parse_segment_filename(filename):
    """
    Parse *filename* according to the format of EUMETCAST EARS AVHRR segment
    filenames (e.g. 'avhrr_20100112_152600_noaa19.hrp' or
    'N19-AVEA-AVHEAR00-NA-NA-20100112152600.000000000Z-1011664').
    
    Returns (platform, satellite number, `datetime.datetime` instance)
    
    """
    import re
    from datetime import datetime
    
    basename = os.path.basename(filename).lower()
    
    # Try to match ~'avhrr_20100112_152600_noaa19.hrp'
    m = re.match('avhrr_(\d{8})_(\d{6})_([a-z]+)(\d*)\..*', basename)
    if m is not None:
        date, time, platform, satnumber = m.groups()
        datetime_s = date + time
    else:
        # Try to match ~'N19-AVEA-AVHEAR00-NA-NA-20100112152600.000000000Z-1011664'
        m = re.match('([a-z])(\d+).*(\d{14})\..*', basename)
        if m is not None:
            short, satnumber, datetime_s = m.groups()
            try:
                platform = _platform_shorts[short]
            except KeyError:
                raise NotImplementedError("Platform %r not understood" % short)
        else:
            raise ValueError("Couldn't parse segment filename %r" % filename)
    
    _datetime = datetime.strptime(datetime_s, "%Y%m%d%H%M%S")
    
    return platform, int(satnumber), _datetime


def cat(files):
    """
    Concatenate *files*
    
    """
    with open(files[0] + '.full', 'wb') as out:
        for file in files:
            with open(file, 'rb') as f:
                out.write(f.read())
    
    return out.name


def prepare_hrpt(filename):
    """
    Prepare HRPT file *filename* for AAPP.
    
    """
    platform, satnumber, datetime = parse_segment_filename(filename)
    
    aapp_filename = ("hrpt16_{platform}-{satnumber}_"
                     "{date.day:02d}-{monthname}-{date.year}_"
                     "{date.hour:02d}:{date.minute:02d}:{date.second:02d}."
                     "xxx_xxxx").format(platform=platform.upper(),
                                        satnumber=satnumber, date=datetime,
                                        monthname=_MONTHS[datetime.month - 1])
    
    aapp_filename = os.path.join(AAPP_incoming, aapp_filename)
    if os.path.exists(aapp_filename):
        print("%r already exists" % aapp_filename)
    else:
        abspath = os.path.abspath(filename)
        print("ln -s %r %r"
              % (abspath, aapp_filename))
        os.symlink(abspath, aapp_filename)


def prepare():
    from optparse import OptionParser
    
    parser = OptionParser(usage="%(progname)s segment [segment ...]")
    parser.add_option('-c', '--concatenate', help="Concatenate segment files",
                      action='store_true')
    opts, filenames = parser.parse_args()
    
    if len(filenames) == 0:
        parser.error("No segment filenames provided")
    if len(filenames) > 1:
        if opts.concatenate:
            print("Concatenating segments")
            filenames = [cat(filenames)]
    
    for filename in filenames:
        prepare_hrpt(filename)


if __name__ == '__main__':
    prepare()
