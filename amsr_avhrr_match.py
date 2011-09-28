"""
Script for matching AMSR-E and AVHRR swaths

"""

import os
import logging
logger = logging.getLogger(__name__)
from runutils import process_scenes


#: Should results be plotted?
_PLOTTING = False

#: Directory for mapper files
MATCH_DIR = os.environ.get('MATCH_DIR', '.')


def process_noaa_scene(satname, orbit, amsr_filename=None, ctype=None):
    from pps_runutil import get_ppsProductArguments
    from pps_basic_configure import AVHRR_DIR, OUTPUT_DIR, AUX_DIR
    
    #argv = [sys.argv[0], 'satproj', satname, orbit]
    # get_ppsProductArguments can't take non-string orbit
    argv = ['', 'satproj', satname, str(orbit)]
    ppsarg, arealist = get_ppsProductArguments(argv) #@UnusedVariable
    avhrr_filename = os.path.join(AVHRR_DIR, ppsarg.files.avhrr)
    if amsr_filename:
        amsr_filenames = [amsr_filename]
    else:
        from amsr_avhrr.match import find_amsr
        amsr_filenames = find_amsr(avhrr_filename)
        logger.debug("Found AMSR-E files: %r" % amsr_filenames)
    
    cpp_filename = os.path.join(OUTPUT_DIR, ppsarg.files.cpp)
    physiography_filename = os.path.join(AUX_DIR, ppsarg.files.physiography)
    if ctype is not None:
        ctype_filename = os.path.join(OUTPUT_DIR, ppsarg.files.pge02)
    else:
        ctype_filename = None
    
    for amsr_filename in amsr_filenames:
        process_case(amsr_filename, avhrr_filename, cpp_filename,
                      physiography_filename, ctype, ctype_filename)


def _match_file(amsr_filename, avhrr_filename):
    return "match--%s--%s.h5" % (os.path.basename(avhrr_filename),
                                 os.path.basename(amsr_filename))

def _fig_base(match_file):
    return match_file.rsplit('.h5', 1)[0] + '--'

def _plot_title(amsr_filename, avhrr_filename):
    "%r\nmatched with\n%r" % (os.path.basename(avhrr_filename),
                              os.path.basename(amsr_filename))


def process_case(amsr_filename, avhrr_filename, cpp_filename,
                  physiography_filename, ctype=None, ctype_filename=None):
    """
    Match, plot, and validate scene defined by the given files.
    
    """
    match_file = _match_file(amsr_filename, avhrr_filename)
    match_path = os.path.join(MATCH_DIR, match_file)
    if os.path.exists(match_path):
        logger.info("Reading match from %r" % match_path)
        from amsr_avhrr.match import MatchMapper
        mapper = MatchMapper.from_file(match_path)
    else:
        logger.info("Matching AMSR-E and AVHRR swaths")
        from amsr_avhrr.match import match
        mapper = match(amsr_filename, avhrr_filename)
        mapper.write(match_path)
        logger.info("Match written to %r" % match_path)
    
    if False:
        logger.warning("Replacing '.h5' extension with '.hdf' in CPP filename, "
                       "due to a bug in CPP.")
        cpp_filename = '.hdf'.join(cpp_filename.rsplit('.h5', 1)) # last '.h5'
    if os.path.exists(cpp_filename):
        compare_lwps(mapper, amsr_filename, cpp_filename,
                     physiography_filename, avhrr_filename, ctype,
                     ctype_filename)
    else:
        logger.warning("No CPP product found")
    
    if _PLOTTING:
        logger.debug("Plotting time difference")
        from amsr_avhrr.util import get_amsr_lonlat
        lon, lat = get_amsr_lonlat(amsr_filename)
        from amsr_avhrr.plotting import plot_array
        fig = plot_array(lon, lat, mapper.time_diff,
                         legend="time difference (s)")
        fig_base = _fig_base(match_file)
        fig.set_size_inches(20, 12)
        fig.suptitle(_plot_title(amsr_filename, avhrr_filename))
        fig.savefig(fig_base + "time_diff.png")


def write_data(data, name, filename, mode=None, attributes=None):
    import h5py
    with h5py.File(filename, mode) as f:
        d = f.create_dataset(name, data=data.compressed())
        if attributes:
            for k, v in attributes.items():
                d.attrs[k] = v


def compare_lwps(mapper, amsr_filename, cpp_filename,
                 physiography_filename, avhrr_filename=None, ctype=None,
                 ctype_filename=None):
    """
    Compare liquid water paths in *amsr_filename* and *cpp_filename*, with
    matching in *mapper*. Sea mask is taken from *physiography_filename*.
    
    If plotting is on, AVHRR lon/lat is read from *avhrr_filename*.
    
    """
    
    from amsr_avhrr.util import get_cpp_lwp
    cpp_lwp_avhrr_proj = get_cpp_lwp(cpp_filename)
    cpp_lwp = mapper(cpp_lwp_avhrr_proj)
    
    from amsr_avhrr.util import get_amsr_lwp, get_amsr_lonlat
    amsr_lwp = get_amsr_lwp(amsr_filename)
    lon, lat = get_amsr_lonlat(amsr_filename) #@UnusedVariable
    
    from epshdf import read_physiography
    landuse =  read_physiography(physiography_filename, 1, 0, 0).landuse
    
    # Sea pixels in AMSR-E swath
    sea = mapper(landuse == 16)
    
    # Select only sea pixels (AMSR-E lwp is only available over sea)
    selection = sea
    restrictions = ['sea']
    
    # Select only pixels with cloud type *ctype*
    if ctype is not None and ctype_filename:
        from epshdf import read_cloudtype
        ctype_obj = read_cloudtype(ctype_filename, 1, 0, 0)
        selection &= (mapper(ctype_obj.cloudtype == ctype))
        restrictions.append('cloud type == %d' % ctype)
    
    from amsr_avhrr.validation import validate_lwp
    lwp_diff, restrictions = validate_lwp(amsr_lwp, cpp_lwp, selection,
                                          restrictions)
    if lwp_diff is None:
        logger.warning("No matches with restrictions: %s" %
                       '; '.join(restrictions))
    else:
        diff_file = (_fig_base(_match_file(amsr_filename, avhrr_filename)) +
                     'lwp_diff.h5')
        write_data(lwp_diff, 'lwp_diff', diff_file, mode='w',
                   attributes={'restrictions': restrictions})
    
    if _PLOTTING:
        title = _plot_title(amsr_filename, avhrr_filename)
        fig_base = _fig_base(_match_file(amsr_filename, avhrr_filename))
        
        logger.debug("Plotting lwp arrays")
        from amsr_avhrr.plotting import imshow_lwps
        fig = imshow_lwps(amsr_lwp, cpp_lwp.mean(axis=-1),
                          mapper.time_diff.mean(axis=-1),
                          sea.mean(axis=-1) > .5)
        fig.suptitle(title)
        fig.set_size_inches(20, 12)
        fig.savefig(fig_base + "lwp_arrays.png")
        
        logger.debug("Plotting lwp swaths")
        from amsr_avhrr.util import get_avhrr_lonlat
        avhrr_lonlat = get_avhrr_lonlat(avhrr_filename)
        from amsr_avhrr.plotting import Field, plot_fields
        fields = [Field(cpp_lwp_avhrr_proj, desc='CPP cwp', *avhrr_lonlat),
                  Field(amsr_lwp, lon, lat, desc='AMSR-E lwp')]
        fig = plot_fields(fields)
        fig.suptitle(title)
        fig.set_size_inches(20, 10)
        fig.savefig(fig_base + "lwp_swaths.png")


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="Usage: %prog [options] satproj "
                                "<satname> <orbit> \n"
                                "or:    %prog [options] CASE [...]\n"
                                "          where CASE ~ 'possible/path/"
                                "noaa19_20110916_0959_03456*'")
    parser.add_option('-p', '--plot', help="Create plots", action='store_true')
    parser.add_option('-v', '--verbose', action='store_true')
    parser.add_option('-d', '--debug', action='store_true',
                      help="Don't ignore errors")
    parser.add_option('-c', '--cloudtype', type='int',
                      help="Screen by CLOUDTYPE (integer value)")
    opts, args = parser.parse_args()
    
    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
        logger.debug("Verbose")
    else:
        logging.basicConfig(level=logging.INFO)
    
    if opts.plot:
        logger.debug("Plotting enabled")
        _PLOTTING = True
    
    # Command line handling
    if args[0] == 'satproj':
        satname, orbit = args[1:]
        process_noaa_scene(satname, orbit, ctype=opts.cloudtype)
    else:
        process_scenes(args, process_noaa_scene, ignore_errors=not opts.debug,
                       ctype=opts.cloudtype)