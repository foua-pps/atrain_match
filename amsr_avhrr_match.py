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

#: Radius of AMSR-E footprint (m)
AMSR_RADIUS = 10e3


def process_noaa_scene(satname, orbit, amsr_filename=None, ctype=None,
                       reff_max=None, lwp_max=None, water=False,
                       n_neighbours=8):
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
                      physiography_filename, ctype, ctype_filename, reff_max,
                      lwp_max, water, n_neighbours)


def _match_file(amsr_filename, avhrr_filename):
    return "match--%s--%s.h5" % (os.path.basename(avhrr_filename),
                                 os.path.basename(amsr_filename))

def _fig_base(match_file):
    return match_file.rsplit('.h5', 1)[0] + '--'

def _plot_title(amsr_filename, avhrr_filename):
    "%r\nmatched with\n%r" % (os.path.basename(avhrr_filename),
                              os.path.basename(amsr_filename))


def process_case(amsr_filename, avhrr_filename, cpp_filename,
                  physiography_filename, ctype=None, ctype_filename=None,
                  reff_max=None, lwp_max=None, water=False, n_neighbours=8):
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
        mapper = match(amsr_filename, avhrr_filename,
                       radius_of_influence=AMSR_RADIUS,
                       n_neighbours=n_neighbours)
        mapper.write(match_path, compression=_COMPRESSION)
        logger.info("Match written to %r" % match_path)
    
    if False:
        logger.warning("Replacing '.h5' extension with '.hdf' in CPP filename, "
                       "due to a bug in CPP.")
        cpp_filename = '.hdf'.join(cpp_filename.rsplit('.h5', 1)) # last '.h5'
    if os.path.exists(cpp_filename):
        compare_lwps(mapper, amsr_filename, cpp_filename,
                     physiography_filename, avhrr_filename, ctype,
                     ctype_filename, reff_max, lwp_max, water)
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


def get_sea(mapper, physiography_filename):
    """
    Get sea map from *physiography_filename*, mapped to *mapper*'s target.
    
    """
    from epshdf import read_physiography
    landuse = read_physiography(physiography_filename, 1, 0, 0).landuse
    
    # Sea pixels in AMSR-E swath
    return mapper(landuse == 16)


def select_pixels(mapper, amsr_lwp, cpp_cwp, sea,
                  lwp_max=None,
                  ctype=None, ctype_filename=None,
                  reff_max=None, cpp_filename=None, water=False):
    """
    Find valid pixels
    
    """
    from amsr_avhrr.util import get_cpp_product
    # Select only sea pixels (AMSR-E lwp is only available over sea)
    selection = sea
    restrictions = ['sea']
    
    # Select only pixels with cloud type *ctype*
    if ctype is not None:
        if not ctype_filename:
            raise ValueError("Need *ctype_filename* for screening on *ctype*")
        from epshdf import read_cloudtype
        ctype_obj = read_cloudtype(ctype_filename, 1, 0, 0)
        selection &= mapper(ctype_obj.cloudtype == ctype)
        restrictions.append('cloud type == %d' % ctype)
    
    if reff_max is not None:
        reff = get_cpp_product(cpp_filename, 'reff')
        selection &= mapper(reff < reff_max)
        restrictions.append('effective radius < %.2g' % reff_max)
    
    if water:
        phase = get_cpp_product(cpp_filename, 'cph')
        selection &= mapper(phase == 1)
        restrictions.append('CPP phase is water')
    
    amsr_lwp_3d = amsr_lwp.reshape(amsr_lwp.shape[0], amsr_lwp.shape[1], 1)
    selection &= 0 <= amsr_lwp_3d # Remove pixels with negative lwp (nodata)
    if lwp_max is None:
        restrictions.append('0 <= AMSR-E lwp')
    else:
        selection &= amsr_lwp_3d < lwp_max
        restrictions.append('0 <= AMSR-E lwp < %.2g' % lwp_max)
    
    selection &= cpp_cwp >= 0 # Remove pixels with negative cwp (nodata)
    restrictions.append('CPP cwp >= 0')
    
    selection.fill_value = False # masked values are never part of selection
    
    return selection, restrictions

def compare_lwps(mapper, amsr_filename, cpp_filename,
                 physiography_filename, avhrr_filename, ctype=None,
                 ctype_filename=None, reff_max=None, lwp_max=None,
                 water=False):
    """
    Compare liquid water paths in *amsr_filename* and *cpp_filename*, with
    matching in *mapper*. Sea mask is taken from *physiography_filename*.
    
    If plotting is on, AVHRR lon/lat is read from *avhrr_filename*.
    
    """
    
    from amsr_avhrr.util import get_cpp_product
    cpp_cwp_avhrr_proj = get_cpp_product(cpp_filename, 'cwp')
    cpp_cwp = mapper(cpp_cwp_avhrr_proj)
    
    from amsr_avhrr.util import get_amsr_lwp, get_amsr_lonlat
    amsr_lwp = get_amsr_lwp(amsr_filename)
    lon, lat = get_amsr_lonlat(amsr_filename) #@UnusedVariable
    
    sea = get_sea(mapper, physiography_filename)
    
    selection, restrictions = select_pixels(mapper, amsr_lwp, cpp_cwp, sea,
                                            lwp_max, ctype, ctype_filename,
                                            reff_max, cpp_filename, water)
    logger.debug("Selected pixels: %s" % '; '.join(restrictions))
    
    from amsr_avhrr.validation import validate_lwp
    lwp_diff = validate_lwp(amsr_lwp, cpp_cwp, selection)
    if lwp_diff is None:
        logger.warning("No matches with restrictions: %s" %
                       '; '.join(restrictions))
    else:
        from amsr_avhrr.util import write_data
        diff_file = (_fig_base(_match_file(amsr_filename, avhrr_filename)) +
                     'lwp_diff.h5')
        write_data(lwp_diff, 'lwp_diff', diff_file, mode='w',
                   attributes={'restrictions': restrictions})
        selection_2d = selection.all(axis=-1).filled(False)
        write_data(cpp_cwp.mean(axis=-1)[selection_2d],
                   'cpp_cwp', diff_file, mode='a')
        write_data(amsr_lwp[selection_2d], 'amsr_lwp', diff_file,
                   mode='a')
        write_data(selection_2d, 'selection', diff_file, mode='a')
        write_data(lon[selection_2d], 'longitudes', diff_file, mode='a')
        write_data(lat[selection_2d], 'latitudes', diff_file, mode='a')
    
    if _PLOTTING:
        title = _plot_title(amsr_filename, avhrr_filename)
        fig_base = _fig_base(_match_file(amsr_filename, avhrr_filename))
        
        logger.debug("Plotting lwp arrays")
        from amsr_avhrr.plotting import imshow_lwps
        fig = imshow_lwps(amsr_lwp, cpp_cwp.mean(axis=-1),
                          mapper.time_diff.mean(axis=-1),
                          sea.mean(axis=-1) > .5, lwp_max=lwp_max)
        fig.suptitle(title)
        fig.set_size_inches(20, 12)
        fig.savefig(fig_base + "lwp_arrays.png")
        
        logger.debug("Plotting lwp swaths")
        from amsr_avhrr.util import get_avhrr_lonlat
        avhrr_lonlat = get_avhrr_lonlat(avhrr_filename)
        from amsr_avhrr.plotting import Field, plot_fields
        fields = [Field(cpp_cwp_avhrr_proj, desc='CPP cwp', *avhrr_lonlat),
                  Field(amsr_lwp, lon, lat, desc='AMSR-E lwp')]
        fig = plot_fields(fields, break_value=lwp_max)
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
                      help="Only include CLOUDTYPE (integer value)")
    parser.add_option('-r', '--reff_max', type='float',
                      help="Screen out effective radii > REFF_MAX")
    parser.add_option('-l', '--lwp_max', type='float',
                      help="Screen out AMSR-E liquid water path > LWP_MAX")
    parser.add_option('-n', '--neighbours', type='int',
                      help="Number of nearest AVHRR neighbours to use")
    parser.add_option('-w', '--water', action='store_true',
                      help="Select only CPP water phase pixels")
    opts, args = parser.parse_args()
    
    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
        logger.debug("Verbose")
    else:
        logging.basicConfig(level=logging.INFO)
    
    if opts.plot:
        logger.debug("Plotting enabled")
        _PLOTTING = True
    
    processing_kwargs = {}
    if opts.cloudtype is not None:
        processing_kwargs['ctype'] = opts.cloudtype
    if opts.reff_max is not None:
        processing_kwargs['reff_max'] = opts.reff_max
    if opts.lwp_max is not None:
        processing_kwargs['lwp_max'] = opts.lwp_max
    if opts.neighbours is not None:
        processing_kwargs['n_neighbours'] = opts.neighbours
    if opts.water:
        processing_kwargs['water'] = True
    
    # Command line handling
    if args[0] == 'satproj':
        satname, orbit = args[1:]
        process_noaa_scene(satname, orbit, **processing_kwargs)
    else:
        process_scenes(args, process_noaa_scene, ignore_errors=not opts.debug,
                       **processing_kwargs)
