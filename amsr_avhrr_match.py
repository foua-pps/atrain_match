"""
Script for matching AMSR-E and AVHRR swaths

"""

import os
import logging
logger = logging.getLogger(__name__)
from runutils import process_scenes
from validate_cph import CPP_PHASE_VALUES_v2012, CPP_PHASE_VALUES
from cloudsat_calipso_avhrr_match import find_files_from_avhrr
from config import PPS_FORMAT_2012_OR_EARLIER
from common import MatchupError


#: Should results be plotted?
_PLOTTING = False

#: Directory for mapper files
MATCH_DIR = os.environ.get('MATCH_DIR', '.')

#Time threshold, i.e. max time diff to be considered as a match
TIME_THR=600.0
#NB! I'm not sure there really is a sorting due to TIME_THR. Do check in
#    the time_diff dataset after running! /Sara Hornquist 2012-06-18

#: Radius of AMSR-E footprint (m)
AMSR_RADIUS = 10e3

#No-data value for lwp and reff
NODATA_CPP = 65535


def process_noaa_scene(avhrr_filename, options, amsr_filename=None, ctype=None,
                       reff_max=None, lwp_max=None, water=False,
                       n_neighbours=8):
    from pps_basic_configure import AVHRR_DIR, OUTPUT_DIR, AUX_DIR

    import pps_arguments 
    #ppsarg = pps_arguments.create_pps_argument(pps_arguments.USE_SUNSATANGLES_GLOB,
    #                                           pps_arguments.SATELLITE_PROJ,
    #                                           satname,orbit)
    #avhrr_filename = os.path.join(AVHRR_DIR, ppsarg.files.avhrr)
    if amsr_filename:
        amsr_filenames = [amsr_filename]
    else:
        from amsr_avhrr.match import find_amsr
        amsr_filenames = find_amsr(avhrr_filename)
        logger.debug("Found AMSR-E files: %r" % amsr_filenames)

    try:
        pps_files = find_files_from_avhrr(avhrr_filename, options)
    except MatchupError:
        #Stop this case, but allow other cases to go on
        logger.warning("Can not find required PPS files for case %s" % \
                       os.path.basename(avhrr_filename))
        tmp_file = os.path.join(MATCH_DIR,
                                "WARNING_missing_PPSfile_for_%s" % \
                                (os.path.basename(avhrr_filename)))
        cmdstr = "touch %s"%(tmp_file)
        os.system(cmdstr)
        return
    
    physiography_filename  = pps_files.physiography
    sunsat_filename  = pps_files.sunsatangles
    cpp_filename = pps_files.cpp    
    #cpp_filename = os.path.join(OUTPUT_DIR, ppsarg.files.cpp)
    #physiography_filename = os.path.join(AUX_DIR, ppsarg.files.physiography)
    if ctype is not None:
        #ctype_filename = os.path.join(OUTPUT_DIR, ppsarg.files.pge02)
        ctype_filename = pps_files.cloudtype
    else:
        ctype_filename = None

    for amsr_filename in amsr_filenames:
        process_case(amsr_filename, avhrr_filename, cpp_filename, sunsat_filename,
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


def process_case(amsr_filename, avhrr_filename, cpp_filename, sunsat_filename,
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
        mapper = match(amsr_filename, avhrr_filename, sunsat_filename,
                       radius_of_influence=AMSR_RADIUS,
                       time_threshold=TIME_THR,
                       n_neighbours=n_neighbours)
        mapper.write(match_path)
        logger.info("Match written to %r" % match_path)
    
    if False:
        logger.warning("Replacing '.h5' extension with '.hdf' in CPP filename, "
                       "due to a bug in CPP.")
        cpp_filename = '.hdf'.join(cpp_filename.rsplit('.h5', 1)) # last '.h5'
    if os.path.exists(cpp_filename):
        compare_lwps(mapper, amsr_filename, cpp_filename, sunsat_filename,
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
        fig_base = _fig_base(match_path)
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
        if PPS_FORMAT_2012_OR_EARLIER:
            reff = get_cpp_product(cpp_filename, 'reff')
        else:
            # Convert from m to mu
            reff_tmp = get_cpp_product(cpp_filename, 'cpp_reff')
            reff = np.where((reff_tmp == NODATA_CPP),
                            NODATA_CPP,
                            1000000.0 * reff_tmp)
        selection &= mapper(reff < reff_max)
        restrictions.append('effective radius < %.2g' % reff_max)
    
    if water:
        if PPS_FORMAT_2012_OR_EARLIER:
            phase = get_cpp_product(cpp_filename, 'cph')
            selection &= mapper(phase == CPP_PHASE_VALUES_v2012['liquid'])
        else:
            phase = get_cpp_product(cpp_filename, 'cpp_phase')
            selection &= mapper(phase == CPP_PHASE_VALUES['liquid'])
        restrictions.append('CPP phase is water')
    
    amsr_lwp_3d = amsr_lwp.reshape(amsr_lwp.shape[0], amsr_lwp.shape[1], 1)
    selection &= 0 <= amsr_lwp_3d # Remove pixels with negative lwp (nodata)
    if lwp_max is None:
        restrictions.append('0 <= AMSR-E lwp')
    else:
        selection &= amsr_lwp_3d < lwp_max
        restrictions.append('0 <= AMSR-E lwp < %.2g' % lwp_max)

    selection &= cpp_cwp >= 0 # Remove pixels with negative cwp (nodata)
    selection &= cpp_cwp != NODATA_CPP # Remove pixels with cwp==nodata
    restrictions.append('CPP cwp >= 0')


    selection.fill_value = False # masked values are never part of selection
    
    return selection, restrictions

def compare_lwps(mapper, amsr_filename, cpp_filename, sunsat_filename,
                 physiography_filename, avhrr_filename, ctype=None,
                 ctype_filename=None, reff_max=None, lwp_max=None,
                 water=False):
    """
    Compare liquid water paths in *amsr_filename* and *cpp_filename*, with
    matching in *mapper*. Sea mask is taken from *physiography_filename*.
    
    If plotting is on, AVHRR lon/lat is read from *avhrr_filename*.
    
    """
    import numpy as np
    from amsr_avhrr.util import get_cpp_product
    #SHq: originally the program read cwp, but we should validate lwp.
    #     They are treated the same way, so it is easy to change. The
    #     variables in this program are still called cwp, anyway.
    #cpp_cwp_avhrr_proj = get_cpp_product(cpp_filename, 'cwp')
    if PPS_FORMAT_2012_OR_EARLIER:
        cpp_cwp_avhrr_proj = get_cpp_product(cpp_filename, 'lwp')
    else:
        cwp_tmp = get_cpp_product(cpp_filename, 'cpp_lwp')
        #Convert from kg/m2 to g/m2
        cpp_cwp_avhrr_proj = np.where((cwp_tmp == NODATA_CPP),
                                      NODATA_CPP,
                                      1000.0 * cwp_tmp)
    cpp_cwp = mapper(cpp_cwp_avhrr_proj)

    from amsr_avhrr.util import get_amsr_lwp, get_amsr_lonlat
    amsr_lwp = get_amsr_lwp(amsr_filename)
    lon, lat = get_amsr_lonlat(amsr_filename) #@UnusedVariable

    #Read time for amsr and avhrr, and create a time_diff
    from amsr_avhrr.util import get_amsr_time, get_avhrr_time
    avhrr_time = get_avhrr_time(sunsat_filename)
    amsr_time = get_amsr_time(amsr_filename)

    amsr_time_expanded=np.ones(amsr_lwp.shape)
    for i in range(amsr_lwp.shape[0]):
        amsr_time_expanded[i]=np.ones(amsr_lwp.shape[1])*amsr_time[i]
    avhrr_time_expanded=np.ones(cpp_cwp_avhrr_proj.shape)
    for i in range(cpp_cwp_avhrr_proj.shape[0]):
        avhrr_time_expanded[i]=np.ones(cpp_cwp_avhrr_proj.shape[1])*avhrr_time[i]
    avhrr_time_remapped = mapper(avhrr_time_expanded)
    avhrr_time_remapped = avhrr_time_remapped.mean(axis=-1)

    time_diff=avhrr_time_remapped-amsr_time_expanded
    #for select_pixels we need an expanded timediff
    time_diff_expanded=np.ones(cpp_cwp.shape)
    for i in range(cpp_cwp.shape[0]):
        for j in range(cpp_cwp.shape[1]):
            time_diff_expanded[i][j]=np.ones(cpp_cwp.shape[2])*time_diff[i][j]
    
        
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
        diff_file = os.path.join(MATCH_DIR,
                                 (_fig_base(_match_file(amsr_filename,
                                                        avhrr_filename)) +
                                  'lwp_diff.h5'))
        write_data(lwp_diff, 'lwp_diff', diff_file, mode='w',
                   attributes={'restrictions': restrictions})
        selection_2d = selection.all(axis=-1).filled(False)
        write_data(cpp_cwp.mean(axis=-1)[selection_2d],
                   'cpp_cwp', diff_file, mode='a')
        write_data(amsr_lwp[selection_2d], 'amsr_lwp', diff_file,
                   mode='a')
        write_data(amsr_time_expanded[selection_2d], 'amsr_time', diff_file,
                   mode='a')
        write_data(avhrr_time_remapped[selection_2d], 'avhrr_time', diff_file,
                   mode='a')
        write_data(time_diff[selection_2d], 'time_diff', diff_file,
                   mode='a')
        
        write_data(selection_2d, 'selection', diff_file, mode='a')
        write_data(lon[selection_2d], 'longitudes', diff_file, mode='a')
        write_data(lat[selection_2d], 'latitudes', diff_file, mode='a')
    
    if _PLOTTING:
        title = _plot_title(amsr_filename, avhrr_filename)
        fig_base = os.path.join(MATCH_DIR,
                                _fig_base(_match_file(amsr_filename,
                                                      avhrr_filename)))
        
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
        avhrr_lonlat = get_avhrr_lonlat(sunsat_filename)
        from amsr_avhrr.plotting import Field, plot_fields
        #SHQ change from cwp to lwp
        #fields = [Field(cpp_cwp_avhrr_proj, desc='CPP cwp', *avhrr_lonlat),
        fields = [Field(cpp_cwp_avhrr_proj, desc='CPP lwp', *avhrr_lonlat),
                  Field(amsr_lwp, lon, lat, desc='AMSR-E lwp')]
        fig = plot_fields(fields, break_value=lwp_max)
        fig.suptitle(title)
        fig.set_size_inches(20, 10)
        fig.savefig(fig_base + "lwp_swaths.png")


if __name__ == '__main__':
    
    import ConfigParser
    CONFIG_PATH = os.environ.get('ATRAINMATCH_CONFIG_DIR', './etc')
    CONF = ConfigParser.ConfigParser()
    CONF.read(os.path.join(CONFIG_PATH, "atrain_match.cfg"))
    OPTIONS = {}
    for option, value in CONF.items('general', raw = True):
        OPTIONS[option] = value
        
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
        #satname, orbit = args[1:]
        #satname = args[1]
        #orbit = int(args[2])
        process_noaa_scene(filename, orbit, **processing_kwargs)
    else:
        process_scenes(args, process_noaa_scene, OPTIONS, ignore_errors=not opts.debug,
                       **processing_kwargs)
