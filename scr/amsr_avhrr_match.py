"""
Script for matching AMSR-E and AVHRR swaths

"""

import os
import logging
logger = logging.getLogger(__name__)
from runutils import process_scenes


#: Should results be plotted?
_PLOTTING = False


def process_noaa_scene(satname, orbit, amsr_filename=None):
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
    
    for amsr_filename in amsr_filenames:
        process_case(amsr_filename, avhrr_filename, cpp_filename,
                      physiography_filename)


def process_case(amsr_filename, avhrr_filename, cpp_filename,
                  physiography_filename):
    """
    Match, plot, and validate scene defined by the given files.
    
    """
    match_file = "match--%s--%s.h5" % (os.path.basename(avhrr_filename),
                                     os.path.basename(amsr_filename))
    if os.path.exists(match_file):
        logger.info("Reading match from %r" % match_file)
        from amsr_avhrr.match import MatchMapper
        mapper = MatchMapper.from_file(match_file)
    else:
        logger.info("Matching AMSR-E and AVHRR swaths")
        from amsr_avhrr.match import match
        mapper = match(amsr_filename, avhrr_filename)
        mapper.write(match_file)
        logger.info("Match written to %r" % match_file)
    
    if _PLOTTING:
        from amsr_avhrr.util import get_amsr_lonlat
        lon, lat = get_amsr_lonlat(amsr_filename)
        title = "%r\nmatched with\n%r" % (os.path.basename(avhrr_filename),
                                          os.path.basename(amsr_filename))
        from amsr_avhrr.plotting import plot_array
        fig = plot_array(lon, lat, mapper.time_diff,
                         legend="time difference (s)")
        fig_base = match_file.rsplit('.h5', 1)[0] + '--'
        fig.set_size_inches(20, 12)
        fig.suptitle(title)
        fig.savefig(fig_base + "time_diff.png")
    
    # Compare liquid water paths
    from amsr_avhrr.util import get_amsr_lwp
    amsr_lwp = get_amsr_lwp(amsr_filename)
    #amsr_lwp[amsr_lwp < 0] = -1
    
    if False:
        logger.warning("Replacing '.h5' extension with '.hdf' in CPP filename, "
                       "due to a bug in CPP.")
        cpp_filename = '.hdf'.join(cpp_filename.rsplit('.h5', 1)) # last '.h5'
    if os.path.exists(cpp_filename):
        from amsr_avhrr.util import get_cpp_lwp
        cpp_lwp_avhrr_proj = get_cpp_lwp(cpp_filename)
        cpp_lwp = mapper(cpp_lwp_avhrr_proj)
        
        from epshdf import read_physiography
        landuse =  read_physiography(physiography_filename, 1, 0, 0).landuse
        
        # Sea pixels in AMSR-E swath
        sea = mapper(landuse == 16)
        
        if _PLOTTING:
            from amsr_avhrr.plotting import imshow_lwps
            fig = imshow_lwps(amsr_lwp, cpp_lwp, mapper.time_diff, sea)
            fig.suptitle(title)
            fig.set_size_inches(20, 12)
            fig.savefig(fig_base + "lwp_arrays.png")
            
            # Plot fields in individual plots
            from amsr_avhrr.util import get_avhrr_lonlat
            avhrr_lonlat = get_avhrr_lonlat(avhrr_filename)
            from amsr_avhrr.plotting import Field, plot_fields
            fields = [Field(cpp_lwp_avhrr_proj, desc='CPP cwp', *avhrr_lonlat),
                      Field(amsr_lwp, lon=lon, lat=lat, desc='AMSR-E lwp')]
            fig = plot_fields(fields)
            fig.suptitle(title)
            fig.set_size_inches(20, 10)
            fig.savefig(fig_base + "lwp_swaths.png")
        
        from amsr_avhrr.match import validate_lwp
        validate_lwp(amsr_lwp, cpp_lwp, sea, plotting=_PLOTTING)
    else:
        logger.warning("No CPP product found")


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="Usage: (progname) [options] satproj "
                                "<satname> <orbit> \n"
                                "or:    (progname) [options] case_filename "
                                "[...]")
    parser.add_option('-p', '--plot', help="Create plots", action='store_true')
    parser.add_option('-v', '--verbose', action='store_true')
    opts, args = parser.parse_args()
    
    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    if opts.plot:
        _PLOTTING = True
    
    # Command line handling
    if args[0] == 'satproj':
        satname, orbit = args[1:]
        process_noaa_scene(satname, orbit)
    else:
        process_scenes(args, process_noaa_scene)