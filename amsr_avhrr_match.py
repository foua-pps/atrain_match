"""
Match AMSR and AVHRR data

"""

from __future__ import with_statement
import numpy as np
from datetime import datetime
import os
from ppshdf_cloudproducts import PpsProduct
import logging
logger = logging.getLogger(__name__)

TAI93 = datetime(1993, 1, 1)


class MatchMapper(object):
    """
    Map arrays from one swath to another.
    
    """
    def __init__(self, rows, cols, mask):
        self.rows = rows
        self.cols = cols
        self.mask = mask
    
    def __call__(self, array):
        """
        Maps *array* to target swath.
        
        """
        return np.ma.array(array[self.rows, self.cols], mask=self.mask)


def get_amsr_lonlat(filename):
    """
    Get (lon, lat) from AMSR file *filename*.
    
    """
    import h5py
    
    with h5py.File(filename, 'r') as f:
        lon = f['Swath1/Geolocation Fields/Longitude'][:]
        lat = f['Swath1/Geolocation Fields/Latitude'][:]
    
    return lon, lat


def get_avhrr_lonlat(filename):
    """
    Get (lon, lat) from AVHRR file *filename* (remapped by PPS).
    
    """
    import pps_io
    
    geo = pps_io.readAvhrrGeoData(filename)
    
    return geo.longitude, geo.latitude


def match_lonlat(source, target, radius_of_influence=1e3):
    """
    Produce a masked array of the same shape as the arrays in *target*, with
    indices of nearest neighbours in *source*. *source* and *target* should be
    tuples (lon, lat) of the source and target swaths, respectively.
    
    Note::
    
        Fastest matching is obtained when *target* has lower resolution than
        *source*.
    
    """
    from pyresample.geometry import SwathDefinition
    from pyresample.kd_tree import get_neighbour_info
    
    source_def = SwathDefinition(*source)
    target_def = SwathDefinition(*target)
    
    valid_in, valid_out, indices, distances = get_neighbour_info( #@UnusedVariable
        source_def, target_def, radius_of_influence, 1)
    
    indices.shape = target_def.shape
    distances.shape = target_def.shape
    
    rows = indices // source_def.shape[0]
    cols = indices % source_def.shape[1]
    mask = distances > radius_of_influence
    
    return MatchMapper(rows, cols, mask)


def get_amsr_time(filename):
    """
    Get time of each scanline in AMSR file *filename*.
    
    Note::
    
        Time is converted from seconds since 1993-01-01 to seconds since
        1970-01-01.
    
    """
    import h5py
    
    with h5py.File(filename) as f:
        sec1993 = f['Swath1/Geolocation Fields/Time']['Time'][:]
    
    from calendar import timegm
    epoch_diff = timegm(TAI93.utctimetuple())
    
    sec1970 = sec1993.round().astype(np.int64) + epoch_diff
    
    return sec1970

def get_avhrr_time(filename):
    """
    Get time of each scanline in AVHRR file *filename*.
    
    """
    import pps_io
    
    geo = pps_io.readAvhrrGeoData(filename)
    
    n_scanlines = geo.longitude.shape[0]
    sec1970 = np.linspace(geo.sec1970_start, geo.sec1970_end, n_scanlines)
    
    return sec1970


def match(amsr_filename, avhrr_filename):
    """
    Find matching indices in AVHRR array for each element in AMSR swath.
    
    Returns three masked arrays (rows, columns, time_difference), where masked-
    out elements had too large differences in either space or time.
    
    """
    from config import RESOLUTION, sec_timeThr
    
    avhrr_lonlat = get_avhrr_lonlat(avhrr_filename)
    amsr_lonlat = get_amsr_lonlat(amsr_filename)
    
    mapper = match_lonlat(avhrr_lonlat, amsr_lonlat, RESOLUTION * 1e3)
    
    avhrr_time = get_avhrr_time(avhrr_filename)
    amsr_time = get_amsr_time(amsr_filename)
    
    time_diff = np.abs(avhrr_time[mapper.rows] - amsr_time.reshape((amsr_time.size, 1)))
    time_diff = np.ma.array(time_diff, mask=mapper.mask.copy())
    
    mapper.mask[time_diff > sec_timeThr] = True
    
    print("Time diff (min, max): %r" % ((time_diff.min(), time_diff.max()),))
    
    return mapper, time_diff


def adjust_lon(lon):
    """
    Analyse longitudes in *lon*, to see whether there is a gap somewhere in the
    middle of the range in *lon*. If there is a gap, shift longitudes up.
    
    Returns (lon_shifted, lon0)
    
    """
    n, bins = np.histogram(lon, bins=50)
    lon0 = -180
    if 0 in n:
        lo, hi = (n == 0).nonzero()[0][[0, -1]] + 1
        lon0 = np.mean((bins[lo], bins[hi]))
    return np.where(lon < lon0, lon + 360, lon), lon0


def plot_array(lon, lat, data, title, legend):
    """
    Make a plot of *data*
    
    """
    import matplotlib.pyplot as pl
    from mpl_toolkits.basemap import Basemap
    
    lon, lon0 = adjust_lon(lon)
    
    fig = pl.figure()
    fig.suptitle(title)
    ax = fig.add_subplot(111)
    
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=lon0, urcrnrlon=lon0 + 360, resolution='c', ax=ax)
    m.drawcoastlines(linewidth=.5)
    dashes = [1, 3]
    labels = [1, 1, 1, 1]
    parallels = range(-80, 81, 20)
    m.drawparallels(parallels, linewidth=.5, dashes=dashes, labels=labels)
    meridian0 = int(np.ceil(lon0 / 45.)) * 45
    meridians = range(meridian0, meridian0 + 360, 45)
    m.drawmeridians(meridians, linewidth=.5, dashes=dashes, labels=labels)
    
    x, y = m(lon, lat)
    mesh = m.pcolormesh(x, y, data)
    cbar = fig.colorbar(mesh)
    cbar.ax.set_ylabel(legend)


def get_amsr_lwp(filename):
    """
    Return liquid water path (lwp) from AMSR-E file *filename*. The units of lwp
    are converted from mm to g m**-2.
    
    """
    import h5py
    
    with h5py.File(filename, 'r') as f:
        lwp_mm = PpsProduct(f['Swath1/Data Fields/High_res_cloud'][:],
                            description='lwp (mm)',
                            gain=f['Swath1/Data Fields/'
                                   'High_res_cloud'].attrs['Scale'],
                            nodata=-9990, scale_up=True)
    
    density = 1e3 # Density of water
    return lwp_mm.array * density


def get_cpp_lwp(filename):
    """
    Return liquid water path (lwp) from PPS CPP file *filename*. Units of the
    returned array are g m**-2.
    
    """
    import h5py
    
    with h5py.File(filename, 'r') as f:
        lwp = PpsProduct(f['cwp'][:], description=f['cwp'].attrs['description'],
                            gain=f['cwp'].attrs['gain'],
                            nodata=f['cwp'].attrs['no_data_value'],
                            scale_up=True)
    
    return lwp.array


def imshow_lwps(amsr_lwp, cpp_lwp, time_diff, sea, title=None):
    """
    Show *amsr_lwp* and *cpp_lwp* side by side. *sea* is used to draw a
    background, and mask out any pixels which are not sea.
    
    """
    import matplotlib.pyplot as pl
    
    from utility_functions import broken_cmap
    from matplotlib.colors import ListedColormap
    
    fig = pl.figure()
    
    vmin = 0
    vmax = max([amsr_lwp.max(), cpp_lwp.max()])
    cmap = broken_cmap(np.array([vmin, vmax]), break_value=170)
    ground_sea_map = ListedColormap(['g', 'b'], name="ground/sea map")
    
    ax = fig.add_subplot(131)
    ax.imshow(np.ma.array(sea, mask=sea.mask + sea), cmap=ground_sea_map)
    im = ax.imshow(np.ma.array(amsr_lwp, mask=~sea),
                   vmin=vmin, vmax=vmax, cmap=cmap)
    cbar = fig.colorbar(im)
    cbar.set_label('g m**-2')
    ax.set_title('AMSR-E lwp')
    
    ax = fig.add_subplot(132, sharex=ax, sharey=ax)
    ax.imshow(np.ma.array(sea, mask=sea.mask + sea), cmap=ground_sea_map)
    im = ax.imshow(np.ma.array(cpp_lwp, mask=~sea),
                   vmin=vmin, vmax=vmax, cmap=cmap)
    cbar = fig.colorbar(im)
    cbar.set_label('g m**-2')
    ax.set_title('PPS CPP cwp')
    
    ax = fig.add_subplot(133, sharex=ax, sharey=ax)
    im = ax.imshow(time_diff)
    cbar = fig.colorbar(im)
    cbar.set_label('s')
    ax.set_title('Time difference')
    
    if title:
        fig.suptitle(title)


def _test():
    import sys
    from pps_runutil import get_ppsProductArguments
    from pps_ioproxy import pps_ioproxy
    from pps_basic_configure import AVHRR_DIR, OUTPUT_DIR
    
    amsr_filename = sys.argv.pop(1)
    ppsarg, arealist = get_ppsProductArguments(sys.argv) #@UnusedVariable
    ioproxy = pps_ioproxy(ppsarg)
    avhrr_filename = os.path.join(AVHRR_DIR, ppsarg.files.avhrr)
    
    mapper, time_diff = match(amsr_filename, avhrr_filename)
    
    lon, lat = get_amsr_lonlat(amsr_filename)
    title = "%r\nmatched with\n%r" % (os.path.basename(avhrr_filename),
                                      os.path.basename(amsr_filename))
    plot_array(lon, lat, time_diff, title=title, legend="time difference (s)")
    
    # Compare liquid water paths
    amsr_lwp = get_amsr_lwp(amsr_filename)
    #amsr_lwp[amsr_lwp < 0] = -1
    
    logger.warning("Replacing '.h5' with '.hdf' in CPP filename, due to a bug "
                   "in CPP.")
    cpp_filename = os.path.join(OUTPUT_DIR, ppsarg.files.cpp.replace('.h5', '.hdf'))
    if os.path.exists(cpp_filename):
        cpp_lwp = mapper(get_cpp_lwp(cpp_filename))
        
        landuse = ioproxy.getLanduseOnly()
        # Sea pixels in AMSR-E swath
        sea = mapper(landuse == 16)
        
        imshow_lwps(amsr_lwp, cpp_lwp, time_diff, sea, title=title)
        
        validate_lwp(amsr_lwp, cpp_lwp, sea)
    else:
        logger.warning("No CPP product found")
    
    import matplotlib.pyplot as pl
    pl.show()


def validate_lwp(amsr_lwp, cpp_lwp, sea, threshold=170):
    """
    Compare liquid water path, lwp, in *amsr_filename* and *cpp_filename* files.
    False pixels in *sea* are masked out. Only values above *threshold* (g
    m**-2) are considered.
    
    """
    amsr_masked = np.ma.array(amsr_lwp, mask=~sea + (amsr_lwp < threshold))
    cpp_masked = np.ma.array(cpp_lwp, mask=~sea)
    
    diff = amsr_masked - cpp_masked
    
    print('=' * 40)
    print("AMSR-E lwp - CPP cwp")
    print("Non-sea pixels screened out")
    print("Pixels where AMSR-E lwp < %r g m**-2 screened out" % threshold)
    print("Number of pixels in comparison: %d" % diff.compressed().size)
    print("bias:    %.4g" % diff.mean())
    print("std:     %.4g" % diff.std())
    print("rel std: %.4g %%" % abs(100. * diff.std() / diff.mean()))
    
    from matplotlib import pyplot as pl
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.hist(diff.compressed())


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    _test()
