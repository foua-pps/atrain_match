"""
Match AMSR and AVHRR data

"""

from __future__ import with_statement
import numpy as np
import matplotlib
matplotlib.use('agg')
from datetime import datetime
import os
import logging
logger = logging.getLogger(__name__)

TAI93 = datetime(1993, 1, 1)


class MatchMapper(object):
    """
    Map arrays from one swath to another.
    
    """
    def __init__(self, rows, cols, pixel_mask, time_diff=None,
                 time_threshold=None):
        self._rows = rows
        self._cols = cols
        self._pixel_mask = pixel_mask
        self._time_diff = time_diff
        self.time_threshold = time_threshold
    
    def __call__(self, array):
        """
        Maps *array* to target swath.
        
        """
        return np.ma.array(array[self.rows, self.cols], mask=self.mask)
    
    @property
    def rows(self):
        return np.ma.array(self._rows, mask=self.mask, fill_value=-1,
                           hard_mask=True)
    
    @property
    def cols(self):
        return np.ma.array(self._cols, mask=self.mask, fill_value=-1,
                           hard_mask=True)
    
    @property
    def time_diff(self):
        """Time difference in seconds"""
        if self._time_diff is None:
            return None
        # Only use pixel mask
        return np.ma.array(self._time_diff, mask=self._pixel_mask,
                           fill_value=np.inf, hard_mask=True)
    
    @time_diff.setter
    def time_diff(self, value):
        self._time_diff = value
    
    @property
    def mask(self):
        if not None in (self.time_diff, self.time_threshold):
            return (self._pixel_mask +
                    (abs(self.time_diff) > self.time_threshold))
        return self._pixel_mask
    
    def write(self, filename):
        """
        Write mapper to hdf5 file *filename*.
        
        """
        import h5py
        
        with h5py.File(filename, 'w') as f:
            f.create_dataset('rows', data=self.rows.filled())
            f.create_dataset('cols', data=self.cols.filled())
            f.create_dataset('pixel_mask', data=self._pixel_mask)
            if self.time_diff is not None:
                f.create_dataset('time_diff', data=self.time_diff.filled())
            if self.time_threshold is not None:
                f.attrs['time_threshold'] = self.time_threshold
    
    @classmethod
    def from_file(cls, filename):
        """
        Create a mapper from contents of *filename*.
        
        """
        import h5py
        
        with h5py.File(filename, 'r') as f:
            rows = f['rows'][:]
            cols = f['cols'][:]
            pixel_mask = f['pixel_mask'][:]
            time_diff = None
            time_threshold = None
            try:
                time_diff = f['time_diff'][:]
            except KeyError:
                pass
            try:
                time_threshold = f.attrs['time_threshold']
            except KeyError:
                pass
        
        return cls(rows=rows, cols=cols, pixel_mask=pixel_mask,
                   time_diff=time_diff, time_threshold=time_threshold)


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
    
    rows = indices // source_def.shape[1]
    cols = indices % source_def.shape[1]
    # Make sure all indices are valid
    rows[rows >= source_def.shape[0]] = -1
    cols[cols >= source_def.shape[1]] = -1
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


def match(amsr_filename, avhrr_filename, radius_of_influence=1e3,
          time_threshold=None):
    """
    Find matching indices in AVHRR array for each element in AMSR swath.
    
    Arguments:
    
        amsr_filename: string
            full path of AMSR-E HDF5 file
        avhrr_filename: string
            full path of AVHRR PPS HDF5 file
        radius_of_influence: float
            radius of influence in meters in pixel-pixel matching (default: 1000 m)
        time_threshold: float
            largest absolute time difference to include in match
    
    Returns:
    
        mapper: `MatchMapper` instance.
    
    """
    avhrr_lonlat = get_avhrr_lonlat(avhrr_filename)
    amsr_lonlat = get_amsr_lonlat(amsr_filename)
    
    mapper = match_lonlat(avhrr_lonlat, amsr_lonlat, radius_of_influence)
    
    avhrr_time = get_avhrr_time(avhrr_filename)
    amsr_time = get_amsr_time(amsr_filename)
    
    time_diff = np.abs(avhrr_time[mapper.rows] - amsr_time.reshape((amsr_time.size, 1)))
    
    mapper.time_diff = time_diff
    mapper.time_threshold = time_threshold
    
    logger.debug("Time diff (min, max): %r" % ((time_diff.min(),
                                                time_diff.max()),))
    
    return mapper


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


def plot_array(lon, lat, data, legend):
    """
    Make a plot of *data*
    
    """
    import matplotlib.pyplot as pl
    from mpl_toolkits.basemap import Basemap
    
    lon, lon0 = adjust_lon(lon)
    
    fig = pl.figure()
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
    
    return fig


class Field(object):
    def __init__(self, data, lon, lat, desc='Field'):
        self.data = data
        self.lon = lon
        self.lat = lat
        self.desc = desc
    
def plot_fields(fields, break_value=170):
    """
    Plot *fields*. Each element in *fields* should be a `Field` instance.
    
    Plots on orthogonal projection, with lon_0, lat_0 taken from mean of lon,
    lat of first element in *fields*.
    
    """
    from mpl_toolkits.basemap import Basemap
    from matplotlib import pyplot as pl
    from utility_functions import broken_cmap
    
    lon_0 = fields[0].lon.mean()
    lat_0 = fields[0].lat.mean()
    
    vmin = 0
    vmax = min(f.data.max() for f in  fields)
    if vmax < vmin:
        vmax = max(break_value, *(f.data.max() for f in  fields))
    cmap = broken_cmap(np.array([vmin, vmax]), break_value=break_value)
    
    fig = pl.figure()
    ax = None
    for ix, f in enumerate(fields):
        ax = fig.add_subplot(len(fields) % 2 + 1, len(fields) // 2 + 1, ix + 1,
                             sharex=ax, sharey=ax)
        m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, lat_ts=0,
                    resolution='c', ax=ax)
        m.drawcoastlines(linewidth=.5, color='g')
        
        step = (f.data.shape[1] // 256) or 1 # don't use full gigantic arrays
        _slice = (slice(None, None, step),) * 2
        x, y = m(f.lon[_slice], f.lat[_slice])
        # Mask out pixels outside projection limb
        #on_map = (x < 1e20) | (y < 1e20)
        #data = np.ma.array(f.data[_slice], mask=~on_map)
        data = f.data[_slice]
        logger.debug("data.shape = %r" % (data.shape,))
        mesh = m.pcolor(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap)
        # pcolormesh is much faster, but I can't get rid of off-projection drawing
        # mesh = m.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap)
        fig.colorbar(mesh)
        ax.set_title(f.desc)
    
    return fig


def get_amsr_lwp(filename):
    """
    Return liquid water path (lwp) from AMSR-E file *filename*. The units of lwp
    are converted from mm to g m**-2.
    
    """
    from ppshdf_cloudproducts import PpsProduct
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
    from ppshdf_cloudproducts import CppProducts
    
    cpp = CppProducts(filename=filename, product_names=['cwp'])
    return cpp.products['cwp'].array


def imshow_lwps(amsr_lwp, cpp_lwp, time_diff, sea, title=None):
    """
    Show *amsr_lwp* and *cpp_lwp* side by side. *sea* is used to draw a
    background, and mask out any pixels which are not sea.
    
    """
    import matplotlib.pyplot as pl
    
    from utility_functions import broken_cmap
    from matplotlib.colors import ListedColormap
    
    fig = pl.figure()
    
    break_value = 170 # g m**-2
    vmin = 0
    vmax = min(amsr_lwp.max(), cpp_lwp.max())
    if vmax < vmin:
        vmax = max(break_value, amsr_lwp.max(), cpp_lwp.max())
    cmap = broken_cmap(np.array([vmin, vmax]), break_value=break_value)
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
    
    return fig


def find_amsr(avhrr_filename):
    """
    Find AMSR-E files matching *avhrr_filename*. Returns a list of file paths.
    
    """
    from file_finders import AmsrFileFinder, PpsFileFinder
    pps_finder = PpsFileFinder()
    parsed = pps_finder.parse(avhrr_filename)
    
    amsr_finder = AmsrFileFinder()
    return amsr_finder.find(parsed['datetime'])


def _test():
    import sys
    from pps_runutil import get_ppsProductArguments
    from pps_ioproxy import pps_ioproxy
    from pps_basic_configure import AVHRR_DIR, OUTPUT_DIR
    
    
    # Command line argument handling
    if sys.argv[1] == 'satproj':
        amsr_filenames = None
    else:
        amsr_filenames = [sys.argv.pop(1)]
    ppsarg, arealist = get_ppsProductArguments(sys.argv) #@UnusedVariable
    ioproxy = pps_ioproxy(ppsarg)
    avhrr_filename = os.path.join(AVHRR_DIR, ppsarg.files.avhrr)
    if not amsr_filenames:
        amsr_filenames = find_amsr(avhrr_filename)
        logger.debug("Found AMSR-E files: %r" % amsr_filenames)
    
    
    for amsr_filename in amsr_filenames:
        match_file = "match--%s--%s.h5" % (os.path.basename(avhrr_filename),
                                         os.path.basename(amsr_filename))
        if os.path.exists(match_file):
            logger.info("Reading match from %r" % match_file)
            mapper = MatchMapper.from_file(match_file)
        else:
            logger.info("Matching AMSR-E and AVHRR swaths")
            mapper = match(amsr_filename, avhrr_filename)
            mapper.write(match_file)
            logger.info("Match written to %r" % match_file)
        
        lon, lat = get_amsr_lonlat(amsr_filename)
        title = "%r\nmatched with\n%r" % (os.path.basename(avhrr_filename),
                                          os.path.basename(amsr_filename))
        fig = plot_array(lon, lat, mapper.time_diff,
                         legend="time difference (s)")
        fig_base = match_file.rsplit('.h5', 1)[0] + '--'
        fig.set_size_inches(20, 12)
        fig.suptitle(title)
        fig.savefig(fig_base + "time_diff.png")
        
        # Compare liquid water paths
        amsr_lwp = get_amsr_lwp(amsr_filename)
        #amsr_lwp[amsr_lwp < 0] = -1
        
        cpp_filename = os.path.join(OUTPUT_DIR, ppsarg.files.cpp)
        if False:
            logger.warning("Replacing '.h5' extension with '.hdf' in CPP filename, "
                           "due to a bug in CPP.")
            cpp_filename = '.hdf'.join(cpp_filename.rsplit('.h5', 1)) # last '.h5'
        if os.path.exists(cpp_filename):
            cpp_lwp_avhrr_proj = get_cpp_lwp(cpp_filename)
            cpp_lwp = mapper(cpp_lwp_avhrr_proj)
            
            landuse = ioproxy.getLanduseOnly()
            # Sea pixels in AMSR-E swath
            sea = mapper(landuse == 16)
            
            fig = imshow_lwps(amsr_lwp, cpp_lwp, mapper.time_diff, sea)
            fig.suptitle(title)
            fig.set_size_inches(20, 12)
            fig.savefig(fig_base + "lwp_arrays.png")
            
            # Plot fields in individual plots
            avhrr_lonlat = get_avhrr_lonlat(avhrr_filename)
            fields = [Field(cpp_lwp_avhrr_proj, desc='CPP cwp', *avhrr_lonlat),
                      Field(amsr_lwp, lon=lon, lat=lat, desc='AMSR-E lwp')]
            fig = plot_fields(fields)
            fig.suptitle(title)
            fig.set_size_inches(20, 10)
            fig.savefig(fig_base + "lwp_swaths.png")
            
            validate_lwp(amsr_lwp, cpp_lwp, sea)
        else:
            logger.warning("No CPP product found")
    
#    import matplotlib.pyplot as pl
#    pl.show()


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
