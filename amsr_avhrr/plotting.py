"""
Plotting functions for AMSR-E AVHRR matching / validation

"""

import numpy as np
import logging
logger = logging.getLogger(__name__)

import matplotlib
try:
    matplotlib.use('agg', warn=False)
except TypeError: # For some reason, I don't always have 'warn'
    matplotlib.use('agg')


def imshow_lwps(amsr_lwp, cpp_lwp, time_diff, sea, title=None, lwp_max=None):
    """
    Show *amsr_lwp* and *cpp_lwp* side by side. *sea* is used to draw a
    background, and mask out any pixels which are not sea.
    
    """
    import matplotlib.pyplot as pl
    
    from utility_functions import broken_cmap_r
    from matplotlib.colors import ListedColormap
    
    # Use average of all AVHRR pixels in AMSR footprint
    if len(cpp_lwp) == 3:
        cpp_lwp = cpp_lwp.mean(axis=-1)
    
    fig = pl.figure()
    
    vmin, vmax, break_value = limits([amsr_lwp, cpp_lwp], lwp_max)
    cmap = broken_cmap_r(np.array([vmin, vmax]), break_value=break_value)
    ground_sea_cmap = ListedColormap(['g', 'b'], name="ground/sea map")
    sea_map = np.ma.array(sea, mask=sea.mask + sea)
    
    ax = fig.add_subplot(131)
    ax.imshow(sea_map, cmap=ground_sea_cmap)
    im = ax.imshow(np.ma.array(amsr_lwp, mask=~sea),
                   vmin=vmin, vmax=vmax, cmap=cmap)
    cbar = fig.colorbar(im)
    cbar.set_label('g m**-2')
    ax.set_title('AMSR-E lwp')
    
    ax = fig.add_subplot(132, sharex=ax, sharey=ax)
    ax.imshow(sea_map, cmap=ground_sea_cmap)
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


class Field(object):
    def __init__(self, data, lon, lat, desc='Field'):
        self.data = data
        self.lon = lon
        self.lat = lat
        self.desc = desc


def limits(arrays, break_value=None):
    """
    Returns data limits for color scaling (vmin, vmax, break_value).
    
    vmax is the smallest of the maxima::
    
    >>> limits([np.arange(13), np.arange(6, 20)], 3)
    (0, 12, 3)
    
    If *break_value* is not given, vmax is returned::
    
    >>> limits([np.arange(17)])
    (0, 16, 16)
    
    vmin is always set to 0::
    
    >>> limits([np.arange(4, 10)])
    (0, 9, 9)
    
    break_values > vmax are retained::
    
    >>> limits([np.arange(4, 10)], 20)
    (0, 9, 20)
    
    When all maxima are below 0, vmax is set to *break_value*::
    
    >>> limits([np.arange(-10, -4)], 8)
    (0, 8, 8)
    
    """
    vmin = 0
    vmax = min(a.max() for a in arrays)
    if break_value is None:
        break_value = vmax
    if vmax < vmin:
        vmax = max(break_value, *(a.max() for a in arrays))
    
    return vmin, vmax, break_value


def plot_fields(fields, break_value=None):
    """
    Plot *fields*. Each element in *fields* should be a `Field` instance.
    
    Plots on orthogonal projection, with lon_0, lat_0 taken from mean of lon,
    lat of first element in *fields*.
    
    """
    from mpl_toolkits.basemap import Basemap
    from matplotlib import pyplot as pl
    from utility_functions import broken_cmap_r
    
    lon_0 = fields[0].lon.mean()
    lat_0 = fields[0].lat.mean()
    
    vmin, vmax, break_value = limits([f.data for f in fields], break_value)
    cmap = broken_cmap_r(vmin=vmin, vmax=vmax, break_value=break_value)
    
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
        #mesh = m.pcolor(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap)
        # pcolormesh is much faster, but I can't get rid of off-projection drawing
        mesh = m.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap)
        fig.colorbar(mesh)
        ax.set_title(f.desc)
        
        m.drawmeridians(range(0, 360, 20), linewidth=.5)
        m.drawparallels(range(-80, 90, 10), linewidth=.5)
        # Mark 70 deg latitudes (cut off for validation)
        m.drawparallels([-70, 70], color='r', dashes=[1, 0], latmax=70)
    
    return fig


def plot_hists(fields, bins=500):
    """
    Plot histograms of *fields* (list of (data, label) tuples).
    
    """
    from matplotlib import pyplot as pl
    fig = pl.figure()
    ax = fig.add_subplot(111)
    
    for data, label in fields:
        values, _bins = np.histogram(data, bins=bins)
        ax.plot((_bins[:-1] + _bins[1:]) / 2., values,
                label=label + ' %d' % values.sum())
    
    ax.legend()
    return fig


def plot_hist(data, **kwargs):
    from matplotlib import pyplot as plt
    
    mean = data.mean()
    median = np.median(data)
    std = data.std()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, bars = ax.hist(data, **kwargs) #@UnusedVariable
    ax.set_ylabel('frequency')
    ax.axvline(mean, label='mean = %.2f' % mean, color='r')
    ax.axvline(median, label='median = %.2f' % median, color='r',
               linestyle='--')
    ax.hlines(n.max() / 2, mean - std, mean + std,
              label='std = %.2f' % std, colors='r', linestyles='dashdot')
    ax.grid()
    ax.legend()
    
    return fig


def density(x, y, bins=None, log=True):
    """
    Create a 2d density plot of correlated values in *x* and *y*.
    
    """
    H, xedges, yedges = np.histogram2d(y, x, bins=bins) #@UnusedVariable
    
    from matplotlib import pyplot as plt
    from matplotlib.colors import LogNorm
    fig = plt.figure()
    ax = fig.add_subplot(111)
    norm = LogNorm(vmin=1, vmax=10**np.ceil(np.log10(H.max()))) if log else None
    ax.plot((0, min(H.shape)), (0, min(H.shape)), 'k', linewidth=.5, alpha=.5)
    im = ax.imshow(H, norm=norm, origin='lower')
    fig.colorbar(im)
    
    return fig


def distribution_map(lon, lat):
    """
    Create a plot of the distribution of pixels on a world map.
    
    """
    from mpl_toolkits.basemap import Basemap
    from matplotlib import pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-90,
                urcrnrlon=180, urcrnrlat=90, resolution='c', ax=ax)
    m.drawcoastlines(linewidth=.5, color='grey')
    meridians = np.arange(-180, 181, 30)
    parallels = np.arange(-90, 91, 30)
    m.drawmeridians(meridians, linewidth=.5, color='grey')
    m.drawparallels(parallels, linewidth=.5, color='grey')
    ax.set_xticks(meridians[::2])
    ax.set_yticks(parallels[::2])
    
    H, LATS, LONS = np.histogram2d(lat.ravel(), lon.ravel(),
                                   bins=[np.arange(-90, 90, 1),
                                         np.arange(-180, 180, 1)])
    im = m.pcolor(LONS, LATS, np.ma.masked_equal(H, 0), edgecolors='none')
    cbar = fig.colorbar(im, orientation='horizontal')
    cbar.set_label('Frequency')
    
    return fig


def basemap2area_def(m, x_size, y_size, **kwargs):
    """Get AreaDefinition instance from Basemap *m*
    
    :Parameters:
    m : `mpl_toolkits.basemap.Basemap` instance
        geometry.AreaDefinition object
    **kwargs: Keyword arguments
        Additional initialization arguments for AreaDefinition
        
    :Returns:
    area_def : `pyresample.geometry.AreaDefinition` instance
    
    """
    from pyresample.geometry import AreaDefinition
    
    d2r = np.pi / 180
    
    llcrnrx, llcrnry = m.projtran(m.llcrnrlon * d2r, m.llcrnrlat * d2r)
    urcrnrx, urcrnry = m.projtran(m.urcrnrlon * d2r, m.urcrnrlat * d2r)
    area_extent = llcrnrx, llcrnry, urcrnrx, urcrnry
    proj_dict = dict(m.projparams)
    proj_dict['a'] = m.rmajor
    proj_dict['b'] = m.rminor
    
    return AreaDefinition('area_id', 'name', 'proj_id', proj_dict, x_size,
                          y_size, area_extent)
