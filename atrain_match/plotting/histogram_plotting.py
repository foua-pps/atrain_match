"""
Plotting functions for AMSR-E IMAGER matching / validation

"""

import numpy as np
import logging
logger = logging.getLogger(__name__)

import matplotlib 
import matplotlib.pyplot as plt
try:
    matplotlib.use('agg', warn=False)
except TypeError: # For some reason, I don't always have 'warn'
    matplotlib.use('agg')


def imshow_lwps(amsr_lwp, cpp_lwp, time_diff, sea, title=None, lwp_max=None):
    """
    Show *amsr_lwp* and *cpp_lwp* side by side. *sea* is used to draw a
    background, and mask out any pixels which are not sea.
    
    """
    import matplotlib.pyplot as plt
    
    from utility_functions import broken_cmap_r
    from matplotlib.colors import ListedColormap

    #Comment: utility_functions is a separate package, created by J.Malm,
    #but which we can not find. If we want to run this part of the code, it
    #might be changed to something from matplotlib instead.
    #By calling match_util_match without '-p', this part of the code is not
    #executed.
    #Sara Hornquist 2015-03-12

    
    # Use average of all IMAGER pixels in AMSR footprint
    if len(cpp_lwp) == 3:
        cpp_lwp = cpp_lwp.mean(axis=-1)
    
    fig = plt.figure()
    
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
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    
    lon, lon0 = adjust_lon(lon)
    
    fig = plt.figure()
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
    from matplotlib import pyplot as plt
    from utility_functions import broken_cmap_r
    
    #Comment: utility_functions is a separate package, created by J.Malm,
    #but which we can not find. If we want to run this part of the code, it
    #might be changed to something from matplotlib instead.
    #broken_cmap_r is supposed to create a colormap with a gradient on one side
    #of a threshold, and another gradient on the other side, and a unique
    #colour on the theshold it self.
    #Sara Hornquist 2015-03-12

    
    lon_0 = fields[0].lon.mean()
    lat_0 = fields[0].lat.mean()
    
    vmin, vmax, break_value = limits([f.data for f in fields], break_value)
    cmap = broken_cmap_r(vmin=vmin, vmax=vmax, break_value=break_value)
    
    fig = plt.figure()
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
    from matplotlib import pyplot as plt
    fig = plt.figure()
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
    from reshaped_files_plotting.plot_ctth_print_bias_and_std_stats import half_sample_mode, my_iqr
    mode = half_sample_mode(data)
    iqr = my_iqr(data)
    q1 = np.percentile(data,25)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, bars = ax.hist(data, **kwargs) #@UnusedVariable
    #plt.hist(data, **kwargs)
    ax.set_ylabel('frequency')
    ax.axvline(mean, label='mean = %.2f' % mean, color='r')
    ax.axvline(median, label='median = %.2f' % median, color='r',
               linestyle='--')
    ax.axvline(mode, label='mode = %.2f' % mode, color='r',
               linestyle=':')
    ax.hlines(n.max() / 2, mean - std, mean + std,
              label='std = %.2f' % std, colors='r', linestyles='dashdot')
    ax.hlines(n.max() / 4, q1, q1 + iqr,
              label='iqr = %.2f' % iqr, colors='r', linestyle=':')
    ax.grid()
    ax.legend()
    print(np.median(data), np.percentile(data,25), np.percentile(data,75))
    
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

def atrain_scatter(fig, ax, x, y, binsize, xymin = None,  xymax=None, vmax=250, 
                   do_colorbar=True, to_km=1.0, ptype='scatter'):
    import copy

    if xymax is None:
        xymax = np.max([np.max(x),np.max(y)])
    if xymin is None:
        xymin = np.min([np.min(x),np.min(y)])
    n_edges = int((xymax-xymin)*1.0/binsize)+1
    edgesx= np.linspace(xymin, xymax, n_edges)
    edgesy= np.linspace(xymin, xymax, n_edges)
    H, xe, ye = np.histogram2d(x, y, bins=[edgesx,edgesy])
    xi = np.searchsorted(edgesx, x)# - edgesx[0])/(n_edges+1)).astype(np.int)               
    yi = np.searchsorted(edgesy, y)#np.floor((y - edgesy[0])/(n_edges+1)).astype(np.int)  #-1?
    xi = xi-1
    yi = yi-1
    #################################
    #Move pixels out side in side ?:
    yi[yi==n_edges] = n_edges-1 
    xi[xi==n_edges] = n_edges-1
        #yi(yi==2)==1 This is what we need!?
    yi[yi==n_edges-1] = n_edges-2
    xi[xi==n_edges-1] = n_edges-2
    #################################
    z = H[xi,yi] 
    idx = z.argsort()
    my_cmap = copy.copy(matplotlib.cm.get_cmap("inferno_r", lut=100))
    cmap_vals = my_cmap(np.arange(100)) #extractvalues as an array
    #print cmap_vals[0]
    cmap_vals[0:5] = cmap_vals[5] #change the first values to less white
    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "new_inferno_r", cmap_vals) 

    if ptype in ['hexbin']:
        b =plt.hexbin(to_km*x,to_km*y, gridsize=binsize,  cmap=my_cmap,  
                      vmin=10, vmax=vmax)
    if ptype in ['hist2d']:
        b =plt.hist2d(to_km*x,to_km*y,bins=binsize,  cmap=my_cmap,  
                      vmin=1, vmax=vmax)
    if ptype in ['scatter']:
        b = plt.scatter(to_km*x[idx], to_km*y[idx], c=z[idx], 
                        edgecolor='', cmap=my_cmap, vmin=1, vmax=vmax, 
                        alpha=1.0, marker='.', edgecolors=None,  rasterized=True)

    ax.set_ylim(xymin,to_km*xymax)
    ax.set_xlim(xymin,to_km*xymax) 
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #if 'pressure' in truth:        
    #    plt.gca().invert_xaxis()
    #    plt.gca().invert_yaxis()
    #    plt.xticks([950, 550, 150])
    #    plt.yticks([950, 550, 150])
    #else:  
    #    plt.xticks([5, 10, 15])
    #    plt.yticks([0, 5, 10, 15])
    if do_colorbar and ptype in ['scatter']:
        cax = fig.add_axes([0.80, 0.65, 0.06, 0.22])
        cbar = fig.colorbar(b, cax=cax, )
    return b    


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
