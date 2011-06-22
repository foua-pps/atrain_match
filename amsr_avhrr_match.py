"""
Match AMSR and AVHRR data

"""

from __future__ import with_statement
import numpy as np
from datetime import datetime
import os

TAI93 = datetime(1993, 1, 1)


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
    
    return (np.ma.array(rows, mask=mask), np.ma.array(cols, mask=mask))


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
    
    rows, cols = match_lonlat(avhrr_lonlat, amsr_lonlat, RESOLUTION * 1e3)
    
    avhrr_time = get_avhrr_time(avhrr_filename)
    amsr_time = get_amsr_time(amsr_filename)
    
    time_diff = np.abs(avhrr_time[rows] - amsr_time.reshape((amsr_time.size, 1)))
    time_diff = np.ma.array(time_diff, mask=rows.mask)
    
    rows.mask[time_diff > sec_timeThr] = True
    if cols.mask is not rows.mask:
        cols.mask[time_diff > sec_timeThr] = True
    
    print("Time diff (min, max): %r" % ((time_diff.min(), time_diff.max()),))
    
    return rows, cols, time_diff


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
    Make a plot of the time difference *time_diff*
    
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
    
    pl.show()


def _test():
    import sys
    amsr_filename = sys.argv[1]
    avhrr_filename = sys.argv[2]
    
    rows, cols, time_diff = match(amsr_filename, avhrr_filename) #@UnusedVariable
    
    lon, lat = get_amsr_lonlat(amsr_filename)
    plot_array(lon, lat, time_diff,
               title="%r\nmatched with\n%r" % (os.path.basename(avhrr_filename),
                                               os.path.basename(amsr_filename)),
               legend="time difference (s)")


if __name__ == '__main__':
    _test()
