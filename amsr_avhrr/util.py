"""
Utilities

"""

import numpy as np

from datetime import datetime
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
    
    density = 1e3 # Density of water [kg m**-3]
    return lwp_mm.array * density # [mm * kg m**-3 = g m**-2]


def get_cpp_product(filename, product):
    """Get *product* from CPP file *filename*."""
    from ppshdf_cloudproducts import CppProducts
    
    cpp = CppProducts(filename=filename, product_names=[product])
    return cpp.products[product].array


def get_diff_data(filenames):
    """
    Read data from diff files *filenames* and return concatenated arrays
    (lwp_diff, cwp, lwp).
    
    """
    import h5py
    lwp_diffs = []
    cwps = []
    lwps = []
    for filename in filenames:
        with h5py.File(filename, 'r') as f:
            data = f['lwp_diff'][:]
            restrictions = f['lwp_diff'].attrs['restrictions']
            if len(lwp_diffs) > 0:
                if not (restrictions == lwp_diffs[-1][-1]).all():
                    raise RuntimeError(
                        "Inconsistent restrictions: %r != %r" % (restrictions, lwp_diffs[-1][-1]))
            lwp_diffs.append((filename, data, restrictions))
            if 'amsr_lwp' in f.keys():
                lwps.append(f['amsr_lwp'][:])
            if 'cpp_cwp' in f.keys():
                cwps.append(f['cpp_cwp'][:])
    
    lwp_diff_array = np.concatenate([data for 
            (filename, data, restrictions) in lwp_diffs])
    
    if not 0 in (len(cwps), len(lwps)):
        cwp_array = np.concatenate(cwps)
        lwp_array = np.concatenate(lwps)
    else:
        cwp_array = None
        lwp_array = None
    
    return lwp_diff_array, restrictions, cwp_array, lwp_array
