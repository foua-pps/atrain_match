# -*- coding: utf-8 -*-

"""
Utilities

"""

import numpy as np
import logging
logger = logging.getLogger(__name__)
from read_cloudproducts_cci import cci_read_prod
import h5py
from datetime import datetime
TAI93 = datetime(1993, 1, 1)

#: h5py compression settings (True, or an integer in range(10))
_COMPRESSION = True


def get_amsr_lonlat(filename):
    """
    Get (lon, lat) from AMSR file *filename*.
    
    """
    
    with h5py.File(filename, 'r') as f:
        lon = f['Swath1/Geolocation Fields/Longitude'][:]
        lat = f['Swath1/Geolocation Fields/Latitude'][:]
        if lon.shape[1]>0:
            print "Warning, expected 1D vectors"
    if f:
        f.close()
    return lon, lat


def get_avhrr_lonlat(filename):
    """
    Get (lon, lat) from AVHRR file *filename* (remapped by PPS).
    
    """
    import pps_io
    
    geo = pps_io.readAvhrrGeoData(filename)
    
    return geo.longitude, geo.latitude

def get_avhrr_lonlat_and_time_cci(filename):
    """
    Get (lon, lat) from CCI AVHRR file .
    
    """    
    geo = cci_read_prod(filename, 'geotime')
    n_scanlines = geo.longitude.shape[0]
    sec1970 = np.linspace(geo.sec1970_start, geo.sec1970_end, n_scanlines)
    return (geo.longitude, geo.latitude), sec1970


def get_amsr_time(filename):
    """
    Get time of each scanline in AMSR file *filename*.
    
    Note::
    
        Time is converted from seconds since 1993-01-01 to seconds since
        1970-01-01.
    
    """
    with h5py.File(filename, 'r') as f:
        sec1993 = f['Swath1/Geolocation Fields/Time']['Time'][:]
    if f:
        f.close()
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
    
    with h5py.File(filename, 'r') as f:
        lwp_mm = PpsProduct(f['Swath1/Data Fields/High_res_cloud'][:],
                            description='lwp (mm)',
                            gain=f['Swath1/Data Fields/'
                                   'High_res_cloud'].attrs['Scale'],
                            nodata=-9990, scale_up=True)
    if f:
        f.close()
    
    density = 1e3 # Density of water [kg m**-3]
    return lwp_mm.array * density # [mm * kg m**-3 = g m**-2]


def get_cpp_product(filename, product):
    """Get *product* from CPP file *filename*."""
    from ppshdf_cloudproducts import CppProducts
    
    try:
        cpp = CppProducts.from_h5(filename, product_names=[product])
    except TypeError:
        logger.warning("Running with old ACPG (< r2_45) -- using old "
                       "CppProducts interface")
        cpp = CppProducts(filename=filename, product_names=[product])
    return cpp.products[product].array

def get_cpp_product_cci(filename, product):
    """Get *product* from CPP file *filename*."""
    cpp_cph = cci_read_prod(filename, 'phase')
    return cpp_cph

def get_diff_data(filenames, aux_fields=None):
    """
    Read data from diff files *filenames* and return concatenated lwp_diff.
    
    Dataset names in *aux_fields* will also be extracted, concatenated and
    returned.
    
    Returns:
    
    (lwp_diff, restrictions, *aux_fields)
    
    """
    lwp_diffs = []
    aux = {}
    if aux_fields is not None:
        for name in aux_fields:
            aux[name] = []
    
    for filename in filenames:
        with h5py.File(filename, 'r') as f:
            data = f['lwp_diff'][:]
            restrictions = f['lwp_diff'].attrs['restrictions']
            if len(lwp_diffs) > 0:
                if not (restrictions == lwp_diffs[-1][-1]).all():
                    raise RuntimeError(
                        "Inconsistent restrictions: %r != %r" % (restrictions, lwp_diffs[-1][-1]))
            lwp_diffs.append((filename, data, restrictions))
            for name in aux_fields:
                aux[name].append(f[name][:])
        if f:
            f.close()


    
    lwp_diff_array = np.concatenate([data for 
            (filename, data, restrictions) in lwp_diffs])
    for name in aux_fields:
        aux[name] = np.concatenate(aux[name])
    
    return [lwp_diff_array, restrictions] + [aux[name] for name in aux_fields]


def write_data(data, name, filename, mode='a', attributes=None):
    """
    Write *data* and any *attributes* (dict) to dataset *name* in HDF5 file
    *filename*. *mode* is the h5py file access mode (default 'a', for append).
    
    """
    with h5py.File(filename, mode) as f:
        if hasattr(data, 'mask'):
            data = data.compressed()
        d = f.create_dataset(name, data=data, compression=_COMPRESSION)
        if attributes:
            for k, v in attributes.items():
                d.attrs[k] = v
    #if f:
    #    f.close()

def write_data_to_open_file(data, name, filehandle,  attributes=None):
    """
    Write *data* and any *attributes* (dict) to dataset *name* in HDF5 file
    *filename*. *mode* is the h5py file access mode (default 'a', for append).
    
    """
    if hasattr(data, 'mask'):
        data = data.compressed()
    d = filehandle.create_dataset(name, data=data, compression=_COMPRESSION)
    if attributes:
        for k, v in attributes.items():
            d.attrs[k] = v


def diff_file_stats(filename):
    """
    Print some statistics for diff file *filename*.
    
    """
    from runutils import parse_scene
    import sys
    write = sys.stdout.write
    write('Scene: %s %s %d' % parse_scene(filename.replace('match--', '')))
    
    lwp_diff, restrictions, lat = get_diff_data([filename], ['latitudes'])
    
    #print("Restrictions: %s" % '; '.join(restrictions))
    write(' ¤ Valid pixels: % 6d' % lwp_diff.size)
    write(' ¤ Latitudes: %5.2f - %5.2f' % (lat.min(), lat.max()))
    print('')
