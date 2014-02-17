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


def get_avhrr_sunsat(filename):
    """
    Get (sunz, satz) from AVHRR file *filename* (remapped by PPS).
    
    """
    import pps_io
    sunsat = pps_io.readSunSatAngles(filename)
    sunz = np.where(np.logical_or(sunsat.sunz.data == sunsat.sunz.no_data,
                                  sunsat.sunz.data == sunsat.sunz.missing_data),
                    sunsat.sunz.no_data,
                    sunsat.sunz.data * sunsat.sunz.gain +
                    sunsat.sunz.intercept)
    satz = np.where(np.logical_or(sunsat.satz.data == sunsat.satz.no_data,
                                  sunsat.satz.data == sunsat.satz.missing_data),
                    sunsat.satz.no_data,
                    sunsat.satz.data * sunsat.satz.gain +
                    sunsat.satz.intercept)
    return sunz, satz

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
    from calipso import createAvhrrTime
    from cloudsat_calipso_avhrr_match import (
        get_satid_datetime_orbit_from_fname)
    import time
    geo = pps_io.readAvhrrGeoData(filename)
    values = get_satid_datetime_orbit_from_fname(filename)
    #datetime=values["date_time"]
    geo_ok = createAvhrrTime(geo, values)

    n_scanlines = geo.longitude.shape[0]
    sec1970 = np.linspace(geo_ok.sec1970_start, geo_ok.sec1970_end, n_scanlines)
    tim1 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(geo_ok.sec1970_start))
    tim2 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(geo_ok.sec1970_end))

    logger.info("Starttime avhrr/viirs: %s, end time: %s"%(tim1, tim2))
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
    
    cpp = CppProducts.from_h5(filename, product_names=[product])
    return cpp.products[product].array

def get_cpp_product_cci(filename, product):
    """Get *product* from CPP file *filename*."""
    cpp_cph = cci_read_prod(filename, 'phase')
    return cpp_cph

def reduce_cpp_lwp_data(lwp, sunsat_filename, dcwp, flag=None):
    """Validate for a sub-set of lwp-data"""
    from config import PPS_FORMAT_2012_OR_EARLIER
    from amsr_avhrr_match import NODATA_CPP, NODATA_CPP_v2012
    
    if PPS_FORMAT_2012_OR_EARLIER:
        nodata = NODATA_CPP_v2012
    else:
        nodata = NODATA_CPP
    return reduce_cpp_data(lwp, nodata, sunsat_filename, dcwp, flag=flag)
   
def reduce_cpp_cph_data(cph, sunsat_filename, dcwp, flag=None):
    """Validate for a sub-set of cph-data"""
    from config import PPS_FORMAT_2012_OR_EARLIER
    from validate_cph import CPP_PHASE_VALUES, CPP_PHASE_VALUES_v2012
    
    if PPS_FORMAT_2012_OR_EARLIER:
        nodata = CPP_PHASE_VALUES_v2012['no_observation']
    else:
        nodata = CPP_PHASE_VALUES['no_data']
    return reduce_cpp_data(cph, nodata, sunsat_filename, dcwp, flag=flag)


def reduce_cpp_data(arr, nodata, sunsat_filename, dcwp, flag=None):
    """
    Validate for a sub-set of the pixels for lwp- or cph-data. Configure
    inside this function which limitations to do, any of:
        sun/satellite angles max and/or min
        remove sunglint
        remove pixels where dcwp is over a limit
        """
        
    from amsr_avhrr_match import NODATA_CPP

    print "Make reduction:"
    #Configure which type of reductions to do
    #  'None' means that do not make a reduction on that quantity
    #   example values: min_angle=72, max_angle=78 (or use just one of them)
    #                   dcwp_limit=2 (it is in kg/m2), no_sunglint=1
    max_angle = None
    min_angle = 72
    dcwp_limit = None
    no_sunglint = None


    #Filter out sunglint, if configured for
    if ((no_sunglint is not None) and (int(no_sunglint) == 1)):
        if (flag is None):
            print "no sunglint flag available, can not filter out those"
        else:
            print "Reduction: sunglint removal"
            from pps_basic_util import evaluate_bit
            from pps_common import SM_FLAG_BIT_SUNGLINT

            arr = np.where(
                evaluate_bit(flag.astype('uint16'), SM_FLAG_BIT_SUNGLINT),
                nodata,
                arr)
    
    #If there is a dcwp-limit, use that one
    if (dcwp_limit is not None):
        if (dcwp == None):
            print "no dcwp data available, can not filter on dcwp"
        else:
            print "Reduction: remove dcwp >= ", dcwp_limit
            dcwp_limit = float(dcwp_limit)
            arr = np.where(np.logical_and(dcwp != NODATA_CPP,
                                          dcwp >= dcwp_limit),
                           nodata,
                           arr)

    #Check for reductions in angles (can be both or either of max/min)
    if ((max_angle is None) and (min_angle is None)):
        #No configuration set to reduce arr due to angles
        return arr

    sunz, satz = get_avhrr_sunsat(sunsat_filename)
    if (max_angle is not None):
        print "Reduction: remove sunsat >= ", max_angle
        max_angle = float(max_angle)
        arr = np.where(np.logical_or(sunz >= max_angle,
                                     satz >= max_angle),
                       nodata,
                       arr)
    if (min_angle is not None):
        print "Reduction: remove sunsat < ", min_angle
        min_angle = float(min_angle)
        arr = np.where(np.logical_and(sunz < min_angle,
                                      satz < min_angle),
                       nodata,
                       arr)
    return arr

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
