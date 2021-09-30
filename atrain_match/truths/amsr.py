# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
import logging
from atrain_match.utils.validate_lwp_util import LWP_THRESHOLD
from atrain_match.libs.extract_imager_along_track import imager_track_from_matched
from atrain_match.utils.common import (ProcessingError, elements_within_range)
from atrain_match.utils.runutils import do_some_logging
from atrain_match.truths.calipso import find_break_points
from atrain_match.matchobject_io import (TruthImagerTrackObject,
                                         AmsrObject)
import atrain_match.config as config
import h5py
import numpy as np
from datetime import datetime
from calendar import timegm

TAI93 = datetime(1993, 1, 1)

logger = logging.getLogger(__name__)


def calculate_amsr_rof(overlap, imager_radius):
    """ Calculate radius of influence for nearest neighbour search. 
        
        100% overlap:  fov_semi_minor_axis_diam / 2 -IMAGER_RADIUS
        50% overlap :  fov_semi_minor_axis_diam / 2 
        25% overlap:  fov_semi_minor_axis_diam / 2  + 0.5*IMAGER_RADIUS
    """
    fov_semi_minor_axis_diam = 7000.
    rof = 0.5*fov_semi_minor_axis_diam + imager_radius - 2*overlap*imager_radius
    return rof


def get_amsr_rof(SETTINGS, imager):
    """ Get radius of influence as a function of overlap and imager pixel size. """
    if SETTINGS['AMSR_IMAGER_RADIUS'] < 0:
        if imager == 'seviri':
            # pixel size 4.8 x 3 km
            imager_radius = 5660.
        elif imager == 'avhrr':
            # pixel size 4.3 x 3 km
            imager_radius = 5325
        elif imager == 'viirs':
            pass
        elif imager == 'modis':
            pass
    return calculate_amsr_rof(SETTINGS['AMSR_OVERLAP'], imager_radius)


def get_amsr(filename, SETTINGS):

    if SETTINGS['AMSR_SENSOR'].lower() in ['amsre', 'amsr-e']:
        AMSR_SENSOR = 'AMSR-E'
        lwp_conversion = 1E3 # Density of water [kg m**-3]
    elif SETTINGS['AMSR_SENSOR'].lower() == 'amsr2_jaxa':
        AMSR_SENSOR = 'AMSR2'
        AMSR2_DATA_SOURCE = 'JAXA'
        lwp_conversion = 1E3 # kg m^-2 to g m^-2
    elif SETTINGS['AMSR_SENSOR'].lower() == 'amsr2_nsidc':
        AMSR_SENSOR = 'AMSR2'
        AMSR2_DATA_SOURCE = 'NSIDC'
        lwp_conversion = 1E3 # kg m^-2 to g m^-2
    else:
        raise Exception('Please specifiy AMSR_SENSOR in the config ' +\
                        ' file. [AMSR2_JAXA, AMSR2_NSIDC, AMSR-E]')

    if AMSR_SENSOR == 'AMSR-E':
        if ".h5" in filename:
            retv = read_amsre_h5(filename)
        else:
            # hdf4 file:
            retv = read_amsre_hdf4(filename)
    elif AMSR_SENSOR == 'AMSR2':
        if AMSR2_DATA_SOURCE == 'JAXA':
            # HDF5 format
            retv = read_amsr2_h5(filename)
        else:
            # HDF-EOS5 format
            retv = read_amsr2_he5(filename)

    n_lat_scans = len(retv.latitude) * 1.0 / (len(retv.sec1993))  # = 242!
    epoch_diff = timegm(TAI93.utctimetuple())
    nadir_sec_1970 = retv.sec1993 + epoch_diff
    retv.sec_1970 = np.repeat(nadir_sec_1970.ravel(), n_lat_scans)
    retv.sec1993 = None
    # convert to g m^-2
    retv.lwp = retv.lwp_mm.ravel() * lwp_conversion

    logger.info("Extract {} lwp between 0 and {} g/m-2".format(AMSR_SENSOR,
                                                               LWP_THRESHOLD))
    use_amsr = np.logical_and(retv.lwp >= 0,
                              retv.lwp < LWP_THRESHOLD)
    retv = retv.extract_elements(idx=use_amsr)
    return retv


def read_amsre_h5(filename):
    retv = AmsrObject()

    with h5py.File(filename, 'r') as f:
        # ravel AMSR-E data to 1 dimension
        retv.longitude = f['Swath1/Geolocation Fields/Longitude'][:].ravel()
        retv.latitude = f['Swath1/Geolocation Fields/Latitude'][:].ravel()
        retv.sec1993 = f['Swath1/Geolocation Fields/Time']['Time'][:]
        # description='lwp (mm)',
        lwp_gain = f['Swath1/Data Fields/High_res_cloud'].attrs['Scale']  # .ravel()
        retv.lwp_mm = f['Swath1/Data Fields/High_res_cloud'][:].ravel() * lwp_gain
    if f:
        f.close()
    return retv
  
  
def read_amsr2_h5(filename):
    retv = AmsrObject()

    with h5py.File(filename, 'r') as f:
        # ravel AMSR-E data to 1 dimension
        retv.longitude = f['Longitude of Observation Point'][:].ravel()
        retv.latitude = f['Latitude of Observation Point'][:].ravel()
        retv.sec1993 = f['Scan Time'][:]
        # description='lwp (mm)',
        lwp_gain = f['Geophysical Data'].attrs['SCALE FACTOR']  # .ravel()
        retv.lwp_mm = f['Geophysical Data'][:].ravel() * lwp_gain
    if f:
        f.close()
    return retv
  

def read_amsr2_he5(filename):
    retv = AmsrObject()
    with h5py.File(filename, 'r') as f:
        geoloc = f['HDFEOS']['SWATHS']['AMSR2_L1R']['Geolocation Fields']
        data = f['HDFEOS']['SWATHS']['AMSR2_L1R']['Data Fields']
      
        retv.longitude = geoloc['Longitude'][:].ravel()
        retv.latitude = geoloc['Latitude'][:].ravel()
        retv.sec1993 = geoloc['tai93time'][:]
        retv.lwp_mm = data['CloudWaterPath'][:].ravel()
    if f:
        f.close()
    return retv
  

def read_amsre_hdf4(filename):
    from pyhdf.SD import SD, SDC
    from pyhdf.HDF import HDF  # HC
    import pyhdf.VS

    retv = AmsrObject()
    h4file = SD(filename, SDC.READ)
    # datasets = h4file.datasets()
    # attributes = h4file.attributes()
    # for idx, attr in enumerate(attributes.keys()):
    #    print idx, attr
    for sds in ["Longitude", "Latitude", "High_res_cloud"]:
        data = h4file.select(sds).get()
        if sds in ["Longitude", "Latitude"]:
            retv.all_arrays[sds.lower()] = data.ravel()
        elif sds in ["High_res_cloud"]:
            lwp_gain = h4file.select(sds).attributes()['Scale']
            retv.all_arrays["lwp_mm"] = data.ravel() * lwp_gain

        # print h4file.select(sds).info()
    h4file = HDF(filename, SDC.READ)
    vs = h4file.vstart()
    data_info_list = vs.vdatainfo()
    # print "1D data compound/Vdata"
    for item in data_info_list:
        # 1D data compound/Vdata
        name = item[0]
        # print name
        if name in ["Time"]:
            data_handle = vs.attach(name)
            data = np.array(data_handle[:])
            retv.all_arrays["sec1993"] = data
            data_handle.detach()
        else:
            pass
            # print name
        # data = np.array(data_handle[:])
        # attrinfo_dic = data_handle.attrinfo()
        # factor = data_handle.findattr('factor')
        # offset = data_handle.findattr('offset')
        # print data_handle.factor
        # data_handle.detach()
    # print data_handle.attrinfo()
    h4file.close()
    # for key in retv.all_arrays.keys():
    #    print key, retv.all_arrays[key]
    return retv


def reshape_amsr(amsrfiles, imager, SETTINGS):
    amsr = get_amsr(amsrfiles[0])
    for i in range(len(amsrfiles) - 1):
        new_amsr = get_amsr(amsrfiles[i + 1])
        amsr_start_all = amsr.sec_1970.ravel()
        amsr_new_all = new_amsr.sec_1970.ravel()
        if not amsr_start_all[0] < amsr_new_all[0]:
            raise ProcessingError("AMSR files are in the wrong order")
        # Concatenate the feature values
        amsr = amsr + new_amsr
    start_break, end_break = find_break_points(amsr, imager, SETTINGS)
    amsr = amsr.extract_elements(starti=start_break,
                                 endi=end_break)
    return amsr


def match_amsr_imager(amsr, cloudproducts, SETTINGS):
    retv = TruthImagerTrackObject(truth='amsr')
    retv.imager_instrument = cloudproducts.instrument.lower()
    retv.amsr = amsr
    if (getattr(cloudproducts.cpp, "cpp_lwp") < 0).all():
        logger.warning("Not matching AMSR-E with scene with no lwp.")
        return None
        # return MatchupError("No imager Lwp.") # if only LWP matching?

    from atrain_match.utils.common import map_imager_distances
    n_neighbours = 8
    if config.RESOLUTION == 5:
        n_neighbours = 5
    
    if SETTINGS['AMSR_RADIUS'] < 0:
        AMSR_RADIUS = get_amsr_rof(SETTINGS, cloudproducts.imager)
    else:
        AMSR_RADIUS = SETTINGS['AMSR_RADIUS']
        
    mapper_and_dist = map_imager_distances(cloudproducts,
                                           amsr.longitude.ravel(),
                                           amsr.latitude.ravel(),
                                           radius_of_influence=AMSR_RADIUS,
                                           n_neighbours=n_neighbours)
    cal, cap = mapper_and_dist["mapper"]
    distances = mapper_and_dist["distances"]
    cal_1 = cal[:, 0]
    cap_1 = cap[:, 0]

    calnan = np.where(cal_1 == config.NODATA, np.nan, cal_1)
    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None
    # check if it is within time limits:
    if len(cloudproducts.time.shape) > 1:
        imager_time_vector = [cloudproducts.time[line, pixel] for line, pixel in zip(cal_1, cap_1)]
        imager_lines_sec_1970 = np.where(cal_1 != config.NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal_1 != config.NODATA, cloudproducts.time[cal_1], np.nan)
    # Find all matching Amsr pixels within +/- sec_timeThr from the IMAGER data
    imager_sunz_vector = np.array([cloudproducts.imager_angles.sunz.data[line, pixel]
                                   for line, pixel in zip(cal_1, cap_1)])
    idx_match = np.logical_and(
        elements_within_range(amsr.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr"]),
        imager_sunz_vector <= 84)  # something larger than 84 (max for lwp)

    if idx_match.sum() == 0:
        logger.warning("No light (sunz<84)  matches in region within time threshold %d s.", SETTINGS["sec_timeThr"])
        return None
    retv.amsr = retv.amsr.extract_elements(idx=idx_match)

    # Amsr line, pixel inside IMAGER swath (one neighbour):
    retv.amsr.imager_linnum = np.repeat(cal_1, idx_match).astype('i')
    retv.amsr.imager_pixnum = np.repeat(cap_1, idx_match).astype('i')
    retv.amsr.imager_linnum_nneigh = np.repeat(cal, idx_match, axis=0)
    retv.amsr.imager_pixnum_nneigh = np.repeat(cap, idx_match, axis=0)
    retv.amsr.imager_amsr_dist = np.repeat(distances, idx_match, axis=0)

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.amsr.sec_1970 - retv.imager.sec_1970

    do_some_logging(retv, amsr)
    logger.debug("Extract imager lwp along track!")

    retv = imager_track_from_matched(retv, SETTINGS,
                                     cloudproducts,
                                     extract_radiances=False,
                                     extract_aux_segments=False,
                                     extract_ctth=False,
                                     extract_ctype=False,
                                     aux_params=['fractionofland', 'landuse'],
                                     extract_some_data_for_x_neighbours=True)
    return retv

# ------------------------------------------------------------------------------
if __name__ == '__main__':

    filename = "/home/a001865/temp/AMSR_E_L2_Ocean_V06_201002011115_A.hdf"
    read_amsr_hdf4(filename)
