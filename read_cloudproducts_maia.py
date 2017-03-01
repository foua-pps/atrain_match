"""
  Use this module to read maia cloudproducts
  2016 SMHI, N.Hakansson a001865
"""

import time
import os
import h5py
import numpy as np
import calendar
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

from read_cloudproducts_and_nwp_pps import (CtypeObj, CtthObj, CmaObj,
                                            createAvhrrTime,
                                            imagerAngObj, imagerGeoObj)
from config import NODATA
ATRAIN_MATCH_NODATA = NODATA
from get_flag_info import get_maia_ct_flag, get_day_night_twilight_info_maia


def get_satid_datetime_orbit_from_fname_maia(avhrr_filename):
    # Get satellite name, time, and orbit number from avhrr_file
    # viiCT_npp_DB_20120817_S035411_E035535_DES_N_La052_Lo-027_00001.h5
    sl_ = os.path.basename(avhrr_filename).split('_')
    date_time = datetime.strptime(sl_[3] + sl_[4], '%Y%m%dS%H%M%S')

    sat_id = sl_[1].lower()
    values = {"satellite": sat_id,
              "date_time": date_time,
              "orbit": "99999",
              "date": date_time.strftime("%Y%m%d"),
              "year": date_time.year,
              "month": "%02d" % (date_time.month),
              "time": date_time.strftime("%H%M"),
              #"basename":sat_id + "_" + date_time.strftime("%Y%m%d_%H%M_99999"),#"20080613002200-ESACCI",
              "ccifilename": avhrr_filename,
              "ppsfilename": None}
    values['basename'] = values["satellite"] + "_" + \
        values["date"] + "_" + values["time"] + "_" + values["orbit"]
    return values


def maia_read_all(filename):
    """Read geolocation, angles info, ctth, and cloudtype
    """
    from runutils import unzip_file

    unzipped = unzip_file(filename)
    if unzipped:
        filename = unzipped

    logger.info("Opening file %s" % (filename))
    with h5py.File(filename, 'r') as maia_h5:
        logger.info("Reading angles ...")
        avhrrAngObj = read_maia_angobj(maia_h5)
        logger.info("Reading cloud type ...")
        # , avhrrAngObj)
        ctype, cma, ctth = read_maia_ctype_cmask_ctth(maia_h5)
        logger.info("Reading longitude, latitude and time ...")
        avhrrGeoObj = read_maia_geoobj(maia_h5, filename)
        logger.info("Reading surface temperature")
        surft = read_maia_surftemp(maia_h5)
        logger.info("Not reading cloud microphysical properties")
        cpp = None
        logger.info("Not reading channel data")
        avhrrObj = None

    if unzipped:
        os.remove(unzipped)

    return avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surft, cpp, cma


def read_maia_ctype_cmask_ctth(maia_h5):
    """Read cloudtype and flag info from filename
    """
    ctth = CtthObj()
    ctype = CtypeObj()
    cma = CmaObj()
    #ctth = CtthObj()
    cma.cma_ext
    maia_ct_bitflag = maia_h5['DATA']['CloudType'].value  # bit-4-8
    maia_ct = get_maia_ct_flag(maia_ct_bitflag)
    maia_ct[maia_ct == 7] = 6
    maia_ct[maia_ct == 9] = 7
    maia_ct[maia_ct == 11] = 8
    maia_ct[maia_ct == 13] = 9
    maia_ct[maia_ct == 15] = 11
    maia_ct[maia_ct == 16] = 12
    maia_ct[maia_ct == 17] = 13
    maia_ct[maia_ct == 18] = 14
    maia_ct[maia_ct == 19] = 10

    ctype.cloudtype = maia_ct
    cma.cma_ext = 0 * maia_ct
    cma.cma_ext[maia_ct == 0] = 255
    cma.cma_ext[maia_ct == 1] = 0
    cma.cma_ext[maia_ct == 2] = 0
    cma.cma_ext[maia_ct == 3] = 3
    cma.cma_ext[maia_ct == 4] = 3
    cma.cma_ext[np.logical_and(maia_ct > 4, maia_ct < 20)] = 1
    #0 -nodata
    # 1-2 clear
    #34 -snow/ice
    # 5-19 cloudy
    ctype.phaseflag = None
    ctype.ct_conditions = None
    ctth.height = 5000.0 * maia_ct
    ctth.temperature = ATRAIN_MATCH_NODATA + 0 * maia_ct
    ctth.pressure = ATRAIN_MATCH_NODATA + 0 * maia_ct
    # ctype.landseaflag =
    return ctype, cma, ctth


def read_maia_angobj(maia_h5):
    """Read angles info from filename
    """
    avhrrAngObj = imagerAngObj()
    avhrrAngObj.satz.data = 0.01 * maia_h5['DATA']['Sat_zenith'].value
    avhrrAngObj.sunz.data = ATRAIN_MATCH_NODATA + 0 * avhrrAngObj.satz.data
    avhrrAngObj.azidiff.data = ATRAIN_MATCH_NODATA + 0 * avhrrAngObj.satz.data
    maia_cma_bitflag = maia_h5['DATA']['CloudMask'].value
    daynight_flags = get_day_night_twilight_info_maia(maia_cma_bitflag)
    (no_qflag, night_flag, twilight_flag,
     day_flag, all_dnt_flag) = daynight_flags
    # Since we dont have sunz angle, set proxy for them using illumination flag
    # This to not include object in reshaped file just for maia!
    avhrrAngObj.sunz.data[night_flag == 1] = 120
    avhrrAngObj.sunz.data[twilight_flag == 1] = 85
    avhrrAngObj.sunz.data[day_flag == 1] = 00
    return avhrrAngObj


def read_maia_surftemp(maia_h5):
    """Read surftemp from file
    """
    surftemp = 273.15 + 0.01 * maia_h5['DATA']['Tsurf'].value
    surftemp[maia_h5['DATA']['Tsurf'].value == -9999] = ATRAIN_MATCH_NODATA
    return surftemp


def read_maia_geoobj(maia_h5, filename):
    """Read geolocation and time info from filename
    """
    GeoObj = imagerGeoObj()
    GeoObj.latitude = 0.0001 * maia_h5['DATA']['Latitude'].value
    GeoObj.longitude = 0.0001 * maia_h5['DATA']['Longitude'].value
    GeoObj.nodata = 999
    GeoObj.longitude[GeoObj.longitude == -9999] = GeoObj.nodata
    GeoObj.latitude[GeoObj.latitude == -9999] = GeoObj.nodata

    # viiCT_npp_GL_20150711_S211124_E211248_ASC_D_La-40_Lo-108_19188.h5
    # viiCT_npp_DB_20120817_S035411_E035535_DES_N_La052_Lo-027_00001.h5
    sl_ = os.path.basename(filename).split('_')
    date_time_start = datetime.strptime(sl_[3] + sl_[4], '%Y%m%dS%H%M%S')
    date_time_end = datetime.strptime(sl_[3] + sl_[5], '%Y%m%dE%H%M%S')

    GeoObj.sec1970_start = calendar.timegm(date_time_start.timetuple())
    GeoObj.sec1970_end = calendar.timegm(date_time_end.timetuple())
    GeoObj.num_of_lines = GeoObj.latitude.shape[0]

    GeoObj = createAvhrrTime(GeoObj, values={}, Trust_sec_1970=True)
    tim1 = time.strftime("%Y%m%d %H:%M",
                         time.gmtime(GeoObj.sec1970_start))
    tim2 = time.strftime("%Y%m%d %H:%M",
                         time.gmtime(GeoObj.sec1970_end))
    logger.info("Starttime: %s, end time: %s" % (tim1, tim2))
    logger.info("Min lon: %f, max lon: %d" % (
        np.min(np.where(
                np.equal(GeoObj.longitude, GeoObj.nodata),
                99999,
                GeoObj.longitude)),
        np.max(GeoObj.longitude)))
    logger.info("Min lat: %d, max lat: %d" % (
        np.min(np.where(
                np.equal(GeoObj.latitude, GeoObj.nodata),
                99999,
                GeoObj.latitude)),
        np.max(GeoObj.latitude)))

    return GeoObj

if __name__ == "__main__":
    pass
