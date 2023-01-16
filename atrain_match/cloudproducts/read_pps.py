#-*- coding: utf-8 -*-
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
"""Reader for PPS cloudproducts and intermediate files."""

from atrain_match.utils.runutils import do_some_geo_obj_logging
import atrain_match.config as config
import numpy as np
import os
import netCDF4
import h5py
import logging
import time
import calendar
from datetime import datetime
logger = logging.getLogger(__name__)
ATRAIN_MATCH_NODATA = config.NODATA
# logger.debug('Just so you know: this module has a logger...')


def get_satid_datetime_orbit_from_fname_pps(imager_filename, as_oldstyle=False):
    """Get sattelite start time, orbit etc from PPS filename."""
    # satname, _datetime, orbit = runutils.parse_scene(imager_filename)
    # returnd orbit as int, loosing leeding zeros, use %05d to get it right.
    # Get satellite name, time, and orbit number from imager_file
    if config.PPS_FORMAT_2012_OR_EARLIER or as_oldstyle:
        sl_ = os.path.basename(imager_filename).split('_')
        date_time = datetime.strptime(sl_[1] + sl_[2], '%Y%m%d%H%M')
        values = {"satellite": sl_[0],
                  "date_time": date_time,
                  "orbit": sl_[3].split('.')[0],
                  "date": sl_[1],
                  "year": date_time.year,
                  "month": "%02d" % (date_time.month),
                  # "lines_lines": sl_[5] + "_" + sl_[6],
                  "lines_lines": "*",
                  "time": sl_[2],
                  "imager_filename": imager_filename}
    else:  # PPS v2014-filenames
        sl_ = os.path.basename(imager_filename).split('_')
        date_time = datetime.strptime(sl_[5], '%Y%m%dT%H%M%S%fZ')
        date_time_end = datetime.strptime(sl_[6].split('Z')[0], '%Y%m%dT%H%M%S%f')

        values = {"satellite": sl_[3],
                  "date_time": date_time,
                  "date_time_end": date_time_end,
                  "orbit": sl_[4],
                  # "date":"%04d%02d%02d"%(date_time.year, dat_time.month, date_time.day),
                  "date": date_time.strftime("%Y%m%d"),
                  "year": date_time.year,
                  "month": "%02d" % (date_time.month),
                  "lines_lines": "*",
                  "time": date_time.strftime("%H%M%S"),
                  "ppsfilename": imager_filename}
    values['basename'] = (values["satellite"] +
                          "_" + values["date"] +
                          "_" + values["time"] +
                          "_" + values["orbit"])
    values["jday"] = date_time.timetuple().tm_yday

    return values


def create_imager_time(obt, values=None):
    """Make matrix with time for each pixel from start and end time."""
    if obt.sec1970_start < 0:
        obt.sec1970_start = calendar.timegm(values['date_time'])
    if obt.sec1970_end < obt.sec1970_start:
        obt.sec1970_end = calendar.timegm(values['date_time_end'])
    obt.time = np.linspace(obt.sec1970_start, obt.sec1970_end, obt.num_of_lines)
    return obt


class AllImagerData(object):
    """Class to hold all imager cloudproduct data."""

    def __init__(self, array_dict=None):
        """Init with data from array_dict, or empty."""
        self.latitude = None
        self.longitude = None
        self.nodata = None  # latlon nodata
        self.sec1970_end = None
        self.sec1970_start = None
        self.time = None
        self.num_of_lines = None
        self.instrument = "imager"
        # self.imager_geo = None
        self.imager_angles = None
        self.imager_channeldata = None
        self.ctype = None
        self.cma = None
        self.ctth = None
        self.unc = None
        self.aux = AuxiliaryObj({'surftemp': None})
        self.cpp = None
        self.nwp_segments = None
        if array_dict is not None:
            self.__dict__.update(array_dict)


class AuxiliaryObj(object):
    """Class to hold auxiliary cloudproduct data."""

    def __init__(self, array_dict):
        """Init with arrays from dictionary array_dict."""
        self.surftemp = None
        self.t500 = None
        self.t700 = None
        self.t850 = None
        self.t950 = None
        self.t1000 = None
        self.t900 = None
        self.t250 = None
        self.t800 = None
        self.psur = None
        self.ptro = None
        self.t2m = None
        self.ciwv = None
        self.text_r06 = None
        self.text_t11 = None
        self.text_t37t12 = None
        self.text_t11t12 = None
        self.text_t37 = None
        self.thr_t11ts_inv = None
        self.thr_t11t37_inv = None
        self.thr_t37t12_inv = None
        self.thr_t11t12_inv = None
        self.thr_t85t11_inv = None
        self.thr_t11ts = None
        self.thr_t11t37 = None
        self.thr_t37t12 = None
        self.thr_t11t12 = None
        self.thr_t85t11 = None
        self.thr_r06 = None
        self.thr_r09 = None
        self.emis1 = None
        self.emis6 = None
        self.emis8 = None
        self.emis9 = None
        self.snowa = None
        self.snowd = None
        self.seaice = None
        self.r37 = None
        self.landuse = None
        self.fractionofland = None
        self.elevation = None
        self.__dict__.update(array_dict)


class SmallDataObject(object):
    """Object to hold data, gain intercept etc.

    This might be unneeded, as we now scale data directly when
    reading.

    TODO:
        Remove

    """

    def __init__(self):
        self.data = None
        self.gain = 1.0
        self.intercept = 0.0
        self.no_data = 0.0
        self.missing_data = 0.0


class ImagerAngObj(object):
    """Class to hold imager angle data."""

    def __init__(self):
        self.satz = SmallDataObject()
        self.sunz = SmallDataObject()
        self.satazimuth = SmallDataObject()
        self.sunazimuth = SmallDataObject()
        self.azidiff = SmallDataObject()


class NewImagerData:
    """Class to hold imager channel data."""

    def __init__(self):
        self.channel = {}


class ImagerChannelData:
    """Class to hold imager channel data for one channel."""

    def __init__(self):
        """Init empty object."""
        self.channel = ""
        self.des = ""
        self.gain = 1.0
        self.intercept = 0.0
        self.SZA_corr_done = False
        self.data = None


class CmaObj:
    """Class to hold imager data for CMA and CMAPROB."""

    def __init__(self):
        """Init all datasets to None."""
        self.cma_ext = None
        self.cma_quality = None
        self.cma_bin = None
        self.cma_prob = None
        self.cma_aflag = None
        self.cma_testlist0 = None  # new in v2018
        self.cma_testlist1 = None  # new in v2018
        self.cma_testlist2 = None  # new in v2018
        self.cma_testlist3 = None  # new in v2018
        self.cma_testlist4 = None  # new in v2018
        self.cma_testlist5 = None  # new in v2018


class CtypeObj:
    """Class to hold imager data for cloud type."""

    def __init__(self):
        """Init all datasets to None."""
        self.cloudtype = None
        self.ct_statusflag = None
        self.ct_quality = None
        self.ct_conditions = None
        self.phaseflag = None  # v2012


class CtthObj:
    """Class to hold imager data for CTTH."""

    def __init__(self):
        """Init all datasets to None."""
        self.height = None
        self.height_corr = None
        self.temperature = None
        self.pressure = None
        self.ctth_statusflag = None


class CppObj:
    """Class to hold imager data for CPP."""

    def __init__(self):
        """Init all datasets to None."""
        self.cpp_cot = None
        self.cpp_cwp = None
        self.cpp_dcot = None
        self.cpp_dcwp = None
        self.cpp_lwp = None
        self.cpp_iwp = None
        self.cpp_quality = None
        self.cpp_status_flag = None
        self.cpp_conditions = None
        self.cpp_phase = None
        self.cpp_phase_extended = None
        self.cpp_reff = None


def read_ctth_h5(filename):
    h5file = h5py.File(filename, 'r')
    ctth = CtthObj()
    ctth.height = h5file['ctth_alti'].value.astype(np.float)
    ctth.temperature = h5file['ctth_tempe'].value.astype(np.float)
    ctth.pressure = h5file['ctth_pres'].value.astype(np.float)
    ctth.ctth_statusflag = h5file['ctth_status_flag'].value
    ctth.h_gain = h5file['ctth_alti'].attrs['scale_factor']
    ctth.h_intercept = h5file['ctth_alti'].attrs['add_offset']
    ctth.t_gain = h5file['ctth_tempe'].attrs['scale_factor']
    ctth.t_intercept = h5file['ctth_tempe'].attrs['add_offset']
    ctth.p_gain = h5file['ctth_pres'].attrs['scale_factor']
    ctth.p_intercept = h5file['ctth_pres'].attrs['add_offset']
    ctth.h_nodata = h5file['ctth_alti'].attrs['_FillValue']
    ctth.t_nodata = h5file['ctth_tempe'].attrs['_FillValue']
    ctth.p_nodata = h5file['ctth_pres'].attrs['_FillValue']
    hmask = ctth.height == ctth.h_nodata
    tmask = ctth.temperature == ctth.t_nodata
    pmask = ctth.pressure == ctth.p_nodata
    ctth.height[~hmask] = ctth.height[~hmask]*ctth.h_gain + ctth.h_intercept
    ctth.pressure[~pmask] = ctth.pressure[~pmask]*ctth.p_gain + ctth.p_intercept
    ctth.temperature[~tmask] = ctth.temperature[~tmask]*ctth.t_gain + ctth.t_intercept
    ctth.height[hmask] = ATRAIN_MATCH_NODATA
    ctth.pressure[pmask] = ATRAIN_MATCH_NODATA
    ctth.temperature[tmask] = ATRAIN_MATCH_NODATA

    logger.debug("min-h: %d, max-h: %d, h_nodata: %d",
                 np.min(ctth.height), np.max(ctth.height), ctth.h_nodata)
    h5file.close()
    return ctth


def read_ctth_nc(filename):
    """Read for PPS CTTH from netcdf file."""
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    ctth = CtthObj()
    ctth.height = pps_nc.variables['ctth_alti'][0, :, :].astype(np.float)
    ctth.temperature = pps_nc.variables['ctth_tempe'][0, :, :].astype(np.float)
    ctth.pressure = pps_nc.variables['ctth_pres'][0, :, :].astype(np.float)
    ctth.ctth_statusflag = pps_nc.variables['ctth_status_flag'][0, :, :]
    # Currently unpacked arrays later in calipso.py
    ctth.h_gain = 1.0
    ctth.h_intercept = 0.0
    ctth.p_gain = 1.0  #
    ctth.p_intercept = 0.0
    ctth.t_gain = 1  #
    ctth.t_intercept = 0.0
    ctth.h_nodata = pps_nc.variables['ctth_alti']._FillValue
    ctth.t_nodata = pps_nc.variables['ctth_tempe']._FillValue
    ctth.p_nodata = pps_nc.variables['ctth_pres']._FillValue
    logger.debug("min-h: %d, max-h: %d, h_nodata: %d",
                 np.min(ctth.height), np.max(ctth.height.data), ctth.h_nodata)
    # already scaled
    if np.ma.is_masked(ctth.height):
        ctth.height.data[ctth.height.mask] = ATRAIN_MATCH_NODATA
        ctth.height = ctth.height.data
    if np.ma.is_masked(ctth.pressure):
        ctth.pressure.data[ctth.pressure.mask] = ATRAIN_MATCH_NODATA
        ctth.pressure = ctth.pressure.data
    if np.ma.is_masked(ctth.temperature):
        ctth.temperature.data[ctth.temperature.mask] = ATRAIN_MATCH_NODATA
        ctth.temperature = ctth.temperature.data
    pps_nc.close()
    return ctth


def read_cloudtype_h5(filename):
    h5file = h5py.File(filename, 'r')
    ctype = CtypeObj()
    ctype.cloudtype = h5file['ct'].value
    ctype.ct_conditions = h5file['ct_conditions'].value
    ctype.ct_statusflag = h5file['ct_status_flag'].value
    ctype.ct_quality = h5file['ct_quality'].value
    h5file.close()
    return ctype


def read_cloudtype_nc(filename):
    """Read for PPS cloud type from netcdf file."""
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    ctype = CtypeObj()
    ctype.cloudtype = pps_nc.variables['ct'][0, :, :]
    ctype.ct_conditions = pps_nc.variables['ct_conditions'][0, :, :]
    ctype.ct_statusflag = pps_nc.variables['ct_status_flag'][0, :, :]
    ctype.ct_quality = pps_nc.variables['ct_quality'][0, :, :]
    pps_nc.close()
    return ctype


def read_cma_h5(filename):
    h5file = h5py.File(filename, 'r')
    cma = CmaObj()
    if 'cma_extended' not in h5file.keys():
        if 'cloud_probability' in h5file.keys():
            logger.error("This CMA-file seem lika a CMAPROB file!")
    cma.cma_ext = h5file['cma_extended'].value
    cma.cma_bin = np.int64(0*cma.cma_ext.copy())
    cma.cma_bin[cma.cma_ext == 1] = 1.0
    cma.cma_bin[cma.cma_ext == 2] = 1.0
    cma.cma_bin[cma.cma_ext < 0] = ATRAIN_MATCH_NODATA
    cma.cma_bin[cma.cma_ext > 10] = ATRAIN_MATCH_NODATA

    # try KeyError 'cma'
    h5file.close()
    return cma


def read_cmaprob_h5(filename, cma):
    h5file = h5py.File(filename, 'r')
    if cma is None:
        cma = CmaObj()
    if 'cma_extended' in h5file.keys():
        if 'cloud_probability' not in h5file.keys():
            logger.error("\n This file looks like a normal CMA-file (not CMAPROB) %s", filename)
    name = 'cloud_probability'
    if name not in h5file.keys():
        logger.info("This CMA-file is old.")
        name = "cmaprob"
    cma.cma_prob = h5file[name].value
    h5file.close()
    return cma


def read_cmaprob_nc(filename, cma):
    """Read for PPS cloud probability from netcdf file."""
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    if cma is None:
        cma = CmaObj()
    if 'cma_extended' in pps_nc.variables.keys():
        if 'cmaprob' not in pps_nc.variables.keys():
            logger.error("\n This file looks not like a CMAPROB-file like a normal CMA-file %s", filename)
    cma.cma_prob = pps_nc.variables['cmaprob'][0, :, :]
    if np.ma.is_masked(cma.cma_prob):
        mask = cma.cma_prob.mask
        cma.cma_prob =  cma.cma_prob.data
        cma.cma_prob[mask] =  ATRAIN_MATCH_NODATA
    pps_nc.close()
    return cma


def read_cma_nc(filename):
    """Read for PPS cloud mask from netcdf file."""
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    cma = CmaObj()
    if 'cma_extended' not in pps_nc.variables.keys():
        if 'cloud_probability' in pps_nc.variables.keys():
            logger.error("Probably you shold set CMAP_PROB_VALIDATION=True!")
    cma.cma_ext = pps_nc.variables['cma_extended'][0, :, :]
    cma.cma_bin = 0*np.int64(cma.cma_ext.copy())
    cma.cma_bin[cma.cma_ext == 1] = 1.0
    cma.cma_bin[cma.cma_ext == 2] = 1.0
    cma.cma_bin[cma.cma_ext < 0] = ATRAIN_MATCH_NODATA
    cma.cma_bin[cma.cma_ext > 10] = ATRAIN_MATCH_NODATA
    cma.cma_quality = pps_nc.variables['cma_quality'][0, :, :]

    # cma.cma_testlist0 = pps_nc.variables['cma_testlist0'][0, :, :]
    # cma.cma_testlist1 = pps_nc.variables['cma_testlist1'][0, :, :]
    # cma.cma_testlist2 = pps_nc.variables['cma_testlist2'][0, :, :]
    # cma.cma_testlist3 = pps_nc.variables['cma_testlist3'][0, :, :]
    # cma.cma_testlist4 = pps_nc.variables['cma_testlist4'][0, :, :]
    # cma.cma_testlist5 = pps_nc.variables['cma_testlist5'][0, :, :]
    for var_name in [
            'cma_aerosol',  # new updated name
            'cma_aerosolflag',
            'cma_dust',
            'cma_testlist0',
            'cma_testlist1',
            'cma_testlist2',
            'cma_testlist3',
            'cma_testlist4',
            'cma_testlist5']:
        if var_name in pps_nc.variables.keys():
            array = pps_nc.variables[var_name][0, :, :]
            atrain_name = var_name
            if var_name == 'cma_aerosol':
                atrain_name = 'cma_aerosolflag'
            if np.ma.is_masked(array):
                mask = array.mask
                data = array.data
                data[mask] = 0
                setattr(cma, atrain_name, data)
            else:
                setattr(cma, atrain_name, array)
        else:
            logger.info("No %s in cma file", var_name)
    pps_nc.close()
    return cma


def read_imager_data_nc(pps_nc):
    """Read for PPS level1c from netcdf file."""
    imager_data = NewImagerData()

    bad = None
    if 'qual_flags' in pps_nc.variables.keys():
        qual = pps_nc.variables['qual_flags'][0, :, :]
        bad = np.sum(qual[:, 1:], axis=1) > 0
        
    for var in pps_nc.variables.keys():
        if 'image' in var or hasattr(pps_nc.variables[var], "sensor"):
            image = pps_nc.variables[var]
            if not hasattr(image, 'id_tag'):
                continue
            id_tag = image.id_tag
            if image.id_tag in ['satzenith', 'sunzenith',
                                'azimuthdiff',
                                'sunazimuth', 'satazimuth']:
                continue
            logger.debug("reading channel %s", id_tag)
            one_channel = ImagerChannelData()
            # channel = image.channel
            data_temporary = image[0, :, :]
            if np.ma.is_masked(one_channel.data):
                one_channel.data = data_temporary.data
                one_channel.data[data_temporary.mask] = ATRAIN_MATCH_NODATA
            else:
                one_channel.data = data_temporary
            if bad is not None and id_tag != 'qual_flags':
                one_channel.data[bad, :] = ATRAIN_MATCH_NODATA
            if hasattr(image, 'description'):
                logger.debug("reading channel %s", image.description)
                one_channel.des = image.description
            one_channel.id_tag = image.id_tag
            if hasattr(pps_nc.variables[var], "sun_zenith_angle_correction_applied"):
                corr_done_attr = pps_nc.variables[var].sun_zenith_angle_correction_applied
                if corr_done_attr.upper() in ["TRUE"]:
                    one_channel.SZA_corr_done = True
            one_channel.gain = 1.0  # data already scaled
            one_channel.intercept = 0.0  # data already scaled
            imager_data.channel[id_tag] = one_channel
            fill_value = pps_nc.variables[var]._FillValue
            imager_data.nodata = fill_value
            imager_data.missingdata = fill_value
            imager_data.no_data = fill_value
            imager_data.missing_data = fill_value
    return imager_data


def read_pps_angobj_nc(pps_nc):
    """Read for PPS level1c from netcdf file."""
    angle_obj = ImagerAngObj()
    for varname in pps_nc.variables.keys():
        this_is = "non_angle_variable"
        if varname in ['satzenith',
                       'sunzenith',
                       'azimuthdiff',
                       'sunazimuth',
                       'satazimuth']:
            this_is = varname
        else:
            # some times we got angels in imager file
            # they are then named imageX as varname
            if hasattr(pps_nc.variables[varname], 'id_tag'):
                this_is = pps_nc.variables[varname].id_tag

        if this_is in['satzenith']:
            angle_obj.satz.data = pps_nc.variables[varname][0, :, :].astype(np.float)
            angle_obj.satz.no_data = pps_nc.variables[varname]._FillValue
            angle_obj.satz.intercept = pps_nc.variables[varname].add_offset
            angle_obj.satz.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['sunzenith']:
            angle_obj.sunz.data = pps_nc.variables[varname][0, :, :].astype(np.float)
            angle_obj.sunz.no_data = pps_nc.variables[varname]._FillValue
            angle_obj.sunz.intercept = pps_nc.variables[varname].add_offset
            angle_obj.sunz.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['azimuthdiff']:
            angle_obj.azidiff.data = pps_nc.variables[varname][0, :, :].astype(np.float)
            angle_obj.azidiff.no_data = pps_nc.variables[varname]._FillValue
            angle_obj.azidiff.intercept = pps_nc.variables[varname].add_offset
            angle_obj.azidiff.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['sunazimuth']:
            angle_obj.sunazimuth.data = pps_nc.variables[varname][0, :, :].astype(np.float)
            angle_obj.sunazimuth.no_data = pps_nc.variables[varname]._FillValue
            angle_obj.sunazimuth.intercept = pps_nc.variables[varname].add_offset
            angle_obj.sunazimuth.gain = pps_nc.variables[varname].scale_factor
        elif this_is in['satazimuth']:
            angle_obj.satazimuth.data = pps_nc.variables[varname][0, :, :].astype(np.float)
            angle_obj.satazimuth.no_data = pps_nc.variables[varname]._FillValue
            angle_obj.satazimuth.intercept = pps_nc.variables[varname].add_offset
            angle_obj.satazimuth.gain = pps_nc.variables[varname].scale_factor

    # already scaled
    if np.ma.is_masked(angle_obj.sunz.data):
        angle_obj.sunz.data[angle_obj.sunz.data.mask] = ATRAIN_MATCH_NODATA
        angle_obj.sunz.data = angle_obj.sunz.data.data
    if np.ma.is_masked(angle_obj.satz.data):
        angle_obj.satz.data[angle_obj.satz.data.mask] = ATRAIN_MATCH_NODATA
        angle_obj.satz.data = angle_obj.satz.data.data
    if np.ma.is_masked(angle_obj.azidiff.data):
        angle_obj.azidiff.data[angle_obj.azidiff.data.mask] = ATRAIN_MATCH_NODATA
        angle_obj.azidiff.data = angle_obj.azidiff.data.data
    if np.ma.is_masked(angle_obj.sunazimuth.data):
        angle_obj.sunazimuth.data[angle_obj.sunazimuth.data.mask] = ATRAIN_MATCH_NODATA
        angle_obj.sunazimuth.data = angle_obj.sunazimuth.data.data
    if np.ma.is_masked(angle_obj.satazimuth.data):
        angle_obj.satazimuth.data[angle_obj.satazimuth.data.mask] = ATRAIN_MATCH_NODATA
        angle_obj.satazimuth.data = angle_obj.satazimuth.data.data
    return angle_obj


def read_pps_angobj_h5(filename):
    """Read angles info from file filename."""
    h5file = h5py.File(filename, 'r')
    angle_obj = ImagerAngObj()

    for var in h5file.keys():
        if 'image' in var:
            image = h5file[var]
            if (image.attrs['description'] == "sun zenith angle" or
                    image.attrs['description'] == "Solar zenith angle"):
                # print "reading sunz"
                angle_obj.sunz.data = image['data'].value.astype(np.float)
                angle_obj.sunz.gain = image['what'].attrs['gain']
                angle_obj.sunz.intercept = image['what'].attrs['offset']
                angle_obj.sunz.no_data = image['what'].attrs['nodata']
                angle_obj.sunz.missing_data = image['what'].attrs['missingdata']
            elif (image.attrs['description'] == "satellite zenith angle" or
                  image.attrs['description'] == "Satellite zenith angle"):
                angle_obj.satz.data = image['data'].value.astype(np.float)
                angle_obj.satz.gain = image['what'].attrs['gain']
                angle_obj.satz.intercept = image['what'].attrs['offset']
                angle_obj.satz.no_data = image['what'].attrs['nodata']
                angle_obj.satz.missing_data = image['what'].attrs['missingdata']
            elif (image.attrs['description'] ==
                  "relative sun-satellite azimuth difference angle" or
                  image.attrs['description'] ==
                  "Relative satellite-sun azimuth angle"):
                angle_obj.azidiff.data = image['data'].value.astype(np.float)
                angle_obj.azidiff.gain = image['what'].attrs['gain']
                angle_obj.azidiff.intercept = image['what'].attrs['offset']
                angle_obj.azidiff.no_data = image['what'].attrs['nodata']
                angle_obj.azidiff.missing_data = image['what'].attrs['missingdata']
    sunzmask = np.logical_or(angle_obj.sunz.data == angle_obj.sunz.no_data,
                             angle_obj.sunz.data == angle_obj.sunz.missing_data)
    satzmask = np.logical_or(angle_obj.satz.data == angle_obj.satz.no_data,
                             angle_obj.satz.data == angle_obj.satz.missing_data)
    diffmask = np.logical_or(angle_obj.azidiff.data == angle_obj.azidiff.no_data,
                             angle_obj.azidiff.data == angle_obj.azidiff.missing_data)
    angle_obj.sunz.data[~sunzmask] = angle_obj.sunz.data[~sunzmask]*angle_obj.sunz.gain + angle_obj.sunz.intercept
    angle_obj.satz.data[~satzmask] = angle_obj.satz.data[~satzmask]*angle_obj.satz.gain + angle_obj.satz.intercept
    angle_obj.azidiff.data[~diffmask] = angle_obj.azidiff.data[~diffmask] * \
        angle_obj.azidiff.gain + angle_obj.azidiff.intercept
    angle_obj.sunz.data[sunzmask] = ATRAIN_MATCH_NODATA
    angle_obj.satz.data[satzmask] = ATRAIN_MATCH_NODATA
    angle_obj.azidiff.data[diffmask] = ATRAIN_MATCH_NODATA
    h5file.close()
    return angle_obj


def read_pps_geoobj_nc(pps_nc):
    """Read geolocation and time info from netcdf file."""
    all_imager_obj = AllImagerData()

    longitude = pps_nc.variables['lon'][::]
    if np.ma.is_masked(longitude):
        all_imager_obj.longitude = pps_nc.variables['lon'][::].data
        all_imager_obj.longitude[longitude.mask] = -999.0
    else:
        all_imager_obj.longitude = longitude
    latitude = pps_nc.variables['lat'][::]
    if np.ma.is_masked(latitude):
        all_imager_obj.latitude = pps_nc.variables['lat'][::].data
        all_imager_obj.latitude[latitude.mask] = -999.0
    else:
        all_imager_obj.latitude = latitude
    all_imager_obj.nodata = pps_nc.variables['lon']._FillValue
    all_imager_obj.num_of_lines = all_imager_obj.latitude.shape[0]
    # import pdb
    # pdb.set_trace()
    time_temp = pps_nc.variables['time'].units  # to 1970 s
    seconds = np.float64(pps_nc.variables['time'][::])  # from PPS often 0
    if 'seconds' in time_temp:
        if 'T' in time_temp:
            time_obj = time.strptime(time_temp, 'seconds since %Y-%m-%dT%H:%M:%S+00:00')
        elif time_temp == u'seconds since 1970-01-01':
            time_obj = time.strptime(time_temp, 'seconds since %Y-%m-%d')
        else:
            time_obj = time.strptime(time_temp, 'seconds since %Y-%m-%d %H:%M:%S.%f +00:00')

        sec_since_1970 = calendar.timegm(time_obj)

        all_imager_obj.sec1970_start = (sec_since_1970 +
                                        np.float64(np.min(pps_nc.variables['time_bnds'][::])) +
                                        seconds)
        all_imager_obj.sec1970_end = (sec_since_1970 +
                                      np.float64(np.max(pps_nc.variables['time_bnds'][::])) +
                                      seconds)
    else:
        try:
            all_imager_obj.sec1970_start = calendar.timegm(time.strptime(pps_nc.variables['lon'].start_time,
                                                                         '%Y-%m-%d %H:%M:%S.%f'))
            all_imager_obj.sec1970_end = calendar.timegm(time.strptime(pps_nc.variables['lon'].end_time,
                                                                       '%Y-%m-%d %H:%M:%S.%f'))

        except:    
            all_imager_obj.sec1970_start = calendar.timegm(time.strptime(pps_nc.start_time,  '%Y-%m-%d %H:%M:%S'))
            all_imager_obj.sec1970_end = calendar.timegm(time.strptime(pps_nc.end_time,  '%Y-%m-%d %H:%M:%S'))
    # print type(all_imager_obj.sec1970_start)
    all_imager_obj.sec1970_start = np.float64(all_imager_obj.sec1970_start)
    all_imager_obj.sec1970_end = np.float64(all_imager_obj.sec1970_end)
    do_some_geo_obj_logging(all_imager_obj)
    return all_imager_obj


def read_pps_geoobj_h5(filename):
    """Read geolocation and time info from filename hdf5."""
    h5file = h5py.File(filename, 'r')
    all_imager_obj = AllImagerData()
    in_fillvalue1 = h5file['where/lon/what'].attrs['nodata']
    in_fillvalue2 = h5file['where/lon/what'].attrs['missingdata']
    all_imager_obj.nodata = -999.0
    gain = h5file['where/lon/what'].attrs['gain']
    intercept = h5file['where/lon/what'].attrs['offset']
    all_imager_obj.longitude = h5file['where/lon']['data'].value*gain + intercept
    all_imager_obj.latitude = h5file['where/lat']['data'].value*gain + intercept

    all_imager_obj.longitude[h5file['where/lon']['data'].value == in_fillvalue1] = all_imager_obj.nodata
    all_imager_obj.latitude[h5file['where/lat']['data'].value == in_fillvalue1] = all_imager_obj.nodata
    all_imager_obj.longitude[h5file['where/lon']['data'].value == in_fillvalue2] = all_imager_obj.nodata
    all_imager_obj.latitude[h5file['where/lat']['data'].value == in_fillvalue2] = all_imager_obj.nodata

    all_imager_obj.num_of_lines = all_imager_obj.latitude.shape[0]
    all_imager_obj.sec1970_start = h5file['how'].attrs['startepochs']
    all_imager_obj.sec1970_end = h5file['how'].attrs['endepochs']
    do_some_geo_obj_logging(all_imager_obj)

    return all_imager_obj


def read_cpp_h5(filename):
    density = 1e3
    h5file = h5py.File(filename, 'r')
    cpp_obj = CppObj()
    for cpp_key in cpp_obj.__dict__.keys():
        data = read_cpp_h5_one_var(h5file, cpp_key)
        if cpp_key in ["cpp_lwp"]:
            logger.debug("Convert from CPP-lwp from kg/m-2 to g/m-2")
            data[data > 0] = density * data[data > 0]
        setattr(cpp_obj, cpp_key, data)
    h5file.close()
    return cpp_obj


def read_cpp_h5_one_var(h5file, cpp_key):
    if cpp_key in h5file.keys():
        logger.debug("Read %s", cpp_key)
        cpp_var_value = h5file[cpp_key].value
        nodata = h5file[cpp_key].attrs['_FillValue']
        if cpp_key in ["cpp_phase", "cpp_phase_extended"]:
            gain = 1.0
            intercept = 0.0
        else:
            gain = h5file[cpp_key].attrs['scale_factor']
            intercept = h5file[cpp_key].attrs['add_offset']

        cpp_data = np.where(cpp_var_value != nodata,
                            cpp_var_value * gain + intercept,
                            ATRAIN_MATCH_NODATA)
        return cpp_data
    else:
        logger.debug("NO %s field, Continue ", cpp_key)
        return None


def read_cpp_nc_one_var(ncFile, cpp_key):
    """Read one CPP dataset from netcdf file."""
    density = 1e3
    if cpp_key in ncFile.variables.keys():
        logger.debug("Read %s", cpp_key)
        cpp_var = ncFile.variables[cpp_key][0, :, :]
        if np.ma.is_masked(cpp_var):
            cpp_data = cpp_var.data.astype(np.float)
            cpp_data[cpp_var.mask] = ATRAIN_MATCH_NODATA
        else:
            cpp_data = cpp_var
        if cpp_key in ["cpp_lwp", "cmic_lwp"]:
            logger.debug("Convert from CPP-lwp from to g/m-2")
            cpp_data[cpp_data > 0] = density * cpp_data[cpp_data > 0]
        return cpp_data
    else:
        logger.debug("NO %s field, Continue ", cpp_key)
        return None


def read_cpp_nc(filename):
    """Read all CPP datasets from netcdf file."""
    pps_nc_cpp = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    cpp_obj = CppObj()
    for cpp_key in cpp_obj.__dict__.keys():
        cpp_key_new = cpp_key
        data = read_cpp_nc_one_var(pps_nc_cpp, cpp_key)
        if data is None:    
            cpp_key_new = cpp_key.replace("cpp", "cmic")
            data = read_cpp_nc_one_var(pps_nc_cpp, cpp_key_new)    
        setattr(cpp_obj, cpp_key, data)
        if data is not None:
            logger.debug("Read cpp_keys: %s, %s", cpp_key, cpp_key_new)
    pps_nc_cpp.close()
    return cpp_obj


def read_nwp_h5(filename, nwp_key):

    import h5py
    h5file = h5py.File(filename, 'r')
    if nwp_key in h5file.keys():
        logger.debug("Read NWP %s", nwp_key)
        value = h5file[nwp_key].value
        gain = h5file[nwp_key].attrs['gain']
        intercept = h5file[nwp_key].attrs['intercept']
        nodat = h5file[nwp_key].attrs['nodata']
        data = np.where(value != nodat, value * gain + intercept, value)
        h5file.close()
        return data
    else:
        logger.debug("NO NWP %s File, Continue", nwp_key)
        return None


def read_etc_nc(ncFile, etc_key):
    """Read a datasets from netcdf file."""
    if etc_key in ncFile.variables.keys():
        logger.debug("Read %s", etc_key)
        nwp_var = ncFile.variables[etc_key][0, :, :]
        if np.ma.is_masked(nwp_var):
            if 'emis' in etc_key:
                # set emissivity 1.0 where we miss data
                nwp_var.data[nwp_var.mask] = 1.0
            nwp_var = nwp_var.data
        return nwp_var
    else:
        logger.debug("NO %s field, Continue ", etc_key)
        return None


def read_segment_data(filename):
    import h5py
    product = {}
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        for attribute in list(h5file.attrs):
            product[attribute] = h5file.attrs[attribute]
        for attribute in list(h5file['satdef'].attrs):
            product[attribute] = h5file['satdef'].attrs[attribute]
        product['colidx'] = h5file['satdef']['colidx'].value
        product['rowidx'] = h5file['satdef']['rowidx'].value
        logger.debug("Read segment info moisture")
        for moist_data in ['moist', 'surfaceMoist']:
            data = h5file['segment_area'][moist_data]
            gain = h5file.attrs['m_gain']
            intercept = h5file.attrs['m_intercept']
            nodata = h5file.attrs['m_nodata']
            # data[data!=nodata] = data[data!=nodata] * (gain) + intercept
            product[moist_data] = data
        logger.debug("Read segment info pressure")
        for pressure_data in ['pressure', 'surfacePressure', 'ptro']:
            # pressure is in Pa in segments file
            data = h5file['segment_area'][pressure_data]
            gain = h5file.attrs['p_gain']
            intercept = h5file.attrs['p_intercept']
            nodata = h5file.attrs['p_nodata']
            # data[data!=nodata] = data[data!=nodata] * (gain/100) + intercept/100 # Pa => hPa
            data_float = np.array(data, dtype=np.float)
            data_float[data != nodata] = data_float[data != nodata] * (gain/100) + intercept/100  # Pa => hPa
            product[pressure_data] = data_float
        logger.debug("Read segment info height")
        for geoheight_data in ['geoheight', 'surfaceGeoHeight']:
            # geo height is in meters in segment file
            data = h5file['segment_area'][geoheight_data]
            gain = h5file.attrs['h_gain']
            intercept = h5file.attrs['h_intercept']
            nodata = h5file.attrs['h_nodata']
            data[data != nodata] = data[data != nodata] * gain + intercept
            product[geoheight_data] = data
        logger.debug("Read segment info temperature")
        for temperature_data in ['temp']:
            # temperature are measured in Kelvin in segmentfile
            data = h5file['segment_area'][temperature_data]
            gain = h5file.attrs['t_gain']
            intercept = h5file.attrs['t_intercept']
            nodata = h5file.attrs['t_nodata']
            data_float = np.array(data, dtype=np.float)
            data_float[data != nodata] = data_float[data != nodata] * gain + intercept
            product[temperature_data] = data_float
        for temperature_data in ['t850', 'ttro', 'surfaceLandTemp', 'surfaceSeaTemp']:
            data = h5file['segment_area'][temperature_data]
            gain = h5file.attrs['t_gain']
            intercept = h5file.attrs['t_intercept']
            nodata = h5file.attrs['t_nodata']
            data_float = np.array(data, dtype=np.float)
            data_float[data != nodata] = data_float[data != nodata] * gain + intercept
            product[temperature_data] = data_float
        for misc_data in ['meanElevation', 'fractionOfLand']:
            product[misc_data] = h5file['segment_area'][misc_data]
        logger.debug("Read segment info brightness temperature")

        for tb_data in ['tb11clfree_sea',
                        'tb12clfree_sea',
                        'tb11clfree_land',
                        'tb12clfree_land',
                        'tb4clfree_sea',
                        'tb5clfree_sea',
                        'tb4clfree_land',
                        'tb5clfree_land'
                        'tb11lowcloud_sea',
                        'tb12lowcloud_sea',
                        'tb11lowcloud_land',
                        'tb12lowcloud_land']:
            try:

                data = h5file['segment_area'][tb_data]
                gain = h5file.attrs['t_gain']
                intercept = h5file.attrs['tb_intercept']
                nodata = h5file.attrs['t_nodata']
                data_float = np.array(data, dtype=np.float)
                data_float[data != nodata] = data_float[data != nodata] * gain + intercept
                name = tb_data
                name = name.replace('4', '11')
                name = name.replace('5', '12')
                product[name] = data_float
            except ValueError:
                pass
        h5file.close()
        return product
    else:
        logger.info("NO segment %s File, Continue", filename)
        return None


def read_thr_h5(filename, h5_obj_type, thr_type):
    import h5py
    product = None
    if thr_type in ["emis1", "emis6", "emis8", "emis9"]:
        if filename is not None:
            h5file = h5py.File(filename, 'r')
            if 1 == 1:  # h5_obj_type in h5file.keys():
                value = h5file[h5_obj_type].value
                gain = h5file.attrs['gain']
                intercept = h5file.attrs['intercept']
                product = value * gain + intercept
                product[product < 0] = 1.0
                product[product > 1.0] = 1.0
                logger.debug("Read EMIS: %s", thr_type)
            else:
                logger.info("ERROR", "Couldn't read %s file, Continue", thr_type)
            h5file.close()
        else:
            logger.debug("NO EMIS %s File, Continue", thr_type)
        return product
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        if h5_obj_type in h5file.keys():
            value = h5file[h5_obj_type].value
            gain = h5file[h5_obj_type].attrs['gain']
            intercept = h5file[h5_obj_type].attrs['intercept']
            product = value * gain + intercept
            logger.debug("Read THR: %s" % (thr_type))
        else:
            logger.error("Could not read %s File, Continue", thr_type)
        h5file.close()
    else:
        logger.debug("NO THR %s File, Continue", thr_type)
    return product


def read_imager_data_h5(filename):
    h5file = h5py.File(filename, 'r')
    imager_data = NewImagerData()
    for var in h5file.keys():
        if 'image' in var:
            image = h5file[var]
            if 'description' not in image.attrs.keys() and "modis" in filename:
                my_description = "MODIS %s" % (image.attrs['channel'])
            else:
                my_description = image.attrs['description']
            logger.debug("reading channel %s", my_description)
            one_channel = ImagerChannelData()
            one_channel.data = image['data'].value
            one_channel.des = my_description
            one_channel.gain = 1.0
            one_channel.intercept = 0.0
            gain = image['what'].attrs['gain']
            intercept = image['what'].attrs['offset']
            imager_data.channel.append(one_channel)
            imager_data.nodata = image['what'].attrs['nodata']
            imager_data.missingdata = image['what'].attrs['missingdata']
            imager_data.no_data = imager_data.nodata
            imager_data.missing_data = imager_data.missingdata
            mask = np.logical_or(one_channel.data == imager_data.no_data,
                                 one_channel.data == imager_data.missing_data)
            one_channel.data[~mask] = one_channel.data[~mask]*gain + intercept
            one_channel.data[mask] = ATRAIN_MATCH_NODATA
    h5file.close()
    return imager_data


def read_all_intermediate_files(pps_files, SETTINGS):
    """Read data sets from pps intermediate files."""

    CTTH_TYPES = SETTINGS['CTTH_TYPES']
    aux_dict = {}
    if pps_files.nnextra is None:
        pass
    else:
        pps_nc_nnextra = netCDF4.Dataset(pps_files.nnextra, 'r', format='NETCDF4')
        for item in pps_nc_nnextra.variables.keys():
            if item[0:2] =='nn':
                aux_dict[item] = read_etc_nc(pps_nc_nnextra, item)
        pps_nc_nnextra.close()
        
    if pps_files.seaice is None:
        pass
    elif '.nc' in pps_files.seaice:
        pps_nc_seaice = netCDF4.Dataset(pps_files.seaice, 'r', format='NETCDF4')
        aux_dict["seaice"] = read_etc_nc(pps_nc_seaice, "seaice")
        pps_nc_seaice.close()
    else:
        logger.info("Not reading PPS seaice data")
    if pps_files.physiography is None:
        logger.info("Not reading PPS physiography data")
    elif '.nc' in pps_files.physiography:
        pps_nc_physiography = netCDF4.Dataset(pps_files.physiography, 'r', format='NETCDF4')
        aux_dict["landuse"] = read_etc_nc(pps_nc_physiography, "landuse")
        aux_dict["fractionofland"] = read_etc_nc(pps_nc_physiography, "fractionofland")
        aux_dict["elevation"] = read_etc_nc(pps_nc_physiography, "elevation")
        pps_nc_physiography.close()
    else:
        logger.info("Not reading PPS physiography data")
    if pps_files.r37 is None:
        pass
    else:
        pps_nc_r37 = netCDF4.Dataset(pps_files.r37, 'r', format='NETCDF4')
        aux_dict["r37_sza_correction_done"] = read_etc_nc(pps_nc_r37, "r37")
        pps_nc_r37.close()
    if pps_files.nwp_tsur is None:
        pass
    elif '.nc' in pps_files.nwp_tsur:
        pps_nc_nwp = netCDF4.Dataset(pps_files.nwp_tsur, 'r', format='NETCDF4')
        aux_dict['surftemp'] = read_etc_nc(pps_nc_nwp, "tsur")
        aux_dict['t500'] = read_etc_nc(pps_nc_nwp, "t500")
        aux_dict['t700'] = read_etc_nc(pps_nc_nwp, "t700")
        aux_dict['t850'] = read_etc_nc(pps_nc_nwp, "t850")
        aux_dict['t950'] = read_etc_nc(pps_nc_nwp, "t950")
        aux_dict['ttro'] = read_etc_nc(pps_nc_nwp, "ttro")
        aux_dict['ciwv'] = read_etc_nc(pps_nc_nwp, "ciwv")
        aux_dict['t1000'] = read_etc_nc(pps_nc_nwp, "t1000")
        aux_dict['t900'] = read_etc_nc(pps_nc_nwp, "t900")
        aux_dict['t800'] = read_etc_nc(pps_nc_nwp, "t800")
        aux_dict['t250'] = read_etc_nc(pps_nc_nwp, "t250")
        aux_dict['ptro'] = read_etc_nc(pps_nc_nwp, "ptro")
        # psur in in Pa, as ctth_pressure is in Pa
        aux_dict['psur'] = read_etc_nc(pps_nc_nwp, "psur")
        aux_dict['t2m'] = read_etc_nc(pps_nc_nwp, "t2m")
        aux_dict['h2m'] = read_etc_nc(pps_nc_nwp, "h2m")
        aux_dict['u10m'] = read_etc_nc(pps_nc_nwp, "u10m")
        aux_dict['v10m'] = read_etc_nc(pps_nc_nwp, "v10m")
        aux_dict['snowa'] = read_etc_nc(pps_nc_nwp, "snowa")
        aux_dict['snowd'] = read_etc_nc(pps_nc_nwp, "snowd")
    else:
        aux_dict['surftemp'] = read_nwp_h5(pps_files.nwp_tsur, "tsur")
        aux_dict['t500'] = read_nwp_h5(pps_files.nwp_t500, "t500")
        aux_dict['t700'] = read_nwp_h5(pps_files.nwp_t700, "t700")
        aux_dict['t850'] = read_nwp_h5(pps_files.nwp_t850, "t850")
        aux_dict['t950'] = read_nwp_h5(pps_files.nwp_t950, "t950")
        aux_dict['ttro'] = read_nwp_h5(pps_files.nwp_ttro, "ttro")
        aux_dict['ciwv'] = read_nwp_h5(pps_files.nwp_ciwv, "ciwv")
    if pps_files.text_t11 is None:
        pass
        logger.info("Not reading PPS texture data")
    elif '.nc' in pps_files.text_t11:
        pps_nc_txt = netCDF4.Dataset(pps_files.text_t11, 'r', format='NETCDF4')
        for ttype in ['r06', 't11', 't37t12', 't37', 't11t12']:
            text_type = 'text_' + ttype
            aux_dict[text_type] = read_etc_nc(pps_nc_txt, ttype)
        pps_nc_txt.close()
    else:
        for ttype in ['r06', 't11', 't37t12', 't37']:
            h5_obj_type = ttype + '_text'
            text_type = 'text_' + ttype
            aux_dict[text_type] = read_thr_h5(getattr(pps_files, text_type),
                                              h5_obj_type, text_type)
    if pps_files.thr_t11ts is None:
        pass
        logger.info("Not reading PPS threshold data")
    elif '.nc' in pps_files.thr_t11ts:
        pps_nc_thr = netCDF4.Dataset(pps_files.thr_t11ts, 'r', format='NETCDF4')
        for nc_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv',
                            't11ts', 't11t37', 't37t12', 't11t12',
                            'r09', 'r06', 't85t11_inv', 't85t11']:
            thr_type = 'thr_' + nc_obj_type
            aux_dict[thr_type] = read_etc_nc(pps_nc_thr, nc_obj_type)
        pps_nc_thr.close()
    else:
        for h5_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv',
                            't11ts', 't11t37', 't37t12', 't11t12',
                            'r09', 'r06', 't85t11_inv', 't85t11']:
            thr_type = 'thr_' + h5_obj_type
            aux_dict[thr_type] = read_thr_h5(getattr(pps_files, thr_type),
                                             h5_obj_type, thr_type)
    if pps_files.emis is None:
        pass
        logger.info("Not reading PPS Emissivity data")
    elif '.nc' in pps_files.emis:
        pps_nc_thr = netCDF4.Dataset(pps_files.emis, 'r', format='NETCDF4')
        for emis_type in ['emis1', "emis6", 'emis8', 'emis9']:
            aux_dict[emis_type] = read_etc_nc(pps_nc_thr, emis_type)
        pps_nc_thr.close()
    else:
        for h5_obj_type in ['emis1', "emis6", 'emis8', 'emis9']:
            emis_type = h5_obj_type
            aux_dict[emis_type] = read_thr_h5(getattr(pps_files, "emis"),
                                              h5_obj_type, emis_type)
    if len(CTTH_TYPES) > 1:
        for ctth_type in CTTH_TYPES[1:]:  # already read first
            aux_dict[ctth_type] = read_ctth_nc(pps_files.ctth[ctth_type])
    aux_obj = AuxiliaryObj(aux_dict)
    return aux_obj


def pps_read_all(pps_files, imager_file, SETTINGS):
    """Read all PPS data and return cloudproducts object."""
    logger.info("Read Imager geolocation data")
    if '.nc' in imager_file:
        pps_nc = netCDF4.Dataset(imager_file, 'r', format='NETCDF4')
        cloudproducts = read_pps_geoobj_nc(pps_nc)
    else:
        # use mpop?
        cloudproducts = read_pps_geoobj_h5(imager_file)
    # create time info for each pixel
    values = get_satid_datetime_orbit_from_fname_pps(imager_file)
    cloudproducts = create_imager_time(cloudproducts, values)
    logger.info("Read sun and satellites angles data")
    if '.nc' in pps_files.sunsatangles:
        pps_nc_ang = netCDF4.Dataset(pps_files.sunsatangles, 'r', format='NETCDF4')
        cloudproducts.imager_angles = read_pps_angobj_nc(pps_nc_ang)
        pps_nc_ang.close()
    else:
        # use mpop?
        cloudproducts.imager_angles = read_pps_angobj_h5(pps_files.sunsatangles)
    logger.info("Read Imager data")
    if '.nc' in imager_file:
        cloudproducts.imager_channeldata = read_imager_data_nc(pps_nc)
        pps_nc.close()
    else:
        cloudproducts.imager_channeldata = read_imager_data_h5(imager_file)
    for imager in ["avhrr", "viirs", "modis", "seviri"]:
        if imager in os.path.basename(imager_file):
            cloudproducts.instrument = imager
    logger.debug("%s, %s, %s", pps_files.cloudtype, pps_files.ctth, pps_files.cma)

    # CPP
    if pps_files.cpp is not None:
        logger.info("Read CPP data")
        if '.nc' in pps_files.cpp:
            cloudproducts.cpp = read_cpp_nc(pps_files.cpp)
        else:
            cloudproducts.cpp = read_cpp_h5(pps_files.cpp)
    # CMA
    if pps_files.cma is not None:
        logger.info("Read PPS Cloud mask")
        logger.debug(pps_files.cma)
        if '.nc' in pps_files.cma:
            cloudproducts.cma = read_cma_nc(pps_files.cma)
        else:
            cloudproducts.cma = read_cma_h5(pps_files.cma)
    # CMAPROB
    if pps_files.cmaprob is not None:
        logger.info("Read PPS Cloud mask prob")
        if '.nc' in pps_files.cmaprob:
            cloudproducts.cma = read_cmaprob_nc(pps_files.cmaprob, cloudproducts.cma)
        else:
            cloudproducts.cma = read_cmaprob_h5(pps_files.cmaprob, cloudproducts.cma)
    # CTYPE
    if pps_files.cloudtype is not None:
        logger.info("Read PPS Cloud type")
        if '.nc' in pps_files.cloudtype:
            cloudproducts.ctype = read_cloudtype_nc(pps_files.cloudtype)
        else:
            cloudproducts.ctype = read_cloudtype_h5(pps_files.cloudtype)
    # CTTH
    CTTH_TYPES = SETTINGS["CTTH_TYPES"]
    if len(pps_files.ctth.keys()) >= 1:
        logger.info("Read PPS CTTH")
        if '.nc' in pps_files.ctth[CTTH_TYPES[0]]:
            # read first ctth as primary one
            cloudproducts.ctth = read_ctth_nc(pps_files.ctth[CTTH_TYPES[0]])
        else:
            cloudproducts.ctth = read_ctth_h5(pps_files.ctth[CTTH_TYPES[0]])

    logger.info("Read PPS full resolution intermediate files")
    cloudproducts.aux = read_all_intermediate_files(pps_files, SETTINGS)

    logger.info("Read PPS NWP segment resolution data")
    cloudproducts.nwp_segments = read_segment_data(getattr(pps_files, 'nwp_segments'))

    return cloudproducts


if __name__ == "__main__":
    pass
