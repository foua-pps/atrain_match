"""
Script to collocate SEVIRI L1B HRIT data with CALIPSO layered data using
atrain_match functionalities.
atrain_match: https://github.com/foua-pps/atrain_match
Author: Daniel Philipp (DWD, 2021)
"""

import argparse
import glob
import os
import xarray as xr
import numpy as np
import satpy
from pyorbital import astronomy
from pyorbital.orbital import get_observer_look
import datetime

# atrain_match functionailities
from atrain_match.utils.common import Cross
from atrain_match.utils.runutils import read_config_info
from atrain_match.matchobject_io import write_truth_imager_match_obj
import atrain_match.config as config
from atrain_match.libs.truth_imager_match import (find_truth_files,
                                    get_additional_calipso_files_if_requested,
                                    get_calipso_matchups)

AM_PATHS, SETTINGS = read_config_info()

# constants
NANOSEC_TO_SEC = 1E-9
FULL_DISK_DIM = 3712
BANDS = ['VIS006', 'VIS008', 'IR_016', 'IR_039', 'WV_062', 'WV_073',
         'IR_087', 'IR_097', 'IR_108', 'IR_120', 'IR_134']


class HRITData(object):
    def __init__(self, array_dict=None):
        """ Class holding imager data as AllImagerData in atrain_match. """
        # time information
        self.sec1970_start = None
        self.sec1970_end = None
        self.time = None

        # instrument name
        self.instrument = None
        self.type = None

        # lat/lon
        self.latitude = None
        self.longitude = None

        # satellite and solar zenith angles
        self.satzen = None
        self.solzen = None
        self.satazi = None

        # sensor measurements
        self.vis006 = None
        self.vis008 = None
        self.ir_016 = None
        self.ir_039 = None
        self.ir_062 = None
        self.ir_073 = None
        self.ir_087 = None
        self.ir_097 = None
        self.ir_108 = None
        self.ir_120 = None
        self.ir_134 = None

        if array_dict is not None:
            self.__dict__.update(array_dict)


def get_satellite_angles(dataset, lons, lats):
    """Compute satellite angles.
    Returns:
        Satellite azimuth angle, Satellite zenith angle in degrees
    """
    sat_lon, sat_lat, sat_alt = satpy.utils.get_satpos(dataset)
    sat_alt *= 0.001

    # Compute angles
    sata, satel = get_observer_look(
        sat_lon,
        sat_lat,
        sat_alt,
        dataset.attrs['start_time'],
        lons, lats, 0)
    satz = 90 - satel

    return sata, satz


def parse_scenesfile_hrit(file):
    """
    Decompose epilogue filename to get information about imager/time.
    Exemplary epilogue filename pattern:
    H-000-MSG4__-MSG4________-_________-EPI______-201803011200-__
    """

    # extract satellite name (MSGX) and timestamp
    filename = os.path.basename(file)
    items = [x for x in filename.split('_') if x != '']
    satname = items[1][1:]
    time = items[-2][1:-1]

    # convert timestamp to datetime object
    date_time = datetime.datetime.strptime(time + '00', '%Y%m%d%H%M%S')

    values = {"satellite": satname.lower(),
              "date_time": date_time,
              "orbit": "99999",
              "date": date_time.strftime("%Y%m%d"),
              "year": date_time.year,
              "month": "%02d" % (date_time.month),
              "time": date_time.strftime("%H%M"),
              "ccifilename": filename,
              "ppsfilename": None}

    # basename for output filename
    values['basename'] = values["satellite"] + "_" + values["date"] + "_" + \
                         values["time"] + "_" + values["orbit"]

    return satname.lower(), time, values


def run_collocation(epi_file, values, odir):
    """ Calculate and write collocations of one slot. """
    truth_files = find_truth_files(date_time = values['date_time'],
                                   AM_PATHS=AM_PATHS,
                                   SETTINGS=SETTINGS,
                                   values=values,
                                   truth='calipso')
    if truth_files is None:
        truth_files = []
    print('Found ', len(truth_files), ' CALIPSO truth file(s).')

    if len(truth_files) > 0:

        extra_files = get_additional_calipso_files_if_requested(
                                            calipso_files=truth_files,
                                            SETTINGS=SETTINGS)
        calipso5km, calipso1km, calipso5km_aerosol = extra_files

        hrit_data, hrit_success = read_hrit(epi_file,
                                            utc_time=values['date_time']
                                            )

        if hrit_success:
            collocations = get_calipso_matchups(truth_files,
                                                values,
                                                hrit_data,
                                                AM_PATHS, SETTINGS,
                                                calipso1km, calipso5km,
                                                calipso5km_aerosol)
        else:
            collocations = None


        if collocations is not None and hrit_success:
            try:
                # save to disk as HDF5 matchup file
                match_file = os.path.join(odir, values['basename'] + '.h5')
                write_truth_imager_match_obj(match_file,
                                             collocations,
                                             SETTINGS,
                                             imager_obj_name='seviri_hrit')
                print('Saved: ', match_file)
            except AttributeError:
                print('No matches found. Skipping.')
        else:
            print('No caliop files in temporal range or HRIT read failed. '
                  'Skipping.')
    else:
        print('No caliop files found. Skipping.')




def read_hrit(epi_file, utc_time):
    """
    Read L1B HRIT data, calculate auxiliary data and put them in
    a structure atrain_match can handle.
    """
    # search for HRIT files in given slot
    directory = os.path.dirname(epi_file)
    hrit_files = glob.glob(os.path.join(directory, 'H-000-*'))

    hrit_success = True

    if len(hrit_files) < 114:
        hrit_success = False

    try:
        # read SEVIRI slot
        scene = satpy.Scene(reader="seviri_l1b_hrit",
                            filenames=hrit_files,
                            reader_kwargs={'calib_mode': 'GSICS'})
    except:
        hrit_success = False
        print('Unable to read HRIT for this slot!')

    if hrit_success:

        scene.load(BANDS)
        # get lon/lat
        lon, lat = scene['VIS006'].attrs['area'].get_lonlats()
        lon[np.fabs(lon) > 360] = np.nan
        lat[np.fabs(lat) > 90] = np.nan

        # calculate satellite zenith and satellite azimuth on the fly
        satazi, satzen = get_satellite_angles(scene['VIS006'], lon, lat)

        # get acquisition time of each pixel (scanline)
        acq_time_1d = scene['VIS006'].acq_time.values
        datetime_1970 = datetime.datetime(1970, 1, 1, 0, 0, 0)
        acq_time_2d = np.ones((FULL_DISK_DIM, FULL_DISK_DIM)) * np.nan
        for t in range(FULL_DISK_DIM):
            if not np.isnat(acq_time_1d[t]):
                element = acq_time_1d[t].tolist() * NANOSEC_TO_SEC
                a = datetime.datetime.utcfromtimestamp(element)
                acq_time_2d[:, t] = (a - datetime_1970).total_seconds()

        # add data to atrain_match structure
        hrit = HRITData()
        hrit.vis006 = scene['VIS006'].values
        hrit.vis008 = scene['VIS008'].values
        hrit.ir_016 = scene['IR_016'].values
        hrit.ir_039 = scene['IR_039'].values
        hrit.ir_062 = scene['WV_062'].values
        hrit.ir_073 = scene['WV_073'].values
        hrit.ir_087 = scene['IR_087'].values
        hrit.ir_097 = scene['IR_097'].values
        hrit.ir_108 = scene['IR_108'].values
        hrit.ir_120 = scene['IR_120'].values
        hrit.ir_134 = scene['IR_134'].values
        hrit.longitude = lon
        hrit.latitude = lat
        hrit.time = acq_time_2d
        hrit.sec1970_start = np.nanmin(acq_time_2d)
        hrit.sec1970_end = np.nanmax(acq_time_2d)
        hrit.solzen = astronomy.sun_zenith_angle(utc_time, lon, lat)
        hrit.satzen = satzen
        hrit.satazi = satazi
        hrit.instrument = 'seviri'
        hrit.type = 'hrit'

        return hrit, hrit_success

    return None, hrit_success


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--hrit_epilogue_filelist', required=True)
    parser.add_argument('--odir', required=True)
    args = parser.parse_args()

    odir = args.odir
    epilogue_list = args.hrit_epilogue_filelist

    print('')

    # iterate through epilogue filelist (one epi file for each slot)
    with open(epilogue_list, 'r') as fh:
        for line in fh:
            if line.rstrip() in "":
                continue
            else:
                satname, time, values = parse_scenesfile_hrit(line)
                print('--------------------------------------------')
                print('Collocating {}'.format(time))
                print('--------------------------------------------\n')
                try:
                    run_collocation(line, values, odir)
                except:
                    print('ERROR IN {}! SKIPPING.'.format(line))

    print('#########################################')
    print('FINISHED SUCCESSFULLY')
    print('#########################################')
