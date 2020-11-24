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

# Author(s)
# Nina Hakansson

import numpy as np
from glob import glob
import os
from atrain_match.matchobject_io import (read_truth_imager_match_obj, 
                                         write_truth_imager_match_obj)

import time

instrument = "avhrr"
truth = "calipso"
version = "v2018"

SATELLITES = ["noaa18", "noaa19", "noaa17", "metopa", "metopb"]
YEAR_LIST = ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015']
BASE_DIR = '/nobackup/smhid14/sm_erjoh/pps/Atrain_match/patchCMSAF_Oct2020_qualityflagcheck/Reshaped_Files'
# BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files_from_kg_temp"
MAKE_EXTRA_CHECK = False

SETTINGS = {"WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE": False}


def remove_doubles(mObj, mObj2):
    non_doubles = [ind for ind in range(len(mObj.calipso.sec_1970))
                   if mObj.calipso.sec_1970[ind] not in mObj2.calipso.sec_1970]

    if MAKE_EXTRA_CHECK:
        # Takes a lot of time ....
        doubles_time_and_id = [(mObj.calipso.sec_1970[ind], mObj.calipso.profile_id[ind, 0])
                               for ind in range(len(mObj.calipso.sec_1970)) if ind not in non_doubles]

        extra_check = [time_and_id for time_and_id in doubles_time_and_id if time_and_id in zip(
            mObj2.calipso.sec_1970, mObj2.calipso.profile_id[:, 0])]
        if len(extra_check) != len(mObj.calipso.sec_1970) - len(non_doubles):
            print("Some points identified as doubles does not have "
                  "the same time and profile id as ther corresponding double"
                  "There is something wrong (in the code) here!")
            raise ValueError

    print("Found {:d} doubles (out of {:d}) in the new object".format(
        len(mObj.calipso.sec_1970) - len(non_doubles), len(mObj.calipso.sec_1970)))
    
    #: All are doublets
    if len(non_doubles) == 0:
        return -1
    else:
        mObj = mObj.extract_elements(idx=np.array(non_doubles))
        return mObj
    # import pdb;pdb.set_trace()
        

if __name__ == '__main__':
    '''
    Creates monthly merged reshape files.
    BASE_DIR is the dir where original reshape files are stored
    Merged results will be in BASE_DIR + _merged_{truth}/{satellite}/
    '''
    tic = time.time()
    files_cant_read = []
    files_all_doubles = []
    match_calipso_merged = None
    for satellite in SATELLITES:
        ROOT_DIR = BASE_DIR + "/{satellite}/5km/%s/%s/*%s%s*_*{truth}*.h5".format(
            satellite=satellite, truth=truth)
        print(ROOT_DIR)
        OUT_DIR_TEMPLATE = BASE_DIR + "_merged_{truth}/{satellite}/".format(
            truth=truth, satellite=satellite)
        outfile_template = "5km_{satellite}_%s%s01_0000_99999_{truth}_{instrument}_match.h5".format(
            satellite=satellite, truth=truth, instrument=instrument)
    
        for year in YEAR_LIST:
            for month in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
                num_n = 0
                files_all_doubles_month = []
                files = sorted(glob(ROOT_DIR % (year, month,
                                                year, month)))
                if len(files) == 0:
                    continue
                print('\n Year = %s, mon = %s \n' %(year, month))
                OUT_DIR = OUT_DIR_TEMPLATE
                if not os.path.exists(OUT_DIR):
                    os.makedirs(OUT_DIR)
                for filename in files:
                    print(os.path.basename(filename))
                    # match_calipso_new=read_truth_imager_match_obj(filename, truth=truth)
                    try:
                        match_calipso_new = read_truth_imager_match_obj(filename, truth=truth)
                    except:
                        print("problem with", os.path.basename(filename))
                        files_cant_read.append(filename)
                        continue
                        # if (match_calipso_new.cloudsat.RVOD_CWC_status is None or
                        # len(match_calipso_new.cloudsat.RVOD_CWC_status) != len(match_calipso_new.avhrr.cpp_lwp)):
                        #   print("Missing RVOD_CWC_status")
                        #   continue
                        
                    print("reading", os.path.basename(filename))
                    if match_calipso_merged is None:
                        match_calipso_merged = match_calipso_new
                    else:
                        match_calipso_new = remove_doubles(match_calipso_new, match_calipso_merged)
                        if match_calipso_new == -1:
                            files_all_doubles.append(filename)
                            files_all_doubles_month.append(filename)
                            continue
                        match_calipso_merged = match_calipso_merged + match_calipso_new
                    num_n += 1
    
                if num_n > 0:
                    filename_merged = outfile_template % (year, month)
                    outfile = os.path.join(OUT_DIR, filename_merged)
    
                    print(len(match_calipso_merged.calipso.sec_1970))
                    write_truth_imager_match_obj(
                        outfile, match_calipso_merged, SETTINGS, imager_obj_name='pps')
                    match_calipso_merged = None
                    print("File with all doublets %s files" %len(files_all_doubles_month))
    
    print("Could not read %s files" %len(files_cant_read))
    print("File with all doublets %s files, total" %len(files_all_doubles))
    print(time.time()-tic)
    
