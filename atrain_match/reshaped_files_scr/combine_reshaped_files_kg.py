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
from matchobject_io import (readTruthImagerMatchObj,
                            writeTruthImagerMatchObj)

import time

instrument = "avhrr"
truth = "calipso"
version = "v2018"

SATELLITES = ["noaa19"]
YEAR_LIST = ["2011"]
BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files_from_kg_temp"
MAKE_EXTRA_CHECK = False

SETTINGS ={"WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE": False}
tic = time.time()

def remove_doubles(mObj, mObj2):
    non_doubles = [ind for ind in range(len(mObj.calipso.sec_1970)) if mObj.calipso.sec_1970[ind] not in mObj2.calipso.sec_1970]

    if MAKE_EXTRA_CHECK:
        #Takes a lot of time .... 
        doubles_time_and_id = [(mObj.calipso.sec_1970[ind], mObj.calipso.profile_id[ind,0]) for ind  in range(len(mObj.calipso.sec_1970)) if ind not in  non_doubles]

        extra_check = [time_and_id for time_and_id in doubles_time_and_id if time_and_id in zip(mObj2.calipso.sec_1970, mObj2.calipso.profile_id[:,0])]
        if len(extra_check) != len(mObj.calipso.sec_1970)-len(non_doubles):
            print("Some points identified as doubles does not have "
                  "the same time and profile id as ther corresponding double"
                  "There is something wrong (in the code) here!")
            raise ValueError
            
                           
    print("Found {:d} doubles (out of {:d}) in the new object".format(
        len(mObj.calipso.sec_1970)-len(non_doubles), len(mObj.calipso.sec_1970)))
    #import pdb;pdb.set_trace()
    mObj = mObj.extract_elements(idx=np.array(non_doubles))
    return mObj   

caObj_merged = None
for satellite in SATELLITES:
    ROOT_DIR = BASE_DIR + "/{satellite}/5km/%s/%s/*%s%s*_*{truth}*.h5".format(
        satellite=satellite, truth=truth)
    print(ROOT_DIR)
    OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_{truth}/{satellite}/".format(
        truth=truth, satellite=satellite)
    outfile_template = "5km_{satellite}_%s%s01_0000_99999_{truth}_{instrument}_match.h5".format(
        satellite=satellite, truth=truth, instrument=instrument)

    for year in YEAR_LIST:
        for month in ["01", "02","03","04","05","06","07","08","09","10","11","12"]:
            OUT_DIR = OUT_DIR_TEMPLATE
            if not os.path.exists(OUT_DIR):
                os.makedirs(OUT_DIR)
            num_n = 0    
            files = sorted(glob(ROOT_DIR%(year, month,
                                          year, month)))
            if len(files)==0:
                continue
            for filename in files:
                print(os.path.basename(filename))
                #caObj_new=readTruthImagerMatchObj(filename, truth=truth) 
                try:
                    caObj_new=readTruthImagerMatchObj(filename, truth=truth)
                except:
                    print("problem with", os.path.basename(filename))
                    continue
                    #if caObj_new.cloudsat.RVOD_CWC_status is None or len(caObj_new.cloudsat.RVOD_CWC_status) != len(caObj_new.avhrr.cpp_lwp):
                    #    print("Missing RVOD_CWC_status")
                    #    continue
                num_n +=1
                print("reading",os.path.basename(filename))
                if caObj_merged is None:
                    caObj_merged =  caObj_new  
                else:    
                    caObj_new  = remove_doubles(caObj_new, caObj_merged)
                    caObj_merged = caObj_merged + caObj_new             
            if num_n>0:
                filename_merged = outfile_template%(year,month)
                outfile = os.path.join(OUT_DIR, filename_merged)
                
                print(len(caObj_merged.calipso.sec_1970))
                writeTruthImagerMatchObj(
                    outfile, caObj_merged, SETTINGS, imager_obj_name = 'pps')   
                caObj_merged = None


print(time.time()-tic)
