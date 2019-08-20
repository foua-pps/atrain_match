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

import numpy as np
from glob import glob
import os
from matchobject_io import (readTruthImagerMatchObj,
                            writeTruthImagerMatchObj)


instrument = "modis"
satellite = "eos2"
truth = "calipso"
version = "v2018"

BASE_DIR = "/my_path/global_{instrument}_{version}_created20181023/".format(instrument=instrument, version=version)
ROOT_DIR = BASE_DIR + "/Reshaped_Files/{satellite}/1km/2010/%s/*2010%s%s_*{truth}*.h5".format(
    satellite=satellite, truth=truth)
OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_{truth}/{satellite}/1km/2010/%s/".format(
    truth=truth,satellite=satellite)
outfile_template = "1km_{satellite}_2010%s%s_%s00_00000_{truth}_{instrument}_match.h5".format(
    satellite=satellite, truth=truth, instrument = instrument)

SETTINGS ={"WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE": False}

caObj_night = None
caObj_day = None
for year in [2010]:
    for month in ["02","03","04","05","06","07","08","09","10","11","12","01"]:
        for day in ["14", "01"]:
            OUT_DIR = OUT_DIR_TEMPLATE%(month)
            if not os.path.exists(OUT_DIR):
                os.makedirs(OUT_DIR)

            files = sorted(glob(ROOT_DIR%(month,month,day)))
            if len(files)==0:
                continue
            num_n = 0  
            num_d = 0
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
                if int(os.path.basename(filename).split('_')[3][0:2])<12:
                    #time between 0000 and 11:59                
                    num_n +=1
                    print("reading",os.path.basename(filename))
                    if caObj_night is None:
                        caObj_night =  caObj_new  
                    else:    
                        caObj_night = caObj_night + caObj_new             
                else :
                    num_d +=1
                    print("reading",os.path.basename(filename))
                    if caObj_day is None:
                        caObj_day =  caObj_new  
                    else:    
                        caObj_day = caObj_day + caObj_new
            if num_n>0:
                filename_night = outfile_template%(month,day,"00")
                outfile = os.path.join(OUT_DIR, filename_night)
                writeTruthImagerMatchObj(
                    outfile, caObj_night, SETTINGS, imager_obj_name = 'pps')   
                caObj_night = None

            if num_d>0:
                filename_day = outfile_template%(month,day,"12")
                outfile = os.path.join(OUT_DIR, filename_day)
                writeTruthImagerMatchObj(
                    outfile, caObj_day, SETTINGS, imager_obj_name = 'pps')   
                caObj_day = None 

