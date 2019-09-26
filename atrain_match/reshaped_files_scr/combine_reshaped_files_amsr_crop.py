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
from matchobject_io import (readAmsrImagerMatchObj,
                            writeAmsrImagerMatchObj,
                            AmsrImagerTrackObject)
from my_dir import ADIR

instrument = "gac"
satellite = "noaa18"
truth = "amsr"
version = "v2014"

BASE_DIR = ADIR + "/DATA_MISC/reshaped_files_validation_2018/global_{instrument}_{version}_created20180927/".format(instrument=instrument, version=version)
ROOT_DIR = BASE_DIR + "/Reshaped_Files/{satellite}/5km/*/*{truth}*.h5".format(
    satellite=satellite, truth=truth)
OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_{truth}_lwp/{satellite}/5km/".format(
    truth=truth,satellite=satellite)
outfile_template = "5km_{satellite}_00000000_0000_00000_{truth}_{instrument}_match.h5".format(
    satellite=satellite, truth=truth, instrument = instrument)

print ROOT_DIR


aObj = AmsrImagerTrackObject()

for year in [2010]:#2012/02","2012/05", "2012/08", "2013/07", "2014/02", "2014/04", "2014/09"]:
    #for month in ["01","02","03","04","05","06","07","08","09","10","11","12"]:
    #for month in ["06"]:
    for month in ["01"]:#, "02","03","04","05","06","07","08","09","10","11","12"]:
        for day in ["01"]:#, "14"]:
            OUT_DIR = OUT_DIR_TEMPLATE
            if not os.path.exists(OUT_DIR):
                os.makedirs(OUT_DIR)

            files = glob(ROOT_DIR)
            if len(files) == 0:
                continue
            num_n = 0
            for filename in files:
                print  os.path.basename(filename)
                try:
                    aObj_new=readAmsrImagerMatchObj(filename)

                except:
                    print "problem with", os.path.basename(filename)
                    continue

                if aObj_new.imager.all_arrays["cpp_lwp"] is None:
                    print "problem with lwp is None", os.path.basename(filename)
                    continue

                from utils.validate_lwp_util import get_lwp_diff_inner
                diff, pps_lwp, amsr_lwp, selection = get_lwp_diff_inner(aObj_new, True)


                if aObj.diff_sec_1970 is not None and len(diff) > 0:
                    aObj.diff_sec_1970 = np.concatenate([aObj.diff_sec_1970, aObj_new.diff_sec_1970[selection]], 0)
                    aObj.imager.cpp_lwp =  np.concatenate([aObj.imager.cpp_lwp, pps_lwp.ravel()], 0)
                    aObj.amsr.lwp = np.concatenate([aObj.amsr.lwp, np.array(amsr_lwp).ravel()],0)
                    for name in ["satz", "sunz", "longitude", "latitude"]:
                        aObj.imager.all_arrays[name] = np.concatenate([aObj.imager.all_arrays[name],
                                                                      aObj_new.imager.all_arrays[name][selection]],
                                                                     0)
                elif len(diff) > 0:
                    aObj.diff_sec_1970 = aObj_new.diff_sec_1970[selection]
                    aObj.imager.cpp_lwp =  pps_lwp.ravel()
                    aObj.amsr.lwp = np.array(amsr_lwp).ravel()
                    for name in ["satz", "sunz", "longitude", "latitude"]:
                        aObj.imager.all_arrays[name] = aObj_new.imager.all_arrays[name][selection]


                num_n +=1
            if num_n > 0:
                filename_night = outfile_template
                outfile = os.path.join(OUT_DIR, filename_night)
                writeAmsrImagerMatchObj(
                    outfile, aObj, imager_obj_name = 'pps')
                aObj = AmsrImagerTrackObject()

