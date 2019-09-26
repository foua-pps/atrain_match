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
import os
from glob import glob
import re
import numpy as np
from scipy import ndimage
from matchobject_io import read_truth_imager_match_obj
from plotting.along_track_plotting import (plot_cal_clsat_geoprof_imager,
                                           drawCalPPSHeightPlot_PrototypePPSHeight  )
from my_dir import ADIR
ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files/"
            "global_modis_14th_created20161108/")
ROOT_DIR_GAC_nn = (ADIR + "/DATA_MISC/reshaped_files/"
                   "ATRAIN_RESULTS_GAC_nn21/Reshaped_Files/noaa18/")
ROOT_DIR_GAC_old = (ADIR + "/DATA_MISC/reshaped_files/"
                    "ATRAIN_RESULTS_GAC_old/Reshaped_Files/noaa18/")
ROOT_DIR_GAC_old = (ADIR + "/DATA_MISC/reshaped_files/"
                    "ATRAIN_RESULTS_GAC_20161222/Reshaped_Files/noaa18/")

ROOT_DIR_GAC_nn_new = (ADIR + "/DATA_MISC/reshaped_files/"
                       "ATRAIN_RESULTS_GAC_nn20161125/Reshaped_Files/noaa18/")
ROOT_DIR_GAC_nn_imager = (ADIR + "/DATA_MISC/reshaped_files/"
                         "ATRAIN_RESULTS_GAC_nnimager_20161202/Reshaped_Files/noaa18/")

BASE_DIR = ADIR + "/DATA_MISC/reshaped_files_validation_2018/"
ROOT_DIR_v2014_GAC = (BASE_DIR + "global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/")
ROOT_DIR_v2018_GAC = (BASE_DIR + "global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/")
re_name = re.compile("_global_(\w+_\w+_\w+)\/")


files = glob(ROOT_DIR_v2018_GAC + "/200*/*cali*h5")
# for ROOT_DIR in [ROOT_DIR_GAC_nn_imager, ROOT_DIR_GAC_nn, ROOT_DIR_GAC_nn_new]:
for filename in files:
    match = re_name.search(ROOT_DIR)
    name = "no_name"
    if match:
        name = match.group(1)
    basename = os.path.basename(filename)
    match_calipso =  read_truth_imager_match_obj(filename)
    match_calipso_OLD =  read_truth_imager_match_obj(filename.replace(ROOT_DIR_v2018_GAC, ROOT_DIR_v2014_GAC))
    height_pps =  match_calipso.imager.all_arrays['ctth_height']
    height_pps_old =  match_calipso_OLD.imager.all_arrays['ctth_height']
    try:
        use = np.logical_or(height_pps > 0, height_pps_old > 0)
    except:
        continue
    # use = np.logical_and(np.abs(match_calipso.imager.all_arrays['latitude'])>70, use)
    xmin = [i for i, x in enumerate(match_calipso.imager.all_arrays['latitude']) if x<-50]
    try:
        xmin = xmin[0]
    except:
        xmin = 0
    drawCalPPSHeightPlot_PrototypePPSHeight(match_calipso.calipso,
                                            use,
                                            height_pps_old + match_calipso.calipso.all_arrays['elevation'],
                                            height_pps + match_calipso.calipso.all_arrays['elevation'],
                                            ADIR + "/PICTURES_FROM_PYTHON/CTTH_LAPSE_RATE_INVESTIGATION/",
                                            "test_plot_file_part_nn_ctth_%s_%s"%(name, basename.split('.h5')[0]),
                                            file_type='png',
                                            xmin=xmin,
                                            xmax=xmin + 2500,
                                            instrument='imager')#, MAXHEIGHT = 18000)

    # plt.show()
    # plt.close("all")
