#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read all COT profiles of Probability of detecting a cloud layer with a fixed COT: POD(cloudy_layer,cot). Then find the COT value for each grid point where POD first exceeds 50 % and plot the resulting map.
"""

import os,sys
import pdb
import h5py
from glob import glob
import numpy as np
import remap_module_KG
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util_KG import (read_POD_CLOUDY_LAYER_afternoon, read_POD_CLOUDY_LAYER_morning)

layer_cots = [0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.55,0.65,0.75,0.85,0.95,1.5,2.5,3.5,4.5]
layer_cot_names = ["0_025","0_075","0_125","0_175","0_225","0_275","0_325","0_375","0_425","0_475","0_55","0_65","0_75","0_85","0_95","1_50","2_50","3_50","4_50"]

N_layers = 19

#BASE_PLOT_DIR = "/nobackup/smhid13/sm_kgkar/atrain_match_CCI_V3/plot_stats/figs_stats_global_cotranges_TCC_POD_final"
BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD_final_night_600km"
BASE_PLOT_DIR_MORNING = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD_morning"

# AFTERNOON SATELLITES
#*********************

# Read POD profiles

filename = BASE_PLOT_DIR + "/pod_cloudy_layer_results.h5"
POD_CLOUDY_LAYER_AFTERNOON, lats, lons, radius_km = read_POD_CLOUDY_LAYER_afternoon(filename)
grid_shape = POD_CLOUDY_LAYER_AFTERNOON['0_025'].size
detection_limit = np.zeros(grid_shape)

for cot in range(N_layers):
    detection_limit = np.where(np.logical_and(np.greater(POD_CLOUDY_LAYER_AFTERNOON[layer_cot_names[cot]], 0.5),np.less(detection_limit, 0.02)),layer_cots[cot], detection_limit)


## # Write results to an hdf5 file

COMPRESS_LVL = 6

detection_limit_filename = BASE_PLOT_DIR + "/detection_limits_layer_cots.h5"

with h5py.File(detection_limit_filename, 'w') as f:
    f.create_dataset('detection_limit', data=detection_limit,
                     compression=COMPRESS_LVL)
    f.create_dataset('latitude', data=lats,
                     compression=COMPRESS_LVL)
    f.create_dataset('longitude', data=lons,
                     compression=COMPRESS_LVL)
    f.create_dataset('fibonacci_radius_km', data=radius_km)
    f.close()

remap_module_KG.remap_and_plot_score_on_several_areas('afternoon_satellites', detection_limit, lats, lons, radius_km, 0.0, 5.0)

sys.exit()

# MORNING SATELLITES
#*******************

# Read POD profiles

filename = BASE_PLOT_DIR_MORNING + "/pod_cloudy_layer_results.h5"
POD_CLOUDY_LAYER_MORNING, lats, lons, radius_km = read_POD_CLOUDY_LAYER_morning(filename)
grid_shape = POD_CLOUDY_LAYER_MORNING['0_025'].size
detection_limit = np.zeros(grid_shape)

for cot in range(N_layers):
    detection_limit = np.where(np.logical_and(np.greater(POD_CLOUDY_LAYER_MORNING[layer_cot_names[cot]], 0.5),np.less(detection_limit, 0.02)),layer_cots[cot], detection_limit)


## # Write results to an hdf5 file

COMPRESS_LVL = 6

detection_limit_filename = BASE_PLOT_DIR_MORNING + "/detection_limits_layer_cots.h5"

with h5py.File(detection_limit_filename, 'w') as f:
    f.create_dataset('detection_limit', data=detection_limit,
                     compression=COMPRESS_LVL)
    f.create_dataset('latitude', data=lats,
                     compression=COMPRESS_LVL)
    f.create_dataset('longitude', data=lons,
                     compression=COMPRESS_LVL)
    f.create_dataset('fibonacci_radius_km', data=radius_km)
    f.close()

remap_module_KG.remap_and_plot_score_on_several_areas('morning_satellites', detection_limit, lats, lons, radius_km, 0.0, 5.0)


