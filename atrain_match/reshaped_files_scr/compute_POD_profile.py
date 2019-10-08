#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read all filtered CALIPSO CFC, POD(cloudy) and FAR(cloudy) fields and compute COT profile of Probability of detecting a cloud layer with a fixed COT: POD(cloudy_layer,cot).
"""

import os
import pdb
import h5py
from glob import glob
import numpy as np
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util_KG import (PerformancePlottingObject,
                                          ppsMatch_Imager_CalipsoObject,write_filtered_TCC_POD_with_FAR_with_N,write_filtered_TCC_POD_without_FAR_without_N,read_filtered_TCC_POD_with_FAR_with_N, read_filtered_TCC_POD_without_FAR_without_N,write_POD_CLOUDY_LAYER)

N_layers = 19

# Don't use the following values on tot_matches! We need the value for each grid point, not the total ones!
#tot_matches_afternoon = 23305814 #Number from previous atrain_match validation - too many???/KG 2007-04-26
#tot_matches_morning = 1849721 #Number from previous atrain_match validation - too many???/KG 2007-04-26
#tot_matches_afternoon = 23191345 # summed self.N
#tot_matches_morning = 1844023 # summed self.N

cots = [0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.0,2.0,3.0,4.0,5.0]
cot_names = ["0_00","0_05","0_10","0_15","0_20","0_25","0_30","0_35","0_40","0_45","0_50","0_60","0_70","0_80","0_90","1_00","2_00","3_00","4_00","5_00"]
layer_cots = [0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.55,0.65,0.75,0.85,0.95,1.5,2.5,3.5,4.5]
layer_cot_names = ["0_025","0_075","0_125","0_175","0_225","0_275","0_325","0_375","0_425","0_475","0_55","0_65","0_75","0_85","0_95","1_50","2_50","3_50","4_50"]

BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD_final_night_600km"
#BASE_PLOT_DIR = "/nobackup/smhid13/sm_kgkar/atrain_match_CCI_V3/plot_stats/figs_stats_global_cotranges_TCC_POD_final"
BASE_PLOT_DIR_MORNING = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD_morning"

# AFTERNOON SATELLITES
#*********************

# Start with unfiltered data and read also FAR(cloudy)

calipso_cfc = {}
PODcloudy = {}
POD_CLOUDY_LAYER_AFTERNOON = {}

pod_filename = BASE_PLOT_DIR + "/stats_filtered_0_00.h5"

calipso_cfc[cot_names[0]], PODcloudy[cot_names[0]], FARcloudy, tot_matches_afternoon, lats, lons, radius_km = read_filtered_TCC_POD_with_FAR_with_N(pod_filename)

# Create POD(cloudy_layer,cot) for all layer cots!

for step in range(N_layers):
    print "Step= ", step
    pod_filename = BASE_PLOT_DIR + "/stats_filtered_%s.h5" % cot_names[step+1]
    calipso_cfc[cot_names[step+1]], PODcloudy[cot_names[step+1]], lats, lons, radius_km = read_filtered_TCC_POD_without_FAR_without_N(pod_filename)

    POD_CLOUDY_LAYER_AFTERNOON[layer_cot_names[step]] = (calipso_cfc[cot_names[step]]*tot_matches_afternoon*PODcloudy[cot_names[step]]-(calipso_cfc[cot_names[step]]*tot_matches_afternoon-(calipso_cfc[cot_names[step]]-calipso_cfc[cot_names[step+1]])*tot_matches_afternoon)*PODcloudy[cot_names[step+1]])/((calipso_cfc[cot_names[step]]-calipso_cfc[cot_names[step+1]])*tot_matches_afternoon)


# Write results to an hdf5 file

pod_cloudy_layer_filename = BASE_PLOT_DIR + "/pod_cloudy_layer_results.h5"
write_POD_CLOUDY_LAYER(pod_cloudy_layer_filename,POD_CLOUDY_LAYER_AFTERNOON, lats, lons, radius_km)

# Write also FAR values to an hdf5 file

COMPRESS_LVL = 6

far_filename = BASE_PLOT_DIR + "/false_alarm_rate_unfiltered.h5"

with h5py.File(far_filename, 'w') as f:
    f.create_dataset('false_alarm_rate', data=FARcloudy,
                     compression=COMPRESS_LVL)
    f.create_dataset('latitude', data=lats,
                     compression=COMPRESS_LVL)
    f.create_dataset('longitude', data=lons,
                     compression=COMPRESS_LVL)
    f.create_dataset('fibonacci_radius_km', data=radius_km)
    f.close()



# MORNING SATELLITES
#*******************

# Start with unfiltered data and read also FAR(cloudy)

## calipso_cfc = {}
## PODcloudy = {}
## POD_CLOUDY_LAYER_MORNING = {}

## pod_filename = BASE_PLOT_DIR_MORNING + "/stats_filtered_0_00.h5"

## calipso_cfc[cot_names[0]], PODcloudy[cot_names[0]], FARcloudy, tot_matches_morning, lats, lons, radius_km = read_filtered_TCC_POD_with_FAR_with_N(pod_filename)

## # Create POD(cloudy_layer,cot) for all layer cots!

## for step in range(N_layers):
##     print "Step= ", step
##     pod_filename = BASE_PLOT_DIR_MORNING + "/stats_filtered_%s.h5" % cot_names[step+1]
##     calipso_cfc[cot_names[step+1]], PODcloudy[cot_names[step+1]], lats, lons, radius_km = read_filtered_TCC_POD_without_FAR_without_N(pod_filename)
##     POD_CLOUDY_LAYER_MORNING[layer_cot_names[step]] = (calipso_cfc[cot_names[step]]*tot_matches_morning*PODcloudy[cot_names[step]]-(calipso_cfc[cot_names[step]]*tot_matches_morning-(calipso_cfc[cot_names[step]]-calipso_cfc[cot_names[step+1]])*tot_matches_morning)*PODcloudy[cot_names[step+1]])/((calipso_cfc[cot_names[step]]-calipso_cfc[cot_names[step+1]])*tot_matches_morning)
## # Write results to a hdf5 file

## pod_cloudy_layer_filename = BASE_PLOT_DIR_MORNING + "/pod_cloudy_layer_results.h5"
## write_POD_CLOUDY_LAYER(pod_cloudy_layer_filename,POD_CLOUDY_LAYER_MORNING, lats, lons, radius_km)

## # Write also FAR values to an hdf5 file

## COMPRESS_LVL = 6

## far_filename = BASE_PLOT_DIR_MORNING + "/false_alarm_rate_unfiltered.h5"

## with h5py.File(far_filename, 'w') as f:
##     f.create_dataset('false_alarm_rate', data=FARcloudy,
##                      compression=COMPRESS_LVL)
##     f.create_dataset('latitude', data=lats,
##                      compression=COMPRESS_LVL)
##     f.create_dataset('longitude', data=lons,
##                      compression=COMPRESS_LVL)
##     f.create_dataset('fibonacci_radius_km', data=radius_km)
##     f.close()

