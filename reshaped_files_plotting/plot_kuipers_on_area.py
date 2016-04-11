#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Read all matched data and make some plotting
"""
import os
from glob import glob
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       get_some_info_from_caobj,
                                       get_detection_stats_on_area_map_grid)
isModis1km = False
isNPP_v2014 = False
isGAC_v2014_morning_sat = False
isGAC_v2014 = True


onlyCirrus=False
isACPGv2012=False
if isModis1km:
    num_files_to_read = 24*6
    isGAC=False
    figure_name = "figure_global_modis_14th_"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_created20160324/"
    files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/*/2010??14_*/*h5")
elif isNPP_v2014:
    num_files_to_read = 30
    isGAC=False
    figure_name = "figure_local_npp"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/sh_reshaped_patch_2014/"
    files = glob(ROOT_DIR + "Reshaped_Files/npp/1km/????/06/arc*/*h5")
elif isGAC_v2014_morning_sat:
    num_files_to_read = 30*3
    isGAC=True
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa17/5km/20??/??/*/*h5")
    files = files + glob(ROOT_DIR + "metop*/5km/20??/??/*/*h5")
    figure_name = "figure_morning_sat_"
elif isGAC_v2014:
    num_files_to_read = 30
    isGAC=True
    figure_name = "figure_"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa18/5km/20??/??/*/*noaa*h5")
    files = files + glob(ROOT_DIR + "noaa19/5km/20??/??/*/*noaa*h5")



#arctic_super_5010_test
#antarctica_test
#'ease_nh_test'
#my_area = 'ease_nh_test'
my_area = 'euro_arctic_test'
my_obj = PerformancePlottingObject()
my_obj.area.set_area(area_name=my_area)
caObj = CalipsoAvhrrTrackObject()

num = 0
for filename in files:
    #print  os.path.basename(filename)
    num +=1
    try :
        caObj_new=readCaliopAvhrrMatchObj(filename)        
    except:
        print "skipping file %s"%(filename)
        continue
    if num >num_files_to_read:
        print "Get info from some %d files!"%(num_files_to_read)
        my_obj = get_some_info_from_caobj(my_obj, caObj, isGAC=isGAC)
        my_obj = get_detection_stats_on_area_map_grid(my_obj)        
        caObj = caObj_new
        num=0
    else:
        caObj = caObj + caObj_new

#Get info from the last files
my_obj =get_some_info_from_caobj(my_obj, caObj, isGAC=isGAC)
my_obj = get_detection_stats_on_area_map_grid(my_obj)  

PLOT_DIR = "/home/a001865/Documents/ATRAIN_MATCH/plots8/"
#Calcualte kuipers
my_obj.area.calculate_kuipers()
my_obj.area.plot_kuipers(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
#Calcualte hitrate
my_obj.area.calculate_hitrate()
my_obj.area.plot_hitrate(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
#Calcualte threat_score
my_obj.area.calculate_threat_score()
my_obj.area.plot_threat_score(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
#Calcualte threat_score_clear
my_obj.area.calculate_threat_score_clear()
my_obj.area.plot_threat_score_clear(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
