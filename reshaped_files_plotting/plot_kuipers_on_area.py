#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Read all matched data and make some plottingxf
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
method = 'KG'

onlyCirrus=False
isACPGv2012=False
if isModis1km:
    num_files_to_read = 24*4
    isGAC=False
    figure_name = "figure_global_modis_14th_%s_"%(method)
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20160419/"
    files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/??/2010??14_*/*h5")
elif isNPP_v2014:
    num_files_to_read = 30
    isGAC=False
    figure_name = "figure_local_npp_%s_"%(method)
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/sh_reshaped_patch_2014/"
    files = glob(ROOT_DIR + "Reshaped_Files/npp/1km/????/06/arc*/*h5")
elif isGAC_v2014_morning_sat:
    num_files_to_read = 30*3
    isGAC=True
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa17/5km/20??/??/*/*h5")
    files = files + glob(ROOT_DIR + "metop*/5km/20??/??/*/*h5")
    figure_name = "figure_morning_sat_%s_"%(method)
    #files = glob(ROOT_DIR + "noaa17/5km/20??/1*/*/*noaa*h5")
elif isGAC_v2014:
    num_files_to_read = 30
    isGAC=True
    figure_name = "figure_%s_"%(method)
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa18/5km/20??/??/*/*noaa*h5")
    files = files + glob(ROOT_DIR + "noaa19/5km/20??/??/*/*noaa*h5")
    #files = glob(ROOT_DIR + "noaa19/5km/2010/09/*/*noaa*h5")

#area definition moved!!
#arctic_super_5010_test
#my_area = 'antarctica_test'
#my_area ='npole'
#'ease_nh_test'
#my_area = 'ease_nh_test'
#my_area = 'ibla_57n20e'
#my_area = 'cea5km_test'
#my_area = 'robinson_test'
#my_area = 'euro_arctic_test'
#my_area = 'ease_world_test'
my_obj = PerformancePlottingObject()
my_obj.area.set_area(radius_km=75)
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
        my_obj = get_some_info_from_caobj(my_obj, caObj, isGAC=isGAC, method=method)
        my_obj = get_detection_stats_on_area_map_grid(my_obj)        
        caObj = caObj_new
        num=0
    else:
        caObj = caObj + caObj_new

#Get info from the last files too
my_obj =get_some_info_from_caobj(my_obj, caObj, isGAC=isGAC, method=method)
my_obj = get_detection_stats_on_area_map_grid(my_obj)  

my_obj.area.PLOT_DIR = "/home/a001865/Documents/ATRAIN_MATCH/plots_modis/"
my_obj.area.figure_name=figure_name
#Calcualte kuipers
my_obj.area.calculate_kuipers()
my_obj.area._remap_score( vmin=0.0, score='Kuipers')
#my_obj.area.plot_kuipers(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
#Calcualte hitrate
my_obj.area.calculate_hitrate()
my_obj.area._remap_score( vmin=0.5, score='Hitrate')
#my_obj.area.plot_hitrate(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
#Calcualte threat_score
my_obj.area.calculate_threat_score()
my_obj.area._remap_score( score='Threat_Score')
#my_obj.area.plot_threat_score(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
#Calcualte threat_score_clear
my_obj.area.calculate_threat_score_clear()
my_obj.area._remap_score( score='Threat_Score_Clear')
#my_obj.area.plot_threat_score_clear(PLOT_DIR=PLOT_DIR, figure_name=figure_name)
my_obj.area.calculate_bias()
my_obj.area._remap_score( vmin=-25, vmax=25, score='Bias')
my_obj.area.calculate_calipso_cfc()
my_obj.area._remap_score( vmin=0, vmax=100.0, score='calipso_cfc')
my_obj.area.calculate_pod_cloudy()
my_obj.area._remap_score( score='PODcloudy',vmin=0.5)
my_obj.area.calculate_pod_clear()
my_obj.area._remap_score( score='PODclear',vmin=0.5)
my_obj.area.calculate_far_cloudy()
my_obj.area._remap_score( score='FARcloudy', vmax=0.5)
my_obj.area.calculate_far_clear()
my_obj.area._remap_score( score='FARclear', vmax=0.5)
