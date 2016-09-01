#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Read all matched data and make some plotting
"""
import os
from glob import glob
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       get_some_info_from_caobj)
isModis1km = True
isNPP_v2014 = False
isGAC_v2014_morning_sat = False
isGAC_v2014 = True
method = 'KG'
DNT="all"

onlyCirrus=False
isACPGv2012=False
if isModis1km:
    num_files_to_read = 24*4
    isGAC=False
    figure_name = "figure_global_modis_14th_%s_dnt_%s_"%(method, DNT)
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20160615/"
    files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/??/2010??14_*/*h5")
elif isNPP_v2014:
    num_files_to_read = 30
    isGAC=False
    figure_name = "figure_local_npp_%s_dnt_%s_"%(method, DNT)
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/sh_reshaped_patch_2014/"
    files = glob(ROOT_DIR + "Reshaped_Files/npp/1km/????/06/arc*/*h5")
elif isGAC_v2014_morning_sat:
    num_files_to_read = 30*3
    isGAC=True
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa17/5km/20??/??/*/*h5")
    files = files + glob(ROOT_DIR + "metop*/5km/20??/??/*/*h5")
    figure_name = "figure_morning_sat_%s_dnt_%s_"%(method, DNT)
    #files = glob(ROOT_DIR + "noaa17/5km/20??/1*/*/*noaa*h5")
elif isGAC_v2014:
    num_files_to_read = 30
    isGAC=True
    figure_name = "figure_%s_dnt_%s_"%(method, DNT)
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa18/5km/20??/??/*/*noaa*h5")
    files = files + glob(ROOT_DIR + "noaa19/5km/20??/??/*/*noaa*h5")
    #files = glob(ROOT_DIR + "noaa19/5km/2010/09/*/*noaa*h5")


pplot_obj = PerformancePlottingObject()
pplot_obj.flattice.set_flattice(radius_km=75)
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
        my_obj = get_some_info_from_caobj(caObj, isGAC=isGAC, method=method, DNT=DNT)
        pplot_obj.add_detection_stats_on_fib_lattice(my_obj)        
        caObj = caObj_new
        num=0
    else:
        caObj = caObj + caObj_new

#Get info from the last files too
my_obj =get_some_info_from_caobj(caObj, isGAC=isGAC, method=method, DNT=DNT)
pplot_obj.add_detection_stats_on_fib_lattice(my_obj)  


pplot_obj.flattice.PLOT_DIR = "/home/a001865/ATRAIN_MATCH_KUIPERS_PLOT/"
pplot_obj.flattice.figure_name=figure_name
pplot_obj.flattice.calculate_lapse_rate()
pplot_obj.flattice.remap_and_plot_score_on_several_areas(
    vmin=-25.0, vmax=0.0, score='lapse_rate', screen_out_valid=True)
#Calcualte scores
#pplot_obj.flattice.calculate_increased_hitrate()
#pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='increased_Hitrate', vmin=-0.05, vmax=0.05)
pplot_obj.flattice.calculate_bias()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-25.0, vmax=25.0, score='Bias', 
                          screen_out_valid=True)
pplot_obj.flattice.calculate_kuipers()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0.0, score='Kuipers')
#Calcualte hitrate
pplot_obj.flattice.calculate_hitrate()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0.5, score='Hitrate')
#Calcualte threat_score
pplot_obj.flattice.calculate_threat_score()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='Threat_Score')
#Calcualte threat_score_clear
pplot_obj.flattice.calculate_threat_score_clear()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='Threat_Score_Clear')
pplot_obj.flattice.calculate_calipso_cfc()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0, vmax=100.0, score='calipso_cfc')
pplot_obj.flattice.calculate_pps_cfc()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0, vmax=100.0, score='pps_cfc')
pplot_obj.flattice.calculate_pod_cloudy()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='PODcloudy', vmin=0.5)
pplot_obj.flattice.calculate_pod_clear()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='PODclear', vmin=0.5)
pplot_obj.flattice.calculate_far_cloudy()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='FARcloudy', vmax=0.5)
pplot_obj.flattice.calculate_far_clear()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='FARclear', vmax=0.5)
pplot_obj.flattice.calculate_RMS()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='RMS',  vmax=50.0, screen_out_valid=True)
