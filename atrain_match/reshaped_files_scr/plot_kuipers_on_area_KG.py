#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Read all matched data and make some plotting
"""
import os
import pdb
from glob import glob
import numpy as np
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util_KG import (PerformancePlottingObject,
                                          ppsMatch_Imager_CalipsoObject,write_filtered_TCC_POD_with_FAR_with_N,write_filtered_TCC_POD_without_FAR_without_N)
cots = [0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.0,2.0,3.0,4.0,5.0]

chosen_cot = 5.0
isGAC_CCI = False
isGAC_CCI_morning = False
isModis1km = False
isModis1km_lvl2 = False
isModis1km_nnctth = False
isNPP_v2014 = False
isGAC_v2014_morning_sat = False
isGAC_v2014 = True
isGAC_CCI_V3 = False
cci_orbits = False
method = 'KG' #Nina or KG or BASIC==no filter
DNT="night" #"all/day/night/twilight"
filter_method = 'no' #no or satz
radius_km = 600 #t.ex 75 250 500
## BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figures_globalplots_morning_all"
#BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges"
#BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_AMT_paper_filtered"
#BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_AMT_paper_unfiltered"
#BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD"
#BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD_final"
BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD_final_night_600km"

pod_filename = BASE_PLOT_DIR + "/stats_filtered_5_00.h5"

onlyCirrus=False
isACPGv2012=False
if isGAC_CCI:
    num_files_to_read = 90
    isGAC=True
    satellites = "cci_noaa18_noaa19"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/cci_reshaped_tgz/"
    files = glob(ROOT_DIR + "noaa19/*/????/??/*/*.h5")
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/cci_reshaped_tgz/"
    files = files + glob(ROOT_DIR + "noaa18/*/????/??/*/*.h5")
elif isGAC_CCI_morning:
    num_files_to_read = 90
    isGAC=True
    satellites = "cci_noaa17_metopa"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/cci_reshaped_tgz/"
    files = glob(ROOT_DIR + "metopa/*/????/??/*/*.h5")
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/cci_reshaped_tgz/"
    files = files + glob(ROOT_DIR + "noaa17/*/????/??/*/*.h5")
elif isModis1km:
    num_files_to_read = 1
    isGAC=False
    satellites = "eos_modis_v2014_2"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20160615/"
    #files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/12/2010??14_*/*h5")
    files = glob(ROOT_DIR + "Reshaped_Files/merged/*.h5")
elif isModis1km_lvl2:
    num_files_to_read = 1
    isGAC=False
    satellites = "eos_modis_lvl2_C6"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20170324/"
    #files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/12/2010??14_*/*h5")
    files = glob(ROOT_DIR + "Reshaped_Files_merged/eos2/1km/2010/*/*.h5")
elif isModis1km_nnctth:
    num_files_to_read = 1
    isGAC=False
    satellites = "eos_modis_v2018"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20170324/"
    #files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/12/2010??14_*/*h5")
    files = glob(ROOT_DIR + "Reshaped_Files_merged/eos2/1km/2010/*/*.h5")
elif isNPP_v2014:
    num_files_to_read = 30
    isGAC=False
    satellites = "npp"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/sh_reshaped_patch_2014/"
    files = glob(ROOT_DIR + "Reshaped_Files/npp/1km/????/06/arc*/*h5")
elif isGAC_v2014_morning_sat:
    num_files_to_read = 20 #1#30*3
    isGAC=True
    ROOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/Reshaped_Files/"
    #ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    #files = glob(ROOT_DIR + "noaa17/5km/20??/??/*/*h5")
    #files = files + glob(ROOT_DIR + "metop*/5km/20??/??/*/*h5")
    files = glob(ROOT_DIR + "metopa/5km/*/*/*/*h5")
    files = files + glob(ROOT_DIR + "metopb/5km/*/*/*/*h5")
    files = files + glob(ROOT_DIR + "noaa17/5km/*/*/*/*h5") 
    satellites = "metopa_metopb_noaa17"
    if cci_orbits:
        satellites = "part_metopa_noaa17"
        files = glob(ROOT_DIR + "merged/noaa17*2006*h5")
        files = files + glob(ROOT_DIR + "merged/noaa17*2007*h5")
        files = files + glob(ROOT_DIR + "merged/metopa*20*0*h5") #07.08.09.10
elif isGAC_v2014:
    num_files_to_read = 20
    isGAC=True
    satellites = "nooa18_nooaa19"
    ROOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/Reshaped_Files/"
    #ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa18/5km/*/*/*/*.h5")
    #files = glob(ROOT_DIR + "noaa18/5km/*/*/*/*.h5")
    files = files + glob(ROOT_DIR + "noaa19/5km/*/*/*/*.h5")
    #pdb.set_trace()
    #files = glob(ROOT_DIR + "noaa18*h5")
    #files = files + glob(ROOT_DIR + "noaa19*h5")
    #files = glob(ROOT_DIR + "merged/noaa19*2014*h5")
    if cci_orbits:
        satellites = "part_noaa18_noaa19"
        files = glob(ROOT_DIR + "merged/noaa18*200*h5") #06,07,08,09
        files = files + glob(ROOT_DIR + "merged/noaa19*20*0*h5") #2009, 2010
elif isGAC_CCI_V3:
    num_files_to_read = 20
    isGAC=True
    satellites = "nooa18_nooaa19"
    ROOT_DIR = "/nobackup/smhid13/sm_kgkar/atrain_match_CCI_V3/Reshaped_Files/"
    #ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa18/5km/*/*/*/*.h5")
    #files = glob(ROOT_DIR + "noaa18/5km/*/*/*/*.h5")
    files = files + glob(ROOT_DIR + "noaa19/5km/*/*/*/*.h5")
    #pdb.set_trace()
    #files = glob(ROOT_DIR + "noaa18*h5")
    #files = files + glob(ROOT_DIR + "noaa19*h5")
    #files = glob(ROOT_DIR + "merged/noaa19*2014*h5")
    if cci_orbits:
        satellites = "part_noaa18_noaa19"
        files = glob(ROOT_DIR + "merged/noaa18*200*h5") #06,07,08,09
        files = files + glob(ROOT_DIR + "merged/noaa19*20*0*h5") #2009, 2010
    
        

pplot_obj = PerformancePlottingObject()
pplot_obj.flattice.set_flattice(radius_km=radius_km)
pplot_obj.flattice.PLOT_DIR = BASE_PLOT_DIR
pplot_obj.flattice.DNT = DNT
pplot_obj.flattice.satellites = satellites
pplot_obj.flattice.filter_method = filter_method
pplot_obj.flattice.cc_method = method
pplot_obj.flattice.cotfilt_value = chosen_cot # Testing adding extra parameter
pplot_obj.flattice.isGAC=isGAC


caObj = CalipsoAvhrrTrackObject()
temp_obj = ppsMatch_Imager_CalipsoObject()
temp_obj.DNT = pplot_obj.flattice.DNT
temp_obj.satellites = pplot_obj.flattice.satellites
temp_obj.filter_method = pplot_obj.flattice.filter_method
temp_obj.cc_method = pplot_obj.flattice.cc_method
temp_obj.cotfilt_value = pplot_obj.flattice.cotfilt_value
temp_obj.isGAC = pplot_obj.flattice.isGAC

num = 0
for filename in files:
    print  os.path.basename(filename)

    num +=1
    try :
        caObj_new=readCaliopAvhrrMatchObj(filename)#, var_to_skip='segment')
        if "modis_lvl2" in satellites:
            caObj_new.avhrr.all_arrays["cloudtype"] = np.where(
                caObj_new.modis.all_arrays["cloud_emissivity"]>100,1,7)
            caObj_new.avhrr.all_arrays["ctth_temperature"] = (
                caObj_new.modis.all_arrays["temperature"]-15000+150)       
            caObj_new.avhrr.all_arrays["ctth_height"] = (
                caObj_new.modis.all_arrays["height"])
        
         
    except:
        #raise
        print "skipping file %s"%(filename)
        continue
    if num_files_to_read==1:
        print "Get info from one file!"
        temp_obj.get_some_info_from_caobj(caObj_new)
        print "Got info, now remap to the lattice"
        pplot_obj.add_detection_stats_on_fib_lattice(temp_obj)     
    elif num >num_files_to_read:
        print "Get info from some %d files!"%(num_files_to_read)
        temp_obj.get_some_info_from_caobj(caObj)
        print "got info, now remap to the lattice"
        pplot_obj.add_detection_stats_on_fib_lattice(temp_obj) 
        print "Got info from some files!"
        caObj = caObj_new
        num=0
    else:
        caObj = caObj + caObj_new

#Get info from the last files too 
if num_files_to_read!=1:
    print "Get info from last files!"
    temp_obj.get_some_info_from_caobj(caObj)
    pplot_obj.add_detection_stats_on_fib_lattice(temp_obj) 

## pplot_obj.flattice.calculate_bias()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-25.0, vmax=25.0, score='Bias',  screen_out_valid=True)
## pplot_obj.flattice.calculate_kuipers()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0.0, score='Kuipers')
#Calcualte hitrate
pplot_obj.flattice.calculate_hitrate()
pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0.5, score='Hitrate')
pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0.0, vmax=10000.0, score='Number_of_matchups')

#start to calculate
if 'modis' in satellites:
    #needs surftemp
    pplot_obj.flattice.calculate_height_bias_lapse()
    pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='lapse_bias_low')
    pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-5000.0, vmax=5000.0, score='lapse_bias_high')
    pplot_obj.flattice.calculate_lapse_rate()
    pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-25.0, vmax=0.0, score='lapse_rate')
    #Calcualte scores #nneds r13
    pplot_obj.flattice.calculate_increased_hitrate()
    pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='increased_Hitrate', vmin=-0.05, vmax=0.05)

## pplot_obj.flattice.calculate_height_mae()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=000.0, vmax=2000.0, score='ctth_mae_low')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=000.0, vmax=5000.0, score='ctth_mae_high')
## pplot_obj.flattice.calculate_height_bias()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='ctth_bias_low')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-5000.0, vmax=5000.0, score='ctth_bias_high')

## pplot_obj.flattice.calculate_temperature_bias()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-10.0, vmax=10.0, score='ctth_bias_temperature_low')
## pplot_obj.flattice.calculate_temperature_bias_t11()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-10.0, vmax=10.0, score='ctth_bias_temperature_low_t11')


## pplot_obj.flattice.calculate_height_bias_type()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='ctth_bias_type_0')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='ctth_bias_type_1')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='ctth_bias_type_2')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='ctth_bias_type_3')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='ctth_bias_type_4')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-2000.0, vmax=2000.0, score='ctth_bias_type_5')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-4000.0, vmax=4000.0, score='ctth_bias_type_6')
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=-4000.0, vmax=4000.0, score='ctth_bias_type_7')






#Calcualte threat_score
## pplot_obj.flattice.calculate_threat_score()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='Threat_Score')
## #Calcualte threat_score_clear
## pplot_obj.flattice.calculate_threat_score_clear()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='Threat_Score_Clear')
pplot_obj.flattice.calculate_calipso_cfc() #Needed for POD-profiles!
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0, vmax=100.0, score='calipso_cfc')
## pplot_obj.flattice.calculate_pps_cfc()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( vmin=0, vmax=100.0, score='pps_cfc')
pplot_obj.flattice.calculate_pod_cloudy() #Needed for POD-profiles!
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='PODcloudy', vmin=0.5)
## pplot_obj.flattice.calculate_pod_clear()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='PODclear', vmin=0.5)
pplot_obj.flattice.calculate_far_cloudy() #Needed as complement to POD-profiles
#pdb.set_trace()
write_filtered_TCC_POD_with_FAR_with_N(pod_filename, pplot_obj.flattice.calipso_cfc, pplot_obj.flattice.PODcloudy, pplot_obj.flattice.FARcloudy, pplot_obj.flattice.N, pplot_obj.flattice.lats, pplot_obj.flattice.lons, radius_km)
#write_filtered_TCC_POD_with_FAR_with_N_withPOD_above5(pod_filename, pplot_obj.flattice.calipso_cfc, pplot_obj.flattice.PODcloudy, pplot_obj.flattice.FARcloudy, pplot_obj.flattice.N, pplot_obj.flattice.lats, pplot_obj.flattice.lons, radius_km)
#write_filtered_TCC_POD_without_FAR_without_N(pod_filename, pplot_obj.flattice.calipso_cfc, pplot_obj.flattice.PODcloudy, pplot_obj.flattice.lats, pplot_obj.flattice.lons, radius_km) #Used for POD-profiles!
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='FARcloudy', vmax=0.5)
## pplot_obj.flattice.calculate_far_clear()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='FARclear', vmax=0.5)
## pplot_obj.flattice.calculate_RMS()
## pplot_obj.flattice.remap_and_plot_score_on_several_areas( score='RMS',  vmax=50.0, screen_out_valid=True)
## pdb.set_trace()

## name = "fig_%s_ccm_%s_%sfilter_dnt_%s_r%skm_"%(
##     pplot_obj.flattice.satellites,
##     pplot_obj.flattice.cc_method,
##     pplot_obj.flattice.filter_method,
##     pplot_obj.flattice.DNT, 
##     pplot_obj.flattice.radius_km)
## print name
## for measure in pplot_obj.flattice.__dict__['all_arrays'].keys():
##     if 'polar' in measure:
##         print measure, "%3.2f"%(
##             pplot_obj.flattice.__dict__['all_arrays'][measure] )
