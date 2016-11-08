
"""Read all matched data and make tables
"""
import os
from glob import glob
import numpy as np
from scipy import ndimage
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)

def get_clear_cloudy_vectors(caObj, use):#, both_have_ct):
    nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
    isCalipsoCloudy = np.logical_and(
        nlay > 0, 
        caObj.calipso.all_arrays['cloud_fraction']>0.5)
    isCalipsoCloudy = np.logical_and(
        isCalipsoCloudy, 
        caObj.calipso.all_arrays['total_optical_depth_5km']>0.15)
    isCalipsoClear = np.logical_and(nlay == 0, meancl<0.01)
    isCalipsoClear = np.logical_and(
        isCalipsoClear, 
        caObj.calipso.all_arrays['total_optical_depth_5km']<0)
    isCloudyPPSorMAIA = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                                       caObj.avhrr.all_arrays['cloudtype']<21) 
    isClearPPSorMAIA = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                      caObj.avhrr.all_arrays['cloudtype']<5)
    nodata = np.sum(caObj.avhrr.all_arrays['cloudtype'][use]>200)
    use = np.logical_and(use, np.logical_or(isCloudyPPSorMAIA, isClearPPSorMAIA))
    use = np.logical_and(use, np.logical_or(isCalipsoCloudy, isCalipsoClear))
    #use = np.logical_and(use, both_have_ct)
    #print isCalipsoCloudy[use].all()
    #print np.sum(isCalipsoCloudy[use])
    #print np.sum(np.logical_and(isCalipsoCloudy[use], 
    #                            isCloudyPPSorMAIA[use]))

    PODcloudy = (np.sum(np.logical_and(isCalipsoCloudy[use], 
                                       isCloudyPPSorMAIA[use]))*1.0
                 /np.sum(isCalipsoCloudy[use]))
    FARcloudy = (np.sum(np.logical_and(isCalipsoClear[use], 
                                       isCloudyPPSorMAIA[use]))*1.0
                 /np.sum(isCloudyPPSorMAIA[use]))
    PODclear = (np.sum(np.logical_and(isCalipsoClear[use], 
                                       isClearPPSorMAIA[use]))*1.0
                 /np.sum(isCalipsoClear[use]))
    Num = np.sum(use)
    part_nodata = nodata*1.0/(nodata+Num)
    print "N: %d POD-cloudy: %3.2f FAR-cloudy: %3.2f POD-clear %3.2f"%(# Part-CT-nodata: %3.2f"%(
        Num, PODcloudy, FARcloudy , PODclear)#, part_nodata)
                 

def print_common_stats(caObjMAIA, caObjPPS):
    pps_sec_1970 = caObjPPS.calipso.all_arrays['sec_1970']
    maia_sec_1970 = caObjMAIA.calipso.all_arrays['sec_1970']
    pps_lat = caObjPPS.calipso.all_arrays['latitude']
    maia_lat = caObjMAIA.calipso.all_arrays['latitude']
    min_time_pps =  np.min(pps_sec_1970)
    max_time_pps =  np.max(pps_sec_1970)
    min_time_maia =  np.min(maia_sec_1970) 
    max_time_maia =  np.max(maia_sec_1970)
    min_lat_pps =  np.min(pps_lat)
    max_lat_pps =  np.max(pps_lat)
    min_lat_maia =  np.min(maia_lat) 
    max_lat_maia =  np.max(maia_lat)

    maia_profile_id = caObjMAIA.calipso.profile_id
    pps_profile_id = caObjPPS.calipso.profile_id

    #maia_ct = caObjMAIA.avhrr.all_arrays['cloudtype']
    pps_bt11 = caObjPPS.avhrr.all_arrays['bt11micron']
    pps_profile_id[pps_bt11<0]=-9
    
    #profile_id_in_both = np.intersect1d(maia_profile_id,pps_profile_id) 

    use_pps = pps_sec_1970>0
    use_maia = maia_sec_1970>0
    use_pps[pps_sec_1970<min_time_maia] = False
    use_pps[pps_sec_1970>max_time_maia] = False
    use_maia[maia_sec_1970<min_time_pps] = False
    use_maia[maia_sec_1970>max_time_pps] = False
    use_pps[pps_lat<min_lat_maia] = False
    use_pps[pps_lat>max_lat_maia] = False
    use_maia[maia_lat<min_lat_pps] = False
    use_maia[maia_lat>max_lat_pps] = False

    use_maia_same_profile = np.array([p_id in pps_profile_id[use_pps] for p_id in maia_profile_id])
    use_pps_same_profile =  np.array([p_id in maia_profile_id[use_maia] for p_id in pps_profile_id])
    use_maia = np.logical_and(use_maia, use_maia_same_profile)
    use_pps = np.logical_and(use_pps, use_pps_same_profile)


    #both_have_ct = np.logical_and(caObjMAIA.avhrr.all_arrays['cloudtype']<20,
    #                              caObjPPS.avhrr.all_arrays['cloudtype']<200)

    print "MAIA |",
    get_clear_cloudy_vectors(caObjMAIA, use_maia)
    print "PPS  |",
    get_clear_cloudy_vectors(caObjPPS, use_pps)#, both_have_ct)


MAIA_ROOT_DIR = ("/home/a001865/VALIDATIONS_TEST/"
                 "VIIRS_npp_pps_for_maia_comparison/Reshaped_Files/"
                 "npp/1km/2012/")
PPS_ROOT_DIR = ("/home/a001865/VALIDATIONS_TEST/"
                "VIIRS_npp_pps_for_maia_comparison_the_pps/Reshaped_Files/"
                "npp/1km/2012/")

for month in ['06','07','08','10']:
    maia_files = glob(MAIA_ROOT_DIR + "%s/*/*h5"%(month))
    pps_files = glob(PPS_ROOT_DIR + "%s/*/*h5"%(month))
    caObjMAIA = CalipsoAvhrrTrackObject()
    #caObjPPS = CalipsoAvhrrTrackObject()
    for filename in maia_files:
     caObjMAIA +=  readCaliopAvhrrMatchObj(filename)  
    for filename in pps_files:
        caObjPPS =  readCaliopAvhrrMatchObj(filename) 
        print "Scene %s"%(os.path.basename(filename))
        print_common_stats(caObjMAIA, caObjPPS)
