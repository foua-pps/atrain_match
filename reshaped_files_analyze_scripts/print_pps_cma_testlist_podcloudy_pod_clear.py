
"""Read all matched data and make tables
"""
import os
from glob import glob
import numpy as np
from scipy import ndimage
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from get_flag_info import get_pixels_where_test_is_passed

def print_common_stats(caObj, use, name_dict):
    nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
    isCalipsoCloudy = np.logical_and(
        nlay > 0, 
        caObj.calipso.all_arrays['cloud_fraction']>0.5)
    isCalipsoCloudy = np.logical_and(
        isCalipsoCloudy, 
        caObj.calipso.all_arrays['total_optical_depth_5km']>0.15)
    #isCalipsoClear = caObj.calipso.all_arrays['cloud_fraction']<0.5
    isCalipsoClear = nlay == 0
    isCalipsoClear = np.logical_and(isCalipsoClear, meancl<0.01)
    isCalipsoClear = np.logical_and(
        isCalipsoClear, 
        caObj.calipso.all_arrays['total_optical_depth_5km']<0)
    isCalipsoSnowIce = np.logical_and(
        isCalipsoClear,
        caObj.calipso.all_arrays['nsidc_surface_type']>50)
    print np.sum(isCalipsoSnowIce)
    isCalipsoNotSnowIce = np.logical_and(
        isCalipsoClear,
        caObj.calipso.all_arrays['nsidc_surface_type']<=0)

    isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                                 caObj.avhrr.all_arrays['cloudtype']<21) 
    isClearPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                 caObj.avhrr.all_arrays['cloudtype']<5)
    gotLight = caObj.avhrr.all_arrays['sunz']<95
    nodata = np.sum(caObj.avhrr.all_arrays['cloudtype'][use]>200)
    use = np.logical_and(use, np.logical_or(isCloudyPPS, isClearPPS))
    use = np.logical_and(use, np.logical_or(isCalipsoCloudy, isCalipsoClear))
    PODcloudy = (np.sum(np.logical_and(isCalipsoCloudy[use], 
                                       isCloudyPPS[use]))*1.0
                 /np.sum(isCalipsoCloudy[use]))
    FARcloudy = (np.sum(np.logical_and(isCalipsoClear[use], 
                                       isCloudyPPS[use]))*1.0
                 /np.sum(isCloudyPPS[use]))
    PODclear = (np.sum(np.logical_and(isCalipsoClear[use], 
                                       isClearPPS[use]))*1.0
                 /np.sum(isCalipsoClear[use]))
    Num = np.sum(use)
    part_nodata = nodata*1.0/(nodata+Num)
    print "N: %d POD-cloudy: %3.2f FAR-cloudy: %3.2f POD-clear %3.2f"%(
        Num, PODcloudy, FARcloudy , PODclear)
    all_pix = use.copy() 
    for var in ['cma_testlist0',
                'cma_testlist1',
                'cma_testlist2',
                'cma_testlist3',
                'cma_testlist4',
                'cma_testlist5',
            ]:
        print var
        for bit_nr in xrange(0,16):
            if bit_nr not in name_dict[var].keys():
                print "not using", var, bit_nr
                continue   
            test_is_on = get_pixels_where_test_is_passed(
                caObj.avhrr.all_arrays[var], bit_nr=bit_nr)
            all_pix = np.logical_and(all_pix,np.equal(test_is_on,False))
            use_this_test = np.logical_and(use, test_is_on)
            #caObj.avhrr.all_arrays[var].
            #print np.sum(test_is_on), np.sum(all_pix), np.sum(use)
            #print np.sum(caObj.avhrr.all_arrays[var]==4)
            #print np.sum(np.logical_and(~test_is_on, 
            #                            caObj.avhrr.all_arrays[var]==4))
            PODcloudy = (np.sum(np.logical_and(isCalipsoCloudy[use_this_test], 
                                               isCloudyPPS[use_this_test]))*1.0
                         /np.sum(isCalipsoCloudy[use]))
            FARcloudy = (np.sum(np.logical_and(isCalipsoClear[use_this_test], 
                                           isCloudyPPS[use_this_test]))*1.0
                         /np.sum(isCloudyPPS[use_this_test]))
            PODclear = (np.sum(np.logical_and(isCalipsoClear[use_this_test], 
                                              isClearPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoClear[use]))
            PODsnow = (np.sum(np.logical_and(isCalipsoSnowIce[use_this_test], 
                                             isClearPPS[use_this_test]))*1.0
                       /np.sum(np.logical_and(isCalipsoSnowIce[use],gotLight[use])))
            FARsnow_not_clouds = (np.sum(np.logical_and(isCalipsoNotSnowIce[use_this_test], 
                                                        isClearPPS[use_this_test]))*1.0
                                  /np.sum(isClearPPS[use_this_test]))
            FARclear = (np.sum(np.logical_and(isCalipsoCloudy[use_this_test], 
                                              isClearPPS[use_this_test]))*1.0
                        /np.sum(isClearPPS[use_this_test]))
            MISSclear = (np.sum(np.logical_and(isCalipsoClear[use_this_test], 
                                              isCloudyPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoClear[use]))
            MISSsnow = (np.sum(np.logical_and(
                np.logical_and(
                    isCalipsoSnowIce[use_this_test],
                    gotLight[use_this_test]), 
                isCloudyPPS[use_this_test]))*1.0
                        /np.sum(isCalipsoClear[use]))
            LOSTsnow = (np.sum(np.logical_and(
                np.logical_and(
                    isCalipsoSnowIce[use_this_test],
                    gotLight[use_this_test]), 
                isCloudyPPS[use_this_test]))*1.0
                        /np.sum(np.logical_and(isCalipsoSnowIce[use],gotLight[use])))

            PRINT_FOR_OUTPUT_ON_SCREEN = True
            test_name = name_dict[var][bit_nr]
            if PRINT_FOR_OUTPUT_ON_SCREEN:
                print "%s test_bit: %s"%(var, str(bit_nr).rjust(2,' ') ),

            if (var ==  'cma_testlist5' or (var == 'cma_testlist4' and bit_nr>1) and 
                'snow' in name_dict[var][bit_nr]):
                print "N: %s POD-clear: %s FAR-clear: %s POD-snow: %s FAR_not_clouds %s LOSTsnow %s %s "%(
                    str(np.sum(use_this_test)).rjust(8,' '), 
                    ("%3.2f"%(PODclear*100)).rjust(5,' '), 
                    ("%3.2f"%(FARclear*100)).rjust(5,' '), 
                    ("%3.2f"%(PODsnow*100)).rjust(5,' '), 
                    ("%3.2f"%(FARsnow_not_clouds*100)).rjust(5,' '),
                    ("%3.2f"%(LOSTsnow*100)).rjust(5,' '),
                    test_name) 
            elif var ==  'cma_testlist5' or (var == 'cma_testlist4' and bit_nr>1):
                print "N: %s POD-clear: %s FAR-clear: %s %s"%(
                    str(np.sum(use_this_test)).rjust(8,' '), 
                    ("%3.2f"%(PODclear*100)).rjust(5,' '), 
                    ("%3.2f"%(FARclear*100)).rjust(5,' '), 
                    test_name)    
            else:    
                print "N: %s POD-cloudy: %s FAR-cloudy: %s MISS-clear: %s MISS-snow: %s %s"%(
                    str(np.sum(use_this_test)).rjust(8,' '), 
                    ("%3.2f"%(PODcloudy*100)).rjust(5,' '), 
                    ("%3.2f"%(FARcloudy*100)).rjust(5,' '), 
                    ("%3.2f"%(MISSclear*100)).rjust(5,' '), 
                    ("%3.2f"%(MISSsnow*100)).rjust(5,' '), 
                    test_name)  

            #print "%s test_bit:%d N: %d POD-cloudy: %3.2f FAR-cloudy: %3.2f POD-clear %3.2f FAR-clear %3.2f"%(
            #    var, bit_nr, Num, PODcloudy, FARcloudy , PODclear, FARclear)

    print "should be zero", np.sum(np.logical_and(isCalipsoClear[all_pix], 
                                                  isCloudyPPS[all_pix])),
    print np.sum(np.logical_and(isCalipsoClear[use], 
                                isCloudyPPS[use]))            

    Num = np.sum(use)
    part_nodata = nodata*1.0/(nodata+Num)
    #print "N: %d POD-cloudy: %3.2f FAR-cloudy: %3.2f POD-clear %3.2f"%(
    #    Num, PODcloudy, FARcloudy , PODclear)
                


ROOT_DIR = ("/home/a001865/DATA_MISC/reshaped_files/"
            "global_modis_14th_created20161108/")
ROOT_DIR_GAC = ("/home/a001865/DATA_MISC/reshaped_files/"
            "ATRAIN_RESULTS_GAC_oldctth/Reshaped_Files/noaa18/")
ROOT_DIR_GAC = ("/home/a001865/DATA_MISC/reshaped_files/"
            "ATRAIN_RESULTS_GAC_hsatz/Reshaped_Files/noaa18/")

files = glob(ROOT_DIR + "Reshaped_Files/merged/modis*h5")
files = glob(ROOT_DIR_GAC + "5km/20??/*/*/*h5")


test_list_file = "/home/a001865/git/acpg_develop/acpg/pges/cloudmask/pps_pge01_cmasktests.h"
TEST_NAMEFILE = open(test_list_file, 'r')
name_dict = {'cma_testlist0': {},
             'cma_testlist1': {},
             'cma_testlist2': {},
             'cma_testlist3': {},
             'cma_testlist4': {},
             'cma_testlist5': {}}
for line in TEST_NAMEFILE:
    if 'define SM_ACMG_' in line:
        list_of_line = line.split(' ')
        (name, list_nr, bit) = (list_of_line[1].replace('SM_ACMG_',''), list_of_line[2].replace(',',''), list_of_line[3])
        var = 'cma_testlist' + list_nr
        name_dict[var][int(bit)] =  name

    
caObjPPS = CalipsoAvhrrTrackObject()
for filename in files:
    caObjPPS +=  readCaliopAvhrrMatchObj(filename) 
    print "Scene %s"%(os.path.basename(filename))
use = caObjPPS.avhrr.all_arrays['bt11micron']>-9
use = caObjPPS.avhrr.all_arrays['sunz']>95
print_common_stats(caObjPPS, use, name_dict)

