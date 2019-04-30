import os
from glob import glob
import re
import numpy as np
from scipy import ndimage
from matchobject_io import readTruthImagerMatchObj
from plotting.along_track_plotting import (drawCalClsatGEOPROFImagerPlot,
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
#for ROOT_DIR in [ROOT_DIR_GAC_nn_imager, ROOT_DIR_GAC_nn,  ROOT_DIR_GAC_nn_new]:
for filename in files:
    match = re_name.search(ROOT_DIR)
    name = "no_name"
    if match:
        name = match.group(1)
    basename = os.path.basename(filename)
    caObj =  readTruthImagerMatchObj(filename) 
    caObj_OLD =  readTruthImagerMatchObj(filename.replace(ROOT_DIR_v2018_GAC, ROOT_DIR_v2014_GAC)) 
    height_pps =  caObj.imager.all_arrays['ctth_height']
    height_pps_old =  caObj_OLD.imager.all_arrays['ctth_height']
    try:
        use = np.logical_or(height_pps>0, height_pps_old>0)
    except:
        continue
    #use = np.logical_and(np.abs(caObj.imager.all_arrays['latitude'])>70, use)
    xmin = [i for i,x in enumerate(caObj.imager.all_arrays['latitude']) if x<-50]
    try:
        xmin = xmin[0]
    except:
        xmin = 0
    drawCalPPSHeightPlot_PrototypePPSHeight(caObj.calipso,
                                            use,
                                            height_pps_old + caObj.calipso.all_arrays['elevation'],
                                            height_pps + caObj.calipso.all_arrays['elevation'],                          
                                            ADIR + "/PICTURES_FROM_PYTHON/CTTH_LAPSE_RATE_INVESTIGATION/",
                                            "test_plot_file_part_nn_ctth_%s_%s"%(name,basename.split('.h5')[0]),
                                            file_type='png',
                                            xmin=xmin,
                                            xmax=xmin+2500,
                                            instrument='imager')#, MAXHEIGHT = 18000)
    
    #plt.show()
    #plt.close("all")
