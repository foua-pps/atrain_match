import os
from glob import glob
import re
import numpy as np
from scipy import ndimage
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)
from cloudsat_calipso_imager_plot import (drawCalClsatGEOPROFImagerPlot,
                                         drawCalPPSHeightPlot_PrototypePPSHeight  ) 
ROOT_DIR = ("/home/a001865/DATA_MISC/reshaped_files/"
            "global_modis_14th_created20161108/")
ROOT_DIR_GAC_nn = ("/home/a001865/DATA_MISC/reshaped_files/"
                   "ATRAIN_RESULTS_GAC_nn21/Reshaped_Files/noaa18/")
ROOT_DIR_GAC_old = ("/home/a001865/DATA_MISC/reshaped_files/"
                    "ATRAIN_RESULTS_GAC_old/Reshaped_Files/noaa18/")
ROOT_DIR_GAC_old = ("/home/a001865/DATA_MISC/reshaped_files/"
                    "ATRAIN_RESULTS_GAC_20161222/Reshaped_Files/noaa18/")

ROOT_DIR_GAC_nn_new = ("/home/a001865/DATA_MISC/reshaped_files/"
                       "ATRAIN_RESULTS_GAC_nn20161125/Reshaped_Files/noaa18/")
ROOT_DIR_GAC_nn_imager = ("/home/a001865/DATA_MISC/reshaped_files/"
                         "ATRAIN_RESULTS_GAC_nnimager_20161202/Reshaped_Files/noaa18/")

re_name = re.compile("_RESULTS_GAC_(\w+)\/")
caObj = CalipsoImagerTrackObject()

#for ROOT_DIR in [ROOT_DIR_GAC_nn_imager, ROOT_DIR_GAC_nn,  ROOT_DIR_GAC_nn_new]:
for ROOT_DIR in [ROOT_DIR_GAC_nn_imager]:
    match = re_name.search(ROOT_DIR)
    name = "no_name"
    if match:
        name = match.group(1)
    files = glob(ROOT_DIR + "5km/2009/*/*/*h5")
    for filename in files:
        basename = os.path.basename(filename)
        print "Scene %s"%(basename)
        caObj =  readCaliopImagerMatchObj(filename) 
        caObj_OLD =  readCaliopImagerMatchObj(filename.replace(ROOT_DIR, ROOT_DIR_GAC_old)) 
        height_pps =  caObj.imager.all_arrays['ctth_height']
        height_pps_old =  caObj_OLD.imager.all_arrays['ctth_height']
        use = np.logical_and(height_pps>0, height_pps_old>0)
        drawCalPPSHeightPlot_PrototypePPSHeight(caObj.calipso,
                                                use,
                                                height_pps_old + caObj.calipso.all_arrays['elevation'],
                                                height_pps + caObj.calipso.all_arrays['elevation'],                          
                                                "/home/a001865/PICTURES_FROM_PYTHON/CTTH_LAPSE_RATE_INVESTIGATION/",
                                                "test_plot_file_part_nn_ctth_%s_%s"%(name,basename.split('.h5')[0]),
                                                file_type='png',
                                                xmin=3000,
                                                xmax=6000,
                                                instrument='imager')
