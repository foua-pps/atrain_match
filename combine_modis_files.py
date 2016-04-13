
import numpy as np
from glob import glob
import os
from matchobject_io import (readCaliopAvhrrMatchObj,
                            writeCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
caObj = CalipsoAvhrrTrackObject()


ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_created20160402/Reshaped_Files/eos2/1km/2010/%s/*/*.h5"

#for month in ["01","02","03","04"]:
for month in ["01","02","03","04","05","06","07","08","09","10","11","12"]:
    files = glob(ROOT_DIR%(month))
    for filename in files:
        print  os.path.basename(filename)
        caObj_new=readCaliopAvhrrMatchObj(filename) 
        #print np.max(caObj_new.avhrr.all_arrays['surftemp'])
        if np.max(caObj_new.avhrr.all_arrays['sunz']>90):
            caObj = caObj + caObj_new
        else:
            print "skipping it"
    writeCaliopAvhrrMatchObj("modis_merged_night_%s.h5"%(month), 
                             caObj, avhrr_obj_name = 'pps')   
    caObj = CalipsoAvhrrTrackObject()    

