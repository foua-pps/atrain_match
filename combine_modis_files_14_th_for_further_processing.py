
import numpy as np
from glob import glob
import os
from matchobject_io import (readCaliopImagerMatchObj,
                            writeCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)


instrument = "seviri"
day = "14"
day_long = "14th"
satellite = "meteosat9"

BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_{instrument}_{day_long}_created20180222".format(
    day_long = day_long,
    instrument=instrument)
ROOT_DIR = BASE_DIR + "/Reshaped_Files/{satellite}/1km/2010/%s/*/*caliop*.h5".format(
    satellite=satellite)
OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_calipso_cbase/{satellite}/1km/2010/%s/".format(
    satellite=satellite)
outfile_template = "1km_%s_{satellite}_2010%s{day}_0000_00000_caliop_{instrument}_match.h5".format(
    satellite=satellite, day=day, instrument = instrument)

print ROOT_DIR


caObj_night = CalipsoImagerTrackObject()
caObj_day = CalipsoImagerTrackObject()

for year in ["2010"]:#2012/02","2012/05", "2012/08", "2013/07", "2014/02", "2014/04", "2014/09"]:
    #for month in ["01","02","03","04","05","06","07","08","09","10","11","12"]:
    #for month in ["06"]:
    for month in ["02","03","04","05","06","07","08","09","10","11","12","01"]:
        OUT_DIR = OUT_DIR_TEMPLATE%(month)
        if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

        files = glob(ROOT_DIR%(month))
        if len(files)==0:
            continue
        num_n = 0  
        num_d = 0
        for filename in files:
            print  os.path.basename(filename)
            caObj_new=readCaliopImagerMatchObj(filename) 
            if np.nanmax(caObj_new.imager.all_arrays['sunz']>=90):
                num_n +=1
                print "reading",os.path.basename(filename)
                caObj_night = caObj_night + caObj_new
            else :
                num_d +=1
                print "reading",os.path.basename(filename)
                caObj_day = caObj_day + caObj_new
        if num_n>0:
            filename_night = outfile_template%("night",month)
            outfile = os.path.join(OUT_DIR, filename_night)
            writeCaliopImagerMatchObj(
                  outfile, caObj_night, imager_obj_name = 'pps')   
            caObj_night = CalipsoImagerTrackObject()    
        if num_d>0:
            filename_day = outfile_template%("day",month)
            outfile = os.path.join(OUT_DIR, filename_day)
            writeCaliopImagerMatchObj(
                outfile, caObj_day, imager_obj_name = 'pps')   
            caObj_day = CalipsoImagerTrackObject()   
