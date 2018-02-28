
import numpy as np
from glob import glob
import os
from matchobject_io import (readAmsrAvhrrMatchObj,
                            writeAmsrAvhrrMatchObj,
                            AmsrAvhrrTrackObject)

instrument = "seviri"
day = "14"
day_long = "14th"
satellite = "meteosat9"

BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_{instrument}_{day_long}_created20180222".format(
    day_long = day_long,
    instrument=instrument)
ROOT_DIR = BASE_DIR + "/Reshaped_Files/{satellite}/1km/2010/%s/*/*amsr*.h5".format(
    satellite=satellite)
OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_amsr/{satellite}/1km/2010/%s/".format(
    satellite=satellite)
outfile_template = "1km_%s_eos2_2010%s{day}_0000_00000_amsr_{instrument}_match.h5".format(
    day=day, instrument = instrument)

#BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20170519"
#ROOT_DIR = BASE_DIR + "/Reshaped_Files/eos2/1km/2010/%s/*/*amsr*.h5"
#OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_amsr/eos2/1km/2010/%s/"
#outfile_template = "1km_%s_eos2_2010%s14_0000_00000_amsr-GEOPROF_modis_match.h5"

clsatObj_night = AmsrAvhrrTrackObject()
clsatObj_day = AmsrAvhrrTrackObject()

for year in ["2010"]:#2012/02","2012/05", "2012/08", "2013/07", "2014/02", "2014/04", "2014/09"]:
    #for month in ["01","02","03","04","05","06","07","08","09","10","11","12"]:
    #for month in ["06"]:
    for month in ["01","02","03","04","05","06","07","08","09","10","11","12"]:
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
            clsatObj_new=readAmsrAvhrrMatchObj(filename) 
            if np.nanmax(clsatObj_new.avhrr.all_arrays['sunz']>=90):
                num_n +=1
                print "reading",os.path.basename(filename)
                clsatObj_night = clsatObj_night + clsatObj_new
            else :
                num_d +=1
                print "reading",os.path.basename(filename)
                clsatObj_day = clsatObj_day + clsatObj_new
        if num_n>0:
            filename_night = outfile_template%("night",month)
            outfile = os.path.join(OUT_DIR, filename_night)
            writeAmsrAvhrrMatchObj(
                  outfile, clsatObj_night)   
            clsatObj_night = AmsrAvhrrTrackObject()    
        if num_d>0:
            filename_day = outfile_template%("day",month)
            outfile = os.path.join(OUT_DIR, filename_day)
            writeAmsrAvhrrMatchObj(
                outfile, clsatObj_day)   
            clsatObj_day = AmsrAvhrrTrackObject()   
