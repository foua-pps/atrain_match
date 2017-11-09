
import numpy as np
from glob import glob
import os
from matchobject_io import (readCaliopAvhrrMatchObj,
                            writeCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)


ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20161108/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"

OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20161108/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_nnviirs_20161205/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_nnviirs_20161205/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_nnviirs_20161205/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_nnviirs_20161205/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_nnmersi2_20161206/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_nnmersi2_20161206/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_AEROSOL/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_AEROSOL/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_AEROSOL/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_AEROSOL/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_AEROSOL_SMOKE/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_AEROSOL_SMOKE/Reshaped_Files/merged/"
ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_AEROSOL_SMOKE/Reshaped_Files/eos2/1km/%s/%s/*/*.h5"
OUT_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_MAY_AEROSOL_SMOKE/Reshaped_Files/merged/"

BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_nnMERSI2"
BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_nnAVHRR_20170315"
BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_nnVIIRS_20170315"

BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20170330"
#BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/ATRAIN_RESULTS_MODIS_NOVEMBER_nnAVHRR_with_gac"
ROOT_DIR = BASE_DIR + "/Reshaped_Files/eos2/1km/2010/%s/*/*caliop*.h5"
OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged/eos2/1km/2010/%s/"
outfile_template = "1km_%s_eos2_2010%s14_0000_00000_caliop_modis_match.h5"

BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_01st_created20170519"
ROOT_DIR = BASE_DIR + "/Reshaped_Files/eos2/1km/2010/%s/*/*caliop*.h5"
OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_calipso_cbase/eos2/1km/2010/%s/"
outfile_template = "1km_%s_eos2_2010%s01_0000_00000_caliop_modis_match.h5"

BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20170519"
ROOT_DIR = BASE_DIR + "/Reshaped_Files/eos2/1km/2010/%s/*/*caliop*.h5"
OUT_DIR_TEMPLATE = BASE_DIR + "/Reshaped_Files_merged_calipso_cbase//eos2/1km/2010/%s/"
outfile_template = "1km_%s_eos2_2010%s14_0000_00000_caliop_modis_match.h5"


caObj_night = CalipsoAvhrrTrackObject()
caObj_day = CalipsoAvhrrTrackObject()

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
            caObj_new=readCaliopAvhrrMatchObj(filename) 
            if np.nanmax(caObj_new.avhrr.all_arrays['sunz']>=90):
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
            writeCaliopAvhrrMatchObj(
                  outfile, caObj_night, avhrr_obj_name = 'pps')   
            caObj_night = CalipsoAvhrrTrackObject()    
        if num_d>0:
            filename_day = outfile_template%("day",month)
            outfile = os.path.join(OUT_DIR, filename_day)
            writeCaliopAvhrrMatchObj(
                outfile, caObj_day, avhrr_obj_name = 'pps')   
            caObj_day = CalipsoAvhrrTrackObject()   
