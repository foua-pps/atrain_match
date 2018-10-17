import os
import re
import netCDF4
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
from read_cloudproducts_and_nwp_pps import (read_ctth_nc, read_pps_angobj_nc)
from read_modis_products import read_modis_h5

def print_data_file_pressure_calipso(out_file_h, reshaped_files):
    out_text = ""
    step = 50#hPa
    pressure_v = np.array(xrange(0,1200,step))
    pressure_v[-1]=1400
    num = len(pressure_v)-1
    num_satz00_10 = np.zeros(num)
    num_satz10_20 = np.zeros(num)
    num_satz20_30 = np.zeros(num)
    num_satz30_40 = np.zeros(num)
    num_satz40_50 = np.zeros(num)
    num_satz50_60 = np.zeros(num)
    num_satz60_70 = np.zeros(num)
    num_satz70_80 = np.zeros(num)
    for val in pressure_v:
        out_text += "%d "%(val + 0.5*step)
    out_text += "\n"

    satz_step = 10
    satz_v = np.array(xrange(0,80,satz_step))
    from matchobject_io import (readCaliopImagerMatchObj,
                                CalipsoImagerTrackObject)

    for filename in reshaped_files:
        print filename
        out_text += "CALIPSO %s\n"%(filename)
        caObj=readCaliopImagerMatchObj(filename)
        ca_pressure = caObj.calipso.all_arrays['layer_top_pressure'][:,0] #hPa nodata -9999
        use_ok = ca_pressure>0
        out_text += "CALIPSO satz00_00 "
        for ind in xrange(0,len(pressure_v)-1):
            h_min = pressure_v[ind]
            h_max = pressure_v[ind+1]
            use_p = np.logical_and(ca_pressure>=h_min, ca_pressure<h_max)
            use_i = np.logical_and(use_ok, use_p)
            out_text += "%d "%( np.sum(use_i)) 
        out_text += "\n"  
            #print out_text              
    out_file_h.write(out_text)
    out_file_h.write("\n")  

def print_data_file_pressure(out_file_h, nn_files, ANGLE_DIR, ctth_label="CTTHold", glob_ctth_label="CTTHold", modisfile_template="/temp/local/local/"):
    out_text = ""
    step = 50#hPa
    pressure_v = np.array(xrange(0,1200,step))
    pressure_v[-1]=1400
    num = len(pressure_v)-1
    num_satz00_10 = np.zeros(num)
    num_satz10_20 = np.zeros(num)
    num_satz20_30 = np.zeros(num)
    num_satz30_40 = np.zeros(num)
    num_satz40_50 = np.zeros(num)
    num_satz50_60 = np.zeros(num)
    num_satz60_70 = np.zeros(num)
    num_satz70_80 = np.zeros(num)
    for val in pressure_v:
        out_text += "%d "%(val + 0.5*step)
    out_text += "\n"

    satz_step = 10
    satz_v = np.array(xrange(0,80,satz_step))

    for filename in nn_files:
        out_text += "%s %s\n"%(ctth_label, filename)
        try:
            if 'C6'  in ctth_label:
                from runutils import parse_scenesfile_v2014
                satellite, datetime_obj = parse_scenesfile_v2014(filename)
                #this is the modis collection 6 data we should read!
                modis_pattern = datetime_obj.strftime(modisfile_template)
                print    modis_pattern
                modis_filename = glob(modis_pattern)
                print modis_filename 
                ctth = read_modis_h5(modis_filename[0])
                ctth.pressure = 100*ctth.pressure #hPa=>Pa
            else:
                ctth = read_ctth_nc(filename.replace( glob_ctth_label, ctth_label))
            anlge_file = filename.replace(glob_ctth_label, "sunsatangles")
            anlge_file = os.path.join(ANGLE_DIR,os.path.basename(anlge_file))
            pps_nc_ang = netCDF4.Dataset(anlge_file, 'r', format='NETCDF4')
            angles = read_pps_angobj_nc(pps_nc_ang)
            pps_nc_ang.close()
        except RuntimeError:
            print "skipping file: \n %s"%(filename.replace( "CTTHold", ctth_label))
            continue
        
        satz = angles.satz.data

        use_ok = np.logical_and(ctth.pressure>=7000,ctth.height<65535) 

        for satz_i in satz_v:
            use_satz = np.logical_and(satz>=satz_i, satz<satz_i+satz_step)
            out_text += "%s satz%d_%d "%(ctth_label, satz_i,satz_i+satz_step)
            for ind in xrange(0,len(pressure_v)-1):
                h_min = 100*pressure_v[ind]
                h_max = 100*pressure_v[ind+1]
                use_p = np.logical_and(ctth.pressure>=h_min, ctth.pressure<h_max)
                use_i = np.logical_and(use_ok, use_p)
                out_text += "%d "%( np.sum(np.logical_and(use_i,use_satz))) 
            out_text += "\n"  
            #print out_text              
    out_file_h.write(out_text)
    out_file_h.write("\n")     

def investigate_nn_ctth_satz():
    month = "02"
    ANGLE_DIR = "/nobackup/smhid13/sm_ninha/pps/modis_netcdf_01st_nn_ctth/import/SUNZ_data/%s/"%(month)
    PLOT_DIR = "/nobackup/smhid13/sm_ninha/atrain_matching/atrain_matching_01st_nnctth/pictures/"
    MODISFILE_TEMPLATE = "/home/a001865/DATA_MISC/atrain_match_testcases/modis/modis_06_file/MYD06_L2.A%Y%j.%H%M.006.*.h5"


    if 2==1:
        for hour in xrange(0,24):
            ROOT_DIR =  "/nobackup/smhid13/sm_ninha/pps/modis_netcdf_01st_nn_ctth/export/2010%s01_%02d/"%(month,hour) 
            out_filename = "/nobackup/smhid13/sm_ninha/atrain_matching/atrain_matching_01st_nnctth/pictures/satz_statistics_modis_%s_h%02d.txt"%(month,hour)
            out_file_h = open(out_filename,'w')
            nn_files =  glob(ROOT_DIR + "/*CTTHold_*.nc")
            for ctth_label in  ["MODIS-C6"]:
                print_data_file_pressure(out_file_h, nn_files,  ANGLE_DIR, ctth_label=ctth_label, glob_ctth_label="CTTHold",  modisfile_template=MODISFILE_TEMPLATE)
            out_file_h.close()

    #ROOT_DIR = "/run/media/a001865/SAMSUNG/A001865/DATA/ctth_satz_compare/"
    #ANGLE_DIR = "/run/media/a001865/SAMSUNG/A001865/DATA/GAC-NOAA18-2009/GAC_CLARA_A2_noaa18_2009_nwp_and_calipso/TEMP_GAC_DATA/"
    ROOT_DIR = "/home/a001865/SAFNWC_PPS/export"
    ANGLE_DIR = "/home/a001865/SAFNWC_PPS/import/ANC_data/remapped"
    MODISFILE_TEMPLATE = "/home/a001865/DATA_MISC/modis/MYD06/MYD06_L2.A%Y%j.%H%M.006.*.h5"
    if False:
        nn_files =  glob(ROOT_DIR + "/*CTTHnnImagerNoRTTOV_*20100114*.nc")
        out_filename = "satz_statistics_and_modis.txt"
        out_file_h = open(out_filename,'w')
        for ctth_label in  ["MODIS-C6", "CTTHold", "CTTHnnImager", "CTTHnnImagerNoRTTOV"]:
            print_data_file_pressure(out_file_h, nn_files,  ANGLE_DIR, 
                                     ctth_label=ctth_label, glob_ctth_label="CTTHnnImagerNoRTTOV",  
                                     modisfile_template=MODISFILE_TEMPLATE)
        out_file_h.close()

    # calipso från reshaped files:    
    if True:
        for month in ['02','04','06','08','10','12']:
            ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_01st_created20180316/Reshaped_Files_merged_calipso_cbase/eos2/1km/2010/%s/"%(month)
            match_files =  glob(ROOT_DIR + "/*caliop*.h5")
            out_filename = "satz_statistics_%s_and_calipso4.txt"%(month)
            out_file_h = open(out_filename,'w')   
            print_data_file_pressure_calipso(out_file_h,match_files)

if __name__ == "__main__":
    investigate_nn_ctth_satz()
