# Program orrb_process_master

# This program supervises the processing of a full month of CloudSat-CALIPSO-PPS matchups and
# generates output statistics and plots for each individual case

import os, string
import sys

import Numeric

STUDIED_MONTH = "200706"
MAIN_DATADIR = "/data/proj_nsc1/safworks/kgkarl/ORR-B-datasets/"
MAIN_RUNDIR = "/data/proj/saf/kgkarl/calipso_cloudsat/scr/orr-b-stat/"
PPS_DIR = "%sppsfiles/%s/" % (MAIN_DATADIR,STUDIED_MONTH)
CALIPSO_DIR = "%scalipso/%s/CLay/" % (MAIN_DATADIR,STUDIED_MONTH)
CLOUDSAT_DIR = "%scloudsat/%s/2B-GEOPROF/" % (MAIN_DATADIR,STUDIED_MONTH)
OUTPUT_DIR = "%sresults/" % MAIN_RUNDIR
INDATA_DIR = "%sindata/" % MAIN_RUNDIR


PLOT_OPTION = 1  # Set to zero if you don't want to generate the plot but just compute statistics

EMISS_FILTERING = 0  # A way to filter out cases with the thinnest topmost CALIPSO layers
ICE_COVER_SEA = 0 # If set, calculations are restricted to ice cover over sea using NSIDC and IGBP data
ICE_FREE_SEA = 0 # If set, calculations are restricted to ice-free sea using NSIDC and IGBP data
SNOW_COVER_LAND = 0 # If set, calculations are restricted to snow over land using NSIDC and IGBP data
SNOW_FREE_LAND = 0 # If set, calculations are restricted to snow-free land using NSIDC and IGBP data
COASTAL_ZONE = 0 # If set, calculations are restricted to coastal regions using NSIDC data (mixed microwave region)


# -----------------------------------------------------
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
	print "INFO","Usage: %s <YYYYMM> "%sys.argv[0]
	sys.exit(-9)
    else:
        STUDIED_MONTH = "%s" % sys.argv[1]
        PPS_DIR = "%sppsfiles/%s/" % (MAIN_DATADIR,STUDIED_MONTH)
        CALIPSO_DIR = "%scalipso/%s/CLay/" % (MAIN_DATADIR,STUDIED_MONTH)
        CLOUDSAT_DIR = "%scloudsat/%s/2B-GEOPROF/" % (MAIN_DATADIR,STUDIED_MONTH)

    # Read all input files from PPS, CloudSat and calipso and place them in lists

    pps_cloudtype_files = open("%scloudtype_list_%s.dat" % (INDATA_DIR, STUDIED_MONTH), "r")
    calipso_files = open("%scalipso_list_%s.dat" % (INDATA_DIR, STUDIED_MONTH), "r")
    cloudsat_files = open("%scloudsat_list_%s.dat" % (INDATA_DIR, STUDIED_MONTH), "r")
 
    cloudtype_list = pps_cloudtype_files.readlines()
    calipso_list = calipso_files.readlines()
    cloudsat_list = cloudsat_files.readlines()
    pps_cloudtype_files.close()
    calipso_files.close()
    cloudsat_files.close()
    
    # Process all cases one by one. Statistics will be produced and written (appended) by program
    # cloudsat_calipso_avhrr_orrb-stat.py

    processing_modes = ["BASIC","EMISSFILT","ICE_COVER_SEA","ICE_FREE_SEA","SNOW_COVER_LAND","SNOW_FREE_LAND","COASTAL_ZONE"]

    for scene in range(len(cloudtype_list)):
        cloudtype_file = "%s%s" % (PPS_DIR,cloudtype_list[scene])
        cloudtype_file = string.split(cloudtype_file,"\n")[0] # Just to remove training carriage return
        namecomps = string.split(cloudtype_list[scene],"cloudtype")
        ctth_file = "%s%s%s" % (PPS_DIR,namecomps[0],"ctth.h5")
        avhrr_file = "%s%s%s" % (PPS_DIR,namecomps[0],"avhrr.h5")
        nwp_tsur_file = "%s%s%s" % (PPS_DIR,namecomps[0],"nwp_tsur.h5")
        calipso_file = "%s%s" % (CALIPSO_DIR,calipso_list[scene])
        calipso_file = string.split(calipso_file,"\n")[0] # Just to remove training carriage return
        cloudsat_file = "%s%s" % (CLOUDSAT_DIR,cloudsat_list[scene])
        cloudsat_file = string.split(cloudsat_file,"\n")[0] # Just to remove training carriage return

        # Run repeatedly with 7 different processing modes:
        #            BASIC: Default run with plotting (actually the only case with plotting!)
        #            EMISSFILT: Run for the case with thinnest clouds removed
        #            ICE_COVER_SEA: Run for the case with only ice covered ocean
        #            ICE_FREE_SEA: Run for the case with only ice free ocean
        #            SNOW_COVER_LAND: Run for the case with only snow covered land
        #            SNOW_FREE_LAND: Run for the case with only snow free land
        #            COASTAL_ZONE: Run for the case with coastal zone

        for mode in processing_modes:        
            cmdstr ="python cloudsat_calipso_avhrr_orrb-stat.py %s %s %s %s %s %s %s" \
                 % (cloudsat_file,calipso_file,cloudtype_file,ctth_file,avhrr_file,nwp_tsur_file, mode)
            os.system(cmdstr)
