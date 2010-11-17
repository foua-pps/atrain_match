
import pdb

import os, string
import sys
from pps_basic_configure import *
from pps_error_messages import *
import numpy.oldnumeric as Numeric

from calipso import *

from cloudsat import *
from cloudsat_cwc import *

from config import AREA, RESHAPE_DIR
#MAIN_DIR = "/data/proj/safworks/adam/calipso_data"
#MAIN_DIR = "/local_disk/calipso_data"
#SUB_DIR = "metop02_calipso_2007spring"
#SUB_DIR = "noaa18_calipso_2007Aug"

#CTYPE_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
#AVHRR_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

#MAIN_DIR = "/local_disk/cloudsat_data"
#CLOUDSAT_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

# -----------------------------------------------------
def getCloudsatAvhrrMatch(avhrrfile,cloudsatfile,ctypefile,ctthfile,surftfile,sunanglefile,cloudsat_type):
    import string,os
    import pps_io
    import epshdf
    import numpy.oldnumeric as Numeric

    #print "Entering cloudsat_avhrr_matchup.py"
    basename = os.path.basename(ctypefile).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_")
    #dirname = os.path.dirname(ctypefile)
    #dirname = os.path.dirname(ctypefile).split("/")[0:2]
    #dirname = string.join(dirname,"/")
    
    savepath = "%s/%s/%skm/%s/%s/%s" %(RESHAPE_DIR, base_sat, RESOLUTION, base_year, base_month, AREA)
    ca_match_file = "%s/%skm_%s_cloudsat-%s_avhrr_match.h5"%(savepath, RESOLUTION, basename, cloudsat_type)
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    print "Match-up file: ",ca_match_file
    #ca_match_file = "%s/matched_files/1km_%s_cloudsat-%s_avhrr_match.h5"%(dirname,basename,cloudsat_type)
    #print "Match-up file: ",ca_match_file
    
    if not os.path.exists(ca_match_file):
        # Read AVHRR lon,lat data
        write_log("INFO","Read AVHRR geolocation data")
        avhrrGeoObj = pps_io.readAvhrrGeoData(avhrrfile)
        
        # Read AVHRR sunsatangels
        write_log("INFO","Read AVHRR Sun -and Satellites Angels data")
        avhrrAngObj = pps_io.readSunSatAngles(sunanglefile) #, withAbsoluteAzimuthAngles=True)
        
        # Read AVHRR data
        write_log("INFO","Read AVHRR data")
        avhrrObj = pps_io.readAvhrrData(avhrrfile)

        # Read PPS Cloud Type data
        write_log("INFO","Read PPS Cloud Type")
        ctype = epshdf.read_cloudtype(ctypefile,1,1,0)
        try:
            ctth = epshdf.read_cloudtop(ctthfile,1,1,1,0,1)
        except:
            ctth = None
        
        # Read remapped NWP Surface temperature data
        write_log("INFO","Read NWP surface temperature")
        nwpinst = epshdf.read_nwpdata(surftfile)
        surft = nwpinst.gain*nwpinst.data.astype('d')+nwpinst.intercept
        
        # --------------------------------------------------------------------
        #write_log("INFO","Read CLOUDSAT GEOPROF and CWC-RWOD data")
        # Read CLOUDSAT Radar data:
        if cloudsat_type == 'GEOPROF':
            write_log("INFO","Read CLOUDSAT GEOPROF data")
            cloudsat = reshapeCloudsat(cloudsatfile,avhrrGeoObj)
            #cloudsat = get_cloudsat(cloudsatfile)
            retv,min_diff,max_diff = match_cloudsat_avhrr(ctypefile,cloudsat,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj)
            writeCloudsatAvhrrMatchObj(ca_match_file,retv,6)
        elif cloudsat_type == 'CWC-RVOD':
            write_log("INFO","Read CLOUDSAT CWC-RWOD data")
            cloudsatcwc = reshapeCloudsat1kmCwc(cloudsatfile,avhrrGeoObj)
            #cloudsatcwc = get_cloudsatCwc(cloudsatfile)
            retv,min_diff,max_diff = match_cloudsatCwc_avhrr(ctypefile,cloudsatcwc,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj)
            writeCloudsatCwcAvhrrMatchObj(ca_match_file,retv,6)   

    else:
        if cloudsat_type == 'GEOPROF':
            retv = readCloudsatAvhrrMatchObj(ca_match_file)
            min_diff = -9.0    # We don't store this information - only extracted during first time of processing
            max_diff = -9.0    # We don't store this information - only extracted during first time of processing
        else:            
            retv = readCloudsatCwcAvhrrMatchObj(ca_match_file)
            min_diff = -9.0    # We don't store this information - only extracted during first time of processing
            max_diff = -9.0    # We don't store this information - only extracted during first time of processing

    return retv,min_diff,max_diff

# -----------------------------------------------------
#if __name__ == "__main__":
#    import sys
#    if len(sys.argv) < 5:
#	write_log("INFO","Usage: %s <cloudsat-hdf5-file> <pps cloudtype file> <pps ctth file> <pps avhrr file>"%sys.argv[0],moduleid=MODULE_ID)
#	sys.exit(-9)
#    else:
#        cloudsatfile = "%s/%s"%(CLOUDSAT_DIR,sys.argv[1])
#        ctypefile = "%s/%s"%(CTYPE_DIR,sys.argv[2])
#        ctthfile = "%s/%s"%(CTYPE_DIR,sys.argv[3])
#        avhrrfile = "%s/%s"%(AVHRR_DIR,sys.argv[4])

#    import string,os
#    import numpy.oldnumeric as Numeric
#    #import rpy
    
#    basename = os.path.basename(ctypefile).split(".h5")[0]
#    basename = string.join(basename.split("_")[0:4],"_")
    
#    sl = string.split(sys.argv[2],"_")
#    platform = sl[0]
#    norbit = string.atoi(sl[3])
#    yyyymmdd = sl[1]

#    caObj = getCloudsatAvhrrMatch(avhrrfile,cloudsatfile,ctypefile,ctthfile)

