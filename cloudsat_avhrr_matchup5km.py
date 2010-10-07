
import pdb

import os, string
import sys
from pps_basic_configure import *
from pps_error_messages import *
import numpy.oldnumeric as Numeric

from calipso5km import *
#from cloudsat5km import *
#from cloudsat5km_cwc import *
from cloudsat_calipso_avhrr_match import AREA5KM, RESHAPE_DIR
#MAIN_DIR = "/data/proj/safworks/adam/calipso_data"
#MAIN_DIR = "/local_disk/calipso_data"
#SUB_DIR = "metop02_calipso_2007spring"
#SUB_DIR = "noaa18_calipso_2007Aug"

#CTYPE_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
#AVHRR_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

#MAIN_DIR = "/local_disk/cloudsat_data"
#CLOUDSAT_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

# -----------------------------------------------------
def getCloudsat5kmAvhrr5kmMatch(avhrrfile5km,cloudsatfile5km,ctypefile5km,ctthfile5km,surftfile5km,sunanglefile5km,cloudsat_type):
    import string,os
    import pps_io
    import epshdf
    import numpy.oldnumeric as Numeric

    #print "Entering cloudsat_avhrr_matchup.py"
    basename = os.path.basename(ctypefile5km).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_")
    #dirname = os.path.dirname(ctypefile)
#    dirname = os.path.dirname(ctypefile5km).split("/")[0:2]
#    dirname = string.join(dirname,"/")
    savepath = "%s/%s/5km/%s/%s/%s" %(RESHAPE_DIR, base_sat, base_year, base_month, AREA5KM)
    ca_match_file5km = "%s/5km_%s_cloudsat-%s_avhrr_match.h5"%(savepath,basename,cloudsat_type)
    if not os.path.exists(savepath):
        os.makedirs(savepath)


    print "Match-up file: ",ca_match_file5km

    if not os.path.exists(ca_match_file5km):
        # Read AVHRR lon,lat data
        write_log("INFO","Read AVHRR geolocation data")
        pdb.set_trace()
        avhrrGeoObj5km = pps_io.readAvhrrGeoData(avhrrfile5km)
        
        # Read AVHRR sunsatangels
        write_log("INFO","Read AVHRR Sun -and Satellites Angels data")
        avhrrAngObj5km = pps_io.readSunSatAngles(sunanglefile5km) #, withAbsoluteAzimuthAngles=True)
        
        # Read AVHRR data
        write_log("INFO","Read AVHRR data")
        avhrrObj5km = pps_io.readAvhrrData(avhrrfile5km)

        # Read PPS Cloud Type data
        write_log("INFO","Read PPS Cloud Type")
        ctype5km = epshdf.read_cloudtype(ctypefile5km,1,1,0)
        try:
            ctth5km = epshdf.read_cloudtop(ctthfile5km,1,1,1,0,1)
        except:
            ctth5km = None
        
        # Read remapped NWP Surface temperature data
        write_log("INFO","Read NWP surface temperature")
        nwpinst5km = epshdf.read_nwpdata(surftfile5km)
        surft5km = nwpinst5km.gain*nwpinst5km.data.astype('d')+nwpinst5km.intercept

        # --------------------------------------------------------------------
        write_log("INFO","Read CLOUDSAT GEOPROF and CWC-RWOD data")
        # Read CLOUDSAT Radar data:
        if cloudsat_type == 'GEOPROF':
            from cloudsat5km import *
            write_log("INFO","Read CLOUDSAT GEOPROF data")
            cloudsat5km = reshapeCloudsat5km(cloudsatfile5km,avhrrGeoObj5km)
            #cloudsat5km = get_cloudsat5km(cloudsatfile5km)
            retv5km,min_diff,max_diff = match_cloudsat5km_avhrr5km(ctypefile5km,cloudsat5km,avhrrGeoObj5km,avhrrObj5km,ctype5km,ctth5km,surft5km,avhrrAngObj5km)
            writeCloudsat5kmAvhrr5kmMatchObj(ca_match_file5km,retv5km,6)
        elif cloudsat_type == 'CWC-RVOD':
            from cloudsat5km_cwc import *
            write_log("INFO","Read CLOUDSAT CWC-RWOD data")
            cloudsat5kmCwc = reshapeCloudsat5kmCwc(cloudsatfile5km,avhrrGeoObj5km)
            #cloudsat5kmCwc = get_cloudsat5kmCwc(cloudsatfile5km)
            retv5km,min_diff,max_diff = match_cloudsat5kmCwc_avhrr5km(ctypefile5km,cloudsat5kmCwc,avhrrGeoObj5km,avhrrObj5km,ctype5km,ctth5km,surft5km,avhrrAngObj5km)
            writeCloudsat5kmCwcAvhrr5kmMatchObj(ca_match_file5km,retv5km,6)   

    else:
        if cloudsat_type == 'GEOPROF':
            from cloudsat5km import *
            retv5km = readCloudsat5kmAvhrr5kmMatchObj(ca_match_file5km)
            min_diff = -9.0    # We don't store this information - only extracted during first time of processing
            max_diff = -9.0    # We don't store this information - only extracted during first time of processing
        else:
            from cloudsat5km_cwc import *       
            retv5km = readCloudsat5kmCwcAvhrr5kmMatchObj(ca_match_file5km)
            min_diff = -9.0    # We don't store this information - only extracted during first time of processing
            max_diff = -9.0    # We don't store this information - only extracted during first time of processing
    ##########################################################################################################################################
    #pdb.set_trace() 
    return retv5km,min_diff,max_diff

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

