#HISTORY CHANGES BY ANKE TETZLAFF#

#080430:
# Could not run in the Arctic using the tmpaid creation
# Used hard coded area 'arctic_super_5010' instead

#080416: Got error message in line 232 and 234:
#         "only rank-0 arrays can be converted to Python scalars";
#         changed start_sec1970 and end_sec1970 to rank0 array 
#         by setting start_sec1970[0] and end_sec1970[0]
import inspect

from pps_basic_configure import *
from pps_error_messages import *

from calipso5km import *
from setup import AREA5KM, SUB_DIR, DATA_DIR

#MAIN_DIR = "/data/proj/safworks/adam/cloudsat_data"
#MAIN_DIR = "/data/proj_nsc1/safworks/kgkarl/ORR-B-datasets/cloudsat/"
#MAIN_DIR = "/data/proj_nsc1/safworks/calipso_cloudsat/data/arctic"
#SUB_DIR = "matchups/"
#SUB_DIR = "metop02_calipso_2007spring"
#SUB_DIR = "noaa18_calipso_cloudsat_2007DEC_KG"
#SUB_DIR = "noaa17_calipso_cloudsat_2007JUN_KG"
#SUB_DIR = "noaa18_calipso_cloudsat_2007JUN_KG"

#SATPOS_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
#EPHE_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
#COMPRESS_LVL = 6
COVERAGE_DIR = "%s/5km/%s"%(SUB_DIR,AREA5KM)
#NLINES=1000
#NLINES=6000
#SWATHWD=2048
#NODATA=-9


class Cloudsat5kmCwcRvodObject:
    def _init_(self):
	# Geolocation
        self.Profile_time=None
        self.TAI_start=None
        self.latitude=None
        self.longitude=None
        self.Height=None
        self.elevation=None
        self.sec_1970=None
        self.avhrr_linnum=None
        self.avhrr_pixnum=None
	# The data# Geolocation
        self.Data_quality=None
        self.Data_targetID=None
        self.RVOD_liq_water_path=None
        self.RVOD_liq_water_path_uncertainty=None
        self.RVOD_ice_water_path=None
        self.RVOD_ice_water_path_uncertainty=None
        self.LO_RVOD_liquid_water_path=None
        self.LO_RVOD_liquid_water_path_uncertainty=None
        self.IO_RVOD_ice_water_path=None
        self.IO_RVOD_ice_water_path_uncertainty=None	
        self.RVOD_liq_water_content=None
        self.RVOD_liq_water_content_uncertainty=None
        self.RVOD_ice_water_content=None
        self.RVOD_ice_water_content_uncertainty=None
        self.LO_RVOD_liquid_water_content=None
        self.LO_RVOD_liquid_water_content_uncertainty=None
        self.IO_RVOD_ice_water_content=None
        self.IO_RVOD_ice_water_content_uncertainty=None	
        self.Temp_min_mixph_K=None
        self.Temp_max_mixph_K=None

class Cloudsat5kmCwcAvhrr5kmTrackObject:
    def __init__(self):
        self.avhrr5km=ppsAvhrr5kmObject()
        self.cloudsat5kmcwc=Cloudsat5kmCwcRvodObject()
        self.diff_sec_1970=None
##########################################################################################################################################   
#pdb.set_trace()

# ----------------------------------------
def readCloudsat5kmCwcAvhrr5kmMatchObj(filename):
    import _pyhl

    retv5kmcwc = Cloudsat5kmCwcAvhrr5kmTrackObject()

    a=_pyhl.read_nodelist(filename)
    a.selectAll()
    a.fetch()

    # Match-Up - time difference:
    retv5kmcwc.diff_sec_1970 = a.getNode("/diff_sec_1970").data()

    # Cloudsat:

    # AVHRR:
    retv5kmcwc.avhrr5km.longitude = a.getNode("/avhrr/longitude").data()
    retv5kmcwc.avhrr5km.latitude = a.getNode("/avhrr/latitude").data()
    retv5kmcwc.avhrr5km.sec_1970 = a.getNode("/avhrr/sec_1970").data()
    retv5kmcwc.avhrr5km.cloudtype = a.getNode("/avhrr/cloudtype").data()
    retv5kmcwc.avhrr5km.ctth_height = a.getNode("/avhrr/ctth_height").data()
    retv5kmcwc.avhrr5km.ctth_pressure = a.getNode("/avhrr/ctth_pressure").data()
    retv5kmcwc.avhrr5km.ctth_temperature = a.getNode("/avhrr/ctth_temperature").data()
    retv5kmcwc.avhrr5km.bt11micron = a.getNode("/avhrr/bt11micron").data()
    retv5kmcwc.avhrr5km.bt12micron = a.getNode("/avhrr/bt12micron").data()
    retv5kmcwc.avhrr5km.surftemp = a.getNode("/avhrr/surftemp").data()
    retv5kmcwc.avhrr5km.satz = a.getNode("/avhrr/satz").data()
    # CWC-RVOD:
    # ====
# Geolocation
    retv5kmcwc.cloudsat5kmcwc.longitude = a.getNode("/cloudsatcwc/longitude").data()
    retv5kmcwc.cloudsat5kmcwc.latitude = a.getNode("/cloudsatcwc/latitude").data()
    retv5kmcwc.cloudsat5kmcwc.elevation = a.getNode("/cloudsatcwc/elevation").data()
    retv5kmcwc.cloudsat5kmcwc.Height = a.getNode("/cloudsatcwc/Height").data()
    retv5kmcwc.cloudsat5kmcwc.avhrr_linnum = a.getNode("/cloudsatcwc/avhrr_linnum").data()
    retv5kmcwc.cloudsat5kmcwc.avhrr_pixnum = a.getNode("/cloudsatcwc/avhrr_pixnum").data()
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    retv5kmcwc.cloudsat5kmcwc.Profile_time = a.getNode("/cloudsatcwc/Profile_time").data()
    #retv5kmcwc.cloudsat5kmcwc.TAI_start = a.getNode("/cloudsatcwc/TAI_start").data()
    retv5kmcwc.cloudsat5kmcwc.sec_1970 = a.getNode("/cloudsatcwc/sec_1970").data()
# The data
    retv5kmcwc.cloudsat5kmcwc.Data_quality = a.getNode("/cloudsatcwc/Data_quality").data()
    retv5kmcwc.cloudsat5kmcwc.Data_targetID = a.getNode("/cloudsatcwc/Data_targetID").data()
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_path = a.getNode("/cloudsatcwc/RVOD_liq_water_path").data()
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_path_uncertainty = a.getNode("/cloudsatcwc/RVOD_liq_water_path_uncertainty").data()
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_path = a.getNode("/cloudsatcwc/RVOD_ice_water_path").data()
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_path_uncertainty = a.getNode("/cloudsatcwc/RVOD_ice_water_path_uncertainty").data()
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_path = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_path").data()
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_path_uncertainty = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_path_uncertainty").data()
    retv5kmcwc.cloudsat5kmcwc.IO_RVOD_ice_water_path = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_path").data()
    retv5kmcwc.cloudsat5kmcwc.IO_RVOD_ice_water_path_uncertainty = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_path_uncertainty").data()

   
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_content = a.getNode("/cloudsatcwc/RVOD_liq_water_content").data()
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_content_uncertainty = a.getNode("/cloudsatcwc/RVOD_liq_water_content_uncertainty").data()
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_content = a.getNode("/cloudsatcwc/RVOD_ice_water_content").data()
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_content_uncertainty = a.getNode("/cloudsatcwc/RVOD_ice_water_content_uncertainty").data()
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_content = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_content").data()
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_content_uncertainty = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_content_uncertainty").data()
    retv5kmcwc.cloudsat5kmcwc.IO_RVOD_ice_water_content = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_content").data()
    retv5kmcwc.cloudsat5kmcwc.IO_RVOD_ice_water_content_uncertainty = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_content_uncertainty").data()
    retv5kmcwc.cloudsat5kmcwc.Temp_min_mixph_K = a.getNode("/cloudsatcwc/Temp_min_mixph_K").data()
    retv5kmcwc.cloudsat5kmcwc.Temp_max_mixph_K = a.getNode("/cloudsatcwc/Temp_max_mixph_K").data()

##########################################################################################################################################
    #pdb.set_trace()

    return retv5kmcwc


# ----------------------------------------
def writeCloudsat5kmCwcAvhrr5kmMatchObj(filename,ca_obj,compress_lvl):
    import _pyhl
    status = -1
    
    a=_pyhl.nodelist()

    shape = [ca_obj.cloudsat5kmcwc.longitude.shape[0]]

    # Match-Up - time difference:
    # ====
    b=_pyhl.node(_pyhl.DATASET_ID,"/diff_sec_1970")
    b.setArrayValue(1,shape,ca_obj.diff_sec_1970,"double",-1)
    a.addNode(b)


    # AVHRR
    # ====
    b=_pyhl.node(_pyhl.GROUP_ID,"/avhrr")
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/longitude")
    shape = ca_obj.avhrr5km.longitude.shape
    b.setArrayValue(1,shape,ca_obj.avhrr5km.longitude,"float",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/latitude")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.latitude,"float",-1)
    a.addNode(b)

    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/sec_1970")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.sec_1970,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_pressure")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.ctth_pressure,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_temperature")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.ctth_temperature,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_height")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.ctth_height,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/cloudtype")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.cloudtype,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/bt11micron")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.bt11micron,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/bt12micron")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.bt12micron,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/surftemp")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.surftemp,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/satz")
    b.setArrayValue(1,shape,ca_obj.avhrr5km.satz,"double",-1)
    a.addNode(b)
    # CWC-RVOD
    # ====
    shapecwc = [ca_obj.cloudsat5kmcwc.longitude.shape[0]]
    shape2dcwc = ca_obj.cloudsat5kmcwc.Height.shape
    shapeTAIcwc = [ca_obj.cloudsat5kmcwc.Temp_min_mixph_K.shape[0]]
    shapecwclong = ca_obj.cloudsat5kmcwc.Data_quality.shape
    shapecwclongtrippel = ca_obj.cloudsat5kmcwc.Profile_time.shape
# Geolocation
    b=_pyhl.node(_pyhl.GROUP_ID,"/cloudsatcwc")
    a.addNode(b)
    
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/longitude")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.longitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/latitude")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.latitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/elevation")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.elevation,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Height")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.Height,"double",-1)
    a.addNode(b)       
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/avhrr_linnum")
    b.setArrayValue(1,shape,ca_obj.cloudsat5kmcwc.avhrr_linnum,"int",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/avhrr_pixnum")
    b.setArrayValue(1,shape,ca_obj.cloudsat5kmcwc.avhrr_pixnum,"int",-1)
    a.addNode(b)
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Profile_time")
    b.setArrayValue(1,shapecwclongtrippel,ca_obj.cloudsat5kmcwc.Profile_time,"double",-1)
    a.addNode(b) 
    #b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/TAI_start")
    #b.setArrayValue(1,shapeTAIcwc,ca_obj.cloudsat5kmcwc.TAI_start,"double",-1)
    #a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/sec_1970")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.sec_1970,"double",-1)
    a.addNode(b)
# The data
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Data_quality")
    b.setArrayValue(1,shapecwclong,ca_obj.cloudsat5kmcwc.Data_quality,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Data_targetID")
    b.setArrayValue(1,shapecwclong,ca_obj.cloudsat5kmcwc.Data_targetID,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.RVOD_liq_water_path,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.RVOD_liq_water_path_uncertainty,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.RVOD_ice_water_path,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.RVOD_ice_water_path_uncertainty,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.LO_RVOD_liquid_water_path,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.LO_RVOD_liquid_water_path_uncertainty,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.IO_RVOD_ice_water_path,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsat5kmcwc.IO_RVOD_ice_water_path_uncertainty,"double",-1)
    a.addNode(b)
##########################################################################################################################################
    #pdb.set_trace()
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.RVOD_liq_water_content,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.RVOD_liq_water_content_uncertainty,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.RVOD_ice_water_content,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.RVOD_ice_water_content_uncertainty,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.LO_RVOD_liquid_water_content,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.LO_RVOD_liquid_water_content_uncertainty,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.IO_RVOD_ice_water_content,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsat5kmcwc.IO_RVOD_ice_water_content_uncertainty,"double",-1)
    a.addNode(b)    

    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Temp_min_mixph_K")
    b.setArrayValue(1,shapeTAIcwc,ca_obj.cloudsat5kmcwc.Temp_min_mixph_K,"double",-1)
    a.addNode(b)

    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Temp_max_mixph_K")
    b.setArrayValue(1,shapeTAIcwc,ca_obj.cloudsat5kmcwc.Temp_max_mixph_K,"double",-1)
    a.addNode(b)

    status = a.write(filename,compress_lvl)

    return status

# -----------------------------------------------------
def select_cloudsat5kmCwc_inside_avhrr(cloudsat5kmObj,cal5km,sec1970_start_end,sec_timeThr):
    # Cloudsat times are already in UTC in sec since 1970
    import numpy

    sec1970_start,sec1970_end = sec1970_start_end
    
    # Select the points inside the avhrr swath:
    # Allowing for sec_timeThr seconds deviation:
    idx_time_okay5km = numpy.logical_and(numpy.greater(\
        cloudsat5kmObj.sec1970,sec1970_start - sec_timeThr),
                                   numpy.less(\
        cloudsat5kmObj.sec1970,sec1970_end + sec_timeThr))
    
    #idx_match5km = numpy.not_equal(cal5km,NODATA)
    idx_match5km = numpy.where(numpy.not_equal(cal5km,NODATA),idx_time_okay5km,False)
    
    return idx_match5km

# -----------------------------------------------------
def get_cloudsat5kmCwc(filenamecwc):
    import _pypps_filters
    import time

    # Read CLOUDSAT Radar data:
    cloudsatcwc = read_cloudsat5kmCwc(filenamecwc)

    # CWC-RVOD
    # ====

    lonCloudsat5kmcwc = cloudsatcwc.longitude[:,1].ravel()
    latCloudsat5kmcwc = cloudsatcwc.latitude[:,1].ravel()
    ndimcwc = lonCloudsat5kmcwc.shape[0]
    # Time in seconds since 1970:
    
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone

    print "CWC-RVOD TAI start time: ",cloudsatcwc.TAI_start
    start_sec1970cwc = cloudsatcwc.TAI_start + dsec
    start_timecwc = time.gmtime(start_sec1970cwc[0])
    end_sec1970cwc = cloudsatcwc.TAI_start + dsec + cloudsatcwc.Profile_time[ndimcwc-1,1]
    end_timecwc = time.gmtime(end_sec1970cwc[0])
    print "CWC-RVOD Start and end times: ",start_timecwc,end_timecwc
    cloudsatcwc.sec1970 = cloudsatcwc.Profile_time[:,1] + start_sec1970cwc

    # Any other possible derived field should come here...

    # --------------------------------------------------------------------

    return cloudsatcwc

# -----------------------------------------------------
def read_cloudsat5kmCwc(filenamecwc):
    import _pyhl
    
    import numpy
    
    # CWC-RVOD
    # ====
    acwc=_pyhl.read_nodelist(filenamecwc)
    bcwc=acwc.getNodeNames()
    acwc.selectAll()
    acwc.fetch()

    retv5kmcwc = Cloudsat5kmCwcRvodObject()
    
    root5kmcwc="/cloudsat"
# Geolocation
    ccwc=acwc.getNode("%s/Geolocation Fields/Longitude"%root5kmcwc)
    retv5kmcwc.longitude=ccwc.data()
    ccwc=acwc.getNode("%s/Geolocation Fields/Latitude"%root5kmcwc)
    retv5kmcwc.latitude=ccwc.data()
    ccwc=acwc.getNode("%s/Geolocation Fields/DEM_elevation"%root5kmcwc)
    retv5kmcwc.elevation=ccwc.data()
    ccwc=acwc.getNode("%s/Geolocation Fields/Height"%root5kmcwc)
    retv5kmcwc.Height=ccwc.data()
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    ccwc=acwc.getNode("%s/Geolocation Fields/Profile_time"%root5kmcwc)
    retv5kmcwc.Profile_time=ccwc.data()
    ccwc=acwc.getNode("%s/Geolocation Fields/TAI_start"%root5kmcwc)
    retv5kmcwc.TAI_start=ccwc.data()
# The data
    ccwc=acwc.getNode("%s/Data Fields/Data_quality"%root5kmcwc)
    retv5kmcwc.Data_quality=ccwc.data()
    ccwc=acwc.getNode("%s/Data Fields/Data_targetID"%root5kmcwc)
    retv5kmcwc.Data_targetID=ccwc.data()
# IWP and LWP		
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_path"%root5kmcwc)
    retv5kmcwc.RVOD_liq_water_path=ccwc.data()
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_path_uncertainty"%root5kmcwc)
    retv5kmcwc.RVOD_liq_water_path_uncertainty=ccwc.data()    
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_path"%root5kmcwc)
    retv5kmcwc.RVOD_ice_water_path=ccwc.data()
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_path_uncertainty"%root5kmcwc)
    retv5kmcwc.RVOD_ice_water_path_uncertainty=ccwc.data()
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_path"%root5kmcwc)
    retv5kmcwc.LO_RVOD_liquid_water_path=ccwc.data()			
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_path_uncertainty"%root5kmcwc)
    retv5kmcwc.LO_RVOD_liquid_water_path_uncertainty=ccwc.data()
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_path"%root5kmcwc)
    retv5kmcwc.IO_RVOD_ice_water_path=ccwc.data()
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_path_uncertainty"%root5kmcwc)
    retv5kmcwc.IO_RVOD_ice_water_path_uncertainty=ccwc.data()
# IWC and LWC
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_content"%root5kmcwc)
    retv5kmcwc.RVOD_liq_water_content=ccwc.data()#.astype('Int16')		
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_content_uncertainty"%root5kmcwc)
    retv5kmcwc.RVOD_liq_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_content"%root5kmcwc)
    retv5kmcwc.RVOD_ice_water_content=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_content_uncertainty"%root5kmcwc)
    retv5kmcwc.RVOD_ice_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_content"%root5kmcwc)
    retv5kmcwc.LO_RVOD_liquid_water_content=ccwc.data()#.astype('Int16')			
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_content_uncertainty"%root5kmcwc)
    retv5kmcwc.LO_RVOD_liquid_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_content"%root5kmcwc)
    retv5kmcwc.IO_RVOD_ice_water_content=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_content_uncertainty"%root5kmcwc)
    retv5kmcwc.IO_RVOD_ice_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/Temp_min_mixph_K"%root5kmcwc)
    retv5kmcwc.Temp_min_mixph_K=ccwc.data()
    ccwc=acwc.getNode("%s/Data Fields/Temp_max_mixph_K"%root5kmcwc)
    retv5kmcwc.Temp_max_mixph_K=ccwc.data()

    return retv5kmcwc

#retv5kmcwc.RVOD_liq_water_content.shape


# --------------------------------------------
def get_cloudsat5kmCwc_avhrr5km_linpix(avhrrIn,lon,lat,clCwcTime):
    #import numpy.oldnumpy as numpy
    import numpy
    tmppcs="tmpproj"
    define_pcs5km(tmppcs, "Plate Caree, central meridian at 15E",
               ['proj=eqc','ellps=bessel', 'lon_0=15'])

    startline=0
    Inside=0
    HasEncounteredMatch=0
    i=0

    if not os.path.exists(COVERAGE_DIR):
        os.makedirs(COVERAGE_DIR)

    while startline < avhrrIn.longitude.shape[0]:
        
        write_log("INFO","Calling get_cloudsatxwx_avhrr_linpix: start-line = ",startline)
        tmpaid = "tmparea_%d"%i
        endline = startline+NLINES
        coverage_filename = "%s/coverage_avhrr_cloudsat_matchup_%s_%.5d_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrrIn.satellite,avhrrIn.orbit,
                             startline,endline,tmpaid)
        write_log("INFO","Coverage filename = ",coverage_filename)
        ##########################################################################################################################################
        #pdb.set_trace()
        cal,cap,ok = get_cloudsat5kmCwc_avhrr5km_linpix_segment(avhrrIn,lon,lat,clCwcTime,
                                                      (startline,endline),
                                                      SWATHWD5km,tmppcs,tmpaid,
                                                      coverage_filename)
        if ok:
            HasEncounteredMatch=1
            write_log("INFO","There was a match...")
            
        #if not ok and HasEncounteredMatch:
        #    write_log("INFO","Data is now empty. Leave the loop...")
        #    break
        
        if(startline==0):
            # First time:
            cloudsat5km_avhrr_line,cloudsat5km_avhrr_pixel = numpy.array(cal),numpy.array(cap)
        else:
            # Merge:
            cloudsat5km_avhrr_line = numpy.where(numpy.equal(cloudsat5km_avhrr_line,-9),cal,cloudsat5km_avhrr_line)
            cloudsat5mk_avhrr_pixel = numpy.where(numpy.equal(cloudsat5km_avhrr_pixel,-9),cap,cloudsat5km_avhrr_pixel)

        startline=startline+NLINES
        i=i+1
        ##########################################################################################################################################
        #pdb.set_trace()
    return cloudsat5km_avhrr_line,cloudsat5km_avhrr_pixel

# --------------------------------------------
def get_cloudsat5kmCwc_avhrr5km_linpix_segment(avhrrIn,lon,lat,timecwc,lines,swath_width,tmppcs,
                                      tmpaid,covfilename):
    #import numpy.oldnumpy as numpy
    import numpy
    import _satproj
    import area,pcs
    import pps_gisdata

    ndim = lon.shape[0]
    
    if avhrrIn.longitude.shape[0] > lines[1]:
        lines_end = lines[1]
    else:
        lines_end = avhrrIn.longitude.shape[0]
    lines_start = lines[0]
    
    write_log("INFO","lines_start,lines_end: ",lines_start,lines_end)
    nlines = lines_end - lines_start
    lonarr = avhrrIn.longitude[lines_start:lines_end,::]
    latarr = avhrrIn.latitude[lines_start:lines_end,::]

    idx_start = lines_start*swath_width
    idx_end   = lines_end*swath_width
    idx = numpy.arange(idx_start,idx_end)
    
    linearr = numpy.divide(idx,swath_width)
    write_log("INFO","Start and end line numbers: ",linearr[0],linearr[idx.shape[0]-1])
    
    linearr = numpy.reshape(linearr,(nlines,swath_width))
    pixelarr = numpy.fmod(idx,swath_width).astype('l')
    pixelarr = numpy.reshape(pixelarr,(nlines,swath_width))

    """
    write_log("INFO","Get bounding box...")
    bounds = getBoundingBox(lonarr,latarr)
    write_log("INFO","Bounding box (lon,lat): ",bounds)
    
    area_x_tup = pps_gisdata.c2s([(bounds[0],bounds[1]),(bounds[2],bounds[3])],tmppcs)
    dimx = int((area_x_tup[1][0] - area_x_tup[0][0])/1000.0 + 1.0)
    dimy = int((area_x_tup[1][1] - area_x_tup[0][1])/1000.0 + 1.0)
    write_log("INFO","X,Y dims: ",dimx,dimy)
    define_longlat_ll(tmpaid, "Temp area def", tmppcs,
                      pcs.d2r((bounds[0],bounds[1])), # lower left corner (lon, lat)
                      (dimx,dimy), 1000)

    areaObj = area.area(tmpaid)
    """
    areaObj = area.area(AREA5KM)

    
    #if not os.path.exists(covfilename):
    #    write_log("INFO","Create Coverage map...")
    #    cov = _satproj.create_coverage(areaObj,lonarr,latarr,1)
    #    print covfilename
    #    writeCoverage(cov,covfilename,"satproj","arctic_super_5010")
    #else:
    #    write_log("INFO","Read the AVHRR-CLOUDSAT matchup coverage from file...")
    #    cov,info = readCoverage(covfilename)
    write_log("INFO","Create Coverage map...")
    cov = _satproj.create_coverage(areaObj,lonarr,latarr,0)
    print covfilename
    writeCoverage5km(cov,covfilename,"satproj",AREA5KM)
    
    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA)
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA)


    write_log("INFO","Go through cloudsat track:")
    cloudsat_avhrr_line = []
    cloudsat_avhrr_pixel = []
    cloudsat_avhrr_line_time = []
    cloudsat_avhrr_pixel_time = []
    kk=0
    for i in range(ndim):
        xy_tup=pps_gisdata.lonlat2xy(AREA5KM,lon[i],lat[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
##        dimx=4500 Should be 5010!!!/KG
##        dimy=4500
        dimx=mapped_line.shape[1]
        dimy=mapped_line.shape[0]
        #if i== 1301:
        #    pdb.set_trace()
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            #if mapped_line[y,x]== 5512:
            #    pdb.set_trace()
            cloudsat_avhrr_line.append(mapped_line[y,x])
            cloudsat_avhrr_pixel.append(mapped_pixel[y,x])
            cloudsat_avhrr_line_time.append(-9)
            cloudsat_avhrr_pixel_time.append(-9)
            
        else:
            cloudsat_avhrr_line.append(-9)
            cloudsat_avhrr_pixel.append(-9)
            cloudsat_avhrr_line_time.append(-9)
            cloudsat_avhrr_pixel_time.append(-9)
            kk=kk+1

    cloudsat_avhrr_line = numpy.array(cloudsat_avhrr_line)
    cloudsat_avhrr_pixel = numpy.array(cloudsat_avhrr_pixel)
    cloudsat_avhrr_line_time = numpy.array(cloudsat_avhrr_line_time)
    cloudsat_avhrr_pixel_time = numpy.array(cloudsat_avhrr_pixel_time)
    # Control the time diference
    match_cloudsat_points = numpy.where(numpy.not_equal(cloudsat_avhrr_line,-9))
    avhrr_time = (cloudsat_avhrr_line[match_cloudsat_points] * DSEC_PER_AVHRR_SCALINE5KM) + avhrrIn.sec1970_start
    cl_time = timecwc[match_cloudsat_points]
    time_diff = avhrr_time-cl_time
    #max_time_diff_allowed = 50*60 #Based on that a lap is 102 min
    max_time_diff_allowed = sec_timeThr
    time_match = numpy.where(abs(time_diff)<max_time_diff_allowed)
    if time_match[0].shape[0]==0:
        x=numpy.repeat(cloudsat_avhrr_line_time,numpy.not_equal(cloudsat_avhrr_line_time,-9))
    else:
        cloudsat_avhrr_line_time[match_cloudsat_points[0][time_match]]= cloudsat_avhrr_line[match_cloudsat_points[0][time_match]]
        cloudsat_avhrr_pixel_time[match_cloudsat_points[0][time_match]] = cloudsat_avhrr_pixel[match_cloudsat_points[0][time_match]]
        x=numpy.repeat(cloudsat_avhrr_line_time,numpy.not_equal(cloudsat_avhrr_line_time,-9)) 

    write_log("INFO","Number of matching points = ",x.shape[0])
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0
    ##########################################################################################################################################
    #pdb.set_trace()
    return cloudsat_avhrr_line_time,cloudsat_avhrr_pixel_time,matchOk

# -----------------------------------------------------
def match_cloudsat5kmCwc_avhrr5km(ctypefile,cloudsat5kmObjcwc,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj):
    #import numpy.oldnumpy as numpy
    import time
    import string
    import numpy
    retv5kmcwc = Cloudsat5kmCwcAvhrr5kmTrackObject()

    lonCloudsat5kmcwc = cloudsat5kmObjcwc.longitude[:,1].ravel()		
    latCloudsat5kmcwc = cloudsat5kmObjcwc.latitude[:,1].ravel()
    timeCloudsat5kmcwc = cloudsat5kmObjcwc.sec1970.ravel()
    ProfileTimeCloudsat5kmcwc = cloudsat5kmObjcwc.Profile_time[:,1].ravel()	
    ndim5kmcwc = lonCloudsat5kmcwc.shape[0]				
    # Add time
    if avhrrGeoObj.sec1970_end < avhrrGeoObj.sec1970_start:
        avhrr5km_sec1970_end = int(DSEC_PER_AVHRR_SCALINE5KM*avhrrObj.num_of_lines+avhrrGeoObj.sec1970_start)
        avhrrGeoObj.sec1970_end = avhrr5km_sec1970_end
    # --------------------------------------------------------------------

    cal5kmcwc,cap5kmcwc=get_cloudsat5kmCwc_avhrr5km_linpix(avhrrGeoObj,lonCloudsat5kmcwc,latCloudsat5kmcwc,timeCloudsat5kmcwc)	# Cwc

    #print len(cal), len(cap), ndim

    idx_cloudsatcwc = numpy.arange(ndim5kmcwc)							# Cwc
    idx_cloudsatcwc_on_avhrr=numpy.repeat(idx_cloudsatcwc,numpy.not_equal(cal5kmcwc,NODATA)) 	# Cwc   
    lon_cloudsatcwc = numpy.repeat(lonCloudsat5kmcwc,numpy.not_equal(cal5kmcwc,NODATA))		# Cwc
    lat_cloudsatcwc = numpy.repeat(latCloudsat5kmcwc,numpy.not_equal(cal5kmcwc,NODATA))		# Cwc
    ProfileTime_cloudsatcwc = numpy.repeat(ProfileTimeCloudsat5kmcwc,numpy.not_equal(cal5kmcwc,NODATA))
    # Cloudsat line,pixel inside AVHRR swath:

    cal_on_avhrrcwc = numpy.repeat(cal5kmcwc,numpy.not_equal(cal5kmcwc,NODATA))	# Cwc
    cap_on_avhrrcwc = numpy.repeat(cap5kmcwc,numpy.not_equal(cal5kmcwc,NODATA))	# Cwc
    timetup_startcwc = time.gmtime(avhrrGeoObj.sec1970_start)		# Cwc
    timetup_endcwc   = time.gmtime(avhrrGeoObj.sec1970_end)		# Cwc

    #sec_timeThrcwc = 60*20 # Allow for 20 minute deviation between AVHRR and CLOUDSAT Cwcmatchup
    
    secTupcwc = avhrrGeoObj.sec1970_start,avhrrGeoObj.sec1970_end
    idx_matchcwc = select_cloudsat5kmCwc_inside_avhrr(cloudsat5kmObjcwc,cal5kmcwc,secTupcwc,sec_timeThr)	########kanske blir fel??????????

    print "Generate all datatypes (lat,lon,cty,ctth,surft) on the cloudsat track!"

    ctype_track = []
    ctth_height_track = []
    ctth_pressure_track = []
    ctth_temperature_track = []
    lon_avhrr_track = []
    lat_avhrr_track = []
    surft_track = []
    bt11micron_track = []
    bt12micron_track = []
    satz_track = []
    for i in range(cal_on_avhrrcwc.shape[0]):
        lat_avhrr_track.append(avhrrGeoObj.latitude[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]])
        lon_avhrr_track.append(avhrrGeoObj.longitude[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]])
        ctype_track.append(ctype.cloudtype[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]])
        surft_track.append(surft[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]])
        if avhrrObj.channels[3].data[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] == avhrrObj.nodata:
            b11 = -9.
        else:
            b11 = avhrrObj.channels[3].data[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] * avhrrObj.channels[3].gain + avhrrObj.channels[3].intercept
        bt11micron_track.append(b11)
        if avhrrObj.channels[4].data[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] == avhrrObj.nodata:
            b12 = -9.
        else:
            b12 = avhrrObj.channels[4].data[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] * avhrrObj.channels[4].gain + avhrrObj.channels[4].intercept
        bt12micron_track.append(b12)
        if ctth == None:
            continue
        if ctth.height[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] == ctth.h_nodata:
            hh = -9.
        else:
            hh = ctth.height[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] * ctth.h_gain + ctth.h_intercept
        ctth_height_track.append(hh)
        if ctth.temperature[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] == ctth.t_nodata:
            tt = -9.
        else:
            tt = ctth.temperature[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] * ctth.t_gain + \
                 ctth.t_intercept
        ctth_temperature_track.append(tt)
        if ctth.pressure[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] == ctth.p_nodata:
            pp = -9.
        else:
            pp = ctth.pressure[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] * ctth.p_gain + ctth.p_intercept
        ctth_pressure_track.append(pp)
        if avhrr5kmAngObj.satz.data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == avhrr5kmAngObj.satz.no_data or avhrr5kmAngObj.satz.data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == avhrr5kmAngObj.satz.missing_data:
            ang = -9
        else:
            ang = avhrr5kmAngObj.satz.data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] * avhrr5kmAngObj.satz.gain + avhrr5kmAngObj.satz.intercept
        satz_track.append(ang)
    retv5kmcwc.avhrr5km.latitude = numpy.array(lat_avhrr_track)
    retv5kmcwc.avhrr5km.longitude = numpy.array(lon_avhrr_track)
    retv5kmcwc.avhrr5km.cloudtype = numpy.array(ctype_track)
    retv5kmcwc.avhrr5km.bt11micron = numpy.array(bt11micron_track)
    retv5kmcwc.avhrr5km.bt12micron = numpy.array(bt12micron_track)
    retv5kmcwc.avhrr5km.satz = numpy.array(satz_track)
    if ctth:
        retv5kmcwc.avhrr5km.ctth_height = numpy.array(ctth_height_track)
        retv5kmcwc.avhrr5km.ctth_pressure = numpy.array(ctth_pressure_track)
        retv5kmcwc.avhrr5km.ctth_temperature = numpy.array(ctth_temperature_track)
    retv5kmcwc.avhrr5km.surftemp = numpy.array(surft_track)

    print "AVHRR-PPS Cloud Type,latitude: shapes = ",\
          retv5kmcwc.avhrr5km.cloudtype.shape,retv5kmcwc.avhrr5km.latitude.shape

# CWC-RVOD
    # ====

    #retv5kmcwc.cloudsat5kmcwc.TAI_start=cloudsat5kmObjcwc.TAI_start
    #retv5kmcwc.cloudsat5kmcwc.Profile_time=numpy.repeat(ProfileTime_cloudsatcwc,idx_matchcwc)
    retv5kmcwc.cloudsat5kmcwc.sec_1970=numpy.repeat(cloudsat5kmObjcwc.sec1970,idx_matchcwc)
    retv5kmcwc.cloudsat5kmcwc.latitude=numpy.repeat(latCloudsat5kmcwc,idx_matchcwc)
    retv5kmcwc.cloudsat5kmcwc.longitude=numpy.repeat(lonCloudsat5kmcwc,idx_matchcwc)

    x = numpy.repeat(cloudsat5kmObjcwc.Height[::,0],idx_matchcwc) # copies the first column
    col_range5kmcwc = cloudsat5kmObjcwc.Height.shape[1]
    for i in range(1,col_range5kmcwc):
        x = numpy.concatenate((x,numpy.repeat(cloudsat5kmObjcwc.Height[::,i],idx_matchcwc))) # Adds all columns to one long column
    N = x.shape[0]/col_range5kmcwc 		# Finds how many columns there should be
    retv5kmcwc.cloudsat5kmcwc.Height = numpy.reshape(x,(col_range5kmcwc,N)).astype('d')	# Reshape the long x to a 125xN matrix 

    retv5kmcwc.cloudsat5kmcwc.elevation = numpy.repeat(\
        cloudsat5kmObjcwc.elevation.ravel(),idx_matchcwc.ravel()).astype('d')
    retv5kmcwc.cloudsat5kmcwc.Data_quality=cloudsat5kmObjcwc.Data_quality
    retv5kmcwc.cloudsat5kmcwc.Data_targetID=cloudsat5kmObjcwc.Data_targetID
    
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_path = numpy.repeat(cloudsat5kmObjcwc.RVOD_liq_water_path,idx_matchcwc).astype('d')
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_path_uncertainty = numpy.repeat(\
        cloudsat5kmObjcwc.RVOD_liq_water_path_uncertainty,idx_matchcwc).astype('d')
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_path = numpy.repeat(cloudsat5kmObjcwc.RVOD_ice_water_path,idx_matchcwc).astype('d')
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_path_uncertainty = numpy.repeat(\
        cloudsat5kmObjcwc.RVOD_ice_water_path_uncertainty,idx_matchcwc).astype('d')
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_path = numpy.repeat(cloudsat5kmObjcwc.LO_RVOD_liquid_water_path,idx_matchcwc).astype('d')
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_path_uncertainty = numpy.repeat(\
        cloudsat5kmObjcwc.LO_RVOD_liquid_water_path_uncertainty,idx_matchcwc).astype('d')
    retv5kmcwc.cloudsat5kmcwc.IO_RVOD_ice_water_path = numpy.repeat(cloudsat5kmObjcwc.IO_RVOD_ice_water_path,idx_matchcwc).astype('d')
    retv5kmcwc.cloudsat5kmcwc.IO_RVOD_ice_water_path_uncertainty = numpy.repeat(\
        cloudsat5kmObjcwc.IO_RVOD_ice_water_path_uncertainty,idx_matchcwc).astype('d')

    x11 = numpy.repeat(cloudsat5kmObjcwc.RVOD_liq_water_content[::,0],idx_matchcwc)
    x12 = numpy.repeat(cloudsat5kmObjcwc.LO_RVOD_liquid_water_content[::,0],idx_matchcwc)
    x13 = numpy.repeat(cloudsat5kmObjcwc.RVOD_ice_water_content[::,0],idx_matchcwc)
    x14 = numpy.repeat(cloudsat5kmObjcwc.IO_RVOD_ice_water_content[::,0],idx_matchcwc)
    for i in range(1,col_range5kmcwc):
        x11 = numpy.concatenate((x11,numpy.repeat(cloudsat5kmObjcwc.RVOD_liq_water_content[::,i],idx_matchcwc))) 
        x12 = numpy.concatenate((x12,numpy.repeat(cloudsat5kmObjcwc.LO_RVOD_liquid_water_content[::,i],idx_matchcwc))) 
        x13 = numpy.concatenate((x13,numpy.repeat(cloudsat5kmObjcwc.RVOD_ice_water_content[::,i],idx_matchcwc))) 
        x14 = numpy.concatenate((x14,numpy.repeat(cloudsat5kmObjcwc.IO_RVOD_ice_water_content[::,i],idx_matchcwc))) 

    N11 = x11.shape[0]/col_range5kmcwc
    N12 = x12.shape[0]/col_range5kmcwc
    N13 = x13.shape[0]/col_range5kmcwc
    N14 = x14.shape[0]/col_range5kmcwc
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_content = numpy.reshape(x11,(col_range5kmcwc,N11)).astype('d')
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_content = numpy.reshape(x12,(col_range5kmcwc,N12)).astype('d')
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_content = numpy.reshape(x13,(col_range5kmcwc,N13)).astype('d')
    retv5kmcwc.cloudsat5kmcwc. IO_RVOD_ice_water_content= numpy.reshape(x14,(col_range5kmcwc,N14)).astype('d')

    x21 = numpy.repeat(cloudsat5kmObjcwc.RVOD_liq_water_content_uncertainty[::,0],idx_matchcwc)
    x22 = numpy.repeat(cloudsat5kmObjcwc.LO_RVOD_liquid_water_content_uncertainty[::,0],idx_matchcwc)
    x23 = numpy.repeat(cloudsat5kmObjcwc.RVOD_ice_water_content_uncertainty[::,0],idx_matchcwc)
    x24 = numpy.repeat(cloudsat5kmObjcwc.IO_RVOD_ice_water_content_uncertainty[::,0],idx_matchcwc)
    for i in range(1,col_range5kmcwc):
        x21 = numpy.concatenate((x21,numpy.repeat(cloudsat5kmObjcwc.RVOD_liq_water_content_uncertainty[::,i],idx_matchcwc))) 
        x22 = numpy.concatenate((x22,numpy.repeat(cloudsat5kmObjcwc.LO_RVOD_liquid_water_content_uncertainty[::,i],idx_matchcwc))) 
        x23 = numpy.concatenate((x23,numpy.repeat(cloudsat5kmObjcwc.RVOD_ice_water_content_uncertainty[::,i],idx_matchcwc))) 
        x24 = numpy.concatenate((x24,numpy.repeat(cloudsat5kmObjcwc.IO_RVOD_ice_water_content_uncertainty[::,i],idx_matchcwc))) 

    N21 = x21.shape[0]/col_range5kmcwc	
    N22 = x22.shape[0]/col_range5kmcwc
    N23 = x23.shape[0]/col_range5kmcwc
    N24 = x24.shape[0]/col_range5kmcwc
    retv5kmcwc.cloudsat5kmcwc.RVOD_liq_water_content_uncertainty = numpy.reshape(x21,(col_range5kmcwc,N21)).astype('d')
    retv5kmcwc.cloudsat5kmcwc.LO_RVOD_liquid_water_content_uncertainty = numpy.reshape(x22,(col_range5kmcwc,N22)).astype('d')
    retv5kmcwc.cloudsat5kmcwc.RVOD_ice_water_content_uncertainty = numpy.reshape(x23,(col_range5kmcwc,N23)).astype('d')
    retv5kmcwc.cloudsat5kmcwc. IO_RVOD_ice_water_content_uncertainty = numpy.reshape(x24,(col_range5kmcwc,N24)).astype('d')

    retv5kmcwc.cloudsat5kmcwc.Temp_min_mixph_K = cloudsat5kmObjcwc.Temp_min_mixph_K
    retv5kmcwc.cloudsat5kmcwc.Temp_max_mixph_K = cloudsat5kmObjcwc.Temp_max_mixph_K
    
    retv5kmcwc.cloudsat5kmcwc.avhrr_linnum = cal_on_avhrrcwc.astype('i')
    retv5kmcwc.cloudsat5kmcwc.avhrr_pixnum = cap_on_avhrrcwc.astype('i')

    retv5kmcwc.avhrr5km.sec_1970 = numpy.add(avhrrGeoObj.sec1970_start,cal_on_avhrrcwc * DSEC_PER_AVHRR_SCALINE5KM)
    retv5kmcwc.diff_sec_1970 = retv5kmcwc.cloudsat5kmcwc.sec_1970 - retv5kmcwc.avhrr5km.sec_1970
    min_diff = numpy.minimum.reduce(retv5kmcwc.diff_sec_1970)
    max_diff = numpy.maximum.reduce(retv5kmcwc.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-cloudsat): ",\
          numpy.maximum.reduce(retv5kmcwc.diff_sec_1970),numpy.minimum.reduce(retv5kmcwc.diff_sec_1970)

    print "AVHRR observation time of first cloudsat-avhrr match: ",\
          time.gmtime(retv5kmcwc.avhrr5km.sec_1970[0])
    print "AVHRR observation time of last cloudsat-avhrr match: ",\
          time.gmtime(retv5kmcwc.avhrr5km.sec_1970[N-1])

    ll = []
    for i in range(ndim5kmcwc):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat5kmcwc[i],latCloudsat5kmcwc[i],idx_matchcwc[i])))

    basename = os.path.basename(ctypefile).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_") 
    
    datapath = "%s/%s/5km/%s/%s/%s" %(DATA_DIR, base_sat, base_year, base_month, AREA5KM)    
    if not os.path.exists(datapath):
        os.makedirs(datapath)    
    fd = open("%s/5km_%s_cloudtype_cloudsat-CWC-RVOD_track2.txt"%(datapath,basename),"w")
    fd.writelines(ll)
    fd.close()

    ll = []
    for i in range(N):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_cloudsatcwc[i],lat_cloudsatcwc[i],0)))
    fd = open("%s/5km_%s_cloudtype_cloudsat-CWC-RVOD_track_excl.txt"%(datapath,basename),"w")
    fd.writelines(ll)
    fd.close()

    return retv5kmcwc,min_diff,max_diff
    
#------------------------------------------------------------------------------------------------------------

def reshapeCloudsat5kmCwc(cloudsatfiles,avhrr):
    import time
    import numpy
    import sys
    
    clsat = Cloudsat5kmCwcRvodObject
    if avhrr.sec1970_end<avhrr.sec1970_start:
        avhrr_end = int(DSEC_PER_AVHRR_SCALINE5KM*avhrr.num_of_lines+avhrr.sec1970_start)
    else:
        avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    
    startCloudsat = get_cloudsat5kmCwc(cloudsatfiles[0])
    startCloudsat.Profile_time = numpy.add(startCloudsat.Profile_time,startCloudsat.TAI_start)
    
    for i in range(len(cloudsatfiles)-1):
        newCloudsat = get_cloudsat5kmCwc(cloudsatfiles[i+1])
        newCloudsat.Profile_time = numpy.add(newCloudsat.Profile_time,newCloudsat.TAI_start)
        
        clsat_start_all = startCloudsat.sec1970.ravel()
        clsat_new_all = newCloudsat.sec1970.ravel()
        
        if not clsat_start_all[0]<clsat_new_all[0]:
            print "cloudsat files are in the wrong order"
            print("Program cloudsat5km.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9)
        clsat_break = numpy.argmin(numpy.abs(clsat_start_all - clsat_new_all[0]))+1

        
        startCloudsat.sec1970 = numpy.concatenate((startCloudsat.sec1970[0:clsat_break],newCloudsat.sec1970))
        startCloudsat.longitude = numpy.concatenate((startCloudsat.longitude[0:clsat_break,:],newCloudsat.longitude))
        startCloudsat.latitude = numpy.concatenate((startCloudsat.latitude[0:clsat_break,:],newCloudsat.latitude))
        startCloudsat.elevation = numpy.concatenate((startCloudsat.elevation[0:clsat_break,:],newCloudsat.elevation))
        startCloudsat.Profile_time = numpy.concatenate((startCloudsat.Profile_time[0:clsat_break,:],newCloudsat.Profile_time))
        startCloudsat.Height = numpy.concatenate((startCloudsat.Height[0:clsat_break,:],newCloudsat.Height))
        startCloudsat.Data_quality = numpy.concatenate((startCloudsat.Data_quality[0:clsat_break,:],newCloudsat.Data_quality))
        startCloudsat.Data_targetID = numpy.concatenate((startCloudsat.Data_targetID[0:clsat_break,:],newCloudsat.Data_targetID))
        startCloudsat.RVOD_liq_water_path = numpy.concatenate((startCloudsat.RVOD_liq_water_path[0:clsat_break,:],newCloudsat.RVOD_liq_water_path))
        startCloudsat.RVOD_liq_water_path_uncertainty = numpy.concatenate((startCloudsat.RVOD_liq_water_path_uncertainty[0:clsat_break,:],newCloudsat.RVOD_liq_water_path_uncertainty))
        startCloudsat.RVOD_ice_water_path = numpy.concatenate((startCloudsat.RVOD_ice_water_path[0:clsat_break,:],newCloudsat.RVOD_ice_water_path))
        startCloudsat.RVOD_ice_water_path_uncertainty = numpy.concatenate((startCloudsat.RVOD_ice_water_path_uncertainty[0:clsat_break,:],newCloudsat.RVOD_ice_water_path_uncertainty))
        startCloudsat.LO_RVOD_liquid_water_path = numpy.concatenate((startCloudsat.LO_RVOD_liquid_water_path[0:clsat_break,:],newCloudsat.LO_RVOD_liquid_water_path))
        startCloudsat.LO_RVOD_liquid_water_path_uncertainty = numpy.concatenate((startCloudsat.LO_RVOD_liquid_water_path_uncertainty[0:clsat_break,:],newCloudsat.LO_RVOD_liquid_water_path_uncertainty))
        startCloudsat.IO_RVOD_ice_water_path = numpy.concatenate((startCloudsat.IO_RVOD_ice_water_path[0:clsat_break,:],newCloudsat.IO_RVOD_ice_water_path))
        startCloudsat.IO_RVOD_ice_water_path_uncertainty = numpy.concatenate((startCloudsat.IO_RVOD_ice_water_path_uncertainty[0:clsat_break,:],newCloudsat.IO_RVOD_ice_water_path_uncertainty))
        startCloudsat.RVOD_liq_water_content = numpy.concatenate((startCloudsat.RVOD_liq_water_content[0:clsat_break,:],newCloudsat.RVOD_liq_water_content))
        startCloudsat.RVOD_liq_water_content_uncertainty = numpy.concatenate((startCloudsat.RVOD_liq_water_content_uncertainty[0:clsat_break,:],newCloudsat.RVOD_liq_water_content_uncertainty))
        startCloudsat.RVOD_ice_water_content = numpy.concatenate((startCloudsat.RVOD_ice_water_content[0:clsat_break,:],newCloudsat.RVOD_ice_water_content))
        startCloudsat.RVOD_ice_water_content_uncertainty = numpy.concatenate((startCloudsat.RVOD_ice_water_content_uncertainty[0:clsat_break,:],newCloudsat.RVOD_ice_water_content_uncertainty))
        startCloudsat.LO_RVOD_liquid_water_content = numpy.concatenate((startCloudsat.LO_RVOD_liquid_water_content[0:clsat_break,:],newCloudsat.LO_RVOD_liquid_water_content))
        startCloudsat.LO_RVOD_liquid_water_content_uncertainty = numpy.concatenate((startCloudsat.LO_RVOD_liquid_water_content_uncertainty[0:clsat_break,:],newCloudsat.LO_RVOD_liquid_water_content_uncertainty))
        startCloudsat.IO_RVOD_ice_water_content = numpy.concatenate((startCloudsat.IO_RVOD_ice_water_content[0:clsat_break,:],newCloudsat.IO_RVOD_ice_water_content))
        startCloudsat.IO_RVOD_ice_water_content_uncertainty = numpy.concatenate((startCloudsat.IO_RVOD_ice_water_content_uncertainty[0:clsat_break,:],newCloudsat.IO_RVOD_ice_water_content_uncertainty))

    start_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_start - sec_timeThr))))-1 # Minus one to get one extra, just to be certain
    end_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain

    clsat.sec1970 = startCloudsat.sec1970[start_break:end_break] 
    clsat.longitude = startCloudsat.longitude[start_break:end_break,:] 
    clsat.latitude = startCloudsat.latitude[start_break:end_break,:] 
    clsat.elevation = startCloudsat.elevation[start_break:end_break,:] 
    clsat.Profile_time = startCloudsat.Profile_time[start_break:end_break,:]  
    clsat.Data_quality = startCloudsat.Data_quality[start_break:end_break,:] 
    clsat.Data_targetID = startCloudsat.Data_targetID[start_break:end_break,:] 
    clsat.Height = startCloudsat.Height[start_break:end_break,:] 
    clsat.RVOD_liq_water_path = startCloudsat.RVOD_liq_water_path[start_break:end_break,:]
    clsat.RVOD_liq_water_path_uncertainty = startCloudsat.RVOD_liq_water_path_uncertainty[start_break:end_break,:]
    clsat.RVOD_ice_water_path = startCloudsat.RVOD_ice_water_path[start_break:end_break,:]
    clsat.RVOD_ice_water_path_uncertainty = startCloudsat.RVOD_ice_water_path_uncertainty[start_break:end_break,:]
    clsat.LO_RVOD_liquid_water_path = startCloudsat.LO_RVOD_liquid_water_path[start_break:end_break,:]
    clsat.LO_RVOD_liquid_water_path_uncertainty = startCloudsat.LO_RVOD_liquid_water_path_uncertainty[start_break:end_break,:]
    clsat.IO_RVOD_ice_water_path = startCloudsat.IO_RVOD_ice_water_path[start_break:end_break,:]
    clsat.IO_RVOD_ice_water_path_uncertainty = startCloudsat.IO_RVOD_ice_water_path_uncertainty[start_break:end_break,:]
    clsat.RVOD_liq_water_content = startCloudsat.RVOD_liq_water_content[start_break:end_break,:]
    clsat.RVOD_liq_water_content_uncertainty = startCloudsat.RVOD_liq_water_content_uncertainty[start_break:end_break,:]
    clsat.RVOD_ice_water_content = startCloudsat.RVOD_ice_water_content[start_break:end_break,:]
    clsat.RVOD_ice_water_content_uncertainty = startCloudsat.RVOD_ice_water_content_uncertainty[start_break:end_break,:]
    clsat.LO_RVOD_liquid_water_content = startCloudsat.LO_RVOD_liquid_water_content[start_break:end_break,:]
    clsat.LO_RVOD_liquid_water_content_uncertainty = startCloudsat.LO_RVOD_liquid_water_content_uncertainty[start_break:end_break,:]
    clsat.IO_RVOD_ice_water_content = startCloudsat.IO_RVOD_ice_water_content[start_break:end_break,:]
    clsat.IO_RVOD_ice_water_content_uncertainty = startCloudsat.IO_RVOD_ice_water_content_uncertainty[start_break:end_break,:]
    clsat.Temp_min_mixph_K = startCloudsat.Temp_min_mixph_K
    clsat.Temp_max_mixph_K = startCloudsat.Temp_max_mixph_K
    if clsat.Profile_time.shape[0] <= 0:
        print("No time match, please try with some other CloudSat files")
        print("Program cloudsat5km.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)
    else:
        clsat.TAI_start = clsat.Profile_time[0,0]
        clsat.Profile_time = clsat.Profile_time - clsat.TAI_start
    
    return clsat
