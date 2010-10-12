#HISTORY CHANGES BY ANKE TETZLAFF#

#080430:
# Could not run in the Arctic using the tmpaid creation
# Used hard coded area 'arctic_super_5010' instead

#080416: Got error message in line 232 and 234:
#         "only rank-0 arrays can be converted to Python scalars";
#         changed start_sec1970 and end_sec1970 to rank0 array 
#         by setting start_sec1970[0] and end_sec1970[0]

from pps_basic_configure import *
from pps_error_messages import *

from calipso import *
from cloudsat_calipso_avhrr_match import AREA1KM, SUB_DIR, DATA_DIR, sec_timeThr

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

COVERAGE_DIR = "%s/1km/%s"%(SUB_DIR,AREA1KM)

#NLINES=1000
#NLINES=6000
#SWATHWD=2048
#NODATA=-9


class CloudsatCwcRvodObject:
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

class CloudsatCwcAvhrrTrackObject:
    def __init__(self):
        self.avhrr=ppsAvhrrObject()
        self.cloudsatcwc=CloudsatCwcRvodObject()
        self.diff_sec_1970=None
##########################################################################################################################################   
#pdb.set_trace()

# ----------------------------------------
def readCloudsatCwcAvhrrMatchObj(filename):
    import _pyhl

    retvcwc = CloudsatCwcAvhrrTrackObject()

    a=_pyhl.read_nodelist(filename)
    a.selectAll()
    a.fetch()

    # Match-Up - time difference:
    retvcwc.diff_sec_1970 = a.getNode("/diff_sec_1970").data()

    # Cloudsat:

    # AVHRR:
    retvcwc.avhrr.longitude = a.getNode("/avhrr/longitude").data()
    retvcwc.avhrr.latitude = a.getNode("/avhrr/latitude").data()
    retvcwc.avhrr.sec_1970 = a.getNode("/avhrr/sec_1970").data()
    retvcwc.avhrr.cloudtype = a.getNode("/avhrr/cloudtype").data()
    retvcwc.avhrr.ctth_height = a.getNode("/avhrr/ctth_height").data()
    retvcwc.avhrr.ctth_pressure = a.getNode("/avhrr/ctth_pressure").data()
    retvcwc.avhrr.ctth_temperature = a.getNode("/avhrr/ctth_temperature").data()
    retvcwc.avhrr.bt11micron = a.getNode("/avhrr/bt11micron").data()
    retvcwc.avhrr.bt12micron = a.getNode("/avhrr/bt12micron").data()
    retvcwc.avhrr.surftemp = a.getNode("/avhrr/surftemp").data()
    retvcwc.avhrr.satz = a.getNode("/avhrr/satz").data()
    # CWC-RVOD:
    # ====
# Geolocation
    retvcwc.cloudsatcwc.longitude = a.getNode("/cloudsatcwc/longitude").data()
    retvcwc.cloudsatcwc.latitude = a.getNode("/cloudsatcwc/latitude").data()
    retvcwc.cloudsatcwc.elevation = a.getNode("/cloudsatcwc/elevation").data()
    retvcwc.cloudsatcwc.Height = a.getNode("/cloudsatcwc/Height").data()
    retvcwc.cloudsatcwc.avhrr_linnum = a.getNode("/cloudsatcwc/avhrr_linnum").data()
    retvcwc.cloudsatcwc.avhrr_pixnum = a.getNode("/cloudsatcwc/avhrr_pixnum").data()
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    retvcwc.cloudsatcwc.Profile_time = a.getNode("/cloudsatcwc/Profile_time").data()
    #retvcwc.cloudsatcwc.TAI_start = a.getNode("/cloudsatcwc/TAI_start").data()
    retvcwc.cloudsatcwc.sec_1970 = a.getNode("/cloudsatcwc/sec_1970").data()
# The data
    retvcwc.cloudsatcwc.Data_quality = a.getNode("/cloudsatcwc/Data_quality").data()
    retvcwc.cloudsatcwc.Data_targetID = a.getNode("/cloudsatcwc/Data_targetID").data()
    retvcwc.cloudsatcwc.RVOD_liq_water_path = a.getNode("/cloudsatcwc/RVOD_liq_water_path").data()
    retvcwc.cloudsatcwc.RVOD_liq_water_path_uncertainty = a.getNode("/cloudsatcwc/RVOD_liq_water_path_uncertainty").data()
    retvcwc.cloudsatcwc.RVOD_ice_water_path = a.getNode("/cloudsatcwc/RVOD_ice_water_path").data()
    retvcwc.cloudsatcwc.RVOD_ice_water_path_uncertainty = a.getNode("/cloudsatcwc/RVOD_ice_water_path_uncertainty").data()
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_path = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_path").data()
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_path_uncertainty = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_path_uncertainty").data()
    retvcwc.cloudsatcwc.IO_RVOD_ice_water_path = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_path").data()
    retvcwc.cloudsatcwc.IO_RVOD_ice_water_path_uncertainty = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_path_uncertainty").data()

   
    retvcwc.cloudsatcwc.RVOD_liq_water_content = a.getNode("/cloudsatcwc/RVOD_liq_water_content").data()
    retvcwc.cloudsatcwc.RVOD_liq_water_content_uncertainty = a.getNode("/cloudsatcwc/RVOD_liq_water_content_uncertainty").data()
    retvcwc.cloudsatcwc.RVOD_ice_water_content = a.getNode("/cloudsatcwc/RVOD_ice_water_content").data()
    retvcwc.cloudsatcwc.RVOD_ice_water_content_uncertainty = a.getNode("/cloudsatcwc/RVOD_ice_water_content_uncertainty").data()
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_content = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_content").data()
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_content_uncertainty = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_content_uncertainty").data()
    retvcwc.cloudsatcwc.IO_RVOD_ice_water_content = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_content").data()
    retvcwc.cloudsatcwc.IO_RVOD_ice_water_content_uncertainty = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_content_uncertainty").data()
    retvcwc.cloudsatcwc.Temp_min_mixph_K = a.getNode("/cloudsatcwc/Temp_min_mixph_K").data()
    retvcwc.cloudsatcwc.Temp_max_mixph_K = a.getNode("/cloudsatcwc/Temp_max_mixph_K").data()

##########################################################################################################################################
    #pdb.set_trace()

    return retvcwc


# ----------------------------------------
def writeCloudsatCwcAvhrrMatchObj(filename,ca_obj,compress_lvl):
    import _pyhl
    status = -1
    
    a=_pyhl.nodelist()

    shape = [ca_obj.cloudsatcwc.longitude.shape[0]]

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
    shape = ca_obj.avhrr.longitude.shape
    b.setArrayValue(1,shape,ca_obj.avhrr.longitude,"float",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/latitude")
    b.setArrayValue(1,shape,ca_obj.avhrr.latitude,"float",-1)
    a.addNode(b)

    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/sec_1970")
    b.setArrayValue(1,shape,ca_obj.avhrr.sec_1970,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_pressure")
    b.setArrayValue(1,shape,ca_obj.avhrr.ctth_pressure,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_temperature")
    b.setArrayValue(1,shape,ca_obj.avhrr.ctth_temperature,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_height")
    b.setArrayValue(1,shape,ca_obj.avhrr.ctth_height,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/cloudtype")
    b.setArrayValue(1,shape,ca_obj.avhrr.cloudtype,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/bt11micron")
    b.setArrayValue(1,shape,ca_obj.avhrr.bt11micron,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/bt12micron")
    b.setArrayValue(1,shape,ca_obj.avhrr.bt12micron,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/surftemp")
    b.setArrayValue(1,shape,ca_obj.avhrr.surftemp,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/satz")
    b.setArrayValue(1,shape,ca_obj.avhrr.satz,"double",-1)
    a.addNode(b)

    # CWC-RVOD
    # ====
    shapecwc = [ca_obj.cloudsatcwc.longitude.shape[0]]
    shape2dcwc = ca_obj.cloudsatcwc.Height.shape
    #shapeTAIcwc = [ca_obj.cloudsatcwc.TAI_start.shape[0]]
    shapeTAIcwc = [ca_obj.cloudsatcwc.Temp_min_mixph_K.shape[0]]
    shapecwclong = [ca_obj.cloudsatcwc.Profile_time.shape[0]]
# Geolocation
    b=_pyhl.node(_pyhl.GROUP_ID,"/cloudsatcwc")
    a.addNode(b)
    
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/longitude")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.longitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/latitude")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.latitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/elevation")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.elevation,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Height")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.Height,"short",-1)
    a.addNode(b)       
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/avhrr_linnum")
    b.setArrayValue(1,shape,ca_obj.cloudsatcwc.avhrr_linnum,"int",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/avhrr_pixnum")
    b.setArrayValue(1,shape,ca_obj.cloudsatcwc.avhrr_pixnum,"int",-1)
    a.addNode(b)
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Profile_time")
    b.setArrayValue(1,shapecwclong,ca_obj.cloudsatcwc.Profile_time,"double",-1)
    a.addNode(b) 
    #b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/TAI_start")
    #b.setArrayValue(1,shapeTAIcwc,ca_obj.cloudsatcwc.TAI_start,"double",-1)
    #a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/sec_1970")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.sec_1970,"double",-1)
    a.addNode(b)
# The data
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Data_quality")
    b.setArrayValue(1,shapecwclong,ca_obj.cloudsatcwc.Data_quality,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Data_targetID")
    b.setArrayValue(1,shapecwclong,ca_obj.cloudsatcwc.Data_targetID,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.RVOD_liq_water_path,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.RVOD_liq_water_path_uncertainty,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.RVOD_ice_water_path,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.RVOD_ice_water_path_uncertainty,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.LO_RVOD_liquid_water_path,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.LO_RVOD_liquid_water_path_uncertainty,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_path")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.IO_RVOD_ice_water_path,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_path_uncertainty")
    b.setArrayValue(1,shapecwc,ca_obj.cloudsatcwc.IO_RVOD_ice_water_path_uncertainty,"uchar",-1)
    a.addNode(b)
##########################################################################################################################################
    #pdb.set_trace()
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.RVOD_liq_water_content,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_liq_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.RVOD_liq_water_content_uncertainty,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.RVOD_ice_water_content,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/RVOD_ice_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.RVOD_ice_water_content_uncertainty,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.LO_RVOD_liquid_water_content,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/LO_RVOD_liquid_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.LO_RVOD_liquid_water_content_uncertainty,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_content")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.IO_RVOD_ice_water_content,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/IO_RVOD_ice_water_content_uncertainty")
    b.setArrayValue(1,shape2dcwc,ca_obj.cloudsatcwc.IO_RVOD_ice_water_content_uncertainty,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Temp_min_mixph_K")
    b.setArrayValue(1,shapeTAIcwc,ca_obj.cloudsatcwc.Temp_min_mixph_K,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Temp_max_mixph_K")
    b.setArrayValue(1,shapeTAIcwc,ca_obj.cloudsatcwc.Temp_max_mixph_K,"double",-1)
    a.addNode(b)

    status = a.write(filename,compress_lvl)

    return status

# -----------------------------------------------------
def select_cloudsatCwc_inside_avhrr(cloudsatObj,cal,sec1970_start_end,sec_timeThr):
    # Cloudsat times are already in UTC in sec since 1970
    import numpy

    sec1970_start,sec1970_end = sec1970_start_end
    
    # Select the points inside the avhrr swath:
    # Allowing for sec_timeThr seconds deviation:
    idx_time_okay = numpy.logical_and(numpy.greater(\
        cloudsatObj.sec1970,sec1970_start - sec_timeThr),
                                   numpy.less(\
        cloudsatObj.sec1970,sec1970_end + sec_timeThr))
    
    #idx_match = numpy.not_equal(cal,NODATA)
    idx_place_okay = numpy.where(numpy.not_equal(cal,NODATA),idx_time_okay,False)
    idx_match = idx_place_okay
    return idx_match

# -----------------------------------------------------
def get_cloudsatCwc(filenamecwc):
    import _pypps_filters
    import time

    # Read CLOUDSAT Radar data:
    cloudsatcwc = read_cloudsatCwc(filenamecwc)
##########################################################################################################################################
    #pdb.set_trace()
    # CWC-RVOD
    # ====

    lonCloudsatcwc = cloudsatcwc.longitude.ravel()
    latCloudsatcwc = cloudsatcwc.latitude.ravel()
    ndimcwc = lonCloudsatcwc.shape[0]
    # Time in seconds since 1970:
    
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone

    print "CWC-RVOD TAI start time: ",cloudsatcwc.TAI_start
    start_sec1970cwc = cloudsatcwc.TAI_start + dsec
    start_timecwc = time.gmtime(start_sec1970cwc[0])
    end_sec1970cwc = cloudsatcwc.TAI_start + dsec + cloudsatcwc.Profile_time[ndimcwc-1]
    end_timecwc = time.gmtime(end_sec1970cwc[0])
    print "CWC-RVOD Start and end times: ",start_timecwc,end_timecwc
    cloudsatcwc.sec1970 = cloudsatcwc.Profile_time + start_sec1970cwc

    # Any other possible derived field should come here...

    # --------------------------------------------------------------------

    return cloudsatcwc

# -----------------------------------------------------
def read_cloudsatCwc(filenamecwc):
    import _pyhl
    import numpy.oldnumeric as Numeric
    import numpy
    
    # CWC-RVOD
    # ====
    acwc=_pyhl.read_nodelist(filenamecwc)
    bcwc=acwc.getNodeNames()
    acwc.selectAll()
    acwc.fetch()

    retvcwc = CloudsatCwcRvodObject()
    rootcwc="/2B-CWC-RVOD"
    
# Geolocation
    ccwc=acwc.getNode("%s/Geolocation Fields/Longitude"%rootcwc)
    retvcwc.longitude=numpy.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Geolocation Fields/Latitude"%rootcwc)
    retvcwc.latitude=numpy.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Geolocation Fields/DEM_elevation"%rootcwc)
    retvcwc.elevation=numpy.fromstring(ccwc.data(),'Int16').astype('Int16')
    ccwc=acwc.getNode("%s/Geolocation Fields/Height"%rootcwc)
    retvcwc.Height=ccwc.data()
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    ccwc=acwc.getNode("%s/Geolocation Fields/Profile_time"%rootcwc)
    retvcwc.Profile_time=numpy.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Geolocation Fields/TAI_start"%rootcwc)
    retvcwc.TAI_start=numpy.fromstring(ccwc.data(),'d').astype('d')
# The data
    ccwc=acwc.getNode("%s/Data Fields/Data_quality"%rootcwc)
    retvcwc.Data_quality=numpy.fromstring(ccwc.data(),'uint8').astype('uint8')
    ccwc=acwc.getNode("%s/Data Fields/Data_targetID"%rootcwc)
    retvcwc.Data_targetID=numpy.fromstring(ccwc.data(),'uint8').astype('uint8')
# IWP and LWP		
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_path"%rootcwc)
    retvcwc.RVOD_liq_water_path=numpy.fromstring(ccwc.data(),'Int16').astype('Int16')			
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_path_uncertainty"%rootcwc)
    retvcwc.RVOD_liq_water_path_uncertainty=numpy.fromstring(ccwc.data(),'uint8').astype('uint8')    
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_path"%rootcwc)
    retvcwc.RVOD_ice_water_path=numpy.fromstring(ccwc.data(),'Int16').astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_path_uncertainty"%rootcwc)
    retvcwc.RVOD_ice_water_path_uncertainty=numpy.fromstring(ccwc.data(),'uint8').astype('uint8')
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_path"%rootcwc)
    retvcwc.LO_RVOD_liquid_water_path=numpy.fromstring(ccwc.data(),'Int16').astype('Int16')			
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_path_uncertainty"%rootcwc)
    retvcwc.LO_RVOD_liquid_water_path_uncertainty=numpy.fromstring(ccwc.data(),'uint8').astype('uint8')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_path"%rootcwc)
    retvcwc.IO_RVOD_ice_water_path=numpy.fromstring(ccwc.data(),'Int16').astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_path_uncertainty"%rootcwc)
    retvcwc.IO_RVOD_ice_water_path_uncertainty=numpy.fromstring(ccwc.data(),'uint8').astype('uint8')
# IWC and LWC
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_content"%rootcwc)
    retvcwc.RVOD_liq_water_content=ccwc.data()#.astype('Int16')		
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_content_uncertainty"%rootcwc)
    retvcwc.RVOD_liq_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_content"%rootcwc)
    retvcwc.RVOD_ice_water_content=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_content_uncertainty"%rootcwc)
    retvcwc.RVOD_ice_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_content"%rootcwc)
    retvcwc.LO_RVOD_liquid_water_content=ccwc.data()#.astype('Int16')			
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_content_uncertainty"%rootcwc)
    retvcwc.LO_RVOD_liquid_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_content"%rootcwc)
    retvcwc.IO_RVOD_ice_water_content=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_content_uncertainty"%rootcwc)
    retvcwc.IO_RVOD_ice_water_content_uncertainty=ccwc.data()#.astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/Temp_min_mixph_K"%rootcwc)
    retvcwc.Temp_min_mixph_K=numpy.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Data Fields/Temp_max_mixph_K"%rootcwc)
    retvcwc.Temp_max_mixph_K=numpy.fromstring(ccwc.data(),'f').astype('d')
##########################################################################################################################################
    #pdb.set_trace()
    return retvcwc

#retvcwc.RVOD_liq_water_content.shape


# --------------------------------------------
def get_cloudsatCwc_avhrr_linpix(avhrrIn,lon,lat):
    #import numpy.oldnumeric as Numeric
    import numpy
    tmppcs="tmpproj"
    define_pcs(tmppcs, "Plate Caree, central meridian at 15E",
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
        coverage_filename = "%s/coverage_avhrr_cloudsat-CWC_RVOD_matchup_%s_%.5d_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrrIn.satellite,avhrrIn.orbit,
                             startline,endline,tmpaid)
        write_log("INFO","Coverage filename = ",coverage_filename)
        cal,cap,ok = get_cloudsatCwc_avhrr_linpix_segment(avhrrIn,lon,lat,
                                                      (startline,endline),
                                                      SWATHWD,tmppcs,tmpaid,
                                                      coverage_filename)
        if ok:
            HasEncounteredMatch=1
            write_log("INFO","There was a match...")
            
        if not ok and HasEncounteredMatch:
            write_log("INFO","Data is now empty. Leave the loop...")
            break
        
        if(startline==0):
            # First time:
            cloudsat_avhrr_line,cloudsat_avhrr_pixel = numpy.array(cal),numpy.array(cap)
        else:
            # Merge:
            cloudsat_avhrr_line = numpy.where(numpy.equal(cloudsat_avhrr_line,-9),cal,cloudsat_avhrr_line)
            cloudsat_avhrr_pixel = numpy.where(numpy.equal(cloudsat_avhrr_pixel,-9),cap,cloudsat_avhrr_pixel)

        startline=startline+NLINES
        i=i+1
        ##########################################################################################################################################
    #pdb.set_trace()
    return cloudsat_avhrr_line,cloudsat_avhrr_pixel

# --------------------------------------------
def get_cloudsatCwc_avhrr_linpix_segment(avhrrIn,lon,lat,lines,swath_width,tmppcs,
                                      tmpaid,covfilename):
    #import numpy.oldnumeric as Numeric
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
    areaObj = area.area(AREA1KM)

    
    if not os.path.exists(covfilename):
        write_log("INFO","Create Coverage map...")
        cov = _satproj.create_coverage(areaObj,lonarr,latarr,1)
        print covfilename
        writeCoverage(cov,covfilename,"satproj",AREA1KM)
    else:
        write_log("INFO","Read the AVHRR-CLOUDSAT matchup coverage from file...")
        cov,info = readCoverage(covfilename)
    
    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA)
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA)


    write_log("INFO","Go through cloudsat track:")
    cloudsat_avhrr_line = []
    cloudsat_avhrr_pixel = []
    for i in range(ndim):
        xy_tup=pps_gisdata.lonlat2xy(AREA1KM,lon[i],lat[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
##        dimx=4500 Should be 5010!!!/KG
##        dimy=4500
        #dimx=5010
        #dimy=5010
        dimx=mapped_line.shape[1]#1002#5010
        dimy=mapped_line.shape[0]#1002#5010
        
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            cloudsat_avhrr_line.append(mapped_line[y,x])
            cloudsat_avhrr_pixel.append(mapped_pixel[y,x])
        else:
            cloudsat_avhrr_line.append(-9)
            cloudsat_avhrr_pixel.append(-9)

    cloudsat_avhrr_line = numpy.array(cloudsat_avhrr_line)
    cloudsat_avhrr_pixel = numpy.array(cloudsat_avhrr_pixel)

    x=numpy.repeat(cloudsat_avhrr_line,numpy.not_equal(cloudsat_avhrr_line,-9))
    write_log("INFO","Number of matching points = ",x.shape[0])
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0

    return cloudsat_avhrr_line,cloudsat_avhrr_pixel,matchOk

# -----------------------------------------------------
def match_cloudsatCwc_avhrr(ctypefile,cloudsatObjcwc,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj):
    #import numpy.oldnumeric as Numeric
    import time
    import string
    import numpy
    retvcwc = CloudsatCwcAvhrrTrackObject()

    lonCloudsatcwc = cloudsatObjcwc.longitude.ravel()		# Cwc
    latCloudsatcwc = cloudsatObjcwc.latitude.ravel()		# Cwc
    ndimcwc = lonCloudsatcwc.shape[0]				# Cwc

    # --------------------------------------------------------------------

    calcwc,capcwc=get_cloudsatCwc_avhrr_linpix(avhrrGeoObj,lonCloudsatcwc,latCloudsatcwc)	# Cwc

    #print len(cal), len(cap), ndim

    idx_cloudsatcwc = numpy.arange(ndimcwc)							# Cwc
    idx_cloudsatcwc_on_avhrr=numpy.repeat(idx_cloudsatcwc,numpy.not_equal(calcwc,NODATA)) 	# Cwc   
    lon_cloudsatcwc = numpy.repeat(lonCloudsatcwc,numpy.not_equal(calcwc,NODATA))		# Cwc
    lat_cloudsatcwc = numpy.repeat(latCloudsatcwc,numpy.not_equal(calcwc,NODATA))		# Cwc
    # Cloudsat line,pixel inside AVHRR swath:

    cal_on_avhrrcwc = numpy.repeat(calcwc,numpy.not_equal(calcwc,NODATA))	# Cwc
    cap_on_avhrrcwc = numpy.repeat(capcwc,numpy.not_equal(calcwc,NODATA))	# Cwc
    timetup_startcwc = time.gmtime(avhrrGeoObj.sec1970_start)		# Cwc
    timetup_endcwc   = time.gmtime(avhrrGeoObj.sec1970_end)		# Cwc

    #sec_timeThrcwc = 60*20 # Allow for 20 minute deviation between AVHRR and CLOUDSAT Cwcmatchup
    secTupcwc = avhrrGeoObj.sec1970_start,avhrrGeoObj.sec1970_end
    idx_matchcwc = select_cloudsatCwc_inside_avhrr(cloudsatObjcwc,calcwc,secTupcwc,sec_timeThr)	########kanske blir fel??????????


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
        if avhrrAngObj.satz.data[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] == avhrrAngObj.satz.no_data or avhrrAngObj.satz.data[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] == avhrrAngObj.satz.missing_data:
            ang = -9
        else:
            ang = avhrrAngObj.satz.data[cal_on_avhrrcwc[i],cap_on_avhrrcwc[i]] * avhrrAngObj.satz.gain + avhrrAngObj.satz.intercept
        satz_track.append(ang)
        
    retvcwc.avhrr.latitude = numpy.array(lat_avhrr_track)
    retvcwc.avhrr.longitude = numpy.array(lon_avhrr_track)
    retvcwc.avhrr.cloudtype = numpy.array(ctype_track)
    retvcwc.avhrr.bt11micron = numpy.array(bt11micron_track)
    retvcwc.avhrr.bt12micron = numpy.array(bt12micron_track)
    retvcwc.avhrr.satz = numpy.array(satz_track)
    if ctth:
        retvcwc.avhrr.ctth_height = numpy.array(ctth_height_track)
        retvcwc.avhrr.ctth_pressure = numpy.array(ctth_pressure_track)
        retvcwc.avhrr.ctth_temperature = numpy.array(ctth_temperature_track)
    retvcwc.avhrr.surftemp = numpy.array(surft_track)

    print "AVHRR-PPS Cloud Type,latitude: shapes = ",\
          retvcwc.avhrr.cloudtype.shape,retvcwc.avhrr.latitude.shape

# CWC-RVOD
    # ====
    retvcwc.cloudsatcwc.TAI_start=cloudsatObjcwc.TAI_start
    retvcwc.cloudsatcwc.Profile_time=cloudsatObjcwc.Profile_time
    retvcwc.cloudsatcwc.sec_1970=numpy.repeat(cloudsatObjcwc.sec1970,idx_matchcwc)
    retvcwc.cloudsatcwc.latitude=numpy.repeat(cloudsatObjcwc.latitude,idx_matchcwc)
    retvcwc.cloudsatcwc.longitude=numpy.repeat(cloudsatObjcwc.longitude,idx_matchcwc)

    x = numpy.repeat(cloudsatObjcwc.Height[::,0],idx_matchcwc) # copies the first column
    for i in range(1,125):
        x = numpy.concatenate((x,numpy.repeat(cloudsatObjcwc.Height[::,i],idx_matchcwc))) # Adds all columns to one long column
    N = x.shape[0]/125 		# Finds how many columns there should be
    retvcwc.cloudsatcwc.Height = numpy.reshape(x,(125,N)).astype('Int16')	# Reshape the long x to a 125xN matrix 

    retvcwc.cloudsatcwc.elevation = numpy.repeat(\
        cloudsatObjcwc.elevation.ravel(),idx_matchcwc.ravel()).astype('Int16')
    retvcwc.cloudsatcwc.Data_quality=cloudsatObjcwc.Data_quality
    retvcwc.cloudsatcwc.Data_targetID=cloudsatObjcwc.Data_targetID
    
    retvcwc.cloudsatcwc.RVOD_liq_water_path = numpy.repeat(cloudsatObjcwc.RVOD_liq_water_path,idx_matchcwc).astype('Int16')
    retvcwc.cloudsatcwc.RVOD_liq_water_path_uncertainty = numpy.repeat(\
        cloudsatObjcwc.RVOD_liq_water_path_uncertainty,idx_matchcwc).astype('Int8')
    retvcwc.cloudsatcwc.RVOD_ice_water_path = numpy.repeat(cloudsatObjcwc.RVOD_ice_water_path,idx_matchcwc).astype('Int16')
    retvcwc.cloudsatcwc.RVOD_ice_water_path_uncertainty = numpy.repeat(\
        cloudsatObjcwc.RVOD_ice_water_path_uncertainty,idx_matchcwc).astype('Int8')
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_path = numpy.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_path,idx_matchcwc).astype('Int16')
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_path_uncertainty = numpy.repeat(\
        cloudsatObjcwc.LO_RVOD_liquid_water_path_uncertainty,idx_matchcwc).astype('Int8')
    retvcwc.cloudsatcwc.IO_RVOD_ice_water_path = numpy.repeat(cloudsatObjcwc.IO_RVOD_ice_water_path,idx_matchcwc).astype('Int16')
    retvcwc.cloudsatcwc.IO_RVOD_ice_water_path_uncertainty = numpy.repeat(\
        cloudsatObjcwc.IO_RVOD_ice_water_path_uncertainty,idx_matchcwc).astype('Int8')

    x11 = numpy.repeat(cloudsatObjcwc.RVOD_liq_water_content[::,0],idx_matchcwc)
    x12 = numpy.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content[::,0],idx_matchcwc)
    x13 = numpy.repeat(cloudsatObjcwc.RVOD_ice_water_content[::,0],idx_matchcwc)
    x14 = numpy.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content[::,0],idx_matchcwc)
    for i in range(1,125):
        x11 = numpy.concatenate((x11,numpy.repeat(cloudsatObjcwc.RVOD_liq_water_content[::,i],idx_matchcwc))) 
        x12 = numpy.concatenate((x12,numpy.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content[::,i],idx_matchcwc))) 
        x13 = numpy.concatenate((x13,numpy.repeat(cloudsatObjcwc.RVOD_ice_water_content[::,i],idx_matchcwc))) 
        x14 = numpy.concatenate((x14,numpy.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content[::,i],idx_matchcwc))) 

    N11 = x11.shape[0]/125 	
    N12 = x12.shape[0]/125 	
    N13 = x13.shape[0]/125 	
    N14 = x14.shape[0]/125 
    retvcwc.cloudsatcwc.RVOD_liq_water_content = numpy.reshape(x11,(125,N11)).astype('Int16')
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_content = numpy.reshape(x12,(125,N12)).astype('Int16')
    retvcwc.cloudsatcwc.RVOD_ice_water_content = numpy.reshape(x13,(125,N13)).astype('Int16')
    retvcwc.cloudsatcwc. IO_RVOD_ice_water_content= numpy.reshape(x14,(125,N14)).astype('Int16')

    x21 = numpy.repeat(cloudsatObjcwc.RVOD_liq_water_content_uncertainty[::,0],idx_matchcwc)
    x22 = numpy.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content_uncertainty[::,0],idx_matchcwc)
    x23 = numpy.repeat(cloudsatObjcwc.RVOD_ice_water_content_uncertainty[::,0],idx_matchcwc)
    x24 = numpy.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content_uncertainty[::,0],idx_matchcwc)
    for i in range(1,125):
        x21 = numpy.concatenate((x21,numpy.repeat(cloudsatObjcwc.RVOD_liq_water_content_uncertainty[::,i],idx_matchcwc))) 
        x22 = numpy.concatenate((x22,numpy.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content_uncertainty[::,i],idx_matchcwc))) 
        x23 = numpy.concatenate((x23,numpy.repeat(cloudsatObjcwc.RVOD_ice_water_content_uncertainty[::,i],idx_matchcwc))) 
        x24 = numpy.concatenate((x24,numpy.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content_uncertainty[::,i],idx_matchcwc))) 

    N21 = x21.shape[0]/125 	
    N22 = x22.shape[0]/125 	
    N23 = x23.shape[0]/125 	
    N24 = x24.shape[0]/125 
    retvcwc.cloudsatcwc.RVOD_liq_water_content_uncertainty = numpy.reshape(x21,(125,N21)).astype('Int8')
    retvcwc.cloudsatcwc.LO_RVOD_liquid_water_content_uncertainty = numpy.reshape(x22,(125,N22)).astype('Int8')
    retvcwc.cloudsatcwc.RVOD_ice_water_content_uncertainty = numpy.reshape(x23,(125,N23)).astype('Int8')
    retvcwc.cloudsatcwc. IO_RVOD_ice_water_content_uncertainty = numpy.reshape(x24,(125,N24)).astype('Int8')

    retvcwc.cloudsatcwc.Temp_min_mixph_K = cloudsatObjcwc.Temp_min_mixph_K
    retvcwc.cloudsatcwc.Temp_max_mixph_K = cloudsatObjcwc.Temp_max_mixph_K
    
    retvcwc.cloudsatcwc.avhrr_linnum = cal_on_avhrrcwc.astype('i')
    retvcwc.cloudsatcwc.avhrr_pixnum = cap_on_avhrrcwc.astype('i')

    retvcwc.avhrr.sec_1970 = numpy.add(avhrrGeoObj.sec1970_start,cal_on_avhrrcwc * DSEC_PER_AVHRR_SCALINE)
    retvcwc.diff_sec_1970 = retvcwc.cloudsatcwc.sec_1970 - retvcwc.avhrr.sec_1970
    min_diff = numpy.minimum.reduce(retvcwc.diff_sec_1970)
    max_diff = numpy.maximum.reduce(retvcwc.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-cloudsat): ",\
          numpy.maximum.reduce(retvcwc.diff_sec_1970),numpy.minimum.reduce(retvcwc.diff_sec_1970)

    print "AVHRR observation time of first cloudsat-avhrr match: ",\
          time.gmtime(retvcwc.avhrr.sec_1970[0])
    print "AVHRR observation time of last cloudsat-avhrr match: ",\
          time.gmtime(retvcwc.avhrr.sec_1970[N-1])

    ll = []
    for i in range(ndimcwc):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsatcwc[i],latCloudsatcwc[i],idx_matchcwc[i])))

    basename = os.path.basename(ctypefile).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_") 
    datapath = "%s/%s/1km/%s/%s/%s" %(DATA_DIR, base_sat, base_year, base_month, AREA1KM)
    if not os.path.exists(datapath):
        os.makedirs(datapath)   
    fd = open("%s/1km_%s_cloudtype_cloudsat-CWC-RVOD_track2.txt"%(datapath,basename),"w")
    fd.writelines(ll)
    fd.close()

    ll = []
    for i in range(N):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_cloudsatcwc[i],lat_cloudsatcwc[i],0)))
    fd = open("%s/1km_%s_cloudtype_cloudsat-CWC-RVOD_track_excl.txt"%(datapath,basename),"w")
    fd.writelines(ll)
    fd.close()

    return retvcwc,min_diff,max_diff
    
#------------------------------------------------------------------------------------------------------------

def reshapeCloudsat1kmCwc(cloudsatfiles,avhrr):
    import time
    import numpy
    import sys
    
    clsat = CloudsatCwcRvodObject
    if avhrr.sec1970_end<avhrr.sec1970_start:
        avhrr_end = int(DSEC_PER_AVHRR_SCALINE5KM*avhrr.num_of_lines+avhrr.sec1970_start)
    else:
        avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    
    startCloudsat = get_cloudsatCwc(cloudsatfiles[0])
    startCloudsat.Profile_time = numpy.add(startCloudsat.Profile_time,startCloudsat.TAI_start)
    
    for i in range(len(cloudsatfiles)-1):
        newCloudsat = get_cloudsatCwc(cloudsatfiles[i+1])
        newCloudsat.Profile_time = numpy.add(newCloudsat.Profile_time,newCloudsat.TAI_start)
        
        clsat_start_all = startCloudsat.sec1970.ravel()
        clsat_new_all = newCloudsat.sec1970.ravel()
        
        if not clsat_start_all[0]<clsat_new_all[0]:
            print "cloudsat files are in the wrong order"
            print("Program cloudsat5km.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9)
        clsat_break = numpy.argmin(numpy.abs(clsat_start_all - clsat_new_all[0]))+1

        startCloudsat.sec1970 = numpy.concatenate((startCloudsat.sec1970[0:clsat_break],newCloudsat.sec1970))
        startCloudsat.longitude = numpy.concatenate((startCloudsat.longitude[0:clsat_break],newCloudsat.longitude))
        startCloudsat.latitude = numpy.concatenate((startCloudsat.latitude[0:clsat_break],newCloudsat.latitude))
        startCloudsat.elevation = numpy.concatenate((startCloudsat.elevation[0:clsat_break],newCloudsat.elevation))
        startCloudsat.Profile_time = numpy.concatenate((startCloudsat.Profile_time[0:clsat_break],newCloudsat.Profile_time))
        startCloudsat.Height = numpy.concatenate((startCloudsat.Height[0:clsat_break,:],newCloudsat.Height))
        startCloudsat.Data_quality = numpy.concatenate((startCloudsat.Data_quality[0:clsat_break],newCloudsat.Data_quality))
        startCloudsat.Data_targetID = numpy.concatenate((startCloudsat.Data_targetID[0:clsat_break],newCloudsat.Data_targetID))
        startCloudsat.RVOD_liq_water_path = numpy.concatenate((startCloudsat.RVOD_liq_water_path[0:clsat_break],newCloudsat.RVOD_liq_water_path))
        startCloudsat.RVOD_liq_water_path_uncertainty = numpy.concatenate((startCloudsat.RVOD_liq_water_path_uncertainty[0:clsat_break],newCloudsat.RVOD_liq_water_path_uncertainty))
        startCloudsat.RVOD_ice_water_path = numpy.concatenate((startCloudsat.RVOD_ice_water_path[0:clsat_break],newCloudsat.RVOD_ice_water_path))
        startCloudsat.RVOD_ice_water_path_uncertainty = numpy.concatenate((startCloudsat.RVOD_ice_water_path_uncertainty[0:clsat_break],newCloudsat.RVOD_ice_water_path_uncertainty))
        startCloudsat.LO_RVOD_liquid_water_path = numpy.concatenate((startCloudsat.LO_RVOD_liquid_water_path[0:clsat_break],newCloudsat.LO_RVOD_liquid_water_path))
        startCloudsat.LO_RVOD_liquid_water_path_uncertainty = numpy.concatenate((startCloudsat.LO_RVOD_liquid_water_path_uncertainty[0:clsat_break],newCloudsat.LO_RVOD_liquid_water_path_uncertainty))
        startCloudsat.IO_RVOD_ice_water_path = numpy.concatenate((startCloudsat.IO_RVOD_ice_water_path[0:clsat_break],newCloudsat.IO_RVOD_ice_water_path))
        startCloudsat.IO_RVOD_ice_water_path_uncertainty = numpy.concatenate((startCloudsat.IO_RVOD_ice_water_path_uncertainty[0:clsat_break],newCloudsat.IO_RVOD_ice_water_path_uncertainty))
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
    clsat.longitude = startCloudsat.longitude[start_break:end_break] 
    clsat.latitude = startCloudsat.latitude[start_break:end_break] 
    clsat.elevation = startCloudsat.elevation[start_break:end_break] 
    clsat.Profile_time = startCloudsat.Profile_time[start_break:end_break]  
    clsat.Data_quality = startCloudsat.Data_quality[start_break:end_break] 
    clsat.Data_targetID = startCloudsat.Data_targetID[start_break:end_break] 
    clsat.Height = startCloudsat.Height[start_break:end_break,:] 
    clsat.RVOD_liq_water_path = startCloudsat.RVOD_liq_water_path[start_break:end_break]
    clsat.RVOD_liq_water_path_uncertainty = startCloudsat.RVOD_liq_water_path_uncertainty[start_break:end_break]
    clsat.RVOD_ice_water_path = startCloudsat.RVOD_ice_water_path[start_break:end_break]
    clsat.RVOD_ice_water_path_uncertainty = startCloudsat.RVOD_ice_water_path_uncertainty[start_break:end_break]
    clsat.LO_RVOD_liquid_water_path = startCloudsat.LO_RVOD_liquid_water_path[start_break:end_break]
    clsat.LO_RVOD_liquid_water_path_uncertainty = startCloudsat.LO_RVOD_liquid_water_path_uncertainty[start_break:end_break]
    clsat.IO_RVOD_ice_water_path = startCloudsat.IO_RVOD_ice_water_path[start_break:end_break]
    clsat.IO_RVOD_ice_water_path_uncertainty = startCloudsat.IO_RVOD_ice_water_path_uncertainty[start_break:end_break]
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
        clsat.TAI_start = clsat.Profile_time[0]
        clsat.Profile_time = clsat.Profile_time - clsat.TAI_start
    
    return clsat
