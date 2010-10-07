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

#MAIN_DIR = "/data/proj/safworks/adam/cloudsat_data"
MAIN_DIR = "/data/proj_nsc1/safworks/kgkarl/ORR-B-datasets/cloudsat/"
#MAIN_DIR = "/data/proj_nsc1/safworks/calipso_cloudsat/data/arctic"
SUB_DIR = "matchups/"
#SUB_DIR = "metop02_calipso_2007spring"
#SUB_DIR = "noaa18_calipso_cloudsat_2007DEC_KG"
#SUB_DIR = "noaa17_calipso_cloudsat_2007JUN_KG"
#SUB_DIR = "noaa18_calipso_cloudsat_2007JUN_KG"

SATPOS_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
EPHE_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
COMPRESS_LVL = 6
COVERAGE_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

#NLINES=1000
NLINES=6000
SWATHWD=2048
NODATA=-9


class CloudsatObject:
    def __init__(self):
        self.longitude=None
        self.latitude=None
        self.avhrr_linnum=None
        self.avhrr_pixnum=None
        self.elevation=None
        self.Profile_time=None
        self.TAI_start=None 
        self.sec_1970=None
        # The data:
        self.CPR_Cloud_mask=None
        self.CPR_Echo_Top=None
        self.Clutter_reduction_flag=None
        self.Data_quality=None
        self.Data_targetID=None
        self.Gaseous_Attenuation=None
        self.MODIS_Cloud_Fraction=None
        self.MODIS_cloud_flag=None
        self.Radar_Reflectivity=None
        self.Height=None
        self.SigmaZero=None
        self.SurfaceHeightBin=None
        self.SurfaceHeightBin_fraction=None
        self.sem_NoiseFloor=None
        self.sem_NoiseFloorVar=None
        self.sem_NoiseGate=None

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

class CloudsatAvhrrTrackObject:
    def __init__(self):
        self.avhrr=ppsAvhrrObject()
        self.cloudsat=CloudsatObject()
        self.cloudsatcwc=CloudsatCwcRvodObject()
        self.diff_sec_1970=None
##########################################################################################################################################   
#pdb.set_trace()

# ----------------------------------------
def readCloudsatAvhrrMatchObj(filename):
    import _pyhl

    retv = CloudsatAvhrrTrackObject()

    a=_pyhl.read_nodelist(filename)
    a.selectAll()
    a.fetch()

    # Match-Up - time difference:
    retv.diff_sec_1970 = a.getNode("/diff_sec_1970").data()

    # Cloudsat:
    retv.cloudsat.longitude = a.getNode("/cloudsat/longitude").data()
    retv.cloudsat.latitude = a.getNode("/cloudsat/latitude").data()
    retv.cloudsat.avhrr_linnum = a.getNode("/cloudsat/avhrr_linnum").data()
    retv.cloudsat.avhrr_pixnum = a.getNode("/cloudsat/avhrr_pixnum").data()
    
    retv.cloudsat.cloud_mask = a.getNode("/cloudsat/cloud_mask").data()
    retv.cloudsat.Radar_Reflectivity = a.getNode("/cloudsat/Radar_Reflectivity").data()
    retv.cloudsat.Height = a.getNode("/cloudsat/Height").data()
    retv.cloudsat.echo_top = a.getNode("/cloudsat/echo_top").data()
    retv.cloudsat.SurfaceHeightBin = a.getNode("/cloudsat/SurfaceHeightBin").data()
    retv.cloudsat.SurfaceHeightBin_fraction = a.getNode("/cloudsat/SurfaceHeightBin_fraction").data()
    
    retv.cloudsat.elevation = a.getNode("/cloudsat/elevation").data()
    retv.cloudsat.sec_1970 = a.getNode("/cloudsat/sec_1970").data()
    retv.cloudsat.MODIS_Cloud_Fraction = a.getNode("/cloudsat/MODIS_Cloud_Fraction").data()
    retv.cloudsat.MODIS_cloud_flag = a.getNode("/cloudsat/MODIS_cloud_flag").data()

    # AVHRR:
    retv.avhrr.longitude = a.getNode("/avhrr/longitude").data()
    retv.avhrr.latitude = a.getNode("/avhrr/latitude").data()
    retv.avhrr.sec_1970 = a.getNode("/avhrr/sec_1970").data()
    retv.avhrr.cloudtype = a.getNode("/avhrr/cloudtype").data()
    retv.avhrr.ctth_height = a.getNode("/avhrr/ctth_height").data()
    retv.avhrr.ctth_pressure = a.getNode("/avhrr/ctth_pressure").data()
    retv.avhrr.ctth_temperature = a.getNode("/avhrr/ctth_temperature").data()
    retv.avhrr.bt11micron = a.getNode("/avhrr/bt11micron").data()
    retv.avhrr.bt12micron = a.getNode("/avhrr/bt12micron").data()
    retv.avhrr.surftemp = a.getNode("/avhrr/surftemp").data()
    
    # CWC-RVOD:
    # ====
# Geolocation
    retv.cloudsatcwc.longitude = a.getNode("/cloudsatcwc/longitude").data()
    retv.cloudsatcwc.latitude = a.getNode("/cloudsatcwc/latitude").data()
    retv.cloudsatcwc.elevation = a.getNode("/cloudsatcwc/elevation").data()
    retv.cloudsatcwc.Height = a.getNode("/cloudsatcwc/Height").data()
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    retv.cloudsatcwc.Profile_time = a.getNode("/cloudsatcwc/Profile_time").data()
    retv.cloudsatcwc.TAI_start = a.getNode("/cloudsatcwc/TAI_start").data()
    retv.cloudsatcwc.sec_1970 = a.getNode("/cloudsatcwc/sec_1970").data()
# The data
    retv.cloudsatcwc.Data_quality = a.getNode("/cloudsatcwc/Data_quality").data()
    retv.cloudsatcwc.Data_targetID = a.getNode("/cloudsatcwc/Data_targetID").data()
    retv.cloudsatcwc.RVOD_liq_water_path = a.getNode("/cloudsatcwc/RVOD_liq_water_path").data()
    retv.cloudsatcwc.RVOD_liq_water_path_uncertainty = a.getNode("/cloudsatcwc/RVOD_liq_water_path_uncertainty").data()
    retv.cloudsatcwc.RVOD_ice_water_path = a.getNode("/cloudsatcwc/RVOD_ice_water_path").data()
    retv.cloudsatcwc.RVOD_ice_water_path_uncertainty = a.getNode("/cloudsatcwc/RVOD_ice_water_path_uncertainty").data()
    retv.cloudsatcwc.LO_RVOD_liquid_water_path = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_path").data()
    retv.cloudsatcwc.LO_RVOD_liquid_water_path_uncertainty = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_path_uncertainty").data()
    retv.cloudsatcwc.IO_RVOD_ice_water_path = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_path").data()
    retv.cloudsatcwc.IO_RVOD_ice_water_path_uncertainty = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_path_uncertainty").data()

   
    retv.cloudsatcwc.RVOD_liq_water_content = a.getNode("/cloudsatcwc/RVOD_liq_water_content").data()
    retv.cloudsatcwc.RVOD_liq_water_content_uncertainty = a.getNode("/cloudsatcwc/RVOD_liq_water_content_uncertainty").data()
    retv.cloudsatcwc.RVOD_ice_water_content = a.getNode("/cloudsatcwc/RVOD_ice_water_content").data()
    retv.cloudsatcwc.RVOD_ice_water_content_uncertainty = a.getNode("/cloudsatcwc/RVOD_ice_water_content_uncertainty").data()
    retv.cloudsatcwc.LO_RVOD_liquid_water_content = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_content").data()
    retv.cloudsatcwc.LO_RVOD_liquid_water_content_uncertainty = a.getNode("/cloudsatcwc/LO_RVOD_liquid_water_content_uncertainty").data()
    retv.cloudsatcwc.IO_RVOD_ice_water_content = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_content").data()
    retv.cloudsatcwc.IO_RVOD_ice_water_content_uncertainty = a.getNode("/cloudsatcwc/IO_RVOD_ice_water_content_uncertainty").data()
    retv.cloudsatcwc.Temp_min_mixph_K = a.getNode("/cloudsatcwc/Temp_min_mixph_K").data()
    retv.cloudsatcwc.Temp_max_mixph_K = a.getNode("/cloudsatcwc/Temp_max_mixph_K").data()

##########################################################################################################################################
    #pdb.set_trace()

    return retv


# ----------------------------------------
def writeCloudsatAvhrrMatchObj(filename,ca_obj,compress_lvl):
    import _pyhl
    status = -1
    
    a=_pyhl.nodelist()

    shape = [ca_obj.cloudsat.longitude.shape[0]]

    # Match-Up - time difference:
    # ====
    b=_pyhl.node(_pyhl.DATASET_ID,"/diff_sec_1970")
    b.setArrayValue(1,shape,ca_obj.diff_sec_1970,"double",-1)
    a.addNode(b)

    # Cloudsat
    # ====
    b=_pyhl.node(_pyhl.GROUP_ID,"/cloudsat")
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/longitude")
    b.setArrayValue(1,shape,ca_obj.cloudsat.longitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/latitude")
    b.setArrayValue(1,shape,ca_obj.cloudsat.latitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/avhrr_linnum")
    b.setArrayValue(1,shape,ca_obj.cloudsat.avhrr_linnum,"int",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/avhrr_pixnum")
    b.setArrayValue(1,shape,ca_obj.cloudsat.avhrr_pixnum,"int",-1)
    a.addNode(b)

    shape2d = ca_obj.cloudsat.cloud_mask.shape
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/cloud_mask")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat.cloud_mask,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/Radar_Reflectivity")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat.Radar_Reflectivity,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/Height")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat.Height,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/echo_top")
    b.setArrayValue(1,shape,ca_obj.cloudsat.echo_top,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/SurfaceHeightBin")
    b.setArrayValue(1,shape,ca_obj.cloudsat.SurfaceHeightBin,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/SurfaceHeightBin_fraction")
    b.setArrayValue(1,shape,ca_obj.cloudsat.SurfaceHeightBin_fraction,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/elevation")
    b.setArrayValue(1,shape,ca_obj.cloudsat.elevation,"short",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/sec_1970")
    b.setArrayValue(1,shape,ca_obj.cloudsat.sec_1970,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/MODIS_Cloud_Fraction")
    b.setArrayValue(1,shape,ca_obj.cloudsat.MODIS_Cloud_Fraction,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/MODIS_cloud_flag")
    b.setArrayValue(1,shape,ca_obj.cloudsat.MODIS_cloud_flag,"uchar",-1)
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

    # CWC-RVOD
    # ====
    shapecwc = [ca_obj.cloudsatcwc.longitude.shape[0]]
    shape2dcwc = ca_obj.cloudsatcwc.Height.shape
    shapeTAIcwc = [ca_obj.cloudsatcwc.TAI_start.shape[0]]
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
# International Atomic Time (TAI) seconds from Jan 1, 1993:
##########################################################################################################################################
    #pdb.set_trace()
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/Profile_time")
    b.setArrayValue(1,shapecwclong,ca_obj.cloudsatcwc.Profile_time,"double",-1)
    a.addNode(b) 
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsatcwc/TAI_start")
    b.setArrayValue(1,shapeTAIcwc,ca_obj.cloudsatcwc.TAI_start,"double",-1)
    a.addNode(b)
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
##########################################################################################################################################
    #pdb.set_trace()
    return status

# -----------------------------------------------------
def select_cloudsat_inside_avhrr(cloudsatObj,cal,sec1970_start_end,sec_timeThr):
    # Cloudsat times are already in UTC in sec since 1970
    import numpy.oldnumeric as Numeric

    sec1970_start,sec1970_end = sec1970_start_end
    
    # Select the points inside the avhrr swath:
    # Allowing for sec_timeThr seconds deviation:
    idx_okay = Numeric.logical_and(Numeric.greater(\
        cloudsatObj.sec1970,sec1970_start - sec_timeThr),
                                   Numeric.less(\
        cloudsatObj.sec1970,sec1970_end + sec_timeThr))
    
    idx_match = Numeric.not_equal(cal,NODATA)
    
    return idx_match

# -----------------------------------------------------
def get_cloudsat(filename,filenamecwc):
    import _pypps_filters
    import numpy.oldnumeric as Numeric
    import time

    # Read CLOUDSAT Radar data:
    cloudsat,cloudsatcwc = read_cloudsat(filename,filenamecwc)

    # GEOPROF
    # ====
    lonCloudsat = cloudsat.longitude.ravel()
    latCloudsat = cloudsat.latitude.ravel()
    ndim = lonCloudsat.shape[0]

    # --------------------------------------------------------------------
    # Time in seconds since 1970:
    
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone

    print "GEOPROF TAI start time: ",cloudsat.TAI_start
    start_sec1970 = cloudsat.TAI_start + dsec
    start_time = time.gmtime(start_sec1970[0])
    end_sec1970 = cloudsat.TAI_start + dsec + cloudsat.Profile_time[ndim-1]
    end_time = time.gmtime(end_sec1970[0])
    print "GEOPROF Start and end times: ",start_time,end_time
    cloudsat.sec1970 = cloudsat.Profile_time + start_sec1970

    # CWC-RVOD
    # ====

    lonCloudsatcwc = cloudsatcwc.longitude.ravel()
    latCloudsatcwc = cloudsatcwc.latitude.ravel()
    ndimcwc = lonCloudsatcwc.shape[0]

    # --------------------------------------------------------------------
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

    return cloudsat,cloudsatcwc

# -----------------------------------------------------
def read_cloudsat(filename,filenamecwc):
    import _pyhl
    import numpy.oldnumeric as Numeric
    
    # GEOPROF
    # ====
    a=_pyhl.read_nodelist(filename)
     
    b=a.getNodeNames()
    a.selectAll()
    a.fetch()

    retv = CloudsatObject()
    root="/2B-GEOPROF"

    d=a.getNode("%s/Geolocation Fields/Longitude"%root)
    retv.longitude=Numeric.fromstring(d.data(),'f').astype('d')
    d=a.getNode("%s/Geolocation Fields/Latitude"%root)
    retv.latitude=Numeric.fromstring(d.data(),'f').astype('d')
    
    # International Atomic Time (TAI) seconds from Jan 1, 1993:
    d=a.getNode("%s/Geolocation Fields/Profile_time"%root)
    retv.Profile_time=Numeric.fromstring(d.data(),'f').astype('d')
    d=a.getNode("%s/Geolocation Fields/TAI_start"%root)
    retv.TAI_start=Numeric.fromstring(d.data(),'d').astype('d')

    d=a.getNode("%s/Geolocation Fields/DEM_elevation"%root)
    retv.elevation=Numeric.fromstring(d.data(),'Int16').astype('Int16')
   
    # The radar data:
    # 2-d arrays of dimension Nx125:
    d=a.getNode("%s/Data Fields/CPR_Cloud_mask"%root)
    retv.CPR_Cloud_mask=d.data()

##########################################################################################################################################
    #pdb.set_trace()
    d=a.getNode("%s/Geolocation Fields/Height"%root)
    retv.Height=d.data()

    d=a.getNode("%s/Data Fields/Radar_Reflectivity"%root)
    retv.Radar_Reflectivity=d.data().astype('Int16')
    d=a.getNode("%s/Data Fields/Gaseous_Attenuation"%root)
    retv.Gaseous_Attenuation=d.data().astype('Int16')

    d=a.getNode("%s/Data Fields/CPR_Echo_Top"%root)
    retv.CPR_Echo_Top=Numeric.fromstring(d.data(),'b').astype('b')
    d=a.getNode("%s/Data Fields/Clutter_reduction_flag"%root)
    retv.Clutter_reduction_flag=Numeric.fromstring(d.data(),'b').astype('b')
    d=a.getNode("%s/Data Fields/Data_quality"%root)
    retv.Data_quality=Numeric.fromstring(d.data(),'b').astype('b')
    d=a.getNode("%s/Data Fields/Data_targetID"%root)
    retv.Data_targetID=Numeric.fromstring(d.data(),'b').astype('b')
    d=a.getNode("%s/Data Fields/MODIS_Cloud_Fraction"%root)
    retv.MODIS_Cloud_Fraction=Numeric.fromstring(d.data(),'b').astype('b')
    d=a.getNode("%s/Data Fields/MODIS_cloud_flag"%root)
    retv.MODIS_cloud_flag=Numeric.fromstring(d.data(),'b').astype('b')
    d=a.getNode("%s/Data Fields/Sigma-Zero"%root)
    retv.SigmaZero=Numeric.fromstring(d.data(),'Int16').astype('Int16')
    d=a.getNode("%s/Data Fields/SurfaceHeightBin"%root)
    retv.SurfaceHeightBin=Numeric.fromstring(d.data(),'b').astype('b')
    d=a.getNode("%s/Data Fields/SurfaceHeightBin_fraction"%root)
    retv.SurfaceHeightBin_fraction=Numeric.fromstring(d.data(),'f').astype('d')
    d=a.getNode("%s/Data Fields/sem_NoiseFloor"%root)
    retv.sem_NoiseFloor=Numeric.fromstring(d.data(),'f').astype('d')
    d=a.getNode("%s/Data Fields/sem_NoiseFloorVar"%root)
    retv.sem_NoiseFloorVar=Numeric.fromstring(d.data(),'f').astype('d')
    d=a.getNode("%s/Data Fields/sem_NoiseGate"%root)
    retv.sem_NoiseGate=Numeric.fromstring(d.data(),'b').astype('b')
    
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
    retvcwc.longitude=Numeric.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Geolocation Fields/Latitude"%rootcwc)
    retvcwc.latitude=Numeric.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Geolocation Fields/DEM_elevation"%rootcwc)
    retvcwc.elevation=Numeric.fromstring(ccwc.data(),'Int16').astype('Int16')
    ccwc=acwc.getNode("%s/Geolocation Fields/Height"%rootcwc)
    retvcwc.Height=ccwc.data()
# International Atomic Time (TAI) seconds from Jan 1, 1993:
    ccwc=acwc.getNode("%s/Geolocation Fields/Profile_time"%rootcwc)
    retvcwc.Profile_time=Numeric.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Geolocation Fields/TAI_start"%rootcwc)
    retvcwc.TAI_start=Numeric.fromstring(ccwc.data(),'d').astype('d')
# The data
    ccwc=acwc.getNode("%s/Data Fields/Data_quality"%rootcwc)
    retvcwc.Data_quality=Numeric.fromstring(ccwc.data(),'b').astype('b')
    ccwc=acwc.getNode("%s/Data Fields/Data_targetID"%rootcwc)
    retvcwc.Data_targetID=Numeric.fromstring(ccwc.data(),'b').astype('b')
# IWP and LWP		
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_path"%rootcwc)
    retvcwc.RVOD_liq_water_path=Numeric.fromstring(ccwc.data(),'Int16').astype('Int16')			
    ccwc=acwc.getNode("%s/Data Fields/RVOD_liq_water_path_uncertainty"%rootcwc)
    retvcwc.RVOD_liq_water_path_uncertainty=Numeric.fromstring(ccwc.data(),'b').astype('b')
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_path"%rootcwc)
    retvcwc.RVOD_ice_water_path=Numeric.fromstring(ccwc.data(),'Int16').astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/RVOD_ice_water_path_uncertainty"%rootcwc)
    retvcwc.RVOD_ice_water_path_uncertainty=Numeric.fromstring(ccwc.data(),'b').astype('b')
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_path"%rootcwc)
    retvcwc.LO_RVOD_liquid_water_path=Numeric.fromstring(ccwc.data(),'Int16').astype('Int16')			
    ccwc=acwc.getNode("%s/Data Fields/LO_RVOD_liquid_water_path_uncertainty"%rootcwc)
    retvcwc.LO_RVOD_liquid_water_path_uncertainty=Numeric.fromstring(ccwc.data(),'b').astype('b')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_path"%rootcwc)
    retvcwc.IO_RVOD_ice_water_path=Numeric.fromstring(ccwc.data(),'Int16').astype('Int16')
    ccwc=acwc.getNode("%s/Data Fields/IO_RVOD_ice_water_path_uncertainty"%rootcwc)
    retvcwc.IO_RVOD_ice_water_path_uncertainty=Numeric.fromstring(ccwc.data(),'b').astype('b')
# IWC and LWC
##########################################################################################################################################
    #pdb.set_trace()
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
    retvcwc.Temp_min_mixph_K=Numeric.fromstring(ccwc.data(),'f').astype('d')
    ccwc=acwc.getNode("%s/Data Fields/Temp_max_mixph_K"%rootcwc)
    retvcwc.Temp_max_mixph_K=Numeric.fromstring(ccwc.data(),'f').astype('d')
##########################################################################################################################################
    #pdb.set_trace()
    return retv,retvcwc

#retvcwc.RVOD_liq_water_content.shape


# --------------------------------------------
def get_cloudsat_avhrr_linpix(avhrrIn,lon,lat):
    import numpy.oldnumeric as Numeric

    tmppcs="tmpproj"
    define_pcs(tmppcs, "Plate Caree, central meridian at 15E",
               ['proj=eqc','ellps=bessel', 'lon_0=15'])

    startline=0

    Inside=0
    HasEncounteredMatch=0
    i=0
    while startline < avhrrIn.longitude.shape[0]:
        
        write_log("INFO","Calling get_cloudsat_avhrr_linpix: start-line = ",startline)
        tmpaid = "tmparea_%d"%i
        endline = startline+NLINES
        coverage_filename = "%s/coverage_avhrr_cloudsat_matchup_%s_%.5d_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrrIn.satellite,avhrrIn.orbit,
                             startline,endline,tmpaid)
        write_log("INFO","Coverage filename = ",coverage_filename)
        cal,cap,ok = get_cloudsat_avhrr_linpix_segment(avhrrIn,lon,lat,
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
            cloudsat_avhrr_line,cloudsat_avhrr_pixel = Numeric.array(cal),Numeric.array(cap)
        else:
            # Merge:
            cloudsat_avhrr_line = Numeric.where(Numeric.equal(cloudsat_avhrr_line,-9),cal,cloudsat_avhrr_line)
            cloudsat_avhrr_pixel = Numeric.where(Numeric.equal(cloudsat_avhrr_pixel,-9),cap,cloudsat_avhrr_pixel)

        startline=startline+NLINES
        i=i+1
        ##########################################################################################################################################
    #pdb.set_trace()
    return cloudsat_avhrr_line,cloudsat_avhrr_pixel

# --------------------------------------------
def get_cloudsat_avhrr_linpix_segment(avhrrIn,lon,lat,lines,swath_width,tmppcs,
                                      tmpaid,covfilename):
    import numpy.oldnumeric as Numeric
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
    idx = Numeric.arange(idx_start,idx_end)
    
    linearr = Numeric.divide(idx,swath_width)
    write_log("INFO","Start and end line numbers: ",linearr[0],linearr[idx.shape[0]-1])
    
    linearr = Numeric.reshape(linearr,(nlines,swath_width))
    pixelarr = Numeric.fmod(idx,swath_width).astype('l')
    pixelarr = Numeric.reshape(pixelarr,(nlines,swath_width))

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
    areaObj = area.area("arctic_super_5010")

    
    if not os.path.exists(covfilename):
        write_log("INFO","Create Coverage map...")
        cov = _satproj.create_coverage(areaObj,lonarr,latarr,1)
        print covfilename
        writeCoverage(cov,covfilename,"satproj","arctic_super_5010")
    else:
        write_log("INFO","Read the AVHRR-CLOUDSAT matchup coverage from file...")
        cov,info = readCoverage(covfilename)
    
    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA)
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA)


    write_log("INFO","Go through cloudsat track:")
    cloudsat_avhrr_line = []
    cloudsat_avhrr_pixel = []
    for i in range(ndim):
        xy_tup=pps_gisdata.lonlat2xy("arctic_super_5010",lon[i],lat[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
##        dimx=4500 Should be 5010!!!/KG
##        dimy=4500
        dimx=5010
        dimy=5010
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            cloudsat_avhrr_line.append(mapped_line[y,x])
            cloudsat_avhrr_pixel.append(mapped_pixel[y,x])
        else:
            cloudsat_avhrr_line.append(-9)
            cloudsat_avhrr_pixel.append(-9)

    cloudsat_avhrr_line = Numeric.array(cloudsat_avhrr_line)
    cloudsat_avhrr_pixel = Numeric.array(cloudsat_avhrr_pixel)

    x=Numeric.repeat(cloudsat_avhrr_line,Numeric.not_equal(cloudsat_avhrr_line,-9))
    write_log("INFO","Number of matching points = ",x.shape[0])
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0

    return cloudsat_avhrr_line,cloudsat_avhrr_pixel,matchOk

# -----------------------------------------------------
def match_cloudsat_avhrr(ctypefile,cloudsatObj,cloudsatObjcwc,avhrrGeoObj,avhrrObj,ctype,ctth,surft):
    import numpy.oldnumeric as Numeric
    import time
    import string
    
    retv = CloudsatAvhrrTrackObject()
    
    lonCloudsat = cloudsatObj.longitude.ravel()
    latCloudsat = cloudsatObj.latitude.ravel()
    ndim = lonCloudsat.shape[0]

    lonCloudsatcwc = cloudsatObjcwc.longitude.ravel()		# Cwc
    latCloudsatcwc = cloudsatObjcwc.latitude.ravel()		# Cwc
    ndimcwc = lonCloudsat.shape[0]				# Cwc

    # --------------------------------------------------------------------

    cal,cap = get_cloudsat_avhrr_linpix(avhrrGeoObj,lonCloudsat,latCloudsat)
    calcwc,capcwc=get_cloudsat_avhrr_linpix(avhrrGeoObj,lonCloudsatcwc,latCloudsatcwc)	# Cwc

    #print len(cal), len(cap), ndim

    idx_cloudsat = Numeric.arange(ndim)
    idx_cloudsat_on_avhrr=Numeric.repeat(idx_cloudsat,Numeric.not_equal(cal,NODATA))    
    lon_cloudsat = Numeric.repeat(lonCloudsat,Numeric.not_equal(cal,NODATA))
    lat_cloudsat = Numeric.repeat(latCloudsat,Numeric.not_equal(cal,NODATA))

    idx_cloudsatcwc = Numeric.arange(ndimcwc)							# Cwc
    idx_cloudsatcwc_on_avhrr=Numeric.repeat(idx_cloudsatcwc,Numeric.not_equal(calcwc,NODATA)) 	# Cwc   
    lon_cloudsatcwc = Numeric.repeat(lonCloudsatcwc,Numeric.not_equal(calcwc,NODATA))		# Cwc
    lat_cloudsatcwc = Numeric.repeat(latCloudsatcwc,Numeric.not_equal(calcwc,NODATA))		# Cwc
    # Cloudsat line,pixel inside AVHRR swath:
    cal_on_avhrr = Numeric.repeat(cal,Numeric.not_equal(cal,NODATA))
    cap_on_avhrr = Numeric.repeat(cap,Numeric.not_equal(cal,NODATA))
    timetup_start = time.gmtime(avhrrGeoObj.sec1970_start)
    timetup_end   = time.gmtime(avhrrGeoObj.sec1970_end)

    cal_on_avhrrcwc = Numeric.repeat(calcwc,Numeric.not_equal(calcwc,NODATA))	# Cwc
    cap_on_avhrrcwc = Numeric.repeat(capcwc,Numeric.not_equal(calcwc,NODATA))	# Cwc
    timetup_startcwc = time.gmtime(avhrrGeoObj.sec1970_start)		# Cwc
    timetup_endcwc   = time.gmtime(avhrrGeoObj.sec1970_end)		# Cwc

    sec_timeThr = 60*20 # Allow for 20 minute deviation between AVHRR and CLOUDSAT matchup
    secTup = avhrrGeoObj.sec1970_start,avhrrGeoObj.sec1970_end
    idx_match = select_cloudsat_inside_avhrr(cloudsatObj,cal,secTup,sec_timeThr)

    sec_timeThrcwc = 60*20 # Allow for 20 minute deviation between AVHRR and CLOUDSAT Cwcmatchup
    secTupcwc = avhrrGeoObj.sec1970_start,avhrrGeoObj.sec1970_end
    idx_matchcwc = select_cloudsat_inside_avhrr(cloudsatObjcwc,calcwc,secTupcwc,sec_timeThrcwc)	########kanske blir fel??????????

    print "Make CPR echo top array..."
    retv.cloudsat.echo_top = Numeric.repeat(cloudsatObj.CPR_Echo_Top,idx_match).astype('b')

    retv.cloudsat.sec_1970 = Numeric.repeat(cloudsatObj.sec1970,idx_match)

    retv.cloudsat.latitude = Numeric.repeat(latCloudsat,idx_match)
    retv.cloudsat.longitude = Numeric.repeat(lonCloudsat,idx_match)

    retv.cloudsat.MODIS_Cloud_Fraction = Numeric.repeat(cloudsatObj.MODIS_Cloud_Fraction,idx_match)
    retv.cloudsat.MODIS_cloud_flag = Numeric.repeat(cloudsatObj.MODIS_cloud_flag,idx_match)
    
    print "cap_on_avhrr.shape: ",cap_on_avhrr.shape
    retv.cloudsat.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.cloudsat.avhrr_pixnum = cap_on_avhrr.astype('i')
    
    print "Make CPR cloud mask array..."
    x = Numeric.repeat(cloudsatObj.CPR_Cloud_mask[::,0],idx_match)
    for i in range(1,125):
        x = Numeric.concatenate((x,Numeric.repeat(cloudsatObj.CPR_Cloud_mask[::,i],idx_match)))
    N = x.shape[0]/125
    retv.cloudsat.cloud_mask = Numeric.reshape(x,(125,N)).astype('b')

    missing_data = -9.9
    print "Make Radar reflectivity array..."
    x = Numeric.repeat(cloudsatObj.Radar_Reflectivity[::,0],idx_match)
    for i in range(1,125):
        x = Numeric.concatenate((x,Numeric.repeat(cloudsatObj.Radar_Reflectivity[::,i],idx_match)))
    N = x.shape[0]/125
    RadarRefl = Numeric.reshape(x,(125,N))
    RadarRefl = Numeric.where(Numeric.less(RadarRefl,0),missing_data,RadarRefl)
    retv.cloudsat.Radar_Reflectivity = RadarRefl.astype('Int16')

    x = Numeric.repeat(cloudsatObj.Height[::,0],idx_match)
    for i in range(1,125):
        x = Numeric.concatenate((x,Numeric.repeat(cloudsatObj.Height[::,i],idx_match)))
    N = x.shape[0]/125
    retv.cloudsat.Height = Numeric.reshape(x,(125,N)).astype('Int16')

    # One-d arrays:
    retv.cloudsat.SurfaceHeightBin = \
                       Numeric.repeat(cloudsatObj.SurfaceHeightBin.ravel(),
                                      idx_match.ravel()).astype('b')
    retv.cloudsat.SurfaceHeightBin_fraction = \
                       Numeric.repeat(cloudsatObj.SurfaceHeightBin_fraction.ravel(),
                                      idx_match.ravel()).astype('d')

    print "Cloudsat observation time of first cloudsat-avhrr match: ",\
          time.gmtime(retv.cloudsat.sec_1970[0])
    print "Cloudsat observation time of last cloudsat-avhrr match: ",\
          time.gmtime(retv.cloudsat.sec_1970[N-1])
    
    # Elevation is given in km's. Convert to meters:
    retv.cloudsat.elevation = Numeric.repeat(\
        cloudsatObj.elevation.ravel(),idx_match.ravel()).astype('Int16')

    retv.avhrr.sec_1970 = Numeric.add(avhrrGeoObj.sec1970_start,
                                      cal_on_avhrr * DSEC_PER_AVHRR_SCALINE)
    retv.diff_sec_1970 = retv.cloudsat.sec_1970 - retv.avhrr.sec_1970
    min_diff = Numeric.minimum.reduce(retv.diff_sec_1970)
    max_diff = Numeric.maximum.reduce(retv.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-cloudsat): ",\
          Numeric.maximum.reduce(retv.diff_sec_1970),Numeric.minimum.reduce(retv.diff_sec_1970)

    print "AVHRR observation time of first cloudsat-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[0])
    print "AVHRR observation time of last cloudsat-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[N-1])

    # Make the latitude and pps cloudtype on the cloudsat track:
    # line and pixel arrays have equal dimensions
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
    for i in range(cal_on_avhrr.shape[0]):
        lat_avhrr_track.append(avhrrGeoObj.latitude[cal_on_avhrr[i],cap_on_avhrr[i]])
        lon_avhrr_track.append(avhrrGeoObj.longitude[cal_on_avhrr[i],cap_on_avhrr[i]])
        ctype_track.append(ctype.cloudtype[cal_on_avhrr[i],cap_on_avhrr[i]])
        surft_track.append(surft[cal_on_avhrr[i],cap_on_avhrr[i]])
        if avhrrObj.channels[3].data[cal_on_avhrr[i],cap_on_avhrr[i]] == avhrrObj.nodata:
            b11 = -9.
        else:
            b11 = avhrrObj.channels[3].data[cal_on_avhrr[i],cap_on_avhrr[i]] * avhrrObj.channels[3].gain + avhrrObj.channels[3].intercept
        bt11micron_track.append(b11)
        if avhrrObj.channels[4].data[cal_on_avhrr[i],cap_on_avhrr[i]] == avhrrObj.nodata:
            b12 = -9.
        else:
            b12 = avhrrObj.channels[4].data[cal_on_avhrr[i],cap_on_avhrr[i]] * avhrrObj.channels[4].gain + avhrrObj.channels[4].intercept
        bt12micron_track.append(b12)
        if ctth == None:
            continue
        if ctth.height[cal_on_avhrr[i],cap_on_avhrr[i]] == ctth.h_nodata:
            hh = -9.
        else:
            hh = ctth.height[cal_on_avhrr[i],cap_on_avhrr[i]] * ctth.h_gain + ctth.h_intercept
        ctth_height_track.append(hh)
        if ctth.temperature[cal_on_avhrr[i],cap_on_avhrr[i]] == ctth.t_nodata:
            tt = -9.
        else:
            tt = ctth.temperature[cal_on_avhrr[i],cap_on_avhrr[i]] * ctth.t_gain + \
                 ctth.t_intercept
        ctth_temperature_track.append(tt)
        if ctth.pressure[cal_on_avhrr[i],cap_on_avhrr[i]] == ctth.p_nodata:
            pp = -9.
        else:
            pp = ctth.pressure[cal_on_avhrr[i],cap_on_avhrr[i]] * ctth.p_gain + ctth.p_intercept
        ctth_pressure_track.append(pp)

    retv.avhrr.latitude = Numeric.array(lat_avhrr_track)
    retv.avhrr.longitude = Numeric.array(lon_avhrr_track)
    retv.avhrr.cloudtype = Numeric.array(ctype_track)
    retv.avhrr.bt11micron = Numeric.array(bt11micron_track)
    retv.avhrr.bt12micron = Numeric.array(bt12micron_track)
    if ctth:
        retv.avhrr.ctth_height = Numeric.array(ctth_height_track)
        retv.avhrr.ctth_pressure = Numeric.array(ctth_pressure_track)
        retv.avhrr.ctth_temperature = Numeric.array(ctth_temperature_track)
    retv.avhrr.surftemp = Numeric.array(surft_track)

    print "AVHRR-PPS Cloud Type,latitude: shapes = ",\
          retv.avhrr.cloudtype.shape,retv.avhrr.latitude.shape

    ll = []
    for i in range(ndim):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],idx_match[i])))
    basename = os.path.basename(ctypefile).split(".h5")[0]
    basename = string.join(basename.split("_")[0:4],"_")    
    fd = open("data/%s_cloudtype_cloudsat_track2.txt"%(basename),"w")
    fd.writelines(ll)
    fd.close()

    ll = []
    for i in range(N):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_cloudsat[i],lat_cloudsat[i],0)))
    fd = open("data/%s_cloudtype_cloudsat_track_excl.txt"%(basename),"w")
    fd.writelines(ll)
    fd.close()
    
# CWC-RVOD
    # ====
    retv.cloudsatcwc.TAI_start=cloudsatObjcwc.TAI_start
    retv.cloudsatcwc.Profile_time=cloudsatObjcwc.Profile_time
    retv.cloudsatcwc.sec_1970=Numeric.repeat(cloudsatObjcwc.sec1970,idx_matchcwc)
    retv.cloudsatcwc.latitude=Numeric.repeat(cloudsatObjcwc.latitude,idx_matchcwc)
    retv.cloudsatcwc.longitude=Numeric.repeat(cloudsatObjcwc.longitude,idx_matchcwc)

    x = Numeric.repeat(cloudsatObjcwc.Height[::,0],idx_matchcwc) # copies the first column
    for i in range(1,125):
        x = Numeric.concatenate((x,Numeric.repeat(cloudsatObj.Height[::,i],idx_matchcwc))) # Adds all columns to one long column
    N = x.shape[0]/125 		# Finds how many columns there should be
    retv.cloudsatcwc.Height = Numeric.reshape(x,(125,N)).astype('Int16')	# Reshape the long x to a 125xN matrix 

    retv.cloudsatcwc.elevation = Numeric.repeat(\
        cloudsatObjcwc.elevation.ravel(),idx_matchcwc.ravel()).astype('Int16')
    retv.cloudsatcwc.Data_quality=cloudsatObjcwc.Data_quality
    retv.cloudsatcwc.Data_targetID=cloudsatObjcwc.Data_targetID
    
    retv.cloudsatcwc.RVOD_liq_water_path = Numeric.repeat(cloudsatObjcwc.RVOD_liq_water_path,idx_matchcwc).astype('Int16')
    retv.cloudsatcwc.RVOD_liq_water_path_uncertainty = Numeric.repeat(\
        cloudsatObjcwc.RVOD_liq_water_path_uncertainty,idx_matchcwc).astype('Int8')
    retv.cloudsatcwc.RVOD_ice_water_path = Numeric.repeat(cloudsatObjcwc.RVOD_ice_water_path,idx_matchcwc).astype('Int16')
    retv.cloudsatcwc.RVOD_ice_water_path_uncertainty = Numeric.repeat(\
        cloudsatObjcwc.RVOD_ice_water_path_uncertainty,idx_matchcwc).astype('Int8')
    retv.cloudsatcwc.LO_RVOD_liquid_water_path = Numeric.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_path,idx_matchcwc).astype('Int16')
    retv.cloudsatcwc.LO_RVOD_liquid_water_path_uncertainty = Numeric.repeat(\
        cloudsatObjcwc.LO_RVOD_liquid_water_path_uncertainty,idx_matchcwc).astype('Int8')
    retv.cloudsatcwc.IO_RVOD_ice_water_path = Numeric.repeat(cloudsatObjcwc.IO_RVOD_ice_water_path,idx_matchcwc).astype('Int16')
    retv.cloudsatcwc.IO_RVOD_ice_water_path_uncertainty = Numeric.repeat(\
        cloudsatObjcwc.IO_RVOD_ice_water_path_uncertainty,idx_matchcwc).astype('Int8')
##########################################################################################################################################
    #pdb.set_trace()
    x11 = Numeric.repeat(cloudsatObjcwc.RVOD_liq_water_content[::,0],idx_matchcwc)
    x12 = Numeric.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content[::,0],idx_matchcwc)
    x13 = Numeric.repeat(cloudsatObjcwc.RVOD_ice_water_content[::,0],idx_matchcwc)
    x14 = Numeric.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content[::,0],idx_matchcwc)
    for i in range(1,125):
        x11 = Numeric.concatenate((x11,Numeric.repeat(cloudsatObjcwc.RVOD_liq_water_content[::,i],idx_matchcwc))) 
        x12 = Numeric.concatenate((x12,Numeric.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content[::,i],idx_matchcwc))) 
        x13 = Numeric.concatenate((x13,Numeric.repeat(cloudsatObjcwc.RVOD_ice_water_content[::,i],idx_matchcwc))) 
        x14 = Numeric.concatenate((x14,Numeric.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content[::,i],idx_matchcwc))) 

    N11 = x11.shape[0]/125 	
    N12 = x12.shape[0]/125 	
    N13 = x13.shape[0]/125 	
    N14 = x14.shape[0]/125 
    retv.cloudsatcwc.RVOD_liq_water_content = Numeric.reshape(x11,(125,N11)).astype('Int16')
    retv.cloudsatcwc.LO_RVOD_liquid_water_content = Numeric.reshape(x12,(125,N12)).astype('Int16')
    retv.cloudsatcwc.RVOD_ice_water_content = Numeric.reshape(x13,(125,N13)).astype('Int16')
    retv.cloudsatcwc. IO_RVOD_ice_water_content= Numeric.reshape(x14,(125,N14)).astype('Int16')

    x21 = Numeric.repeat(cloudsatObjcwc.RVOD_liq_water_content_uncertainty[::,0],idx_matchcwc)
    x22 = Numeric.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content_uncertainty[::,0],idx_matchcwc)
    x23 = Numeric.repeat(cloudsatObjcwc.RVOD_ice_water_content_uncertainty[::,0],idx_matchcwc)
    x24 = Numeric.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content_uncertainty[::,0],idx_matchcwc)
    for i in range(1,125):
        x21 = Numeric.concatenate((x21,Numeric.repeat(cloudsatObjcwc.RVOD_liq_water_content_uncertainty[::,i],idx_matchcwc))) 
        x22 = Numeric.concatenate((x22,Numeric.repeat(cloudsatObjcwc.LO_RVOD_liquid_water_content_uncertainty[::,i],idx_matchcwc))) 
        x23 = Numeric.concatenate((x23,Numeric.repeat(cloudsatObjcwc.RVOD_ice_water_content_uncertainty[::,i],idx_matchcwc))) 
        x24 = Numeric.concatenate((x24,Numeric.repeat(cloudsatObjcwc.IO_RVOD_ice_water_content_uncertainty[::,i],idx_matchcwc))) 

    N21 = x21.shape[0]/125 	
    N22 = x22.shape[0]/125 	
    N23 = x23.shape[0]/125 	
    N24 = x24.shape[0]/125 
    retv.cloudsatcwc.RVOD_liq_water_content_uncertainty = Numeric.reshape(x21,(125,N21)).astype('Int8')
    retv.cloudsatcwc.LO_RVOD_liquid_water_content_uncertainty = Numeric.reshape(x22,(125,N22)).astype('Int8')
    retv.cloudsatcwc.RVOD_ice_water_content_uncertainty = Numeric.reshape(x23,(125,N23)).astype('Int8')
    retv.cloudsatcwc. IO_RVOD_ice_water_content_uncertainty = Numeric.reshape(x24,(125,N24)).astype('Int8')

    retv.cloudsatcwc.Temp_min_mixph_K = cloudsatObjcwc.Temp_min_mixph_K
    retv.cloudsatcwc.Temp_max_mixph_K = cloudsatObjcwc.Temp_max_mixph_K

##########################################################################################################################################
    #pdb.set_trace()
    return retv,min_diff,max_diff


# -----------------------------------------------------
if __name__ == "__main__":
    # Testing:
    import string
    import epshdf
    import pps_io
    import numpy.oldnumeric as Numeric
    
    CLOUDSAT_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

    cloudsatfile = "%s/2007151082929_05796_CS_2B-GEOPROF_GRANULE_P_R04_E02.h5"%(CLOUDSAT_DIR)

    # --------------------------------------------------------------------
    write_log("INFO","Read CLOUDSAT data")
    # Read CLOUDSAT Radar data:
    cloudsat = get_cloudsat(cloudsatfile)
    #cloudsat = read_cloudsat(cloudsatfile)

    lonCloudsat = cloudsat.longitude.ravel()
    latCloudsat = cloudsat.latitude.ravel()

    # Test:
    ndim = lonCloudsat.ravel().shape[0]
    idx_match = Numeric.zeros((ndim,),'b')
    idx_match[0:10] = 1

    x = Numeric.repeat(cloudsat.Height[::,0],idx_match)
    for i in range(1,125):
        x = Numeric.concatenate((x,Numeric.repeat(cloudsat.Height[::,i],idx_match)))
    N = x.shape[0]/125
    cloudsat.Height = Numeric.reshape(x,(125,N))
