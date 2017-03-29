import logging
logging.basicConfig(level=logging.info)
logger = logging.getLogger(__name__)

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
