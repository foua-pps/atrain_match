import ppshdf_cloudproducts
import ppshdf_helpers
from epshdf import (SafRegion, 
                    CloudType, 
                    SunSatAngleData)
import netCDF4	
import numpy as np
import pps_io
import calendar
import datetime
from pps_error_messages import write_log
import time
#import h5py
#from mpop.satin.nwcsaf_pps import NwcSafPpsChannel

def daysafter4713bc_to_sec1970(bcdate_array):
    #import pps_time_util #@UnresolvedImport
    
    bcdate=np.min(bcdate_array)
    ddays=np.floor(bcdate)
    dseconds=np.floor(24*60*60*(bcdate-ddays))
    ddays = ddays-2433283

    time_delta=datetime.timedelta(days=ddays, seconds=dseconds)
    the_time = datetime.datetime(1950,1,1,12,0)+time_delta
    sec_1970=calendar.timegm(the_time.timetuple())
    sec_1970_array = 24*60*60*(bcdate_array - bcdate) + sec_1970
    return sec_1970_array

def cci_read_ctth(filename):
    #filename="20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
    #filename_pps="noaa18_20080613_0022_99999_satproj_00000_13793_cloudtype.h5"
    #ct=h5py.File(filename_pps)
    #ct.attrs['sec_1970']

    #cci_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    #dir(cci_nc)
    #cci_nc.variables
    #cci_nc.variables['cth'].add_offset
    #cci_nc.variables['cth'][:] 
    #cci_nc.variables['cth'].scale_factor
    write_log("INFO", "Opening file %s"%(filename))
    cci_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    write_log("INFO", "Reading ctth ...")
    ctth = read_cci_ctth(cci_nc)
    write_log("INFO", "Reading cloud type ...")
    ctype = read_cci_ctype(cci_nc)
    write_log("INFO", "Reading longitude, latitude and time ...")
    avhrrGeoObj = read_cci_geoobj(cci_nc)
    write_log("INFO", "Reading angles ...")
    avhrrAngObj = read_cci_angobj(cci_nc)
    write_log("INFO", "Not reading surface temperature")
    surft = None
    write_log("INFO", "Not reading cloud liquid water path")
    cppLwp = None
    write_log("INFO", "Not reading cloud phase")
    cppCph = None
    write_log("INFO", "Not reading channel data")
    avhrrObj = None    
    return avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surft, cppLwp, cppCph

def read_cci_ctype(cci_nc):
    #epshdf.read_cloudtype(pps_files.cloudtype, 1, 1, 0)
    #cci_nc.variables['cc_total'][:])
    #cci_nc.variables['ls_flag'][:])
    ctype=CloudType()
    ctype.region=SafRegion()
    #ctype.region.xsize=10
    #ctype.region.ysize=10
    ctype.des="This is a cci cloudtype"
    ctype.cloudtype_des="This is the cloudtype description"
    ctype.qualityflag_des="This is the qualityflag description"
    ctype.phaseflag_des="This is the phaseflag description"
    ctype.sec_1970="???"
    ctype.satellite_id="????"

    ctype.cloudtype_lut.append("This is the first cloudtype lut")
    ctype.cloudtype_lut.append("This is the second cloudtype lut")
    ctype.qualityflag_lut.append("This is the first qualityflag lut")
    ctype.qualityflag_lut.append("This is the second qualityflag lut")
    ctype.phaseflag_lut.append("This is the first phaseflag lut")
    ctype.phaseflag_lut.append("This is the second phaseflag lut")

    #pps 1: cloudfree land 2:cloudfree sea
    #cci lsflag 0:sea 1:land
    ctype.cloudtype = 2 - cci_nc.variables['lsflag'][::]
    #if cloud cover over 0.5 set ctype to 9: Medium level cumiliform cloud
    ctype.cloudtype = np.where(cci_nc.variables['cc_total'][::]>0.5,9,ctype.cloudtype)
    ctype.qualityflag= ctype.cloudtype
    ctype.phaseflag=None
    return ctype

def read_cci_angobj(cci_nc):
    #pps_io.readSunSatAngles(pps_files.sunsatangles) 
    avhrrAngObj = SunSatAngleData()
    avhrrAngObj.satz.data=(
        cci_nc.variables['satellite_zenith_view_no1'][::] 
        - avhrrAngObj.satz.intercept )/avhrrAngObj.satz.gain
    avhrrAngObj.sunz.data=(
        cci_nc.variables['solar_zenith_view_no1'][::]
        -avhrrAngObj.sunz.intercept )/avhrrAngObj.sunz.gain
    avhrrAngObj.azidiff= None #cci_nc.variables['rel_azimuth_view_no1'][:]-AngObj.satz.gain 
    return avhrrAngObj

def read_cci_geoobj(cci_nc):
    #avhrrGeoObj = pps_io.readAvhrrGeoData(avhrr_file)
    avhrrGeoObj = pps_io.GeoLocationData()
    #cci_nc.variables['lon'].add_offset
    avhrrGeoObj.longitude=cci_nc.variables['lon'][::] 
    #cci_nc.variables['lon'].scale_factor
    #cci_nc.variables['lat'].add_offset
    avhrrGeoObj.latitude = cci_nc.variables['lat'][::] 
    #cci_nc.variables['lat'].scale_factor
    #for pps these are calculated in calipso.py, but better read them
    # from file because time are already available on arrays in teh netcdf files
    #dsec = calendar.timegm((1993,1,1,0,0,0,0,0,0)) #TAI to UTC
    time_temp= daysafter4713bc_to_sec1970(cci_nc.variables['time'][::])
    avhrrGeoObj.time = time_temp[::]#[:,0] #time_temp[:,0]

    avhrrGeoObj.sec1970_start = np.min(avhrrGeoObj.time)  
    avhrrGeoObj.sec1970_end = np.max(avhrrGeoObj.time)

    tim1=time.strftime("%Y %m %d %H:%M", time.gmtime(avhrrGeoObj.sec1970_start))
    tim2=time.strftime("%Y %m %d %H:%M", time.gmtime(avhrrGeoObj.sec1970_end))
    write_log("INFO", "Starttime: %s, end time: %s"%(tim1,tim2))
    #t=datetime.datetime(2008,06,13,00,22)
    #avhrrGeoObj.sec1970_start=calendar.timegm(t.timetuple())
    #avhrrGeoObj.sec1970_end=avhrrGeoObj.sec1970_start+110*60
    #avhrrGeoObj.time=None
    #print avhrrGeoObj.sec1970_start
    #print avhrrGeoObj.sec1970_end
    return  avhrrGeoObj

def read_cci_ctth(cci_nc):
    ctth = ppshdf_cloudproducts.CloudTop()
    ctth.region = ppshdf_helpers.SafRegion() 
    #ctth.region.xsize=10
    #ctth.region.ysize=10
    cth=cci_nc.variables['cth'][::]
    cth_data=np.where(cth.mask,-999,cth.data)
    #*cci_nc.variables['cth'].scale_factor+cci_nc.variables['cth'].add_offset)
    
    ctth.h_gain = 200.0
    ctth.h_intercept = 0.0
    ctth.h_nodata = 255
    ctth.h_valid_max = 254
    ctth.h_valid_min = 0
    data=(1/ctth.h_gain)*(1000*cth_data-ctth.h_intercept) #*1000 CCI data in km
    ctth.height = np.where(cth.mask,ctth.h_nodata,data)
    
    
    ctt=cci_nc.variables['ctt'][::]
    ctt_data=np.where(ctt.mask,-999,ctt.data)
    #*cci_nc.variables['cth'].scale_factor+cci_nc.variables['cth'].add_offset)
    ctth.t_gain = 1.0
    ctth.t_intercept = 100
    ctth.t_nodata = 255
    ctth.t_valid_max = 254
    ctth.t_valid_min = 0
    data=(1/ctth.t_gain)*(ctt_data-ctth.t_intercept)
    ctth.temperature =  np.where(ctt.mask,ctth.t_nodata,data)
    
    ctp=cci_nc.variables['ctp'][::]
    ctp_data=np.where(ctp.mask,-999,ctp.data)#Upscaled
    #*cci_nc.variables['cth'].scale_factor+cci_nc.variables['cth'].add_offset)
    ctth.p_gain = 25.0
    ctth.p_intercept = 0.0
    ctth.p_nodata = 255
    ctth.p_valid_max = 254
    ctth.p_valid_min = 0
    data=(1/ctth.p_gain)*(ctp_data-ctth.p_intercept)
    ctth.pressure =  np.where(ctp.mask,ctth.p_nodata,data)


    ctth.des = "CLOUD CCI Cloud Top Temperature & Height"
    ctth.ctt_des = "This is the ctt_des description"
    ctth.cth_des = "This is the cth_des description"
    ctth.ctp_des = "This is the ctp_des description"
    ctth.cloudiness_des = "This is the cloudiness description"
    ctth.processingflag_des = "This is the processingflag description"
    ctth.sec_1970 = "Is this used in atrain_match???"
    ctth.satellite_id = "???"

    ctth.processingflag_lut = []

    ctth.c_nodata = 255
    ctth.cloudiness = None
    return ctth

if __name__ == "__main__":
    filename="20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
    pps_objects = cci_read_ctth(filename)
    print pps_objects
