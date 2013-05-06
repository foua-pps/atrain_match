import ppshdf_cloudproducts
import ppshdf_helpers
from epshdf import *
import netCDF4	
#from mpop.satin.nwcsaf_pps import NwcSafPpsChannel
import numpy as np
#import h5py
import pps_io
import calendar
import datetime

def daysafter4713bc_to_sec1970(bcdate_array):
    #import pps_time_util #@UnresolvedImport
    import datetime
    bcdate=np.min(bcdate_array)
    ddays=np.floor(bcdate)
    dseconds=np.floor(24*60*60*(bcdate-ddays))
    ddays = ddays-2433283

    time_delta=datetime.timedelta(days=ddays, seconds=dseconds)
    the_time = datetime.datetime(1950,1,1,12,0)+time_delta
    sec_1970=calendar.timegm(the_time.timetuple())
    sec_1970_array = 24*60*60*(bcdate_array - bcdate) + sec_1970
    return sec_1970_array


filename="20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
def cci_read_ctth(filename):
    filename="20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
    #filename_pps="noaa18_20080613_0022_99999_satproj_00000_13793_cloudtype.h5"
    #ct=h5py.File(filename_pps)
    #ct.attrs['sec_1970']

    #ctype_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    #dir(ctype_nc)
    #ctype_nc.variables
    #ctype_nc.variables['cth'].add_offset
    #ctype_nc.variables['cth'][:] 
    #ctype_nc.variables['cth'].scale_factor

    ctype_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    ctth = read_cci_ctth(ctype_nc)
    ctype = read_cci_ctype(ctype_nc)
    avhrrGeoObj = read_cci_geoobj(ctype_nc)
    avhrrAngObj = read_cci_angobj(ctype_nc)
    surft = None
    cppLwp = None
    cppCph = None
    avhrrObj = None
    print ctth
    print ctype
    print avhrrGeoObj
    print avhrrAngObj
    
    return avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, surft, cppLwp, cppCph

def read_cci_ctype(ctype_nc):
    #epshdf.read_cloudtype(pps_files.cloudtype, 1, 1, 0)
    #ctype_nc.variables['cc_total'][:])
    #ctype_nc.variables['ls_flag'][:])
    a=CloudType()
    a.region=SafRegion()
    #a.region.xsize=10
    #a.region.ysize=10
    a.des="This is a cci cloudtype"
    a.cloudtype_des="This is the cloudtype description"
    a.qualityflag_des="This is the qualityflag description"
    a.phaseflag_des="This is the phaseflag description"
    a.sec_1970="???"
    a.satellite_id="????"


    a.cloudtype_lut.append("This is the first cloudtype lut")
    a.cloudtype_lut.append("This is the second cloudtype lut")
    a.qualityflag_lut.append("This is the first qualityflag lut")
    a.qualityflag_lut.append("This is the second qualityflag lut")
    a.phaseflag_lut.append("This is the first phaseflag lut")
    a.phaseflag_lut.append("This is the second phaseflag lut")

    #pps 1: cloudfree land 2:cloudfree sea
    #cci lsflag 0:sea 1:land
    a.cloudtype = 2 - ctype_nc.variables['lsflag'][0:10000,::]
    #if cloud cover over 0.5 set ctype to 9: Medium level cumiliform cloud
    a.cloudtype = np.where(ctype_nc.variables['cc_total'][0:10000,::]>0.5,9,a.cloudtype)
    a.qualityflag= a.cloudtype
    a.phaseflag=None
    return a

def read_cci_angobj(ctype_nc):
    #pps_io.readSunSatAngles(pps_files.sunsatangles) 
    avhrrAngObj = SunSatAngleData()
    avhrrAngObj.satz.data=(
        ctype_nc.variables['satellite_zenith_view_no1'][0:10000,::] 
        - avhrrAngObj.satz.intercept )/avhrrAngObj.satz.gain
    avhrrAngObj.sunz.data=(
        ctype_nc.variables['solar_zenith_view_no1'][0:10000,::]
        -avhrrAngObj.sunz.intercept )/avhrrAngObj.sunz.gain
    avhrrAngObj.azidiff= None #ctype_nc.variables['rel_azimuth_view_no1'][:]-AngObj.satz.gain 
    return avhrrAngObj

def read_cci_geoobj(ctype_nc):
    #avhrrGeoObj = pps_io.readAvhrrGeoData(avhrr_file)
    avhrrGeoObj = pps_io.GeoLocationData()
#    ctype_nc.variables['lon'].add_offset
    avhrrGeoObj.longitude=ctype_nc.variables['lon'][0:10000,::] 
#    ctype_nc.variables['lon'].scale_factor
#    ctype_nc.variables['lat'].add_offset
    avhrrGeoObj.latitude = ctype_nc.variables['lat'][0:10000,::] 
    #avhrrGeoObj.latitude = np.where(avhrrGeoObj.latitude<0,avhrrGeoObj.latitude+180,avhrrGeoObj.latitude)
#    ctype_nc.variables['lat'].scale_factor
    #for pps these are calculated in calipso.py, but better read them
    # from file because time are already available on arrays in teh netcdf files
    dsec = calendar.timegm((1993,1,1,0,0,0,0,0,0)) #TAI to UTC
    time_temp= daysafter4713bc_to_sec1970(ctype_nc.variables['time'][0:10000,::])
    avhrrGeoObj.time = time_temp[0:10000,::]#[:,0] #time_temp[:,0]

    avhrrGeoObj.sec1970_start = np.min(avhrrGeoObj.time)  
    avhrrGeoObj.sec1970_end = np.max(avhrrGeoObj.time)
    #avhrrGeoObj.time = np.linspace(avhrrGeoObj.sec1970_start,avhrrGeoObj.sec1970_end, time_temp.shape[0])
    import time
    print time.gmtime(avhrrGeoObj.sec1970_start)
    print time.gmtime(avhrrGeoObj.sec1970_end)
    #t=datetime.datetime(2008,06,13,00,22)
    #avhrrGeoObj.sec1970_start=calendar.timegm(t.timetuple())
    #avhrrGeoObj.sec1970_end=avhrrGeoObj.sec1970_start+110*60
    #avhrrGeoObj.time=None
    #print avhrrGeoObj.sec1970_start
    #print avhrrGeoObj.sec1970_end
    return  avhrrGeoObj

def read_cci_ctth(ctype_nc):
    a = ppshdf_cloudproducts.CloudTop()
    a.region = ppshdf_helpers.SafRegion() 
    #a.region.xsize=10
    #a.region.ysize=10
    cth=ctype_nc.variables['cth'][0:10000,::]
    cth_data=np.where(cth.mask,-999,cth.data)#*ctype_nc.variables['cth'].scale_factor+ctype_nc.variables['cth'].add_offset)
    
    a.h_gain = 200.0
    a.h_intercept = 0.0
    a.h_nodata = 255
    a.h_valid_max = 254
    a.h_valid_min = 0
    data=(1/a.h_gain)*(1000*cth_data-a.h_intercept) #*1000 CCI data is in km
    a.height = np.where(cth.mask,a.h_nodata,data)
    
    
    ctt=ctype_nc.variables['ctt'][0:10000,::]
    ctt_data=np.where(ctt.mask,-999,ctt.data)#*ctype_nc.variables['cth'].scale_factor+ctype_nc.variables['cth'].add_offset)
    a.t_gain = 1.0
    a.t_intercept = 100
    a.t_nodata = 255
    a.t_valid_max = 254
    a.t_valid_min = 0
    data=(1/a.t_gain)*(ctt_data-a.t_intercept)
    a.temperature =  np.where(ctt.mask,a.t_nodata,data)
    
    ctp=ctype_nc.variables['ctp'][0:10000,::]
    ctp_data=np.where(ctp.mask,-999,ctp.data)#*ctype_nc.variables['cth'].scale_factor+ctype_nc.variables['cth'].add_offset)
    a.p_gain = 25.0
    a.p_intercept = 0.0
    a.p_nodata = 255
    a.p_valid_max = 254
    a.p_valid_min = 0
    data=(1/a.p_gain)*(ctp_data-a.p_intercept)
    a.pressure =  np.where(ctp.mask,a.p_nodata,data)


    a.des = "CLOUD CCI Cloud Top Temperature & Height"
    a.ctt_des = "This is the ctt_des description"
    a.cth_des = "This is the cth_des description"
    a.ctp_des = "This is the ctp_des description"
    a.cloudiness_des = "This is the cloudiness description"
    a.processingflag_des = "This is the processingflag description"
    a.sec_1970 = "Is this used in atrain_match???"
    a.satellite_id = "???"

    a.processingflag_lut = []

    a.c_nodata = 255
    a.cloudiness = None
    return a

if __name__ == "__main__":
    filename="20080613002200-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv1.0.nc"
    pps_objects = cci_read_ctth(filename)
    print pps_objects
