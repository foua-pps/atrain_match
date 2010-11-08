#HISTORY CHANGES BY ANKE TETZLAFF#

#080430:
# Could not run in the Arctic using the tmpaid creation
# Used hard coded area 'arctic_super_1002_5km' instead

#080416: Got error message in line 232 and 234:
#         "only rank-0 arrays can be converted to Python scalars";
#         changed start_sec1970 and end_sec1970 to rank0 array 
#         by setting start_sec1970[0] and end_sec1970[0]
import inspect

from pps_basic_configure import *
from pps_error_messages import *

from calipso5km import *
from config import AREA5KM, SUB_DIR, DATA_DIR

#MAIN_DIR = "/data/proj/safworks/adam/cloudsat5km_data"
#MAIN_DIR = "/data/proj_nsc1/safworks/kgkarl/ORR-B-datasets/cloudsat5km/"
#MAIN_DIR = "/data/proj_nsc1/safworks/calipso_cloudsat5km/data/arctic"
#SUB_DIR = "matchups/"
#SUB_DIR = "metop02_calipso_2007spring"
#SUB_DIR = "noaa18_calipso_cloudsat5km_2007DEC_KG"
#SUB_DIR = "noaa17_calipso_cloudsat5km_2007JUN_KG"
#SUB_DIR = "noaa18_calipso_cloudsat5km_2007JUN_KG"

#SATPOS_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
#EPHE_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
#COMPRESS_LVL = 6

COVERAGE_DIR = "%s/5km/%s"%(SUB_DIR,AREA5KM)
#NLINES=1000
#NLINES=6000
#SWATHWD=2048
#SWATHWD5km=409
#NODATA=-9


class Cloudsat5kmObject:
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

class Cloudsat5kmAvhrr5kmTrackObject:
    def __init__(self):
        self.avhrr5km=ppsAvhrr5kmObject()
        self.cloudsat5km=Cloudsat5kmObject()
        self.diff_sec_1970=None

# ----------------------------------------
def readCloudsat5kmAvhrr5kmMatchObj(filename5km):
    import _pyhl

    retv5km = Cloudsat5kmAvhrr5kmTrackObject()

    a=_pyhl.read_nodelist(filename5km)
    a.selectAll()
    a.fetch()

    # Match-Up - time difference:
    retv5km.diff_sec_1970 = a.getNode("/diff_sec_1970").data()

    # cloudsat5km:
    retv5km.cloudsat5km.longitude = a.getNode("/cloudsat/longitude").data()
    retv5km.cloudsat5km.latitude = a.getNode("/cloudsat/latitude").data()
    retv5km.cloudsat5km.avhrr_linnum = a.getNode("/cloudsat/avhrr_linnum").data()
    retv5km.cloudsat5km.avhrr_pixnum = a.getNode("/cloudsat/avhrr_pixnum").data()
    
    retv5km.cloudsat5km.cloud_mask = a.getNode("/cloudsat/cloud_mask").data()
    retv5km.cloudsat5km.Radar_Reflectivity = a.getNode("/cloudsat/Radar_Reflectivity").data()
    retv5km.cloudsat5km.Height = a.getNode("/cloudsat/Height").data()
    retv5km.cloudsat5km.echo_top = a.getNode("/cloudsat/echo_top").data()
    retv5km.cloudsat5km.SurfaceHeightBin = a.getNode("/cloudsat/SurfaceHeightBin").data()
    retv5km.cloudsat5km.SurfaceHeightBin_fraction = a.getNode("/cloudsat/SurfaceHeightBin_fraction").data()
    
    retv5km.cloudsat5km.elevation = a.getNode("/cloudsat/elevation").data()
    retv5km.cloudsat5km.sec_1970 = a.getNode("/cloudsat/sec_1970").data()
    retv5km.cloudsat5km.MODIS_Cloud_Fraction = a.getNode("/cloudsat/MODIS_Cloud_Fraction").data()
    retv5km.cloudsat5km.MODIS_cloud_flag = a.getNode("/cloudsat/MODIS_cloud_flag").data()

    # AVHRR:
    retv5km.avhrr5km.longitude = a.getNode("/avhrr/longitude").data()
    retv5km.avhrr5km.latitude = a.getNode("/avhrr/latitude").data()
    retv5km.avhrr5km.sec_1970 = a.getNode("/avhrr/sec_1970").data()
    retv5km.avhrr5km.cloudtype = a.getNode("/avhrr/cloudtype").data()
    retv5km.avhrr5km.ctth_height = a.getNode("/avhrr/ctth_height").data()
    retv5km.avhrr5km.ctth_pressure = a.getNode("/avhrr/ctth_pressure").data()
    retv5km.avhrr5km.ctth_temperature = a.getNode("/avhrr/ctth_temperature").data()
    retv5km.avhrr5km.bt11micron = a.getNode("/avhrr/bt11micron").data()
    retv5km.avhrr5km.bt12micron = a.getNode("/avhrr/bt12micron").data()
    retv5km.avhrr5km.surftemp = a.getNode("/avhrr/surftemp").data()
    retv5km.avhrr5km.satz = a.getNode("/avhrr/satz").data()

    return retv5km


# ----------------------------------------
def writeCloudsat5kmAvhrr5kmMatchObj(filename5km,ca_obj,compress_lvl):
    import _pyhl
    status = -1

    a=_pyhl.nodelist()

    shape = [ca_obj.cloudsat5km.longitude.shape[0]]

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
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.longitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/latitude")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.latitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/avhrr_linnum")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.avhrr_linnum,"int",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/avhrr_pixnum")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.avhrr_pixnum,"int",-1)
    a.addNode(b)

    shape2d = ca_obj.cloudsat5km.cloud_mask.shape
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/cloud_mask")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat5km.cloud_mask,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/Radar_Reflectivity")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat5km.Radar_Reflectivity,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/Height")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat5km.Height,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/echo_top")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.echo_top,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/SurfaceHeightBin")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.SurfaceHeightBin,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/SurfaceHeightBin_fraction")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.SurfaceHeightBin_fraction,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/elevation")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.elevation,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/sec_1970")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.sec_1970,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/MODIS_Cloud_Fraction")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.MODIS_Cloud_Fraction,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/MODIS_cloud_flag")
    b.setArrayValue(1,shape,ca_obj.cloudsat5km.MODIS_cloud_flag,"double",-1)
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

    status = a.write(filename5km,compress_lvl)

    return status

# -----------------------------------------------------
def select_cloudsat5km_inside_avhrr5km(cloudsat5kmObj,cal5km,sec1970_start_end,sec_timeThr):
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
def get_cloudsat5km(filename5km):
    import _pypps_filters
    #import numpy
    import time

    # Read CLOUDSAT Radar data:
    cloudsat5km = read_cloudsat5km(filename5km)

    # GEOPROF
    # ====
    lonCloudsat5km = cloudsat5km.longitude[:,1].ravel()
    latCloudsat5km = cloudsat5km.latitude[:,1].ravel()
    ndim5km = lonCloudsat5km.shape[0]

    # --------------------------------------------------------------------
    # Time in seconds since 1970:
    
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone

    print "GEOPROF TAI start time: ",cloudsat5km.TAI_start
    start_sec1970 = cloudsat5km.TAI_start + dsec
    start_time = time.gmtime(start_sec1970[0])
    end_sec1970 = cloudsat5km.TAI_start + dsec + cloudsat5km.Profile_time[:,1][ndim5km-1]
    end_time = time.gmtime(end_sec1970[0])
    print "GEOPROF Start and end times: ",start_time,end_time
    cloudsat5km.sec1970 = cloudsat5km.Profile_time[:,1] + start_sec1970

#    ClTime=CloudSatStartEndObject()
#    ClTime.CloudSatStart = start_sec1970
#    ClTime.CloudSatEnd = end_sec1970
#    filename_ClTime = './testdata/5km/matched_files/ClTime/Temp.h5'
#    if not os.path.exists('./testdata/5km/matched_files/ClTime/'):
#        os.makedirs('./testdata/5km/matched_files/ClTime/')
#    
#    writeCloudSatStartEnd(filename_ClTime,ClTime)

    return cloudsat5km

# -----------------------------------------------------
def read_cloudsat5km(filename5km):
    import _pyhl
    import numpy

    # GEOPROF
    # ====
    a=_pyhl.read_nodelist(filename5km)
     
    b=a.getNodeNames()
    a.selectAll()
    a.fetch()

    retv5km = Cloudsat5kmObject()
    root="/cloudsat"

    d=a.getNode("%s/Geolocation Fields/Longitude"%root)
    retv5km.longitude=d.data().astype('d')
    d=a.getNode("%s/Geolocation Fields/Latitude"%root)
    retv5km.latitude=d.data().astype('d')
    
    # International Atomic Time (TAI) seconds from Jan 1, 1993:
    d=a.getNode("%s/Geolocation Fields/Profile_time"%root)
    retv5km.Profile_time=d.data().astype('d')
    d=a.getNode("%s/Geolocation Fields/TAI_start"%root)
    retv5km.TAI_start=d.data().astype('d')

    d=a.getNode("%s/Geolocation Fields/DEM_elevation"%root)
    retv5km.elevation=d.data().astype('d')
   
    # The radar data:
    # 2-d arrays of dimension Nx125:
    d=a.getNode("%s/Data Fields/CPR_Cloud_mask"%root)
    retv5km.CPR_Cloud_mask=d.data()

    d=a.getNode("%s/Geolocation Fields/Height"%root)
    retv5km.Height=d.data()

    d=a.getNode("%s/Data Fields/Radar_Reflectivity"%root)
    retv5km.Radar_Reflectivity=d.data().astype('d')
    d=a.getNode("%s/Data Fields/Gaseous_Attenuation"%root)
    retv5km.Gaseous_Attenuation=d.data().astype('d')

    d=a.getNode("%s/Data Fields/CPR_Echo_Top"%root)
    retv5km.CPR_Echo_Top=d.data().astype('d')
    d=a.getNode("%s/Data Fields/Clutter_reduction_flag"%root)
    retv5km.Clutter_reduction_flag=d.data().astype('d')
    d=a.getNode("%s/Data Fields/Data_quality"%root)
    retv5km.Data_quality=d.data().astype('d')
    d=a.getNode("%s/Data Fields/Data_targetID"%root)
    retv5km.Data_targetID=d.data().astype('d')
    d=a.getNode("%s/Data Fields/MODIS_Cloud_Fraction"%root)
    retv5km.MODIS_Cloud_Fraction=d.data().astype('d')
    d=a.getNode("%s/Data Fields/MODIS_cloud_flag"%root)
    retv5km.MODIS_cloud_flag=d.data().astype('d')
    d=a.getNode("%s/Data Fields/Sigma-Zero"%root)
    retv5km.SigmaZero=d.data().astype('d')
    d=a.getNode("%s/Data Fields/SurfaceHeightBin"%root)
    retv5km.SurfaceHeightBin=d.data().astype('d')
    d=a.getNode("%s/Data Fields/SurfaceHeightBin_fraction"%root)
    retv5km.SurfaceHeightBin_fraction=d.data().astype('d')
    d=a.getNode("%s/Data Fields/sem_NoiseFloor"%root)
    retv5km.sem_NoiseFloor=d.data().astype('d')
    d=a.getNode("%s/Data Fields/sem_NoiseFloorVar"%root)
    retv5km.sem_NoiseFloorVar=d.data().astype('d')
    d=a.getNode("%s/Data Fields/sem_NoiseGate"%root)
    retv5km.sem_NoiseGate=d.data().astype('d')

    return retv5km
# --------------------------------------------
def get_cloudsat5km_avhrr5km_linpix(avhrr5kmIn,avhrrname,lon,lat,clTime):
    import numpy

    tmppcs="tmpproj"
    define_pcs5km(tmppcs, "Plate Caree, central meridian at 15E",
               ['proj=eqc','ellps=bessel', 'lon_0=15'])

    orbittime =  os.path.basename(avhrrname).split("_")[1:3]
    
    startline=0
    Inside=0
    HasEncounteredMatch=0
    i=0
        
    if not os.path.exists(COVERAGE_DIR):
        os.makedirs(COVERAGE_DIR)
    while startline < avhrr5kmIn.longitude.shape[0]:
        
        write_log("INFO","Calling get_cloudsat_avhrr_linpix: start-line = ",startline)
        tmpaid = "tmparea_%d"%i
        endline = startline+NLINES
        coverage_filename5km = "%s/coverage_avhrr_cloudsat_matchup_%s_%s_%s_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrr5kmIn.satellite,orbittime[0],orbittime[1],
                             startline,endline,tmpaid)

        write_log("INFO","Coverage filename5km = ",coverage_filename5km)
        cal,cap,ok = get_cloudsat5km_avhrr5km_linpix_segment(avhrr5kmIn,lon,lat,clTime,
                                                      (startline,endline),
                                                      SWATHWD5km,tmppcs,tmpaid,
                                                      coverage_filename5km)
        if ok:
            HasEncounteredMatch=1
            write_log("INFO","There was a match...")
            
        #if not ok and HasEncounteredMatch:
        #    write_log("INFO","Data is now empty. Leave the loop...")
        #    break
        
        if(startline==0):
            # First time:
            cloudsat_avhrr5km_line,cloudsat_avhrr5km_pixel = numpy.array(cal),numpy.array(cap)
        else:
            # Merge:
            cloudsat_avhrr5km_line = numpy.where(numpy.equal(cloudsat_avhrr5km_line,-9),cal,cloudsat_avhrr5km_line)
            cloudsat_avhrr5km_pixel = numpy.where(numpy.equal(cloudsat_avhrr5km_pixel,-9),cap,cloudsat_avhrr5km_pixel)

        startline=startline+NLINES
        i=i+1
##########################################################################################################################################
        #pdb.set_trace() 
    return cloudsat_avhrr5km_line,cloudsat_avhrr5km_pixel

# --------------------------------------------
def get_cloudsat5km_avhrr5km_linpix_segment(avhrr5kmIn,lon5km,lat5km,time5km,lines,swath_width5km,tmppcs,
                                      tmpaid,covfilename5km):
    import numpy
    import _satproj
    import area,pcs
    import pps_gisdata

    ndim5km= lon5km.shape[0]
    
    if avhrr5kmIn.longitude.shape[0] > lines[1]:
        lines_end = lines[1]
    else:
        lines_end = avhrr5kmIn.longitude.shape[0]
        
    lines_start = lines[0]

    write_log("INFO","lines_start,lines_end: ",lines_start,lines_end)
    nlines = lines_end - lines_start
    lonarr5km = avhrr5kmIn.longitude[lines_start:lines_end,::]
    latarr5km = avhrr5kmIn.latitude[lines_start:lines_end,::]

    idx_start = lines_start*swath_width5km
    idx_end   = lines_end*swath_width5km
    idx5km = numpy.arange(idx_start,idx_end)
    
    linearr5km = numpy.divide(idx5km,swath_width5km)
    write_log("INFO","Start and end line numbers: ",linearr5km[0],linearr5km[idx5km.shape[0]-1])

    linearr5km = numpy.reshape(linearr5km,(nlines,swath_width5km))
    pixelarr5km = numpy.fmod(idx5km,swath_width5km).astype('l')
    pixelarr5km = numpy.reshape(pixelarr5km,(nlines,swath_width5km))

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

    
#    if not os.path.exists(covfilename5km):
#        write_log("INFO","Create Coverage map...")
#        cov5km= _satproj.create_coverage(areaObj,lonarr5km,latarr5km,1)
#        print covfilename5km
#        writeCoverage5km(cov5km,covfilename5km,"satproj",AREA5KM)
#    else:
#        write_log("INFO","Read the AVHRR-CLOUDSAT matchup coverage from file...")
#        cov5km,info = readCoverage5km(covfilename5km)
    
    write_log("INFO","Create Coverage map...")
    cov5km= _satproj.create_coverage(areaObj,lonarr5km,latarr5km,0)
    print covfilename5km
    writeCoverage5km(cov5km,covfilename5km,"satproj",AREA5KM)

    mapped_line = _satproj.project(cov5km.coverage,cov5km.rowidx,cov5km.colidx,linearr5km,NODATA)
    mapped_pixel = _satproj.project(cov5km.coverage,cov5km.rowidx,cov5km.colidx,pixelarr5km,NODATA)


    write_log("INFO","Go through cloudsat track:")
    cloudsat5km_avhrr5km_line = []
    cloudsat5km_avhrr5km_pixel = []
    cloudsat5km_avhrr5km_line_time = []
    cloudsat5km_avhrr5km_pixel_time = []
    kk=0
    for i in range(ndim5km):
        xy_tup=pps_gisdata.lonlat2xy(AREA5KM,lon5km[i],lat5km[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
##        dimx=4500 Should be 5010!!!/KG
##        dimy=4500
        dimx=mapped_line.shape[1]#1002#5010
        dimy=mapped_line.shape[0]#1002#5010
        #if i== 1301:
        #        pdb.set_trace()
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            #if mapped_line[y,x]== 5512:
            #    pdb.set_trace()
            cloudsat5km_avhrr5km_line.append(mapped_line[y,x])
            cloudsat5km_avhrr5km_pixel.append(mapped_pixel[y,x])
            cloudsat5km_avhrr5km_line_time.append(-9)
            cloudsat5km_avhrr5km_pixel_time.append(-9)
        else:
            
            cloudsat5km_avhrr5km_line.append(-9)
            cloudsat5km_avhrr5km_pixel.append(-9)
            cloudsat5km_avhrr5km_line_time.append(-9)
            cloudsat5km_avhrr5km_pixel_time.append(-9)
            kk=kk+1
            
    cloudsat5km_avhrr5km_line = numpy.array(cloudsat5km_avhrr5km_line)
    cloudsat5km_avhrr5km_pixel = numpy.array(cloudsat5km_avhrr5km_pixel)
    cloudsat5km_avhrr5km_line_time = numpy.array(cloudsat5km_avhrr5km_line_time)
    cloudsat5km_avhrr5km_pixel_time = numpy.array(cloudsat5km_avhrr5km_pixel_time)
    # Control the time diference
    match_cloudsat_points = numpy.where(numpy.not_equal(cloudsat5km_avhrr5km_line,-9))
    avhrr_time = (cloudsat5km_avhrr5km_line[match_cloudsat_points] * DSEC_PER_AVHRR_SCALINE5KM) + avhrr5kmIn.sec1970_start
    cl_time = time5km[match_cloudsat_points]
    time_diff = avhrr_time-cl_time
    #max_time_diff_allowed = 50*60 #Based on that a lap is 102 min
    max_time_diff_allowed = sec_timeThr
    time_match = numpy.where(abs(time_diff)<max_time_diff_allowed)
    if time_match[0].shape[0]==0:
        x=numpy.repeat(cloudsat5km_avhrr5km_line_time,numpy.not_equal(cloudsat5km_avhrr5km_line_time,-9))
    else:
        cloudsat5km_avhrr5km_line_time[match_cloudsat_points[0][time_match]]= cloudsat5km_avhrr5km_line[match_cloudsat_points[0][time_match]]
        cloudsat5km_avhrr5km_pixel_time[match_cloudsat_points[0][time_match]] = cloudsat5km_avhrr5km_pixel[match_cloudsat_points[0][time_match]]
        x=numpy.repeat(cloudsat5km_avhrr5km_line_time,numpy.not_equal(cloudsat5km_avhrr5km_line_time,-9)) 
        
    write_log("INFO","Number of matching points = ",x.shape[0])
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0
    ##########################################################################################################################################
    #pdb.set_trace()
    return cloudsat5km_avhrr5km_line_time,cloudsat5km_avhrr5km_pixel_time,matchOk

# -----------------------------------------------------
def match_cloudsat5km_avhrr5km(ctypefile5km,cloudsat5kmObj,avhrr5kmGeoObj,avhrr5kmObj,ctype5km,ctth5km,surft5km,avhrr5kmAngObj):
    import numpy
    import time
    import string
    
    retv5km = Cloudsat5kmAvhrr5kmTrackObject()

    lonCloudsat5km = cloudsat5kmObj.longitude[:,1].ravel()
    latCloudsat5km = cloudsat5kmObj.latitude[:,1].ravel()
    timeCloudsat5km = cloudsat5kmObj.sec1970.ravel()
    ndim5km = lonCloudsat5km.shape[0]

    # Add time
    if avhrr5kmGeoObj.sec1970_end<avhrr5kmGeoObj.sec1970_start:
        avhrr5km_sec1970_end = int(DSEC_PER_AVHRR_SCALINE5KM*avhrr5kmObj.num_of_lines+avhrr5kmGeoObj.sec1970_start)
        avhrr5kmGeoObj.sec1970_end = avhrr5km_sec1970_end

    # --------------------------------------------------------------------

    cal5km,cap5km = get_cloudsat5km_avhrr5km_linpix(avhrr5kmGeoObj,ctypefile5km,lonCloudsat5km,latCloudsat5km,timeCloudsat5km)

    #print len(cal), len(cap), ndim

    idx_cloudsat5km = numpy.arange(ndim5km)
    idx_cloudsat5km_on_avhrr5km = numpy.repeat(idx_cloudsat5km,numpy.not_equal(cal5km,NODATA))    
    lon_cloudsat5km = numpy.repeat(lonCloudsat5km,numpy.not_equal(cal5km,NODATA))
    lat_cloudsat5km = numpy.repeat(latCloudsat5km,numpy.not_equal(cal5km,NODATA))

    # Cloudsat line,pixel inside AVHRR swath:
    cal_on_avhrr5km = numpy.repeat(cal5km,numpy.not_equal(cal5km,NODATA))
    cap_on_avhrr5km = numpy.repeat(cap5km,numpy.not_equal(cal5km,NODATA))
    timetup_start = time.gmtime(avhrr5kmGeoObj.sec1970_start)
    timetup_end   = time.gmtime(avhrr5kmGeoObj.sec1970_end)

    #sec_timeThr = 60*20 # Allow for 20 minute deviation between AVHRR and CLOUDSAT matchup
    secTup = avhrr5kmGeoObj.sec1970_start,avhrr5kmGeoObj.sec1970_end
    idx_match5km = select_cloudsat5km_inside_avhrr5km(cloudsat5kmObj,cal5km,secTup,sec_timeThr)

    print "Make CPR echo top array..."
    retv5km.cloudsat5km.echo_top = numpy.repeat(cloudsat5kmObj.CPR_Echo_Top,idx_match5km).astype('d')#.astype('b')

    retv5km.cloudsat5km.sec_1970 = numpy.repeat(cloudsat5kmObj.sec1970,idx_match5km)

    retv5km.cloudsat5km.latitude = numpy.repeat(latCloudsat5km,idx_match5km)
    retv5km.cloudsat5km.longitude = numpy.repeat(lonCloudsat5km,idx_match5km)

    retv5km.cloudsat5km.MODIS_Cloud_Fraction = numpy.repeat(cloudsat5kmObj.MODIS_Cloud_Fraction,idx_match5km)
    retv5km.cloudsat5km.MODIS_cloud_flag = numpy.repeat(cloudsat5kmObj.MODIS_cloud_flag,idx_match5km)
    
    print "cap_on_avhrr.shape: ",cap_on_avhrr5km.shape
    retv5km.cloudsat5km.avhrr_linnum = cal_on_avhrr5km.astype('i')
    retv5km.cloudsat5km.avhrr_pixnum = cap_on_avhrr5km.astype('i')
    
    print "Make CPR cloud mask array..."
    x = numpy.repeat(cloudsat5kmObj.CPR_Cloud_mask[::,0],idx_match5km)
    for i in range(1,125):
        x = numpy.concatenate((x,numpy.repeat(cloudsat5kmObj.CPR_Cloud_mask[::,i],idx_match5km)))
    N = x.shape[0]/125
    retv5km.cloudsat5km.cloud_mask = numpy.reshape(x,(125,N)).astype('d')

    missing_data = -9.9
    print "Make Radar reflectivity array..."
    x = numpy.repeat(cloudsat5kmObj.Radar_Reflectivity[::,0],idx_match5km)
    for i in range(1,125):
        x = numpy.concatenate((x,numpy.repeat(cloudsat5kmObj.Radar_Reflectivity[::,i],idx_match5km)))
    N = x.shape[0]/125
    RadarRefl = numpy.reshape(x,(125,N))
    RadarRefl = numpy.where(numpy.less(RadarRefl,0),missing_data,RadarRefl)
    retv5km.cloudsat5km.Radar_Reflectivity = RadarRefl.astype('d')#.astype('Int16')

    x = numpy.repeat(cloudsat5kmObj.Height[::,0],idx_match5km)
    for i in range(1,125):
        x = numpy.concatenate((x,numpy.repeat(cloudsat5kmObj.Height[::,i],idx_match5km)))
    N = x.shape[0]/125
    retv5km.cloudsat5km.Height = numpy.reshape(x,(125,N)).astype('d')#.astype('Int16')

    # One-d arrays:
    retv5km.cloudsat5km.SurfaceHeightBin = \
                       numpy.repeat(cloudsat5kmObj.SurfaceHeightBin.ravel(),
                                      idx_match5km.ravel()).astype('d')#.astype('b')
    retv5km.cloudsat5km.SurfaceHeightBin_fraction = \
                       numpy.repeat(cloudsat5kmObj.SurfaceHeightBin_fraction.ravel(),
                                      idx_match5km.ravel()).astype('d')

    print "cloudsat5km observation time of first cloudsat5km-avhrr match: ",\
          time.gmtime(retv5km.cloudsat5km.sec_1970[0])
    print "cloudsat5km observation time of last cloudsat5km-avhrr match: ",\
          time.gmtime(retv5km.cloudsat5km.sec_1970[N-1])

    # Elevation is given in km's. Convert to meters:
    retv5km.cloudsat5km.elevation = numpy.repeat(\
        cloudsat5kmObj.elevation.ravel(),idx_match5km.ravel()).astype('d')#.astype('Int16')

    retv5km.avhrr5km.sec_1970 = numpy.add(avhrr5kmGeoObj.sec1970_start,
                                      cal_on_avhrr5km * DSEC_PER_AVHRR_SCALINE5KM)
    retv5km.diff_sec_1970 = retv5km.cloudsat5km.sec_1970 - retv5km.avhrr5km.sec_1970
    min_diff = numpy.minimum.reduce(retv5km.diff_sec_1970)
    max_diff = numpy.maximum.reduce(retv5km.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-cloudsat5km): ",\
          numpy.maximum.reduce(retv5km.diff_sec_1970),numpy.minimum.reduce(retv5km.diff_sec_1970)

    print "AVHRR observation time of first cloudsat5km-avhrr match: ",\
          time.gmtime(retv5km.avhrr5km.sec_1970[0])
    print "AVHRR observation time of last cloudsat5km-avhrr match: ",\
          time.gmtime(retv5km.avhrr5km.sec_1970[N-1])

    # Make the latitude and pps cloudtype on the cloudsat5km track:
    # line and pixel arrays have equal dimensions
    print "Generate all datatypes (lat,lon,cty,ctth,surft) on the cloudsat5km track!"

    ctype_track = []
    ctth_height_track = []
    ctth_pressure_track = []
    ctth_temperature_track = []
    lon_avhrr5km_track = []
    lat_avhrr5km_track = []
    surft_track = []
    bt11micron_track = []
    bt12micron_track = []
    satz_track = []
    for i in range(cal_on_avhrr5km.shape[0]):
        lat_avhrr5km_track.append(avhrr5kmGeoObj.latitude[cal_on_avhrr5km[i],cap_on_avhrr5km[i]])
        lon_avhrr5km_track.append(avhrr5kmGeoObj.longitude[cal_on_avhrr5km[i],cap_on_avhrr5km[i]])
        ctype_track.append(ctype5km.cloudtype[cal_on_avhrr5km[i],cap_on_avhrr5km[i]])
        surft_track.append(surft5km[cal_on_avhrr5km[i],cap_on_avhrr5km[i]])
        if avhrr5kmObj.channels[3].data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == avhrr5kmObj.nodata:
            b11 = -9.
        else:
            b11 = avhrr5kmObj.channels[3].data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] * avhrr5kmObj.channels[3].gain + avhrr5kmObj.channels[3].intercept
        bt11micron_track.append(b11)
        if avhrr5kmObj.channels[4].data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == avhrr5kmObj.nodata:
            b12 = -9.
        else:
            b12 = avhrr5kmObj.channels[4].data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] * avhrr5kmObj.channels[4].gain + avhrr5kmObj.channels[4].intercept
        bt12micron_track.append(b12)
        if ctth5km == None:
            continue
        if ctth5km.height[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == ctth5km.h_nodata:
            hh = -9.
        else:
            hh = ctth5km.height[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] * ctth5km.h_gain + ctth5km.h_intercept
        ctth_height_track.append(hh)
        if ctth5km.temperature[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == ctth5km.t_nodata:
            tt = -9.
        else:
            tt = ctth5km.temperature[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] * ctth5km.t_gain + \
                 ctth5km.t_intercept
        ctth_temperature_track.append(tt)
        if ctth5km.pressure[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == ctth5km.p_nodata:
            pp = -9.
        else:
            pp = ctth5km.pressure[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] * ctth5km.p_gain + ctth5km.p_intercept
        ctth_pressure_track.append(pp)
        if avhrr5kmAngObj.satz.data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == avhrr5kmAngObj.satz.no_data or avhrr5kmAngObj.satz.data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] == avhrr5kmAngObj.satz.missing_data:
            ang = -9
        else:
            ang = avhrr5kmAngObj.satz.data[cal_on_avhrr5km[i],cap_on_avhrr5km[i]] * avhrr5kmAngObj.satz.gain + avhrr5kmAngObj.satz.intercept
        satz_track.append(ang)
        
    retv5km.avhrr5km.latitude = numpy.array(lat_avhrr5km_track)
    retv5km.avhrr5km.longitude = numpy.array(lon_avhrr5km_track)
    retv5km.avhrr5km.cloudtype = numpy.array(ctype_track)
    retv5km.avhrr5km.bt11micron = numpy.array(bt11micron_track)
    retv5km.avhrr5km.bt12micron = numpy.array(bt12micron_track)
    retv5km.avhrr5km.satz = numpy.array(satz_track)
    if ctth5km:
        retv5km.avhrr5km.ctth_height = numpy.array(ctth_height_track)
        retv5km.avhrr5km.ctth_pressure = numpy.array(ctth_pressure_track)
        retv5km.avhrr5km.ctth_temperature = numpy.array(ctth_temperature_track)
    retv5km.avhrr5km.surftemp = numpy.array(surft_track)

    print "AVHRR-PPS Cloud Type,latitude: shapes = ",\
          retv5km.avhrr5km.cloudtype.shape,retv5km.avhrr5km.latitude.shape

    ll = []
    for i in range(ndim5km):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat5km[i],latCloudsat5km[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat5km[i],latCloudsat5km[i],idx_match5km[i])))
    basename = os.path.basename(ctypefile5km).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_")
    
    datapath = "%s/%s/5km/%s/%s/%s" %(DATA_DIR, base_sat, base_year, base_month, AREA5KM)
    if not os.path.exists(datapath):
        os.makedirs(datapath) 
           
    fd = open("%s/5km_%s_cloudtype_cloudsat-GEOPROF_track2.txt"%(datapath,basename),"w")
    fd.writelines(ll)
    fd.close()

    ll = []
    for i in range(N):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_cloudsat5km[i],lat_cloudsat5km[i],0)))
    fd = open("%s/5km_%s_cloudtype_cloudsat-GEOPROF_track_excl.txt"%(datapath,basename),"w")
    fd.writelines(ll)
    fd.close()

    return retv5km,min_diff,max_diff

#------------------------------------------------------------------------------------------

def reshapeCloudsat5km(cloudsatfiles,avhrr):
    import time
    import numpy
    import sys
    
    clsat = Cloudsat5kmObject()
    if avhrr.sec1970_end<avhrr.sec1970_start:
        avhrr_end = int(DSEC_PER_AVHRR_SCALINE5KM*avhrr.num_of_lines+avhrr.sec1970_start)    
    else:
        avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    
    startCloudsat = get_cloudsat5km(cloudsatfiles[0])
    startCloudsat.Profile_time = numpy.add(startCloudsat.Profile_time,startCloudsat.TAI_start)
    
    for i in range(len(cloudsatfiles)-1):
        newCloudsat = get_cloudsat5km(cloudsatfiles[i+1])
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
        startCloudsat.CPR_Cloud_mask = numpy.concatenate((startCloudsat.CPR_Cloud_mask[0:clsat_break,:],newCloudsat.CPR_Cloud_mask))
        startCloudsat.CPR_Echo_Top = numpy.concatenate((startCloudsat.CPR_Echo_Top[0:clsat_break,:],newCloudsat.CPR_Echo_Top))
        startCloudsat.Clutter_reduction_flag = numpy.concatenate((startCloudsat.Clutter_reduction_flag[0:clsat_break,:],newCloudsat.Clutter_reduction_flag))
        startCloudsat.Data_quality = numpy.concatenate((startCloudsat.Data_quality[0:clsat_break,:],newCloudsat.Data_quality))
        startCloudsat.Data_targetID = numpy.concatenate((startCloudsat.Data_targetID[0:clsat_break,:],newCloudsat.Data_targetID))
        startCloudsat.Gaseous_Attenuation = numpy.concatenate((startCloudsat.Gaseous_Attenuation[0:clsat_break,:],newCloudsat.Gaseous_Attenuation))
        startCloudsat.MODIS_Cloud_Fraction = numpy.concatenate((startCloudsat.MODIS_Cloud_Fraction[0:clsat_break,:],newCloudsat.MODIS_Cloud_Fraction))
        startCloudsat.MODIS_cloud_flag = numpy.concatenate((startCloudsat.MODIS_cloud_flag[0:clsat_break,:],newCloudsat.MODIS_cloud_flag))
        startCloudsat.Radar_Reflectivity = numpy.concatenate((startCloudsat.Radar_Reflectivity[0:clsat_break,:],newCloudsat.Radar_Reflectivity))
        startCloudsat.Height = numpy.concatenate((startCloudsat.Height[0:clsat_break,:],newCloudsat.Height))
        startCloudsat.SigmaZero = numpy.concatenate((startCloudsat.SigmaZero[0:clsat_break,:],newCloudsat.SigmaZero))
        startCloudsat.SurfaceHeightBin = numpy.concatenate((startCloudsat.SurfaceHeightBin[0:clsat_break,:],newCloudsat.SurfaceHeightBin))
        startCloudsat.SurfaceHeightBin_fraction = numpy.concatenate((startCloudsat.SurfaceHeightBin_fraction[0:clsat_break,:],newCloudsat.SurfaceHeightBin_fraction))
        startCloudsat.sem_NoiseFloor = numpy.concatenate((startCloudsat.sem_NoiseFloor[0:clsat_break,:],newCloudsat.sem_NoiseFloor))
        
        startCloudsat.sem_NoiseFloorVar = numpy.concatenate((startCloudsat.sem_NoiseFloorVar[0:clsat_break,:],newCloudsat.sem_NoiseFloorVar))
        startCloudsat.sem_NoiseGate = numpy.concatenate((startCloudsat.sem_NoiseGate[0:clsat_break,:],newCloudsat.sem_NoiseGate))
        

    start_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_start - sec_timeThr))))-1 # Minus one to get one extra, just to be certain
    end_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain

    clsat.sec1970 = startCloudsat.sec1970[start_break:end_break] 
    clsat.longitude = startCloudsat.longitude[start_break:end_break,:] 
    clsat.latitude = startCloudsat.latitude[start_break:end_break,:] 
    clsat.elevation = startCloudsat.elevation[start_break:end_break,:] 
    clsat.Profile_time = startCloudsat.Profile_time[start_break:end_break,:] 
    clsat.CPR_Cloud_mask = startCloudsat.CPR_Cloud_mask[start_break:end_break,:] 
    clsat.CPR_Echo_Top = startCloudsat.CPR_Echo_Top[start_break:end_break,:] 
    clsat.Clutter_reduction_flag = startCloudsat.Clutter_reduction_flag[start_break:end_break,:]  
    clsat.Data_quality = startCloudsat.Data_quality[start_break:end_break,:] 
    clsat.Data_targetID = startCloudsat.Data_targetID[start_break:end_break,:] 
    clsat.Gaseous_Attenuation = startCloudsat.Gaseous_Attenuation[start_break:end_break,:] 
    clsat.MODIS_Cloud_Fraction = startCloudsat.MODIS_Cloud_Fraction[start_break:end_break,:] 
    clsat.MODIS_cloud_flag = startCloudsat.MODIS_cloud_flag[start_break:end_break,:] 
    clsat.Radar_Reflectivity = startCloudsat.Radar_Reflectivity[start_break:end_break,:] 
    clsat.Height = startCloudsat.Height[start_break:end_break,:] 
    clsat.SigmaZero = startCloudsat.SigmaZero[start_break:end_break,:] 
    clsat.SurfaceHeightBin = startCloudsat.SurfaceHeightBin[start_break:end_break,:] 
    clsat.SurfaceHeightBin_fraction = startCloudsat.SurfaceHeightBin_fraction[start_break:end_break,:] 
    clsat.sem_NoiseFloor = startCloudsat.sem_NoiseFloor[start_break:end_break,:] 
    clsat.sem_NoiseFloorVar = startCloudsat.sem_NoiseFloorVar[start_break:end_break,:] 
    clsat.sem_NoiseGate = startCloudsat.sem_NoiseGate[start_break:end_break,:] 
    
    if clsat.Profile_time.shape[0] <= 0:
        print("No time match, please try with some other CloudSat files")
        print("Program cloudsat5km.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)
    else:
        clsat.TAI_start = clsat.Profile_time[0,0]
        clsat.Profile_time = clsat.Profile_time - clsat.TAI_start
    
    return clsat
    
    
    
# -----------------------------------------------------
if __name__ == "__main__":
    # Testing:
    import string
    import epshdf
    import pps_io
    import numpy
    
    CLOUDSAT_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

    cloudsatfile = "%s/2007151082929_05796_CS_2B-GEOPROF_GRANULE_P_R04_E02.h5"%(CLOUDSAT_DIR)

    # --------------------------------------------------------------------
    write_log("INFO","Read CLOUDSAT data")
    # Read CLOUDSAT Radar data:
    cloudsat5km = get_cloudsat5km(cloudsat5kmfile)
    #cloudsat5km = read_cloudsat5km(cloudsat5kmfile)

    loncloudsat5km = cloudsat5km.longitude.ravel()
    latcloudsat5km = cloudsat5km.latitude.ravel()

    # Test:
    ndim = loncloudsat5km.ravel().shape[0]
    idx_match = numpy.zeros((ndim,),'b')
    idx_match[0:10] = 1

    x = numpy.repeat(cloudsat5km.Height[::,0],idx_match)
    for i in range(1,125):
        x = numpy.concatenate((x,numpy.repeat(cloudsat5km.Height[::,i],idx_match)))
    N = x.shape[0]/125
    cloudsat5km.Height = numpy.reshape(x,(125,N))
