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
from config import AREA, SUB_DIR, DATA_DIR, sec_timeThr, DSEC_PER_AVHRR_SCALINE
from common import MatchupError, elements_within_range

COVERAGE_DIR = "%s/%skm/%s"%(SUB_DIR,RESOLUTION,AREA)



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

class CloudsatAvhrrTrackObject:
    def __init__(self):
        self.avhrr=ppsAvhrrObject()
        self.cloudsat=CloudsatObject()
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
    retv.avhrr.satz = a.getNode("/avhrr/satz").data()
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
    b.setArrayValue(1,shape2d,ca_obj.cloudsat.cloud_mask,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/Radar_Reflectivity")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat.Radar_Reflectivity,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/Height")
    b.setArrayValue(1,shape2d,ca_obj.cloudsat.Height,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/echo_top")
    b.setArrayValue(1,shape,ca_obj.cloudsat.echo_top,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/SurfaceHeightBin")
    b.setArrayValue(1,shape,ca_obj.cloudsat.SurfaceHeightBin,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/SurfaceHeightBin_fraction")
    b.setArrayValue(1,shape,ca_obj.cloudsat.SurfaceHeightBin_fraction,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/elevation")
    b.setArrayValue(1,shape,ca_obj.cloudsat.elevation,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/sec_1970")
    b.setArrayValue(1,shape,ca_obj.cloudsat.sec_1970,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/MODIS_Cloud_Fraction")
    b.setArrayValue(1,shape,ca_obj.cloudsat.MODIS_Cloud_Fraction.astype('double'),"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/cloudsat/MODIS_cloud_flag")
    b.setArrayValue(1,shape,ca_obj.cloudsat.MODIS_cloud_flag.astype('double'),"double",-1)
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
    
    status = a.write(filename,compress_lvl)
##########################################################################################################################################
    #pdb.set_trace()
    return status

# -----------------------------------------------------
def select_cloudsat_inside_avhrr(cloudsatObj,cal,sec1970_start_end,sec_timeThr):
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
def get_cloudsat(filename):
    import _pypps_filters
    #import numpy
    import time

    # Read CLOUDSAT Radar data:
    cloudsat = read_cloudsat(filename)

    if RESOLUTION == 1:
        lonCloudsat = cloudsat.longitude.ravel()
        latCloudsat = cloudsat.latitude.ravel()
        timeCloudsat = cloudsat.Profile_time.ravel()
    elif RESOLUTION == 5:
        lonCloudsat = cloudsat.longitude[:,1].ravel()
        latCloudsat = cloudsat.latitude[:,1].ravel()
        timeCloudsat = cloudsat.Profile_time[:,1].ravel()
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
    cloudsat.sec1970 = timeCloudsat + start_sec1970

   # --------------------------------------------------------------------

    return cloudsat

# -----------------------------------------------------
def read_cloudsat(filename):
    import _pyhl
    import numpy
    
    # GEOPROF
    # ====
    a=_pyhl.read_nodelist(filename)
     
    b=a.getNodeNames()
    a.selectAll()
    a.fetch()
    retv = CloudsatObject()
    if RESOLUTION == 1:
        root="/2B-GEOPROF"
        d=a.getNode("%s/Geolocation Fields/Longitude"%root)    
        retv.longitude=numpy.fromstring(d.data(),'f').astype('d')
        d=a.getNode("%s/Geolocation Fields/Latitude"%root)
        retv.latitude=numpy.fromstring(d.data(),'f').astype('d')    
        # International Atomic Time (TAI) seconds from Jan 1, 1993:
        d=a.getNode("%s/Geolocation Fields/Profile_time"%root)
        retv.Profile_time=numpy.fromstring(d.data(),'f').astype('d')
        d=a.getNode("%s/Geolocation Fields/TAI_start"%root)
        retv.TAI_start=numpy.fromstring(d.data(),'d').astype('d')
        d=a.getNode("%s/Geolocation Fields/DEM_elevation"%root)
        retv.elevation=numpy.fromstring(d.data(),'Int16').astype('Int16')  
        # The radar data:
        # 2-d arrays of dimension Nx125:
        d=a.getNode("%s/Data Fields/CPR_Cloud_mask"%root)
        retv.CPR_Cloud_mask=d.data()
        d=a.getNode("%s/Geolocation Fields/Height"%root)
        retv.Height=d.data()
        d=a.getNode("%s/Data Fields/Radar_Reflectivity"%root)
        retv.Radar_Reflectivity=d.data().astype('Int16')
        d=a.getNode("%s/Data Fields/Gaseous_Attenuation"%root)
        retv.Gaseous_Attenuation=d.data().astype('Int16')
        d=a.getNode("%s/Data Fields/CPR_Echo_Top"%root)
        retv.CPR_Echo_Top=numpy.fromstring(d.data(),'b').astype('b')
        d=a.getNode("%s/Data Fields/Clutter_reduction_flag"%root)
        retv.Clutter_reduction_flag=numpy.fromstring(d.data(),'b').astype('b')
        d=a.getNode("%s/Data Fields/Data_quality"%root)
        retv.Data_quality=numpy.fromstring(d.data(),'b').astype('b')
        d=a.getNode("%s/Data Fields/Data_targetID"%root)
        retv.Data_targetID=numpy.fromstring(d.data(),'b').astype('b')
        d=a.getNode("%s/Data Fields/MODIS_Cloud_Fraction"%root)
        retv.MODIS_Cloud_Fraction=numpy.fromstring(d.data(),'b').astype('b')
        d=a.getNode("%s/Data Fields/MODIS_cloud_flag"%root)
        retv.MODIS_cloud_flag=numpy.fromstring(d.data(),'b').astype('b')
        d=a.getNode("%s/Data Fields/Sigma-Zero"%root)
        retv.SigmaZero=numpy.fromstring(d.data(),'Int16').astype('Int16')
        d=a.getNode("%s/Data Fields/SurfaceHeightBin"%root)
        retv.SurfaceHeightBin=numpy.fromstring(d.data(),'b').astype('b')
        d=a.getNode("%s/Data Fields/SurfaceHeightBin_fraction"%root)
        retv.SurfaceHeightBin_fraction=numpy.fromstring(d.data(),'f').astype('d')
        d=a.getNode("%s/Data Fields/sem_NoiseFloor"%root)
        retv.sem_NoiseFloor=numpy.fromstring(d.data(),'f').astype('d')
        d=a.getNode("%s/Data Fields/sem_NoiseFloorVar"%root)
        retv.sem_NoiseFloorVar=numpy.fromstring(d.data(),'f').astype('d')
        d=a.getNode("%s/Data Fields/sem_NoiseGate"%root)
        retv.sem_NoiseGate=numpy.fromstring(d.data(),'b').astype('b')
    elif RESOLUTION == 5:
        root="/cloudsat"
        d=a.getNode("%s/Geolocation Fields/Longitude"%root)
        retv.longitude=d.data().astype('d')
        d=a.getNode("%s/Geolocation Fields/Latitude"%root)
        retv.latitude=d.data().astype('d')    
        # International Atomic Time (TAI) seconds from Jan 1, 1993:
        d=a.getNode("%s/Geolocation Fields/Profile_time"%root)
        retv.Profile_time=d.data().astype('d')
        d=a.getNode("%s/Geolocation Fields/TAI_start"%root)
        retv.TAI_start=d.data().astype('d')
        d=a.getNode("%s/Geolocation Fields/DEM_elevation"%root)
        retv.elevation=d.data().astype('d')   
        # The radar data:
        # 2-d arrays of dimension Nx125:
        d=a.getNode("%s/Data Fields/CPR_Cloud_mask"%root)
        retv.CPR_Cloud_mask=d.data()
        d=a.getNode("%s/Geolocation Fields/Height"%root)
        retv.Height=d.data()
        d=a.getNode("%s/Data Fields/Radar_Reflectivity"%root)
        retv.Radar_Reflectivity=d.data().astype('d')
        d=a.getNode("%s/Data Fields/Gaseous_Attenuation"%root)
        retv.Gaseous_Attenuation=d.data().astype('d')
        d=a.getNode("%s/Data Fields/CPR_Echo_Top"%root)
        retv.CPR_Echo_Top=d.data().astype('d')
        d=a.getNode("%s/Data Fields/Clutter_reduction_flag"%root)
        retv.Clutter_reduction_flag=d.data().astype('d')
        d=a.getNode("%s/Data Fields/Data_quality"%root)
        retv.Data_quality=d.data().astype('d')
        d=a.getNode("%s/Data Fields/Data_targetID"%root)
        retv.Data_targetID=d.data().astype('d')
        d=a.getNode("%s/Data Fields/MODIS_Cloud_Fraction"%root)
        retv.MODIS_Cloud_Fraction=d.data().astype('d')
        d=a.getNode("%s/Data Fields/MODIS_cloud_flag"%root)
        retv.MODIS_cloud_flag=d.data().astype('d')
        d=a.getNode("%s/Data Fields/Sigma-Zero"%root)
        retv.SigmaZero=d.data().astype('d')
        d=a.getNode("%s/Data Fields/SurfaceHeightBin"%root)
        retv.SurfaceHeightBin=d.data().astype('d')
        d=a.getNode("%s/Data Fields/SurfaceHeightBin_fraction"%root)
        retv.SurfaceHeightBin_fraction=d.data().astype('d')
        d=a.getNode("%s/Data Fields/sem_NoiseFloor"%root)
        retv.sem_NoiseFloor=d.data().astype('d')
        d=a.getNode("%s/Data Fields/sem_NoiseFloorVar"%root)
        retv.sem_NoiseFloorVar=d.data().astype('d')
        d=a.getNode("%s/Data Fields/sem_NoiseGate"%root)
        retv.sem_NoiseGate=d.data().astype('d')
    return retv
# --------------------------------------------
def get_cloudsat_avhrr_linpix(avhrrIn,avhrrname,lon,lat,clTime):
    import numpy

    tmppcs="tmpproj"
    define_pcs(tmppcs, "Plate Caree, central meridian at 15E",
               ['proj=eqc','ellps=bessel', 'lon_0=15'])
    orbittime =  os.path.basename(avhrrname).split("_")[1:3]
    
    startline=0
    Inside=0
    HasEncounteredMatch=0
    i=0    
    if not os.path.exists(COVERAGE_DIR):
        os.makedirs(COVERAGE_DIR)
    while startline < avhrrIn.longitude.shape[0]:
        
        write_log("INFO","Calling get_cloudsat_avhrr_linpix: start-line = ",startline)
        tmpaid = "tmparea_%d"%i
        endline = startline+NLINES
        coverage_filename = "%s/coverage_avhrr_cloudsat_matchup_%s_%s_%s_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrrIn.satellite,orbittime[0],orbittime[1],
                             startline,endline,tmpaid)
        write_log("INFO","Coverage filename = ",coverage_filename)
        cal,cap,ok = get_cloudsat_avhrr_linpix_segment(avhrrIn,lon,lat,clTime,
                                                      (startline,endline),
                                                      SWATHWD,tmppcs,tmpaid,
                                                      coverage_filename)
        if ok:
            HasEncounteredMatch=1
            write_log("INFO","There was a match...")

        # Do not like this one /Erik    
        #if not ok and HasEncounteredMatch:
        #    write_log("INFO","Data is now empty. Leave the loop...")
        #    break
        
        if(startline==0):
            # First time:
            cloudsat_avhrr_line,cloudsat_avhrr_pixel = numpy.array(cal),numpy.array(cap)
        else:
            # Merge:
            cloudsat_avhrr_line = numpy.where(numpy.equal(cloudsat_avhrr_line,-9),cal,cloudsat_avhrr_line)
            cloudsat_avhrr_pixel = numpy.where(numpy.equal(cloudsat_avhrr_pixel,-9),cap,cloudsat_avhrr_pixel)

        startline=startline+NLINES
        i=i+1

    return cloudsat_avhrr_line,cloudsat_avhrr_pixel

# --------------------------------------------
def get_cloudsat_avhrr_linpix_segment(avhrrIn,lon,lat,cltime,lines,swath_width,tmppcs,
                                      tmpaid,covfilename):
    import numpy
    import _satproj
    import area,pcs
    import pps_gisdata
    import numpy
    
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
    areaObj = area.area(AREA)

    # Dont use thise one /Erik
    #if not os.path.exists(covfilename):
    #    write_log("INFO","Create Coverage map...")
    #    cov = _satproj.create_coverage(areaObj,lonarr,latarr,1)
    #    print covfilename
    #    writeCoverage(cov,covfilename,"satproj",AREA1KM)
    #else:
    #    write_log("INFO","Read the AVHRR-CLOUDSAT matchup coverage from file...")
    #    cov,info = readCoverage(covfilename)
    # Do like this instead
    write_log("INFO","Create Coverage map...")
    cov = _satproj.create_coverage(areaObj,lonarr,latarr,1)
    print covfilename
    writeCoverage(cov,covfilename,"satproj",AREA)
    
    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA)
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA)


    write_log("INFO","Go through cloudsat track:")
    cloudsat_avhrr_line = []
    cloudsat_avhrr_pixel = []
    cloudsat_avhrr_line_time = []
    cloudsat_avhrr_pixel_time = []
    for i in range(ndim):
        xy_tup=pps_gisdata.lonlat2xy(AREA,lon[i],lat[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
        dimx=mapped_line.shape[1]#1002#5010
        dimy=mapped_line.shape[0]#1002#5010
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            cloudsat_avhrr_line.append(mapped_line[y,x])
            cloudsat_avhrr_pixel.append(mapped_pixel[y,x])
            cloudsat_avhrr_line_time.append(-9)
            cloudsat_avhrr_pixel_time.append(-9)
        else:
            cloudsat_avhrr_line.append(-9)
            cloudsat_avhrr_pixel.append(-9)
            cloudsat_avhrr_line_time.append(-9)
            cloudsat_avhrr_pixel_time.append(-9)
            
    cloudsat_avhrr_line = numpy.array(cloudsat_avhrr_line)
    cloudsat_avhrr_pixel = numpy.array(cloudsat_avhrr_pixel)
    cloudsat_avhrr_line_time = numpy.array(cloudsat_avhrr_line_time)
    cloudsat_avhrr_pixel_time = numpy.array(cloudsat_avhrr_pixel_time)
    # Control the time diference
    match_cloudsat_points = numpy.where(numpy.not_equal(cloudsat_avhrr_line,-9))
    avhrr_time = (cloudsat_avhrr_line[match_cloudsat_points] * DSEC_PER_AVHRR_SCALINE) + avhrrIn.sec1970_start
    cl_time = cltime[match_cloudsat_points]
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

    return cloudsat_avhrr_line_time, cloudsat_avhrr_pixel_time, matchOk

# -----------------------------------------------------
def match_cloudsat_avhrr(ctypefile,cloudsatObj,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj):
    import numpy
    import time
    import string
    import numpy
    
    retv = CloudsatAvhrrTrackObject()
    
    if RESOLUTION ==1:
        lonCloudsat = cloudsatObj.longitude.ravel()
        latCloudsat = cloudsatObj.latitude.ravel()
    
    elif RESOLUTION ==5:
        lonCloudsat = cloudsatObj.longitude[:,1].ravel()
        latCloudsat = cloudsatObj.latitude[:,1].ravel()
        
    timeCloudsat = cloudsatObj.sec1970.ravel()
    ndim = lonCloudsat.shape[0]
# Add time
    if avhrrGeoObj.sec1970_end<avhrrGeoObj.sec1970_start:
        avhrr_sec1970_end = int(DSEC_PER_AVHRR_SCALINE*avhrrObj.num_of_lines+avhrrGeoObj.sec1970_start)
        avhrrGeoObj.sec1970_end = avhrr_sec1970_end
    # --------------------------------------------------------------------

    cal,cap = get_cloudsat_avhrr_linpix(avhrrGeoObj,ctypefile,lonCloudsat,latCloudsat,timeCloudsat)
    
    calnan = numpy.where(cal == NODATA, numpy.nan, cal)
    if (~numpy.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
    avhrr_lines_sec_1970 = calnan * DSEC_PER_AVHRR_SCALINE + avhrrGeoObj.sec1970_start

    # Find all matching Cloudsat pixels within +/- sec_timeThr from the AVHRR data
    idx_match = elements_within_range(cloudsatObj.sec1970, avhrr_lines_sec_1970, sec_timeThr)
    #            numpy.logical_and(cloudsatObj.sec1970 > avhrr_lines_sec_1970 - sec_timeThr,
    #                              cloudsatObj.sec1970 < avhrr_lines_sec_1970 + sec_timeThr)
    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)

    lon_cloudsat = numpy.repeat(lonCloudsat, idx_match)
    lat_cloudsat = numpy.repeat(latCloudsat, idx_match)

    # Cloudsat line,pixel inside AVHRR swath:
    cal_on_avhrr = numpy.repeat(cal, idx_match)
    cap_on_avhrr = numpy.repeat(cap, idx_match)

    print "Make CPR echo top array..."
    retv.cloudsat.echo_top = numpy.repeat(cloudsatObj.CPR_Echo_Top,idx_match).astype('d')

    retv.cloudsat.sec_1970 = numpy.repeat(cloudsatObj.sec1970,idx_match)

    retv.cloudsat.latitude = numpy.repeat(latCloudsat,idx_match)
    retv.cloudsat.longitude = numpy.repeat(lonCloudsat,idx_match)

    retv.cloudsat.MODIS_Cloud_Fraction = numpy.repeat(cloudsatObj.MODIS_Cloud_Fraction,idx_match)
    retv.cloudsat.MODIS_cloud_flag = numpy.repeat(cloudsatObj.MODIS_cloud_flag,idx_match)
    print "cap_on_avhrr.shape: ",cap_on_avhrr.shape
    retv.cloudsat.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.cloudsat.avhrr_pixnum = cap_on_avhrr.astype('i')
    
    missing_data = -9.9
    col_dim = cloudsatObj.CPR_Cloud_mask.shape[1]
    x_ccm = numpy.repeat(cloudsatObj.CPR_Cloud_mask[::,0],idx_match)
    x_rr = numpy.repeat(cloudsatObj.Radar_Reflectivity[::,0],idx_match)
    x_H = numpy.repeat(cloudsatObj.Height[::,0],idx_match)
    
    for i in range(1,col_dim):
        x_ccm = numpy.concatenate((x_ccm,numpy.repeat(cloudsatObj.CPR_Cloud_mask[::,i],idx_match)))
        x_rr = numpy.concatenate((x_rr,numpy.repeat(cloudsatObj.Radar_Reflectivity[::,i],idx_match)))
        x_H = numpy.concatenate((x_H,numpy.repeat(cloudsatObj.Height[::,i],idx_match)))
    
    N_ccm = x_ccm.shape[0]/col_dim
    retv.cloudsat.cloud_mask = numpy.reshape(x_ccm,(col_dim,N_ccm)).astype('d')    
    N_rr = x_rr.shape[0]/col_dim
    RadarRefl = numpy.reshape(x_rr,(col_dim,N_rr))
    RadarRefl = numpy.where(numpy.less(RadarRefl,0),missing_data,RadarRefl)
    retv.cloudsat.Radar_Reflectivity = RadarRefl.astype('d')
    N_H = x_H.shape[0]/col_dim
    retv.cloudsat.Height = numpy.reshape(x_H,(col_dim,N_H)).astype('d')

    # One-d arrays:
    retv.cloudsat.SurfaceHeightBin = \
                       numpy.repeat(cloudsatObj.SurfaceHeightBin.ravel(),
                                      idx_match.ravel()).astype('d')
    retv.cloudsat.SurfaceHeightBin_fraction = \
                       numpy.repeat(cloudsatObj.SurfaceHeightBin_fraction.ravel(),
                                      idx_match.ravel()).astype('d')

    print "Cloudsat observation time of first cloudsat-avhrr match: ",\
          time.gmtime(retv.cloudsat.sec_1970[0])
    print "Cloudsat observation time of last cloudsat-avhrr match: ",\
          time.gmtime(retv.cloudsat.sec_1970[N_H-1])
    
    # Elevation is given in km's. Convert to meters:
    retv.cloudsat.elevation = numpy.repeat(\
        cloudsatObj.elevation.ravel(),idx_match.ravel()).astype('d')

    retv.avhrr.sec_1970 = numpy.add(avhrrGeoObj.sec1970_start,
                                      cal_on_avhrr * DSEC_PER_AVHRR_SCALINE)
    retv.diff_sec_1970 = retv.cloudsat.sec_1970 - retv.avhrr.sec_1970
    min_diff = numpy.minimum.reduce(retv.diff_sec_1970)
    max_diff = numpy.maximum.reduce(retv.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-cloudsat): ",\
          numpy.maximum.reduce(retv.diff_sec_1970),numpy.minimum.reduce(retv.diff_sec_1970)

    print "AVHRR observation time of first cloudsat-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[0])
    print "AVHRR observation time of last cloudsat-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[N_H-1])

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
    satz_track = []
    # maby skould take consideration of missing data also, as for satz.
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
        if avhrrAngObj.satz.data[cal_on_avhrr[i],cap_on_avhrr[i]] == avhrrAngObj.satz.no_data or avhrrAngObj.satz.data[cal_on_avhrr[i],cap_on_avhrr[i]] == avhrrAngObj.satz.missing_data:
            ang = -9
        else:
            ang = avhrrAngObj.satz.data[cal_on_avhrr[i],cap_on_avhrr[i]] * avhrrAngObj.satz.gain + avhrrAngObj.satz.intercept
        satz_track.append(ang)
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

    retv.avhrr.latitude = numpy.array(lat_avhrr_track)
    retv.avhrr.longitude = numpy.array(lon_avhrr_track)
    retv.avhrr.cloudtype = numpy.array(ctype_track)
    retv.avhrr.surftemp = numpy.array(surft_track)
    retv.avhrr.bt11micron = numpy.array(bt11micron_track)
    retv.avhrr.bt12micron = numpy.array(bt12micron_track)
    retv.avhrr.satz = numpy.array(satz_track)
    if ctth:
        retv.avhrr.ctth_height = numpy.array(ctth_height_track)
        retv.avhrr.ctth_pressure = numpy.array(ctth_pressure_track)
        retv.avhrr.ctth_temperature = numpy.array(ctth_temperature_track)

    print "AVHRR-PPS Cloud Type,latitude: shapes = ",\
          retv.avhrr.cloudtype.shape,retv.avhrr.latitude.shape

    ll = []
    for i in range(ndim):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],idx_match[i])))
    basename = os.path.basename(ctypefile).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_")
    
    datapath = "%s/%s/%skm/%s/%s/%s" %(DATA_DIR, base_sat, RESOLUTION, base_year, base_month, AREA)
    if not os.path.exists(datapath):
        os.makedirs(datapath)   
    fd = open("%s/%skm_%s_cloudtype_cloudsat-GEOPROF_track2.txt"%(datapath, RESOLUTION, basename),"w")
    fd.writelines(ll)
    fd.close()
    ll = []
    for i in range(N_H):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_cloudsat[i],lat_cloudsat[i],0)))
    fd = open("%s/%skm_%s_cloudtype_cloudsat-GEOPROF_track_excl.txt"%(datapath, RESOLUTION, basename),"w")
    fd.writelines(ll)
    fd.close()
    
    return retv,min_diff,max_diff

#------------------------------------------------------------------------------------------

def reshapeCloudsat(cloudsatfiles,avhrr):
    import time
    import numpy
    import sys
    
    clsat = CloudsatObject()
    if avhrr.sec1970_end<avhrr.sec1970_start:
        avhrr_end = int(DSEC_PER_AVHRR_SCALINE*avhrr.num_of_lines+avhrr.sec1970_start)
    else:
        avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    
    startCloudsat = get_cloudsat(cloudsatfiles[0])
    startCloudsat.Profile_time = numpy.add(startCloudsat.Profile_time,startCloudsat.TAI_start)
    
    for i in range(len(cloudsatfiles)-1):
        newCloudsat = get_cloudsat(cloudsatfiles[i+1])
        newCloudsat.Profile_time = numpy.add(newCloudsat.Profile_time,newCloudsat.TAI_start)
        
        clsat_start_all = startCloudsat.sec1970.ravel()
        clsat_new_all = newCloudsat.sec1970.ravel()
        
        if not clsat_start_all[0]<clsat_new_all[0]:
            print "cloudsat files are in the wrong order"
            print("Program cloudsat.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9)
        clsat_break = numpy.argmin(numpy.abs(clsat_start_all - clsat_new_all[0]))+1
        
        startCloudsat.sec1970 = numpy.concatenate((startCloudsat.sec1970[0:clsat_break],newCloudsat.sec1970))
        if RESOLUTION ==1:
            startCloudsat.longitude = numpy.concatenate((startCloudsat.longitude[0:clsat_break],newCloudsat.longitude))
            startCloudsat.latitude = numpy.concatenate((startCloudsat.latitude[0:clsat_break],newCloudsat.latitude))
            startCloudsat.elevation = numpy.concatenate((startCloudsat.elevation[0:clsat_break],newCloudsat.elevation))
            startCloudsat.Profile_time = numpy.concatenate((startCloudsat.Profile_time[0:clsat_break],newCloudsat.Profile_time))
            startCloudsat.CPR_Echo_Top = numpy.concatenate((startCloudsat.CPR_Echo_Top[0:clsat_break],newCloudsat.CPR_Echo_Top))
            startCloudsat.Clutter_reduction_flag = numpy.concatenate((startCloudsat.Clutter_reduction_flag[0:clsat_break],newCloudsat.Clutter_reduction_flag))
            startCloudsat.Data_quality = numpy.concatenate((startCloudsat.Data_quality[0:clsat_break],newCloudsat.Data_quality))
            startCloudsat.CPR_Cloud_mask = numpy.concatenate((startCloudsat.CPR_Cloud_mask[0:clsat_break,:],newCloudsat.CPR_Cloud_mask))        
            startCloudsat.Data_targetID = numpy.concatenate((startCloudsat.Data_targetID[0:clsat_break],newCloudsat.Data_targetID))
            startCloudsat.Gaseous_Attenuation = numpy.concatenate((startCloudsat.Gaseous_Attenuation[0:clsat_break,:],newCloudsat.Gaseous_Attenuation))
            startCloudsat.MODIS_Cloud_Fraction = numpy.concatenate((startCloudsat.MODIS_Cloud_Fraction[0:clsat_break],newCloudsat.MODIS_Cloud_Fraction))
            startCloudsat.MODIS_cloud_flag = numpy.concatenate((startCloudsat.MODIS_cloud_flag[0:clsat_break],newCloudsat.MODIS_cloud_flag))
            startCloudsat.Radar_Reflectivity = numpy.concatenate((startCloudsat.Radar_Reflectivity[0:clsat_break,:],newCloudsat.Radar_Reflectivity))
            startCloudsat.Height = numpy.concatenate((startCloudsat.Height[0:clsat_break,:],newCloudsat.Height))
            startCloudsat.SigmaZero = numpy.concatenate((startCloudsat.SigmaZero[0:clsat_break],newCloudsat.SigmaZero))
            startCloudsat.SurfaceHeightBin = numpy.concatenate((startCloudsat.SurfaceHeightBin[0:clsat_break],newCloudsat.SurfaceHeightBin))
            startCloudsat.SurfaceHeightBin_fraction = numpy.concatenate((startCloudsat.SurfaceHeightBin_fraction[0:clsat_break],newCloudsat.SurfaceHeightBin_fraction))
            startCloudsat.sem_NoiseFloor = numpy.concatenate((startCloudsat.sem_NoiseFloor[0:clsat_break],newCloudsat.sem_NoiseFloor))        
            startCloudsat.sem_NoiseFloorVar = numpy.concatenate((startCloudsat.sem_NoiseFloorVar[0:clsat_break],newCloudsat.sem_NoiseFloorVar))
            startCloudsat.sem_NoiseGate = numpy.concatenate((startCloudsat.sem_NoiseGate[0:clsat_break],newCloudsat.sem_NoiseGate))
            start_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_start - sec_timeThr))))-1 # Minus one to get one extra, just to be certain
            end_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain
    
            clsat.sec1970 = startCloudsat.sec1970[start_break:end_break] 
            clsat.longitude = startCloudsat.longitude[start_break:end_break] 
            clsat.latitude = startCloudsat.latitude[start_break:end_break] 
            clsat.elevation = startCloudsat.elevation[start_break:end_break] 
            clsat.Profile_time = startCloudsat.Profile_time[start_break:end_break] 
            clsat.CPR_Cloud_mask = startCloudsat.CPR_Cloud_mask[start_break:end_break,:] 
            clsat.CPR_Echo_Top = startCloudsat.CPR_Echo_Top[start_break:end_break] 
            clsat.Clutter_reduction_flag = startCloudsat.Clutter_reduction_flag[start_break:end_break]  
            clsat.Data_quality = startCloudsat.Data_quality[start_break:end_break] 
            clsat.Data_targetID = startCloudsat.Data_targetID[start_break:end_break] 
            clsat.Gaseous_Attenuation = startCloudsat.Gaseous_Attenuation[start_break:end_break,:] 
            clsat.MODIS_Cloud_Fraction = startCloudsat.MODIS_Cloud_Fraction[start_break:end_break] 
            clsat.MODIS_cloud_flag = startCloudsat.MODIS_cloud_flag[start_break:end_break] 
            clsat.Radar_Reflectivity = startCloudsat.Radar_Reflectivity[start_break:end_break,:] 
            clsat.Height = startCloudsat.Height[start_break:end_break,:] 
            clsat.SigmaZero = startCloudsat.SigmaZero[start_break:end_break] 
            clsat.SurfaceHeightBin = startCloudsat.SurfaceHeightBin[start_break:end_break] 
            clsat.SurfaceHeightBin_fraction = startCloudsat.SurfaceHeightBin_fraction[start_break:end_break] 
            clsat.sem_NoiseFloor = startCloudsat.sem_NoiseFloor[start_break:end_break] 
            clsat.sem_NoiseFloorVar = startCloudsat.sem_NoiseFloorVar[start_break:end_break] 
            clsat.sem_NoiseGate = startCloudsat.sem_NoiseGate[start_break:end_break] 
            if clsat.Profile_time.shape[0] <= 0:
                print("No time match, please try with some other CloudSat files")
                print("Program cloudsat.py at line %i" %(inspect.currentframe().f_lineno+1))
                sys.exit(-9)
            else:
                clsat.TAI_start = clsat.Profile_time[0]
                clsat.Profile_time = clsat.Profile_time - clsat.TAI_start
        
        elif RESOLUTION ==5:
            startCloudsat.longitude = numpy.concatenate((startCloudsat.longitude[0:clsat_break,:],newCloudsat.longitude))
            startCloudsat.latitude = numpy.concatenate((startCloudsat.latitude[0:clsat_break,:],newCloudsat.latitude))
            startCloudsat.elevation = numpy.concatenate((startCloudsat.elevation[0:clsat_break,:],newCloudsat.elevation)) 
            startCloudsat.Profile_time = numpy.concatenate((startCloudsat.Profile_time[0:clsat_break,:],newCloudsat.Profile_time))  
            startCloudsat.CPR_Echo_Top = numpy.concatenate((startCloudsat.CPR_Echo_Top[0:clsat_break,:],newCloudsat.CPR_Echo_Top))
            startCloudsat.Clutter_reduction_flag = numpy.concatenate((startCloudsat.Clutter_reduction_flag[0:clsat_break,:],newCloudsat.Clutter_reduction_flag))
            startCloudsat.Data_quality = numpy.concatenate((startCloudsat.Data_quality[0:clsat_break,:],newCloudsat.Data_quality))
            startCloudsat.CPR_Cloud_mask = numpy.concatenate((startCloudsat.CPR_Cloud_mask[0:clsat_break,:],newCloudsat.CPR_Cloud_mask))
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
                print("Program cloudsat.py at line %i" %(inspect.currentframe().f_lineno+1))
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
    cloudsat = get_cloudsat(cloudsatfile)
    #cloudsat = read_cloudsat(cloudsatfile)

    lonCloudsat = cloudsat.longitude.ravel()
    latCloudsat = cloudsat.latitude.ravel()

    # Test:
    ndim = lonCloudsat.ravel().shape[0]
    idx_match = numpy.zeros((ndim,),'b')
    idx_match[0:10] = 1

    x = numpy.repeat(cloudsat.Height[::,0],idx_match)
    for i in range(1,125):
        x = numpy.concatenate((x,numpy.repeat(cloudsat.Height[::,i],idx_match)))
    N = x.shape[0]/125
    cloudsat.Height = numpy.reshape(x,(125,N))
