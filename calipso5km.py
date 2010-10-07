
import pdb
import inspect

from pps_basic_configure import *
from pps_error_messages import *

#DSEC_PER_AVHRR_SCALINE = 0.1667 # Full scan period, i.e. the time interval between two consecutive lines (sec)
DSEC_PER_AVHRR_SCALINE5KM = 1.0/6*4 # A "work for the time being" solution
from cloudsat_calipso_avhrr_match import AREA5KM, SUB_DIR, DATA_DIR, sec_timeThr

#AREA = "cea5km_test" #
#AREA = "arctic_super_1002_5km"
#MAIN_DIR = "/data/proj/safworks/adam/calipso_data"
#MAIN_DIR = "/data/proj_nsc1/safworks/calipso_cloudsat/data/arctic/"
#MAIN_DIR = "/data/proj_nsc1/safworks/kgkarl/ORR-B-datasets/calipso5km/"
#SUB_DIR = "metop02_calipso_2007spring"
#SUB_DIR = "matchups/"
#SUB_DIR = "noaa18_calipso_cloudsat_2007DEC_KG"
#SUB_DIR = "noaa17_calipso_cloudsat_2007JUN_KG"
#SUB_DIR = "noaa18_calipso_cloudsat_2007JUN_KG"

#SATPOS_DIR = "%s%s"%(MAIN_DIR,SUB_DIR)
#EPHE_DIR = "%s%s"%(MAIN_DIR,SUB_DIR)
COMPRESS_LVL = 6
COVERAGE_DIR = "%s/5km/%s"%(SUB_DIR,AREA5KM)

#NLINES=1000
NLINES=6000
#SWATHWD=2048
SWATHWD5km=409
NODATA=-9
        
class ppsAvhrr5kmObject:
    def __init__(self):
        self.longitude=None
        self.latitude=None
        self.sec_1970=None
        self.ctth_height=None
        self.ctth_pressure=None
        self.ctth_temperature=None
        self.cloudtype=None
        self.surftemp=None
        self.bt11micron=None
        self.bt12micron=None  
        self.satz=None
         
class Calipso5kmObject:
    def __init__(self):
        self.longitude=None
        self.latitude=None
        self.avhrr_linnum=None
        self.avhrr_pixnum=None
        self.cloud_fraction=None
        self.cloud_top_profile=None
        self.cloud_base_profile=None
        self.cloud_mid_temperature=None
        self.number_of_layers_found=None
        self.igbp=None
        self.nsidc=None
        self.elevation=None
        self.time=None
        self.utc_time=None 
        self.sec_1970=None
        self.feature_classification_flags=None
        self.day_night_flag=None
        self.optical_depth=None
        self.optical_depth_uncertainty=None
#class CloudSatStartEndObject:
#    def __init__(self):
#        self.CloudSatStart = None
#        self.CloudSatEnd = None
        
class Calipso5kmAvhrr5kmTrackObject:
    def __init__(self):
        self.avhrr5km=ppsAvhrr5kmObject()
        self.calipso5km=Calipso5kmObject()
        self.diff_sec_1970=None


class area_interface:
    pass

# ----------------------------------------
def readCaliop5kmAvhrr5kmMatchObj(filename5km):
    import _pyhl

    retv5km = Calipso5kmAvhrr5kmTrackObject()

    a=_pyhl.read_nodelist(filename5km)
    a.selectAll()
    a.fetch()

    # Match-Up - time difference:
    retv5km.diff_sec_1970 = a.getNode("/diff_sec_1970").data()

    # Calipso:
    retv5km.calipso5km.longitude = a.getNode("/calipso/longitude").data()
    retv5km.calipso5km.latitude = a.getNode("/calipso/latitude").data()
    retv5km.calipso5km.avhrr_linnum = a.getNode("/calipso/avhrr_linnum").data()
    retv5km.calipso5km.avhrr_pixnum = a.getNode("/calipso/avhrr_pixnum").data()
    
    retv5km.calipso5km.cloud_fraction = a.getNode("/calipso/cloud_fraction").data()
    retv5km.calipso5km.cloud_top_profile = a.getNode("/calipso/cloud_top_profile").data()
    retv5km.calipso5km.cloud_base_profile = a.getNode("/calipso/cloud_base_profile").data()
    retv5km.calipso5km.cloud_mid_temperature = a.getNode("/calipso/cloud_mid_temperature").data()
    retv5km.calipso5km.elevation = a.getNode("/calipso/elevation").data()
    retv5km.calipso5km.number_of_layers_found = a.getNode("/calipso/number_of_layers_found").data()
    retv5km.calipso5km.optical_depth = a.getNode("/calipso/optical_depth").data()
    retv5km.calipso5km.optical_depth_uncertainty = a.getNode("/calipso/optical_depth_uncertainty").data()
    try:
        retv5km.calipso5km.feature_classification_flags = a.getNode("/calipso/feature_classification_flags").data()
    except:
        print "No feature_classification_flags array in file!"
        pass
    
    retv5km.calipso5km.igbp = a.getNode("/calipso/igbp").data()
    retv5km.calipso5km.nsidc = a.getNode("/calipso/nsidc").data()
    #retv5km.calipso5km.utc_time = a.getNode("/calipso/utc_time").data()
    retv5km.calipso5km.sec_1970 = a.getNode("/calipso/sec_1970").data()

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
def writeCaliop5kmAvhrr5kmMatchObj(filename5km,ca_obj5km,compress_lvl):
    import _pyhl
    status = -1
    
    a=_pyhl.nodelist()

    shape = [ca_obj5km.calipso5km.longitude.shape[0]]

    # Match-Up - time difference:
    # ====
    b=_pyhl.node(_pyhl.DATASET_ID,"/diff_sec_1970")
    b.setArrayValue(1,shape,ca_obj5km.diff_sec_1970,"double",-1)
    a.addNode(b)

    # Calipso
    # ====
    b=_pyhl.node(_pyhl.GROUP_ID,"/calipso")
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/longitude")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.longitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/latitude")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.latitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/avhrr_linnum")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.avhrr_linnum,"int",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/avhrr_pixnum")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.avhrr_pixnum,"int",-1)
    a.addNode(b)

    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_fraction")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.cloud_fraction,"double",-1)
    a.addNode(b)

    shape2d = ca_obj5km.calipso5km.cloud_top_profile.shape
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_top_profile")
    b.setArrayValue(1,shape2d,ca_obj5km.calipso5km.cloud_top_profile,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_base_profile")
    b.setArrayValue(1,shape2d,ca_obj5km.calipso5km.cloud_base_profile,"double",-1)
    a.addNode(b)

    #b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_mid_profile")
    #b.setArrayValue(1,shape2d,ca_obj5km.calipso5km.cloud_mid_profile,"double",-1)
    #a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_mid_temperature")
    b.setArrayValue(1,shape2d,ca_obj5km.calipso5km.cloud_mid_temperature,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/number_of_layers_found")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.number_of_layers_found,"int",-1)
    a.addNode(b)    
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/igbp")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.igbp,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/nsidc")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.nsidc,"uchar",-1)
    a.addNode(b)        
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/optical_depth")
    b.setArrayValue(1,shape2d,ca_obj5km.calipso5km.optical_depth,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/optical_depth_uncertainty")
    b.setArrayValue(1,shape2d,ca_obj5km.calipso5km.optical_depth_uncertainty,"double",-1)
    a.addNode(b)
       
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/elevation")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.elevation,"double",-1)
    a.addNode(b)
    #b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/time")
    #b.setArrayValue(1,shape,ca_obj5km.calipso5km.time,"double",-1)
    #a.addNode(b)
    #b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/utc_time")
    #b.setArrayValue(1,shape,ca_obj5km.calipso5km.utc_time,"double",-1)
    #a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/sec_1970")
    b.setArrayValue(1,shape,ca_obj5km.calipso5km.sec_1970,"double",-1)
    a.addNode(b)
    try:
        b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/feature_classification_flags")

        print ca_obj5km.calipso5km.feature_classification_flags.shape
        #print ca_obj5km.calipso5km.feature_classification_flags.typecode()
        #print ca_obj5km.calipso5km.feature_classification_flags.itemsize()
        b.setArrayValue(1,shape2d,ca_obj5km.calipso5km.feature_classification_flags,"int",-1)
        a.addNode(b)
    except:
        pass

    # AVHRR
    # ====
    b=_pyhl.node(_pyhl.GROUP_ID,"/avhrr")
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/longitude")
    shape = ca_obj5km.avhrr5km.longitude.shape
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.longitude,"float",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/latitude")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.latitude,"float",-1)
    a.addNode(b)

    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/sec_1970")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.sec_1970,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_pressure")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.ctth_pressure,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_temperature")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.ctth_temperature,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/ctth_height")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.ctth_height,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/cloudtype")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.cloudtype,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/bt11micron")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.bt11micron,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/bt12micron")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.bt12micron,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/surftemp")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.surftemp,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/satz")
    b.setArrayValue(1,shape,ca_obj5km.avhrr5km.satz,"double",-1)
    a.addNode(b)

    status = a.write(filename5km,compress_lvl)
    return status

# ----------------------------------------
class SatProjCov5km:
    def __init__(self):
        self.coverage=None
        self.colidx=None
        self.rowidx=None    

# ------------------------------------------------------------------
def writeCoverage5km(covIn,filename5km,inAid,outAid):
    import _pyhl

    a=_pyhl.nodelist()

    b=_pyhl.node(_pyhl.GROUP_ID,"/info")
    a.addNode(b)
    b=_pyhl.node(_pyhl.ATTRIBUTE_ID,"/info/description")
    b.setScalarValue(-1,"MSG coverage from area %s on to area %s"%(inAid,outAid),"string",-1)
    a.addNode(b)
    
    shape=[covIn.coverage.shape[0],covIn.coverage.shape[1]]
    b=_pyhl.node(_pyhl.DATASET_ID,"/coverage")
    b.setArrayValue(1,shape,covIn.coverage,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/rowidx")
    b.setArrayValue(1,shape,covIn.rowidx,"ushort",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/colidx")
    b.setArrayValue(1,shape,covIn.colidx,"ushort",-1)
    a.addNode(b)

    a.write(filename5km,COMPRESS_LVL)
    
    return

# ------------------------------------------------------------------
def readCoverage5km(filename5km):
    import _pyhl

    a=_pyhl.read_nodelist(filename5km)
    b=a.getNodeNames()
    
    a.selectNode("/info/description")
    a.selectNode("/coverage")
    a.selectNode("/rowidx")
    a.selectNode("/colidx")
    a.fetch()

    info={}
    c=a.getNode("/info/description");
    d=c.data()
    info["description"]=d

    c=a.getNode("/coverage")
    coverage=c.data()
    c=a.getNode("/rowidx")
    rowidx=c.data()
    c=a.getNode("/colidx")
    colidx=c.data()
    
    retv5km = SatProjCov5km()
    retv5km.coverage = coverage.astype('Int8')
    retv5km.rowidx = rowidx.astype('Int16')
    retv5km.colidx = colidx.astype('Int16')

    return retv5km,info
# ------------------------------------------------------------------    
#def writeCloudSatStartEnd(filename5km,StartEnd_Obj):
#    import _pyhl
#    import pdb
#    a=_pyhl.nodelist()
#    
#    shape = [StartEnd_Obj.CloudSatStart.shape[0]]
#    b=_pyhl.node(_pyhl.DATASET_ID,"/StartTime")
#    ##########################################################################################################################################
#    #pdb.set_trace()
#    b.setArrayValue(1,shape,StartEnd_Obj.CloudSatStart,"double",-1)
#    a.addNode(b)
#    b=_pyhl.node(_pyhl.DATASET_ID,"/EndTime")
#    b.setArrayValue(1,shape,StartEnd_Obj.CloudSatEnd,"double",-1)
#    a.addNode(b)
#    a.write(filename5km,COMPRESS_LVL)
#    
#    return
# ------------------------------------------------------------------
#def readCloudSatStartEnd(filename5km):
#    import _pyhl
#
#    retv5km = CloudSatStartEndObject()
#
#    a=_pyhl.read_nodelist(filename5km)
#    a.selectAll()
#    a.fetch()
#    
#    retv5km.CloudSatStart = a.getNode("/StartTime").data()
#    retv5km.CloudSatEnd = a.getNode("/EndTime").data()
#    
#    return retv5km
# -----------------------------------------------------
def define_pcs5km(id,name,definition):
    import pcs
    p = pcs.usgs(name,definition)
    pcs.register(id,p)

# -----------------------------------------------------
def define_longlat_ll5km(id, name, pcs_id, ll_ll, size, scale):
    import pcs
    import area
    a = area_interface()
    a.name = name
    a.pcs = pcs.pcs(pcs_id)
    x, y = a.pcs.proj(ll_ll)
    a.extent = (x, y, x + scale * size[0], y + scale * size[1])
    a.xsize = size[0]
    a.ysize = size[1]
    area.register(id, a)


# -----------------------------------------------------
def sec1970_to_julianday5km(sec1970):
    import pps_time_util
    import time
    
    year,month,day,hour,minutes,sec,dummy,dummy,dummy = time.gmtime(sec1970)
    jday_1950 = int(pps_time_util.getJulianDay(year,month,day) - pps_time_util.getJulianDay(1950,1,1))
    jday = jday_1950 + (hour+minutes/60.0+sec/3600)/24.0

    return jday

# -----------------------------------------------------
def avhrr_linepix_from_lonlat_aapp5km(lon,lat,avhrrObj,platform,norbit,yyyymmdd):
    import numpy.oldnumeric as Numeric
    import CreateAngles
    import _py_linepix_lonlat
    
    ndim = lon.shape[0]
    lin = Numeric.zeros((ndim,),'d')
    pix = Numeric.zeros((ndim,),'d')

    if platform.find("metop") >= 0:
        file_satpos = "%s/satpos_M%.2d_%s.txt"%(SATPOS_DIR,string.atoi(platform.split("metop")[1]),yyyymmdd)
    else:
        file_satpos = "%s/satpos_%s_%s.txt"%(SATPOS_DIR,platform,yyyymmdd)

    file_ephe = "%s/ephe_%s.txt"%(EPHE_DIR,yyyymmdd)
    print file_ephe
    
    start_jday,end_jday,sec1970_start,sec1970_end = CreateAngles.get_acqtime(file_ephe,norbit)
    write_log("INFO","From ephemeris file: platform,norbit,start_jday,end_jday = ",platform,norbit,start_jday,end_jday)

    start_jday = sec1970_to_julianday(avhrrObj.sec1970_start)
    end_jday   = sec1970_to_julianday(avhrrObj.sec1970_end)
    write_log("INFO","From AVHRR file: platform,norbit,start_jday,end_jday = ",platform,norbit,start_jday,end_jday)

    #attitude = CreateAngles.get_attitude(platform,norbit,"tle")
    #attitude = avhrrObj.attitude_error["yaw"],avhrrObj.attitude_error["roll"],avhrrObj.attitude_error["pitch"]
    attitude = (0.,0.,0.)
    write_log("INFO","YAW,ROLL,PITCH:",attitude[0],attitude[1],attitude[2])

    start_end_times = (start_jday,end_jday)

    # Get the avhrr line,pixel arrays matching input lon,lat.
    # Those points where the time of the avhrr pixel is too far away from the time of the lon,lat arrays
    # will be set to NODATA (line=-1,pixel=-1):
    # (Not yet implemented the time constraints)
    this = _py_linepix_lonlat.getLinePixFromLonLat(platform,file_satpos,lon,lat,lin,pix,
                                                   attitude,start_end_times)

    return lin,pix

# --------------------------------------------
def getBoundingBox(lon,lat):
    import numpy.oldnumeric as Numeric
    maxlon = Numeric.maximum.reduce(lon.ravel())
    minlon = Numeric.minimum.reduce(lon.ravel())
    maxlat = Numeric.maximum.reduce(lat.ravel())
    minlat = Numeric.minimum.reduce(lat.ravel())

    return minlon,minlat,maxlon,maxlat

# --------------------------------------------
def get_calipso5km_avhrr5km_linpix(avhrrIn,avhrrname,lon,lat,caTime):
    #import numpy.oldnumeric as Numeric
    import numpy
    tmppcs="tmpproj"
    define_pcs5km(tmppcs, "Plate Caree, central meridian at 15E",
               ['proj=eqc','ellps=bessel', 'lon_0=15'])


    startline=0

    """
    # Test code:
    tmpaid="tmparea"
    mask = get_calipso_avhrr_linpix(avhrr,lon,lat,(startline,startline+NLINES),SWATHWD5km,tmppcs,tmpaid)
    import Image
    dimy,dimx = mask.shape
    that = Image.fromstring("L",(dimx,dimy),mask.tostring())
    that.save("./yt.png")
    that.thumbnail((dimx/8,dimy/8))
    that.save("./yt_thumbnail.png")
    """
    orbittime =  os.path.basename(avhrrname).split("_")[1:3]
    
    if not os.path.exists(COVERAGE_DIR):
        os.makedirs(COVERAGE_DIR)
    Inside=0
    HasEncounteredMatch=0
    i=0
    while startline < avhrrIn.longitude.shape[0]:
        write_log("INFO","Calling get_calipso_avhrr_linpix: start-line = ",startline)
        tmpaid = "tmparea_%d"%i
        endline = startline+NLINES
        coverage_filename = "%s/coverage_avhrr_caliop_matchup_%s_%s_%s_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrrIn.satellite,orbittime[0],orbittime[1],
                             startline,endline,tmpaid)
        
        write_log("INFO","Coverage filename = ",coverage_filename)
        cal,cap,ok = get_calipso5km_avhrr5km_linpix_segment(avhrrIn,lon,lat,caTime,
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
            calipso_avhrr_line,calipso_avhrr_pixel = numpy.array(cal),numpy.array(cap)
        else:
            # Merge:
            calipso_avhrr_line = numpy.where(numpy.equal(calipso_avhrr_line,-9),cal,calipso_avhrr_line)
            calipso_avhrr_pixel = numpy.where(numpy.equal(calipso_avhrr_pixel,-9),cap,calipso_avhrr_pixel)

        startline=startline+NLINES
        i=i+1

    return calipso_avhrr_line,calipso_avhrr_pixel

# --------------------------------------------
def get_calipso5km_avhrr5km_linpix_segment(avhrrIn,lon,lat,time5km,lines,swath_width,tmppcs,
                                     tmpaid,covfilename5km):
    #import numpy.oldnumeric as Numeric
    import _satproj
    import area,pcs
    import pps_gisdata
    import numpy
    import time
    ndim = lon.shape[0]
    
    if avhrrIn.longitude.shape[0] > lines[1]:
        lines_end = lines[1]
    else:
        lines_end = avhrrIn.longitude.shape[0]
    lines_start = lines[0]

    write_log("INFO","lines_start,lines_end: ",lines_start,lines_end)
    nlines = lines_end - lines_start # Calculates number of lines
    lonarr = avhrrIn.longitude[lines_start:lines_end,::]    # Get the concerned longitude points
    latarr = avhrrIn.latitude[lines_start:lines_end,::]     # Get the concerned latitude points

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
    """

    areaObj = area.area(AREA5KM)#("arctic_super_5010")
 
#    if not os.path.exists(covfilename5km):
#        write_log("INFO","Create Coverage map...")
#        cov = _satproj.create_coverage(areaObj,lonarr,latarr,1)
#        print covfilename5km
#        writeCoverage5km(cov,covfilename5km,"satproj",AREA)
#    else:
#        write_log("INFO","Read the AVHRR-CALIOP matchup coverage from file...")
#        cov,info = readCoverage5km(covfilename5km)
    
    write_log("INFO","Create Coverage map...")
    cov = _satproj.create_coverage(areaObj,lonarr,latarr,0)
    print covfilename5km
    writeCoverage5km(cov,covfilename5km,"satproj",AREA5KM)

    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA)
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA)


    write_log("INFO","Go through calipso track:")
    calipso_avhrr_line = []
    calipso_avhrr_pixel = []
    calipso_avhrr_line_time = []
    calipso_avhrr_pixel_time = []
    for i in range(ndim):
        xy_tup=pps_gisdata.lonlat2xy(AREA5KM,lon[i],lat[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
        
        dimx=mapped_line.shape[1]#1002#5010
        dimy=mapped_line.shape[0]#1002#5010
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            calipso_avhrr_line.append(mapped_line[y,x])
            calipso_avhrr_pixel.append(mapped_pixel[y,x])
            calipso_avhrr_line_time.append(-9)
            calipso_avhrr_pixel_time.append(-9)
        else:
            calipso_avhrr_line.append(-9)
            calipso_avhrr_pixel.append(-9)
            calipso_avhrr_line_time.append(-9)
            calipso_avhrr_pixel_time.append(-9)

    calipso_avhrr_line = numpy.array(calipso_avhrr_line)
    calipso_avhrr_pixel = numpy.array(calipso_avhrr_pixel)
    calipso_avhrr_line_time = numpy.array(calipso_avhrr_line_time)
    calipso_avhrr_pixel_time = numpy.array(calipso_avhrr_pixel_time)
    
    # Control the time diference
    match_calipso_points = numpy.where(numpy.not_equal(calipso_avhrr_line,-9))
    avhrr_time = (calipso_avhrr_line[match_calipso_points] * DSEC_PER_AVHRR_SCALINE5KM) + avhrrIn.sec1970_start
    cal_time = time5km[match_calipso_points]
    time_diff = avhrr_time-cal_time
    #max_time_diff_allowed = 50*60 #Based on that a lap is 102 min
    max_time_diff_allowed = sec_timeThr
    time_match = numpy.where(abs(time_diff)<max_time_diff_allowed)
    if time_match[0].shape[0]==0:             
        x=numpy.repeat(calipso_avhrr_line_time,numpy.not_equal(calipso_avhrr_line_time,-9))
    else:
        calipso_avhrr_line_time[match_calipso_points[0][time_match]] = calipso_avhrr_line[match_calipso_points[0][time_match]]
        calipso_avhrr_pixel_time[match_calipso_points[0][time_match]] = calipso_avhrr_pixel[match_calipso_points[0][time_match]]
        x=numpy.repeat(calipso_avhrr_line_time,numpy.not_equal(calipso_avhrr_line_time,-9))
    
    write_log("INFO","Number of matching points = ",x.shape[0])
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0

    return calipso_avhrr_line_time,calipso_avhrr_pixel_time,matchOk
                    
# -----------------------------------------------------
def match_calipso5km_avhrr5km(ctypefile5km,calipso5kmObj,avhrr5kmGeoObj,avhrr5kmObj,ctype5km,ctth5km,surft5km,avhrr5kmAngObj):
    #import numpy.oldnumeric as Numeric
    import time
    import string
    import numpy
    import copy
 
    # --------------------------------------------------------------------
    retv5km = Calipso5kmAvhrr5kmTrackObject()
    
    # Convert from TAI time to UTC in seconds since 1970:
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    
    lonCalipso5km = calipso5kmObj.longitude[:,1].ravel()
    latCalipso5km = calipso5kmObj.latitude[:,1].ravel()
    timeCalipso5km = calipso5kmObj.time[:,1].ravel() + dsec

    ndim5km = lonCalipso5km.shape[0]

    if avhrr5kmGeoObj.sec1970_end<avhrr5kmGeoObj.sec1970_start:
        avhrr5km_sec1970_end = int(DSEC_PER_AVHRR_SCALINE5KM*avhrr5kmObj.num_of_lines+avhrr5kmGeoObj.sec1970_start)
        avhrr5kmGeoObj.sec1970_end = avhrr5km_sec1970_end
    
    # --------------------------------------------------------------------
    cal5km,cap5km = get_calipso5km_avhrr5km_linpix(avhrr5kmGeoObj,ctypefile5km,lonCalipso5km,latCalipso5km,timeCalipso5km)

    # --------------------------------------------------------------------    
    
    idx_calipso5km = numpy.arange(ndim5km)
    idx_calipso5km_on_avhrr5km=numpy.repeat(idx_calipso5km,numpy.not_equal(cal5km,NODATA))
    
    lon_calipso5km = numpy.repeat(lonCalipso5km,numpy.not_equal(cal5km,NODATA))
    lat_calipso5km = numpy.repeat(latCalipso5km,numpy.not_equal(cal5km,NODATA))
    # Calipso line,pixel inside AVHRR swath:
    cal_on_avhrr5km = numpy.repeat(cal5km,numpy.not_equal(cal5km,NODATA))
    cap_on_avhrr5km = numpy.repeat(cap5km,numpy.not_equal(cal5km,NODATA))

    timetup_start5km = time.gmtime(avhrr5kmGeoObj.sec1970_start)
    timetup_end5km = time.gmtime(avhrr5kmGeoObj.sec1970_end)

    # Convert from TAI time to UTC in seconds since 1970:
    #dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone

    print "Start and end times: ",time.gmtime(calipso5kmObj.time[0][1] + dsec),time.gmtime(calipso5kmObj.time[ndim5km-1][1] + dsec)
    #print "Start and end times: ",time.gmtime(calipsoObj.time[0][0] + dsec),time.gmtime(calipsoObj.time[ndim-1][0] + dsec)
    #sec_timeThr = 60*20 # Allow for 20 minute deviation between AVHRR and CALIPSO matchup

    # --------------------------------------------------------------------
    secTup5km = avhrr5kmGeoObj.sec1970_start,avhrr5kmGeoObj.sec1970_end
    
    idx_match5km = select_calipso5km_inside_avhrr5km(calipso5kmObj,cal5km,dsec,secTup5km,sec_timeThr)

    retv5km.calipso5km.sec_1970 = numpy.repeat(calipso5kmObj.time[::,1] + dsec,idx_match5km)

    retv5km.calipso5km.cloud_fraction = numpy.repeat(calipso5kmObj.cloud_fraction,idx_match5km)
    retv5km.calipso5km.latitude = numpy.repeat(latCalipso5km,idx_match5km)
    retv5km.calipso5km.longitude = numpy.repeat(lonCalipso5km,idx_match5km)

     # --------------------------------------------------------------------
    
    print "cap_on_avhrr.shape: ",cap_on_avhrr5km.shape
    retv5km.calipso5km.avhrr_linnum = cal_on_avhrr5km.astype('i')
    retv5km.calipso5km.avhrr_pixnum = cap_on_avhrr5km.astype('i')

    #print "Concatenate arrays..."
    #x = Numeric.concatenate((idx_match,idx_match))
    #for i in range(2,10):
    #    x = Numeric.concatenate((x,idx_match))
    #idx_match_2d = Numeric.reshape(x,(ndim,10))

    print "Make cloud top and base arrays..."
    missing_data = -9.9
    #cloud_top = Numeric.repeat(calipsoObj.cloud_top_profile.flat,idx_match_2d.flat)
    #cloud_top = Numeric.where(Numeric.less(cloud_top,0),missing_data,cloud_top)
    #N = cloud_top.flat.shape[0]/10
    #cloud_top = Numeric.reshape(cloud_top,(N,10))
    
    x = numpy.repeat(calipso5kmObj.feature_classification_flags[::,0],idx_match5km)
    for i in range(1,10):
        x = numpy.concatenate(\
            (x,numpy.repeat(calipso5kmObj.feature_classification_flags[::,i],idx_match5km)))
    N = x.shape[0]/10
    retv5km.calipso5km.feature_classification_flags = numpy.reshape(x,(10,N)).astype('i')

    x = numpy.repeat(calipso5kmObj.cloud_top_profile[::,0],idx_match5km)
    for i in range(1,10):
        x = numpy.concatenate(\
            (x,numpy.repeat(calipso5kmObj.cloud_top_profile[::,i],idx_match5km)))
    N = x.shape[0]/10
    retv5km.calipso5km.cloud_top_profile = numpy.reshape(x,(10,N)).astype('d')

    print "Calipso observation time of first calipso-avhrr match: ",\
          time.gmtime(retv5km.calipso5km.sec_1970[0])
    print "Calipso observation time of last calipso-avhrr match: ",\
          time.gmtime(retv5km.calipso5km.sec_1970[N-1])
    
    x = numpy.repeat(calipso5kmObj.cloud_base_profile[::,0],idx_match5km)
    for i in range(1,10):
        x = numpy.concatenate(\
            (x,numpy.repeat(calipso5kmObj.cloud_base_profile[::,i],idx_match5km)))
    N = x.shape[0]/10
    retv5km.calipso5km.cloud_base_profile = numpy.reshape(x,(10,N)).astype('d')

    x = numpy.repeat(calipso5kmObj.cloud_mid_temperature[::,0],idx_match5km)
    for i in range(1,10):
        x = numpy.concatenate(\
            (x,numpy.repeat(calipso5kmObj.cloud_mid_temperature[::,i],idx_match5km)))
    N = x.shape[0]/10
    retv5km.calipso5km.cloud_mid_temperature = numpy.reshape(x,(10,N)).astype('d')
    
    x = numpy.repeat(calipso5kmObj.optical_depth[::,0],idx_match5km)
    for i in range(1,10):
        x = numpy.concatenate(\
            (x,numpy.repeat(calipso5kmObj.optical_depth[::,i],idx_match5km)))
    N = x.shape[0]/10
    retv5km.calipso5km.optical_depth = numpy.reshape(x,(10,N)).astype('d')
    
    x = numpy.repeat(calipso5kmObj.optical_depth_uncertainty[::,0],idx_match5km)
    for i in range(1,10):
        x = numpy.concatenate(\
            (x,numpy.repeat(calipso5kmObj.optical_depth_uncertainty[::,i],idx_match5km)))
    N = x.shape[0]/10
    retv5km.calipso5km.optical_depth_uncertainty = numpy.reshape(x,(10,N)).astype('d')

    #cloud_mid_temp = Numeric.repeat(calipsoObj.cloud_mid_temperature.flat,idx_match_2d.flat)
    #cloud_mid_temp = Numeric.where(Numeric.less(cloud_mid_temp,0),missing_data,cloud_mid_temp)
    #cloud_mid_temp = Numeric.reshape(cloud_mid_temp,(N,10))
    #retv.calipso5km.cloud_mid_temperature = cloud_mid_temp
    
    # IGBP Land Cover:
    retv5km.calipso5km.igbp = numpy.repeat(calipso5kmObj.igbp.ravel(),idx_match5km.ravel())

    # NSIDC Ice and Snow Cover:
    retv5km.calipso5km.nsidc = numpy.repeat(calipso5kmObj.nsidc.ravel(),idx_match5km.ravel())

    # Elevation is given in km's. Convert to meters:
    retv5km.calipso5km.elevation = numpy.repeat(calipso5kmObj.elevation[::,2].ravel()*1000.0,idx_match5km.ravel()).astype('d')
                                            
    retv5km.calipso5km.number_of_layers_found = numpy.repeat(\
        calipso5kmObj.number_of_layers_found.ravel(),idx_match5km.ravel()).astype('i')

    retv5km.avhrr5km.sec_1970 = numpy.add(avhrr5kmGeoObj.sec1970_start,\
                                      cal_on_avhrr5km * DSEC_PER_AVHRR_SCALINE5KM)
    retv5km.diff_sec_1970 = retv5km.calipso5km.sec_1970 - retv5km.avhrr5km.sec_1970
    min_diff = numpy.minimum.reduce(retv5km.diff_sec_1970)
    max_diff = numpy.maximum.reduce(retv5km.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-calipso): ",\
          numpy.maximum.reduce(retv5km.diff_sec_1970),numpy.minimum.reduce(retv5km.diff_sec_1970)

    print "AVHRR observation time of first calipso-avhrr match: ",\
          time.gmtime(retv5km.avhrr5km.sec_1970[0])
    print "AVHRR observation time of last calipso-avhrr match: ",\
          time.gmtime(retv5km.avhrr5km.sec_1970[N-1])

    # Make the latitude and pps cloudtype on the calipso track:
    # line and pixel arrays have equal dimensions
    print "Generate the latitude,cloudtype tracks!"

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
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso[i],latCalipso[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso5km[i],latCalipso5km[i],idx_match5km[i])))
    basename = os.path.basename(ctypefile5km).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_")
    
    datapath = "%s/%s/5km/%s/%s/%s" %(DATA_DIR, base_sat, base_year, base_month, AREA5KM)
    if not os.path.exists(datapath):
        os.makedirs(datapath)
    fd = open("%s/5km_%s_cloudtype_calipso_track2.txt"%(datapath, basename),"w")
    fd.writelines(ll)
    fd.close()

    ll = []
    for i in range(N):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_calipso5km[i],lat_calipso5km[i],0)))
    fd = open("%s/5km_%s_cloudtype_calipso_track_excl.txt"%(datapath, basename),"w")
    fd.writelines(ll)
    fd.close()
    
    # CALIOP Maximum cloud top in km:
    max_cloud_top_calipso = numpy.maximum.reduce(retv5km.calipso5km.cloud_top_profile.ravel())
    print "max_cloud_top_calipso: ",max_cloud_top_calipso

    return retv5km,min_diff,max_diff

# -----------------------------------------------------
def select_calipso5km_inside_avhrr5km(calipso5kmObj,cal5km,dsec,sec1970_start_end,sec_timeThr):
    #import numpy.oldnumeric as Numeric
    import numpy

    sec1970_start,sec1970_end = sec1970_start_end
    
    # Select the points inside the avhrr swath:
    # Allowing for sec_timeThr seconds deviation:
    idx_time_okay5km = numpy.logical_and(numpy.greater(\
        calipso5kmObj.time[:,1],sec1970_start - dsec - sec_timeThr),
                                   numpy.less(\
        calipso5kmObj.time[:,1],sec1970_end - dsec   + sec_timeThr))
        
    idx_place_okay5km = numpy.where(numpy.not_equal(cal5km,NODATA),idx_time_okay5km,False)
    idx_match5km = idx_place_okay5km
    #idx_match5km1 = numpy.logical_and(numpy.not_equal(cal5km,NODATA),idx_time_okay5km)
    #test=[]
    #test1=[]
    #for i in range(4208): #Dont "hard coda" thise
    #    if idx_match5km[i] == idx_match5km1[i]:
    #        test.append(0)
    #    else:
    #        test.append(1)
    #    if idx_time_okay5km[i] == True:
    #        test1.append(0)
    #    else:
    #        test1.append(1)
    #idx_match5km = Numeric.logical_and(Numeric.greater(lin,0),idx_okay5km[::,0])
    #idx_match5km = Numeric.logical_and(idx_match,Numeric.logical_and(Numeric.greater(pix,0),Numeric.less_equal(pix,2048)))
    #print "Number of matches: ",Numeric.repeat(idx_match,idx_match).shape[0]

    # Get the PPS Cloud Types matching CALIPSO:
    #line = Numeric.repeat(lin,idx_match)
    #line = Numeric.floor(line+0.5).astype('i')
    #pixel = Numeric.repeat(pix,idx_match)
    #pixel = Numeric.floor(pixel+0.5).astype('i')
    #print "Number of matches: ",line.shape[0]

    return idx_match5km

# -----------------------------------------------------
def get_calipso5km(filename5km):
    import _pypps_filters
    import numpy.oldnumeric as Numeric
    import numpy

    # Read CALIPSO Lidar (CALIOP) data:        

    calipso5km = read_calipso5km(filename5km)

    caliop5km_max_height = numpy.ones(calipso5km.cloud_top_profile[::,0].shape)*-9
    for i in range(10):
        caliop5km_max_height = numpy.maximum(caliop5km_max_height,calipso5km.cloud_top_profile[::,i] * 1000.)
    
    calipso5km_clmask = Numeric.greater(caliop5km_max_height,0).astype('d')

    calipso5km.cloud_fraction = calipso5km_clmask.copy()

    return calipso5km

# -----------------------------------------------------
def read_calipso5km(filename5km):
    import _pyhl
    
    a=_pyhl.read_nodelist(filename5km)
    b=a.getNodeNames()
    a.selectAll()
    a.fetch()
    
    retv5km = Calipso5kmObject()

    c=a.getNode("/Longitude")
    retv5km.longitude=c.data().astype('d')
    c=a.getNode("/Latitude")
    retv5km.latitude=c.data().astype('d')
    c=a.getNode("/Profile_Time") # Internatiopnal Atomic Time (TAI) seconds from Jan 1, 1993
    retv5km.time=c.data()
    c=a.getNode("/Profile_UTC_Time") # TAI time converted to UTC and stored in format yymmdd.fffffff    
    retv5km.utc_time=c.data()

    c=a.getNode("/Feature_Classification_Flags")
    retv5km.feature_classification_flags=c.data().astype('i')
    c=a.getNode("/Layer_Top_Altitude")
    retv5km.cloud_top_profile=c.data()
    c=a.getNode("/Layer_Base_Altitude")
    retv5km.cloud_base_profile=c.data()
    c=a.getNode("/Number_Layers_Found")
    retv5km.number_of_layers_found=c.data()
    #c=a.getNode("/closest_calipso_cloud_fraction")
    #retv5km.cloud_fraction=c.data()
    c=a.getNode("/Midlayer_Temperature")
    retv5km.cloud_mid_temperature=c.data()

    c=a.getNode("/Day_Night_Flag")
    retv5km.day_night_flag=c.data()
    c=a.getNode("/DEM_Surface_Elevation")
    retv5km.elevation=c.data()
    c=a.getNode("/IGBP_Surface_Type")
    retv5km.igbp=c.data()
    c=a.getNode("/NSIDC_Surface_Type")
    retv5km.nsidc=c.data()
    c=a.getNode("/Feature_Optical_Depth_532")
    retv5km.optical_depth=c.data()
    c=a.getNode("/Feature_Optical_Depth_Uncertainty_532")
    retv5km.optical_depth_uncertainty=c.data()
    ##########################################################################################################################################
    #pdb.set_trace()
    return retv5km
# -----------------------------------------------------
def reshapeCalipso5km(calipsofiles,avhrr):
    import time
    import numpy
    import sys
    
    cal= Calipso5kmObject()
    if avhrr.sec1970_end<avhrr.sec1970_start:
        avhrr_end = int(DSEC_PER_AVHRR_SCALINE5KM*avhrr.num_of_lines+avhrr.sec1970_start)    
    else:
        avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
     
#    filenameClTime = './testdata/5km/matched_files/ClTime/Temp.h5'
#    clSat = readCloudSatStartEnd(filenameClTime)
#    clSat.CloudSatStart
#    clSat.CloudSatEnd
    startCalipso = get_calipso5km(calipsofiles[0])
    for i in range(len(calipsofiles)-1):
        newCalipso = get_calipso5km(calipsofiles[i+1])
        cal_start_all = startCalipso.time[:,1]+dsec
        cal_new_all = newCalipso.time[:,1]+dsec
        
        if not cal_start_all[0]<cal_new_all[0]:
            print "calipso files are in the wrong order"
            print("Program calipso5km.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9)
        
        cal_break = numpy.argmin(numpy.abs(cal_start_all - cal_new_all[0]))+1
        
        startCalipso.time = numpy.concatenate((startCalipso.time[0:cal_break,:],newCalipso.time))
        startCalipso.longitude = numpy.concatenate((startCalipso.longitude[0:cal_break,:],newCalipso.longitude))
        startCalipso.latitude = numpy.concatenate((startCalipso.latitude[0:cal_break,:],newCalipso.latitude))
        startCalipso.cloud_fraction = numpy.concatenate((startCalipso.cloud_fraction[0:cal_break],newCalipso.cloud_fraction))
        startCalipso.cloud_top_profile = numpy.concatenate((startCalipso.cloud_top_profile[0:cal_break,:],newCalipso.cloud_top_profile))
        startCalipso.cloud_base_profile = numpy.concatenate((startCalipso.cloud_base_profile[0:cal_break,:],newCalipso.cloud_base_profile))
        startCalipso.cloud_mid_temperature = numpy.concatenate((startCalipso.cloud_mid_temperature[0:cal_break,:],newCalipso.cloud_mid_temperature))
        startCalipso.number_of_layers_found = numpy.concatenate((startCalipso.number_of_layers_found[0:cal_break,:],newCalipso.number_of_layers_found))
        startCalipso.igbp = numpy.concatenate((startCalipso.igbp[0:cal_break,:],newCalipso.igbp))
        startCalipso.nsidc = numpy.concatenate((startCalipso.nsidc[0:cal_break,:],newCalipso.nsidc))
        startCalipso.elevation = numpy.concatenate((startCalipso.elevation[0:cal_break,:],newCalipso.elevation))
        startCalipso.utc_time = numpy.concatenate((startCalipso.utc_time[0:cal_break,:],newCalipso.utc_time))
        startCalipso.feature_classification_flags = numpy.concatenate((startCalipso.feature_classification_flags[0:cal_break,:],newCalipso.feature_classification_flags))
        startCalipso.day_night_flag = numpy.concatenate((startCalipso.day_night_flag[0:cal_break,:],newCalipso.day_night_flag))
        startCalipso.optical_depth = numpy.concatenate((startCalipso.optical_depth[0:cal_break,:],newCalipso.optical_depth))
        startCalipso.optical_depth_uncertainty = numpy.concatenate((startCalipso.optical_depth_uncertainty[0:cal_break,:],newCalipso.optical_depth_uncertainty))

    start_break = numpy.argmin((numpy.abs((startCalipso.time[:,1] + dsec) - (avhrr_start - sec_timeThr))))-1 # Minus one to get one extra, just to be certain
    end_break = numpy.argmin((numpy.abs((startCalipso.time[:,1] + dsec) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain    

    cal.time = startCalipso.time[start_break:end_break,:]    
    cal.longitude = startCalipso.longitude[start_break:end_break,:]
    cal.latitude = startCalipso.latitude[start_break:end_break,:]
    cal.cloud_fraction = startCalipso.cloud_fraction[start_break:end_break]
    cal.cloud_top_profile = startCalipso.cloud_top_profile[start_break:end_break,:]
    cal.cloud_base_profile=startCalipso.cloud_base_profile[start_break:end_break,:]
    cal.cloud_mid_temperature=startCalipso.cloud_mid_temperature[start_break:end_break,:]
    cal.number_of_layers_found=startCalipso.number_of_layers_found[start_break:end_break,:]
    cal.igbp=startCalipso.igbp[start_break:end_break,:]
    cal.nsidc=startCalipso.nsidc[start_break:end_break,:]
    cal.elevation=startCalipso.elevation[start_break:end_break,:]
    cal.utc_time=startCalipso.utc_time[start_break:end_break,:]
    cal.feature_classification_flags=startCalipso.feature_classification_flags[start_break:end_break,:]
    cal.day_night_flag=startCalipso.day_night_flag[start_break:end_break,:]
    cal.optical_depth=startCalipso.optical_depth[start_break:end_break,:]
    cal.optical_depth_uncertainty=startCalipso.optical_depth_uncertainty[start_break:end_break,:]
    
    if cal.time.shape[0] <= 0:
        print("No time match, please try with some other Calipso files")
        print("Program calipso5km.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)
        
    ##########################################################################################################################################
    #pdb.set_trace()  
    return cal
# -----------------------------------------------------
if __name__ == "__main__":
    # Testing:
    import string
    import epshdf
    import pps_io
    import numpy.oldnumeric as Numeric
    
    MAIN_DIR = "/local_disk/calipso_data"
    SUB_DIR = "noaa18_calipso_2007Aug"

    #PPS_DIR = "/local_disk/data/export"
    PPS_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
    AVHRR_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)
    CALIPSO_DIR = "%s/%s"%(MAIN_DIR,SUB_DIR)

    calipsofile = "%s/CAL_LID_L2_01kmCLay-Prov-V1-20.2007-08-24T10-54-14ZD.h5"%(CALIPSO_DIR)
    ctypefile = "%s/noaa18_20070824_1121_11649_satproj_00000_05012_cloudtype.h5"%(PPS_DIR)
    ctthfile = "%s/noaa18_20070824_1121_11649_satproj_00000_05012_ctth.h5"%(PPS_DIR)
    avhrrfile = "%s/noaa18_20070824_1121_11649_satproj_00000_05012_avhrr.h5"%(AVHRR_DIR)

    sl = string.split(os.path.basename(ctypefile),"_")
    platform = sl[0]
    norbit = string.atoi(sl[3])
    yyyymmdd = sl[1]

    # Read AVHRR lon,lat data
    write_log("INFO","Read AVHRR geolocation data")
    avhrrGeoObj = pps_io.readAvhrrGeoData(avhrrfile)

    # Read PPS Cloud Type data
    write_log("INFO","Read PPS Cloud Type")
    ctype = epshdf.read_cloudtype(ctypefile,1,1,0)
    ctth = epshdf.read_cloudtop(ctthfile,1,1,1,0,1)
    
    # --------------------------------------------------------------------
    write_log("INFO","Read CALIPSO data")
    # Read CALIPSO Lidar (CALIOP) data:
    calipso = get_calipso(calipsofile)

    lonCalipso = calipso.longitude.ravel()
    latCalipso = calipso.latitude.ravel()

    # Calculations with AAPP in ERROR!!! Fixme, Ad 2007-09-19
    #lin,pix = avhrr_linepix_from_lonlat_aapp(lonCalipso,latCalipso,avhrrGeoObj,platform,norbit,yyyymmdd)


    import numpy.oldnumeric as Numeric
    caliop_height = []
    caliop_base = []
    caliop_max_height = Numeric.ones(calipso.cloud_top_profile[::,0].shape)*-9
    for i in range(10):
        hh = Numeric.where(Numeric.greater(calipso.cloud_top_profile[::,i],-9),
                           calipso.cloud_top_profile[::,i] * 1000.,-9)
        caliop_max_height = Numeric.maximum(caliop_max_height,
                                            calipso.cloud_top_profile[::,i] * 1000.)
        caliop_height.append(hh)
        bb = Numeric.where(Numeric.greater(calipso.cloud_base_profile[::,i],-9),
                           calipso.cloud_base_profile[::,i] * 1000.,-9)
        caliop_base.append(bb)

    x = Numeric.repeat(calipso.number_of_layers_found.ravel(),
                       Numeric.greater(calipso.number_of_layers_found.ravel(),0))
    print "Number of points with more than 0 layers: ",x.shape[0]
    
    cal_data_ok = Numeric.greater(caliop_max_height,-9.)

    # Testing...
    import calipso_avhrr_matchup
    caObj = calipso_avhrr_matchup.getCaliopAvhrrMatch(avhrrfile,calipsofile,ctypefile,ctthfile)
    import time
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    print "Original: ",calipso.time[16203,0]+dsec
    print "Matchup:  ",caObj.calipso.sec_1970[3421]
    print calipso.cloud_top_profile[16203]
    print caObj.calipso.cloud_top_profile[::,3421]
    
    
