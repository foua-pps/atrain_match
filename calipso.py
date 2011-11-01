
import pdb
import inspect
import os
from pps_basic_configure import *
from pps_error_messages import * #@UnusedWildImport

from config import AREA, SUB_DIR, DATA_DIR, sec_timeThr, RESOLUTION, COMPRESS_LVL, NLINES, SWATHWD, NODATA
from common import MatchupError, elements_within_range

COVERAGE_DIR = "%s/%ikm/%s"%(SUB_DIR,RESOLUTION,AREA)
class DataObject(object):
    """
    Class to handle data objects with several arrays.
    
    """
    def __getattr__(self, name):
        try:
            return self.all_arrays[name]
        except KeyError:
            raise AttributeError("%s instance has no attribute '%s'" % (self.__class__.__name__, name))
    
    def __setattr__(self, name, value):
        if name == 'all_arrays':
            object.__setattr__(self, name, value)
        else:
            self.all_arrays[name] = value
            
class ppsAvhrrObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'sec_1970': None,
            'ctth_height': None,
            'ctth_pressure': None,
            'ctth_temperature': None,
            'cloudtype': None,
            'surftemp': None,
            'bt11micron': None,
            'bt12micron': None,
            'satz': None,
            'surftemp': None
            }
        
class CalipsoObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'avhrr_linnum': None,
            'avhrr_pixnum': None,
            'cloud_fraction': None,
            'cloud_top_profile': None,
            'cloud_base_profile': None,
            'cloud_mid_temperature': None,
            'number_of_layers_found': None,
            'igbp': None,
            'nsidc': None,
            'elevation': None,
            'time': None,
            'utc_time': None, 
            'sec_1970': None,
            'feature_classification_flags': None,
            'day_night_flag': None,
            'optical_depth': None,
            'optical_depth_uncertainty': None
            }
class CalipsoAvhrrTrackObject:
    def __init__(self):
        self.avhrr=ppsAvhrrObject()
        self.calipso=CalipsoObject()
        self.diff_sec_1970=None

class area_interface:
    pass

class SatProjCov:
    def __init__(self):
        self.coverage=None
        self.colidx=None
        self.rowidx=None 
# ----------------------------------------
def readCaliopAvhrrMatchObj(filename):
    import _pyhl #@UnresolvedImport

    retv = CalipsoAvhrrTrackObject()

    a=_pyhl.read_nodelist(filename)
    a.selectAll()
    a.fetch()

    # Match-Up - time difference:
    retv.diff_sec_1970 = a.getNode("/diff_sec_1970").data()

    # Calipso:
    retv.calipso.longitude = a.getNode("/calipso/longitude").data()
    retv.calipso.latitude = a.getNode("/calipso/latitude").data()
    retv.calipso.avhrr_linnum = a.getNode("/calipso/avhrr_linnum").data()
    retv.calipso.avhrr_pixnum = a.getNode("/calipso/avhrr_pixnum").data()
    
    retv.calipso.cloud_fraction = a.getNode("/calipso/cloud_fraction").data()
    retv.calipso.cloud_top_profile = a.getNode("/calipso/cloud_top_profile").data()
    retv.calipso.cloud_base_profile = a.getNode("/calipso/cloud_base_profile").data()
    retv.calipso.cloud_mid_temperature = a.getNode("/calipso/cloud_mid_temperature").data()
    retv.calipso.elevation = a.getNode("/calipso/elevation").data()
    retv.calipso.number_of_layers_found = a.getNode("/calipso/number_of_layers_found").data()
    try:
        retv.calipso.feature_classification_flags = a.getNode("/calipso/feature_classification_flags").data()
    except:
        print "No feature_classification_flags array in file!"
        pass
    
    retv.calipso.igbp = a.getNode("/calipso/igbp").data()
    retv.calipso.nsidc = a.getNode("/calipso/nsidc").data()
    #retv.calipso.utc_time = a.getNode("/calipso/utc_time").data()
    retv.calipso.sec_1970 = a.getNode("/calipso/sec_1970").data()
    if RESOLUTION == 5:
        retv.calipso.optical_depth = a.getNode("/calipso/optical_depth").data()
        retv.calipso.optical_depth_uncertainty = a.getNode("/calipso/optical_depth_uncertainty").data()
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
    retv.avhrr.satz = a.getNode("/avhrr/satz").data()
    try:
        retv.avhrr.surftemp = a.getNode("/avhrr/surftemp").data()
    except IOError:
        print("No surfttemp. Continue")
    return retv

# ----------------------------------------
def writeCaliopAvhrrMatchObj(filename,ca_obj,compress_lvl):

    import _pyhl #@UnresolvedImport
#    status = -1
    
    a=_pyhl.nodelist()

    shape = [ca_obj.calipso.longitude.shape[0]]

    # Match-Up - time difference:
    # ====
    b=_pyhl.node(_pyhl.DATASET_ID,"/diff_sec_1970")
    b.setArrayValue(1,shape,ca_obj.diff_sec_1970,"double",-1)
    a.addNode(b)
    
    # Calipso
    # ====
    b=_pyhl.node(_pyhl.GROUP_ID,"/calipso")
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/longitude")
    b.setArrayValue(1,shape,ca_obj.calipso.longitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/latitude")
    b.setArrayValue(1,shape,ca_obj.calipso.latitude,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/avhrr_linnum")
    b.setArrayValue(1,shape,ca_obj.calipso.avhrr_linnum,"int",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/avhrr_pixnum")
    b.setArrayValue(1,shape,ca_obj.calipso.avhrr_pixnum,"int",-1)
    a.addNode(b)

    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_fraction")
    b.setArrayValue(1,shape,ca_obj.calipso.cloud_fraction,"double",-1)
    a.addNode(b)

    shape2d = ca_obj.calipso.cloud_top_profile.shape
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_top_profile")
    b.setArrayValue(1,shape2d,ca_obj.calipso.cloud_top_profile,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_base_profile")
    b.setArrayValue(1,shape2d,ca_obj.calipso.cloud_base_profile,"double",-1)
    a.addNode(b)

    #b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_mid_profile")
    #b.setArrayValue(1,shape2d,ca_obj.calipso.cloud_mid_profile,"double",-1)
    #a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/cloud_mid_temperature")
    b.setArrayValue(1,shape2d,ca_obj.calipso.cloud_mid_temperature,"double",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/number_of_layers_found")
    b.setArrayValue(1,shape,ca_obj.calipso.number_of_layers_found,"int",-1)
    a.addNode(b)    
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/igbp")
    b.setArrayValue(1,shape,ca_obj.calipso.igbp,"uchar",-1)
    a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/nsidc")
    b.setArrayValue(1,shape,ca_obj.calipso.nsidc,"uchar",-1)
    a.addNode(b)
    
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/elevation")
    b.setArrayValue(1,shape,ca_obj.calipso.elevation,"double",-1)
    a.addNode(b)
    #b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/time")
    #b.setArrayValue(1,shape,ca_obj.calipso.time,"double",-1)
    #a.addNode(b)
    #b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/utc_time")
    #b.setArrayValue(1,shape,ca_obj.calipso.utc_time,"double",-1)
    #a.addNode(b)
    b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/sec_1970")
    b.setArrayValue(1,shape,ca_obj.calipso.sec_1970,"double",-1)
    a.addNode(b)
    try:
        b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/feature_classification_flags")
        print ca_obj.calipso.feature_classification_flags.shape
        #print ca_obj.calipso.feature_classification_flags.typecode()
        #print ca_obj.calipso.feature_classification_flags.itemsize()
        b.setArrayValue(1,shape2d,ca_obj.calipso.feature_classification_flags,"int",-1)
        a.addNode(b)
    except:
        pdb.set_trace()
        pass
    
    if RESOLUTION == 5:
        b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/optical_depth")
        b.setArrayValue(1,shape2d,ca_obj.calipso.optical_depth,"double",-1)
        a.addNode(b)
        b=_pyhl.node(_pyhl.DATASET_ID,"/calipso/optical_depth_uncertainty")
        b.setArrayValue(1,shape2d,ca_obj.calipso.optical_depth_uncertainty,"double",-1)
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
    try:        
        b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/surftemp")
        b.setArrayValue(1,shape,ca_obj.avhrr.surftemp,"double",-1)
        a.addNode(b)
    except AttributeError:
        print("No surftemp. Continue")
        
    b=_pyhl.node(_pyhl.DATASET_ID,"/avhrr/satz")
    b.setArrayValue(1,shape,ca_obj.avhrr.satz,"double",-1)
    a.addNode(b)
    status = a.write(filename,compress_lvl)

    return status

# ------------------------------------------------------------------
def writeCoverage(covIn,filename,inAid,outAid):
    import _pyhl #@UnresolvedImport
        
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

    a.write(filename,COMPRESS_LVL)
    
    return

# ------------------------------------------------------------------
def readCoverage(filename):
    import _pyhl #@UnresolvedImport #@UnresolvedImport

    a=_pyhl.read_nodelist(filename)
#    b=a.getNodeNames()
    
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
    
    retv = SatProjCov()
    retv.coverage = coverage.astype('Int8')
    retv.rowidx = rowidx.astype('Int16')
    retv.colidx = colidx.astype('Int16')

    return retv,info

# -----------------------------------------------------
def define_pcs(id,name,definition):
    import pcs
    p = pcs.usgs(name,definition)
    pcs.register(id,p)

# -----------------------------------------------------
def define_longlat_ll(id, name, pcs_id, ll_ll, size, scale):
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
def sec1970_to_julianday(sec1970):
    import pps_time_util
    import time
    
    year,month,day,hour,minutes,sec,dummy,dummy,dummy = time.gmtime(sec1970)
    jday_1950 = int(pps_time_util.getJulianDay(year,month,day) - pps_time_util.getJulianDay(1950,1,1))
    jday = jday_1950 + (hour+minutes/60.0+sec/3600)/24.0

    return jday

# -----------------------------------------------------
def avhrr_linepix_from_lonlat_aapp(lon,lat,avhrrObj,platform,norbit,yyyymmdd):
    import numpy
    import CreateAngles
    import _py_linepix_lonlat #@UnresolvedImport
    
    ndim = lon.shape[0]
    lin = numpy.zeros((ndim,),'d')
    pix = numpy.zeros((ndim,),'d')

    if platform.find("metop") >= 0:
        file_satpos = "%s/satpos_M%.2d_%s.txt"%(SATPOS_DIR,string.atoi(platform.split("metop")[1]),yyyymmdd) #@UndefinedVariable
    else:
        file_satpos = "%s/satpos_%s_%s.txt"%(SATPOS_DIR,platform,yyyymmdd) #@UndefinedVariable

    file_ephe = "%s/ephe_%s.txt"%(EPHE_DIR,yyyymmdd) #@UndefinedVariable
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
    import numpy
    maxlon = numpy.maximum.reduce(lon.ravel())
    minlon = numpy.minimum.reduce(lon.ravel())
    maxlat = numpy.maximum.reduce(lat.ravel())
    minlat = numpy.minimum.reduce(lat.ravel())

    return minlon,minlat,maxlon,maxlat

# --------------------------------------------
def get_calipso_avhrr_linpix(avhrrIn,avhrrname,lon,lat,caTime):
    import numpy

    tmppcs="tmpproj"
    define_pcs(tmppcs, "Plate Caree, central meridian at 15E",
               ['proj=eqc','ellps=bessel', 'lon_0=15'])


    startline=0

    """
    # Test code:
    tmpaid="tmparea"
    mask = get_calipso_avhrr_linpix(avhrr,lon,lat,(startline,startline+NLINES),SWATHWD,tmppcs,tmpaid)
    import Image
    dimy,dimx = mask.shape
    that = Image.fromstring("L",(dimx,dimy),mask.tostring())
    that.save("./yt.png")
    that.thumbnail((dimx/8,dimy/8))
    that.save("./yt_thumbnail.png")
    """
    orbittime =  os.path.basename(avhrrname).split("_")[1:3]
#    Inside=0
#    HasEncounteredMatch=0
    i=0
    if not os.path.exists(COVERAGE_DIR):
        os.makedirs(COVERAGE_DIR)
    while startline < avhrrIn.longitude.shape[0]:
        write_log("INFO","Calling get_calipso_avhrr_linpix: start-line = ",startline)
        tmpaid = "tmparea_%d"%i
        endline = startline+NLINES
        coverage_filename = "%s/coverage_avhrr_caliop_matchup_%s_%s_%s_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrrIn.satellite,orbittime[0],orbittime[1],
                             startline,endline,tmpaid)
        write_log("INFO","Coverage filename = ",coverage_filename)
        cal,cap,ok = get_calipso_avhrr_linpix_segment(avhrrIn,lon,lat,caTime,
                                                      (startline,endline),
                                                      SWATHWD,tmppcs,tmpaid,
                                                      coverage_filename)
        if ok:
#            HasEncounteredMatch=1
            write_log("INFO","There was a match...")
        # Do not like this one /Erik    
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
def get_calipso_avhrr_linpix_segment(avhrrIn,lon,lat,catime,lines,swath_width,tmppcs,
                                     tmpaid,covfilename):
    import numpy
    import _satproj
    import area
#    import pcs
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
    """
    areaObj = area.area(AREA)
    # Should this one be used? /Erik
#    if not os.path.exists(covfilename):
#        write_log("INFO","Create Coverage map...")
#        cov = _satproj.create_coverage(areaObj,lonarr,latarr,1)
#        writeCoverage(cov,covfilename,"satproj",AREA)
#    else:
#        write_log("INFO","Read the AVHRR-CALIOP matchup coverage from file...")
#        cov,info = readCoverage(covfilename)
    # Do like this instead /Erik
    write_log("INFO","Create Coverage map...")
    lonarr = lonarr.astype('float64')
    latarr = latarr.astype('float64')
    
    cov = _satproj.create_coverage(areaObj,lonarr,latarr,0) #@UndefinedVariable
    writeCoverage(cov,covfilename,"satproj",AREA)

    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA) #@UndefinedVariable
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA) #@UndefinedVariable

    write_log("INFO","Go through calipso track:")
    calipso_avhrr_line = []
    calipso_avhrr_pixel = []
#    calipso_avhrr_line_time = []
#    calipso_avhrr_pixel_time = []
    for i in range(ndim):
        xy_tup=pps_gisdata.lonlat2xy(AREA,lon[i],lat[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
##        if(x < 4500 and x >= 0 and y >= 0 and y < 4500): Should be 5010!!!/KG
        dimx=mapped_line.shape[1]#1002#5010
        dimy=mapped_line.shape[0]#1002#5010
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            calipso_avhrr_line.append(mapped_line[y,x])
            calipso_avhrr_pixel.append(mapped_pixel[y,x])
#            calipso_avhrr_line_time.append(-9)
#            calipso_avhrr_pixel_time.append(-9)
        else:
            calipso_avhrr_line.append(-9)
            calipso_avhrr_pixel.append(-9)
#            calipso_avhrr_line_time.append(-9)
#            calipso_avhrr_pixel_time.append(-9)

    calipso_avhrr_line = numpy.array(calipso_avhrr_line)
    calipso_avhrr_pixel = numpy.array(calipso_avhrr_pixel)
#    calipso_avhrr_line_time = numpy.array(calipso_avhrr_line_time)
#    calipso_avhrr_pixel_time = numpy.array(calipso_avhrr_pixel_time)
    # Control the time diference
#    match_calipso_points = numpy.where(numpy.not_equal(calipso_avhrr_line,-9))
#    avhrr_time = (calipso_avhrr_line[match_calipso_points] * DSEC_PER_AVHRR_SCALINE) + avhrrIn.sec1970_start
#    cal_time = catime[match_calipso_points]
#    time_diff = avhrr_time-cal_time
    #max_time_diff_allowed = 50*60 #Based on that a lap is 102 min
#    max_time_diff_allowed = sec_timeThr
#    time_match = numpy.where(abs(time_diff)<max_time_diff_allowed)
#    if time_match[0].shape[0]==0:             
#        x=numpy.repeat(calipso_avhrr_line_time,numpy.not_equal(calipso_avhrr_line_time,-9))
#    else:
#        calipso_avhrr_line_time[match_calipso_points[0][time_match]] = calipso_avhrr_line[match_calipso_points[0][time_match]]
#        calipso_avhrr_pixel_time[match_calipso_points[0][time_match]] = calipso_avhrr_pixel[match_calipso_points[0][time_match]]
#        x=numpy.repeat(calipso_avhrr_line_time,numpy.not_equal(calipso_avhrr_line_time,-9))
    x = numpy.repeat(calipso_avhrr_line, numpy.not_equal(calipso_avhrr_line, -9))
    write_log("INFO","Number of matching points = ",x.shape[0])
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0

    return calipso_avhrr_line, calipso_avhrr_pixel, matchOk

#----------------------------------------------------------------------------------------------------
def createAvhrrTime(Obt, filename):
    import os #@Reimport
    import numpy
    from config import DSEC_PER_AVHRR_SCALINE

    filename = os.path.basename(filename)

    if Obt.sec1970_end < Obt.sec1970_start:
        """
        In some GAC edition the end time is negative. If so this if statement 
        tries to estimate the endtime depending on the start time plus number of 
        scanlines multiplied with the estimate scan time for the instrument. 
        This estimation is not that correct but what to do?
        """
        Obt.sec1970_end = int(DSEC_PER_AVHRR_SCALINE * Obt.num_of_lines + Obt.sec1970_start)
    
    if filename.split('_')[-3] != '00000':
        """
        This if statement takes care of a bug in start and end time, 
        that occurs when a file is cute at midnight
        """
        
        import calendar, time
        timediff = Obt.sec1970_end - Obt.sec1970_start
        old_start = time.gmtime(Obt.sec1970_start + (24 * 3600)) # Adds 24 h to get the next day in new start
        new_start = calendar.timegm(time.strptime('%i %i %i' %(old_start.tm_year, \
                                                               old_start.tm_mon, \
                                                               old_start.tm_mday), \
                                                               '%Y %m %d'))
        Obt.sec1970_start = new_start
        Obt.sec1970_end = new_start + timediff
    
    
    Obt.time = numpy.linspace(Obt.sec1970_start, Obt.sec1970_end, Obt.num_of_lines)
    
    return Obt

# -----------------------------------------------------
def match_calipso_avhrr(ctypefile,calipsoObj,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj):

    import time
    import string
    import numpy
    
    retv = CalipsoAvhrrTrackObject()
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone # Convert from TAI time to UTC in seconds since 1970
    if RESOLUTION == 1:
        lonCalipso = calipsoObj.longitude.ravel()
        latCalipso = calipsoObj.latitude.ravel()
        timeCalipso = calipsoObj.time[::,0].ravel() + dsec
        elevationCalipso = calipsoObj.elevation.ravel()
    if RESOLUTION == 5:
        # Use [:,1] Since 5km data has start, center, and end for each pixel
        lonCalipso = calipsoObj.longitude[:,1].ravel()
        latCalipso = calipsoObj.latitude[:,1].ravel()
        timeCalipso = calipsoObj.time[:,1].ravel() + dsec
        elevationCalipso = calipsoObj.elevation[::,2].ravel()
    
    ndim = lonCalipso.shape[0]
    
    # --------------------------------------------------------------------

    cal,cap = get_calipso_avhrr_linpix(avhrrGeoObj,ctypefile,lonCalipso,latCalipso,timeCalipso)
    calnan = numpy.where(cal == NODATA, numpy.nan, cal)
    if (~numpy.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
    avhrrGeoObj = createAvhrrTime(avhrrGeoObj, ctypefile)
    avhrr_lines_sec_1970 = numpy.where(cal != NODATA, avhrrGeoObj.time[cal], numpy.nan)
    
#    avhrr_lines_sec_1970 = calnan * DSEC_PER_AVHRR_SCALINE + avhrrGeoObj.sec1970_start
    # Find all matching Cloudsat pixels within +/- sec_timeThr from the AVHRR data
    idx_match = elements_within_range(timeCalipso, avhrr_lines_sec_1970, sec_timeThr) 

    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)
    
    lon_calipso = numpy.repeat(lonCalipso, idx_match)
    lat_calipso = numpy.repeat(latCalipso, idx_match)
    # Calipso line,pixel inside AVHRR swath:
    cal_on_avhrr = numpy.repeat(cal, idx_match)
    cap_on_avhrr = numpy.repeat(cap, idx_match)

    print "Start and end times: ",time.gmtime(timeCalipso[0]),time.gmtime(timeCalipso[ndim-1])
    
    retv.calipso.sec_1970 = numpy.repeat(timeCalipso,idx_match)

    retv.calipso.cloud_fraction = numpy.repeat(calipsoObj.cloud_fraction,idx_match)
    retv.calipso.latitude = numpy.repeat(latCalipso,idx_match)
    retv.calipso.longitude = numpy.repeat(lonCalipso,idx_match)
    
    print "cap_on_avhrr.shape: ",cap_on_avhrr.shape
    retv.calipso.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.calipso.avhrr_pixnum = cap_on_avhrr.astype('i')
    
    #print "Concatenate arrays..."
    #x = numpy.concatenate((idx_match,idx_match))
    #for i in range(2,10):
    #    x = numpy.concatenate((x,idx_match))
    #idx_match_2d = numpy.reshape(x,(ndim,10))

    print "Make cloud top and base arrays..."
#    missing_data = -9.9
    #cloud_top = numpy.repeat(calipsoObj.cloud_top_profile.flat,idx_match_2d.flat)
    #cloud_top = numpy.where(numpy.less(cloud_top,0),missing_data,cloud_top)
    #N = cloud_top.flat.shape[0]/10
    #cloud_top = numpy.reshape(cloud_top,(N,10))
    
    x_fcf = numpy.repeat(calipsoObj.feature_classification_flags[::,0],idx_match)
    x_ctp = numpy.repeat(calipsoObj.cloud_top_profile[::,0],idx_match)
    x_cbp = numpy.repeat(calipsoObj.cloud_base_profile[::,0],idx_match)
    x_cmt = numpy.repeat(calipsoObj.cloud_mid_temperature[::,0],idx_match)
    
    col_dim = calipsoObj.cloud_mid_temperature.shape[1]
    for i in range(1,col_dim):
        x_fcf = numpy.concatenate(\
            (x_fcf,numpy.repeat(calipsoObj.feature_classification_flags[::,i],idx_match)))
        x_ctp = numpy.concatenate(\
            (x_ctp,numpy.repeat(calipsoObj.cloud_top_profile[::,i],idx_match)))
        x_cbp = numpy.concatenate(\
            (x_cbp,numpy.repeat(calipsoObj.cloud_base_profile[::,i],idx_match)))
        x_cmt = numpy.concatenate(\
            (x_cmt,numpy.repeat(calipsoObj.cloud_mid_temperature[::,i],idx_match)))
    N_fcf = x_fcf.shape[0]/col_dim
    retv.calipso.feature_classification_flags = numpy.reshape(x_fcf,(col_dim,N_fcf)).astype('i')
    N_ctp = x_ctp.shape[0]/col_dim
    retv.calipso.cloud_top_profile = numpy.reshape(x_ctp,(col_dim,N_ctp)).astype('d')
    N_cbp = x_cbp.shape[0]/col_dim
    retv.calipso.cloud_base_profile = numpy.reshape(x_cbp,(col_dim,N_cbp)).astype('d')
    N_cmt = x_cmt.shape[0]/col_dim
    retv.calipso.cloud_mid_temperature = numpy.reshape(x_cmt,(col_dim,N_cmt)).astype('d')
    if RESOLUTION == 5:
        x_od = numpy.repeat(calipsoObj.optical_depth[::,0],idx_match)
        x_odu = numpy.repeat(calipsoObj.optical_depth_uncertainty[::,0],idx_match)
        for i in range(1, col_dim):
            x_od = numpy.concatenate(\
                (x_od,numpy.repeat(calipsoObj.optical_depth[::,i],idx_match)))
            x_odu = numpy.concatenate(\
                (x_odu,numpy.repeat(calipsoObj.optical_depth_uncertainty[::,i],idx_match)))
        N_od = x_od.shape[0]/col_dim
        retv.calipso.optical_depth = numpy.reshape(x_od,(col_dim,N_od)).astype('d')
        N_odu = x_odu.shape[0]/col_dim
        retv.calipso.optical_depth_uncertainty = numpy.reshape(x_odu,(col_dim,N_odu)).astype('d')

    #cloud_mid_temp = numpy.repeat(calipsoObj.cloud_mid_temperature.flat,idx_match_2d.flat)
    #cloud_mid_temp = numpy.where(numpy.less(cloud_mid_temp,0),missing_data,cloud_mid_temp)
    #cloud_mid_temp = numpy.reshape(cloud_mid_temp,(N,10))
    #retv.calipso.cloud_mid_temperature = cloud_mid_temp
    
    # IGBP Land Cover:
    retv.calipso.igbp = numpy.repeat(calipsoObj.igbp.ravel(),idx_match.ravel())

    # NSIDC Ice and Snow Cover:
    retv.calipso.nsidc = numpy.repeat(calipsoObj.nsidc.ravel(),idx_match.ravel())

    # Elevation is given in km's. Convert to meters:
    retv.calipso.elevation = numpy.repeat(elevationCalipso.ravel()*1000.0,
                                            idx_match.ravel()).astype('d')

    retv.calipso.number_of_layers_found = numpy.repeat(\
        calipsoObj.number_of_layers_found.ravel(),idx_match.ravel()).astype('i')
    
    # Time
    retv.avhrr.sec_1970 = avhrrGeoObj.time[cal_on_avhrr]
    retv.diff_sec_1970 = retv.calipso.sec_1970 - retv.avhrr.sec_1970

    min_diff = numpy.minimum.reduce(retv.diff_sec_1970)
    max_diff = numpy.maximum.reduce(retv.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-calipso): ",\
          numpy.maximum.reduce(retv.diff_sec_1970),numpy.minimum.reduce(retv.diff_sec_1970)

    print "AVHRR observation time of first calipso-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[0])
    print "AVHRR observation time of last calipso-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[N_cmt-1])

    # Make the latitude and pps cloudtype on the calipso track:
    # line and pixel arrays have equal dimensions
    print "Generate the latitude,cloudtype tracks!"

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
    #TODO: maby skould take consideration of missing data also, as for satz.
    for i in range(cal_on_avhrr.shape[0]):
        lat_avhrr_track.append(avhrrGeoObj.latitude[cal_on_avhrr[i],cap_on_avhrr[i]])
        lon_avhrr_track.append(avhrrGeoObj.longitude[cal_on_avhrr[i],cap_on_avhrr[i]])
        ctype_track.append(ctype.cloudtype[cal_on_avhrr[i],cap_on_avhrr[i]])
        if surft != None:
            print("No surface temperatur")
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
        if avhrrAngObj.satz.data[cal_on_avhrr[i],cap_on_avhrr[i]] == avhrrAngObj.satz.no_data or avhrrAngObj.satz.data[cal_on_avhrr[i],cap_on_avhrr[i]] == avhrrAngObj.satz.missing_data:
            ang = -9
        else:
            ang = avhrrAngObj.satz.data[cal_on_avhrr[i],cap_on_avhrr[i]] * avhrrAngObj.satz.gain + avhrrAngObj.satz.intercept
        satz_track.append(ang)
    retv.avhrr.latitude = numpy.array(lat_avhrr_track)
    retv.avhrr.longitude = numpy.array(lon_avhrr_track)
    retv.avhrr.cloudtype = numpy.array(ctype_track)
    retv.avhrr.bt11micron = numpy.array(bt11micron_track)
    retv.avhrr.bt12micron = numpy.array(bt12micron_track)
    retv.avhrr.satz = numpy.array(satz_track)
    if ctth:
        retv.avhrr.ctth_height = numpy.array(ctth_height_track)
        retv.avhrr.ctth_pressure = numpy.array(ctth_pressure_track)
        retv.avhrr.ctth_temperature = numpy.array(ctth_temperature_track)
    if surft != None:
        retv.avhrr.surftemp = numpy.array(surft_track)

    print "AVHRR-PPS Cloud Type,latitude: shapes = ",\
          retv.avhrr.cloudtype.shape,retv.avhrr.latitude.shape

    ll = []
    for i in range(ndim):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso[i],latCalipso[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso[i],latCalipso[i],idx_match[i])))
    basename = os.path.basename(ctypefile).split(".h5")[0]
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_")
    datapath = "%s/%s/1km/%s/%s/%s" %(DATA_DIR, base_sat, base_year, base_month, AREA)
    if not os.path.exists(datapath):
        os.makedirs(datapath)
        
    fd = open("%s/%skm_%s_cloudtype_calipso_track2.txt"%(datapath,RESOLUTION,basename),"w")
    fd.writelines(ll)
    fd.close()

    ll = []
    for i in range(N_cmt):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_calipso[i],lat_calipso[i],0)))
    fd = open("%s/%skm_%s_cloudtype_calipso_track_excl.txt"%(datapath,RESOLUTION,basename),"w")
    fd.writelines(ll)
    fd.close()
    
    # CALIOP Maximum cloud top in km:
    max_cloud_top_calipso = numpy.maximum.reduce(retv.calipso.cloud_top_profile.ravel())
    print "max_cloud_top_calipso: ",max_cloud_top_calipso

    return retv,min_diff,max_diff

# -----------------------------------------------------
def select_calipso_inside_avhrr(calipsoObj,cal,dsec,sec1970_start_end,sec_timeThr):
    import numpy

    sec1970_start,sec1970_end = sec1970_start_end
    
    # Select the points inside the avhrr swath:
    # Allowing for sec_timeThr seconds deviation:
    if RESOLUTION == 1:
        idx_time_okay = numpy.logical_and(numpy.greater(\
            calipsoObj.time[:,0],sec1970_start - dsec - sec_timeThr),
                                   numpy.less(\
            calipsoObj.time[:,0],sec1970_end - dsec   + sec_timeThr))
    elif RESOLUTION == 5:
        idx_time_okay = numpy.logical_and(numpy.greater(\
            calipsoObj.time[:,1],sec1970_start - dsec - sec_timeThr),
                                   numpy.less(\
            calipsoObj.time[:,1],sec1970_end - dsec   + sec_timeThr)) ##########################################################################################################################################   
    #pdb.set_trace()
    #idx_match = numpy.not_equal(cal,NODATA)        
    idx_place_okay = numpy.where(numpy.not_equal(cal,NODATA),idx_time_okay,False)
    idx_match = idx_place_okay
    
    #idx_match = numpy.logical_and(numpy.greater(lin,0),idx_okay[::,0])
    #idx_match = numpy.logical_and(idx_match,numpy.logical_and(numpy.greater(pix,0),numpy.less_equal(pix,2048)))
    #print "Number of matches: ",numpy.repeat(idx_match,idx_match).shape[0]

    # Get the PPS Cloud Types matching CALIPSO:
    #line = numpy.repeat(lin,idx_match)
    #line = numpy.floor(line+0.5).astype('i')
    #pixel = numpy.repeat(pix,idx_match)
    #pixel = numpy.floor(pixel+0.5).astype('i')
    #print "Number of matches: ",line.shape[0]

    return idx_match

# -----------------------------------------------------
def get_calipso(filename):
    import _pypps_filters
    import numpy
    # Read CALIPSO Lidar (CALIOP) data:
    calipso = read_calipso(filename)
    ndim = calipso.longitude.shape[0]
    if RESOLUTION == 1:
        lonCalipso = calipso.longitude.ravel()
        latCalipso = calipso.latitude.ravel()
        # --------------------------------------------------------------------
        # Derive the calipso cloud fraction using the 
        # cloud height:       
        winsz=3
        caliop_max_height = numpy.ones(calipso.cloud_top_profile[::,0].shape)*-9
        for i in range(calipso.cloud_top_profile.shape[1]):
            caliop_max_height = numpy.maximum(caliop_max_height,
                                                calipso.cloud_top_profile[::,i] * 1000.)
    
        #calipso_clmask = numpy.greater(calipso.cloud_base_profile[::,0],0).astype('d')
        calipso_clmask = numpy.greater(caliop_max_height,0).astype('d')
        cloud_mask = numpy.concatenate((calipso_clmask,calipso_clmask))
        for idx in range(2,winsz):
            cloud_mask = numpy.concatenate((cloud_mask,calipso_clmask))
        cloud_mask = numpy.reshape(cloud_mask,(winsz,ndim)).astype('d')
    
        calipso.cloud_fraction=numpy.zeros((winsz,ndim),'d')
        _pypps_filters.texture(cloud_mask,calipso.cloud_fraction,winsz,"mean") #@UndefinedVariable
        calipso.cloud_fraction = calipso.cloud_fraction[winsz/2,::]
    
    elif RESOLUTION == 5:
#        lonCalipso = calipso.longitude[:,1].ravel()
#        latCalipso = calipso.latitude[:,1].ravel()
        calipso.cloud_fraction = numpy.where(calipso.cloud_top_profile[:,0] > 0, 1, 0).astype('d')
    
    return calipso

# -----------------------------------------------------
def read_calipso(filename):
    import _pyhl #@UnresolvedImport

    a=_pyhl.read_nodelist(filename)
#    b=a.getNodeNames()
    a.selectAll()
    a.fetch()

    retv = CalipsoObject()

    c=a.getNode("/Longitude")
    retv.longitude=c.data().astype('d')
    c=a.getNode("/Latitude")
    retv.latitude=c.data().astype('d')
    c=a.getNode("/Profile_Time") # Internatiopnal Atomic Time (TAI) seconds from Jan 1, 1993
    retv.time=c.data()
    c=a.getNode("/Profile_UTC_Time") # TAI time converted to UTC and stored in format yymmdd.fffffff    
    retv.utc_time=c.data()

    c=a.getNode("/Feature_Classification_Flags")
    retv.feature_classification_flags=c.data().astype('uint16')
    c=a.getNode("/Layer_Top_Altitude")
    retv.cloud_top_profile=c.data()
    c=a.getNode("/Layer_Base_Altitude")
    retv.cloud_base_profile=c.data()
    c=a.getNode("/Number_Layers_Found")
    retv.number_of_layers_found=c.data()
    #c=a.getNode("/closest_calipso_cloud_fraction")
    #retv.cloud_fraction=c.data()
    c=a.getNode("/Midlayer_Temperature")
    retv.cloud_mid_temperature=c.data()

    c=a.getNode("/Day_Night_Flag")
    retv.day_night_flag=c.data()
    c=a.getNode("/DEM_Surface_Elevation")
    retv.elevation=c.data()
    c=a.getNode("/IGBP_Surface_Type")
    retv.igbp=c.data()
    c=a.getNode("/NSIDC_Surface_Type")
    retv.nsidc=c.data()
    if RESOLUTION == 5:
        c=a.getNode("/Feature_Optical_Depth_532")
        retv.optical_depth=c.data()
        c=a.getNode("/Feature_Optical_Depth_Uncertainty_532")
        retv.optical_depth_uncertainty=c.data()
    return retv
# -----------------------------------------------------
def reshapeCalipso(calipsofiles, avhrr, avhrrfilename):
    import time
    import numpy
    import sys
    
    cal= CalipsoObject()
    avhrr = createAvhrrTime(avhrr, avhrrfilename)
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    
    startCalipso = get_calipso(calipsofiles[0])
    for i in range(len(calipsofiles)-1):
        newCalipso = get_calipso(calipsofiles[i+1])
        if RESOLUTION == 1:
            cal_start_all = startCalipso.time[:,0]+dsec
            cal_new_all = newCalipso.time[:,0]+dsec
        elif RESOLUTION == 5:
            cal_start_all = startCalipso.time[:,1]+dsec
            cal_new_all = newCalipso.time[:,1]+dsec
        if not cal_start_all[0]<cal_new_all[0]:
            print "calipso files are in the wrong order"
            print("Program calipso.py at line %i" %(inspect.currentframe().f_lineno+1))
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
        if RESOLUTION == 5:
            startCalipso.optical_depth = numpy.concatenate((startCalipso.optical_depth[0:cal_break,:],newCalipso.optical_depth))
            startCalipso.optical_depth_uncertainty = numpy.concatenate((startCalipso.optical_depth_uncertainty[0:cal_break,:],newCalipso.optical_depth_uncertainty))

    if RESOLUTION == 1:
        start_break = numpy.argmin((numpy.abs((startCalipso.time[:,0] + dsec) - (avhrr_start - sec_timeThr))))-1 # Minus one to get one extra, just to be certain
        end_break = numpy.argmin((numpy.abs((startCalipso.time[:,0] + dsec) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain    
    if RESOLUTION == 5:
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
    if RESOLUTION == 5:        
        cal.optical_depth=startCalipso.optical_depth[start_break:end_break,:]
        cal.optical_depth_uncertainty=startCalipso.optical_depth_uncertainty[start_break:end_break,:]
    
    if cal.time.shape[0] <= 0:
        print("No time match, please try with some other Calipso files")
        print("Program calipso.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)

    return cal

# -----------------------------------------------------
if __name__ == "__main__":
    # Testing:
    import string
    import epshdf
    import pps_io
    import numpy
    
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

    caliop_height = []
    caliop_base = []
    caliop_max_height = numpy.ones(calipso.cloud_top_profile[::,0].shape)*-9
    for i in range(10):
        hh = numpy.where(numpy.greater(calipso.cloud_top_profile[::,i],-9),
                           calipso.cloud_top_profile[::,i] * 1000.,-9)
        caliop_max_height = numpy.maximum(caliop_max_height,
                                            calipso.cloud_top_profile[::,i] * 1000.)
        caliop_height.append(hh)
        bb = numpy.where(numpy.greater(calipso.cloud_base_profile[::,i],-9),
                           calipso.cloud_base_profile[::,i] * 1000.,-9)
        caliop_base.append(bb)

    x = numpy.repeat(calipso.number_of_layers_found.ravel(),
                       numpy.greater(calipso.number_of_layers_found.ravel(),0))
    print "Number of points with more than 0 layers: ",x.shape[0]
    
    cal_data_ok = numpy.greater(caliop_max_height,-9.)

    # Testing...
    import calipso_avhrr_matchup
    caObj = calipso_avhrr_matchup.getCaliopAvhrrMatch(avhrrfile,calipsofile,ctypefile,ctthfile)
    import time
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    print "Original: ",calipso.time[16203,0]+dsec
    print "Matchup:  ",caObj.calipso.sec_1970[3421]
    print calipso.cloud_top_profile[16203]
    print caObj.calipso.cloud_top_profile[::,3421]
    
    
