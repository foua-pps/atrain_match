import pdb #@UnusedImport

import inspect #@UnusedImport
import os #@UnusedImport
import numpy as np
from pps_basic_configure import *
from pps_error_messages import write_log

from config import (AREA, _validation_results_dir, 
                    sec_timeThr, COMPRESS_LVL, RESOLUTION,
                    NLINES, SWATHWD, NODATA) #@UnusedImport
from common import MatchupError, elements_within_range #@UnusedImport
from config import RESOLUTION as resolution
from config import (OPTICAL_DETECTION_LIMIT,
                    EXCLUDE_CALIPSO_PIXEL_IF_TOTAL_OPTICAL_THICKNESS_TO_LOW,
                    EXCLUDE_ALL_MULTILAYER,
                    EXCLUDE_MULTILAYER_IF_TOO_THIN_TOP_LAYER,
                    EXCLUDE_GEOMETRICALLY_THICK,
                    PPS_VALIDATION,
                    IMAGER_INSTRUMENT)
import time as tm

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
            'ctth_opaque': None,  # True if opaque retrieval was applied
            'cloudtype': None,
            'cloudtype_qflag': None,
            'surftemp': None,
            'bt11micron': None,
            'bt12micron': None,
            'satz': None,
            'lwp': None
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
            'optical_depth_uncertainty': None,
            'single_shot_cloud_cleared_fraction': None,
            'Horizontal_Averaging': None
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
    import h5py #@UnresolvedImport
    
    retv = CalipsoAvhrrTrackObject()
    
    h5file = h5py.File(filename, 'r')
    for group, data_obj in [(h5file['/calipso'], retv.calipso),
                            (h5file['/avhrr'], retv.avhrr)]:
        for dataset in group.keys():        
            if dataset in data_obj.all_arrays.keys():
                data_obj.all_arrays[dataset] = group[dataset].value

#    duplicate_names(retv.calipso)
    
    retv.diff_sec_1970 = h5file['diff_sec_1970'].value

    h5file.close()

    return retv

# ----------------------------------------
def writeCaliopAvhrrMatchObj(filename, ca_obj):
    """
    Write *ca_obj* to *filename*.
    
    """
    from common import write_match_objects
    groups = {'calipso': ca_obj.calipso.all_arrays,
              'avhrr': ca_obj.avhrr.all_arrays}
    write_match_objects(filename, ca_obj.diff_sec_1970, groups)
    
#    h5file.create_group('avhrr')
#    for arname, value in ca_obj.avhrr.all_arrays.items():
#        if value == None or value == []:
#            continue
#        else:
#            h5file.create_dataset(('avhrr' + '/' + arname), data = value)
#    h5file.close()

    status = 1
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
    import pcs #@UnresolvedImport
    p = pcs.usgs(name,definition)
    pcs.register(id,p)

# -----------------------------------------------------
def define_longlat_ll(id, name, pcs_id, ll_ll, size, scale):
    import pcs #@UnresolvedImport
    import area #@UnresolvedImport
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
    #import pps_time_util #@UnresolvedImport
    import time as tm
    import datetime
    year,month,day,hour,minutes,sec,tm_wday,tm_yday,tm_isdst = tm.gmtime(sec1970)
    daysdelta=datetime.datetime(year,month,day,00,0) - datetime.datetime(1950,1,1,00,0)
    jday50 = daysdelta.days
    #jday50 is the same as jday_1950
    #jday_1950 = int(pps_time_util.getJulianDay(year,month,day) - pps_time_util.getJulianDay(1950,1,1))
    jday = jday50 + (hour+minutes/60.0+sec/3600)/24.0
    if not jday==tm_yday:
        print "Is this (%f) really the julian day wanted?"%(jday)
        print "The day of the year is: (%d)"%(tm_yday)
        print "And if it days since 1 januari 1950, i would suggest:( %d)"%(jday50)
 
    return jday

# -----------------------------------------------------
def avhrr_linepix_from_lonlat_aapp(lon,lat,avhrrObj,platform,norbit,yyyymmdd):
    import CreateAngles #@UnresolvedImport
    import _py_linepix_lonlat #@UnresolvedImport
    
    ndim = lon.shape[0]
    lin = np.zeros((ndim,), 'd')
    pix = np.zeros((ndim,), 'd')

    if platform.find("metop") >= 0:
        file_satpos = "%s/satpos_M%.2d_%s.txt"%(SATPOS_DIR,string.atoi(platform.split("metop")[1]),yyyymmdd) #@UndefinedVariable
    else:
        file_satpos = "%s/satpos_%s_%s.txt"%(SATPOS_DIR,platform,yyyymmdd) #@UndefinedVariable

    file_ephe = "%s/ephe_%s.txt"%(EPHE_DIR,yyyymmdd) #@UndefinedVariable
    
    start_jday, end_jday, sec1970_start, sec1970_end = CreateAngles.get_acqtime(file_ephe,norbit)
    write_log("INFO","From ephemeris file: platform,norbit,start_jday,end_jday = ",platform,norbit,start_jday,end_jday) #@UndefinedVariable

    start_jday = sec1970_to_julianday(avhrrObj.sec1970_start)
    end_jday   = sec1970_to_julianday(avhrrObj.sec1970_end)
    write_log("INFO","From AVHRR file: platform,norbit,start_jday,end_jday = ",platform,norbit,start_jday,end_jday) #@UndefinedVariable

    #attitude = CreateAngles.get_attitude(platform,norbit,"tle")
    #attitude = avhrrObj.attitude_error["yaw"],avhrrObj.attitude_error["roll"],avhrrObj.attitude_error["pitch"]
    attitude = (0.,0.,0.)
    write_log("INFO","YAW,ROLL,PITCH:",attitude[0],attitude[1],attitude[2]) #@UndefinedVariable

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
    maxlon = np.maximum.reduce(lon.ravel())
    minlon = np.minimum.reduce(lon.ravel())
    maxlat = np.maximum.reduce(lat.ravel())
    minlat = np.minimum.reduce(lat.ravel())

    return minlon,minlat,maxlon,maxlat

# --------------------------------------------

def get_calipso_avhrr_linpix(avhrrIn, values, lon, lat, caTime, options):

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
    #orbittime =  os.path.basename(avhrrname).split("_")[1:3]
    coverage_dir  = options['coverage_dir'].format(resolution=str(RESOLUTION),
                                                   area=AREA,
                                                   val_dir=_validation_results_dir
                                                   )
    i=0
    if not os.path.exists(coverage_dir):
        os.makedirs(coverage_dir)
    while startline < avhrrIn.longitude.shape[0]:
        write_log("INFO","Calling get_calipso_avhrr_linpix start-line = ",startline)
        endline = startline + NLINES
        tmpaid = "tmparea_%d" %(i)
        coverage_filename =  options['coverage_filename'].format(
            satellite=avhrrIn.satellite,
            tmpaid=tmpaid,
            startline="%.5d"%(startline),
            endline="%.5d" %(endline),
            date=values["date"],#orbittime[0],
            time=values["time"],#orbittime[1],
            atrain_sat="calipso")

        write_log("INFO","Coverage filename = ",coverage_filename) #@UndefinedVariable
        coverage_file = os.path.join(coverage_dir, coverage_filename)
        cal,cap,ok = get_calipso_avhrr_linpix_segment(avhrrIn,lon,lat,caTime,
                                                      (startline,endline),
                                                      SWATHWD,tmppcs,tmpaid,
                                                      coverage_file)
        if ok:
#            HasEncounteredMatch=1
            write_log("INFO","There was a match...") #@UndefinedVariable
        # Do not like this one /Erik    
        #if not ok and HasEncounteredMatch:
        #    write_log("INFO","Data is now empty. Leave the loop...")
        #    break
        if(startline==0):
            # First time:
            calipso_avhrr_line,calipso_avhrr_pixel = np.array(cal),np.array(cap)
        else:
            # Merge:
            calipso_avhrr_line = np.where(np.equal(calipso_avhrr_line,-9),cal,calipso_avhrr_line)
            calipso_avhrr_pixel = np.where(np.equal(calipso_avhrr_pixel,-9),cap,calipso_avhrr_pixel)

        startline=startline+NLINES
        i=i+1

    return calipso_avhrr_line,calipso_avhrr_pixel

# --------------------------------------------
def get_calipso_avhrr_linpix_segment(avhrrIn,lon,lat,catime,lines,swath_width,tmppcs,
                                     tmpaid,covfilename):
    import _satproj #@UnresolvedImport
    import area #@UnresolvedImport
    import pps_gisdata #@UnresolvedImport

    ndim = lon.shape[0]

    if avhrrIn.longitude.shape[0] > lines[1]:
        lines_end = lines[1]
    else:
        lines_end = avhrrIn.longitude.shape[0]
    lines_start = lines[0]
    
    write_log("INFO","lines_start,lines_end: ",lines_start,lines_end) #@UndefinedVariable
    nlines = lines_end - lines_start
    lonarr = avhrrIn.longitude[lines_start:lines_end,::]
    latarr = avhrrIn.latitude[lines_start:lines_end,::]

    idx_start = lines_start*swath_width
    idx_end   = lines_end*swath_width
    idx = np.arange(idx_start,idx_end)
    
    linearr = np.divide(idx,swath_width)
    write_log("INFO","Start and end line numbers: ",linearr[0],linearr[idx.shape[0]-1]) #@UndefinedVariable
    
    linearr = np.reshape(linearr,(nlines,swath_width))
    pixelarr = np.fmod(idx,swath_width).astype('l')
    pixelarr = np.reshape(pixelarr,(nlines,swath_width))

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
    write_log("INFO","Create Coverage map...") #@UndefinedVariable
    lonarr = lonarr.astype('float64')
    latarr = latarr.astype('float64')
    
    cov = _satproj.create_coverage(areaObj,lonarr,latarr,0) #@UndefinedVariable
    writeCoverage(cov,covfilename,"satproj",AREA)

    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA) #@UndefinedVariable
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA) #@UndefinedVariable

    write_log("INFO","Go through calipso track:") #@UndefinedVariable
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
    calipso_avhrr_line = np.array(calipso_avhrr_line)
    calipso_avhrr_pixel = np.array(calipso_avhrr_pixel)
#    calipso_avhrr_line_time = np.array(calipso_avhrr_line_time)
#    calipso_avhrr_pixel_time = np.array(calipso_avhrr_pixel_time)
    # Control the time diference
#    match_calipso_points = np.where(np.not_equal(calipso_avhrr_line,-9))
#    avhrr_time = (calipso_avhrr_line[match_calipso_points] * DSEC_PER_AVHRR_SCALINE) + avhrrIn.sec1970_start
#    cal_time = catime[match_calipso_points]
#    time_diff = avhrr_time-cal_time
    #max_time_diff_allowed = 50*60 #Based on that a lap is 102 min
#    max_time_diff_allowed = sec_timeThr
#    time_match = np.where(abs(time_diff)<max_time_diff_allowed)
#    if time_match[0].shape[0]==0:             
#        x=np.repeat(calipso_avhrr_line_time,np.not_equal(calipso_avhrr_line_time,-9))
#    else:
#        calipso_avhrr_line_time[match_calipso_points[0][time_match]] = calipso_avhrr_line[match_calipso_points[0][time_match]]
#        calipso_avhrr_pixel_time[match_calipso_points[0][time_match]] = calipso_avhrr_pixel[match_calipso_points[0][time_match]]
#        x=np.repeat(calipso_avhrr_line_time,np.not_equal(calipso_avhrr_line_time,-9))
    x = np.repeat(calipso_avhrr_line, np.not_equal(calipso_avhrr_line, -9))
    write_log("INFO","Number of matching points = ",x.shape[0]) #@UndefinedVariable
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0

    return calipso_avhrr_line, calipso_avhrr_pixel, matchOk

#-----------------------------------------------------------------------------
def createAvhrrTime(Obt, values):
    import os #@Reimport
    from config import DSEC_PER_AVHRR_SCALINE
    #import time
    #from datetime import datetime
    import calendar

    #filename = os.path.basename(filename)
    # Ex.: npp_20120827_2236_04321_satproj_00000_04607_cloudtype.h5
    if IMAGER_INSTRUMENT == 'viirs':
    #if filename.split('_')[0] == 'npp':
        if Obt.sec1970_start < 0: #10800
            write_log("WARNING", 
                      "NPP start time negative! " + str(Obt.sec1970_start))
            datetime=values["date_time"]
            Obt.sec1970_start = calendar.timegm(datetime.timetuple())
            #Obt.sec1970_start = calendar.timegm((year, mon, day, hour, mins, sec)) + hundredSec
        num_of_scan = Obt.num_of_lines / 16.
        #if (Obt.sec1970_end - Obt.sec1970_start) / (num_of_scan) > 2:
        #    pdb.set_trace()
       #linetime = np.linspace(1, 10, 20)
       #test = np.apply_along_axis(np.multiply,  0, np.ones([20, 16]), linetime).reshape(30)        
        linetime = np.linspace(Obt.sec1970_start, Obt.sec1970_end, num_of_scan)
        Obt.time = np.apply_along_axis(np.multiply,  0, np.ones([num_of_scan, 16]), linetime).reshape(Obt.num_of_lines)
 

    else:
        if Obt.sec1970_end < Obt.sec1970_start:
            """
            In some GAC edition the end time is negative. If so this if statement 
            tries to estimate the endtime depending on the start time plus number of 
            scanlines multiplied with the estimate scan time for the instrument. 
            This estimation is not that correct but what to do?
            """
            Obt.sec1970_end = int(DSEC_PER_AVHRR_SCALINE * Obt.num_of_lines + Obt.sec1970_start)
        
        if values["ppsfilename"].split('_')[-3] != '00000':
            """
            This if statement takes care of a bug in start and end time, 
            that occurs when a file is cut at midnight
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
        Obt.time = np.linspace(Obt.sec1970_start, Obt.sec1970_end, Obt.num_of_lines)
    
    return Obt

#---------------------------------------------------------------------------
def avhrr_track_from_matched(obt, GeoObj, dataObj, AngObj, 
                             surft, ctth, ctype, 
                             row_matched, col_matched, 
                             avhrrLwp=None, avhrrCph=None):
    ctype_track = []
    ctype_qflag_track = []
    ctth_height_track = []
    ctth_pressure_track = []
    ctth_temperature_track = []
    ctth_opaque_track = []
    lon_avhrr_track = []
    lat_avhrr_track = []
    surft_track = []
    bt11micron_track = []
    bt12micron_track = []
    satz_track = []
    lwp_track = []
    cph_track = []

    #idx = [x in range(row_matched.shape[0])]
    lat_avhrr_track = [GeoObj.latitude[row_matched[idx], col_matched[idx]] 
                     for idx in range(row_matched.shape[0])]
    lon_avhrr_track = [GeoObj.longitude[row_matched[idx], col_matched[idx]]
                     for idx in range(row_matched.shape[0])]
    ctype_track = [ctype.cloudtype[row_matched[idx], col_matched[idx]]
                 for idx in range(row_matched.shape[0])]
    ctype_qflag_track = [ctype.qualityflag[row_matched[idx], col_matched[idx]]
                        for idx in range(row_matched.shape[0])]
    if surft != None:
        surft_track = [surft[row_matched[idx], col_matched[idx]]
                       for idx in range(row_matched.shape[0])]
    #b11   
    if dataObj != None:
        temp = [dataObj.channels[3].data[row_matched[idx], col_matched[idx]]
                for idx in range(row_matched.shape[0])] 
        b11_temp =  [(dataObj.channels[3].data[row_matched[idx], 
                                               col_matched[idx]] * 
                      dataObj.channels[3].gain + 
                      dataObj.channels[3].intercept)       
                     for idx in range(row_matched.shape[0])]
        bt11micron_track = np.where(
            np.logical_or(
                np.equal(temp, dataObj.nodata),
                np.equal(temp, dataObj.missing_data)),
            -9, b11_temp)
    #b12
        temp = [dataObj.channels[4].data[row_matched[idx], col_matched[idx]]
                for idx in range(row_matched.shape[0])]
        b12_temp = [(dataObj.channels[4].data[row_matched[idx], 
                                              col_matched[idx]] * 
                     dataObj.channels[4].gain + 
                     dataObj.channels[4].intercept)
                    for idx in range(row_matched.shape[0])]
        bt12micron_track = np.where(
            np.logical_or(
                np.equal(temp, dataObj.nodata),
                np.equal(temp, dataObj.missing_data)),
            -9, b12_temp)

    temp = [AngObj.satz.data[row_matched[idx], col_matched[idx]] 
            for idx in range(row_matched.shape[0])]
    sats_temp = [(AngObj.satz.data[row_matched[idx], col_matched[idx]] * 
                  AngObj.satz.gain + AngObj.satz.intercept)
                 for idx in range(row_matched.shape[0])]
    satz_track = np.where(
        np.logical_or(
            np.equal(temp, AngObj.satz.no_data),
            np.equal(temp, AngObj.satz.missing_data)),
        -9, sats_temp)
    if ctth == None:
        write_log('INFO', "Not extracting ctth")
        pass
    else:
        write_log('INFO', "Extracting ctth along track ")
        temp = [ctth.height[row_matched[idx], col_matched[idx]]
                for idx in range(row_matched.shape[0])]
        hh_temp = [(ctth.height[row_matched[idx], col_matched[idx]] * 
                    ctth.h_gain + ctth.h_intercept)
                   for idx in range(row_matched.shape[0])]
        ctth_height_track = np.where(np.equal(temp, ctth.h_nodata), 
                                     -9, hh_temp)
        temp = [ctth.temperature[row_matched[idx], col_matched[idx]]
                for idx in range(row_matched.shape[0])]
        tt_temp = [(ctth.temperature[row_matched[idx], col_matched[idx]] * 
                    ctth.t_gain + ctth.t_intercept)
                   for idx in range(row_matched.shape[0])]
        ctth_temperature_track = np.where(np.equal(temp, ctth.t_nodata),
                                          -9,tt_temp)
        temp = [ctth.pressure[row_matched[idx], col_matched[idx]]
                for idx in range(row_matched.shape[0])]
        pp_temp = [(ctth.pressure[row_matched[idx], col_matched[idx]] * 
                   ctth.p_gain + ctth.p_intercept)
                   for idx in range(row_matched.shape[0])]
        ctth_pressure_track = np.where(np.equal(temp, ctth.p_nodata), 
                                      -9, pp_temp)
        if (PPS_VALIDATION ):
            is_opaque = np.bitwise_and(np.right_shift(ctth.processingflag, 3), 1)
            ctth_opaque_track = [is_opaque[row_matched[idx], col_matched[idx]]
                                 for idx in range(row_matched.shape[0])]
    #: TODO Do not use fix value -1 but instead something.no_data
    if avhrrLwp != None:
        lwp_temp = [avhrrLwp[row_matched[idx], col_matched[idx]]
                    for idx in range(row_matched.shape[0])]
        lwp_track = np.where(np.equal(lwp_temp, -1), -9, lwp_temp)
    if avhrrCph != None:
        cph_temp = [avhrrCph[row_matched[idx], col_matched[idx]]
                    for idx in range(row_matched.shape[0])]
        cph_track = np.where(np.equal(cph_temp, -1), -9, cph_temp)


    obt.avhrr.latitude = np.array(lat_avhrr_track)
    obt.avhrr.longitude = np.array(lon_avhrr_track)
    obt.avhrr.cloudtype = np.array(ctype_track)
    obt.avhrr.cloudtype_qflag = np.array(ctype_qflag_track)
    if dataObj !=None:
        obt.avhrr.bt11micron = np.array(bt11micron_track)
        obt.avhrr.bt12micron = np.array(bt12micron_track)
    obt.avhrr.satz = np.array(satz_track)
    if ctth:
        obt.avhrr.ctth_height = np.array(ctth_height_track)
        obt.avhrr.ctth_pressure = np.array(ctth_pressure_track)
        obt.avhrr.ctth_temperature = np.array(ctth_temperature_track)
        if (PPS_VALIDATION):
            obt.avhrr.ctth_opaque = np.array(ctth_opaque_track)
    if surft != None:
        obt.avhrr.surftemp = np.array(surft_track)
    if avhrrLwp != None:
        obt.avhrr.lwp = np.array(lwp_track)
    if avhrrCph != None:
        obt.avhrr.cph = np.array(cph_track)
    return obt

# -----------------------------------------------------------------
def match_calipso_avhrr(values, 
                        calipsoObj, imagerGeoObj, imagerObj, 
                        ctype, ctth, cppCph, surft,
                        avhrrAngObj, options, res=resolution):

    import time
    import string
    
    retv = CalipsoAvhrrTrackObject()
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone # Convert from TAI time to UTC in seconds since 1970
    if res == 1:
        lonCalipso = calipsoObj.longitude.ravel()
        latCalipso = calipsoObj.latitude.ravel()
        timeCalipso = calipsoObj.time[::,0].ravel() + dsec
        elevationCalipso = calipsoObj.elevation.ravel()
    if res == 5:
        # Use [:,1] Since 5km data has start, center, and end for each pixel
        lonCalipso = calipsoObj.longitude[:,1].ravel()
        latCalipso = calipsoObj.latitude[:,1].ravel()
        timeCalipso = calipsoObj.time[:,1].ravel() + dsec
        elevationCalipso = calipsoObj.elevation[::,2].ravel()
    
    ndim = lonCalipso.shape[0]
    
    # --------------------------------------------------------------------
    cal,cap = get_calipso_avhrr_linpix(imagerGeoObj,values,lonCalipso,latCalipso,timeCalipso, options)
    # This function (match_calipso_avhrr) could use the MatchMapper object
    # created in map_avhrr() to make things a lot simpler... See usage in
    # amsr_avhrr_match.py
    #from common import map_avhrr
    #cal, cap = map_avhrr(imagerGeoObj, lonCalipso.ravel(), latCalipso.ravel(),
    #                     radius_of_influence=res * .7 * 1e3) # somewhat larger than radius...
    calnan = np.where(cal == NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
    if (PPS_VALIDATION):
        #CCIcloud already have time as array.
        imagerGeoObj = createAvhrrTime(imagerGeoObj, values)
    
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal,cap)]
        avhrr_lines_sec_1970 = np.where(cal != NODATA, imager_time_vector, np.nan)
    else:
        avhrr_lines_sec_1970 = np.where(cal != NODATA, imagerGeoObj.time[cal], np.nan)

#    avhrr_lines_sec_1970 = calnan * DSEC_PER_AVHRR_SCALINE + imagerGeoObj.sec1970_start
    # Find all matching Calipso pixels within +/- sec_timeThr from the AVHRR data
#    pdb.set_trace()
    idx_match = elements_within_range(timeCalipso, avhrr_lines_sec_1970, sec_timeThr) 
    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)
    
    lon_calipso = np.repeat(lonCalipso, idx_match)
    lat_calipso = np.repeat(latCalipso, idx_match)
    # Calipso line,pixel inside AVHRR swath:
    cal_on_avhrr = np.repeat(cal, idx_match)
    cap_on_avhrr = np.repeat(cap, idx_match)
    write_log('INFO', "Start and end times: ",
              time.gmtime(timeCalipso[0]),
              time.gmtime(timeCalipso[ndim-1]))
    
    retv.calipso.sec_1970 = np.repeat(timeCalipso,idx_match)

    retv.calipso.cloud_fraction = np.repeat(calipsoObj.cloud_fraction,idx_match)
    retv.calipso.latitude = np.repeat(latCalipso,idx_match)
    retv.calipso.longitude = np.repeat(lonCalipso,idx_match)
    #for p in range(20):
    #    print "cal time", time.gmtime(timeCalipso[p])
    #    print "cal lat", retv.calipso.latitude[p]
    #    print "cal lon", retv.calipso.longitude[p]
    #    #print "cci time", time.gmtime(avhrr_lines_sec_1970[p])
    #    print "cci lat", imagerGeoObj.latitude[cal_on_avhrr[p],cap_on_avhrr[p]]
    #    print "cci lat", imagerGeoObj.longitude[cal_on_avhrr[p],cap_on_avhrr[p]]

   
    write_log('INFO',"cap_on_avhrr.shape: ",cap_on_avhrr.shape)
    retv.calipso.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.calipso.avhrr_pixnum = cap_on_avhrr.astype('i')
    
    #print "Concatenate arrays..."
    #x = np.concatenate((idx_match,idx_match))
    #for i in range(2,10):
    #    x = np.concatenate((x,idx_match))
    #idx_match_2d = np.reshape(x,(ndim,10))

    write_log('INFO', "Make cloud top and base arrays...")
#    missing_data = -9.9
    #cloud_top = np.repeat(calipsoObj.cloud_top_profile.flat,idx_match_2d.flat)
    #cloud_top = np.where(np.less(cloud_top,0),missing_data,cloud_top)
    #N = cloud_top.flat.shape[0]/10
    #cloud_top = np.reshape(cloud_top,(N,10))
    
    x_fcf = np.repeat(calipsoObj.feature_classification_flags[::,0],idx_match)
    x_ctp = np.repeat(calipsoObj.cloud_top_profile[::,0],idx_match)
    x_cbp = np.repeat(calipsoObj.cloud_base_profile[::,0],idx_match)
    x_cmt = np.repeat(calipsoObj.cloud_mid_temperature[::,0],idx_match)
    
    col_dim = calipsoObj.cloud_mid_temperature.shape[1]
    for i in range(1,col_dim):
        x_fcf = np.concatenate(\
            (x_fcf,np.repeat(calipsoObj.feature_classification_flags[::,i],idx_match)))
        x_ctp = np.concatenate(\
            (x_ctp,np.repeat(calipsoObj.cloud_top_profile[::,i],idx_match)))
        x_cbp = np.concatenate(\
            (x_cbp,np.repeat(calipsoObj.cloud_base_profile[::,i],idx_match)))
        x_cmt = np.concatenate(\
            (x_cmt,np.repeat(calipsoObj.cloud_mid_temperature[::,i],idx_match)))
    N_fcf = x_fcf.shape[0]/col_dim
    retv.calipso.feature_classification_flags = np.reshape(x_fcf,(col_dim,N_fcf)).astype('i')
    N_ctp = x_ctp.shape[0]/col_dim
    retv.calipso.cloud_top_profile = np.reshape(x_ctp,(col_dim,N_ctp)).astype('d')
    N_cbp = x_cbp.shape[0]/col_dim
    retv.calipso.cloud_base_profile = np.reshape(x_cbp,(col_dim,N_cbp)).astype('d')
    N_cmt = x_cmt.shape[0]/col_dim
    retv.calipso.cloud_mid_temperature = np.reshape(x_cmt,(col_dim,N_cmt)).astype('d')
    if res == 5:
        x_od = np.repeat(calipsoObj.optical_depth[::,0],idx_match)
        x_odu = np.repeat(calipsoObj.optical_depth_uncertainty[::,0],idx_match)
        x_ss = np.repeat(calipsoObj.single_shot_cloud_cleared_fraction[::,0],idx_match)
        for i in range(1, col_dim):
            x_od = np.concatenate(\
                (x_od,np.repeat(calipsoObj.optical_depth[::,i],idx_match)))
            x_odu = np.concatenate(\
                (x_odu,np.repeat(calipsoObj.optical_depth_uncertainty[::,i],idx_match)))
            x_ss = np.concatenate(\
                (x_ss,np.repeat(calipsoObj.single_shot_cloud_cleared_fraction[::,i],idx_match)))
        N_od = x_od.shape[0]/col_dim
        retv.calipso.optical_depth = np.reshape(x_od,(col_dim,N_od)).astype('d')
        N_odu = x_odu.shape[0]/col_dim
        retv.calipso.optical_depth_uncertainty = np.reshape(x_odu,(col_dim,N_odu)).astype('d')
        N_ss = x_ss.shape[0]/col_dim
        retv.calipso.single_shot_cloud_cleared_fraction = np.reshape(x_ss,(col_dim,N_ss)).astype('d')
    #cloud_mid_temp = np.repeat(calipsoObj.cloud_mid_temperature.flat,idx_match_2d.flat)
    #cloud_mid_temp = np.where(np.less(cloud_mid_temp,0),missing_data,cloud_mid_temp)
    #cloud_mid_temp = np.reshape(cloud_mid_temp,(N,10))
    #retv.calipso.cloud_mid_temperature = cloud_mid_temp
    
    # IGBP Land Cover:
    retv.calipso.igbp = np.repeat(calipsoObj.igbp.ravel(),idx_match.ravel())

    # NSIDC Ice and Snow Cover:
    retv.calipso.nsidc = np.repeat(calipsoObj.nsidc.ravel(),idx_match.ravel())

    # Elevation is given in km's. Convert to meters:
    retv.calipso.elevation = np.repeat(elevationCalipso.ravel()*1000.0,
                                            idx_match.ravel()).astype('d')

    retv.calipso.number_of_layers_found = np.repeat(\
        calipsoObj.number_of_layers_found.ravel(),idx_match.ravel()).astype('i')
    
    # Time
    if len(imagerGeoObj.time.shape)>1:
        retv.avhrr.sec_1970= [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal_on_avhrr,cap_on_avhrr)]
    else:
        retv.avhrr.sec_1970 = imagerGeoObj.time[cal_on_avhrr]
    retv.diff_sec_1970 = retv.calipso.sec_1970 - retv.avhrr.sec_1970

    min_diff = np.minimum.reduce(retv.diff_sec_1970)
    max_diff = np.maximum.reduce(retv.diff_sec_1970)
    write_log('INFO', "Maximum and minimum time differences in sec (avhrr-calipso): ",
          np.maximum.reduce(retv.diff_sec_1970),np.minimum.reduce(retv.diff_sec_1970))

    write_log('INFO', "AVHRR observation time of first calipso-avhrr match: ",
          time.gmtime(retv.avhrr.sec_1970[0]))
    write_log('INFO', "AVHRR observation time of last calipso-avhrr match: ",
          time.gmtime(retv.avhrr.sec_1970[N_cmt-1]))

    # Make the latitude and pps cloudtype on the calipso track:
    # line and pixel arrays have equal dimensions
    write_log('INFO', "Generate the latitude,cloudtype tracks!")
    
    # -------------------------------------------------------------------------
    # Pick out the data from the track from AVHRR
    retv = avhrr_track_from_matched(retv, imagerGeoObj, imagerObj, avhrrAngObj, 
                                    surft, ctth, ctype, cal_on_avhrr, 
                                    cap_on_avhrr, avhrrCph=cppCph)
    # -------------------------------------------------------------------------    

# for arname, value in retv1.avhrr.all_arrays.items():
#        if arname not in retv.avhrr.all_arrays.keys():
#            pdb.set_trace()
#        print arname
#        if value == None:
#            if retv.avhrr.all_arrays[arname] == None:
#                continue
#            else:
#                pdb.set_trace()
#        if (value==retv.avhrr.all_arrays[arname]).all() != True:
#            pdb.set_trace()
#        if (value.dtype==retv.avhrr.all_arrays[arname].dtype) != True:
#            pdb.set_trace()
#        if (value.shape==retv.avhrr.all_arrays[arname].shape) != True:
#            pdb.set_trace()
#    for arname in retv.avhrr.all_arrays.keys():
#        if arname not in retv1.avhrr.all_arrays.keys():
#            pdb.set_trace()
#    print('all avhrr correct')
    write_log('INFO', "AVHRR-PPS Cloud Type,latitude: shapes = ",
          retv.avhrr.cloudtype.shape,retv.avhrr.latitude.shape)
    ll = []
    for i in range(ndim):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso[i],latCalipso[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso[i],latCalipso[i],idx_match[i])))
    #basename = os.path.basename(ctypefile).split(".h5")[0]
    #values={"satellite":basename.split("_")[-8]}
    #values["year"] = str(basename.split("_")[-7][0:4])
    #values["month"] = str(basename.split("_")[-7][4:6])
    #values["basename"] = string.join(basename.split("_")[0:4],"_")
    data_path = options['data_dir'].format(val_dir=_validation_results_dir, 
                                           satellite=values["satellite"],
                                           resolution=str(RESOLUTION),
                                           year=values["year"],
                                           month=values["month"],
                                           area=AREA)                                            
    if not os.path.exists(data_path):
        write_log('INFO', "Creating datadir: %s"%(data_path ))
        os.makedirs(data_path)
    data_file = options['data_file'].format(resolution=str(RESOLUTION),
                                            basename=values["basename"],
                                            atrain_sat="calipso",
                                            track="track2")
    filename = data_path +  data_file        
    fd = open(filename,"w")
    fd.writelines(ll)
    fd.close()
    ll = []
    for i in range(N_cmt):
        ll.append(("%7.3f  %7.3f  %d\n"%(lon_calipso[i],lat_calipso[i],0)))
    data_file = options['data_file'].format(resolution=str(RESOLUTION),
                                            basename=values["basename"],
                                            atrain_sat="calipso",
                                            track="track_excl")
    filename = data_path + data_file 
    fd = open(filename,"w")
    fd.writelines(ll)
    fd.close()    
    # CALIOP Maximum cloud top in km:
    max_cloud_top_calipso = np.maximum.reduce(retv.calipso.cloud_top_profile.ravel())
    write_log('INFO', "max_cloud_top_calipso: ",max_cloud_top_calipso)
    return retv,min_diff,max_diff

#===============================================================================
# # -----------------------------------------------------
# def select_calipso_inside_avhrr(calipsoObj,cal,dsec,sec1970_start_end,sec_timeThr):
#    import numpy as np
# 
#    sec1970_start,sec1970_end = sec1970_start_end
#    
#    # Select the points inside the avhrr swath:
#    # Allowing for sec_timeThr seconds deviation:
#    if RESOLUTION == 1:
#        idx_time_okay = np.logical_and(np.greater(\
#            calipsoObj.time[:,0],sec1970_start - dsec - sec_timeThr),
#                                   np.less(\
#            calipsoObj.time[:,0],sec1970_end - dsec   + sec_timeThr))
#    elif RESOLUTION == 5:
#        idx_time_okay = np.logical_and(np.greater(\
#            calipsoObj.time[:,1],sec1970_start - dsec - sec_timeThr),
#                                   np.less(\
#            calipsoObj.time[:,1],sec1970_end - dsec   + sec_timeThr)) ##########################################################################################################################################   
#    #pdb.set_trace()
#    #idx_match = np.not_equal(cal,NODATA)        
#    idx_place_okay = np.where(np.not_equal(cal,NODATA),idx_time_okay,False)
#    idx_match = idx_place_okay
#    
#    #idx_match = np.logical_and(np.greater(lin,0),idx_okay[::,0])
#    #idx_match = np.logical_and(idx_match,np.logical_and(np.greater(pix,0),np.less_equal(pix,2048)))
#    #print "Number of matches: ",np.repeat(idx_match,idx_match).shape[0]
# 
#    # Get the PPS Cloud Types matching CALIPSO:
#    #line = np.repeat(lin,idx_match)
#    #line = np.floor(line+0.5).astype('i')
#    #pixel = np.repeat(pix,idx_match)
#    #pixel = np.floor(pixel+0.5).astype('i')
#    #print "Number of matches: ",line.shape[0]
# 
#    return idx_match
#===============================================================================

# -----------------------------------------------------
def get_calipso(filename, res):
    import _pypps_filters
    # Read CALIPSO Lidar (CALIOP) data:
    clobj = read_calipso(filename, res)
    if res == 1:
        lon = clobj.longitude.ravel()
        ndim = lon.shape[0]
        # --------------------------------------------------------------------
        # Derive the calipso cloud fraction using the 
        # cloud height:       
        winsz = 3
        max_height = np.ones(clobj.cloud_top_profile[::, 0].shape) * -9
        for idx in range(clobj.cloud_top_profile.shape[1]):
            max_height = np.maximum(max_height,
                                    clobj.cloud_top_profile[::, idx] * 1000.)
    
        calipso_clmask = np.greater(max_height, 0).astype('d')
        cloud_mask = np.vstack([ calipso_clmask for idx in range(winsz) ])
        #cloud_mask = np.concatenate((calipso_clmask,calipso_clmask))
        #for idx in range(2,winsz): #@UnusedVariable
        #    cloud_mask = np.concatenate((cloud_mask,calipso_clmask))
        #cloud_mask = np.reshape(cloud_mask,(winsz,ndim)).astype('d')
    
        clobj.cloud_fraction = np.zeros((winsz, ndim), 'd')
        cloud_mask_nodata = -1
        _pypps_filters.texture(cloud_mask, clobj.cloud_fraction,
                               winsz, "mean", cloud_mask_nodata)
        clobj.cloud_fraction = clobj.cloud_fraction[winsz/2, ::]
    
    elif res == 5:

#        lonCalipso = calipso.longitude[:,1].ravel()
#        latCalipso = calipso.latitude[:,1].ravel()
        clobj.cloud_fraction = np.where(clobj.cloud_top_profile[:,0] > 0, 1, 0).astype('d')
        # Strange - this will give 0 cloud fraction in points with no data, wouldn't it????/KG

    
    return clobj

# -----------------------------------------------------
def read_calipso(filename, res):
    
    import _pyhl #@UnresolvedImport
    import h5py #@UnresolvedImport
    

#    if res == 5:
#        h5file = h5py.File(filename, 'r')
#        pdb.set_trace()
#        h5file['Horizontal_Averaging']
#        h5file.close()

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
    if res == 5:
        write_log('INFO', "calipso-file %s" % filename)
        c=a.getNode("/Feature_Optical_Depth_532")
        retv.optical_depth=c.data()
        c=a.getNode("/Feature_Optical_Depth_Uncertainty_532")
        retv.optical_depth_uncertainty=c.data()
        c=a.getNode("/Single_Shot_Cloud_Cleared_Fraction")
        retv.single_shot_cloud_cleared_fraction=c.data()

    return retv

# -----------------------------------------------------
def reshapeCalipso(calipsofiles, avhrr, values, timereshape = True, res=resolution):
    import time
    import sys
    
    cal= CalipsoObject()
    if (PPS_VALIDATION):
        avhrr = createAvhrrTime(avhrr, values)
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start

    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    startCalipso = get_calipso(calipsofiles[0], res)
    # Concatenate the data from the different files
    for i in range(len(calipsofiles) - 1):
        newCalipso = get_calipso(calipsofiles[i + 1], res)
        if res == 1:
            cal_start_all = startCalipso.time[:,0] + dsec
            cal_new_all = newCalipso.time[:,0] + dsec
        elif res == 5:
            cal_start_all = startCalipso.time[:,1] + dsec
            cal_new_all = newCalipso.time[:,1] + dsec
        
        if not cal_start_all[0] < cal_new_all[0]:
            write_log('INFO', "calipso files are in the wrong order")
            print("Program calipso.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9)
            
        cal_break = np.argmin(np.abs(cal_start_all - cal_new_all[0])) + 1
        # Concatenate the feature values
        #arname = array name from calipsoObj
        for arname, value in startCalipso.all_arrays.items(): 
            if value != None:
                if value.size != 1:
                    startCalipso.all_arrays[arname] = np.concatenate((value[0:cal_break,...],newCalipso.all_arrays[arname]))

    # Finds Break point
    if res == 1:
        start_break = np.argmin((np.abs((startCalipso.time[:,0] + dsec) - (avhrr_start - sec_timeThr))))
        end_break = np.argmin((np.abs((startCalipso.time[:,0] + dsec) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain    
    if res == 5:
        start_break = np.argmin((np.abs((startCalipso.time[:,1] + dsec) - (avhrr_start - sec_timeThr))))
        end_break = np.argmin((np.abs((startCalipso.time[:,1] + dsec) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain 
    if start_break != 0:
        start_break = start_break - 1 # Minus one to get one extra, just to be certain
    
    if timereshape == True:
        # Cute the feature values
        #arnameca = array name from calipsoObj
        for arnameca, valueca in startCalipso.all_arrays.items(): 
            if valueca != None:
                if valueca.size != 1:
                    cal.all_arrays[arnameca] = valueca[start_break:end_break,...]
                else:
                    cal.all_arrays[arnameca] = valueca
    else:
        cal = startCalipso
        
    if cal.time.shape[0] <= 0:
        write_log('INFO',("No time match, please try with some other CloudSat files"))
        print("Program calipso.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)  
    return cal, start_break, end_break

#****************************************************

def add1kmTo5km(Obj1, Obj5, start_break, end_break):
    retv = CalipsoObject()
    # First check if length of 5 km and 1 km arrays correspond (i.e. 1 km array = 5 times longer array)
    # Here we check the middle time (index 1) out of the three time values given (start, mid, end) for 5 km data
    #pdb.set_trace()
    if (Obj5.utc_time[:,1] == Obj1.utc_time[2::5]).sum() != Obj5.utc_time.shape[0]:
                              
        print("length mismatch")
        pdb.set_trace()

    #First making a preliminary check of the differences in fraction of cloudy calipso columns in 1 km and 5 km data.

    #pdb.set_trace()
    cfc_5km = 0
    cfc_1km = 0
    len_5km = Obj5.utc_time.shape[0]
    len_1km = Obj5.utc_time.shape[0]*5
    for i in range(len_5km):
        if Obj5.number_of_layers_found[i] > 0:
            cfc_5km = cfc_5km + 1
    for i in range(len_1km):
        if Obj1.number_of_layers_found[i] > 0:
            cfc_1km = cfc_1km + 1


    print "*****CHECKING CLOUD FREQUENCY DIFFERENCES IN 1KM AND 5KM DATASETS:"
    print " "
    print "Number of 5 km FOVS: ", len_5km
    print "Number of cloudy 5 km FOVS:", cfc_5km
    print "Cloudy fraction 5 km: ", float(cfc_5km)/float(len_5km)
    print "Number of 1 km FOVS: ", len_1km
    print "Number of cloudy 1 km FOVS:", cfc_1km 
    print "Cloudy fraction 1 km: ", float(cfc_1km)/float(len_1km)
    print " "
    #pdb.set_trace()    

    # Now calculate the cloud fraction in 5 km data from 1 km data (discretized to 0.0, 0.2, 0.4, 0.6, 0.8 and 1.0).
    
    # In addition, if there are cloud layers in 5 km data but nothing in 1 km data, set cloud fraction to 1.0.
    # This latter case represents when very thin cloud layers are being detected over longer distances
    
    # Finally, if there are cloud layers in 1 km data but not in 5 km data, add a layer to 5 km data and set corresponding
    # COT to 1.0. Cloud base and cloud top for this layer is calculated as averages from original levels (max height for
    # top and min height for base if there are more than one layer).This is a pragmatic solution to take care of a
    # weakness or bug in the CALIPSO retrieval of clouds below 4 km


    for i in range(Obj5.utc_time.shape[0]):
        cfc = 0.0
        for j in range(5):
            if Obj1.number_of_layers_found[i*5+j] > 0:
                cfc = cfc + 0.2000
        if ((Obj5.number_of_layers_found[i] > 0) and (cfc < 0.1)):
            cfc = 1.0
        if ((cfc > 0.1) and (Obj5.number_of_layers_found[i] == 0)): #Add missing layer due to CALIPSO processing bug
            cloudtop_sum = 0.0
            cloudbase_sum = 0.0
            cloud_layers = 0
            feature_array_list = []
            for j in range(5):
                if Obj1.number_of_layers_found[i*5+j] != 0:
                    for k in range(Obj1.number_of_layers_found[i*5+j]):
                        cloudtop_sum = cloudtop_sum + Obj1.cloud_top_profile[i,k]
                        cloudbase_sum = cloudbase_sum + Obj1.cloud_base_profile[i,k]
                        cloud_layers = cloud_layers + 1
                        feature_array_list.append(Obj1.feature_classification_flags[i, k])
            Obj5.number_of_layers_found[i] = 1
            Obj5.cloud_top_profile[i, 0] = cloudtop_sum/cloud_layers
            Obj5.cloud_base_profile[i, 0] = cloudbase_sum/cloud_layers
            Obj5.optical_depth[i, 0] = 1.0 #Just put it safely away from the thinnest cloud layers - the best we can do!
            # Obj5.feature_classification_flags[i, 0] = 22218 if assuming like below:
            # cloud, low quality, water phase, low quality, low broken cumulus, confident, 1 km horizontal averaging)
            feature_array = np.asarray(feature_array_list)
            Obj5.feature_classification_flags[i, 0] = np.median(feature_array[:]) # However, let's take the median value
            Obj5.single_shot_cloud_cleared_fraction[i] = 0.0 # Just put any value, we will not use it! 
            
        if Obj5.cloud_fraction[i] >= 0.0:
            Obj5.cloud_fraction[i]=cfc

    # Cute the feature values
    #arnameca = array name from calipsoObj
    for arnameca, valueca in Obj5.all_arrays.items(): 
        if valueca != None:
            if valueca.size != 1:
                retv.all_arrays[arnameca] = valueca[start_break:end_break,...]
            else:
                retv.all_arrays[arnameca] = valueca
    return retv
    
def use5km_remove_thin_clouds_from_1km(Obj1, Obj5, start_break, end_break):
    retv = CalipsoObject()
    if (Obj5.utc_time[:,1] == Obj1.utc_time[2::5]).sum() != Obj5.utc_time.shape[0]:
        write_log('WARNING', "length mismatch")
        pdb.set_trace()
    for pixel in range(Obj5.utc_time.shape[0]):    
        cloud_max_top = np.max(Obj5.cloud_top_profile[pixel, 0:10])
        if cloud_max_top ==-9999:
            continue
        else:
            cloud_top_max = int(round(1000*cloud_max_top))
        height_profile = 0.001*np.array(range(cloud_top_max, -1, -1))
        optical_thickness = np.zeros(height_profile.shape)
        for lay in range(Obj5.number_of_layers_found[pixel]): 
            cloud_at_these_height_index = np.logical_and(
                Obj5.cloud_top_profile[pixel, lay]>= height_profile, 
                height_profile>=Obj5.cloud_base_profile[pixel, lay])
            eye_this_cloud = np.where(cloud_at_these_height_index ,  1, 0)
            number_of_cloud_boxes = sum(eye_this_cloud)         
            if number_of_cloud_boxes == 0:
                write_log('WARNING', "cloud has no depth!!")
             
            optical_thickness_this_layer = (
                eye_this_cloud*
                Obj5.optical_depth[pixel, lay]*
                1.0/number_of_cloud_boxes)             
            if abs(np.sum(optical_thickness_this_layer) - 
                   Obj5.optical_depth[pixel, lay])>0.001:
                write_log('WARNING', "The sum of the optical thickness profile is "
                       "not the same as total optical thickness of the cloud!!")
             
            optical_thickness = optical_thickness + optical_thickness_this_layer

        optical_thickness_profile = np.cumsum(optical_thickness)
        ok_and_higher_heights = np.where(
            optical_thickness_profile <= OPTICAL_DETECTION_LIMIT, 
            height_profile, cloud_max_top)
        height_limit1 = np.min(ok_and_higher_heights)
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ###Rolles suggestion, sort out pixels where the comparison is still bad
        ### Tested, but not yet decided to use, set sort_put var to True to use
        sort_out_pixels_that_we_have_no_good_truth_for = False
        cloud_tops = []
        for pixel_1km in range(pixel*5, pixel*5+5, 1):
            cloud_top_max = np.max(Obj1.cloud_top_profile[pixel_1km, :])
            if cloud_top_max > 0:
                cloud_tops.append(cloud_top_max)
        if len (cloud_tops)>0:
            optical_depth_var_approx = max(cloud_tops)-min(cloud_tops)
        else:
            optical_depth_var_approx = -9
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        for pixel_1km in range(pixel*5, pixel*5+5, 1):                          
            for lay in range(Obj1.number_of_layers_found[pixel_1km]-1, -1, -1):
                #print optical_depth_var_approx
                if  optical_depth_var_approx > 0.5 and sort_out_pixels_that_we_have_no_good_truth_for:
                    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Obj1.cloud_top_profile[pixel_1km, lay] = -9999 
                    Obj1.cloud_base_profile[pixel_1km, lay] = -9999 
                    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                #Remove all layers of clouds if total optical thickness is to low:    
                elif (np.max(optical_thickness_profile)< OPTICAL_DETECTION_LIMIT and 
                      EXCLUDE_CALIPSO_PIXEL_IF_TOTAL_OPTICAL_THICKNESS_TO_LOW):
                    Obj1.cloud_top_profile[pixel_1km, lay] = -9999 
                    Obj1.cloud_base_profile[pixel_1km, lay] = -9999  
                elif   (Obj5.optical_depth[pixel, 0]< OPTICAL_DETECTION_LIMIT and  #top layer very thin
                        np.max(optical_thickness_profile)> Obj5.optical_depth[pixel, 0]  and #is multilayer
                        EXCLUDE_MULTILAYER_IF_TOO_THIN_TOP_LAYER):
                    #relative thin top layer, total optical thickness thicker than limit
                    Obj1.cloud_top_profile[pixel_1km, lay] = -9999 
                    Obj1.cloud_base_profile[pixel_1km, lay] = -9999   
                elif   (np.max(optical_thickness_profile)> Obj5.optical_depth[pixel, 0]  and #is multilayer
                        EXCLUDE_ALL_MULTILAYER):
                    Obj1.cloud_top_profile[pixel_1km, lay] = -9999 
                    Obj1.cloud_base_profile[pixel_1km, lay] = -9999   
                elif   (Obj1.cloud_top_profile[pixel_1km, 0]-Obj1.cloud_base_profile[pixel_1km, 0]>1  and #is multilayer
                        EXCLUDE_GEOMETRICALLY_THICK):
                    #relative thin top layer, total optical thickness thicker than limit
                    Obj1.cloud_top_profile[pixel_1km, lay] = -9999 
                    Obj1.cloud_base_profile[pixel_1km, lay] = -9999       
                elif height_limit1 < Obj1.cloud_top_profile[pixel_1km, lay]:
                    #cut cloud at limit or at base of cloud
                    Obj1.cloud_top_profile[pixel_1km, lay] =  max(
                        height_limit1, 
                        Obj1.cloud_base_profile[pixel_1km, lay]+0.1)
                              


#save removed clouds and heights so they can be plotted (yellow in figure as clouds calipso sees but npp can't see)

    for arnameca, valueca in Obj1.all_arrays.items(): 
        if valueca != None:
            if valueca.size != 1:
                retv.all_arrays[arnameca] = valueca[start_break:end_break,...]
            else:
                retv.all_arrays[arnameca] = valueca
    return retv
    
# -----------------------------------------------------
if __name__ == "__main__":
    # Testing:
    import string
    import epshdf #@UnresolvedImport
    import pps_io #@UnresolvedImport
    import calipso_avhrr_matchup #@UnresolvedImport
    import time
    
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
    write_log("INFO","Read AVHRR geolocation data") #@UndefinedVariable
    avhrrGeoObj = pps_io.readAvhrrGeoData(avhrrfile)

    # Read PPS Cloud Type data
    write_log("INFO","Read PPS Cloud Type") #@UndefinedVariable
    ctype = epshdf.read_cloudtype(ctypefile,1,1,0)
    ctth = epshdf.read_cloudtop(ctthfile,1,1,1,0,1)
    
    # --------------------------------------------------------------------
    write_log("INFO","Read CALIPSO data") #@UndefinedVariable
    # Read CALIPSO Lidar (CALIOP) data:
    calipso = get_calipso(calipsofile)

    lonCalipso = calipso.longitude.ravel()
    latCalipso = calipso.latitude.ravel()

    # Calculations with AAPP in ERROR!!! Fixme, Ad 2007-09-19
    #lin,pix = avhrr_linepix_from_lonlat_aapp(lonCalipso,latCalipso,avhrrGeoObj,platform,norbit,yyyymmdd)

    caliop_height = []
    caliop_base = []
    caliop_max_height = np.ones(calipso.cloud_top_profile[::,0].shape)*-9
    for i in range(10):
        hh = np.where(np.greater(calipso.cloud_top_profile[::,i],-9),
                           calipso.cloud_top_profile[::,i] * 1000.,-9)
        caliop_max_height = np.maximum(caliop_max_height,
                                            calipso.cloud_top_profile[::,i] * 1000.)
        caliop_height.append(hh)
        bb = np.where(np.greater(calipso.cloud_base_profile[::,i],-9),
                           calipso.cloud_base_profile[::,i] * 1000.,-9)
        caliop_base.append(bb)

    x = np.repeat(calipso.number_of_layers_found.ravel(),
                       np.greater(calipso.number_of_layers_found.ravel(),0))
    print "Number of points with more than 0 layers: ",x.shape[0]
    
    cal_data_ok = np.greater(caliop_max_height,-9.)

    # Testing...
    caObj = calipso_avhrr_matchup.getCaliopAvhrrMatch(avhrrfile,calipsofile,ctypefile,ctthfile)
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    print "Original: ",calipso.time[16203,0]+dsec
    print "Matchup:  ",caObj.calipso.sec_1970[3421]
    print calipso.cloud_top_profile[16203]
    print caObj.calipso.cloud_top_profile[::,3421]
    
    
