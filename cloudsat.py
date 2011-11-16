#HISTORY CHANGES BY ANKE TETZLAFF#

#080430:
# Could not run in the Arctic using the tmpaid creation
# Used hard coded area 'arctic_super_5010' instead

#080416: Got error message in line 232 and 234:
#         "only rank-0 arrays can be converted to Python scalars";
#         changed start_sec1970 and end_sec1970 to rank0 array 
#         by setting start_sec1970[0] and end_sec1970[0]

import pdb
from pps_error_messages import write_log

#from calipso import * #@UnusedWildImport
from config import AREA, SUB_DIR, DATA_DIR, sec_timeThr, RESOLUTION, \
    NODATA, NLINES, SWATHWD
from common import MatchupError, elements_within_range
from calipso import DataObject, ppsAvhrrObject, define_pcs, writeCoverage,\
    createAvhrrTime, avhrr_track_from_matched

COVERAGE_DIR = "%s/%skm/%s"%(SUB_DIR,RESOLUTION,AREA)
 
class CloudsatObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
                            'longitude': None,
                            'latitude': None,
                            'avhrr_linnum': None,
                            'avhrr_pixnum': None,
                            'elevation': None,
                            'Profile_time': None,
                            'sec_1970': None,
                            'sec1970': None,
                            'TAI_start': None,
                            'Temp_min_mixph_K': None,
                            'Temp_max_mixph_K': None,
                            # The data:
                            'CPR_Cloud_mask': None,
                            'CPR_Echo_Top': None,
                            'Clutter_reduction_flag': None,
                            'Data_quality': None,
                            'Data_targetID': None,
                            'Gaseous_Attenuation': None,
                            'MODIS_Cloud_Fraction': None,
                            'MODIS_cloud_flag': None,
                            'Radar_Reflectivity': None,
                            'Height': None,
                            'SigmaZero': None,
                            'SurfaceHeightBin': None,
                            'SurfaceHeightBin_fraction': None,
                            'sem_NoiseFloor': None,
                            'sem_NoiseFloorVar': None,
                            'sem_NoiseGate': None,
                            'RVOD_liq_water_path': None,
                            'RVOD_liq_water_path_uncertainty': None,
                            'RVOD_ice_water_path': None,
                            'RVOD_ice_water_path_uncertainty': None,
                            'LO_RVOD_liquid_water_path': None,
                            'LO_RVOD_liquid_water_path_uncertainty': None,
                            'IO_RVOD_ice_water_path': None,
                            'IO_RVOD_ice_water_path_uncertainty': None,
                            'RVOD_liq_water_content': None,
                            'RVOD_liq_water_content_uncertainty': None,
                            'RVOD_ice_water_content': None,
                            'RVOD_ice_water_content_uncertainty': None,
                            'LO_RVOD_liquid_water_content': None,
                            'LO_RVOD_liquid_water_content_uncertainty': None,
                            'IO_RVOD_ice_water_content': None,
                            'IO_RVOD_ice_water_content_uncertainty': None
                           }


class CloudsatAvhrrTrackObject:
    def __init__(self):
        self.avhrr=ppsAvhrrObject()
        self.cloudsat=CloudsatObject()
        self.diff_sec_1970=None


def duplicate_names(cloudsatObj):
    # For some reason these ones have two names
    cloudsatObj.echo_top = cloudsatObj.CPR_Echo_Top
    cloudsatObj.sec_1970 = cloudsatObj.sec1970
    cloudsatObj.cloud_mask = cloudsatObj.CPR_Cloud_mask

# ----------------------------------------
def readCloudsatAvhrrMatchObj(filename):
    import h5py #@UnresolvedImport
    
    retv = CloudsatAvhrrTrackObject()
    
    h5file = h5py.File(filename, 'r')
    for group, data_obj in [(h5file['/cloudsat'], retv.cloudsat),
                            (h5file['/avhrr'], retv.avhrr)]:
        for dataset in group.keys():        
            if dataset in data_obj.all_arrays.keys():
                data_obj.all_arrays[dataset] = group[dataset].value

    duplicate_names(retv.cloudsat)
    
    retv.diff_sec_1970 = h5file['diff_sec_1970'].value

    h5file.close()

    return retv


# ----------------------------------------
def writeCloudsatAvhrrMatchObj(filename,cl_obj):
    from common import write_match_objects
    groups = {'cloudsat': cl_obj.cloudsat.all_arrays,
              'avhrr': cl_obj.avhrr.all_arrays}
    write_match_objects(filename, cl_obj.diff_sec_1970, groups)
    
    status = 1
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
#    import _pypps_filters
    #import numpy
    import time

    # Read CLOUDSAT Radar data:
    cloudsat = read_cloudsat(filename)

    if RESOLUTION == 1:
        lonCloudsat = cloudsat.longitude.ravel()
#        latCloudsat = cloudsat.latitude.ravel()
        timeCloudsat = cloudsat.Profile_time.ravel()
    elif RESOLUTION == 5:
        lonCloudsat = cloudsat.longitude[:,1].ravel()
#        latCloudsat = cloudsat.latitude[:,1].ravel()
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
    import h5py #@UnresolvedImport
    from config import CLOUDSAT_TYPE

    def get_data(dataset):
        type_name = dataset.value.dtype.names
        try:
            data = dataset.value[type_name[0]]
        except TypeError:
            data = dataset.value
        # Convert 1-dimensional matrices to 1-d arrays
        if len(data.shape) == 2:
            if data.shape[1] == 1:
                return data[:, 0]
            elif data.shape[0] == 1:
                return data[0, :]
        return data

    retv = CloudsatObject()
    h5file = h5py.File(filename, 'r')
    root="2B-" + CLOUDSAT_TYPE
    for group in ['Geolocation Fields', 'Data Fields']:
        tempG = h5file["%s/%s" % (root, group)]
        for dataset in tempG.keys():
            if dataset in retv.all_arrays.keys():
                retv.all_arrays[dataset] = get_data(tempG[dataset])
            elif dataset.lower() in retv.all_arrays.keys():
                retv.all_arrays[dataset.lower()] = get_data(tempG[dataset])
            elif dataset == 'DEM_elevation':
                retv.all_arrays['elevation'] = get_data(tempG[dataset])
            elif dataset == 'Sigma-Zero':           
                retv.all_arrays['SigmaZero'] = get_data(tempG[dataset])
    h5file.close()

    return retv
# --------------------------------------------
def get_cloudsat_avhrr_linpix(avhrrIn,avhrrname,lon,lat,clTime):
    import numpy
    import os
    tmppcs="tmpproj"
    define_pcs(tmppcs, "Plate Caree, central meridian at 15E",
               ['proj=eqc','ellps=bessel', 'lon_0=15'])
    orbittime =  os.path.basename(avhrrname).split("_")[1:3]
    
    startline=0
#    Inside=0
#    HasEncounteredMatch=0
    i=0    
    if not os.path.exists(COVERAGE_DIR):
        os.makedirs(COVERAGE_DIR)
    while startline < avhrrIn.longitude.shape[0]:
        
        write_log("INFO","Calling get_cloudsat_avhrr_linpix: start-line = ",startline)
        tmpaid = "tmparea_%d"%i
        endline = startline + NLINES
        coverage_filename = "%s/coverage_avhrr_cloudsat_matchup_%s_%s_%s_%.5d_%.5d_%s.h5"%\
                            (COVERAGE_DIR,avhrrIn.satellite,orbittime[0],orbittime[1],
                             startline,endline,tmpaid)
        write_log("INFO","Coverage filename = ",coverage_filename)
        cal,cap,ok = get_cloudsat_avhrr_linpix_segment(avhrrIn,lon,lat,clTime,
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
    import area
    import pps_gisdata
    
    ndim = lon.shape[0]
    
    if avhrrIn.longitude.shape[0] > lines[1]:
        lines_end = lines[1]
    else:
        lines_end = avhrrIn.longitude.shape[0]
    lines_start = lines[0]
    # TODO: Do not use swath_with. Use Shape[1] instead. /Erik
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
    cov = _satproj.create_coverage(areaObj,lonarr,latarr,1) #@UndefinedVariable
    pdb.set_trace()
    print covfilename
    writeCoverage(cov,covfilename,"satproj",AREA)
    mapped_line = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,linearr,NODATA) #@UndefinedVariable
    mapped_pixel = _satproj.project(cov.coverage,cov.rowidx,cov.colidx,pixelarr,NODATA) #@UndefinedVariable
    
    pdb.set_trace()
    write_log("INFO","Go through cloudsat track:")
    cloudsat_avhrr_line = []
    cloudsat_avhrr_pixel = []
#    cloudsat_avhrr_line_time = []
#    cloudsat_avhrr_pixel_time = []
    for i in range(ndim):
        xy_tup=pps_gisdata.lonlat2xy(AREA,lon[i],lat[i])
        x,y=int(xy_tup[0]+0.5),int(xy_tup[1]+0.5)
        dimx=mapped_line.shape[1]#1002#5010
        dimy=mapped_line.shape[0]#1002#5010
        if(x < dimx and x >= 0 and y >= 0 and y < dimy):
            cloudsat_avhrr_line.append(mapped_line[y,x])
            cloudsat_avhrr_pixel.append(mapped_pixel[y,x])
#            cloudsat_avhrr_line_time.append(-9)
#            cloudsat_avhrr_pixel_time.append(-9)
        else:
            cloudsat_avhrr_line.append(-9)
            cloudsat_avhrr_pixel.append(-9)
#            cloudsat_avhrr_line_time.append(-9)
#            cloudsat_avhrr_pixel_time.append(-9)
    cloudsat_avhrr_line = numpy.array(cloudsat_avhrr_line)
    cloudsat_avhrr_pixel = numpy.array(cloudsat_avhrr_pixel)
    x=numpy.repeat(cloudsat_avhrr_line, numpy.not_equal(cloudsat_avhrr_line,-9))
#    cloudsat_avhrr_line_time = numpy.array(cloudsat_avhrr_line_time)
#    cloudsat_avhrr_pixel_time = numpy.array(cloudsat_avhrr_pixel_time)
    # Control the time diference
#    match_cloudsat_points = numpy.where(numpy.not_equal(cloudsat_avhrr_line,-9))
#    avhrr_time = (cloudsat_avhrr_line[match_cloudsat_points] * DSEC_PER_AVHRR_SCALINE) + avhrrIn.sec1970_start
#    cl_time = cltime[match_cloudsat_points]
#    time_diff = avhrr_time-cl_time
    #max_time_diff_allowed = 50*60 #Based on that a lap is 102 min
#    max_time_diff_allowed = sec_timeThr
#    time_match = numpy.where(abs(time_diff)<max_time_diff_allowed)
#    if time_match[0].shape[0]==0:
#        x=numpy.repeat(cloudsat_avhrr_line_time,numpy.not_equal(cloudsat_avhrr_line_time,-9))
#    else:
#        cloudsat_avhrr_line_time[match_cloudsat_points[0][time_match]]= cloudsat_avhrr_line[match_cloudsat_points[0][time_match]]
#        cloudsat_avhrr_pixel_time[match_cloudsat_points[0][time_match]] = cloudsat_avhrr_pixel[match_cloudsat_points[0][time_match]]
#        x=numpy.repeat(cloudsat_avhrr_line_time,numpy.not_equal(cloudsat_avhrr_line_time,-9)) 
    write_log("INFO","Number of matching points = ",x.shape[0])
    if x.shape[0] > 0:
        matchOk = 1
    else:
        matchOk = 0
    pdb.set_trace()
    return cloudsat_avhrr_line, cloudsat_avhrr_pixel, matchOk

# -----------------------------------------------------
def match_cloudsat_avhrr(ctypefile,cloudsatObj,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj):
    import numpy
    import time
    import string
#    import sys
#    import inspect
    import os
    retv = CloudsatAvhrrTrackObject()

    if RESOLUTION ==1:
        lonCloudsat = cloudsatObj.longitude.ravel()
        latCloudsat = cloudsatObj.latitude.ravel()
    
    elif RESOLUTION ==5:
        lonCloudsat = cloudsatObj.longitude[:,1].ravel()
        latCloudsat = cloudsatObj.latitude[:,1].ravel()
        
    timeCloudsat = cloudsatObj.sec1970.ravel()
    ndim = lonCloudsat.shape[0]

    # --------------------------------------------------------------------

    #cal,cap = get_cloudsat_avhrr_linpix(avhrrGeoObj,ctypefile,lonCloudsat,latCloudsat,timeCloudsat)
    from common import map_avhrr
    cal, cap = map_avhrr(avhrrGeoObj, lonCloudsat, latCloudsat,
                         radius_of_influence=RESOLUTION * .7 * 1e3)
    
    calnan = numpy.where(cal == NODATA, numpy.nan, cal)
    if (~numpy.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
    avhrrGeoObj = createAvhrrTime(avhrrGeoObj, ctypefile)
    avhrr_lines_sec_1970 = numpy.where(cal != NODATA, avhrrGeoObj.time[cal], numpy.nan)
    # Find all matching Cloudsat pixels within +/- sec_timeThr from the AVHRR data
    idx_match = elements_within_range(cloudsatObj.sec1970, avhrr_lines_sec_1970, sec_timeThr)
    #            numpy.logical_and(cloudsatObj.sec1970 > avhrr_lines_sec_1970 - sec_timeThr,
    #                              cloudsatObj.sec1970 < avhrr_lines_sec_1970 + sec_timeThr)
    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)  
    
    duplicate_names(cloudsatObj)
    
    #arnamecl = array name from cloudsatObj
    for arnamecl, value in cloudsatObj.all_arrays.items(): 
        if value != None:
            if value.ndim == 0:
                retv.cloudsat.all_arrays[arnamecl] = value.copy()
            elif value.ndim == 1:
                retv.cloudsat.all_arrays[arnamecl] = value.copy()[idx_match,:].astype('d')
            elif value.ndim == 2:
                temp = value.copy()[idx_match,:].astype('d')
                if arnamecl == 'Radar_Reflectivity':
                    temp = numpy.where(numpy.less(temp,0),-9.9,temp)
                retv.cloudsat.all_arrays[arnamecl] = temp.transpose()
    
    # Special because in 5km lon and lat is 2dim and shold just be 1dim
    retv.cloudsat.longitude = numpy.repeat(lonCloudsat, idx_match).astype('d')
    retv.cloudsat.latitude = numpy.repeat(latCloudsat, idx_match).astype('d')
    
    # Cloudsat line,pixel inside AVHRR swath:
    cal_on_avhrr = numpy.repeat(cal, idx_match)
    cap_on_avhrr = numpy.repeat(cap, idx_match)
    retv.cloudsat.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.cloudsat.avhrr_pixnum = cap_on_avhrr.astype('i')
        
    N =    retv.cloudsat.sec_1970.shape[0] 
       
    print "Cloudsat observation time of first cloudsat-avhrr match: ",\
        time.gmtime(retv.cloudsat.sec_1970[0])
    print "Cloudsat observation time of last cloudsat-avhrr match: ",\
        time.gmtime(retv.cloudsat.sec_1970[N-1])
    # Time
    retv.avhrr.sec_1970 = avhrrGeoObj.time[cal_on_avhrr]    
    retv.diff_sec_1970 = retv.cloudsat.sec_1970 - retv.avhrr.sec_1970
    
    min_diff = numpy.minimum.reduce(retv.diff_sec_1970)
    max_diff = numpy.maximum.reduce(retv.diff_sec_1970)
    print "Maximum and minimum time differences in sec (avhrr-cloudsat): ",\
          numpy.maximum.reduce(retv.diff_sec_1970),numpy.minimum.reduce(retv.diff_sec_1970)

    print "AVHRR observation time of first cloudsat-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[0])
    print "AVHRR observation time of last cloudsat-avhrr match: ",\
          time.gmtime(retv.avhrr.sec_1970[N-1])               
                
    # Make the latitude and pps cloudtype on the cloudsat track:
    # line and pixel arrays have equal dimensions
    print "Generate all datatypes (lat,lon,cty,ctth,surft) on the cloudsat track!"
    retv = avhrr_track_from_matched(retv, avhrrGeoObj, avhrrObj, avhrrAngObj, \
                                    surft, ctth, ctype, cal_on_avhrr, cap_on_avhrr)

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
    for i in range(N):
        ll.append(("%7.3f  %7.3f  %d\n"%(retv.cloudsat.longitude[i],retv.cloudsat.latitude[i],0)))
    fd = open("%s/%skm_%s_cloudtype_cloudsat-GEOPROF_track_excl.txt"%(datapath, RESOLUTION, basename),"w")
    fd.writelines(ll)
    fd.close()

    return retv,min_diff,max_diff

#------------------------------------------------------------------------------------------

def reshapeCloudsat(cloudsatfiles, avhrr, avhrrfilename):
#    import time
    import numpy
    import sys
    import inspect
    clsat = CloudsatObject()
    avhrr = createAvhrrTime(avhrr, avhrrfilename)
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    
#    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    
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
        # Concatenate the feature values
        #arname = array name from cloudsatObj
        for arname, value in startCloudsat.all_arrays.items(): 
            if value != None:
                if value.size != 1:
                    startCloudsat.all_arrays[arname] = numpy.concatenate((value[0:clsat_break,...],newCloudsat.all_arrays[arname]))
                
    # Finds Break point
    start_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_start - sec_timeThr))))
    if start_break != 0:
        start_break = start_break - 1 # Minus one to get one extra, just to be certain
    end_break = numpy.argmin((numpy.abs((startCloudsat.sec1970) - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain

    # Cute the feature values
    #arnamecl = array name from cloudsatObj
    for arnamecl, valuecl in startCloudsat.all_arrays.items(): 
        if valuecl != None:
            if valuecl.size != 1:
                clsat.all_arrays[arnamecl] = valuecl[start_break:end_break,...]
            else:
                clsat.all_arrays[arnamecl] = valuecl

    if clsat.Profile_time.shape[0] <= 0:
        print("No time match, please try with some other CloudSat files")
        print("Program cloudsat.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)
    else:
        if clsat.Profile_time.ndim == 1:
            clsat.TAI_start = clsat.Profile_time[0]
        else:
            clsat.TAI_start = clsat.Profile_time[0,0]
        clsat.Profile_time = clsat.Profile_time - clsat.TAI_start        
        clsat.TAI_start = numpy.asarray(clsat.TAI_start)

    return clsat
    
    
    
# -----------------------------------------------------
if __name__ == "__main__":
    # Testing:
#    import string
#    import epshdf
#    import pps_io
    import numpy
    
    from config import CLOUDSAT_DIR
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
