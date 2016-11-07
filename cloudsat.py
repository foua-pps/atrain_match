#HISTORY CHANGES BY ANKE TETZLAFF#

#080430:
# Could not run in the Arctic using the tmpaid creation
# Used hard coded area 'arctic_super_5010' instead

#080416: Got error message in line 232 and 234:
#         "only rank-0 arrays can be converted to Python scalars";
#         changed start_sec1970 and end_sec1970 to rank0 array 
#         by setting start_sec1970[0] and end_sec1970[0]

import pdb #@UnusedImport
import logging
logger = logging.getLogger(__name__)
from matchobject_io import (DataObject,
                            ppsAvhrrObject)                            

#from calipso import * #@UnusedWildImport
from config import (AREA, sec_timeThr, RESOLUTION,
                    NODATA, NLINES, SWATHWD, 
                    _validation_results_dir)
from common import (MatchupError, 
                    elements_within_range)
from extract_imager_along_track import avhrr_track_from_matched


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
                            'IO_RVOD_ice_water_content_uncertainty': None,
                            'RVOD_CWC_status': None
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

# -----------------------------------------------------
def match_cloudsat_avhrr(ctypefile,cloudsatObj,avhrrGeoObj,avhrrObj,ctype,ctth,surft,avhrrAngObj, avhrrLwp,options):
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

    #cal,cap = get_cloudsat_avhrr_linpix(avhrrGeoObj,ctypefile,lonCloudsat,latCloudsat,timeCloudsat, options)

    # This function (match_calipso_avhrr) could use the MatchMapper object
    # created in map_avhrr() to make things a lot simpler... See usage in
    # amsr_avhrr_match.py
    #Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    from common import map_avhrr
    cal, cap = map_avhrr(imagerGeoObj, lonCalipso.ravel(), latCalipso.ravel(),
                         radius_of_influence=RESOLUTION*0.7*1000.0) # somewhat larger than radius...

    calnan = numpy.where(cal == NODATA, numpy.nan, cal)
    if (~numpy.isnan(calnan)).sum() == 0:
        raise MatchupError("No matches within region.")
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
                if value.size != 1:
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
                                    surft, ctth, ctype, cal_on_avhrr, cap_on_avhrr, avhrrLwp)

    print "AVHRR-PPS Cloud Type,latitude: shapes = ",\
          retv.avhrr.cloudtype.shape,retv.avhrr.latitude.shape

    ll = []
    for i in range(ndim):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCloudsat[i],latCloudsat[i],idx_match[i])))
    basename = os.path.basename(ctypefile).split(".h5")[0]
    values={"satellite":basename.split("_")[-8]}
    values["year"] = str(basename.split("_")[-7][0:4])
    values["month"] = str(basename.split("_")[-7][4:6])
    values["basename"] = string.join(basename.split("_")[0:4],"_")
    data_path = options['data_dir'].format(val_dir=_validation_results_dir, 
                                           satellite=values["satellite"],
                                           resolution=str(RESOLUTION),
                                           year=values["year"],
                                           month=values["month"],
                                           area=AREA)                                           
    if not os.path.exists(data_path):
        print "Creating datadir: %s"%(data_path )
        os.makedirs(data_path)
    data_file = options['data_file'].format(resolution=str(RESOLUTION),
                                            basename=values["basename"],
                                            atrain_sat="cloudsat-GEOPROF",
                                            track="track2")
    filename = data_path + data_file    
    fd = open(filename,"w")
    fd.writelines(ll)
    fd.close()
    ll = []
    for i in range(N):
        ll.append(("%7.3f  %7.3f  %d\n"%(retv.cloudsat.longitude[i],retv.cloudsat.latitude[i],0)))
    data_file = options['data_file'].format(resolution=str(RESOLUTION),
                                            basename=values["basename"],
                                            atrain_sat="cloudsat-GEOPROF",
                                            track="track_excl")
    filename = data_path + data_file 
    fd = open(filename,"w")
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
    logger.info("Read CLOUDSAT data")
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
