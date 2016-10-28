import pdb #@UnusedImport

import inspect #@UnusedImport
import os #@UnusedImport
import numpy as np

import logging
logger = logging.getLogger(__name__)

from config import (AREA, _validation_results_dir, 
                    sec_timeThr, COMPRESS_LVL, RESOLUTION,
                    NLINES, SWATHWD, NODATA,
                    DO_WRITE_COVERAGE, DO_WRITE_DATA) #@UnusedImport
from common import MatchupError, TimeMatchError, elements_within_range
from config import RESOLUTION as resolution
from config import (OPTICAL_DETECTION_LIMIT,
                    OPTICAL_LIMIT_CLOUD_TOP,
                    ALSO_USE_5KM_FILES,
                    USE_5KM_FILES_TO_FILTER_CALIPSO_DATA,
                    PPS_VALIDATION,
                    IMAGER_INSTRUMENT,
                    PPS_FORMAT_2012_OR_EARLIER)
EXCLUDE_GEOMETRICALLY_THICK = False
import time as tm



from matchobject_io import (CalipsoAvhrrTrackObject,
                            CalipsoObject)


class area_interface:
    pass

class SatProjCov:
    def __init__(self):
        self.coverage=None
        self.colidx=None
        self.rowidx=None 

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

def createAvhrrTime(Obt, values):
    import os #@Reimport
    from config import DSEC_PER_AVHRR_SCALINE
    #import time
    #from datetime import datetime
    import calendar
    import time
    #filename = os.path.basename(filename)
    # Ex.: npp_20120827_2236_04321_satproj_00000_04607_cloudtype.h5
    if IMAGER_INSTRUMENT == 'viirs':
    #if filename.split('_')[0] == 'npp':
        if Obt.sec1970_start < 0: #10800
            logger.warning( 
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

        logger.info("NPP start time :  %s", time.gmtime(Obt.sec1970_start))
        logger.info("NPP end time : %s", time.gmtime(Obt.sec1970_end))
 
    else:
        if Obt.sec1970_end < Obt.sec1970_start:
            """
            In some GAC edition the end time is negative. If so this if statement 
            tries to estimate the endtime depending on the start time plus number of 
            scanlines multiplied with the estimate scan time for the instrument. 
            This estimation is not that correct but what to do?
            """
            Obt.sec1970_end = int(DSEC_PER_AVHRR_SCALINE * Obt.num_of_lines + Obt.sec1970_start)

        datetime=values["date_time"]
        sec1970_start_filename = calendar.timegm(datetime.timetuple())
        diff_filename_infile_time = sec1970_start_filename-Obt.sec1970_start
        diff_hours= abs( diff_filename_infile_time/3600.0  )
        if (diff_hours<13):
            logger.info("Time in file and filename do agree. Difference  %d hours."%diff_hours)
        if (diff_hours>13):
            """
            This if statement takes care of a bug in start and end time, 
            that occurs when a file is cut at midnight
            Former condition needed line number in file name:
            if ((values["ppsfilename"].split('_')[-3] != '00000' and PPS_FORMAT_2012_OR_EARLIER) or
            (values["ppsfilename"].split('_')[-2] != '00000' and not PPS_FORMAT_2012_OR_EARLIER)):
            Now instead check if we aer more than 13 hours off. 
            If we are this is probably the problem, do the check and make sure results are fine afterwards.
            """
            logger.warning("Time in file and filename do not agree! Difference  %d hours.", diff_hours)
            import calendar, time
            timediff = Obt.sec1970_end - Obt.sec1970_start
            old_start = time.gmtime(Obt.sec1970_start + (24 * 3600)) # Adds 24 h to get the next day in new start
            new_start = calendar.timegm(time.strptime('%i %i %i' %(old_start.tm_year, \
                                                                   old_start.tm_mon, \
                                                                   old_start.tm_mday), \
                                                                   '%Y %m %d'))
            Obt.sec1970_start = new_start
            Obt.sec1970_end = new_start + timediff
            diff_filename_infile_time = sec1970_start_filename-Obt.sec1970_start
            diff_hours= abs( diff_filename_infile_time/3600.0)
            if (diff_hours>20):
                logger.error("Time in file and filename do not agree! Difference  %d hours.", diff_hours)
                raise TimeMatchError("Time in file and filename do not agree.")        
        Obt.time = np.linspace(Obt.sec1970_start, Obt.sec1970_end, Obt.num_of_lines)
    return Obt

def calipso_track_from_matched(retv_calipso, calipso, idx_match):
    # Calipso line,pixel inside AVHRR swath:
    for arnameca, valueca in calipso.all_arrays.items(): 
        if valueca != None:
            if valueca.size != 1:
                the_values = valueca[idx_match,...]
                shape_of_data = the_values.shape
                if len(shape_of_data)==1 or shape_of_data[1]==1:
                    #traditionally atrain_match expect arrays to be of chspe (n,) not (n,1)
                    #This is what repeat returns so let it to as before:
                    retv_calipso.all_arrays[arnameca] = np.repeat(np.array(the_values.ravel()),1)
                else:
                    retv_calipso.all_arrays[arnameca] = the_values
            else:
                retv_calipso.all_arrays[arnameca] = valueca
    return retv_calipso

def match_calipso_avhrr(values, 
                        calipsoObj, calipsoObjAerosol, 
                        imagerGeoObj, imagerObj, 
                        ctype, cma, ctth, cppCph, nwp_obj,
                        avhrrAngObj, nwp_segments, options, res=resolution):

    import time
    import string
    
    retv = CalipsoAvhrrTrackObject()
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone # Convert from TAI time to UTC in seconds since 1970
    if res == 1:
        lonCalipso = calipsoObj.longitude.ravel()
        latCalipso = calipsoObj.latitude.ravel()
        timeCalipso_tai = calipsoObj.profile_time_tai[::,0].ravel()
        timeCalipso = calipsoObj.profile_time_tai[::,0].ravel() + dsec
        timeCalipso_utc = calipsoObj.profile_utc_time[::,0].ravel()
        elevationCalipso = calipsoObj.dem_surface_elevation.ravel()
    if res == 5:
        # Use [:,1] Since 5km data has start, center, and end for each pixel
        lonCalipso = calipsoObj.longitude[:,1].ravel()
        latCalipso = calipsoObj.latitude[:,1].ravel()
        timeCalipso_tai = calipsoObj.profile_time_tai[:,1].ravel()
        timeCalipso = calipsoObj.profile_time_tai[:,1].ravel() + dsec
        timeCalipso_utc = calipsoObj.profile_utc_time[:,1].ravel()
        elevationCalipso = calipsoObj.dem_surface_elevation[::,2].ravel()
    
    ndim = lonCalipso.shape[0]
    
    # --------------------------------------------------------------------
    #cal,cap = get_calipso_avhrr_linpix(imagerGeoObj,values,lonCalipso,latCalipso,timeCalipso, options)
    # This function (match_calipso_avhrr) could use the MatchMapper object
    # created in map_avhrr() to make things a lot simpler... See usage in
    # amsr_avhrr_match.py
    #Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    from common import map_avhrr
    cal, cap = map_avhrr(imagerGeoObj, lonCalipso.ravel(), latCalipso.ravel(),
                         radius_of_influence=RESOLUTION*0.7*1000.0) # somewhat larger than radius...
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

    """
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = (9,8))
    ax = fig.add_subplot(111)
    plt.plot(timeCalipso, 'b.')
    plt.plot(avhrr_lines_sec_1970, 'g.')
    plt.show()
    fig.savefig("nina_test_fil.png")
    """

    idx_match = elements_within_range(timeCalipso, avhrr_lines_sec_1970, sec_timeThr) 
    if idx_match.sum() == 0:
        raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)
    retv.calipso = calipso_track_from_matched(retv.calipso, calipsoObj, idx_match)
    cal_on_avhrr = np.repeat(cal, idx_match)
    cap_on_avhrr = np.repeat(cap, idx_match)
    retv.calipso.avhrr_linnum = cal_on_avhrr.astype('i')
    retv.calipso.avhrr_pixnum = cap_on_avhrr.astype('i')
    logger.info("cap_on_avhrr.shape: %s",cap_on_avhrr.shape)

    lon_calipso = np.repeat(lonCalipso, idx_match)
    lat_calipso = np.repeat(latCalipso, idx_match)
    # Calipso line,pixel inside AVHRR swath:
    cal_on_avhrr = np.repeat(cal, idx_match)
    cap_on_avhrr = np.repeat(cap, idx_match)
    logger.info("Start and end times: %s %s",
              time.gmtime(timeCalipso[0]),
              time.gmtime(timeCalipso[ndim-1]))
    
    retv.calipso.sec_1970 = np.repeat(timeCalipso,idx_match)
    retv.calipso.latitude = np.repeat(latCalipso,idx_match)
    retv.calipso.longitude = np.repeat(lonCalipso,idx_match)
    retv.calipso.profile_time_tai = np.repeat(timeCalipso_tai,idx_match)

    # Elevation is given in km's. Convert to meters:
    retv.calipso.elevation = np.repeat(elevationCalipso.ravel()*1000.0,
                                            idx_match.ravel()).astype('d')
    # Time
    if len(imagerGeoObj.time.shape)>1:
        retv.avhrr.sec_1970= [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal_on_avhrr,cap_on_avhrr)]
    else:
        retv.avhrr.sec_1970 = imagerGeoObj.time[cal_on_avhrr]
    retv.diff_sec_1970 = retv.calipso.sec_1970 - retv.avhrr.sec_1970

    min_diff = np.minimum.reduce(retv.diff_sec_1970)
    max_diff = np.maximum.reduce(retv.diff_sec_1970)
    logger.info("Maximum and minimum time differences in sec (avhrr-calipso): %d %d",
          np.maximum.reduce(retv.diff_sec_1970),np.minimum.reduce(retv.diff_sec_1970))
    logger.info("AVHRR observation time of first calipso-avhrr match: %s",
          time.gmtime(retv.avhrr.sec_1970[0]))
    logger.info("AVHRR observation time of last calipso-avhrr match: %s",
          time.gmtime(retv.avhrr.sec_1970[-1]))

    # Make the latitude and pps cloudtype on the calipso track:
    # line and pixel arrays have equal dimensions
    logger.info("Generate the latitude,cloudtype tracks!")
    # -------------------------------------------------------------------------
    # Pick out the data from the track from AVHRR
    from extract_imager_along_track import avhrr_track_from_matched
    retv = avhrr_track_from_matched(retv, imagerGeoObj, imagerObj, avhrrAngObj, 
                                    nwp_obj, ctth, ctype, cma,  cal_on_avhrr, 
                                    cap_on_avhrr, avhrrCph=cppCph, 
                                    nwp_segments=nwp_segments)

    if calipsoObjAerosol is not None:
        retv.calipso_aerosol = calipso_track_from_matched(retv.calipso_aerosol, calipsoObjAerosol, idx_match)

    # -------------------------------------------------------------------------    
    logger.info("AVHRR-PPS Cloud Type,latitude: shapes = %s %s",
          retv.avhrr.cloudtype.shape,retv.avhrr.latitude.shape)
    ll = []
    for i in range(ndim):        
        #ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso[i],latCalipso[i],0)))
        ll.append(("%7.3f  %7.3f  %d\n"%(lonCalipso[i],latCalipso[i],idx_match[i])))
    max_cloud_top_calipso = np.maximum.reduce(retv.calipso.layer_top_altitude.ravel())
    logger.info("max_cloud_top_calipso: %2.1f",max_cloud_top_calipso)
    return retv,min_diff,max_diff

def get_calipso(filename, res, ALAY=False):
    from scipy import ndimage
    # Read CALIPSO Lidar (CALIOP) data:
    clobj = read_calipso(filename, res, ALAY=ALAY)
    if res == 1 and not ALAY:
        lon = clobj.longitude.ravel()
        ndim = lon.shape[0]
        # --------------------------------------------------------------------
        # Derive the calipso cloud fraction using the 
        # cloud height:       
        winsz = 3 #Means a winsz x sinsz KERNEL is used.
        max_height = np.ones(clobj.layer_top_altitude[::, 0].shape) * -9
        for idx in range(clobj.layer_top_altitude.shape[1]):
            max_height = np.maximum(max_height,
                                    clobj.layer_top_altitude[::, idx] * 1000.)
    
        calipso_clmask = np.greater(max_height, 0).astype('d')
        clobj.cloud_fraction = calipso_clmask 
        ##############################################################
        # Replace _pypps_filter (mean over array) with function from scipy.
        # This filtering of single clear/cloud pixels is questionable.
        # Minor investigation (45 scenes npp), shows small decrease in results if removed.
        cloud_fraction_temp =  ndimage.filters.uniform_filter1d(calipso_clmask, size=winsz)
        #don't use filter to set cloudy pixels to clear
        #clobj.cloud_fraction = np.where(
        #    np.logical_and(clobj.cloud_fraction>1,
        #                   cloud_fraction_temp<1.5/winsz),
        #    0,clobj.cloud_fraction)
        clobj.cloud_fraction = np.where(
            np.logical_and(clobj.cloud_fraction<0,
                           cloud_fraction_temp>((winsz-1.5)/winsz)),            
            1,clobj.cloud_fraction)
       ##############################################################
    elif res == 5 and  not ALAY:
        clobj.cloud_fraction = np.where(clobj.layer_top_altitude[:,0] > 0, 1, 0).astype('d')
        # Strange - this will give 0 cloud fraction in points with no data, wouldn't it????/KG
    return clobj


def read_calipso(filename, res, ALAY=False):
    import h5py     
    print "reading file", filename
    scip_these_larger_variables_until_needed = {
        # if any of these are needed just rempve them from the dictionary!
        "Spacecraft_Position": True, #3D-variable
        #2-D variable with second dimension larger than 40:
        "Attenuated_Backscatter_Statistics_1064" : True,
        "Attenuated_Backscatter_Statistics_532" : True,
        "Attenuated_Total_Color_Ratio_Statistics" : True,
        "Volume_Depolarization_Ratio_Statistics" : True,
        "Particulate_Depolarization_Ratio_Statistics" : True,
        "Cirrus_Shape_Parameter" : True,
        "Cirrus_Shape_Parameter_Invalid_Points" : True,
        "Cirrus_Shape_Parameter_Uncertainty" : True
        }
    atrain_match_names = {
        #Add "_tai" to the profile_time name to not forget it is tai time
        "Profile_Time": "profile_time_tai"}
    retv = CalipsoObject()
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        for dataset in h5file.keys():
            if dataset in "metadata_t":
                continue
            if dataset in scip_these_larger_variables_until_needed.keys():
                continue
            name = dataset.lower()
            if dataset in atrain_match_names.keys():
                name = atrain_match_names[dataset]
            data = h5file[dataset].value
            data = np.array(data)
            setattr(retv, name, data) 
        h5file.close()
    return retv  

def discardCalipsoFilesOutsideTimeRange(calipsofiles_list, avhrr, values, res=resolution, ALAY=False):
    import sys
    import time
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    if (PPS_VALIDATION):
        avhrr = createAvhrrTime(avhrr, values)
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    calipso_within_time_range = []
    for i in range(len(calipsofiles_list)):
        current_file = calipsofiles_list[i]
        newCalipso = get_calipso(current_file, res, ALAY=ALAY)
        if res == 1:
            cal_new_all = newCalipso.profile_time_tai[:,0] + dsec
        elif res == 5:
            cal_new_all = newCalipso.profile_time_tai[:,1] + dsec 
        if cal_new_all[0]>avhrr_end or  cal_new_all[-1]<avhrr_start:
            pass
            #print "skipping file %s outside time_limits"%(current_file)
        else:
            print "keeping file %s inside time_limits"%(current_file)
            calipso_within_time_range.append(current_file)
    return calipso_within_time_range

def reshapeCalipso(calipsofiles, res=resolution, ALAY=False):
    import sys
    #concatenate and reshape calipso files
    cal= CalipsoObject()
    startCalipso = get_calipso(calipsofiles[0], res, ALAY=ALAY)
    # Concatenate the data from the different files
    for i in range(len(calipsofiles) - 1):
        newCalipso = get_calipso(calipsofiles[i + 1], res, ALAY=ALAY)
        if res == 1:
            cal_start_all = startCalipso.profile_time_tai[:,0] 
            cal_new_all = newCalipso.profile_time_tai[:,0] 
        elif res == 5:
            cal_start_all = startCalipso.profile_time_tai[:,1] 
            cal_new_all = newCalipso.profile_time_tai[:,1] 
        if not cal_start_all[0] < cal_new_all[0]:
            logger.info("calipso files are in the wrong order")
            print("Program calipso.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9)            
        cal_break = len(cal_start_all) # np.argmin(np.abs(cal_start_all - cal_new_all[0])) + 1
        print "same number", cal_break, len(cal_start_all)
        # Concatenate the feature values
        #arname = array name from calipsoObj
        for arname, value in startCalipso.all_arrays.items(): 
            if np.size(value)>1 or value != None:
                if value.size != 1:
                    startCalipso.all_arrays[arname] = np.concatenate((value[0:cal_break,...], 
                                                                      newCalipso.all_arrays[arname])) 
    cal = startCalipso        
    if cal.profile_time_tai.shape[0] <= 0:
        logger.info(("No time match, please try with some other Calipso files"))
        print("Program calipso.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)  
    return cal

def find_break_points(startCalipso, avhrr, values, res=resolution):
    """
    Find the start and end point where calipso and avhrr matches is within 
    time limits.
    """
    import time
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    if (PPS_VALIDATION):
        avhrr = createAvhrrTime(avhrr, values)
    avhrr_end = avhrr.sec1970_end
    avhrr_start = avhrr.sec1970_start
    # Finds Break point
    if res == 1:
        start_break = np.argmin((np.abs((startCalipso.profile_time_tai[:,0] + dsec) 
                                        - (avhrr_start - sec_timeThr))))
        end_break = np.argmin((np.abs((startCalipso.profile_time_tai[:,0] + dsec) 
                                      - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain    
    if res == 5:
        start_break = np.argmin((np.abs((startCalipso.profile_time_tai[:,1] + dsec) 
                                        - (avhrr_start - sec_timeThr))))
        end_break = np.argmin((np.abs((startCalipso.profile_time_tai[:,1] + dsec) 
                                      - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain 
    if start_break != 0:
        start_break = start_break - 1 # Minus one to get one extra, just to be certain
    return start_break, end_break

def time_reshape_calipso(startCalipso, avhrr=None, values=None, 
                         start_break=None, end_break=None):
    """
    Cut the calipso data at the point where matches with avhrr is within 
    time limits.
    """
    if start_break==None or end_break==None:
        if avhrr==None or values ==None:
            logger.info(("Call function with either avhrr and values or"
                              "start_brak and end_break!"))
            print("Program calipso.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit(-9) 
        start_break, end_break = find_break_points(startCalipso, avhrr, 
                                                   values, res=resolution)
    # Cut the feature values
    #arnameca = array name from calipsoObj
    cal= CalipsoObject()
    for arnameca, valueca in startCalipso.all_arrays.items(): 
        if valueca != None:
            if valueca.size != 1:
                cal.all_arrays[arnameca] = valueca[start_break:end_break,...]
            else:
                cal.all_arrays[arnameca] = valueca
    return cal

def adjust5kmTo1kmresolution(calipso5km):
    logger.info("Repeat 5km calipso data to fit 1km resoluiton")
    calipso= CalipsoObject()
    for arname, value in calipso5km.all_arrays.items(): 
        if np.size(value)>1 or value != None:
            if value.size != 1:
                the_shape = value.shape
                new_values = np.repeat(value,5,axis=0)
                calipso.all_arrays[arname] = new_values                    
    return calipso 

def add5kmVariablesTo1kmresolution(calipso1km, calipso5km):
    logger.info("Repeat 5km calipso data to fit 1km resoluiton")
    for variable_5km  in [ "column_optical_depth_aerosols_1064",
                           "column_optical_depth_aerosols_532",
                           "column_optical_depth_aerosols_uncertainty_1064",
                           "column_optical_depth_aerosols_uncertainty_532",
                           "column_optical_depth_cloud_532",
                           "column_optical_depth_cloud_uncertainty_532",
                           #"feature_optical_depth_532",
                           #"feature_optical_depth_uncertainty_532"
                       ]:                
        data = getattr(calipso5km, variable_5km)
        new_data = np.repeat(data, 5, axis=0)    
        setattr(calipso1km, variable_5km +"_5km", new_data)      
    for variable_5km  in ["feature_optical_depth_532_top_layer_5km",
                          "total_optical_depth_5km"]:                
        data = getattr(calipso5km, variable_5km)
        new_data = np.repeat(data, 5, axis=0)    
        setattr(calipso1km, variable_5km, new_data)        
    return calipso1km 

def add1kmTo5km(Obj1, Obj5, start_break, end_break):
    retv = CalipsoObject()
    # First check if length of 5 km and 1 km arrays correspond (i.e. 1 km array = 5 times longer array)
    # Here we check the middle time (index 1) out of the three time values given (start, mid, end) for 5 km data
    #pdb.set_trace()
    if (Obj5.profile_utc_time[:,1] == Obj1.profile_utc_time[2::5]).sum() != Obj5.profile_utc_time.shape[0]:
                              
        print("length mismatch")
        pdb.set_trace()

    #First making a preliminary check of the differences in fraction of cloudy calipso columns in 1 km and 5 km data.
    #pdb.set_trace()
    cfc_5km = 0
    cfc_1km = 0
    len_5km = Obj5.profile_utc_time.shape[0]
    len_1km = Obj5.profile_utc_time.shape[0]*5
    for i in range(len_5km):
        if Obj5.number_layers_found[i] > 0:
            cfc_5km = cfc_5km + 1
    for i in range(len_1km):
        if Obj1.number_layers_found[i] > 0:
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
    # COT to 1.0. Cloud base and cloud tp for this layer is calculated as averages from original levels (max height for
    # top and min height for base if there are more than one layer).This is a pragmatic solution to take care of a
    # weakness or bug in the CALIPSO retrieval of clouds below 4 km

    for i in range(Obj5.profile_utc_time.shape[0]):
        cfc = 0.0
        for j in range(5):
            if Obj1.number_layers_found[i*5+j] > 0:
                cfc = cfc + 0.2000
        if (Obj5.number_layers_found[i] > 0):
            Obj5.cloud_fraction[i] = 1.0
        if ((cfc > 0.1) and (Obj5.number_layers_found[i] == 0)): #Add missing layer due to CALIPSO processing bug
            cloudtop_sum = 0.0
            cloudbase_sum = 0.0
            cloud_layers = 0
            feature_array_list = []
            for j in range(5):
                if Obj1.number_layers_found[i*5+j] != 0:
                    for k in range(Obj1.number_layers_found[i*5+j]):
                        cloudtop_sum = cloudtop_sum + Obj1.layer_top_altitude[i,k]
                        cloudbase_sum = cloudbase_sum + Obj1.layer_base_altitude[i,k]
                        cloud_layers = cloud_layers + 1
                        feature_array_list.append(Obj1.feature_classification_flags[i, k])
            Obj5.number_layers_found[i] = 1
            Obj5.layer_top_altitude[i, 0] = cloudtop_sum/cloud_layers
            Obj5.layer_base_altitude[i, 0] = cloudbase_sum/cloud_layers
            Obj5.feature_optical_depth_532[i, 0] = 1.0 #Just put it safely away from the thinnest cloud layers - the best we can do!
            # Obj5.feature_classification_flags[i, 0] = 22218 if assuming like below:
            # cloud, low quality, water phase, low quality, low broken cumulus, confident, 1 km horizontal averaging)
            feature_array = np.asarray(feature_array_list)
            Obj5.feature_classification_flags[i, 0] = np.median(feature_array[:]) # However, let's take the median value
            Obj5.single_shot_cloud_cleared_fraction[i] = 0.0 # Just put any value, we will not use it!
            Obj5.cloud_fraction[i] = cfc
    # Cute the feature values
    #arnameca = array name from calipsoObj
    for arnameca, valueca in Obj5.all_arrays.items(): 
        if valueca != None:
            if valueca.size != 1:
                retv.all_arrays[arnameca] = valueca[start_break:end_break,...]
            else:
                retv.all_arrays[arnameca] = valueca
    return retv

def detection_height_from_5km_data(Obj1, Obj5, start_break, end_break):
    retv = CalipsoObject()
    if (Obj5.profile_utc_time[:,1] == Obj1.profile_utc_time[2::5]).sum() != Obj5.profile_utc_time.shape[0]:
        logger.warning("length mismatch")
        #pdb.set_trace()
    Obj1.detection_height_5km = np.ones(Obj1.number_layers_found.shape)*-9                 
    for pixel in range(Obj5.profile_utc_time.shape[0]):
        top = Obj5.layer_top_altitude[pixel, 0]
        base = Obj5.layer_base_altitude[pixel, 0]
        opt_th = Obj5.feature_optical_depth_532[pixel, 0]        
        if base==-9999 or top==-9999 or opt_th==-9999: 
            #can not calculate detection height without data!
            continue            
        pixel_1km_first = 5*pixel
        need_only_to_use_one_layer = False
        if (Obj5.number_layers_found[pixel]==1 or
            (base>=np.max(Obj5.layer_top_altitude[pixel, 1:10])) and
            opt_th>=OPTICAL_LIMIT_CLOUD_TOP):
            need_only_to_use_one_layer = True
        if  need_only_to_use_one_layer:
            #only have one layer or top layer is completely above other layers and thick
            if opt_th <= OPTICAL_LIMIT_CLOUD_TOP:
                #top layer too thin use base of it             
                Obj1.detection_height_5km[pixel_1km_first:pixel_1km_first+5] = base
            else:     
                # filter top layer
                Obj1.detection_height_5km[pixel_1km_first:pixel_1km_first+5] = (
                    base + (top-base)*(opt_th - OPTICAL_LIMIT_CLOUD_TOP)*1.0/opt_th)    
        elif   Obj1.total_optical_depth_5km[pixel_1km_first]<0 <= OPTICAL_LIMIT_CLOUD_TOP:
            bases = Obj5.layer_base_altitude[pixel, 0:10]
            min_base = np.min(bases[bases!=-9999])
            Obj1.detection_height_5km[pixel_1km_first:pixel_1km_first+5] = min_base
        
        else: 
            cloud_max_top = np.max(Obj5.layer_top_altitude[pixel, 0:10])
            cloud_top_max = int(round(1000*cloud_max_top))          
            height_profile = 0.001*np.array(range(cloud_top_max, -1, -1))
            optical_thickness = np.zeros(height_profile.shape)
            for lay in range(Obj5.number_layers_found[pixel]):

            #dont use layers with negative top or base or optical_thickness values
                top = Obj5.layer_top_altitude[pixel, lay]
                base = Obj5.layer_base_altitude[pixel, lay]
                opt_th = Obj5.feature_optical_depth_532[pixel, lay]    
                if (top!=-9999 and base!=-9999 and opt_th!=-9999):
                    cloud_at_these_height_index = np.logical_and(
                        top >= height_profile, 
                        height_profile>=base)
                    eye_this_cloud = np.where(cloud_at_these_height_index ,  1, 0)
                    number_of_cloud_boxes = sum(eye_this_cloud)         
                    if number_of_cloud_boxes == 0 and top>0:
                        cloud_at_these_height_index = np.logical_and(
                            np.logical_and(
                                top  >= height_profile-0.01, 
                                base >=height_profile-0.01),
                            np.logical_and(
                                top  <= height_profile+0.01, 
                                base <= height_profile+0.01))
                        logger.info("Cloud top %.2f base: %.2f "%(top, base))
                        logger.info(" Using height_profile %2.f %2.f"%(
                                np.min(height_profile[cloud_at_these_height_index]),
                            np.max(height_profile[cloud_at_these_height_index])))
                        eye_this_cloud = np.where(cloud_at_these_height_index ,  1, 0)
                        number_of_cloud_boxes = sum(eye_this_cloud)         
                    if number_of_cloud_boxes == 0:
                        logger.warning("cloud has no depth!!")
                    optical_thickness_this_layer = (
                        eye_this_cloud*opt_th*1.0/number_of_cloud_boxes)             
                    if abs(np.sum(optical_thickness_this_layer) - opt_th)>0.001:
                        logger.warning("The sum of the optical thickness profile is "
                                  "not the same as total optical thickness of the cloud!!")             
                    optical_thickness = optical_thickness + optical_thickness_this_layer
            optical_thickness_profile = np.cumsum(optical_thickness)
            ok_and_higher_heights = np.where(
                optical_thickness_profile <= OPTICAL_LIMIT_CLOUD_TOP, 
                height_profile, cloud_max_top)
            height_limit1 = np.min(ok_and_higher_heights)  
            Obj1.detection_height_5km[pixel_1km_first:pixel_1km_first+5] = height_limit1  

    for arnameca, valueca in Obj1.all_arrays.items(): 
        if valueca != None:
            if valueca.size != 1:
                retv.all_arrays[arnameca] = valueca[start_break:end_break,...]
            else:
                retv.all_arrays[arnameca] = valueca
    return retv    

if __name__ == "__main__":
    # Testing:
    import string
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
    logger.info("Read AVHRR geolocation data") #@UndefinedVariable
    avhrrGeoObj = pps_io.readAvhrrGeoData(avhrrfile)

    # Read PPS Cloud Type data
    logger.info("Read PPS Cloud Type") #@UndefinedVariable
    #ctype = read_cloudtype
    #ctth = read_cloudtop
    
    # --------------------------------------------------------------------
    logger.info("Read CALIPSO data") #@UndefinedVariable
    # Read CALIPSO Lidar (CALIOP) data:
    calipso = get_calipso(calipsofile)

    lonCalipso = calipso.longitude.ravel()
    latCalipso = calipso.latitude.ravel()

    # Calculations with AAPP in ERROR!!! Fixme, Ad 2007-09-19
    #lin,pix = avhrr_linepix_from_lonlat_aapp(lonCalipso,latCalipso,avhrrGeoObj,platform,norbit,yyyymmdd)

    caliop_height = []
    caliop_base = []
    caliop_max_height = np.ones(calipso.layer_top_altitude[::,0].shape)*-9
    for i in range(10):
        hh = np.where(np.greater(calipso.layer_top_altitude[::,i],-9),
                           calipso.layer_top_altitude[::,i] * 1000.,-9)
        caliop_max_height = np.maximum(caliop_max_height,
                                            calipso.layer_top_altitude[::,i] * 1000.)
        caliop_height.append(hh)
        bb = np.where(np.greater(calipso.layer_base_altitude[::,i],-9),
                           calipso.layer_base_altitude[::,i] * 1000.,-9)
        caliop_base.append(bb)

    x = np.repeat(calipso.number_layers_found.ravel(),
                       np.greater(calipso.number_layers_found.ravel(),0))
    print "Number of points with more than 0 layers: ",x.shape[0]
    
    cal_data_ok = np.greater(caliop_max_height,-9.)

    # Testing...
    caObj = calipso_avhrr_matchup.getCaliopAvhrrMatch(avhrrfile,calipsofile,ctypefile,ctthfile)
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    print "Original: ",calipso.profile_time_tai[16203,0]+dsec
    print "Matchup:  ",caObj.calipso.sec_1970[3421]
    print calipso.layer_top_altitude[16203]
    print caObj.calipso.layer_top_altitude[3421,::]
    
    
