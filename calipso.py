import pdb
import os
import numpy as np
import logging
logger = logging.getLogger(__name__)

from common import (MatchupError, TimeMatchError, 
                    InputError, ProcessingError, 
                    elements_within_range)
from config import (_validation_results_dir, 
                    sec_timeThr, 
                    COMPRESS_LVL, 
                    RESOLUTION,
                    NODATA,
                    ALSO_USE_5KM_FILES,
                    PPS_VALIDATION,
                    CALIPSO_version4,
                    CALIPSO_version3,
                    CALIPSO_CLOUDY_MIN_CFC,
                    IMAGER_INSTRUMENT,
                    PPS_FORMAT_2012_OR_EARLIER,
                    CALIPSO_REQUIRED)

import time as tm
from matchobject_io import (CalipsoAvhrrTrackObject,
                            CalipsoObject)
from runutils import do_some_logging

def add_validation_ctth_calipso(calipso):
    calipso.validation_height = calipso.layer_top_altitude[:,0].copy()
    calipso.validation_height[calipso.validation_height>=0] = (
        1000.0*calipso.validation_height[calipso.validation_height>=0])
    calipso.validation_height[calipso.validation_height<0] = -9
    return calipso

def calipso_track_from_matched(retv_calipso, calipso, idx_match):
    # Calipso line,pixel inside AVHRR swath:
    for arnameca, valueca in calipso.all_arrays.items(): 
        if valueca is not None:
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
                        caObj, caObjAerosol, 
                        imagerGeoObj, imagerObj, 
                        ctype, cma, ctth, cpp, nwp_obj,
                        avhrrAngObj, nwp_segments, options, res=RESOLUTION):

    from common import map_avhrr
    retv = CalipsoAvhrrTrackObject()
    
    #Nina 20150313 Swithcing to mapping without area as in cpp. Following suggestion from Jakob
    cal, cap = map_avhrr(imagerGeoObj, 
                         caObj.longitude.ravel(),
                         caObj.latitude.ravel(),
                         radius_of_influence=RESOLUTION*0.7*1000.0) # somewhat larger than radius...
    #warn if no matches
    calnan = np.where(cal == NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        if CALIPSO_REQUIRED:
            raise MatchupError("No matches within region.")
        else:
            logger.warning("No matches within region.")
            return None
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal,cap)]
        imager_lines_sec_1970 = np.where(cal != NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != NODATA, imagerGeoObj.time[cal], np.nan)
    idx_match = elements_within_range(caObj.sec_1970, imager_lines_sec_1970, sec_timeThr) 
    if idx_match.sum() == 0:
        if CALIPSO_REQUIRED:
            raise MatchupError("No matches in region within time threshold %d s." % sec_timeThr)
        else:
            logger.warning("No matches in region within time threshold %d s.", sec_timeThr)
            return None    
    retv.calipso = calipso_track_from_matched(retv.calipso, caObj, idx_match)
    
    # Calipso line,pixel inside AVHRR swath:
    retv.calipso.imager_linnum = np.repeat(cal, idx_match).astype('i')
    retv.calipso.imager_pixnum = np.repeat(cap, idx_match).astype('i')

    # Imager time
    retv.avhrr.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.calipso.sec_1970 - retv.avhrr.sec_1970
    do_some_logging(retv, caObj)
    logger.debug("Generate the latitude,cloudtype tracks!")
    from extract_imager_along_track import avhrr_track_from_matched
    retv = avhrr_track_from_matched(retv, imagerGeoObj, imagerObj, avhrrAngObj, 
                                    nwp_obj, ctth, ctype, cma, cpp=cpp, 
                                    nwp_segments=nwp_segments)
    if caObjAerosol is not None:
        retv.calipso_aerosol = calipso_track_from_matched(retv.calipso_aerosol, caObjAerosol, 
                                                          idx_match)
    max_cloud_top_calipso = np.maximum.reduce(retv.calipso.layer_top_altitude.ravel())
    logger.debug("max_cloud_top_calipso: %2.1f",max_cloud_top_calipso)
    return retv

def get_calipso(filename, res, ALAY=False):
    from scipy import ndimage
    # Read CALIPSO Lidar (CALIOP) data:
    cal = read_calipso(filename)
    cal = add_validation_ctth_calipso(cal)
    if res == 1 and not ALAY:
        lon = cal.longitude.ravel()
        # --------------------------------------------------------------------
        # Derive the calipso cloud fraction using the 
        # cloud height:       
        winsz = 3 #Means a winsz x winsz KERNEL is used.   
        calipso_clmask = np.greater_equal(cal.validation_height, 0).astype('d')
        cal.cloud_fraction = calipso_clmask 
        ##############################################################
        # This filtering of single clear/cloud pixels is questionable.
        # Minor investigation (45 scenes npp), shows small decrease in results if removed.
        cloud_fraction_temp =  ndimage.filters.uniform_filter1d(calipso_clmask*1.0, size=winsz)
        #don't use filter to set cloudy pixels to clear
        #If winsz=3: 1clear 2cloudy => cfc = 0.66
        #   winsz=3; 2clear 1cloudy => cfc = 0.33
        cal.cloud_fraction = np.where(
            np.logical_and(cal.cloud_fraction<1.0,
                           cloud_fraction_temp>0.01),            
            cloud_fraction_temp,cal.cloud_fraction)
       ##############################################################
    elif res == 5 and  not ALAY:
        cal.cloud_fraction = np.where(cal.layer_top_altitude[:,0] > 0, 1, 0).astype('d')
        # Strange - this will give 0 cloud fraction in points with no data, wouldn't it????/KG
    return cal


#READING DATA FROM CALIOP:
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
    #Use version 3 "nsidc_surface_type" as name for sea/ice info
    "Snow_Ice_Surface_Type": "nsidc_surface_type",
    #Add "_tai" to the profile_time name to not forget it is tai time           
    "Profile_Time": "profile_time_tai"}

def rearrange_calipso_the_single_shot_info(retv, singleshotdata):

    # Extract number of cloudy single shots (max 15)
    # plus average cloud base and top
    # in 5 km FOV
    name = "ssNumber_Layers_Found"
    data = singleshotdata[name]
    data = np.array(data)
    data_reshaped_15 = data.reshape(-1,15)
    single_shot_cloud_cleared_array = np.sum(data_reshaped_15==0,axis=1) #Number of clear
    single_shot_cloud_cleared_array = 15-single_shot_cloud_cleared_array #Number of cloudy
    name = "Number_cloudy_single_shots" # New name used here
    #pdb.set_trace()
    setattr(retv, name, np.array(single_shot_cloud_cleared_array))
    #We need also average cloud top and cloud base for single_shot clouds
    data = singleshotdata["ssLayer_Base_Altitude"]
    data = np.array(data)
    data_reshaped_5 = data.reshape(-1,5)
    base_array = data_reshaped_5[:,0]
    base_array = base_array.reshape(-1,15)
    base_array = np.where(base_array>0, base_array, 0.0)
    base_mean = np.where(single_shot_cloud_cleared_array>0, 
                         np.divide(np.sum(base_array,axis=1), single_shot_cloud_cleared_array), 
                         -9.0) #Calculate average cloud base
    name = "Average_cloud_base_single_shots"
    setattr(retv, name, base_mean)
    data = singleshotdata["ssLayer_Top_Altitude"]
    data = np.array(data)
    data_reshaped_5 = data.reshape(-1,5)
    top_array = data_reshaped_5[:,0]
    top_array = top_array.reshape(-1,15)
    top_array = np.where(top_array>0, top_array, 0.0)
    top_mean = np.where(single_shot_cloud_cleared_array>0, 
                        np.divide(np.sum(top_array,axis=1), single_shot_cloud_cleared_array),
                        -9.0) #Calculate average cloud top
    name = "Average_cloud_top_single_shots"
    #pdb.set_trace()
    setattr(retv, name, top_mean)
    #extract less information for 1km matching
    if RESOLUTION == 1:
        # create 1km  singel shot cfc
        data = singleshotdata["ssNumber_Layers_Found"]
        data = np.array(data).reshape(-1,3)
        data_reshaped_3 = data.reshape(-1,3)
        single_shot_num_clear_array = np.sum(data_reshaped_3==0,axis=1) #Number of clear
        cfc_single_shots = (3-single_shot_num_clear_array)/3.0 #Number of cloudy
        name = "cfc_single_shots_1km_from_5km_file" # New name used here
        #pdb.set_trace()
        setattr(retv, name, cfc_single_shots)
    return retv

def read_calipso(filename):
    logger.debug("Reading file %s", filename)
    retv = CalipsoObject()
    if filename is not None:
        if "hdf" in filename:
            retv = read_calipso_hdf4(filename, retv)
        else:
            retv = read_calipso_h5(filename, retv)
    #Adopt some variables     
    dsec = tm.mktime((1993,1,1,0,0,0,0,0,0)) - tm.timezone
    #1km
    if retv.profile_time_tai.shape == retv.number_layers_found.shape: 
        retv.latitude = retv.latitude[:,0]
        retv.longitude = retv.longitude [:,0]
        retv.profile_time_tai = retv.profile_time_tai[:,0]
        # Elevation is given in km's. Convert to meters:
        retv.elevation = retv.dem_surface_elevation[:,0]*1000.0 
    #5km
    else:
        retv.latitude = retv.latitude[:,1]
        retv.longitude = retv.longitude[:,1] 
        retv.profile_time_tai = retv.profile_time_tai[:,1]
        # Elevation is given in km's. Convert to meters:
        retv.elevation = retv.dem_surface_elevation[:,2]*1000.0    
    setattr(retv, "sec_1970", retv.profile_time_tai + dsec)
    return retv        

def read_calipso_hdf4(filename, retv):
    from pyhdf.SD import SD, SDC
    from pyhdf.HDF import HDF, HC
    import pyhdf.VS 
    def convert_data(data):
        if len(data.shape) == 2:
            if data.shape[1] == 1:
                return data[:, 0]
            elif data.shape[0] == 1:
                return data[0, :]
        return data
    if filename is not None:
        h4file = SD(filename, SDC.READ)
        datasets = h4file.datasets()
        attributes = h4file.attributes()
        singleshotdata = {}
        for idx, dataset in enumerate(datasets.keys()):
            #non-goups
            if dataset in scip_these_larger_variables_until_needed.keys():        
                continue
            elif dataset[0:8] == 'Surface_':
                continue
            if dataset in ["ssNumber_Layers_Found", 
                           "ssLayer_Base_Altitude", 
                           "ssLayer_Top_Altitude"]:
                singleshotdata[dataset] = h4file.select(dataset).get()
            if dataset[0:2] == "ss":
                #already saved temporarly what we need
                continue            
            name = dataset.lower()
            #print idx, dataset
            if dataset in atrain_match_names.keys():
                name = atrain_match_names[dataset]
            data = np.array(h4file.select(dataset).get())
            setattr(retv, name, data) 
        if "ssNumber_Layers_Found" in singleshotdata.keys():
            # Extract number of cloudy single shots (max 15)
            # plus average cloud base and top
            # in 5 km FOV
            logger.info("Reading single shot information")
            retv = rearrange_calipso_the_single_shot_info(
                retv,
                singleshotdata)
    return retv 

def read_calipso_h5(filename, retv):
    import h5py
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        if "Single_Shot_Detection" in h5file.keys():
            # Extract number of cloudy single shots (max 15)
            # plus average cloud base and top
            # in 5 km FOV
            logger.info("Reading single shot information")
            retv = rearrange_calipso_the_single_shot_info(
                retv, 
                {"ssNumber_Layers_Found": h5file["Single_Shot_Detection/ssNumber_Layers_Found"].value,
                 "ssLayer_Base_Altitude": h5file["Single_Shot_Detection/ssLayer_Base_Altitude"].value,
                 "ssLayer_Top_Altitude": h5file["Single_Shot_Detection/ssLayer_Top_Altitude"].value})
        for dataset in h5file.keys():
            if dataset in ["Lidar_Surface_Detection", # New group V4
                           "Single_Shot_Detection", #Handeled above
                           "metadata_t", 
                           "metadata"]: 
                #skip all groups
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

def discardCalipsoFilesOutsideTimeRange(calipsofiles_list, avhrrGeoObj, values, res=RESOLUTION, ALAY=False):
    avhrr_end = avhrrGeoObj.sec1970_end
    avhrr_start = avhrrGeoObj.sec1970_start
    calipso_within_time_range = []
    for current_file in calipsofiles_list:
        newCalipso = get_calipso(current_file, res, ALAY=ALAY)
        cal_new_all = newCalipso.sec_1970
        if cal_new_all[0]>avhrr_end + sec_timeThr or  cal_new_all[-1] + sec_timeThr<avhrr_start:
            pass
            #print "skipping file %s outside time_limits"%(current_file)
        else:
            logger.debug("Keeping file %s inside time_limits", os.path.basename(current_file))
            calipso_within_time_range.append(current_file)
    return calipso_within_time_range

def reshapeCalipso(calipsofiles, res=RESOLUTION, ALAY=False):
    #concatenate and reshape calipso files
    startCalipso = get_calipso(calipsofiles[0], res, ALAY=ALAY)
    # Concatenate the data from the different files
    for i in range(len(calipsofiles) - 1):
        newCalipso = get_calipso(calipsofiles[i + 1], res, ALAY=ALAY)
        cal_start_all = startCalipso.profile_time_tai 
        cal_new_all = newCalipso.profile_time_tai      
        if not cal_start_all[0] < cal_new_all[0]:
            raise ProcessingError("Calipso files are in the wrong order!")
        # Concatenate the feature values
        for arname, value in startCalipso.all_arrays.items(): 
            if value is not None:
                startCalipso.all_arrays[arname] = np.concatenate((value[0:,...], 
                                                                  newCalipso.all_arrays[arname]))          
    return startCalipso 

def find_break_points(startCalipso, avhrrGeoObj):
    """
    Find the start and end point where calipso and avhrr matches is within 
    time limits.
    """
    avhrr_end = avhrrGeoObj.sec1970_end
    avhrr_start = avhrrGeoObj.sec1970_start
    # Finds Break point
    start_break = np.argmin((np.abs((startCalipso.sec_1970) 
                                    - (avhrr_start - sec_timeThr))))
    end_break = np.argmin((np.abs((startCalipso.sec_1970) 
                                  - (avhrr_end + sec_timeThr)))) + 2    # Plus two to get one extra, just to be certain    
    if start_break != 0:
        start_break = start_break - 1 # Minus one to get one extra, just to be certain
    return start_break, end_break

def time_reshape_calipso(startCalipso,
                         start_break, end_break):
    """
    Cut the calipso data at the point where matches with avhrr is within 
    time limits.
    """
    # Cut the feature values
    #arnameca = array name from caObj
    cal = CalipsoObject()
    for arnameca, valueca in startCalipso.all_arrays.items(): 
        if valueca is not None:
            if valueca.size != 1:
                cal.all_arrays[arnameca] = valueca[start_break:end_break,...]
            else:
                cal.all_arrays[arnameca] = valueca
    return cal

def adjust5kmTo1kmresolution(calipso5km):
    logger.debug("Repeat 5km calipso data to fit 1km resoluiton")
    calipso= CalipsoObject()
    for arname, value in calipso5km.all_arrays.items(): 
        if value is not None:
            new_values = np.repeat(value, 5, axis=0)
            calipso.all_arrays[arname] = new_values                    
    return calipso 

def add5kmVariablesTo1kmresolution(calipso1km, calipso5km):
    logger.debug("Repeat 5km calipso data to fit 1km resoluiton")
    for variable in ["cfc_single_shots_1km_from_5km_file"]:
        if hasattr(calipso5km, variable):
            cfc_single_shot_1km_from_5km_file = getattr(calipso5km, variable)
            setattr(calipso1km, variable, cfc_single_shot_1km_from_5km_file)
    for variable_5km  in [ 
            "column_optical_depth_tropospheric_aerosols_1064",
            "column_optical_depth_tropospheric_aerosols_532",
            "column_optical_depth_tropospheric_aerosols_uncertainty_1064",
            "column_optical_depth_tropospheric_aerosols_uncertainty_532",
            "column_optical_depth_cloud_532",
            "column_optical_depth_cloud_uncertainty_532",
            #"feature_optical_depth_532",
            #"feature_optical_depth_uncertainty_532"
    ]:                
        if not hasattr(calipso5km, variable_5km) and CALIPSO_version3:
            #in version 3 names were different
            variable_5km = variable_5km.replace('tropospheric_','') 
        data = getattr(calipso5km, variable_5km)
        new_data = np.repeat(data, 5, axis=0)    
        setattr(calipso1km, variable_5km +"_5km", new_data)      
    for variable_5km  in ["feature_optical_depth_532_top_layer_5km",
                          "total_optical_depth_5km"]:                
        data = getattr(calipso5km, variable_5km)
        new_data = np.repeat(data, 5, axis=0)    
        setattr(calipso1km, variable_5km, new_data) 
    #for variable_5km  in ["layer_top_altitude", "layer_base_altitude"]:     
    #    data = getattr(calipso5km, variable_5km)#[:,0]*1000
    #    data[data<0] =-9
    #    new_data = np.repeat(data, 5, axis=0)    
    #    setattr(calipso1km, variable_5km + "_5km", new_data)  
    #cfc_5km = np.repeat(calipso5km.cloud_fraction, 5, axis=0)     
    #number_layers_found_5km = np.repeat(calipso5km.number_layers_found,
    # 5, axis=0).ravel() 
    isCloudOnlyIn5km = np.logical_and(calipso1km.cloud_fraction<0.1,
                                      calipso1km.total_optical_depth_5km>0)
    #Give clouds in only 5km  cloud fraction equal to 0.15
    #these clouds are either: 
    #1) thin or 
    #2) other pixels of the 5 1km pixles are cloudy
    #In case 1 we might want to have them as clouds in case 2 we do not.
    #For now just set all of them to 0.15.
    calipso1km.cloud_fraction[
        isCloudOnlyIn5km] = 0.15 
    if CALIPSO_version4:
        isCloudOnlyIn300m = np.logical_and(calipso1km.cloud_fraction<=0,
                                           cfc_single_shot_1km_from_5km_file>0)    
        calipso1km.cloud_fraction[
            isCloudOnlyIn300m] = cfc_single_shot_1km_from_5km_file[isCloudOnlyIn300m]
    return calipso1km 

def add1kmTo5km(Obj1, Obj5):
    # First check if length of 5 km and 1 km arrays correspond 
    # (i.e. 1 km array = 5 times longer array)
    # Here we check the middle time (index 1) out of the three time values 
    # given (start, mid, end) for 5 km data
    if (Obj5.profile_utc_time[:,1] == Obj1.profile_utc_time[2::5]).sum() != Obj5.profile_utc_time.shape[0]:
                              
        print("length mismatch")
        pdb.set_trace()

    #First making a preliminary check of the differences in fraction of cloudy
    # calipso columns in 1 km and 5 km data.
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
    # Now calculate the cloud fraction in 5 km data from 1 km data 
    # (discretized to 0.0, 0.2, 0.4, 0.6, 0.8 and 1.0).
    # In addition, if there are cloud layers in 5 km data but nothing in 1 km 
    # data, set cloud fraction to 1.0.
    # This latter case represents when very thin cloud layers are being 
    # detected over longer distances
    # Finally, if there are cloud layers in 1 km data but not in 5 km data, 
    # add a layer to 5 km data and set corresponding
    # COT to 1.0. Cloud base and cloud tp for this layer is calculated as 
    # averages from original levels (max height for
    # top and min height for base if there are more than one layer).
    # This is a pragmatic solution to take care of a
    # weakness or bug in the CALIPSO retrieval of clouds below 4 km
    Obj5.number_layers_found_1km = 1.0*Obj5.number_layers_found.copy()
    for i in range(Obj5.profile_utc_time.shape[0]):
        cfc = 0.0
        for j in range(5):
            if Obj1.number_layers_found[i*5+j] > 0:
                cfc = cfc + 0.2000
        Obj5.number_layers_found_1km[i] = cfc
        if cfc == 1.0:
            cfc = 0.99 # Just to be able to track the case 
                       # when no cloud layers existed in 5 km data
                       # but all 1 km FOVs were cloudy /KG 20170213  
        if (Obj5.number_layers_found[i] > 0):
            Obj5.cloud_fraction[i] = 1.0
        if ((cfc > 0.1) and (Obj5.number_layers_found[i] == 0)): 
            #Add missing layer due to CALIPSO processing bug
            cloudtop_sum = 0.0
            cloudbase_sum = 0.0
            cloud_layers = 0
            feature_array_list = []
            for j in range(5):
                if Obj1.number_layers_found[i*5+j] != 0:
                    for k in range(Obj1.number_layers_found[i*5+j]):
                        cloudtop_sum = (cloudtop_sum + 
                                        Obj1.layer_top_altitude[i,k])
                        cloudbase_sum = (cloudbase_sum + 
                                         Obj1.layer_base_altitude[i,k])
                        cloud_layers = cloud_layers + 1
                        feature_array_list.append(
                            Obj1.feature_classification_flags[i, k])
            Obj5.number_layers_found[i] = 1
            Obj5.layer_top_altitude[i, 0] = cloudtop_sum/cloud_layers
            Obj5.layer_base_altitude[i, 0] = cloudbase_sum/cloud_layers
            Obj5.feature_optical_depth_532[i, 0] = 1.0 #Just put it safely 
            # away from the thinnest cloud layers - the best we can do!
            # Obj5.feature_classification_flags[i, 0] = 22218 
            # if assuming like below:
            # cloud, low quality, water phase, low quality, low broken cumulus
            # , confident, 1 km horizontal averaging)
            feature_array = np.asarray(feature_array_list)
            Obj5.feature_classification_flags[i, 0] = np.median(
                feature_array[:]) # However, let's take the median value
            Obj5.single_shot_cloud_cleared_fraction[i] = 0.0 #Not used later
            Obj5.cloud_fraction[i] = cfc
    return Obj5


def addSingleShotTo5km(Obj5): # Valid only for CALIPSO-CALIOP version 4.10
    logger.info("Making use of new Single Shot cloud cleared information")
    cfc_5km = 0
    cfc_Single_Shot = 0
    len_5km = Obj5.profile_utc_time.shape[0]
    len_single_shot = Obj5.Number_cloudy_single_shots.shape[0]
    cloudy_5km = np.greater(Obj5.number_layers_found,0)
    cloudy_singleshot = np.greater_equal(Obj5.Number_cloudy_single_shots*1.0, 
                                         15.0*CALIPSO_CLOUDY_MIN_CFC)
    if len_5km == len_single_shot:
        cfc_5km = np.sum(cloudy_5km)
        cfc_Single_Shot = np.sum(cloudy_singleshot)
        cfc_combined = np.sum(np.logical_or(cloudy_5km.ravel(), 
                                            cloudy_singleshot.ravel())) 
    else:
        raise InputError("Different array length: 5km and Single Shot data!")
    logger.info(
        "*CHECKING CLOUD FREQUENCY DIFFERENCES SINGLE SHOT AND 5KM DATASETS:\n"
        " \t \t Number of 5 km FOVS:                 %d\n"
        " \t \t Number of cloudy 5 km FOVS:          %d\n"
        " \t \t Cloudy fraction 5 km:             %3.2f\n"
        " \t \t Number of Single shot FOVS:          %d\n"
        " \t \t Number of cloudy Single shot FOVS:   %d\n"
        " \t \t Cloudy fraction Single Shots:     %3.2f\n"
        " \t \t Number of Combined FOVS:             %d\n"
        " \t \t Number of cloudy combined FOVS:      %d\n"
        " \t \t Cloudy fraction combined:         %3.2f\n",
            len_5km,
            cfc_5km,
            float(cfc_5km)/float(len_5km),
            len_single_shot,
            cfc_Single_Shot,
            float(cfc_Single_Shot)/float(len_single_shot),
            len_single_shot,
            cfc_combined,
            float(cfc_combined)/float(len_single_shot))
    # Now calculate the cloud fraction in 5 km data from single shot data 
    # (fraction of 15 FOVs declared cloudy).
    # In addition, if there are cloud layers in 5 km data but nothing in 
    # single shot data, set cloud fraction to 1.0.
    # This latter case represents when very thin cloud layers are being 
    # detected over longer distances.
    # Finally, if there are cloud layers in single shot data but not in 5 km 
    # data, add a layer to 5 km data and set corresponding
    # COT to 1.0. Cloud base and cloud tp for this layer is calculated as 
    # averages from original levels (max height for
    # top and min height for base if there are more than one layer).
    # This is a pragmatic solution to take care of a
    # weakness or bug in the (CALIPSO-Version3?) retrieval of clouds below 4 km
    cfc = Obj5.Number_cloudy_single_shots.ravel()/15.0
    layers_found_5km = Obj5.number_layers_found.ravel()
    cfc[np.equal(cfc, 1.0)] = 0.99
    cfc[np.greater(layers_found_5km, 0)] = 1.0
    update_these = np.logical_and(np.greater(cfc, 0.001),
                                  np.equal(layers_found_5km, 0))
    Obj5.number_layers_found[update_these] = 1

    Obj5.layer_top_altitude[update_these, 0] = (
        Obj5.Average_cloud_top_single_shots[update_these].ravel()) 
    # Averaged from single shot data
    Obj5.layer_base_altitude[update_these, 0] = (
        Obj5.Average_cloud_base_single_shots[update_these].ravel()) 
    # Averaged from single shot data
    Obj5.feature_optical_depth_532[update_these, 0] = 1.0 #not very thin 
    Obj5.cloud_fraction = cfc
    return Obj5


if __name__ == "__main__":
    # Testing:
    pass
 
