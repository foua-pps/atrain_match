import pdb
import os
import numpy as np
import time as tm
import logging
logger = logging.getLogger(__name__)


import config
from utils.common import (MatchupError, TimeMatchError, 
                    InputError, ProcessingError, 
                    elements_within_range)

from matchobject_io import (CalipsoImagerTrackObject,
                            CalipsoObject)
from utils.runutils import do_some_logging

def add_validation_ctth_calipso(calipso):
    calipso.validation_height = calipso.layer_top_altitude[:,0].copy()
    calipso.validation_height[calipso.validation_height>=0] = (
        1000.0*calipso.validation_height[calipso.validation_height>=0])
    calipso.validation_height[calipso.validation_height<0] = -9
    return calipso

def calipso_track_from_matched(retv_calipso, calipso, idx_match):
    # Calipso line,pixel inside IMAGER swath:
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


def match_calipso_imager(values, 
                        caObj, caObjAerosol, 
                        imagerGeoObj, imagerObj, 
                        ctype, cma, ctth, cpp, nwp_obj,
                        imagerAngObj, nwp_segments, SETTINGS, res=config.RESOLUTION):

    from utils.common import map_imager
    retv = CalipsoImagerTrackObject()
    retv.imager_instrument = imagerObj.instrument.lower()

    cal, cap = map_imager(imagerGeoObj, 
                         caObj.longitude.ravel(),
                         caObj.latitude.ravel(),
                         radius_of_influence=config.RESOLUTION*0.7*1000.0) # somewhat larger than radius...
    #warn if no matches
    calnan = np.where(cal == config.NODATA, np.nan, cal)
    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal,cap)]
        imager_lines_sec_1970 = np.where(cal != config.NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal != config.NODATA, imagerGeoObj.time[cal], np.nan)
    idx_match = elements_within_range(caObj.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr"]) 
    if idx_match.sum() == 0:
        logger.warning("No matches in region within time threshold %d s.", SETTINGS["sec_timeThr"])
        return None    
    retv.calipso = calipso_track_from_matched(retv.calipso, caObj, idx_match)
    
    # Calipso line,pixel inside IMAGER swath:
    retv.calipso.imager_linnum = np.repeat(cal, idx_match).astype('i')
    retv.calipso.imager_pixnum = np.repeat(cap, idx_match).astype('i')

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.calipso.sec_1970 - retv.imager.sec_1970
    do_some_logging(retv, caObj)
    logger.debug("Generate the latitude,cloudtype tracks!")
    from libs.extract_imager_along_track import imager_track_from_matched
    retv = imager_track_from_matched(retv, SETTINGS, 
                                     imagerGeoObj, imagerObj, imagerAngObj, 
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
        #Se low cloudfraction on clear 1km pixels that might be cloudy.
        #don't use filter to set cloudy pixels to clear
        #If winsz=3: 1clear 2cloudy => cfc = 0.066
        #   winsz=3; 2clear 1cloudy => cfc = 0.033
        cal.cloud_fraction = np.where(
            np.logical_and(cal.cloud_fraction<1.0,
                           cloud_fraction_temp>0.01),            
            0.01*cloud_fraction_temp,cal.cloud_fraction)
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
    if config.RESOLUTION == 1:
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

def discardCalipsoFilesOutsideTimeRange(calipsofiles_list, imagerGeoObj, values, SETTINGS, res=config.RESOLUTION, ALAY=False):
    imager_end = imagerGeoObj.sec1970_end
    imager_start = imagerGeoObj.sec1970_start
    calipso_within_time_range = []
    for current_file in calipsofiles_list:
        newCalipso = get_calipso(current_file, res, ALAY=ALAY)
        cal_new_all = newCalipso.sec_1970
        if cal_new_all[0]>imager_end + SETTINGS["sec_timeThr"] or  cal_new_all[-1] + SETTINGS["sec_timeThr"]<imager_start:
            pass
            #print "skipping file %s outside time_limits"%(current_file)
        else:
            logger.debug("Keeping file %s inside time_limits", os.path.basename(current_file))
            calipso_within_time_range.append(current_file)
    return calipso_within_time_range

def reshapeCalipso(calipsofiles, res=config.RESOLUTION, ALAY=False):
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

def find_break_points(startCalipso, imagerGeoObj, SETTINGS):
    """
    Find the start and end point where calipso and imager matches is within 
    time limits.
    """
    imager_end = imagerGeoObj.sec1970_end
    imager_start = imagerGeoObj.sec1970_start
    # Finds Break point
    start_break = np.argmin((np.abs((startCalipso.sec_1970) 
                                    - (imager_start - SETTINGS["sec_timeThr"]))))
    end_break = np.argmin((np.abs((startCalipso.sec_1970) 
                                  - (imager_end + SETTINGS["sec_timeThr"])))) + 2    # Plus two to get one extra, just to be certain    
    if start_break != 0:
        start_break = start_break - 1 # Minus one to get one extra, just to be certain
    return start_break, end_break

def time_reshape_calipso(startCalipso,
                         start_break, end_break):
    """
    Cut the calipso data at the point where matches with imager is within 
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

def add5kmVariablesTo1kmresolution(calipso1km, calipso5km, CALIPSO_version):
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
        if CALIPSO_version == 3:
            if not hasattr(calipso5km, variable_5km) or getattr(calipso5km, variable_5km) is None:
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
    if CALIPSO_version == 4:
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


def addSingleShotTo5km(Obj5, SETTINGS): # Valid only for CALIPSO-CALIOP version 4.10
    logger.info("Making use of new Single Shot cloud cleared information")
    cfc_5km = 0
    cfc_Single_Shot = 0
    len_5km = Obj5.profile_utc_time.shape[0]
    len_single_shot = Obj5.Number_cloudy_single_shots.shape[0]
    cloudy_5km = np.greater(Obj5.number_layers_found,0)
    cloudy_singleshot = np.greater_equal(Obj5.Number_cloudy_single_shots*1.0, 
                                         15.0*SETTINGS['CALIPSO_CLOUDY_MIN_CFC'])
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

def CalipsoCloudOpticalDepth_new(calipso, min_optical_depth, use_old_method=False,
                                 limit_ctop=0.2):
    new_cloud_top = np.ones(calipso.layer_top_altitude.shape)*config.NODATA*1.0
    new_cloud_base = np.ones(calipso.layer_base_altitude.shape)*config.NODATA*1.0
    new_cloud_top_pressure = np.ones(calipso.layer_top_altitude.shape)*config.NODATA*1.0
    new_cloud_base_pressure = np.ones(calipso.layer_base_altitude.shape)*config.NODATA*1.0
    distance_down_in_cloud_we_see = np.ones(calipso.layer_base_altitude.shape)*config.NODATA*1.0
    new_cloud_fraction = np.zeros(calipso.cloud_fraction.shape)
    new_fcf = np.ones(calipso.feature_classification_flags.shape).astype(calipso.feature_classification_flags.dtype)
    depthsum = 0.0*new_cloud_fraction.copy()
    already_detected = depthsum > 99999 #None from the start! All False
    already_detected[calipso.layer_top_altitude[:,0]<0] = True # don't bother with clear pixels
    #to_thin_to_see_at_all = calipso.total_optical_depth_5km < min_optical_depth
    N10 = new_cloud_top.shape[1]
    for layer_j in xrange(N10):
        #Filter OPTICAL_LIMIT_CLOUD_TOP down in EACH layer
        #Notice even if top of layer i is always above layer i-1,
        #it is NOT the case that base of layer i is always above layer i-1!
        # We could have the situation:
        # layer 0 base at 2km top at 0km optical thickness 0.1
        # layer 1 base at 8km top at 9km optical thickness 10
        # Or the situation:
        # layer 0 base at 2km top at 0km optical thickness 0.5
        # layer 1 base at 8km top at 9km optical thickness 0.5
        this_layer_optical_depth = calipso.feature_optical_depth_532[:, layer_j].ravel()
        cloud_base = calipso.layer_base_altitude[:,layer_j]
        cloud_top = calipso.layer_top_altitude[:,layer_j]
        geometrical_distance_cloud = cloud_top-cloud_base
        geometrical_distance_cloud[cloud_base<0] = 0
        geometrical_distance_cloud[cloud_top<0] = 0                                         
        fraction_into_cloud = (limit_ctop*1.0)/this_layer_optical_depth
        fraction_into_cloud[fraction_into_cloud<0] = 0
        fraction_into_cloud[fraction_into_cloud>1.0] = 1.0
        if use_old_method:
            fraction_into_cloud = 0.5
        #For KGs method set fraction_into_cloud = 0.5 always
        distance_down_in_cloud_we_see[:,layer_j] = geometrical_distance_cloud*fraction_into_cloud 
    for layer_j in xrange(N10):
        total_optical_depth_layers_above = depthsum.copy()
        this_layer_optical_depth = calipso.feature_optical_depth_532[:, layer_j].ravel()
        depthsum += this_layer_optical_depth 
        #this_is_the_top_most_layer_we_detect
        update = np.logical_and(
            depthsum >= min_optical_depth,
            np.equal(already_detected,False))
        new_cloud_top[update, 0:(N10-layer_j)] = (calipso.layer_top_altitude[update, layer_j:] -
                                                  distance_down_in_cloud_we_see[update, layer_j:])
        new_cloud_base[update, 0:(N10-layer_j)] = calipso.layer_base_altitude[update, layer_j:]
        new_cloud_top_pressure[update, 0:(N10-layer_j)] = (calipso.layer_top_pressure[update, layer_j:])
        new_cloud_base_pressure[update, 0:(N10-layer_j)] = calipso.layer_base_pressure[update, layer_j:]
        new_fcf[update, 0:(N10-layer_j)] = calipso.feature_classification_flags[update, layer_j:]
        new_cloud_fraction[update] = calipso.cloud_fraction[update]
        already_detected[update] = True
    if use_old_method:
        pass
    else:    
        for layer_j in xrange(N10):
            filtered_to_low = new_cloud_top[:,0]< new_cloud_top[:,layer_j]
            new_cloud_top[filtered_to_low,0] = new_cloud_top[filtered_to_low,layer_j]
        
    new_validation_height = new_cloud_top[:,0].copy()
    new_validation_height[new_validation_height>=0] = new_validation_height[new_validation_height>=0]*1000
    new_validation_height[new_validation_height<0] = -9
    return (new_cloud_top, new_cloud_base, new_cloud_fraction, new_fcf, new_validation_height, 
            new_cloud_top_pressure, new_cloud_base_pressure)


def check_total_optical_depth_and_warn(caObj):
    obj = caObj.calipso
    if  (obj.total_optical_depth_5km is not None and 
         (obj.total_optical_depth_5km < obj.feature_optical_depth_532_top_layer_5km).any()):
        badPix=np.less(obj.total_optical_depth_5km+0.001, 
                       obj.feature_optical_depth_532_top_layer_5km)
        diff=obj.total_optical_depth_5km- obj.feature_optical_depth_532_top_layer_5km
        print "warning", len(obj.total_optical_depth_5km)  
        print len(obj.total_optical_depth_5km[badPix])
        print obj.total_optical_depth_5km[badPix]
        print obj.feature_optical_depth_532_top_layer_5km[badPix] 
        print diff[badPix]
        print obj.number_layers_found[badPix]
        if obj.detection_height_5km is not None:
            print obj.detection_height_5km[badPix] 
        print np.where(badPix)
        print obj.layer_top_altitude[0,badPix]
        print obj.layer_base_altitude[0,badPix]

def CalipsoOpticalDepthHeightFiltering1km(CaObj):
    #Should all layers be updated???
    new_cloud_tops = np.where(
        CaObj.calipso.layer_top_altitude[:,0]*1000 > CaObj.calipso.detection_height_5km,
        CaObj.calipso.detection_height_5km,
        CaObj.calipso.layer_top_altitude[:,0]*1000)
    new_cloud_tops = np.where(
        CaObj.calipso.layer_base_altitude[:,0]*1000 > CaObj.calipso.detection_height_5km,
        CaObj.calipso.layer_base_altitude[:,0]*1000,
        new_cloud_tops)
    clouds_to_update = np.logical_and(
        CaObj.calipso.layer_top_altitude[:,0]*1000>CaObj.calipso.detection_height_5km,
        np.not_equal(CaObj.calipso.detection_height_5km, -9))
    return np.where(clouds_to_update,
                    new_cloud_tops,
                    CaObj.calipso.validation_height)


def detection_height_from_5km_data(Obj1, Obj5, limit_ctop=0.2):
    if (Obj5.profile_utc_time[:,1] == Obj1.profile_utc_time[2::5]).sum() != Obj5.profile_utc_time.shape[0]:
        logger.warning("length mismatch")
    retv = CalipsoCloudOpticalDepth_new(Obj5, 0.0, 
                                        limit_ctop=limit_ctop)
    cloud_top5km, dummy, dummy, dummy, detection_height, dummy, dummy = retv
    #cloud_base1km = Obj1.layer_base_altitude[:,0]
    #cloud_top5km = np.repeat(cloud_top5km[:,0], 5, axis=0) 
    detection_height1km = np.repeat(detection_height, 5, axis=0) 
    #detection_height1km[cloud_top5km<cloud_base1km] = -9
    Obj1.detection_height_5km = detection_height1km
    return Obj1    

def CalipsoOpticalDepthSetThinToClearFiltering1km(CaObj, SETTINGS):
    isThin_clouds = np.logical_and(
        CaObj.calipso.total_optical_depth_5km < SETTINGS['OPTICAL_DETECTION_LIMIT'],
        np.not_equal(CaObj.calipso.total_optical_depth_5km,-9))
    #>0.0 important. Some clouds are missing in 5km data set but present in 1km data set!
    set_to_clear = np.logical_and(
        CaObj.calipso.number_layers_found>0,
        isThin_clouds)
    cloud_fraction =  CaObj.calipso.cloud_fraction
    validation_height = CaObj.calipso.validation_height
    cloud_fraction[set_to_clear] = 0.00001
    validation_height[set_to_clear] = -9
    return cloud_fraction, validation_height



if __name__ == "__main__":
    # Testing:
    pass
 
