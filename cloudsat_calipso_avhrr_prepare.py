# Program cloudsat_calipso_avhrr_prepare.py
# python cloudsat_calipso_avhrr_prepare.py 
import numpy as np
from config import OPTICAL_DETECTION_LIMIT
#Nina 20160314, Removed obsolete function CalipsoCloudFraction.
#If reintroduced remember that the calipso layer variables where before
#tranposed! So there might be needed change of indexing 
# i.e top layer before where in [0,:] and  top layer is now [:,0]

def CalipsoCloudOpticalDepth(cloud_top, cloud_base, optical_depth, cloud_fraction, fcf, min_optical_depth):
    from config import MIN_OPTICAL_DEPTH

    new_cloud_top = np.ones(cloud_top.shape,'d')*np.min(cloud_top)
    new_cloud_base = np.ones(cloud_base.shape,'d')*np.min(cloud_base)
    new_cloud_fraction = np.zeros(cloud_fraction.shape,'d')
    new_fcf = np.ones(fcf.shape).astype(fcf.dtype)
    
    for pixel_i in range(optical_depth.shape[0]):

        depthsum = 0 #Used to sum the optical_depth
        for layer_j in range(optical_depth.shape[1]):
            # Just stops the for loop when there are no more valid value 
            if optical_depth[pixel_i, layer_j] < 0:
                break
            else:
                depthsum = depthsum + optical_depth[pixel_i, layer_j]
            
                # Removes the cloud values for all pixels that have a optical depth (integrated from the top) below MIN_OPTICAL_DEPTH and moves the first valid value to the first column and so on.
                #if depthsum >= MIN_OPTICAL_DEPTh:
                if depthsum >= min_optical_depth:
                    new_cloud_top[pixel_i, 0:(optical_depth.shape[1]-layer_j)] = cloud_top[pixel_i, layer_j:]
                    new_cloud_base[pixel_i, 0:(optical_depth.shape[1]-layer_j)] = cloud_base[pixel_i, layer_j:]
                    # new_cloud_fraction[pixel_i] = 1
                    new_cloud_fraction[pixel_i] = cloud_fraction[pixel_i] # Let's still trust in what is seen in 1 km data/KG
                    new_fcf[pixel_i, 0:(fcf.shape[1]-layer_j)] = fcf[pixel_i, layer_j:]
                    #extrra to get a cloud top that corresponds better to avhrr
                    for k in range(new_cloud_top.shape[1]):
                        if new_cloud_top[pixel_i, k] < 0:
                            break
##                        elif new_cloud_top[pixel_i,k] - new_cloud_base[pixel_i,k] > 2:                  #Let's skip this option
##                            new_cloud_top[pixel_i,k] = new_cloud_base[pixel_i,k] + \                    #Not well motivated!/KG
##                            ((new_cloud_top[pixel_i,k] - new_cloud_base[pixel_i,k]) * 3/4.)
                        else:
                            new_cloud_top[pixel_i, k] = new_cloud_base[pixel_i, k] + \
                            ((new_cloud_top[pixel_i, k] - new_cloud_base[pixel_i, k]) * 1/2.)
                    break
    new_validation_height = new_cloud_top[:,0].copy()
    new_validation_height[new_validation_height>=0] *= 1000
    new_validation_height[new_validation_height<0] = -9
    return new_cloud_top, new_cloud_base, new_cloud_fraction, new_fcf, new_validation_height

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
        CaObj.calipso.layer_base_altitude[:,0] > CaObj.calipso.detection_height_5km,
        CaObj.calipso.layer_base_altitude[:,0]+0.1,
        CaObj.calipso.detection_height_5km)
    clouds_to_update = np.logical_and(
        CaObj.calipso.layer_top_altitude[:,0]>CaObj.calipso.detection_height_5km,
        np.not_equal(CaObj.calipso.detection_height_5km, -9))
    CaObj.calipso.validation_height = np.where(
        clouds_to_update,
        new_cloud_tops,
        CaObj.calipso.validation_height)
    return CaObj

def CalipsoOpticalDepthSetThinToClearFiltering1km(CaObj):
    isThin_clouds = np.logical_and(CaObj.calipso.total_optical_depth_5km < OPTICAL_DETECTION_LIMIT,
                                   np.not_equal(CaObj.calipso.total_optical_depth_5km,-9))
    #>0.0 important. Some clouds are missing in 5km data set but present in 1km data set!
    set_to_clear = np.logical_and(
        CaObj.calipso.number_layers_found>0,
        isThin_clouds)
    CaObj.calipso.cloud_fraction[set_to_clear] = 0.00001
    CaObj.calipso.validation_height[set_to_clear] = -9

    return CaObj
 






















































