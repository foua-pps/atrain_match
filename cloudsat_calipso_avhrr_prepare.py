# Program cloudsat_calipso_avhrr_prepare.py
# python cloudsat_calipso_avhrr_prepare.py 
import numpy as np
from config import OPTICAL_DETECTION_LIMIT
def CalipsoCloudFraction(calipsoObj):
     # This function is obsolete! Now we merge 1 km and 5 km data instead!!
    # Then we take cloud fraction from mainly 1 km, only from 5 km (=1.0) for thinnest clouds!/KG
    new_cloud_fraction = np.zeros(calipsoObj.calipso.cloud_fraction.shape, 'd')

    fcf = calipsoObj.calipso.feature_classification_flags
    new_fcf = np.ones(fcf.shape).astype(fcf.dtype)

    cloud_top = calipsoObj.calipso.cloud_top_profile
    new_cloud_top = np.ones(cloud_top.shape, 'd')*np.min(cloud_top)

    cloud_base = calipsoObj.calipso.cloud_base_profile
    new_cloud_base = np.ones(cloud_base.shape, 'd')*np.min(cloud_base)

    optical_depth = calipsoObj.calipso.optical_depth
    new_optical_depth = np.ones(optical_depth.shape, 'd')*np.min(optical_depth)

    ssccf = calipsoObj.calipso.single_shot_cloud_cleared_fraction
    new_ssccf = np.ones(ssccf.shape,'d')*np.min(ssccf)
    for i in range(ssccf.shape[1]):
        ind = ((ssccf[:,i]>=0) & (ssccf[:,i]<=0.5))
        
        new_cloud_top[:,i] = np.where(ind, cloud_top[:,i], new_cloud_top[:,i])
        new_cloud_base[:,i] = np.where(ind, cloud_base[:,i], new_cloud_base[:,i])
        new_fcf[:,i] = np.where(ind, fcf[:,i], new_fcf[:,i])
        new_ssccf[:,i] = np.where(ind, ssccf[:,i], new_ssccf[:,i])
        if ind.any():
            new_cloud_fraction[i] = cloud_fraction[i]
            new_fcf[:,i] = fcf[:,i]
            new_cloud_top[:,i] = cloud_top[:,i]
            new_cloud_base[:,i] = cloud_base[:,i]
            new_optical_depth[:,i] = optical_depth[:,i]
            new_ssccf[:,i] = ssccf[:,i]
        
    return new_cloud_top, new_cloud_base, new_optical_depth, new_cloud_fraction, new_fcf, new_ssccf


def NinaTestarMedelCloudBaseAndTop(cloud_top, cloud_base):
    #Maybe we should use middle of cloud as cloud top. To get better than nonsens comparison for high (optical thin) clouds. Take middel only for clouds the are certainly very thin, deeper than 1km!
    new_cloud_top = np.where((cloud_top - cloud_base)>1, 0.5*(cloud_top + cloud_base), cloud_top)
    return new_cloud_top


def CloudsatCloudOpticalDepth(cloud_top, cloud_base, optical_depth, cloud_fraction, fcf, min_optical_depth):
    from config import MIN_OPTICAL_DEPTH

    new_cloud_top = np.ones(cloud_top.shape,'d')*np.min(cloud_top)
    new_cloud_base = np.ones(cloud_base.shape,'d')*np.min(cloud_base)
    new_cloud_fraction = np.zeros(cloud_fraction.shape,'d')
    new_fcf = np.ones(fcf.shape).astype(fcf.dtype)
    
    for i in range(optical_depth.shape[1]):

        depthsum = 0 #Used to sum the optical_depth
        for j in range(optical_depth.shape[0]):
            # Just stops the for loop when there are no more valid value 
            if optical_depth[j,i] < 0:
                break
            else:
                depthsum = depthsum + optical_depth[j,i]
            
                # Removes the cloud values for all pixels that have a optical depth (integrated from the top) below MIN_OPTICAL_DEPTH and moves the first valid value to the first column and so on.
                #if depthsum >= MIN_OPTICAL_DEPTh:
                if depthsum >= min_optical_depth:
                    new_cloud_top[0:(optical_depth.shape[0]-j),i] = cloud_top[j:,i]
                    new_cloud_base[0:(optical_depth.shape[0]-j),i] = cloud_base[j:,i]
                    # new_cloud_fraction[i] = 1
                    new_cloud_fraction[i] = cloud_fraction[i] # Let's still trust in what is seen in 1 km data/KG
                    new_fcf[0:(fcf.shape[0]-j),i] = fcf[j:,i]
                    #extrra to get a cloud top that corresponds better to avhrr
                    for k in range(new_cloud_top.shape[0]):
                        if new_cloud_top[k, i] < 0:
                            break
##                        elif new_cloud_top[k, i] - new_cloud_base[k, i] > 2:                  #Let's skip this option
##                            new_cloud_top[k, i] = new_cloud_base[k, i] + \                    #Not well motivated!/KG
##                            ((new_cloud_top[k, i] - new_cloud_base[k, i]) * 3/4.)
                        else:
                            new_cloud_top[k, i] = new_cloud_base[k, i] + \
                            ((new_cloud_top[k, i] - new_cloud_base[k, i]) * 1/2.)
                    break

    return new_cloud_top, new_cloud_base, new_cloud_fraction, new_fcf

def check_total_optical_depth_and_warn(caObj):
    obj = caObj.calipso
    if  (obj.total_optical_depth_5km is not None and 
         (obj.total_optical_depth_5km < obj.optical_depth_top_layer5km).any()):
        badPix=np.less(obj.total_optical_depth_5km+0.001, 
                       obj.optical_depth_top_layer5km)
        diff=obj.total_optical_depth_5km- obj.optical_depth_top_layer5km
        print "warning", len(obj.total_optical_depth_5km)  
        print len(obj.total_optical_depth_5km[badPix])
        print obj.total_optical_depth_5km[badPix]
        print obj.optical_depth_top_layer5km[badPix] 
        print diff[badPix]
        print obj.number_of_layers_found[badPix]
        if obj.detection_height_5km is not None:
            print obj.detection_height_5km[badPix] 
        print np.where(badPix)
        print obj.cloud_top_profile[0,badPix]
        print obj.cloud_base_profile[0,badPix]

def CalipsoOpticalDepthHeightFiltering1km(CaObj):
    #print CaObj.calipso.cloud_top_profile[0,:]
    new_cloud_tops = np.where(
        CaObj.calipso.cloud_base_profile[0,:] > CaObj.calipso.detection_height_5km,
        CaObj.calipso.cloud_base_profile[0,:]+0.1,
        CaObj.calipso.detection_height_5km)
    clouds_to_update = np.logical_and(
        CaObj.calipso.cloud_top_profile[0,:]>CaObj.calipso.detection_height_5km,
        np.not_equal(CaObj.calipso.detection_height_5km, -9))
    CaObj.calipso.cloud_top_profile[0,:] = np.where(
        clouds_to_update,
        new_cloud_tops,
        CaObj.calipso.cloud_top_profile[0,:])
    #print CaObj.calipso.cloud_top_profile[1,:]
    #print CaObj.calipso.cloud_top_profile[0,:]
    #print CaObj.calipso.detection_height_5km
    return CaObj

def CalipsoOpticalDepthSetThinToClearFiltering1km(CaObj):
    isThin_clouds = np.logical_and(CaObj.calipso.total_optical_depth_5km < OPTICAL_DETECTION_LIMIT,
                                   np.not_equal(CaObj.calipso.total_optical_depth_5km,-9))
    #>0.0 important. Some clouds are missing in 5km data set but present in 1km data set!
    set_to_clear = np.logical_and(
        CaObj.calipso.number_of_layers_found>0,
        isThin_clouds)
    #print set_to_clear.shape
    for lay in xrange(CaObj.calipso.cloud_base_profile.shape[0]):
        #print lay
        CaObj.calipso.cloud_top_profile[lay,:] = np.where(
            set_to_clear,
            -9999,
            CaObj.calipso.cloud_top_profile[lay,:])
        CaObj.calipso.cloud_base_profile[lay,:] = np.where(
            set_to_clear,
            -9999,
            CaObj.calipso.cloud_base_profile[lay,:])
    CaObj.calipso.cloud_fraction = np.where(
        set_to_clear,
        0.5,
        CaObj.calipso.cloud_fraction)

    #print CaObj.calipso.cloud_base_profile.shape      
    return CaObj
 






















































