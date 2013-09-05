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

def CloudsatOpticalDepthHeightFiltering1km(CaObj):
    #print CaObj.calipso.cloud_top_profile[0,:]
    new_cloud_tops = np.where(
        CaObj.calipso.cloud_base_profile[0,:]>CaObj.calipso.detection_height_5km,
        CaObj.calipso.cloud_base_profile[0,:]+0.1,
        CaObj.calipso.detection_height_5km)
    CaObj.calipso.cloud_top_profile[0,:] = np.where(
        CaObj.calipso.cloud_top_profile[0,:]>CaObj.calipso.detection_height_5km,
        new_cloud_tops,
        CaObj.calipso.cloud_top_profile[0,:])
    #print CaObj.calipso.cloud_top_profile[1,:]
    #print CaObj.calipso.cloud_top_profile[0,:]
    #print CaObj.calipso.detection_height_5km
    return CaObj

def CalipsoOpticalDepthSeToClearFiltering1km(CaObj):
    thin_clouds = np.where(CaObj.calipso.total_optical_depth_5km < OPTICAL_DETECTION_LIMIT,1,0) 
    set_to_clear = np.where(
        np.logical_and(
            CaObj.calipso.number_of_layers_found>0,
            CaObj.calipso.total_optical_depth_5km < OPTICAL_DETECTION_LIMIT),
        1,0)
    #print set_to_clear.shape
    for lay in xrange(CaObj.calipso.cloud_base_profile.shape[0]):
        #print lay
        CaObj.calipso.cloud_top_profile[lay,:] = np.where(
            set_to_clear==1,
            -9999,
            CaObj.calipso.cloud_top_profile[lay,:])
        CaObj.calipso.cloud_base_profile[lay:] = np.where(
            set_to_clear==1,
            -9999,
            CaObj.calipso.cloud_base_profile[lay:])
        CaObj.calipso.cloud_fraction = np.where(
            set_to_clear,
            0.1,
            CaObj.calipso.cloud_fraction)

    #print CaObj.calipso.cloud_base_profile.shape      
    return CaObj
 
#---------------------------------------------------------------------
def CloudsatAvhrrSatz(clsatObj):
    import np
    from config import AZIMUTH_RANGE
    
    # CloudSat:
    clsat_satz_ok = np.where((AZIMUTH_RANGE[0] <= clsatObj.avhrr.satz) * \
                                (clsatObj.avhrr.satz <= AZIMUTH_RANGE[1]))
    
    clsat_satz_ok_bools = (AZIMUTH_RANGE[0] <= clsatObj.avhrr.satz) == \
                                (clsatObj.avhrr.satz <= AZIMUTH_RANGE[1])
    
    n_clsat_satz_ok = len(clsat_satz_ok[0])
    n_clsat_satz_ok_bools = len(clsat_satz_ok_bools)
    clsatObj.cloudsat.longitude = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.longitude,np.nan)
    clsatObj.cloudsat.latitude = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.latitude,np.nan)
    clsatObj.cloudsat.avhrr_linnum = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.avhrr_linnum,np.nan)
    clsatObj.cloudsat.avhrr_pixnum = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.avhrr_pixnum,np.nan)

    for i in range(clsatObj.cloudsat.cloud_mask.shape[0]):
        clsatObj.cloudsat.cloud_mask[i,:] = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.cloud_mask[i,:],np.nan)
        clsatObj.cloudsat.Radar_Reflectivity[i,:] = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.Radar_Reflectivity[i,:],np.nan)
        clsatObj.cloudsat.Height[i,:] = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.Height[i,:],np.nan)
    clsatObj.cloudsat.echo_top = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.echo_top,np.nan)
    clsatObj.cloudsat.SurfaceHeightBin = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.SurfaceHeightBin,np.nan)
    clsatObj.cloudsat.SurfaceHeightBin_fraction = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.SurfaceHeightBin_fraction,np.nan)
    
    clsatObj.cloudsat.elevation = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.elevation,np.nan)
    clsatObj.cloudsat.sec_1970 = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.sec_1970,np.nan)
    clsatObj.cloudsat.MODIS_Cloud_Fraction = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.MODIS_Cloud_Fraction,np.nan)
    clsatObj.cloudsat.MODIS_cloud_flag = np.where(clsat_satz_ok_bools,clsatObj.cloudsat.MODIS_cloud_flag,np.nan)

    # AVHRR - Cloudsat:
    clsatObj.avhrr.longitude = np.where(clsat_satz_ok_bools,clsatObj.avhrr.longitude,np.nan)
    clsatObj.avhrr.latitude = np.where(clsat_satz_ok_bools,clsatObj.avhrr.latitude,np.nan)
    clsatObj.avhrr.sec_1970 = np.where(clsat_satz_ok_bools,clsatObj.avhrr.sec_1970,np.nan)
    clsatObj.avhrr.cloudtype = np.where(clsat_satz_ok_bools,clsatObj.avhrr.cloudtype,np.nan)
    clsatObj.avhrr.ctth_height = np.where(clsat_satz_ok_bools,clsatObj.avhrr.ctth_height,np.nan)
    clsatObj.avhrr.ctth_pressure = np.where(clsat_satz_ok_bools,clsatObj.avhrr.ctth_pressure,np.nan)
    clsatObj.avhrr.ctth_temperature = np.where(clsat_satz_ok_bools,clsatObj.avhrr.ctth_temperature,np.nan)
    clsatObj.avhrr.bt11micron = np.where(clsat_satz_ok_bools,clsatObj.avhrr.bt11micron,np.nan)
    clsatObj.avhrr.bt12micron = np.where(clsat_satz_ok_bools,clsatObj.avhrr.bt12micron,np.nan)
    clsatObj.avhrr.surftemp = np.where(clsat_satz_ok_bools,clsatObj.avhrr.surftemp,np.nan)
    clsatObj.avhrr.satz = np.where(clsat_satz_ok_bools,clsatObj.avhrr.satz,np.nan)
    return clsatObj

def CalipsoAvhrrSatz(caObj):
    from config import AZIMUTH_RANGE
    import pdb
    # Calipso:
    ca_satz_ok = np.where((AZIMUTH_RANGE[0] <= caObj.avhrr.satz) * \
                                (caObj.avhrr.satz <= AZIMUTH_RANGE[1]))
    ca_satz_ok_bools = (AZIMUTH_RANGE[0] <= caObj.avhrr.satz) == \
                                (caObj.avhrr.satz <= AZIMUTH_RANGE[1])                         
    n_ca_satz_ok = len(ca_satz_ok[0])
    n_ca_satz_ok_bools = len(ca_satz_ok_bools)
    caObj.calipso.longitude = np.where(ca_satz_ok_bools,caObj.calipso.longitude,np.nan)
    caObj.calipso.latitude = np.where(ca_satz_ok_bools,caObj.calipso.latitude,np.nan)
    caObj.calipso.avhrr_linnum = np.where(ca_satz_ok_bools,caObj.calipso.avhrr_linnum,np.nan)
    caObj.calipso.avhrr_pixnum = np.where(ca_satz_ok_bools,caObj.calipso.avhrr_pixnum,np.nan)
    
    caObj.calipso.cloud_fraction = np.where(ca_satz_ok_bools,caObj.calipso.cloud_fraction,np.nan)
    caObj.calipso.elevation = np.where(ca_satz_ok_bools,caObj.calipso.elevation,np.nan)
    caObj.calipso.number_of_layers_found = np.where(ca_satz_ok_bools,caObj.calipso.number_of_layers_found,np.nan)
    
    caObj.calipso.igbp = np.where(ca_satz_ok_bools,caObj.calipso.igbp,np.nan)
    caObj.calipso.nsidc = np.where(ca_satz_ok_bools,caObj.calipso.nsidc  ,np.nan)  
    caObj.calipso.sec_1970 = np.where(ca_satz_ok_bools,caObj.calipso.sec_1970,np.nan)  

    for j in range(caObj.calipso.cloud_top_profile.shape[0]):
        caObj.calipso.cloud_top_profile[j,:] = np.where(ca_satz_ok_bools, caObj.calipso.cloud_top_profile[j,:], np.nan)
        caObj.calipso.cloud_base_profile[j,:] = np.where(ca_satz_ok_bools,caObj.calipso.cloud_base_profile[j,:], np.nan)
        caObj.calipso.cloud_mid_temperature[j,:] = np.where(ca_satz_ok_bools, caObj.calipso.cloud_mid_temperature[j,:], np.nan)
        try:
            caObj.calipso.feature_classification_flags[j,:] = np.where(ca_satz_ok_bools,caObj.calipso.feature_classification_flags[j,:],np.nan)
        except:
            print "No feature_classification_flags array in file!"
            pass
    # AVHRR - Calipso:
    caObj.avhrr.longitude = np.where(ca_satz_ok_bools,caObj.avhrr.longitude,np.nan)
    caObj.avhrr.latitude = np.where(ca_satz_ok_bools,caObj.avhrr.latitude,np.nan)
    caObj.avhrr.sec_1970 = np.where(ca_satz_ok_bools,caObj.avhrr.sec_1970,np.nan)
    caObj.avhrr.cloudtype = np.where(ca_satz_ok_bools,caObj.avhrr.cloudtype,np.nan)
    caObj.avhrr.ctth_height = np.where(ca_satz_ok_bools,caObj.avhrr.ctth_height,np.nan)
    caObj.avhrr.ctth_pressure = np.where(ca_satz_ok_bools,caObj.avhrr.ctth_pressure,np.nan)
    caObj.avhrr.ctth_temperature = np.where(ca_satz_ok_bools,caObj.avhrr.ctth_temperature,np.nan)
    caObj.avhrr.bt11micron = np.where(ca_satz_ok_bools,caObj.avhrr.bt11micron,np.nan)
    caObj.avhrr.bt12micron = np.where(ca_satz_ok_bools,caObj.avhrr.bt12micron,np.nan)
    caObj.avhrr.satz = np.where(ca_satz_ok_bools,caObj.avhrr.satz,np.nan)

    if caObj.avhrr.surftemp == None:        
        print("No Surftemp. Continue")
    else:
        caObj.avhrr.surftemp = np.where(ca_satz_ok_bools,caObj.avhrr.surftemp,np.nan)

    return caObj





















































