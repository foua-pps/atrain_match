# Program cloudsat_calipso_avhrr_prepare.py
# python cloudsat_calipso_avhrr_prepare.py 


def CalipsoCloudFraction(cloud_top, cloud_base, optical_depth, cloud_fraction, fcf, ssccf):
    import numpy
    import pdb

    # This function is obsolete! Now we merge 1 km and 5 km data instead!!
    # Then we take cloud fraction from mainly 1 km, only from 5 km (=1.0) for thinnest clouds!/KG

    new_cloud_top = numpy.ones(cloud_top.shape,'d')*numpy.min(cloud_top)
    new_cloud_base = numpy.ones(cloud_base.shape,'d')*numpy.min(cloud_base)
    new_cloud_fraction = numpy.zeros(cloud_fraction.shape,'d')
    new_fcf = numpy.ones(fcf.shape).astype(fcf.dtype)
    new_ssccf = numpy.ones(ssccf.shape,'d')*numpy.min(ssccf)
    for i in range(ssccf.shape[1]):
        ind = ((ssccf[:,i]>=0) & (ssccf[:,i]<=0.5))
        
        new_cloud_top[:,i] = numpy.where(ind, cloud_top[:,i], new_cloud_top[:,i])
        new_cloud_base[:,i] = numpy.where(ind, cloud_base[:,i], new_cloud_base[:,i])
        new_fcf[:,i] = numpy.where(ind, fcf[:,i], new_fcf[:,i])
        new_ssccf[:,i] = numpy.where(ind, ssccf[:,i], new_ssccf[:,i])
        if ind.any():
            new_cloud_fraction[i] = cloud_fraction[i]
       
    return new_cloud_top, new_cloud_base, new_optical_depth, new_cloud_fraction, new_fcf, new_ssccf

def CloudsatCloudOpticalDepth(cloud_top, cloud_base, optical_depth, cloud_fraction, fcf, min_optical_depth):
    import numpy
    import pdb
    # from config import MIN_OPTICAL_DEPTH

    new_cloud_top = numpy.ones(cloud_top.shape,'d')*numpy.min(cloud_top)
    new_cloud_base = numpy.ones(cloud_base.shape,'d')*numpy.min(cloud_base)
    new_cloud_fraction = numpy.zeros(cloud_fraction.shape,'d')
    new_fcf = numpy.ones(fcf.shape).astype(fcf.dtype)
    
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
                
#---------------------------------------------------------------------
def CloudsatAvhrrSatz(clsatObj):
    import numpy
    from config import AZIMUTH_RANGE
    
    # CloudSat:
    clsat_satz_ok = numpy.where((AZIMUTH_RANGE[0] <= clsatObj.avhrr.satz) * \
                                (clsatObj.avhrr.satz <= AZIMUTH_RANGE[1]))
    
    clsat_satz_ok_bools = (AZIMUTH_RANGE[0] <= clsatObj.avhrr.satz) == \
                                (clsatObj.avhrr.satz <= AZIMUTH_RANGE[1])
    
    n_clsat_satz_ok = len(clsat_satz_ok[0])
    n_clsat_satz_ok_bools = len(clsat_satz_ok_bools)
    clsatObj.cloudsat.longitude = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.longitude,numpy.nan)
    clsatObj.cloudsat.latitude = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.latitude,numpy.nan)
    clsatObj.cloudsat.avhrr_linnum = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.avhrr_linnum,numpy.nan)
    clsatObj.cloudsat.avhrr_pixnum = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.avhrr_pixnum,numpy.nan)

    for i in range(clsatObj.cloudsat.cloud_mask.shape[0]):
        clsatObj.cloudsat.cloud_mask[i,:] = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.cloud_mask[i,:],numpy.nan)
        clsatObj.cloudsat.Radar_Reflectivity[i,:] = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.Radar_Reflectivity[i,:],numpy.nan)
        clsatObj.cloudsat.Height[i,:] = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.Height[i,:],numpy.nan)
    clsatObj.cloudsat.echo_top = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.echo_top,numpy.nan)
    clsatObj.cloudsat.SurfaceHeightBin = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.SurfaceHeightBin,numpy.nan)
    clsatObj.cloudsat.SurfaceHeightBin_fraction = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.SurfaceHeightBin_fraction,numpy.nan)
    
    clsatObj.cloudsat.elevation = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.elevation,numpy.nan)
    clsatObj.cloudsat.sec_1970 = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.sec_1970,numpy.nan)
    clsatObj.cloudsat.MODIS_Cloud_Fraction = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.MODIS_Cloud_Fraction,numpy.nan)
    clsatObj.cloudsat.MODIS_cloud_flag = numpy.where(clsat_satz_ok_bools,clsatObj.cloudsat.MODIS_cloud_flag,numpy.nan)

    # AVHRR - Cloudsat:
    clsatObj.avhrr.longitude = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.longitude,numpy.nan)
    clsatObj.avhrr.latitude = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.latitude,numpy.nan)
    clsatObj.avhrr.sec_1970 = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.sec_1970,numpy.nan)
    clsatObj.avhrr.cloudtype = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.cloudtype,numpy.nan)
    clsatObj.avhrr.ctth_height = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.ctth_height,numpy.nan)
    clsatObj.avhrr.ctth_pressure = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.ctth_pressure,numpy.nan)
    clsatObj.avhrr.ctth_temperature = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.ctth_temperature,numpy.nan)
    clsatObj.avhrr.bt11micron = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.bt11micron,numpy.nan)
    clsatObj.avhrr.bt12micron = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.bt12micron,numpy.nan)
    clsatObj.avhrr.surftemp = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.surftemp,numpy.nan)
    clsatObj.avhrr.satz = numpy.where(clsat_satz_ok_bools,clsatObj.avhrr.satz,numpy.nan)
    return clsatObj

def CalipsoAvhrrSatz(caObj):
    import numpy
    from config import AZIMUTH_RANGE
    import pdb
    # Calipso:
    ca_satz_ok = numpy.where((AZIMUTH_RANGE[0] <= caObj.avhrr.satz) * \
                                (caObj.avhrr.satz <= AZIMUTH_RANGE[1]))
    ca_satz_ok_bools = (AZIMUTH_RANGE[0] <= caObj.avhrr.satz) == \
                                (caObj.avhrr.satz <= AZIMUTH_RANGE[1])                         
    n_ca_satz_ok = len(ca_satz_ok[0])
    n_ca_satz_ok_bools = len(ca_satz_ok_bools)
    caObj.calipso.longitude = numpy.where(ca_satz_ok_bools,caObj.calipso.longitude,numpy.nan)
    caObj.calipso.latitude = numpy.where(ca_satz_ok_bools,caObj.calipso.latitude,numpy.nan)
    caObj.calipso.avhrr_linnum = numpy.where(ca_satz_ok_bools,caObj.calipso.avhrr_linnum,numpy.nan)
    caObj.calipso.avhrr_pixnum = numpy.where(ca_satz_ok_bools,caObj.calipso.avhrr_pixnum,numpy.nan)
    
    caObj.calipso.cloud_fraction = numpy.where(ca_satz_ok_bools,caObj.calipso.cloud_fraction,numpy.nan)
    caObj.calipso.elevation = numpy.where(ca_satz_ok_bools,caObj.calipso.elevation,numpy.nan)
    caObj.calipso.number_of_layers_found = numpy.where(ca_satz_ok_bools,caObj.calipso.number_of_layers_found,numpy.nan)
    
    caObj.calipso.igbp = numpy.where(ca_satz_ok_bools,caObj.calipso.igbp,numpy.nan)
    caObj.calipso.nsidc = numpy.where(ca_satz_ok_bools,caObj.calipso.nsidc  ,numpy.nan)  
    caObj.calipso.sec_1970 = numpy.where(ca_satz_ok_bools,caObj.calipso.sec_1970,numpy.nan)  

    for j in range(caObj.calipso.cloud_top_profile.shape[0]):
        caObj.calipso.cloud_top_profile[j,:] = numpy.where(ca_satz_ok_bools, caObj.calipso.cloud_top_profile[j,:], numpy.nan)
        caObj.calipso.cloud_base_profile[j,:] = numpy.where(ca_satz_ok_bools,caObj.calipso.cloud_base_profile[j,:], numpy.nan)
        caObj.calipso.cloud_mid_temperature[j,:] = numpy.where(ca_satz_ok_bools, caObj.calipso.cloud_mid_temperature[j,:], numpy.nan)
        try:
            caObj.calipso.feature_classification_flags[j,:] = numpy.where(ca_satz_ok_bools,caObj.calipso.feature_classification_flags[j,:],numpy.nan)
        except:
            print "No feature_classification_flags array in file!"
            pass
    # AVHRR - Calipso:
    caObj.avhrr.longitude = numpy.where(ca_satz_ok_bools,caObj.avhrr.longitude,numpy.nan)
    caObj.avhrr.latitude = numpy.where(ca_satz_ok_bools,caObj.avhrr.latitude,numpy.nan)
    caObj.avhrr.sec_1970 = numpy.where(ca_satz_ok_bools,caObj.avhrr.sec_1970,numpy.nan)
    caObj.avhrr.cloudtype = numpy.where(ca_satz_ok_bools,caObj.avhrr.cloudtype,numpy.nan)
    caObj.avhrr.ctth_height = numpy.where(ca_satz_ok_bools,caObj.avhrr.ctth_height,numpy.nan)
    caObj.avhrr.ctth_pressure = numpy.where(ca_satz_ok_bools,caObj.avhrr.ctth_pressure,numpy.nan)
    caObj.avhrr.ctth_temperature = numpy.where(ca_satz_ok_bools,caObj.avhrr.ctth_temperature,numpy.nan)
    caObj.avhrr.bt11micron = numpy.where(ca_satz_ok_bools,caObj.avhrr.bt11micron,numpy.nan)
    caObj.avhrr.bt12micron = numpy.where(ca_satz_ok_bools,caObj.avhrr.bt12micron,numpy.nan)
    caObj.avhrr.satz = numpy.where(ca_satz_ok_bools,caObj.avhrr.satz,numpy.nan)

    if caObj.avhrr.surftemp == None:        
        print("No Surftemp. Continue")
    else:
        caObj.avhrr.surftemp = numpy.where(ca_satz_ok_bools,caObj.avhrr.surftemp,numpy.nan)

    return caObj





















































