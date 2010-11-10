# Program cloudsat_calipso_avhrr_prepare.py
# python cloudsat_calipso_avhrr_prepare.py 








def CloudsatCloudOpticalDepth(cloud_top, cloud_base, optical_depth):
    import numpy
    import pdb
    
    
    new_cloud_top = numpy.ones(cloud_top.shape,'d')*numpy.min(cloud_top)
    new_cloud_base = numpy.ones(cloud_base.shape,'d')*numpy.min(cloud_base)
    
    
    for i in range(optical_depth.shape[1]):
        depthsum = 0 #Used to sum the optical_depth        
        for j in range(optical_depth.shape[0]):
            # Just stops the for loop when there are no more valid value 
            if optical_depth[j,i] < 0: 
                break
            
            depthsum = depthsum + optical_depth[j,i]
            # Removes the cloud values for all pixels that have a optical depth (integrated from the top) below one and moves the first valid value to the first column and so on.
            if depthsum >=1:
                new_cloud_top[0:(optical_depth.shape[0]-j),i] = cloud_top[j:,i]
                new_cloud_base[0:(optical_depth.shape[0]-j),i] = cloud_base[j:,i]                               
                break
               
    return new_cloud_top,new_cloud_base      
                
#---------------------------------------------------------------------
               
def CloudsatCalipsoAvhrrSatz1km(clsatObj,caObj):
    import pdb
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
    caObj.avhrr.surftemp = numpy.where(ca_satz_ok_bools,caObj.avhrr.surftemp,numpy.nan)
    caObj.avhrr.satz = numpy.where(ca_satz_ok_bools,caObj.avhrr.satz,numpy.nan)

    return clsatObj, caObj





















































