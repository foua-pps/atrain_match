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
    from cloudsat_calipso_avhrr_match import AZIMUTH_THR
    # CloudSat:
    clsat_satz_ok = numpy.where(clsatObj.avhrr.satz<=AZIMUTH_THR)
    n_clsat_satz_ok = len(clsat_satz_ok[0])
    clsatObj.cloudsat.longitude = clsatObj.cloudsat.longitude[clsat_satz_ok]
    clsatObj.cloudsat.latitude = clsatObj.cloudsat.latitude[clsat_satz_ok]
    clsatObj.cloudsat.avhrr_linnum = clsatObj.cloudsat.avhrr_linnum[clsat_satz_ok]
    clsatObj.cloudsat.avhrr_pixnum = clsatObj.cloudsat.avhrr_pixnum[clsat_satz_ok]
  
    clsatObj.cloudsat.cloud_mask = numpy.reshape(clsatObj.cloudsat.cloud_mask[:,clsat_satz_ok],(-1,n_clsat_satz_ok))
    clsatObj.cloudsat.Radar_Reflectivity = numpy.reshape(clsatObj.cloudsat.Radar_Reflectivity[:,clsat_satz_ok],(-1,n_clsat_satz_ok))
    clsatObj.cloudsat.Height = numpy.reshape(clsatObj.cloudsat.Height[:,clsat_satz_ok],(-1,n_clsat_satz_ok))
    clsatObj.cloudsat.echo_top = clsatObj.cloudsat.echo_top[clsat_satz_ok]
    clsatObj.cloudsat.SurfaceHeightBin = clsatObj.cloudsat.SurfaceHeightBin[clsat_satz_ok]
    clsatObj.cloudsat.SurfaceHeightBin_fraction = clsatObj.cloudsat.SurfaceHeightBin_fraction[clsat_satz_ok]
    
    clsatObj.cloudsat.elevation = clsatObj.cloudsat.elevation[clsat_satz_ok]
    clsatObj.cloudsat.sec_1970 = clsatObj.cloudsat.sec_1970[clsat_satz_ok]
    clsatObj.cloudsat.MODIS_Cloud_Fraction = clsatObj.cloudsat.MODIS_Cloud_Fraction[clsat_satz_ok]
    clsatObj.cloudsat.MODIS_cloud_flag = clsatObj.cloudsat.MODIS_cloud_flag[clsat_satz_ok]

    # AVHRR - Cloudsat:
    clsatObj.avhrr.longitude = clsatObj.avhrr.longitude[clsat_satz_ok]
    clsatObj.avhrr.latitude = clsatObj.avhrr.latitude[clsat_satz_ok]
    clsatObj.avhrr.sec_1970 = clsatObj.avhrr.sec_1970[clsat_satz_ok]
    clsatObj.avhrr.cloudtype = clsatObj.avhrr.cloudtype[clsat_satz_ok]
    clsatObj.avhrr.ctth_height = clsatObj.avhrr.ctth_height[clsat_satz_ok]
    clsatObj.avhrr.ctth_pressure = clsatObj.avhrr.ctth_pressure[clsat_satz_ok]
    clsatObj.avhrr.ctth_temperature = clsatObj.avhrr.ctth_temperature[clsat_satz_ok]
    clsatObj.avhrr.bt11micron = clsatObj.avhrr.bt11micron[clsat_satz_ok]
    clsatObj.avhrr.bt12micron = clsatObj.avhrr.bt12micron[clsat_satz_ok]
    clsatObj.avhrr.surftemp = clsatObj.avhrr.surftemp[clsat_satz_ok]
    clsatObj.avhrr.satz = clsatObj.avhrr.satz[clsat_satz_ok]   
    
    # Calipso:
    ca_satz_ok = numpy.where(caObj.avhrr.satz<=AZIMUTH_THR)
    n_ca_satz_ok = len(ca_satz_ok[0])
    caObj.calipso.longitude = caObj.calipso.longitude[ca_satz_ok]
    caObj.calipso.latitude = caObj.calipso.latitude[ca_satz_ok]
    caObj.calipso.avhrr_linnum = caObj.calipso.avhrr_linnum[ca_satz_ok]
    caObj.calipso.avhrr_pixnum = caObj.calipso.avhrr_pixnum[ca_satz_ok]
    
    caObj.calipso.cloud_fraction = caObj.calipso.cloud_fraction[ca_satz_ok]
    caObj.calipso.cloud_top_profile = numpy.reshape(caObj.calipso.cloud_top_profile[:,ca_satz_ok],(-1,n_ca_satz_ok))
    caObj.calipso.cloud_base_profile = numpy.reshape(caObj.calipso.cloud_base_profile[:,ca_satz_ok],(-1,n_ca_satz_ok))
    caObj.calipso.cloud_mid_temperature = numpy.reshape(caObj.calipso.cloud_mid_temperature[:,ca_satz_ok],(-1,n_ca_satz_ok))
    caObj.calipso.elevation = caObj.calipso.elevation[ca_satz_ok]
    caObj.calipso.number_of_layers_found = caObj.calipso.number_of_layers_found[ca_satz_ok]
    try:
        caObj.calipso.feature_classification_flags = numpy.reshape(caObj.calipso.feature_classification_flags[:,ca_satz_ok],(-1,n_ca_satz_ok))
    except:
        print "No feature_classification_flags array in file!"
        pass
    
    caObj.calipso.igbp = caObj.calipso.igbp[ca_satz_ok]
    caObj.calipso.nsidc = caObj.calipso.nsidc  [ca_satz_ok]  
    caObj.calipso.sec_1970 = caObj.calipso.sec_1970[ca_satz_ok]

    # AVHRR - Calipso:
    caObj.avhrr.longitude = caObj.avhrr.longitude[ca_satz_ok]
    caObj.avhrr.latitude = caObj.avhrr.latitude[ca_satz_ok]
    caObj.avhrr.sec_1970 = caObj.avhrr.sec_1970[ca_satz_ok]
    caObj.avhrr.cloudtype = caObj.avhrr.cloudtype[ca_satz_ok]
    caObj.avhrr.ctth_height = caObj.avhrr.ctth_height[ca_satz_ok]
    caObj.avhrr.ctth_pressure = caObj.avhrr.ctth_pressure[ca_satz_ok]
    caObj.avhrr.ctth_temperature = caObj.avhrr.ctth_temperature[ca_satz_ok]
    caObj.avhrr.bt11micron = caObj.avhrr.bt11micron[ca_satz_ok]
    caObj.avhrr.bt12micron = caObj.avhrr.bt12micron[ca_satz_ok]
    caObj.avhrr.surftemp = caObj.avhrr.surftemp[ca_satz_ok]
    caObj.avhrr.satz = caObj.avhrr.satz[ca_satz_ok]

    return clsatObj, caObj





















































