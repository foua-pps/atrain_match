#Program cloudsat_calipso_avhrr_statistics.py
import config

def CalculateStatistics(mode, clsatObj, statfile, caObj, cal_MODIS_cflag,
                        cal_vert_feature, avhrr_ctth_csat_ok, data_ok,
                        cal_data_ok, avhrr_ctth_cal_ok, caliop_max_height,
                        process_calipso_ok):
    import Scientific.Statistics
    import numpy
#    import numpy.oldnumpy as Numeric
    import pdb
    # First prepare possible subsetting of CALIOP datasets according to NSIDC and IGBP surface types

    if mode == "EMISSFILT":
        emissfilt_calipso_ok = process_calipso_ok 

    if mode is 'ICE_COVER_SEA':
##         cal_subset = numpy.logical_and(numpy.logical_and(numpy.less(caObj.calipso.nsidc,100),numpy.greater(caObj.calipso.nsidc,10)),numpy.equal(caObj.calipso.igbp,17))
        # Unfortunately, the above formulation used for ORR-B excluded the case when ice-cover was exactly 100 %!!! Very embarrassing!/KG
        cal_subset = numpy.logical_and(numpy.logical_and(numpy.less_equal(caObj.calipso.nsidc,100),numpy.greater(caObj.calipso.nsidc,10)),numpy.equal(caObj.calipso.igbp,17))
    elif mode is 'ICE_FREE_SEA':
        cal_subset = numpy.logical_and(numpy.equal(caObj.calipso.nsidc,0),numpy.equal(caObj.calipso.igbp,17))
    elif mode is 'SNOW_COVER_LAND':
        cal_subset = numpy.logical_and(numpy.logical_and(numpy.less(caObj.calipso.nsidc,104),numpy.greater(caObj.calipso.nsidc,10)),numpy.not_equal(caObj.calipso.igbp,17))
        # Notice that some uncertainty remains about the meaning of IGBP category 15 = "snow and ice". Can this possibly include also the Arctic ice sheet? We hope that it is not!!! However, if it is, the whole classification here might be wrong since this will affect also the definition of IGBP category 17./KG 
    elif mode is 'SNOW_FREE_LAND':
        cal_subset = numpy.logical_and(numpy.equal(caObj.calipso.nsidc,0),numpy.not_equal(caObj.calipso.igbp,17))
    elif mode is 'COASTAL_ZONE':
        cal_subset = numpy.equal(caObj.calipso.nsidc,255)
        
        
    # CLOUD MASK EVALUATION
    #=======================
    
    # CORRELATION CLOUD MASK: CLOUDSAT - AVHRR

    #print "------------------------------------"
    #print "STATISTICS CLOUD MASK: CLOUDSAT - AVHRR" 
    
    
    dummy=clsatObj.cloudsat.latitude.shape[0]
    pixel_position=numpy.arange(dummy)
    cloudsat_cloud_mask=clsatObj.cloudsat.cloud_mask
    cloudsat_cloud_mask=numpy.greater_equal(cloudsat_cloud_mask, config.CLOUDSAT_CLOUDY_THR)
    
    cloudsat_cloud_fraction=numpy.zeros(len(pixel_position))
    
    
    sum_cloudsat_cloud_mask=sum(cloudsat_cloud_mask)
        
    for idx in range (len(pixel_position)):
        if sum_cloudsat_cloud_mask[idx] > 2: # requires at least two cloudy bins
            cloudsat_cloud_fraction[idx]=1
                                        
    cloudsat_clear =  numpy.less(cloudsat_cloud_fraction,1)
    cloudsat_cloudy = numpy.greater_equal(cloudsat_cloud_fraction,1)
    pps_clear = numpy.logical_and(numpy.less_equal(clsatObj.avhrr.cloudtype,4),numpy.greater(clsatObj.avhrr.cloudtype,0))
    pps_cloudy = numpy.logical_and(numpy.greater(clsatObj.avhrr.cloudtype,4),numpy.less(clsatObj.avhrr.cloudtype,20))

    n_clear_clear = numpy.repeat(pps_clear,numpy.logical_and(cloudsat_clear,pps_clear)).shape[0]
    n_cloudy_cloudy = numpy.repeat(pps_cloudy,numpy.logical_and(cloudsat_cloudy,pps_cloudy)).shape[0]
    n_clear_cloudy = numpy.repeat(pps_cloudy,numpy.logical_and(cloudsat_clear,pps_cloudy)).shape[0]
    n_cloudy_clear = numpy.repeat(pps_clear,numpy.logical_and(cloudsat_cloudy,pps_clear)).shape[0]
    
    nclear = numpy.repeat(cloudsat_clear,cloudsat_clear).shape[0]
    ncloudy = numpy.repeat(cloudsat_cloudy,cloudsat_cloudy).shape[0]
    
    #print "Number of clear points (CLOUDSAT): ",nclear
    #print "Number of cloudy points (CLOUDSAT): ",ncloudy
    #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
    #print "Cloudsat-Clear PPS-Cloudy, Cloudsat-Cloudy PPS-Clear",n_clear_cloudy,n_cloudy_clear
    
    if ncloudy > 0:
        pod_cloudy = float(n_cloudy_cloudy)/ncloudy
        far_cloudy = float(n_cloudy_clear)/ncloudy     #Actually not correct computation of FAR!
                                                        #Should be divided by number of PPS cloudy!!!
    else:
        pod_cloudy = -9.0
        far_cloudy = -9.0
    if nclear > 0:
        pod_clear = float(n_clear_clear)/nclear
        far_clear = float(n_clear_cloudy)/nclear
    else:
        pod_clear = -9.0
        far_clear = -9.0
    
    #print "POD-Cloudy: ",pod_cloudy
    #print "POD-Clear: ",pod_clear
    #print "FAR-Cloudy: ",far_cloudy
    #print "FAR-Clear: ",far_clear    
    #print "-----------------------------------------"
    mean_cloudsat=((n_clear_clear+n_clear_cloudy)*0.0 + (n_cloudy_clear+n_cloudy_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
    mean_pps=((n_clear_clear+n_cloudy_clear)*0.0 + (n_cloudy_cloudy+n_clear_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
    bias=mean_pps-mean_cloudsat
    statfile.write("CLOUD MASK CLOUDSAT-PPS TABLE: %s %s %s %s \n" % (n_clear_clear,n_clear_cloudy,n_cloudy_clear,n_cloudy_cloudy))
    statfile.write("CLOUD MASK CLOUDSAT-PPS PROB: %f %f %f %f %f \n" % (pod_cloudy,pod_clear,far_cloudy,far_clear,bias))
    

    # CORRELATION CLOUD MASK: CLOUDSAT - MODIS

    #print "------------------------------------"
    #print "STATISTICS CLOUD MASK: CLOUDSAT - MODIS" 

    modis_clear = numpy.logical_or(numpy.equal(clsatObj.cloudsat.MODIS_cloud_flag,1),
                                        numpy.equal(clsatObj.cloudsat.MODIS_cloud_flag,0))
    modis_cloudy = numpy.logical_or(numpy.equal(clsatObj.cloudsat.MODIS_cloud_flag,3),
                                        numpy.equal(clsatObj.cloudsat.MODIS_cloud_flag,2))

    n_clear_clear = numpy.repeat(modis_clear,numpy.logical_and(cloudsat_clear,modis_clear)).shape[0]
    n_cloudy_cloudy = numpy.repeat(modis_cloudy,numpy.logical_and(cloudsat_cloudy,modis_cloudy)).shape[0]
    n_clear_cloudy = numpy.repeat(modis_cloudy,numpy.logical_and(cloudsat_clear,modis_cloudy)).shape[0]
    n_cloudy_clear = numpy.repeat(modis_clear,numpy.logical_and(cloudsat_cloudy,modis_clear)).shape[0]
    
    nclear = numpy.repeat(cloudsat_clear,cloudsat_clear).shape[0]
    ncloudy = numpy.repeat(cloudsat_cloudy,cloudsat_cloudy).shape[0]
    
    #print "Number of clear points (CLOUDSAT): ",nclear
    #print "Number of cloudy points (CLOUDSAT): ",ncloudy
    #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
    #print "Cloudsat-Clear MODIS-Cloudy, Cloudsat-Cloudy MODIS-Clear",n_clear_cloudy,n_cloudy_clear
    
    if ncloudy > 0:
        pod_cloudy = float(n_cloudy_cloudy)/ncloudy
        far_cloudy = float(n_cloudy_clear)/ncloudy
    else:
        pod_cloudy = -9.0
        far_cloudy = -9.0
    if nclear > 0:
        pod_clear = float(n_clear_clear)/nclear
        far_clear = float(n_clear_cloudy)/nclear
    else:
        pod_clear = -9.0
        far_clear = -9.0
    
    #print "POD-Cloudy: ",pod_cloudy
    #print "POD-Clear: ",pod_clear
    #print "FAR-Cloudy: ",far_cloudy
    #print "FAR-Clear: ",far_clear    
    #print "-----------------------------------------"
    mean_cloudsat=((n_clear_clear+n_clear_cloudy)*0.0 + (n_cloudy_clear+n_cloudy_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
    mean_modis=((n_clear_clear+n_cloudy_clear)*0.0 + (n_cloudy_cloudy+n_clear_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
    bias=mean_modis-mean_cloudsat
    statfile.write("CLOUD MASK CLOUDSAT-MODIS TABLE: %s %s %s %s \n" % (n_clear_clear,n_clear_cloudy,n_cloudy_clear,n_cloudy_cloudy))
    statfile.write("CLOUD MASK CLOUDSAT-MODIS PROB: %f %f %f %f %f \n" % (pod_cloudy,pod_clear,far_cloudy,far_clear,bias))
    
    
    # CORRELATION CLOUD MASK: CALIOP - AVHRR

    if mode is 'EMISSFILT':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),emissfilt_calipso_ok)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),emissfilt_calipso_ok)
        pps_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),emissfilt_calipso_ok)
        pps_cloudy = numpy.logical_and(numpy.greater(caObj.avhrr.cloudtype,4),emissfilt_calipso_ok)

        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - AVHRR - Emissivity filtered limit= ", EMISS_LIMIT
    elif mode is 'ICE_COVER_SEA':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        pps_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),cal_subset)
        pps_cloudy = numpy.logical_and(numpy.greater(caObj.avhrr.cloudtype,4),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - AVHRR - Exclusively over ocean ice:  "
    elif mode is 'ICE_FREE_SEA':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        pps_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),cal_subset)
        pps_cloudy = numpy.logical_and(numpy.greater(caObj.avhrr.cloudtype,4),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - AVHRR - Exclusively over ice-free ocean:  "
    elif mode is 'SNOW_COVER_LAND':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        pps_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),cal_subset)
        pps_cloudy = numpy.logical_and(numpy.greater(caObj.avhrr.cloudtype,4),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - AVHRR - Exclusively over snow-covered land:  "
    elif mode is 'SNOW_FREE_LAND':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        pps_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),cal_subset)
        pps_cloudy = numpy.logical_and(numpy.greater(caObj.avhrr.cloudtype,4),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - AVHRR - Exclusively over snow-free land:  "
    elif mode is 'COASTAL_ZONE':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        pps_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),cal_subset)
        pps_cloudy = numpy.logical_and(numpy.greater(caObj.avhrr.cloudtype,4),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - AVHRR - Exclusively over coastal zone:  "
    else:
        calipso_clear = numpy.less(caObj.calipso.cloud_fraction,0.34)
        calipso_cloudy = numpy.greater(caObj.calipso.cloud_fraction,0.66)
        pps_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0))
        pps_cloudy = numpy.logical_and(numpy.greater(caObj.avhrr.cloudtype,4),numpy.less(caObj.avhrr.cloudtype,20))
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - AVHRR"

    n_clear_clear = numpy.repeat(pps_clear,numpy.logical_and(calipso_clear,pps_clear)).shape[0]
    n_cloudy_cloudy = numpy.repeat(pps_cloudy,numpy.logical_and(calipso_cloudy,pps_cloudy)).shape[0]
    n_clear_cloudy = numpy.repeat(pps_cloudy,numpy.logical_and(calipso_clear,pps_cloudy)).shape[0]
    n_cloudy_clear = numpy.repeat(pps_clear,numpy.logical_and(calipso_cloudy,pps_clear)).shape[0]
    nclear = numpy.repeat(calipso_clear,calipso_clear).shape[0]
    ncloudy = numpy.repeat(calipso_cloudy,calipso_cloudy).shape[0]
    
    
    #print "Number of clear points (CALIPSO): ",nclear
    #print "Number of cloudy points (CALIPSO): ",ncloudy
    #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
    #print "Calipso-Clear PPS-Cloudy, Calipso-Cloudy PPS-Clear",n_clear_cloudy,n_cloudy_clear
    
    if ncloudy > 0:
        pod_cloudy = float(n_cloudy_cloudy)/ncloudy
        far_cloudy = float(n_cloudy_clear)/ncloudy
    else:
        pod_cloudy = -9.0
        far_cloudy = -9.0
    if nclear > 0:
        pod_clear = float(n_clear_clear)/nclear
        far_clear = float(n_clear_cloudy)/nclear
    else:
        pod_clear = -9.0
        far_clear = -9.0
    
    #print "POD-Cloudy: ",pod_cloudy
    #print "POD-Clear: ",pod_clear
    #print "FAR-Cloudy: ",far_cloudy
    #print "FAR-Clear: ",far_clear    
    #print "-----------------------------------------"
    if (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy) > 0:
        mean_caliop=((n_clear_clear+n_clear_cloudy)*0.0 + (n_cloudy_clear+n_cloudy_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        mean_pps=((n_clear_clear+n_cloudy_clear)*0.0 + (n_cloudy_cloudy+n_clear_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        bias=mean_pps-mean_caliop
    else:
        bias = -9.0
    statfile.write("CLOUD MASK CALIOP-PPS TABLE: %s %s %s %s \n" % (n_clear_clear,n_clear_cloudy,n_cloudy_clear,n_cloudy_cloudy))
    statfile.write("CLOUD MASK CALIOP-PPS PROB: %f %f %f %f %f \n" % (pod_cloudy,pod_clear,far_cloudy,far_clear,bias))

    # CORRELATION CLOUD MASK: CALIOP - MODIS

    if mode is 'EMISSFILT':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),emissfilt_calipso_ok)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),emissfilt_calipso_ok)
        modis_clear = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,1),
                                            numpy.equal(cal_MODIS_cflag,0)),emissfilt_calipso_ok)
        modis_cloudy = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,3),
                                            numpy.equal(cal_MODIS_cflag,2)),emissfilt_calipso_ok)

        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - MODIS - Emissivity filtered limit= ", EMISS_LIMIT
    elif mode is 'ICE_COVER_SEA':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        modis_clear = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,1),
                                            numpy.equal(cal_MODIS_cflag,0)),cal_subset)
        modis_cloudy = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,3),
                                            numpy.equal(cal_MODIS_cflag,2)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - MODIS - Exclusively over ocean ice:  "
    elif mode is 'ICE_FREE_SEA':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        modis_clear = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,1),
                                            numpy.equal(cal_MODIS_cflag,0)),cal_subset)
        modis_cloudy = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,3),
                                            numpy.equal(cal_MODIS_cflag,2)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - MODIS - Exclusively over ice-free ocean:  "
    elif mode is 'SNOW_COVER_LAND':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        modis_clear = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,1),
                                            numpy.equal(cal_MODIS_cflag,0)),cal_subset)
        modis_cloudy = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,3),
                                            numpy.equal(cal_MODIS_cflag,2)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - MODIS - Exclusively over snow-covered land:  "
    elif mode is 'SNOW_FREE_LAND':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        modis_clear = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,1),
                                            numpy.equal(cal_MODIS_cflag,0)),cal_subset)
        modis_cloudy = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,3),
                                            numpy.equal(cal_MODIS_cflag,2)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - MODIS - Exclusively over snow-free land:  "
    elif mode is 'COASTAL_ZONE':
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        calipso_cloudy = numpy.logical_and(numpy.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
        modis_clear = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,1),
                                            numpy.equal(cal_MODIS_cflag,0)),cal_subset)
        modis_cloudy = numpy.logical_and(numpy.logical_or(numpy.equal(cal_MODIS_cflag,3),
                                            numpy.equal(cal_MODIS_cflag,2)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - MODIS - Exclusively over coastal zone:  "
    else:
        calipso_clear = numpy.less(caObj.calipso.cloud_fraction,0.34)
        calipso_cloudy = numpy.greater(caObj.calipso.cloud_fraction,0.66)
        modis_clear = numpy.logical_or(numpy.equal(cal_MODIS_cflag,1),
                                            numpy.equal(cal_MODIS_cflag,0))
        modis_cloudy = numpy.logical_or(numpy.equal(cal_MODIS_cflag,3),
                                            numpy.equal(cal_MODIS_cflag,2))
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CALIOP - MODIS"

    n_clear_clear = numpy.repeat(modis_clear,numpy.logical_and(calipso_clear,modis_clear)).shape[0]
    n_cloudy_cloudy = numpy.repeat(modis_cloudy,numpy.logical_and(calipso_cloudy,modis_cloudy)).shape[0]
    n_clear_cloudy = numpy.repeat(modis_cloudy,numpy.logical_and(calipso_clear,modis_cloudy)).shape[0]
    n_cloudy_clear = numpy.repeat(modis_clear,numpy.logical_and(calipso_cloudy,modis_clear)).shape[0]
    nclear = numpy.repeat(calipso_clear,calipso_clear).shape[0]
    ncloudy = numpy.repeat(calipso_cloudy,calipso_cloudy).shape[0]
        
    #print "Number of clear points (CALIPSO): ",nclear
    #print "Number of cloudy points (CALIPSO): ",ncloudy
    #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
    #print "Calipso-Clear MODIS-Cloudy, Calipso-Cloudy MODIS-Clear",n_clear_cloudy,n_cloudy_clear
    
    if ncloudy > 0:
        pod_cloudy = float(n_cloudy_cloudy)/ncloudy
        far_cloudy = float(n_cloudy_clear)/ncloudy
    else:
        pod_cloudy = -9.0
        far_cloudy = -9.0
    if nclear > 0:
        pod_clear = float(n_clear_clear)/nclear
        far_clear = float(n_clear_cloudy)/nclear
    else:
        pod_clear = -9.0
        far_clear = -9.0
    
    #print "POD-Cloudy: ",pod_cloudy
    #print "POD-Clear: ",pod_clear
    #print "FAR-Cloudy: ",far_cloudy
    #print "FAR-Clear: ",far_clear    
    #print "-----------------------------------------"
    if (n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy) > 0:
        mean_caliop=((n_clear_clear+n_clear_cloudy)*0.0 + (n_cloudy_clear+n_cloudy_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        mean_modis=((n_clear_clear+n_cloudy_clear)*0.0 + (n_cloudy_cloudy+n_clear_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        bias=mean_modis-mean_caliop
    else:
        bias=-9.0
    statfile.write("CLOUD MASK CALIOP-MODIS TABLE: %s %s %s %s \n" % (n_clear_clear,n_clear_cloudy,n_cloudy_clear,n_cloudy_cloudy))
    statfile.write("CLOUD MASK CALIOP-MODIS PROB: %f %f %f %f %f \n" % (pod_cloudy,pod_clear,far_cloudy,far_clear,bias))
        

    # CLOUD TYPE EVALUATION - Based exclusively on CALIPSO data (Vertical Feature Mask)
    # =======================


    if mode is 'EMISSFILT':
        calipso_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3)),emissfilt_calipso_ok)
        calipso_medium = numpy.logical_and(numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5)),emissfilt_calipso_ok)
        calipso_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7)),emissfilt_calipso_ok)
        avhrr_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,5),numpy.less_equal(caObj.avhrr.cloudtype,8)),emissfilt_calipso_ok)
        avhrr_medium = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,9),numpy.less_equal(caObj.avhrr.cloudtype,10)),emissfilt_calipso_ok)
        avhrr_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,11),numpy.less_equal(caObj.avhrr.cloudtype,18)),emissfilt_calipso_ok)
        avhrr_frac = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,19),numpy.less_equal(caObj.avhrr.cloudtype,19)),emissfilt_calipso_ok)
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),emissfilt_calipso_ok)
        avhrr_clear = numpy.logical_and(numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0)),emissfilt_calipso_ok)
        #print "STATISTICS CLOUD TYPE: CALIOP - AVHRR - Emissivity filtered limit= ", EMISS_LIMIT
    elif mode is 'ICE_COVER_SEA':
        calipso_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3)),cal_subset)
        calipso_medium = numpy.logical_and(numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5)),cal_subset)
        calipso_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7)),cal_subset)
        avhrr_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,5),numpy.less_equal(caObj.avhrr.cloudtype,8)),cal_subset)
        avhrr_medium = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,9),numpy.less_equal(caObj.avhrr.cloudtype,10)),cal_subset)
        avhrr_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,11),numpy.less_equal(caObj.avhrr.cloudtype,18)),cal_subset)
        avhrr_frac = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,19),numpy.less_equal(caObj.avhrr.cloudtype,19)),cal_subset)
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        avhrr_clear = numpy.logical_and(numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD TYPE: CALIOP - AVHRR - Exclusively over ocean ice:  "
    elif mode is 'ICE_FREE_SEA':
        calipso_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3)),cal_subset)
        calipso_medium = numpy.logical_and(numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5)),cal_subset)
        calipso_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7)),cal_subset)
        avhrr_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,5),numpy.less_equal(caObj.avhrr.cloudtype,8)),cal_subset)
        avhrr_medium = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,9),numpy.less_equal(caObj.avhrr.cloudtype,10)),cal_subset)
        avhrr_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,11),numpy.less_equal(caObj.avhrr.cloudtype,18)),cal_subset)
        avhrr_frac = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,19),numpy.less_equal(caObj.avhrr.cloudtype,19)),cal_subset)
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        avhrr_clear = numpy.logical_and(numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD TYPE: CALIOP - AVHRR - Exclusively over ice-free ocean:  "
    elif mode is 'SNOW_COVER_LAND':
        calipso_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3)),cal_subset)
        calipso_medium = numpy.logical_and(numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5)),cal_subset)
        calipso_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7)),cal_subset)
        avhrr_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,5),numpy.less_equal(caObj.avhrr.cloudtype,8)),cal_subset)
        avhrr_medium = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,9),numpy.less_equal(caObj.avhrr.cloudtype,10)),cal_subset)
        avhrr_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,11),numpy.less_equal(caObj.avhrr.cloudtype,18)),cal_subset)
        avhrr_frac = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,19),numpy.less_equal(caObj.avhrr.cloudtype,19)),cal_subset)
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        avhrr_clear = numpy.logical_and(numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD TYPE: CALIOP - AVHRR - Exclusively over snow-covered land:  "
    elif mode is 'SNOW_FREE_LAND':
        calipso_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3)),cal_subset)
        calipso_medium = numpy.logical_and(numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5)),cal_subset)
        calipso_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7)),cal_subset)
        avhrr_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,5),numpy.less_equal(caObj.avhrr.cloudtype,8)),cal_subset)
        avhrr_medium = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,9),numpy.less_equal(caObj.avhrr.cloudtype,10)),cal_subset)
        avhrr_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,11),numpy.less_equal(caObj.avhrr.cloudtype,18)),cal_subset)
        avhrr_frac = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,19),numpy.less_equal(caObj.avhrr.cloudtype,19)),cal_subset)
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        avhrr_clear = numpy.logical_and(numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD TYPE: CALIOP - AVHRR - Exclusively over snow-free land:  "
    elif mode is 'COASTAL_ZONE':
        calipso_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3)),cal_subset)
        calipso_medium = numpy.logical_and(numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5)),cal_subset)
        calipso_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7)),cal_subset)
        avhrr_low = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,5),numpy.less_equal(caObj.avhrr.cloudtype,8)),cal_subset)
        avhrr_medium = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,9),numpy.less_equal(caObj.avhrr.cloudtype,10)),cal_subset)
        avhrr_high = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,11),numpy.less_equal(caObj.avhrr.cloudtype,18)),cal_subset)
        avhrr_frac = numpy.logical_and(numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,19),numpy.less_equal(caObj.avhrr.cloudtype,19)),cal_subset)
        calipso_clear = numpy.logical_and(numpy.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
        avhrr_clear = numpy.logical_and(numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0)),cal_subset)
        #print "------------------------------------"
        #print "STATISTICS CLOUD TYPE: CALIOP - AVHRR - Exclusively over coastal zone:  "
    else:       
        calipso_low = numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3))
        calipso_medium = numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5))
        calipso_high = numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7))
        avhrr_low = numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,5),numpy.less_equal(caObj.avhrr.cloudtype,8))
        avhrr_medium = numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,9),numpy.less_equal(caObj.avhrr.cloudtype,10))
        avhrr_high = numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,11),numpy.less_equal(caObj.avhrr.cloudtype,18))
        avhrr_frac = numpy.logical_and(numpy.greater_equal(caObj.avhrr.cloudtype,19),numpy.less_equal(caObj.avhrr.cloudtype,19))
        calipso_clear = numpy.less(caObj.calipso.cloud_fraction,0.34)
        avhrr_clear = numpy.logical_and(numpy.less_equal(caObj.avhrr.cloudtype,4),numpy.greater(caObj.avhrr.cloudtype,0))

        #print "STATISTICS CLOUD TYPE: CALIOP - AVHRR"

    # Notice that we have unfortunately changed order in notation compared to cloud mask
    # Here the PPS category is mentioned first and then the CALIOP category 

    n_low_low = numpy.repeat(avhrr_low,numpy.logical_and(calipso_low,avhrr_low)).shape[0]
    n_low_medium = numpy.repeat(avhrr_low,numpy.logical_and(calipso_medium,avhrr_low)).shape[0]
    n_low_high = numpy.repeat(avhrr_low,numpy.logical_and(calipso_high,avhrr_low)).shape[0]
    n_medium_low = numpy.repeat(avhrr_medium,numpy.logical_and(calipso_low,avhrr_medium)).shape[0]
    n_medium_medium = numpy.repeat(avhrr_medium,numpy.logical_and(calipso_medium,avhrr_medium)).shape[0]
    n_medium_high = numpy.repeat(avhrr_medium,numpy.logical_and(calipso_high,avhrr_medium)).shape[0]
    n_high_low = numpy.repeat(avhrr_high,numpy.logical_and(calipso_low,avhrr_high)).shape[0]
    n_high_medium = numpy.repeat(avhrr_high,numpy.logical_and(calipso_medium,avhrr_high)).shape[0]
    n_high_high = numpy.repeat(avhrr_high,numpy.logical_and(calipso_high,avhrr_high)).shape[0]
    n_frac_low = numpy.repeat(avhrr_frac,numpy.logical_and(calipso_low,avhrr_frac)).shape[0]
    n_frac_medium = numpy.repeat(avhrr_frac,numpy.logical_and(calipso_medium,avhrr_frac)).shape[0]
    n_frac_high = numpy.repeat(avhrr_frac,numpy.logical_and(calipso_high,avhrr_frac)).shape[0]

    nlow = numpy.repeat(calipso_low,calipso_low).shape[0]
    nmedium = numpy.repeat(calipso_medium,calipso_medium).shape[0]
    nhigh = numpy.repeat(calipso_high,calipso_high).shape[0]
        
    n_clear_low = numpy.repeat(avhrr_clear,numpy.logical_and(calipso_low,avhrr_clear)).shape[0]
    n_clear_medium = numpy.repeat(avhrr_clear,numpy.logical_and(calipso_medium,avhrr_clear)).shape[0]
    n_clear_high = numpy.repeat(avhrr_clear,numpy.logical_and(calipso_high,avhrr_clear)).shape[0]
    n_low_clear = numpy.repeat(avhrr_low,numpy.logical_and(calipso_clear,avhrr_low)).shape[0]
    n_medium_clear = numpy.repeat(avhrr_medium,numpy.logical_and(calipso_clear,avhrr_medium)).shape[0]
    n_high_clear = numpy.repeat(avhrr_high,numpy.logical_and(calipso_clear,avhrr_high)).shape[0]
    n_frac_clear = numpy.repeat(avhrr_frac,numpy.logical_and(calipso_clear,avhrr_frac)).shape[0]
    
    #print "Number of points with low clouds (CALIPSO): ",nlow
    #print "Number of points with medium clouds (CALIPSO): ",nmedium
    #print "Number of points with high clouds (CALIPSO): ",nhigh
    #print "CALIPSO low-AVHRR low,CALIPSO low-AVHRR medium, CALIPSO low-AVHRR high",n_low_low,n_medium_low,n_high_low
    #print "CALIPSO medium-AVHRR low,CALIPSO medium-AVHRR medium, CALIPSO medium-AVHRR high",n_low_medium,n_medium_medium,n_high_medium
    #print "CALIPSO high-AVHRR low,CALIPSO high-AVHRR medium, CALIPSO high-AVHRR high",n_low_high,n_medium_high,n_high_high
    #print "CALIPSO high-AVHRR frac,CALIPSO medium-AVHRR frac, CALIPSO low-AVHRR frac",n_frac_high,n_frac_medium,n_frac_low
    
    if (n_low_low+n_medium_low+n_high_low+n_frac_low) > 0:
        pod_low = float(n_low_low)/(n_low_low+n_medium_low+n_high_low+n_frac_low)
        far_low = float(n_medium_low+n_high_low+n_frac_low)/(n_low_low+n_medium_low+n_high_low+n_frac_low)
    else:
        pod_low = -9.0
        far_low = -9.0
    if (n_low_medium+n_medium_medium+n_high_medium+n_frac_medium) > 0:
        pod_medium = float(n_medium_medium)/(n_low_medium+n_medium_medium+n_high_medium+n_frac_medium)
        far_medium = float(n_low_medium+n_high_medium+n_frac_medium)/(n_low_medium+n_medium_medium+n_high_medium+n_frac_medium)
    else:
        pod_medium =-9.0
        far_medium =-9.0
    if (n_low_high+n_medium_high+n_high_high+n_frac_high) > 0:
        pod_high = float(n_high_high)/(n_low_high+n_medium_high+n_high_high+n_frac_high)
        far_high = float(n_low_high+n_medium_high+n_frac_high)/(n_low_high+n_medium_high+n_high_high+n_frac_high)
    else:
        pod_high =-9.0
        far_high =-9.0

    #print "POD-Low: ",pod_low
    #print "POD-Medium: ",pod_medium
    #print "POD-High: ",pod_high
    #print "FAR-Low: ",far_low
    #print "FAR-Medium: ",far_medium    
    #print "FAR-High: ",far_high    
    #print "-----------------------------------------"
    statfile.write("CLOUD TYPE CALIOP-PPS TABLE: %s %s %s %s %s %s %s %s %s %s %s %s \n" % (n_low_low,n_low_medium,n_low_high,n_medium_low,n_medium_medium,n_medium_high,n_high_low,n_high_medium,n_high_high,n_frac_low,n_frac_medium,n_frac_high))
    statfile.write("CLOUD TYPE CALIOP-PPS PROB: %f %f %f %f %f %f \n" % (pod_low,pod_medium,pod_high,far_low,far_medium,far_high))
    statfile.write("CLOUD TYPE CALIOP-PPS TABLE MISSED: %s %s %s %s %s %s %s \n" % (n_clear_low,n_clear_medium,n_clear_high,n_low_clear,n_medium_clear,n_high_clear,n_frac_clear))
            

    # CLOUD TOP EVALUATION
    #=======================
    
    # CORRELATION: CLOUDSAT - AVHRR HEIGHT

    #print "STATISTICS CLOUD TOP HEIGHT: CLOUDSAT - AVHRR"

    #if mode not in config.PLOT_MODES:
    if True: # TODO: The above if statement seems to have lost its meaning...
        dummy=clsatObj.cloudsat.latitude.shape[0]
        pixel_position_plain=numpy.arange(dummy)
        clsat_max_height = numpy.repeat(clsatObj.cloudsat.Height[124,::],data_ok)
#        clsat_max_height = numpy.zeros(dummy, 'f')

        for i in range(125):
            height = clsatObj.cloudsat.Height[i,::]
            cmask_ok = clsatObj.cloudsat.cloud_mask[i,::]
            for idx in range(len(pixel_position_plain)):
                #nidx = int(cmask_ok[idx]+0.5)-1
                nidx = int(cmask_ok[idx]+0.5)
                #if nidx < 0:
                if nidx == 0:
                    continue
                #print idx,nidx,colors[nidx]
                #if height[idx] < 0 or height[idx] > MAXHEIGHT:
                if height[idx] < 240*4: # or height[idx] > config.MAXHEIGHT:
                    continue
                base_height = height[idx]-120
                top_height = height[idx]+120
                if nidx >= int(config.CLOUDSAT_CLOUDY_THR):
#                if nidx >= 20:
                    clsat_max_height[idx] = max(clsat_max_height[idx],top_height)  

    okarr = numpy.logical_and(numpy.greater(avhrr_ctth_csat_ok,0.0),data_ok)
    okarr = numpy.logical_and(okarr,numpy.greater(clsat_max_height,0.0))
    clsat_max_height = numpy.repeat(clsat_max_height[::],okarr)
    avhrr_height = numpy.repeat(avhrr_ctth_csat_ok[::],okarr)

    if len(avhrr_height) > 0:
        if len(avhrr_height) > 20:
            corr_cloudsat_avhrr = Scientific.Statistics.correlation(clsat_max_height,avhrr_height)
        else:
            corr_cloudsat_avhrr = -99.0
        #print "Correlation: cloudsat-avhrr = ",corr_cloudsat_avhrr
        diff = avhrr_height-clsat_max_height
        bias = Scientific.Statistics.mean(diff)
        #print "Mean difference cloudsat-avhrr = ",bias
        diff_squared = diff*diff
        RMS_difference = numpy.sqrt(Scientific.Statistics.mean(diff_squared))
        #print "RMS difference cloudsat-avhrr = ",RMS_difference
        diff_squared = (diff-bias)*(diff-bias)
        RMS_difference = numpy.sqrt(Scientific.Statistics.mean(diff_squared))
        #print "Bias-corrected RMS difference cloudsat-avhrr = ",RMS_difference
        #print "Number of matchups: ", len(avhrr_height)
        #print "Number of failing PPS CTTHs: ", len(numpy.repeat(clsatObj.avhrr.cloudtype[::],numpy.greater(clsatObj.avhrr.cloudtype,4)))-len(avhrr_height)
        #print

        statfile.write("CLOUD HEIGHT CLOUDSAT: %f %f %f %s %f \n" % (corr_cloudsat_avhrr,bias,RMS_difference,len(avhrr_height),sum(diff_squared)))
    else:
        statfile.write("CLOUD HEIGHT CLOUDSAT: -9.0 -9.0 -9.0 0 -9.0 \n")
                    
    # CORRELATION: CALIOP - AVHRR HEIGHT
    # FIRST TOTAL FIGURES


    if mode is 'EMISSFILT':
        okcaliop = numpy.logical_and(cal_data_ok,
                                        numpy.logical_and(numpy.greater(avhrr_ctth_cal_ok[::],0),emissfilt_calipso_ok))
        #print "-----------------------------------------"
        #print "STATISTICS CLOUD TOP HEIGHT: CALIOP - AVHRR - Emissivity filtered limit= ", EMISS_LIMIT
    elif mode is 'ICE_COVER_SEA':
        okcaliop = numpy.logical_and(cal_data_ok,
                                        numpy.logical_and(numpy.greater(avhrr_ctth_cal_ok[::],0),cal_subset))
        #print "------------------------------------"
        #print "STATISTICS CLOUD TOP HEIGHT: CALIOP - AVHRR - Exclusively over ocean ice:  "
    elif mode is 'ICE_FREE_SEA':
        okcaliop = numpy.logical_and(cal_data_ok,
                                        numpy.logical_and(numpy.greater(avhrr_ctth_cal_ok[::],0),cal_subset))
        #print "------------------------------------"
        #print "STATISTICS CLOUD TOP HEIGHT: CALIOP - AVHRR - Exclusively over ice-free ocean:  "
    elif mode is 'SNOW_COVER_LAND':
        okcaliop = numpy.logical_and(cal_data_ok,
                                        numpy.logical_and(numpy.greater(avhrr_ctth_cal_ok[::],0),cal_subset))
        #print "------------------------------------"
        #print "STATISTICS CLOUD TOP HEIGHT: CALIOP - AVHRR - Exclusively over snow-covered land:  "
    elif mode is 'SNOW_FREE_LAND':
        okcaliop = numpy.logical_and(cal_data_ok,
                                        numpy.logical_and(numpy.greater(avhrr_ctth_cal_ok[::],0),cal_subset))
        #print "------------------------------------"
        #print "STATISTICS CLOUD TOP HEIGHT: CALIOP - AVHRR - Exclusively over snow-free land:  "
    elif mode is 'COASTAL_ZONE':
        okcaliop = numpy.logical_and(cal_data_ok,
                                        numpy.logical_and(numpy.greater(avhrr_ctth_cal_ok[::],0),cal_subset))
        #print "------------------------------------"
        #print "STATISTICS CLOUD TOP HEIGHT: CALIOP - AVHRR - Exclusively over coastal zone:  "
    else:
        okcaliop = numpy.logical_and(cal_data_ok,
                                        numpy.greater(avhrr_ctth_cal_ok[::],0))
        #print "-----------------------------------------"
        #print "STATISTICS CLOUD TOP HEIGHT: CALIOP - AVHRR"

    #print "ALL CLOUDS:"
    avhrr_height_work = numpy.repeat(avhrr_ctth_cal_ok[::],okcaliop)
    caliop_max_height_work = numpy.repeat(caliop_max_height[::],okcaliop)
    if len(caliop_max_height_work) > 0:
        if len(avhrr_height_work) > 20:
            corr_caliop_avhrr = Scientific.Statistics.correlation(caliop_max_height_work,avhrr_height_work)
        else:
            corr_caliop_avhrr = -99.0
        diff = avhrr_height_work-caliop_max_height_work
        bias = Scientific.Statistics.mean(diff)
        diff_squared = diff*diff
        RMS_difference = numpy.sqrt(Scientific.Statistics.mean(diff_squared))
        diff_squared_biascorr = (diff-bias)*(diff-bias)
        RMS_difference_biascorr = numpy.sqrt(Scientific.Statistics.mean(diff_squared_biascorr))
    else:
        corr_caliop_avhrr = -9.0
        bias = -9.0
        RMS_difference = -9.0
        RMS_difference_biascorr = -9.0
        diff_squared_biascorr = numpy.array([-9.0])
        
    #print "Correlation: caliop-avhrr = ",corr_caliop_avhrr
    #print "Mean difference caliop-avhrr = ", bias
    #print "RMS difference caliop-avhrr = ",RMS_difference
    #print "Bias-corrected RMS difference caliop-avhrr = ",RMS_difference_biascorr
    #print "Number of matchups: ", len(avhrr_height_work)
    #print "Number of failing PPS CTTHs: ", len(numpy.repeat(caObj.avhrr.cloudtype[::],numpy.greater(caObj.avhrr.cloudtype,4)))-len(avhrr_height_work)
    statfile.write("CLOUD HEIGHT CALIOP ALL: %f %f %f %s %f \n" % (corr_caliop_avhrr,bias,RMS_difference,len(avhrr_height_work),sum(diff_squared_biascorr)))
    # THEN FOR LOW CLOUDS (VERTICAL FEATURE MASK CATEGORIES 0-3)
    cal_low_ok = numpy.logical_and(numpy.greater_equal(cal_vert_feature[::],0),numpy.less_equal(cal_vert_feature[::],3))
    cal_low_ok = numpy.logical_and(cal_low_ok,okcaliop)
    avhrr_height_work = numpy.repeat(avhrr_ctth_cal_ok[::],cal_low_ok)
    caliop_max_height_work = numpy.repeat(caliop_max_height[::],cal_low_ok)
    #print "LOW CLOUDS ( >680 hPa):"
    if len(caliop_max_height_work) > 0:
        if len(avhrr_height_work) > 20:
            corr_caliop_avhrr = Scientific.Statistics.correlation(caliop_max_height_work,avhrr_height_work)
        else:
            corr_caliop_avhrr = -99.0
        diff = avhrr_height_work-caliop_max_height_work
        bias = Scientific.Statistics.mean(diff)
        diff_squared = diff*diff
        RMS_difference = numpy.sqrt(Scientific.Statistics.mean(diff_squared))
        diff_squared_biascorr = (diff-bias)*(diff-bias)
        RMS_difference_biascorr = numpy.sqrt(Scientific.Statistics.mean(diff_squared_biascorr))
    else:
        corr_caliop_avhrr = -9.0
        bias = -9.0
        RMS_difference = -9.0
        RMS_difference_biascorr = -9.0
        diff_squared_biascorr = numpy.array([-9.0])
        
    #print "Correlation: caliop-avhrr = ",corr_caliop_avhrr
    #print "Mean difference caliop-avhrr = ", bias
    #print "RMS difference caliop-avhrr = ",RMS_difference
    #print "Bias-corrected RMS difference caliop-avhrr = ",RMS_difference_biascorr
    #print "Number of matchups: ", len(avhrr_height_work)
    #print "Number of failing PPS CTTHs: ", len(numpy.repeat(caObj.avhrr.cloudtype[::],numpy.greater(caObj.avhrr.cloudtype,4)))-len(avhrr_height_work)
    statfile.write("CLOUD HEIGHT CALIOP LOW: %f %f %f %s %f \n" % (corr_caliop_avhrr,bias,RMS_difference,len(avhrr_height_work),sum(diff_squared_biascorr)))
        
    # THEN FOR MEDIUM CLOUDS (VERTICAL FEATURE MASK CATEGORIES 4-5)
    cal_mid_ok = numpy.logical_and(numpy.greater(cal_vert_feature[::],3),numpy.less_equal(cal_vert_feature[::],5))
    cal_mid_ok = numpy.logical_and(cal_mid_ok,okcaliop)
    avhrr_height_work = numpy.repeat(avhrr_ctth_cal_ok[::],cal_mid_ok)
    caliop_max_height_work = numpy.repeat(caliop_max_height[::],cal_mid_ok)
    #print "MEDIUM CLOUDS ( 440-680 hPa):"
    if len(caliop_max_height_work) > 0:
        if len(avhrr_height_work) > 20:
            corr_caliop_avhrr = Scientific.Statistics.correlation(caliop_max_height_work,avhrr_height_work)
        else:
            corr_caliop_avhrr = -99.0
        diff = avhrr_height_work-caliop_max_height_work
        bias = Scientific.Statistics.mean(diff)
        diff_squared = diff*diff
        RMS_difference = numpy.sqrt(Scientific.Statistics.mean(diff_squared))
        diff_squared_biascorr = (diff-bias)*(diff-bias)
        RMS_difference_biascorr = numpy.sqrt(Scientific.Statistics.mean(diff_squared_biascorr))
    else:
        corr_caliop_avhrr = -9.0
        bias = -9.0
        RMS_difference = -9.0
        RMS_difference_biascorr = -9.0
        diff_squared_biascorr = numpy.array([-9.0])
        
    #print "Correlation: caliop-avhrr = ",corr_caliop_avhrr
    #print "Mean difference caliop-avhrr = ", bias
    #print "RMS difference caliop-avhrr = ",RMS_difference
    #print "Bias-corrected RMS difference caliop-avhrr = ",RMS_difference_biascorr
    #print "Number of matchups: ", len(avhrr_height_work)
    #print "Number of failing PPS CTTHs: ", len(numpy.repeat(caObj.avhrr.cloudtype[::],numpy.greater(caObj.avhrr.cloudtype,4)))-len(avhrr_height_work)
    statfile.write("CLOUD HEIGHT CALIOP MEDIUM: %f %f %f %s %f \n" % (corr_caliop_avhrr,bias,RMS_difference,len(avhrr_height_work),sum(diff_squared_biascorr)))
            
    # FINALLY FOR HIGH CLOUDS (VERTICAL FEATURE MASK CATEGORIES 6-7)
    cal_high_ok = numpy.logical_and(numpy.greater(cal_vert_feature[::],5),numpy.less_equal(cal_vert_feature[::],7))
    cal_high_ok = numpy.logical_and(cal_high_ok,okcaliop)
    avhrr_height_work = numpy.repeat(avhrr_ctth_cal_ok[::],cal_high_ok)
    caliop_max_height_work = numpy.repeat(caliop_max_height[::],cal_high_ok)
    #print "HIGH CLOUDS ( <440 hPa):"
    if len(caliop_max_height_work) > 0:
        if len(avhrr_height_work) > 20:
            corr_caliop_avhrr = Scientific.Statistics.correlation(caliop_max_height_work,avhrr_height_work)
        else:
            corr_caliop_avhrr = -99.0
        diff = avhrr_height_work-caliop_max_height_work
        bias = Scientific.Statistics.mean(diff)
        diff_squared = diff*diff
        RMS_difference = numpy.sqrt(Scientific.Statistics.mean(diff_squared))
        diff_squared_biascorr = (diff-bias)*(diff-bias)
        RMS_difference_biascorr = numpy.sqrt(Scientific.Statistics.mean(diff_squared_biascorr))
    else:
        corr_caliop_avhrr = -9.0
        bias = -9.0
        RMS_difference = -9.0
        RMS_difference_biascorr = -9.0
        diff_squared_biascorr = numpy.array([-9.0])
        
    #print "Correlation: caliop-avhrr = ",corr_caliop_avhrr
    #print "Mean difference caliop-avhrr = ", bias
    #print "RMS difference caliop-avhrr = ",RMS_difference
    #print "Bias-corrected RMS difference caliop-avhrr = ",RMS_difference_biascorr
    #print "Number of matchups: ", len(avhrr_height_work)
    #print "Number of failing PPS CTTHs: ", len(numpy.repeat(caObj.avhrr.cloudtype[::],numpy.greater(caObj.avhrr.cloudtype,4)))-len(avhrr_height_work)
    statfile.write("CLOUD HEIGHT CALIOP HIGH: %f %f %f %s %f \n" % (corr_caliop_avhrr,bias,RMS_difference,len(avhrr_height_work),sum(diff_squared_biascorr)))

    statfile.close()
    
