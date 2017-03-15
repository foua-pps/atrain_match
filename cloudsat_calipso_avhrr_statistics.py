#Program cloudsat_calipso_avhrr_statistics.py
import config
import sys
import logging
logger = logging.getLogger(__name__)
import numpy as np
from get_flag_info import (get_semi_opaque_info_pps2014, 
                           get_semi_opaque_info_pps2012,
                           get_sunglint_info_pps2014,
                           get_high_terrain_info_pps2014,
                           get_mountin_info_pps2014,
                           get_inversion_info_pps2014,
                           get_land_coast_sea_info_pps2014,
                           get_land_coast_sea_info_pps2012,
                           get_land_coast_sea_info_cci2014,
                           get_ice_info_pps2014,
                           get_ice_info_pps2012,
                           get_day_night_twilight_info_pps2014,
                           get_day_night_twilight_info_pps2012,
                           get_day_night_twilight_info_cci2014,
                           get_day_night_twilight_info_maia,
                           get_sunglint_info_pps2012,
                           get_mountin_info_pps2012,
                           get_inversion_info_pps2012)

def calculate_ctth_stats(okcaliop, avhrr_ctth_cal_ok,caliop_max_height):
    avhrr_height_work = np.repeat(avhrr_ctth_cal_ok[::],okcaliop)
    caliop_max_height_work = np.repeat(caliop_max_height[::],okcaliop)
    if len(caliop_max_height_work) > 0:
        if len(avhrr_height_work) > 20:
            corr_caliop_avhrr = np.corrcoef(caliop_max_height_work,
                                            avhrr_height_work)[0,1]
        else:
            corr_caliop_avhrr = -99.0
        diff = avhrr_height_work-caliop_max_height_work
        bias = np.mean(diff)
        diff_squared = diff*diff
        RMS_difference = np.sqrt(np.mean(diff_squared))
        diff_squared_biascorr = (diff-bias)*(diff-bias)
#        RMS_difference_biascorr = np.sqrt(np.mean(diff_squared_biascorr))
    else:
        corr_caliop_avhrr = -9.0
        bias = -9.0
        RMS_difference = -9.0
#        RMS_difference_biascorr = -9.0
        diff_squared_biascorr = np.array([-9.0])
    #return (corr_caliop_avhrr,bias,RMS_difference,avhrr_height_work,diff_squared_biascorr)
    return "%f %f %f %s %f "%(corr_caliop_avhrr,bias,RMS_difference,len(avhrr_height_work),sum(diff_squared_biascorr))

def get_subset_for_mode(caObj, mode):
  # First prepare possible subsetting of CALIOP datasets according to NSIDC
    # and IGBP surface types  
    if mode == 'ICE_COVER_SEA':
##         cal_subset = np.logical_and(np.logical_and(np.less(caObj.calipso.nsidc_surface_type,100),np.greater(caObj.calipso.nsidc_surface_type,10)),np.equal(caObj.calipso.igbp_surface_type,17))
        # Unfortunately, the above formulation used for ORR-B excluded the case
        # when ice-cover was exactly 100 %!!! Very embarrassing!/KG
        cal_subset = np.logical_and(
            np.logical_and(np.less_equal(caObj.calipso.nsidc_surface_type,100),
                           np.greater(caObj.calipso.nsidc_surface_type,10)),
            np.equal(caObj.calipso.igbp_surface_type,17))
    elif mode == 'ICE_FREE_SEA':
        cal_subset = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),np.equal(caObj.calipso.igbp_surface_type,17))
    elif mode == 'SNOW_COVER_LAND':
        cal_subset = np.logical_and(
            np.logical_and(np.less(caObj.calipso.nsidc_surface_type,104),
                           np.greater(caObj.calipso.nsidc_surface_type,10)),
            np.not_equal(caObj.calipso.igbp_surface_type,17))
        # Notice that some uncertainty remains about the meaning of IGBP category 15 = "snow and ice". Can this possibly include also the Arctic ice sheet? We hope that it is not!!! However, if it is, the whole classification here might be wrong since this will affect also the definition of IGBP category 17./KG 
    elif mode == 'SNOW_FREE_LAND':
        cal_subset = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                    np.not_equal(caObj.calipso.igbp_surface_type,17))
    elif mode == 'COASTAL_ZONE':
        cal_subset = np.equal(caObj.calipso.nsidc_surface_type,255)
    
    elif mode == 'TROPIC_ZONE':
        cal_subset = np.abs(caObj.calipso.latitude) <= 10
    elif mode == 'TROPIC_ZONE_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                        np.not_equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.abs(caObj.calipso.latitude) <= 10
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'TROPIC_ZONE_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                        np.equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.abs(caObj.calipso.latitude) <= 10
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    
    elif mode == 'SUB_TROPIC_ZONE':
        cal_subset = np.logical_and((np.abs(caObj.calipso.latitude) > 10), 
                                    (np.abs(caObj.calipso.latitude) <= 45))    
    elif mode == 'SUB_TROPIC_ZONE_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and((np.abs(caObj.calipso.latitude) > 10), 
                                        (np.abs(caObj.calipso.latitude) <= 45))
        cal_subset_area = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                         np.not_equal(caObj.calipso.igbp_surface_type,17))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'SUB_TROPIC_ZONE_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and((np.abs(caObj.calipso.latitude) > 10), 
                                        (np.abs(caObj.calipso.latitude) <= 45))
        cal_subset_area = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                         np.equal(caObj.calipso.igbp_surface_type,17))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    
    elif mode == 'HIGH-LATITUDES':
        cal_subset = np.logical_and((np.abs(caObj.calipso.latitude) > 45), 
                                    (np.abs(caObj.calipso.latitude) <= 75))
    elif mode == 'HIGH-LATITUDES_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                        np.not_equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.logical_and((np.abs(caObj.calipso.latitude) > 45), 
                                         (np.abs(caObj.calipso.latitude) <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'HIGH-LATITUDES_SNOW_COVER_LAND':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less(caObj.calipso.nsidc_surface_type,104),
                           np.greater(caObj.calipso.nsidc_surface_type,10)),
            np.not_equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.logical_and((np.abs(caObj.calipso.latitude) > 45), (np.abs(caObj.calipso.latitude) <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'HIGH-LATITUDES_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                        np.equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.logical_and((np.abs(caObj.calipso.latitude) > 45), 
                                         (np.abs(caObj.calipso.latitude) <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'HIGH-LATITUDES_ICE_COVER_SEA':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less_equal(caObj.calipso.nsidc_surface_type,100),
                           np.greater(caObj.calipso.nsidc_surface_type,10)),
            np.equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.logical_and((np.abs(caObj.calipso.latitude) > 45), 
                                         (np.abs(caObj.calipso.latitude) <= 75))
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    
    
    elif mode == 'POLAR':
        cal_subset = np.abs(caObj.calipso.latitude) > 75
    elif mode == 'POLAR_SNOW_FREE_LAND':
        cal_subset_lat = np.logical_and(np.equal(caObj.calipso.nsidc_surface_type,0),
                                        np.not_equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.abs(caObj.calipso.latitude) > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'POLAR_SNOW_COVER_LAND':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less(caObj.calipso.nsidc_surface_type,104),
                           np.greater(caObj.calipso.nsidc_surface_type,10)),
            np.not_equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.abs(caObj.calipso.latitude) > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'POLAR_ICE_FREE_SEA':
        cal_subset_lat = np.logical_and(
            np.equal(caObj.calipso.nsidc_surface_type,0),np.equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.abs(caObj.calipso.latitude) > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    elif mode == 'POLAR_ICE_COVER_SEA':
        cal_subset_lat = np.logical_and(
            np.logical_and(np.less_equal(caObj.calipso.nsidc_surface_type,100),
                           np.greater(caObj.calipso.nsidc_surface_type,10)),
            np.equal(caObj.calipso.igbp_surface_type,17))
        cal_subset_area = np.abs(caObj.calipso.latitude) > 75
        cal_subset = np.logical_and(cal_subset_lat, cal_subset_area)
    
    elif mode == 'BASIC':
        cal_subset = np.bool_(np.ones(caObj.calipso.igbp_surface_type.shape))
    elif mode == 'OPTICAL_DEPTH':
        cal_subset = np.bool_(np.ones(caObj.calipso.igbp_surface_type.shape))
    elif mode == 'STANDARD':
        cal_subset = np.bool_(np.ones(caObj.calipso.igbp_surface_type.shape))
    elif mode == 'OPTICAL_DEPTH_THIN_IS_CLEAR':
        cal_subset = np.bool_(np.ones(caObj.calipso.igbp_surface_type.shape))
    else:
        print('The mode %s is not added in statistic file' %mode)
        sys.exit()
    return cal_subset     

def get_day_night_info(caObj):
    daynight_flags = None
    if config.CCI_CLOUD_VALIDATION or config.MAIA_CLOUD_VALIDATION:
        daynight_flags = get_day_night_twilight_info_cci2014(
        caObj.avhrr.sunz)
    if config.PPS_VALIDATION and  hasattr(caObj.avhrr, 'cloudtype_qflag'):
        if caObj.avhrr.cloudtype_qflag is not None:
            daynight_flags = get_day_night_twilight_info_pps2012(
                caObj.avhrr.cloudtype_qflag)
    if config.PPS_VALIDATION and  hasattr(caObj.avhrr, 'cloudtype_conditions'):
        if caObj.avhrr.cloudtype_conditions is not None:
            daynight_flags = get_day_night_twilight_info_pps2014(
                caObj.avhrr.cloudtype_conditions)     
    if config.PPS_VALIDATION and daynight_flags == None:
        daynight_flags = get_day_night_twilight_info_cci2014(
        caObj.avhrr.sunz)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag) = daynight_flags
    if (no_qflag.sum() + night_flag.sum() + twilight_flag.sum() + day_flag.sum()) != caObj.calipso.longitude.size:          
        print('something wrong with quality flags. It does not sum up. See beginning of statistic file')
        sys.exit()
    return daynight_flags
    
def get_semi_opaque_info(caObj):
    semi_flag = None    
    opaque_flag = None
    if hasattr(caObj.avhrr, 'cloudtype_qflag'):
        #print caObj.avhrr.ctth_opaque
        if caObj.avhrr.ctth_opaque is not None:
            semi_flag, opaque_flag = get_semi_opaque_info_pps2012(
                caObj.avhrr.ctth_opaque) 
    if hasattr(caObj.avhrr, 'cloudtype_conditions'):
        if caObj.avhrr.ctth_status is not None:
            #print caObj.avhrr.ctth_status
            semi_flag, opaque_flag = get_semi_opaque_info_pps2014(
                caObj.avhrr.ctth_status)
    return semi_flag, opaque_flag


def print_cloudsat_stats(clsatObj, statfile):
    if clsatObj != None:
        cloudsat_cloud_mask=clsatObj.cloudsat.cloud_mask
        cloudsat_cloud_mask=np.greater_equal(cloudsat_cloud_mask, 
                                             config.CLOUDSAT_CLOUDY_THR)
        cloudsat_cloud_fraction=np.zeros(clsatObj.cloudsat.latitude.shape[0])
        sum_cloudsat_cloud_mask = sum(cloudsat_cloud_mask)
        cloudsat_cloud_fraction[sum_cloudsat_cloud_mask > 2] = 1 # requires at least two cloudy bins
        print "Warning what version of pps cloudtype is this!!"        
        cloudsat_clear =  np.less(cloudsat_cloud_fraction,1)
        cloudsat_cloudy = np.greater_equal(cloudsat_cloud_fraction,1)
        pps_clear = np.logical_and(np.less_equal(clsatObj.avhrr.cloudtype,4),
                                   np.greater(clsatObj.avhrr.cloudtype,0))
        pps_cloudy = np.logical_and(np.greater(clsatObj.avhrr.cloudtype,4),
                                    np.less(clsatObj.avhrr.cloudtype,20))

        n_clear_clear = np.repeat(
            pps_clear, np.logical_and(cloudsat_clear,pps_clear)).shape[0]
        n_cloudy_cloudy = np.repeat(
            pps_cloudy,np.logical_and(cloudsat_cloudy,pps_cloudy)).shape[0]
        n_clear_cloudy = np.repeat(
            pps_cloudy,np.logical_and(cloudsat_clear,pps_cloudy)).shape[0]
        n_cloudy_clear = np.repeat(
            pps_clear,np.logical_and(cloudsat_cloudy,pps_clear)).shape[0]
        
        nclear = np.repeat(cloudsat_clear,cloudsat_clear).shape[0]
        ncloudy = np.repeat(cloudsat_cloudy,cloudsat_cloudy).shape[0]
        ncloudy_pps = n_cloudy_cloudy+n_clear_cloudy
        nclear_pps = n_cloudy_clear+n_clear_clear
    
        #print "Number of clear points (CLOUDSAT): ",nclear
        #print "Number of cloudy points (CLOUDSAT): ",ncloudy
        #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
        #print "Cloudsat-Clear PPS-Cloudy, Cloudsat-Cloudy PPS-Clear",n_clear_cloudy,n_cloudy_clear
    
        if ncloudy > 0:
            pod_cloudy = float(n_cloudy_cloudy)/ncloudy
        else:
            pod_cloudy = -9.0
        if ncloudy_pps > 0:
            far_cloudy = float(n_clear_cloudy)/ncloudy_pps
        else:
            far_cloudy = -9.0

        if nclear > 0:
            pod_clear = float(n_clear_clear)/nclear
        else:
            pod_clear = -9.0
        if nclear_pps > 0:
            far_clear = float(n_cloudy_clear)/nclear_pps
        else:
            far_clear = -9.0
    
        mean_cloudsat=((n_clear_clear+n_clear_cloudy)*0.0 + (n_cloudy_clear+n_cloudy_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        mean_pps=((n_clear_clear+n_cloudy_clear)*0.0 + (n_cloudy_cloudy+n_clear_cloudy)*1.0)/(n_clear_clear+n_clear_cloudy+n_cloudy_clear+n_cloudy_cloudy)
        bias=mean_pps-mean_cloudsat
        statfile.write("CLOUD MASK CLOUDSAT-PPS TABLE: %s %s %s %s \n" % (n_clear_clear,n_clear_cloudy,n_cloudy_clear,n_cloudy_cloudy))
        statfile.write("CLOUD MASK CLOUDSAT-PPS PROB: %f %f %f %f %f \n" % (pod_cloudy,pod_clear,far_cloudy,far_clear,bias))
    else:
        statfile.write('No CloudSat \n')
        statfile.write('No CloudSat \n')

def print_cloudsat_modis_stats(clsatObj, statfile):
    if clsatObj != None:
        # CORRELATION CLOUD MASK: CLOUDSAT - MODIS
    
        #print "------------------------------------"
        #print "STATISTICS CLOUD MASK: CLOUDSAT - MODIS" 
    
        modis_clear = np.logical_or(
            np.equal(clsatObj.cloudsat.MODIS_cloud_flag,1),
            np.equal(clsatObj.cloudsat.MODIS_cloud_flag,0))
        modis_cloudy = np.logical_or(
            np.equal(clsatObj.cloudsat.MODIS_cloud_flag,3),
            np.equal(clsatObj.cloudsat.MODIS_cloud_flag,2))
        cloudsat_cloud_mask=clsatObj.cloudsat.cloud_mask
        cloudsat_cloud_mask=np.greater_equal(cloudsat_cloud_mask, 
                                             config.CLOUDSAT_CLOUDY_THR)
        cloudsat_cloud_fraction=np.zeros(clsatObj.cloudsat.latitude.shape[0])
        sum_cloudsat_cloud_mask = sum(cloudsat_cloud_mask)
        cloudsat_cloud_fraction[sum_cloudsat_cloud_mask > 2] = 1 # requires at least two cloudy bins
        print "Warning what version of pps cloudtype is this!!"        
        cloudsat_clear =  np.less(cloudsat_cloud_fraction,1)
        cloudsat_cloudy = np.greater_equal(cloudsat_cloud_fraction,1)

        n_clear_clear = np.repeat(
            modis_clear, np.logical_and(cloudsat_clear,modis_clear)).shape[0]
        n_cloudy_cloudy = np.repeat(
            modis_cloudy, np.logical_and(cloudsat_cloudy,modis_cloudy)).shape[0]
        n_clear_cloudy = np.repeat(
            modis_cloudy, np.logical_and(cloudsat_clear,modis_cloudy)).shape[0]
        n_cloudy_clear = np.repeat(
            modis_clear, np.logical_and(cloudsat_cloudy,modis_clear)).shape[0]
        
        nclear = np.repeat(cloudsat_clear,cloudsat_clear).shape[0]
        ncloudy = np.repeat(cloudsat_cloudy,cloudsat_cloudy).shape[0]
        ncloudy_modis = n_cloudy_cloudy+n_clear_cloudy
        nclear_modis = n_cloudy_clear+n_clear_clear
        
        #print "Number of clear points (CLOUDSAT): ",nclear
        #print "Number of cloudy points (CLOUDSAT): ",ncloudy
        #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
        #print "Cloudsat-Clear MODIS-Cloudy, Cloudsat-Cloudy MODIS-Clear",n_clear_cloudy,n_cloudy_clear
        
        if ncloudy > 0:
            pod_cloudy = float(n_cloudy_cloudy)/ncloudy
        else:
            pod_cloudy = -9.0
        if ncloudy_modis > 0:
            far_cloudy = float(n_clear_cloudy)/ncloudy_modis
        else:
            far_cloudy = -9.0

        if nclear > 0:
            pod_clear = float(n_clear_clear)/nclear
        else:
            pod_clear = -9.0
        if nclear_modis > 0:
            far_clear = float(n_cloudy_clear)/nclear_modis
        else:
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
    else:
        statfile.write('No CloudSat \n')
        statfile.write('No CloudSat \n') 


def print_calipso_cmask_stats(caObj, statfile, cal_subset):
    # CLOUD MASK EVALUATION
    #=======================
    
    # CORRELATION CLOUD MASK: CLOUDSAT - AVHRR

    #print "------------------------------------"
    #print "STATISTICS CLOUD MASK: CLOUDSAT - AVHRR" 
    
       
    # CORRELATION CLOUD MASK: CALIOP - AVHRR

    calipso_clear = np.logical_and(
        np.less(caObj.calipso.cloud_fraction,config.CALIPSO_CLEAR_MAX_CFC),cal_subset)
    calipso_cloudy = np.logical_and(
        np.greater(caObj.calipso.cloud_fraction,config.CALIPSO_CLOUDY_MIN_CFC),cal_subset)
    
        
    # For the combined 1km + 5km dataset cloud_fraction can only have values (0.0, 0.2, 0.4, 0.6, 0.8, 1.0). So the threshold should
    # really be set to 0.4, i.e., at least two 1 km columns should be cloudy!. 
    
    pps_clear = np.logical_and(np.logical_and(np.less_equal(caObj.avhrr.cloudtype,4),np.greater(caObj.avhrr.cloudtype,0)),cal_subset)
    pps_cloudy = np.logical_and(np.logical_and(np.greater(caObj.avhrr.cloudtype,4),np.less(caObj.avhrr.cloudtype,20)),cal_subset)

    #print "------------------------------------"
    #print "STATISTICS CLOUD MASK: CALIOP - AVHRR"
    
    n_clear_clear = np.repeat(
        pps_clear,np.logical_and(calipso_clear,pps_clear)).shape[0]
    n_cloudy_cloudy = np.repeat(
        pps_cloudy,np.logical_and(calipso_cloudy,pps_cloudy)).shape[0]
    n_clear_cloudy = np.repeat(
        pps_cloudy,np.logical_and(calipso_clear,pps_cloudy)).shape[0]
    n_cloudy_clear = np.repeat(
        pps_clear,np.logical_and(calipso_cloudy,pps_clear)).shape[0]
    nclear = np.repeat(calipso_clear,calipso_clear).shape[0]
    ncloudy = np.repeat(calipso_cloudy,calipso_cloudy).shape[0]
    ncloudy_pps = n_cloudy_cloudy+n_clear_cloudy
    nclear_pps = n_cloudy_clear+n_clear_clear
    
    
    #print "Number of clear points (CALIPSO): ",nclear
    #print "Number of cloudy points (CALIPSO): ",ncloudy
    #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
    #print "Calipso-Clear PPS-Cloudy, Calipso-Cloudy PPS-Clear",n_clear_cloudy,n_cloudy_clear
    
    if ncloudy > 0:
        pod_cloudy = float(n_cloudy_cloudy)/ncloudy
    else:
        pod_cloudy = -9.0
    if ncloudy_pps > 0:
        far_cloudy = float(n_clear_cloudy)/ncloudy_pps
    else:
        far_cloudy = -9.0
    if nclear > 0:
        pod_clear = float(n_clear_clear)/nclear
    else:
        pod_clear = -9.0
    if nclear_pps > 0:
        far_clear = float(n_cloudy_clear)/nclear_pps
    else:
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


def print_calipso_modis_stats(caObj, statfile, cal_subset, cal_MODIS_cflag):    
    # CORRELATION CLOUD MASK: CALIOP - MODIS
    calipso_clear = np.logical_and(np.less(caObj.calipso.cloud_fraction,0.34),
                                   cal_subset)
    calipso_cloudy = np.logical_and(
        np.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
    if cal_MODIS_cflag != None:
        modis_clear = np.logical_and(
            np.logical_or(np.equal(cal_MODIS_cflag,1),
                          np.equal(cal_MODIS_cflag,0)),cal_subset)
        modis_cloudy = np.logical_and(
            np.logical_or(np.equal(cal_MODIS_cflag,3),
                          np.equal(cal_MODIS_cflag,2)),cal_subset)

    if cal_MODIS_cflag != None:
        n_clear_clear = np.repeat(
            modis_clear,
            np.logical_and(calipso_clear,modis_clear)).shape[0]
        n_cloudy_cloudy = np.repeat(
            modis_cloudy,
            np.logical_and(calipso_cloudy,modis_cloudy)).shape[0]
        n_clear_cloudy = np.repeat(
            modis_cloudy,
            np.logical_and(calipso_clear,modis_cloudy)).shape[0]
        n_cloudy_clear = np.repeat(
            modis_clear,
            np.logical_and(calipso_cloudy,modis_clear)).shape[0]
        nclear = np.repeat(calipso_clear,calipso_clear).shape[0]
        ncloudy = np.repeat(calipso_cloudy,calipso_cloudy).shape[0]
        ncloudy_modis = n_cloudy_cloudy+n_clear_cloudy
        nclear_modis = n_cloudy_clear+n_clear_clear
            
        #print "Number of clear points (CALIPSO): ",nclear
        #print "Number of cloudy points (CALIPSO): ",ncloudy
        #print "Clear-Clear,Cloudy-Cloudy",n_clear_clear,n_cloudy_cloudy
        #print "Calipso-Clear MODIS-Cloudy, Calipso-Cloudy MODIS-Clear",n_clear_cloudy,n_cloudy_clear
        
        if ncloudy > 0:
            pod_cloudy = float(n_cloudy_cloudy)/ncloudy
        else:
            pod_cloudy = -9.0
        if ncloudy_modis > 0:
            far_cloudy = float(n_clear_cloudy)/ncloudy_modis
        else:
            far_cloudy = -9.0

        if nclear > 0:
            pod_clear = float(n_clear_clear)/nclear
        else:
            pod_clear = -9.0
        if nclear_modis > 0:
            far_clear = float(n_cloudy_clear)/nclear_modis
        else:
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
    else:
        statfile.write("No CloudSat \n")
        statfile.write("No CloudSat \n")
    
def print_calipso_stats_ctype(caObj, statfile, cal_subset, cal_vert_feature):

    # CLOUD TYPE EVALUATION - Based exclusively on CALIPSO data (Vertical Feature Mask)
    # =======================
    calipso_low = np.logical_and(
        np.logical_and(np.greater_equal(cal_vert_feature[::],0),
                       np.less_equal(cal_vert_feature[::],3)),
        cal_subset)
    calipso_medium = np.logical_and(
        np.logical_and(np.greater(cal_vert_feature[::],3),
                       np.less_equal(cal_vert_feature[::],5)),
        cal_subset)
    calipso_high = np.logical_and(
        np.logical_and(np.greater(cal_vert_feature[::],5),
                       np.less_equal(cal_vert_feature[::],7)),
        cal_subset)

    if config.CCI_CLOUD_VALIDATION or caObj.avhrr.cloudtype_conditions is not None:
        if config.CCI_CLOUD_VALIDATION:
            logger.info("CCI validation ctype validation not useful!")
        else:    
            logger.info("Assuming cloudtype structure from pps v2014")
        avhrr_low = np.logical_and(
            np.logical_and(np.greater_equal(caObj.avhrr.cloudtype,5),
                           np.less_equal(caObj.avhrr.cloudtype,6)),
            cal_subset)
        avhrr_medium = np.logical_and(
            np.equal(caObj.avhrr.cloudtype,7), cal_subset)
        avhrr_high_op = np.logical_and(
            np.logical_and(np.greater_equal(caObj.avhrr.cloudtype,8),
                           np.less_equal(caObj.avhrr.cloudtype,9)),
            cal_subset)
        avhrr_high_semi = np.logical_and(
            np.logical_and(np.greater_equal(caObj.avhrr.cloudtype,11),
                           np.less_equal(caObj.avhrr.cloudtype,15)),
            cal_subset)
        avhrr_high = np.logical_or(avhrr_high_op,avhrr_high_semi)
        avhrr_frac = np.logical_and(np.equal(caObj.avhrr.cloudtype,10), 
                                    cal_subset)

    else:
        logger.info("Assuming cloudtype structure from pps v2012")
        avhrr_low = np.logical_and(
            np.logical_and(np.greater_equal(caObj.avhrr.cloudtype,5),
                           np.less_equal(caObj.avhrr.cloudtype,8)),
            cal_subset)
        avhrr_medium = np.logical_and(
            np.logical_and(np.greater_equal(caObj.avhrr.cloudtype,9),
                           np.less_equal(caObj.avhrr.cloudtype,10)),
            cal_subset)
        avhrr_high = np.logical_and(
            np.logical_and(np.greater_equal(caObj.avhrr.cloudtype,11),
                           np.less_equal(caObj.avhrr.cloudtype,18)),
            cal_subset)
        avhrr_frac = np.logical_and(
            np.logical_and(np.greater_equal(caObj.avhrr.cloudtype,19),
                           np.less_equal(caObj.avhrr.cloudtype,19)),
            cal_subset)




    calipso_clear = np.logical_and(
        np.less(caObj.calipso.cloud_fraction,0.34),cal_subset)
    calipso_cloudy = np.logical_and(
        np.greater(caObj.calipso.cloud_fraction,0.66),cal_subset)
    avhrr_clear = np.logical_and(
        np.logical_and(np.less_equal(caObj.avhrr.cloudtype,4),
                       np.greater(caObj.avhrr.cloudtype,0)),
        cal_subset)
    
    
    # Notice that we have unfortunately changed order in notation compared to cloud mask
    # Here the PPS category is mentioned first and then the CALIOP category 

    n_low_low = np.repeat(
        avhrr_low,
        np.logical_and(calipso_low,avhrr_low)).shape[0]
    n_low_medium = np.repeat(
        avhrr_low,
        np.logical_and(calipso_medium,avhrr_low)).shape[0]
    n_low_high = np.repeat(
        avhrr_low,
        np.logical_and(calipso_high,avhrr_low)).shape[0]
    n_medium_low = np.repeat(
        avhrr_medium,
        np.logical_and(calipso_low,avhrr_medium)).shape[0]
    n_medium_medium = np.repeat(
        avhrr_medium,
        np.logical_and(calipso_medium,avhrr_medium)).shape[0]
    n_medium_high = np.repeat(
        avhrr_medium,
        np.logical_and(calipso_high,avhrr_medium)).shape[0]
    n_high_low = np.repeat(
        avhrr_high, 
        np.logical_and(calipso_low,avhrr_high)).shape[0]
    n_high_medium = np.repeat(
        avhrr_high,
        np.logical_and(calipso_medium,avhrr_high)).shape[0]
    n_high_high = np.repeat(
        avhrr_high,
        np.logical_and(calipso_high,avhrr_high)).shape[0]
    n_frac_low = np.repeat(
        avhrr_frac,
        np.logical_and(calipso_low,avhrr_frac)).shape[0]
    n_frac_medium = np.repeat(
        avhrr_frac,
        np.logical_and(calipso_medium,avhrr_frac)).shape[0]
    n_frac_high = np.repeat(
        avhrr_frac,
        np.logical_and(calipso_high,avhrr_frac)).shape[0]

#    nlow = np.repeat(calipso_low,calipso_low).shape[0]
#    nmedium = np.repeat(calipso_medium,calipso_medium).shape[0]
#    nhigh = np.repeat(calipso_high,calipso_high).shape[0]
        
    n_clear_low = np.repeat(
        avhrr_clear,
        np.logical_and(calipso_low,avhrr_clear)).shape[0]
    n_clear_medium = np.repeat(
        avhrr_clear,
        np.logical_and(calipso_medium,avhrr_clear)).shape[0]
    n_clear_high = np.repeat(
        avhrr_clear,
        np.logical_and(calipso_high,avhrr_clear)).shape[0]
    n_low_clear = np.repeat(
        avhrr_low,
        np.logical_and(calipso_clear,avhrr_low)).shape[0]
    n_medium_clear = np.repeat(
        avhrr_medium,
        np.logical_and(calipso_clear,avhrr_medium)).shape[0]
    n_high_clear = np.repeat(
        avhrr_high,
        np.logical_and(calipso_clear,avhrr_high)).shape[0]
    n_frac_clear = np.repeat(
        avhrr_frac,
        np.logical_and(calipso_clear,avhrr_frac)).shape[0]
    
    #print "Number of points with low clouds (CALIPSO): ",nlow
    #print "Number of points with medium clouds (CALIPSO): ",nmedium
    #print "Number of points with high clouds (CALIPSO): ",nhigh
    #print "CALIPSO low-AVHRR low,CALIPSO low-AVHRR medium, CALIPSO low-AVHRR high",n_low_low,n_medium_low,n_high_low
    #print "CALIPSO medium-AVHRR low,CALIPSO medium-AVHRR medium, CALIPSO medium-AVHRR high",n_low_medium,n_medium_medium,n_high_medium
    #print "CALIPSO high-AVHRR low,CALIPSO high-AVHRR medium, CALIPSO high-AVHRR high",n_low_high,n_medium_high,n_high_high
    #print "CALIPSO high-AVHRR frac,CALIPSO medium-AVHRR frac, CALIPSO low-AVHRR frac",n_frac_high,n_frac_medium,n_frac_low
    
    if (n_low_low+n_medium_low+n_high_low+n_frac_low) > 0:
        pod_low = float(n_low_low + n_frac_low)/(n_low_low+n_medium_low+n_high_low+n_frac_low)
        far_low = float(n_medium_low+n_high_low)/(n_low_low+n_medium_low+n_high_low+n_frac_low)
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
            

def print_cloudsat_stats_ctop(clsatObj, statfile):

    # CLOUD TOP EVALUATION
    #=======================
    if clsatObj != None: 
        avhrr_ctth_csat_ok = clsatObj.avhrr.avhrr_ctth_csat_ok
        # CORRELATION: CLOUDSAT - AVHRR HEIGHT
    
        #print "STATISTICS CLOUD TOP HEIGHT: CLOUDSAT - AVHRR"
        clsat_max_height = -9 + 0*np.zeros(clsatObj.cloudsat.latitude.shape)
    
        for i in range(125):
            height = clsatObj.cloudsat.Height[:,i]
            cmask_ok = clsatObj.cloudsat.cloud_mask[:,i]
            top_height = height+120
            #top_height[height<240*4] = -9999 #Do not use not sure why these are not used Nina 20170317
            is_cloudy = cmask_ok > config.CLOUDSAT_CLOUDY_THR
            top_height[~is_cloudy] = -9999
            clsat_max_height[clsat_max_height<top_height] =  top_height[clsat_max_height<top_height] 

        okarr = np.greater_equal(avhrr_ctth_csat_ok,0.0)
        okarr = np.logical_and(okarr,np.greater_equal(clsat_max_height,0.0))
        clsat_max_height = np.repeat(clsat_max_height[::],okarr)
        avhrr_height = np.repeat(avhrr_ctth_csat_ok[::],okarr)
    
        if len(avhrr_height) > 0:
            if len(avhrr_height) > 20:
                corr_cloudsat_avhrr = np.corrcoef(clsat_max_height,avhrr_height)[0,1]
            else:
                corr_cloudsat_avhrr = -99.0
            #print "Correlation: cloudsat-avhrr = ",corr_cloudsat_avhrr
            diff = avhrr_height-clsat_max_height
            bias = np.mean(diff)
            #print "Mean difference cloudsat-avhrr = ",bias
            diff_squared = diff*diff
            RMS_difference = np.sqrt(np.mean(diff_squared))
            #print "RMS difference cloudsat-avhrr = ",RMS_difference
            diff_squared = (diff-bias)*(diff-bias)
            RMS_difference = np.sqrt(np.mean(diff_squared))
            #print "Bias-corrected RMS difference cloudsat-avhrr = ",RMS_difference
            #print "Number of matchups: ", len(avhrr_height)
            #print "Number of failing PPS CTTHs: ", len(np.repeat(clsatObj.avhrr.cloudtype[::],np.greater(clsatObj.avhrr.cloudtype,4)))-len(avhrr_height)
            #print
    
            statfile.write("CLOUD HEIGHT CLOUDSAT: %f %f %f %s %f \n" % (corr_cloudsat_avhrr,bias,RMS_difference,len(avhrr_height),sum(diff_squared)))
        else:
            statfile.write("CLOUD HEIGHT CLOUDSAT: -9.0 -9.0 -9.0 0 -9.0 \n")
    else:
        statfile.write('No CloudSat \n')


def print_height_all_low_medium_high(NAME, okcaliop,  statfile, 
                                     cal_vert_feature, avhrr_ctth_cal_ok,
                                     caliop_max_height):
    out_stats = calculate_ctth_stats(okcaliop,avhrr_ctth_cal_ok,
                                     caliop_max_height)
    statfile.write("CLOUD HEIGHT %s ALL: %s\n" %(NAME, out_stats))
#    statfile.write("CLOUD HEIGHT CALIOP ALL: 

    # THEN FOR LOW CLOUDS (VERTICAL FEATURE MASK CATEGORIES 0-3)
    cal_low_ok = np.logical_and(np.greater_equal(cal_vert_feature[::],0),
                                np.less_equal(cal_vert_feature[::],3))
    cal_low_ok = np.logical_and(cal_low_ok, okcaliop) 
    out_stats =calculate_ctth_stats(cal_low_ok,avhrr_ctth_cal_ok,
                                    caliop_max_height)   
    statfile.write("CLOUD HEIGHT %s LOW: %s \n" % (NAME, out_stats))
        
    # THEN FOR MEDIUM CLOUDS (VERTICAL FEATURE MASK CATEGORIES 4-5)
    cal_mid_ok = np.logical_and(np.greater(cal_vert_feature[::],3),
                                np.less_equal(cal_vert_feature[::],5))
    cal_mid_ok = np.logical_and(cal_mid_ok,okcaliop)
    out_stats = calculate_ctth_stats(cal_mid_ok, avhrr_ctth_cal_ok,
                                     caliop_max_height)   
    statfile.write("CLOUD HEIGHT %s MEDIUM: %s \n" % (NAME, out_stats))
            
    # FINALLY FOR HIGH CLOUDS (VERTICAL FEATURE MASK CATEGORIES 6-7)
    cal_high_ok = np.logical_and(np.greater(cal_vert_feature[::],5),
                                 np.less_equal(cal_vert_feature[::],7))
    cal_high_ok = np.logical_and(cal_high_ok,okcaliop)
    out_stats = calculate_ctth_stats(cal_high_ok, avhrr_ctth_cal_ok,
                                     caliop_max_height) 
    statfile.write("CLOUD HEIGHT %s HIGH: %s \n" % (NAME, out_stats))

def print_calipso_stats_ctop(caObj, statfile, cal_subset, cal_vert_feature, 
                             semi_flag, opaque_flag):

    # CORRELATION: CALIOP - AVHRR HEIGHT
    # FIRST TOTAL FIGURES
    avhrr_ctth_cal_ok = caObj.avhrr.avhrr_ctth_cal_ok  
    cal_data_ok = np.greater(caObj.calipso.layer_top_altitude[:,0],-9.)
    caliop_max_height = caObj.calipso.layer_top_altitude[:,0].copy()
    caliop_max_height[caliop_max_height>=0] *= 1000
    caliop_max_height[caliop_max_height<0] = -9

    okcaliop = np.logical_and(
        cal_data_ok, 
        np.logical_and(np.greater_equal(avhrr_ctth_cal_ok[::],0),cal_subset))
    #print "ALL CLOUDS:" 
    print_height_all_low_medium_high("CALIOP", 
                                     okcaliop, 
                                     statfile, cal_vert_feature, 
                                     avhrr_ctth_cal_ok, caliop_max_height)
    

    if config.COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC:
        statfile.write("CLOUD HEIGHT SINGLE-LAYER\n")
        okcaliop_single = np.logical_and(
            okcaliop, 
            np.equal(caObj.calipso.number_layers_found,1)) 
        print_height_all_low_medium_high("CALIOP-SINGLE-LAYER", 
                                         okcaliop_single,
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)

    if (config.COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC and
        (config.ALSO_USE_5KM_FILES or config.RESOLUTION==5)): 
        statfile.write("CLOUD HEIGHT SINGLE-LAYER, NOT THIN\n")
        lim=2*config.OPTICAL_DETECTION_LIMIT
        okcaliop_single_not_thinnest = np.logical_and(
            okcaliop_single, 
            np.greater_equal(caObj.calipso.feature_optical_depth_532_top_layer_5km,lim))
        print_height_all_low_medium_high("CALIOP-SINGLE-LAYER>%f"%(lim), 
                                         okcaliop_single_not_thinnest,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)
        
        statfile.write("CLOUD HEIGHT NOT VERY THIN TOP LAYER\n")
        lim=config.OPTICAL_DETECTION_LIMIT
        okcaliop_not_thinnest_top_layer = np.logical_and(
            okcaliop, 
            np.greater_equal(caObj.calipso.feature_optical_depth_532_top_layer_5km,lim))
        print_height_all_low_medium_high("CALIOP-TOP-LAYER>%f"%(lim), 
                                         okcaliop_not_thinnest_top_layer,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)

        lim=config.OPTICAL_DETECTION_LIMIT
        statfile.write("CLOUD HEIGHT VERY THIN TOP LAYER\n")
        okcaliop_thinnest_top_layer = np.logical_and(
            okcaliop, 
            np.less_equal(caObj.calipso.feature_optical_depth_532_top_layer_5km,lim))
        print_height_all_low_medium_high("CALIOP-TOP-LAYER<=%f"%(lim), 
                                         okcaliop_thinnest_top_layer,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)
           
    if (config.COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE and 
        semi_flag is not None or
        opaque_flag is not None):
            
        #Opaque stats
        statfile.write("CLOUD HEIGHT OPAQUE\n")
        okcaliop_opaque_pps = np.logical_and(okcaliop, opaque_flag)
        print_height_all_low_medium_high("CALIOP-OPAQUE", 
                                         okcaliop_opaque_pps,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)

        #Semi-transparent stats
        statfile.write("CLOUD HEIGHT SEMI\n")
        okcaliop_semi_pps = np.logical_and(okcaliop, semi_flag)
        print_height_all_low_medium_high("CALIOP-SEMI", 
                                         okcaliop_semi_pps,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)


    if (config.COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC and
        (config.ALSO_USE_5KM_FILES or config.RESOLUTION==5) and
        config.COMPILE_RESULTS_SEPARATELY_FOR_SEMI_AND_OPAQUE and 
        semi_flag is not None or
        opaque_flag is not None):
            
        #Thin top layer
        #Opaque stats
        lim=config.OPTICAL_DETECTION_LIMIT
        statfile.write("CLOUD HEIGHT OPAQUE TOP LAYER VERY NOT-THIN\n")
        okcaliop_opaque_pps_not_thin = np.logical_and(
            np.logical_and(okcaliop, opaque_flag),
            np.greater(caObj.calipso.feature_optical_depth_532_top_layer_5km,lim))
        print_height_all_low_medium_high("CALIOP-OPAQUE-TOP-LAYER>%f"%(lim), 
                                         okcaliop_opaque_pps_not_thin,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)
        #Semi-transparent stats
        statfile.write("CLOUD HEIGHT SEMI TOP LAYER VERY NOT-THIN\n")
        okcaliop_semi_pps_not_thin = np.logical_and(
            np.logical_and(okcaliop, semi_flag),
            np.greater(caObj.calipso.feature_optical_depth_532_top_layer_5km,lim))
        print_height_all_low_medium_high("CALIOP-SEMI-TOP-LAYER>%f"%(lim), 
                                         okcaliop_semi_pps_not_thin,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)
        #Not thin top layer
        #Opaque stats
        lim=config.OPTICAL_DETECTION_LIMIT
        statfile.write("CLOUD HEIGHT OPAQUE TOP LAYER VERY THIN\n")
        okcaliop_opaque_pps_thin = np.logical_and(
            np.logical_and(okcaliop, opaque_flag),
            np.less_equal(caObj.calipso.feature_optical_depth_532_top_layer_5km,lim))
        print_height_all_low_medium_high("CALIOP-OPAQUE-TOP-LAYER<=%f"%(lim), 
                                         okcaliop_opaque_pps_thin,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)
        #Semi-transparent stats
        statfile.write("CLOUD HEIGHT SEMI TOP LAYER VERY THIN\n")
        okcaliop_semi_pps_thin = np.logical_and(
            np.logical_and(okcaliop, semi_flag),
            np.less_equal(caObj.calipso.feature_optical_depth_532_top_layer_5km,lim))
        print_height_all_low_medium_high("CALIOP-SEMI-TOP-LAYER<=%f"%(lim), 
                                         okcaliop_semi_pps_thin,  
                                         statfile, cal_vert_feature, 
                                         avhrr_ctth_cal_ok, caliop_max_height)

def CalculateStatistics(mode, clsatObj, statfile, caObj, cal_MODIS_cflag,
                        dnt_flag = None): 
    import sys
 
    # Extract CALIOP Vertical Feature Classification (bits 10-12) from 16 bit
    # representation for topmost cloud layer #Nina 20140120 this is cloud type nit vert feature !!
    cal_vert_feature = np.ones(caObj.calipso.layer_top_altitude[::,0].shape)*-9
    feature_array = 4*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],11),1) + 2*np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],10),1) + np.bitwise_and(np.right_shift(caObj.calipso.feature_classification_flags[::,0],9),1)
    cal_vert_feature = np.where(np.not_equal(caObj.calipso.feature_classification_flags[::,0],1),feature_array[::],cal_vert_feature[::])  



    cal_subset = get_subset_for_mode(caObj, mode)
    semi_flag, opaque_flag = get_semi_opaque_info(caObj)
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag) = get_day_night_info(caObj)

    if dnt_flag == None:
        print('dnt_flag = %s' %'NO DNT FLAG -> ALL PIXELS')
        cal_subset = np.logical_and(cal_subset, all_dnt_flag)
    elif dnt_flag.upper() == 'DAY':
        print('dnt_flag = %s' %dnt_flag.upper())
        cal_subset = np.logical_and(cal_subset, day_flag)
    elif dnt_flag.upper() == 'NIGHT':
        print('dnt_flag = %s' %dnt_flag.upper())
        cal_subset = np.logical_and(cal_subset, night_flag)
    elif dnt_flag.upper() == 'TWILIGHT':
        print('dnt_flag = %s' %dnt_flag.upper())
        cal_subset = np.logical_and(cal_subset, twilight_flag)
    else:
        print('dnt_flag = %s' %dnt_flag.upper())
        print('statistic calculation is not prepared for this dnt_flag')
        sys.exit()

    print_cloudsat_stats(clsatObj, statfile)
    print_cloudsat_modis_stats(clsatObj, statfile)
    print_calipso_cmask_stats(caObj, statfile, cal_subset)
    print_calipso_modis_stats(caObj, statfile, cal_subset, cal_MODIS_cflag)
    print_calipso_stats_ctype(caObj, statfile, cal_subset, cal_vert_feature)
    print_cloudsat_stats_ctop(clsatObj, statfile)
    print_calipso_stats_ctop(caObj, statfile, cal_subset, cal_vert_feature, 
                             semi_flag, opaque_flag)
   


    statfile.close()
    
