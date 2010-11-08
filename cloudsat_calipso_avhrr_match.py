
# Program cloudsat_calipso_avhrr_match.py

# This program is used to process and output statistics for the inter-comparison of AVHRR PPS
# results and CloudSat/CALIPSO observations. It may be run
# repeatedly and supervised by program cloudsat_calipso_process_master.py.

# FOR ERIK: Let's wait with this and only run initially on individual cases. Start the program
# first with the following command: python cloudsat_calipso_avhrr_match.py.
# You will then get a help statement on usage telling you what input parameters are required.

# This particular version of Adam's original CloudSat/CALIPSO matchup and analysis program has been complemented
# with the following:

#  * A method to calculate cloud emissivities for the uppermost CALIPSO cloud layer. With the
#    use of parameters EMISS_FILTERING, EMISS_MIN_HEIGHT and EMISS_LIMIT the thinnest uppermost
#    CALIPSO cloud layers can be analysed and the entire column can be disregarded if the
#    emissivity falls below the EMISS_LIMIT value.
#    Cloud emissivities Ec are calculated as follows:
#
#                     Ec = (I-Iclear)/(B(Tc)-Iclear)
#    where
#        I = Measured radiance in AVHRR channel 4 (11 micron)
#            To be calculated as the Planck radiance for the associated brightness temperature
#        Iclear = Estimated radiance in cloud free situations
#                 To be calculate as the Planck radiance for the NWP-analysed surface temperature
#                 (i.e., neglecting further atmospheric contributions)
#        B(Tc) = Planck radiance for the uppermost cloud layer using CALIPSO mid-layer temperatures

#  * Adjusted scales between CloudSat and CALIPSO datasets. The previous assumption that both datasets
#    had 1 km resolution resulted in that datasets went out of phase for distances longer than about
#    1000 km. An empirical scale factor (CLOUDSAT_TRACK_RESOLUTION) of 1.076 is used to get the most
#    optimal match.

#  * AVHRR cloud top height datasets have been recalculated to heights above mean sea level using
#    CloudSat and CALIPSO elevation data

#  * The MODIS cloud flag has been added to the extracted CALIPSO dataset. This enables direct
#    comparisons to the MODIS cloud mask! Consequently, corresponding MODIS Cloud Mask statistics
#    are calculated and printed.

#  * The Vertical Feature Mask parameter in the CALIPSO dataset has been used to subdivide results
#    into three cloud groups: Low, Medium and High. This has enabled an evaluation of PPS Cloud Type
#    results and a further sub-division of Cloud Top Height results

#  * The National Snow and Ice Data Center (NSIDC) ice and snow mapping results have been added
#    to the extracted Calipso parameters. Together with the IGBP land use classification it is then
#    possible to isolate the study to focus on one of the following categories:

#        ICE_COVER_SEA, ICE_FREE_SEA, SNOW_COVER_LAND, SNOW_FREE_LAND or COASTAL_ZONE

#  RUNNING INSTRUCTIONS:

# The program is capable of running in a wide range of modes (according to description above). These
# various modes are selected by enabling (disabling) the following parameters:

# PLOT_OPTION, EMISS_FILTERING, ICE_COVER_SEA, ICE_FREE_SEA, SNOW_COVER_LAND, SNOW_FREE_LAND, COASTAL_ZONE

# However, notice that only one mode can be chosen for each run. The only exception is PLOT_OPTION
# (i.e., the generation of a PNG plot) which can be combined with EMISS_FILTERING. This also means
# that processing of individual surface categories only generates statistics and not any plots.

# Input data has to be supplied at directories defined by MAIN_DIR and SUB_DIR parameters below.
 
# Every exection of the program prints out statistics for PPS Cloud Mask, Cloud Type and Cloud Top Height
# directly on the screen.

# Don't forget to set all parameters needed for standard ACPG/AHAMAP execution since the matchup
# software uses parts of the ACPG/AHAMAP software.

# FOR ERIK: This can be done by the following command:
#  source /data/proj/saf/kgkarl/GAC_PPS/new_runscripts/source_me_cmsaf

# Dependencies: For a successful run of the program the following supporting python modules must be
#               available in the default run directory:

#               cloudsat.py
#               calipso.py
#               calipso_avhrr_matchup.py
#               cloudsat_avhrr_matchup.py
#               radiance_tb_tables_kgtest.py

#               For full consistency make sure that MAIN_DIR and SUB_DIR parameters are the same also
#               in modules cloudsat.py and calipso.py.

# Output data files: Main results are generally written in the directory MAIN_DIR/SUB_DIR but plotting
#                    results are stored at ./Plot and temporary results at ./Data directories. Thus,
#                    make sure that these directories exist as subdirectories at the default run directory.
#                   This is now made automatic /Erik

# Finally, notice that the matching of the PPS, CloudSat and CALIPSO datasets have been calculated using the
# fix area arctic_super_5010 defined over the Arctic region in Lambert Azimuthal Equal Area projection. Thus,
# for matching data to other regions please modify modules calipso.py and cloudsat.py and replace area
# arctic_super_5010 with the desired area.
# This is now made below /Erik

# /KG March 2010

import os, string
import sys
import numpy.oldnumeric as Numeric
import rpy
import numpy
import inspect
from radiance_tb_tables_kgtest import * #Just use the brightness temperature to radiance conversion/KG

from pps_basic_configure import *
from pps_error_messages import *

import config

from cloudsat_calipso_avhrr_statistics import *
from cloudsat_calipso_avhrr_plot import *
from trajectory_plot import *
from cloudsat_calipso_avhrr_prepare import *

from cloudsat_avhrr_matchup import *
from calipso_avhrr_matchup import *
from cloudsat_avhrr_matchup5km import *
from calipso_avhrr_matchup5km import *

# -----------------------------------------------------

def run(cloudsatfile, calipsofile, ctypefile, ctthfile, avhrrfile, surftfile, sunanglefile, process_mode, Resolution):
    write_log('INFO', "Case: %s" % ctypefile)
    write_log('INFO', "Process mode: %s" % process_mode)
    
    cloudsat_type = os.path.basename(cloudsatfile[0]).split(".h5")[0]
    cloudsat_type = string.join(cloudsat_type.split("_")[-5].split("-")[1:],"-")
        
    basename = os.path.basename(ctypefile).split(".h5")[0]  # delar vid h5 
    base_sat = basename.split("_")[-8]
    base_year = basename.split("_")[-7][0:4]
    base_month = basename.split("_")[-7][4:6]
    basename = string.join(basename.split("_")[0:4],"_") 
    
    #cloudsattype = os.path.basename(cloudsatfile[0]).split(".h5")[0]
    #cloudsattype = cloudsattype.split("_")[4].split("-")[1]
    
    sl = string.split(basename,"_")
    platform = sl[0]
    #print "platform: ", platform
    if platform == "noaa17":
        noaa_number=17
    elif platform == "noaa18":
        noaa_number=18
    elif platform == "metop02":
        noaa_number=2 # Poor man's solution!
    elif platform == 'noaa19':
        noaa_number = 19
    else:
        raise NotImplementedError("Support for satellite %s is not yet implemented." % platform)
        
    norbit = string.atoi(sl[3]) #verkar inte vara string langre?????
    yyyymmdd = sl[1]


    # Now fetch all the datasets for the section of the AREA where all
    # three datasets match. Also get maximum and minimum time differences to AVHRR (in seconds)

    if int(Resolution[0])==1:
        AREA=AREA1KM
 
#        if int(Resolution[2])==1:
#            cloudsat_type = 'GEOPROF'
#        elif int(Resolution[2])==2:
#            cloudsat_type = 'CWC-RVOD'
#        else:
#            write_log("INFO","Define Cloudsat type")
#            print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
#            sys.exit(-9)  
        
        clsatObj,clsat_min_diff,clsat_max_diff = getCloudsatAvhrrMatch(avhrrfile,cloudsatfile,ctypefile,ctthfile,surftfile,sunanglefile,cloudsat_type)
        caObj,ca_min_diff,ca_max_diff = getCaliopAvhrrMatch(avhrrfile,calipsofile,ctypefile,ctthfile,surftfile,sunanglefile)
        clsatObj, caObj = CloudsatCalipsoAvhrrSatz1km(clsatObj,caObj)
        if cloudsat_type=='GEOPROF':
            cllon = clsatObj.cloudsat.longitude.copy()
            cllat = clsatObj.cloudsat.latitude.copy()
        elif cloudsat_type=='CWC-RVOD':
            cllon = clsatObj.cloudsatcwc.longitude.copy()
            cllat = clsatObj.cloudsatcwc.latitude.copy()
        calon = caObj.calipso.longitude.copy()
        calat = caObj.calipso.latitude.copy()
        avhrlon = caObj.avhrr.longitude.copy()
        avhrlat = caObj.avhrr.latitude.copy()
    elif int(Resolution[0])==5:
        AREA=AREA5KM
                
        clsatObj,clsat_min_diff,clsat_max_diff = getCloudsat5kmAvhrr5kmMatch(avhrrfile,cloudsatfile,ctypefile,ctthfile,surftfile,sunanglefile,cloudsat_type)
        caObj,ca_min_diff,ca_max_diff = getCaliop5kmAvhrr5kmMatch(avhrrfile,calipsofile,ctypefile,ctthfile,surftfile,sunanglefile)           
        #caObj,ca_min_diff,ca_max_diff = getCaliop5kmAvhrr5kmMatch(avhrrfile,calipsofile1,calipsofile2,calipsofile3,ctypefile,ctthfile,surftfile,int(Resolution[-1]))
        if cloudsat_type=='GEOPROF':
            cllon = clsatObj.cloudsat5km.longitude.copy()
            cllat = clsatObj.cloudsat5km.latitude.copy()
        elif cloudsat_type=='CWC-RVOD':
            cllon = clsatObj.cloudsat5kmcwc.longitude.copy()
            cllat = clsatObj.cloudsat5kmcwc.latitude.copy()
        calon = caObj.calipso5km.longitude.copy()
        calat = caObj.calipso5km.latitude.copy()
        avhrlon = caObj.avhrr5km.longitude.copy()
        avhrlat = caObj.avhrr5km.latitude.copy()
    else:
        write_log("INFO","Define resolution")
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9) 
    
    ##########################################################################################################################################   
    #pdb.set_trace()
    # Issue a warning if startpoint or endpoint latitude of CloudSat and CALIPSO differ by more than 0.1 degrees
    # Furthermore, if startpoint differs it means that CALIPSO data is not available for the first part of the matchup
    # cross section. This means that we must find the first corresponding CloudSat point to this CALIPSO start point
    # before we can do the plotting (statistics calculations are not affected). Consequently, set the CALIPSO_DISPLACED
    # flag and find correct startpoint just before starting the plotting of CALIPSO data!

    ### 1 KM DATA GEOPROF ###
    if int(Resolution[0])==1 and cloudsat_type=='GEOPROF':   
        CALIPSO_DISPLACED = 0
        latdiff = abs(clsatObj.cloudsat.latitude[0] - caObj.calipso.latitude[0])
        print "latdiff: ", latdiff
        if latdiff > 0.1:
            write_log('INFO', "CloudSat/CALIPSO startpoint differ by %f degrees." % latdiff)
            write_log('INFO', "Cloudsat start lon, lat: %f, %f" % \
                      (clsatObj.cloudsat.longitude[0], clsatObj.cloudsat.latitude[0]))
            write_log('INFO', "CALIPSO start lon, lat: %f, %f" % \
                      (caObj.calipso.longitude[0], caObj.calipso.latitude[0]))
            CALIPSO_DISPLACED = 1
            for j in range(clsatObj.cloudsat.latitude.shape[0]):
                if (abs(clsatObj.cloudsat.latitude[j] - caObj.calipso.latitude[0]) < 0.05) and\
                       (abs(clsatObj.cloudsat.longitude[j] - caObj.calipso.longitude[0]) < 0.1):
                    calipso_displacement=int(j*CLOUDSAT_TRACK_RESOLUTION)
                    write_log('INFO', "CALIPSO_DISPLACEMENT: %d" % calipso_displacement)
                    break
        # First make sure that PPS cloud top heights are converted to height above sea level
        # just as CloudSat and CALIPSO heights are defined. Use corresponding DEM data.            
        elevation = Numeric.where(Numeric.greater(clsatObj.cloudsat.elevation,0),
                            clsatObj.cloudsat.elevation,-9)			# If clsatObj.cloudsat.elevation is bigger then 0 elevation(i,j)=clsatObj.cloudsat.elevation(i,j) else the value =-9

###        elevationcwc = Numeric.where(Numeric.greater(clsatObj.cloudsatcwc.elevation,0),
###                            clsatObj.cloudsatcwc.elevation,-9)

###        data_okcwc = Numeric.ones(clsatObj.cloudsatcwc.elevation.shape,'b')

        data_ok = Numeric.ones(clsatObj.cloudsat.elevation.shape,'b')
        print "Length of CLOUDSAT array: ", len(data_ok)
        lat_ok = Numeric.repeat(clsatObj.cloudsat.latitude[::],data_ok)
        avhrr_ctth_csat_ok = Numeric.repeat(clsatObj.avhrr.ctth_height[::],data_ok)
        avhrr_ctth_csat_ok = Numeric.where(Numeric.greater(avhrr_ctth_csat_ok,0.0),avhrr_ctth_csat_ok[::]+elevation*1.0,avhrr_ctth_csat_ok)

        # CALIPSO

        cal_elevation = Numeric.where(Numeric.greater(caObj.calipso.elevation,0),
                                    caObj.calipso.elevation,-9)
        cal_data_ok = Numeric.ones(caObj.calipso.elevation.shape,'b')
        print "Length of CALIOP array: ", len(cal_data_ok)
        avhrr_ctth_cal_ok = Numeric.repeat(caObj.avhrr.ctth_height[::],cal_data_ok)
        avhrr_ctth_cal_ok = Numeric.where(Numeric.greater(avhrr_ctth_cal_ok,0.0),avhrr_ctth_cal_ok[::]+cal_elevation,avhrr_ctth_cal_ok)                    
        
        if (len(cal_data_ok) == 0) or (len(data_ok) == 0):
            print "Processing stopped: Zero lenght of matching arrays!"
            print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit()
        else:
            # If everything is OK, now create filename for statistics output file and open it for writing.
            # Notice that more than one file
            # (but maximum 2) can be created for one particular noaa orbit.
            
            resultpath = "%s/%s/%ikm/%s/%s/%s/%s" % (setup.RESULT_DIR, base_sat,int(Resolution[0]), base_year, base_month,AREA,process_mode)
            if not os.path.exists(resultpath):
                os.makedirs(resultpath)
            statfilename = "%s/%ikm_%s_cloudsat_calipso_avhrr_stat.dat" % (resultpath,int(Resolution[0]),basename)

            statfile = open(statfilename,"w")
            if process_mode == "BASIC":
                statfile.write("CloudSat min and max time diff: %f %f \n" %(clsat_min_diff,clsat_max_diff))
                statfile.write("CALIPSO min and max time diff: %f %f \n" %(ca_min_diff,ca_max_diff))
            else:
                statfile.write("CloudSat min and max time diff: See results for BASIC! \n")
                statfile.write("CALIPSO min and max time diff: See results for BASIC! \n")
            statfile.write("Start-Stop-Length Cloudsat: %f %f %f %f %s \n" %(clsatObj.cloudsat.latitude[0],clsatObj.cloudsat.longitude[0],clsatObj.cloudsat.latitude[len(clsatObj.cloudsat.latitude)-1],clsatObj.cloudsat.longitude[len(clsatObj.cloudsat.latitude)-1],len(data_ok)))
            statfile.write("Start-Stop-Length CALIPSO: %f %f %f %f %s \n" %(caObj.calipso.latitude[0],caObj.calipso.longitude[0],caObj.calipso.latitude[len(caObj.calipso.latitude)-1],caObj.calipso.longitude[len(caObj.calipso.latitude)-1],len(cal_data_ok)))
        # Extract CALIOP Vertical Feature Classification (bits 10-12) from 16 bit representation
        # for topmost cloud layer
        cal_vert_feature = Numeric.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9
        feature_array = 4*Numeric.bitwise_and(Numeric.right_shift(caObj.calipso.feature_classification_flags[0,::],11),1) + 2*Numeric.bitwise_and(Numeric.right_shift(caObj.calipso.feature_classification_flags[0,::],10),1) + Numeric.bitwise_and(Numeric.right_shift(caObj.calipso.feature_classification_flags[0,::],9),1)

        cal_vert_feature = Numeric.where(Numeric.not_equal(caObj.calipso.feature_classification_flags[0,::],1),feature_array[::],cal_vert_feature[::])   
        # Prepare for plotting, cloud emissivity and statistics calculations
        
        midlayer_temp_kelvin = caObj.calipso.cloud_mid_temperature + 273.15
        caliop_toplay_thickness = Numeric.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9

        caliop_max_height_midlaytemp = Numeric.where(Numeric.greater(caObj.calipso.cloud_top_profile[0,::],-9),midlayer_temp_kelvin[0,::],-9)
        thickness = (caObj.calipso.cloud_top_profile[0,::]-caObj.calipso.cloud_base_profile[0,::])*1000.
        caliop_toplay_thickness = Numeric.where(Numeric.greater(caObj.calipso.cloud_top_profile[0,::],-9),thickness,caliop_toplay_thickness)

        caliop_height = []
        caliop_base = []
        caliop_max_height = Numeric.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9
        
        for i in range(10):
            hh = Numeric.where(Numeric.greater(caObj.calipso.cloud_top_profile[i,::],-9),
                                caObj.calipso.cloud_top_profile[i,::] * 1000.,-9)
                                            
            caliop_max_height = Numeric.maximum(caliop_max_height,
                                                caObj.calipso.cloud_top_profile[i,::] * 1000.)
            # This is actually unnecessary - we know that layer 1 is always the highest layer!!
            # However, arrays caliop_height and caliop_base are needed later for plotting/ KG

            caliop_height.append(hh)
            bb = Numeric.where(Numeric.greater(caObj.calipso.cloud_base_profile[i,::],-9),
                                caObj.calipso.cloud_base_profile[i,::] * 1000.,-9)
            caliop_base.append(bb)
            thickness = hh - bb
            
        x = Numeric.repeat(caObj.calipso.number_of_layers_found.ravel(),
                            Numeric.greater(caObj.calipso.number_of_layers_found.ravel(),0))
        #print "Number of points with more than 0 layers: ",x.shape[0]

        cal_data_ok = Numeric.greater(caliop_max_height,-9.)
        cal_lat_ok = Numeric.repeat(caObj.calipso.latitude[::],cal_data_ok)
        cal_avhrr_ctth_ok = Numeric.repeat(avhrr_ctth_cal_ok[::],cal_data_ok)
        cal_topmidlay_temp_ok = Numeric.repeat(caliop_max_height_midlaytemp[::],cal_data_ok)
        cal_toplay_thickness_ok = Numeric.repeat(caliop_toplay_thickness[::],cal_data_ok)
        cal_maxheight_ok = Numeric.repeat(caliop_max_height[::],cal_data_ok)
        cal_surftemp_ok = Numeric.repeat(caObj.avhrr.surftemp[::],cal_data_ok)
        cal_bt11temp_ok = Numeric.repeat(caObj.avhrr.bt11micron[::],cal_data_ok)

        # Transfer CloudSat MODIS cloud flag to CALIPSO representation

        cal_MODIS_cflag = Numeric.zeros(len(cal_data_ok),'b')
        for i in range(len(cal_data_ok)):
            if not CALIPSO_DISPLACED: 
                cloudsat_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)
                #print len(clsatObj.cloudsat.MODIS_cloud_flag),cloudsat_index
                if len(clsatObj.cloudsat.MODIS_cloud_flag) > (cloudsat_index):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
            else:
                cloudsat_index = int((i+calipso_displacement)/CLOUDSAT_TRACK_RESOLUTION + 0.5)
                if len(clsatObj.cloudsat.MODIS_cloud_flag) > (cloudsat_index-1):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
                    
        # Now calculate cloud emissivity Ec for topmost CALIOP and CloudSat layer

        #    Ec = Numeric.zeros(cal_lat_ok.shape[0],'d')
        Ec = Numeric.zeros(cal_data_ok.shape[0],'d')
        dum1,cwnum,dum2=get_central_wavenumber(noaa_number,273.15)
        
        #    for i in range(cal_lat_ok.shape[0]):
        for i in range(cal_data_ok.shape[0]):
            if cal_data_ok[i]:
                I=tb2radiance_using_central_wavenumber_klm(cwnum,caObj.avhrr.bt11micron[i],noaa_number,'4')
                Iclear=tb2radiance_using_central_wavenumber_klm(cwnum,caObj.avhrr.surftemp[i],noaa_number,'4')
                if I > Iclear:
                    Iclear = I   # Just avoiding too much mismatch between forecasted and real surface temps
                BTc=tb2radiance_using_central_wavenumber_klm(cwnum,caliop_max_height_midlaytemp[i],noaa_number,'4')
                if BTc < Iclear:
                    Ec[i]=(I-Iclear)/(BTc-Iclear)
                else:
                    Ec[i]=-9.0 # Give up on all temperature inversion cases!
            else:
                Ec[i]=-9.0                   
    
        if process_mode is 'EMISSFILT': #Apply only when cloud heights are above EMISS_MIN_HEIGHT
            #print "Emissivity filtering applied!"
            caliop_min_height_ok = Numeric.greater(caliop_max_height, setup.EMISS_MIN_HEIGHT)
            emissfilt_calipso_ok = Numeric.logical_or(Numeric.logical_and(Numeric.greater(Ec, setup.EMISS_LIMIT),caliop_min_height_ok),Numeric.logical_or(Numeric.equal(caliop_max_height,-9.),Numeric.less_equal(caliop_max_height, setup.EMISS_MIN_HEIGHT)))
    ##########################################################################################################################################
    #pdb.set_trace()
    ### 1 KM DATA CWC-RVOD ###                       
    elif int(Resolution[0])==1 and cloudsat_type=='CWC-RVOD':   
        elevationcwc = Numeric.where(Numeric.greater(clsatObj.cloudsatcwc.elevation,0),
                            clsatObj.cloudsatcwc.elevation,-9)

        data_okcwc = Numeric.ones(clsatObj.cloudsatcwc.elevation.shape,'b')
    ### 5 KM DATA ###
    elif int(Resolution[0])==5 and cloudsat_type=='GEOPROF':
        ##########################################################################################################################################
        #pdb.set_trace()
        CALIPSO_DISPLACED = 0
        latdiff = abs(clsatObj.cloudsat5km.latitude[0] - caObj.calipso5km.latitude[0])
        print "latdiff: ", latdiff
        if latdiff > 0.1:
            print "WARNING: CloudSat/CALIPSO startpoint differ by %f degrees!" % latdiff
            print clsatObj.cloudsat5km.latitude[0],clsatObj.cloudsat5km.longitude[0],caObj.calipso5km.latitude[0],caObj.calipso5km.longitude[0]
            CALIPSO_DISPLACED = 1
            for j in range(clsatObj.cloudsat5km.latitude.shape[0]):
                if (abs(clsatObj.cloudsat5km.latitude[j] - caObj.calipso5km.latitude[0]) < 0.05) and\
                       (abs(clsatObj.cloudsat5km.longitude[j] - caObj.calipso5km.longitude[0]) < 0.1):
                    calipso_displacement=int(j*CLOUDSAT5KM_TRACK_RESOLUTION)
                    print "CALIPSO_DISPLACEMENT: ", calipso_displacement
                    break       
        # First make sure that PPS cloud top heights are converted to height above sea level
        # just as CloudSat and CALIPSO heights are defined. Use corresponding DEM data.            
        elevation = Numeric.where(Numeric.greater(clsatObj.cloudsat5km.elevation,0),
                    clsatObj.cloudsat5km.elevation,-9)			# If clsatObj.cloudsat.elevation is bigger then 0 elevation(i,j)=clsatObj.cloudsat.elevation(i,j) else the value =-9

        #elevationcwc = Numeric.where(Numeric.greater(clsatObj.cloudsatcwc.elevation,0),
        #            clsatObj.cloudsatcwc.elevation,-9)

        #data_okcwc = Numeric.ones(clsatObj.cloudsatcwc.elevation.shape,'b')

        data_ok = Numeric.ones(clsatObj.cloudsat5km.elevation.shape,'b')
        print "Length of CLOUDSAT array: ", len(data_ok)
        lat_ok = Numeric.repeat(clsatObj.cloudsat5km.latitude[::],data_ok)
        avhrr_ctth_csat_ok = Numeric.repeat(clsatObj.avhrr5km.ctth_height[::],data_ok)
        avhrr_ctth_csat_ok = Numeric.where(Numeric.greater(avhrr_ctth_csat_ok,0.0),avhrr_ctth_csat_ok[::]+elevation*1.0,avhrr_ctth_csat_ok)

        # CALIPSO

        cal_elevation = Numeric.where(Numeric.greater(caObj.calipso5km.elevation,0),
                    caObj.calipso5km.elevation,-9)
        cal_data_ok = Numeric.ones(caObj.calipso5km.elevation.shape,'b')
        print "Length of CALIOP array: ", len(cal_data_ok)
        avhrr_ctth_cal_ok = Numeric.repeat(caObj.avhrr5km.ctth_height[::],cal_data_ok)
        avhrr_ctth_cal_ok = Numeric.where(Numeric.greater(avhrr_ctth_cal_ok,0.0),avhrr_ctth_cal_ok[::]+cal_elevation,avhrr_ctth_cal_ok) 
        ##########################################################################################################################################
        #pdb.set_trace()
        if (len(cal_data_ok) == 0) or (len(data_ok) == 0):
            print "Processing stopped: Zero lenght of matching arrays!"
            print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit()
        else:
            # If everything is OK, now create filename for statistics output file and open it for writing.
            # Notice that more than one file
            # (but maximum 2) can be created for one particular noaa orbit.
            resultpath = "%s/%s/%ikm/%s/%s/%s/%s" % (setup.RESULT_DIR, base_sat,int(Resolution[0]), base_year, base_month,AREA,process_mode)
            if not os.path.exists(resultpath):
                os.makedirs(resultpath)
            statfilename = "%s/5km_%s_cloudsat_calipso_avhrr_stat.dat" % (resultpath,basename)

            statfile = open(statfilename,"w")
            if process_mode == "BASIC":
                statfile.write("CloudSat min and max time diff: %f %f \n" %(clsat_min_diff,clsat_max_diff))
                statfile.write("CALIPSO min and max time diff: %f %f \n" %(ca_min_diff,ca_max_diff))
            else:
                statfile.write("CloudSat min and max time diff: See results for BASIC! \n")
                statfile.write("CALIPSO min and max time diff: See results for BASIC! \n")
            statfile.write("Start-Stop-Length Cloudsat: %f %f %f %f %s \n" %(clsatObj.cloudsat5km.latitude[0],clsatObj.cloudsat5km.longitude[0],clsatObj.cloudsat5km.latitude[len(clsatObj.cloudsat5km.latitude)-1],clsatObj.cloudsat5km.longitude[len(clsatObj.cloudsat5km.latitude)-1],len(data_ok)))
            statfile.write("Start-Stop-Length CALIPSO: %f %f %f %f %s \n" %(caObj.calipso5km.latitude[0],caObj.calipso5km.longitude[0],caObj.calipso5km.latitude[len(caObj.calipso5km.latitude)-1],caObj.calipso5km.longitude[len(caObj.calipso5km.latitude)-1],len(cal_data_ok)))
        
        # Extract CALIOP Vertical Feature Classification (bits 10-12) from 16 bit representation
        # for topmost cloud layer
        
        cal_vert_feature = Numeric.ones(caObj.calipso5km.cloud_top_profile[0,::].shape)*-9
        feature_array = 4*Numeric.bitwise_and(Numeric.right_shift(caObj.calipso5km.feature_classification_flags[0,::],11),1) + 2*Numeric.bitwise_and(Numeric.right_shift(caObj.calipso5km.feature_classification_flags[0,::],10),1) + Numeric.bitwise_and(Numeric.right_shift(caObj.calipso5km.feature_classification_flags[0,::],9),1)

        cal_vert_feature = Numeric.where(Numeric.not_equal(caObj.calipso5km.feature_classification_flags[0,::],1),feature_array[::],cal_vert_feature[::])
        
        # Prepare for plotting, cloud emissivity and statistics calculations
                
        midlayer_temp_kelvin = caObj.calipso5km.cloud_mid_temperature + 273.15
        caliop_toplay_thickness = Numeric.ones(caObj.calipso5km.cloud_top_profile[0,::].shape)*-9

        caliop_max_height_midlaytemp = Numeric.where(Numeric.greater(caObj.calipso5km.cloud_top_profile[0,::],-9),midlayer_temp_kelvin[0,::],-9)
        thickness = (caObj.calipso5km.cloud_top_profile[0,::]-caObj.calipso5km.cloud_base_profile[0,::])*1000.
        caliop_toplay_thickness = Numeric.where(Numeric.greater(caObj.calipso5km.cloud_top_profile[0,::],-9),thickness,caliop_toplay_thickness)
        ##########################################################################################################################################   
        #pdb.set_trace()
        (new_cloud_top,new_cloud_base) = CloudsatCloudOpticalDepth(caObj.calipso5km.cloud_top_profile, caObj.calipso5km.cloud_base_profile,caObj.calipso5km.optical_depth)
        caObj.calipso5km.cloud_top_profile = new_cloud_top
        caObj.calipso5km.cloud_base_profile = new_cloud_base
        caliop_height = []
        caliop_base = []
        caliop_max_height = Numeric.ones(caObj.calipso5km.cloud_top_profile[0,::].shape)*-9
        
        for i in range(10):
            hh = Numeric.where(Numeric.greater(caObj.calipso5km.cloud_top_profile[i,::],-9),
                                caObj.calipso5km.cloud_top_profile[i,::] * 1000.,-9)
                                            
            caliop_max_height = Numeric.maximum(caliop_max_height,
                                                caObj.calipso5km.cloud_top_profile[i,::] * 1000.)
            # This is actually unnecessary - we know that layer 1 is always the highest layer!!
            # However, arrays caliop_height and caliop_base are needed later for plotting/ KG

            caliop_height.append(hh)
            bb = Numeric.where(Numeric.greater(caObj.calipso5km.cloud_base_profile[i,::],-9),
                                caObj.calipso5km.cloud_base_profile[i,::] * 1000.,-9)
            caliop_base.append(bb)
            thickness = hh - bb
            
        x = Numeric.repeat(caObj.calipso5km.number_of_layers_found.ravel(),
                            Numeric.greater(caObj.calipso5km.number_of_layers_found.ravel(),0))
        #print "Number of points with more than 0 layers: ",x.shape[0]
                        ##########################################################################################################################################   
        #pdb.set_trace()    
        cal_data_ok = Numeric.greater(caliop_max_height,-9.)
        cal_lat_ok = Numeric.repeat(caObj.calipso5km.latitude[::],cal_data_ok)
        cal_avhrr_ctth_ok = Numeric.repeat(avhrr_ctth_cal_ok[::],cal_data_ok)
        cal_topmidlay_temp_ok = Numeric.repeat(caliop_max_height_midlaytemp[::],cal_data_ok)
        cal_toplay_thickness_ok = Numeric.repeat(caliop_toplay_thickness[::],cal_data_ok)
        cal_maxheight_ok = Numeric.repeat(caliop_max_height[::],cal_data_ok)
        cal_surftemp_ok = Numeric.repeat(caObj.avhrr5km.surftemp[::],cal_data_ok)
        cal_bt11temp_ok = Numeric.repeat(caObj.avhrr5km.bt11micron[::],cal_data_ok)

        # Transfer CloudSat MODIS cloud flag to CALIPSO representation

        cal_MODIS_cflag = Numeric.zeros(len(cal_data_ok),'b')
        
        for i in range(len(cal_data_ok)):
            if not CALIPSO_DISPLACED: 
                cloudsat_index = int(i/CLOUDSAT5KM_TRACK_RESOLUTION + 0.5)
                #print len(clsatObj.cloudsat.MODIS_cloud_flag),cloudsat_index
                if len(clsatObj.cloudsat5km.MODIS_cloud_flag) > (cloudsat_index):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat5km.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
            else:
                cloudsat_index = int((i+calipso_displacement)/CLOUDSAT5KM_TRACK_RESOLUTION + 0.5)
                if len(clsatObj.cloudsat5km.MODIS_cloud_flag) > (cloudsat_index-1):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat5km.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
        
        # Now calculate cloud emissivity Ec for topmost CALIOP and CloudSat layer

        #    Ec = Numeric.zeros(cal_lat_ok.shape[0],'d')
        Ec = Numeric.zeros(cal_data_ok.shape[0],'d')
        dum1,cwnum,dum2=get_central_wavenumber(noaa_number,273.15)
        #    for i in range(cal_lat_ok.shape[0]):
        for i in range(cal_data_ok.shape[0]):
            if cal_data_ok[i]:
                I=tb2radiance_using_central_wavenumber_klm(cwnum,caObj.avhrr5km.bt11micron[i],noaa_number,'4')
                Iclear=tb2radiance_using_central_wavenumber_klm(cwnum,caObj.avhrr5km.surftemp[i],noaa_number,'4')
                if I > Iclear:
                    Iclear = I   # Just avoiding too much mismatch between forecasted and real surface temps
                BTc=tb2radiance_using_central_wavenumber_klm(cwnum,caliop_max_height_midlaytemp[i],noaa_number,'4')
                if BTc < Iclear:
                    Ec[i]=(I-Iclear)/(BTc-Iclear)
                else:
                    Ec[i]=-9.0 # Give up on all temperature inversion cases!
            else:
               Ec[i]=-9.0 
  
        if process_mode is 'EMISSFILT': #Apply only when cloud heights are above EMISS_MIN_HEIGHT
            #print "Emissivity filtering applied!"
            caliop_min_height_ok = Numeric.greater(caliop_max_height, setup.EMISS_MIN_HEIGHT)
            emissfilt_calipso_ok = Numeric.logical_or(Numeric.logical_and(Numeric.greater(Ec, setup.EMISS_LIMIT),caliop_min_height_ok),Numeric.logical_or(Numeric.equal(caliop_max_height,-9.),Numeric.less_equal(caliop_max_height, setup.EMISS_MIN_HEIGHT)))         
            
    ### 5 KM DATA CWC-RVOD ###                       
    elif int(Resolution[0])==5 and cloudsat_type=='CWC-RVOD':   
        elevationcwc = Numeric.where(Numeric.greater(clsatObj.cloudsat5kmcwc.elevation,0),
                            clsatObj.cloudsat5kmcwc.elevation,-9)

        data_okcwc = Numeric.ones(clsatObj.cloudsat5kmcwc.elevation.shape,'b')   
##    print "lat,lon start CALIPSO: ", caObj.calipso.latitude[0],caObj.calipso.longitude[0]
##    latdiff = abs(clsatObj.cloudsat.latitude[len(clsatObj.cloudsat.latitude)-1] - caObj.calipso.latitude[len(caObj.calipso.latitude)-1])
##    if latdiff > 0.1:
##        print "WARNING: CloudSat/CALIPSO endpoint differ by %f degrees!" % latdiff
##    print "lat,lon end Cloudsat: ", clsatObj.cloudsat.latitude[len(clsatObj.cloudsat.latitude)-1],clsatObj.cloudsat.longitude[len(clsatObj.cloudsat.latitude)-1]
##    print "lat,lon end CALIPSO: ", caObj.calipso.latitude[len(caObj.calipso.latitude)-1],caObj.calipso.longitude[len(caObj.calipso.latitude)-1]



        # Cloud emissivities Ec are calculated as follows:
        
        #                     Ec = (I-Iclear)/(B(Tc)-Iclear)


#    size_test=clsatObj.cloudsat.cloud_mask.shape
#    cloud=clsatObj.cloudsat.cloud_mask
#    tot=Numeric.zeros((size_test[0],size_test[1]),'d')
#    tot_cero=Numeric.zeros((size_test[0],size_test[1]),'d')
#    cloud_cero=Numeric.zeros((size_test[0],size_test[1]),'d')
#    result=Numeric.zeros((size_test[0],size_test[1]),'d')
#    for row in range(size_test[0]):
#        for col in range(size_test[1]):
#            tot[row,col]=clsatObj.cloudsatcwc.RVOD_ice_water_content[row,col]+clsatObj.cloudsatcwc.RVOD_liq_water_content[row,col] 
#	    if tot[row,col]==0:	        
#		tot_cero[row,col]=0
#	    else:
#		tot_cero[row,col]=1
#	    if cloud[row,col]<=20:
#		cloud_cero[row,col]=0
#	    else:
#                cloud_cero[row,col]=5
#
#	    result[row,col]=tot_cero[row,col]-cloud_cero[row,col]
#        
#	    if result[row,col]==-4:
#		result[row,col]=0
#    if max(max(result))==0 and min(min(result))==0:
#	print "IWC + LWC contains same pixlar as those with cloud"
#    else:
#        if max(max(result))!=0:
#	    print "IWC + LWC contains pixlar that are not represented in cloud produkt"
#	if min(min(result))==-5:
#	    print "Cloud produkt contains pixlar that are not represented in IWC + LWC"

 ##########################################################################################################################################   
    #pdb.set_trace()

    #==============================================================================================
    #Draw plot
    if process_mode in setup.PLOT_MODES:
    ##########################################################################################################################################   
        #pdb.set_trace()
        plotpath = "%s/%s/%ikm/%s/%s/%s" %(setup.PLOT_DIR, base_sat, int(Resolution[0]), base_year, base_month, AREA)
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
            
        #trajectorypath = "%s/trajectory_plot/%ikm/%s" %(MAIN_RUNDIR,int(Resolution[0]),AREA)         
        trajectorypath = "%s/trajectory_plot" %(plotpath)
        if not os.path.exists(trajectorypath):
                os.makedirs(trajectorypath)
        trajectoryname = "%s/%skm_%s_trajectory" %(trajectorypath,int(Resolution[0]),basename)
          ##########################################################################################################################################   
        #pdb.set_trace()     
        # To make it possible to use the same function call to drawCalClsatGEOPROFAvhrr*kmPlot
        # in any processing mode:
        if 'emissfilt_calipso_ok' not in locals():
            emissfilt_calipso_ok = None
        
        if int(Resolution[0])==1 and cloudsat_type=='GEOPROF': 
            drawCalClsatAvhrrPlotSATZ(cllat, clsatObj.avhrr.satz, caObj.avhrr.satz, plotpath, basename, Resolution[0])
            drawCalClsatGEOPROFAvhrr1kmPlot(clsatObj, caObj, elevation, data_ok,
                                            CALIPSO_DISPLACED, caliop_base,
                                            caliop_height, cal_data_ok,
                                            avhrr_ctth_cal_ok, plotpath,
                                            basename, process_mode, emissfilt_calipso_ok)
            drawCalClsatAvhrrPlotTimeDiff(cllat, clsatObj.diff_sec_1970, caObj.diff_sec_1970, plotpath, basename, Resolution[0])
            plotSatelliteTrajectory(cllon,cllat,calon,calat,avhrlon,avhrlat,trajectoryname)
            
        elif int(Resolution[0])==1 and cloudsat_type=='CWC-RVOD':
            drawCalClsatAvhrrPlotTimeDiff(cllat, clsatObj.diff_sec_1970, caObj.diff_sec_1970, plotpath, basename, Resolution[0])
            phase='LW'  
            drawCalClsatCWCAvhrr1kmPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase) #caObj, CALIPSO_DISPLACED, caliop_base, caliop_height, cal_data_ok, avhrr_ctth_cal_ok)
            phase='IW'  
            drawCalClsatCWCAvhrr1kmPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase) #caObj, CALIPSO_DISPLACED, caliop_base, caliop_height, cal_data_ok, avhrr_ctth_cal_ok)
        elif int(Resolution[0])==5 and cloudsat_type=='GEOPROF':
            drawCalClsatAvhrrPlotSATZ(cllat, clsatObj.avhrr5km.satz, caObj.avhrr5km.satz, plotpath, basename, Resolution[0])
            drawCalClsatGEOPROFAvhrr5kmPlot(clsatObj, caObj, elevation, data_ok,
                                            CALIPSO_DISPLACED, caliop_base,
                                            caliop_height, cal_data_ok,
                                            avhrr_ctth_cal_ok, plotpath,
                                            basename, process_mode, emissfilt_calipso_ok)
            drawCalClsatAvhrrPlotTimeDiff(cllat, clsatObj.diff_sec_1970, caObj.diff_sec_1970, plotpath, basename, Resolution[0])
            plotSatelliteTrajectory(cllon,cllat,calon,calat,avhrlon,avhrlat,trajectoryname)          
        elif int(Resolution[0])==5 and cloudsat_type=='CWC-RVOD':
            phase='LW'
            drawCalClsatCWCAvhrr5kmPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase)
            phase='IW'
            drawCalClsatCWCAvhrr5kmPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase)
    #================================================================================================
    #Calculate Statistics
    if cloudsat_type=='GEOPROF':
        if process_mode is 'EMISSFILT':
            process_calipso_ok = emissfilt_calipso_ok
        else:
            process_calipso_ok = 0

        CalculateStatistics(process_mode, clsatObj, statfile, caObj, cal_MODIS_cflag,
                            cal_vert_feature, avhrr_ctth_csat_ok, data_ok,
                            cal_data_ok, avhrr_ctth_cal_ok, caliop_max_height,
                            process_calipso_ok, Resolution)
         

    

    

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 17:
        write_log("INFO","Usage: %s <3*cloudsat-hdf5-file> <6*calipso-hdf5-file> <pps cloudtype file> <pps ctth file> <pps avhrr file> <pps surftemp file> <pps sunsatangle file> <processing_mode> <resolution>" %sys.argv[0],moduleid=MODULE_ID)
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9)

    cloudsatfile = sys.argv[1:4]
    calipsofile = sys.argv[4:-7]
    ctypefile = sys.argv[-7]
    ctthfile = sys.argv[-6]
    avhrrfile = sys.argv[-5]
    surftfile = sys.argv[-4]
    sunanglefile = sys.argv[-3]
    process_mode = sys.argv[-2]
    Resolution = sys.argv[-1]
    
    run(cloudsatfile, calipsofile, ctypefile, ctthfile, avhrrfile, surftfile, sunanglefile, process_mode, Resolution)