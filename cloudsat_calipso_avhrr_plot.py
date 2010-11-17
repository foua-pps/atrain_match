# -*- coding: utf-8 -*-
#Program cloudsat_calipso_avhrr_plot.py

from config import MAIN_DIR, SUB_DIR, MAIN_RUNDIR, CLOUDSAT_TRACK_RESOLUTION, CLOUDSAT_CLOUDY_THR, CLOUDSAT5KM_TRACK_RESOLUTION, AZIMUTH_RANGE, RESOLUTION
import pdb
import os, string
import sys
import rpy
import numpy

def format_title(title):
    """
    Format *title*, possibly adding satellite azimuth range (depending on range).
    
    Returns formatted title.
    """
    if AZIMUTH_RANGE[0] > 0 or AZIMUTH_RANGE[1] < 90:
        title = "%s (satz range [%.1f, %.1f] deg)" % (title, AZIMUTH_RANGE[0], AZIMUTH_RANGE[1])
    
    return title
# -----------------------------------------------------
def drawCalClsatGEOPROFAvhrrPlot(clsatObj_cloudsat, clsatObj_avhrr, caObj_calipso, caObj_avhrr, elevation, data_ok,
                                    CALIPSO_DISPLACED, caliop_base,
                                    caliop_height, cal_data_ok,
                                    avhrr_ctth_cal_ok, plotpath, basename,
                                    mode, emissfilt_calipso_ok=None, file_type='eps'):
    
    # Prepare for Avhrr
    #Second version for avhrr ctth plotting based on data along CALIPSO track
    if mode is 'EMISSFILT':
        cal_plot_data_ok = numpy.logical_and(emissfilt_calipso_ok,numpy.greater(avhrr_ctth_cal_ok[::],0))
    else:
        cal_plot_data_ok = numpy.greater(avhrr_ctth_cal_ok[::],0)
    dummy=caObj_calipso.latitude.shape[0]
    pixel_position=numpy.arange(dummy)
    pixel_position_ok = numpy.repeat(pixel_position[::],cal_plot_data_ok)
    avhrr_ctth_ok = numpy.repeat(avhrr_ctth_cal_ok[::],cal_plot_data_ok)    
    
    
    # Calculates Hihest Cloud Top   
    MAXHEIGHT_CLOUDSAT = numpy.nanmax(numpy.where(clsatObj_cloudsat.cloud_mask>=CLOUDSAT_CLOUDY_THR,clsatObj_cloudsat.Height,numpy.nan))+120
    MAXHEIGHT_CALIPSO = numpy.nanmax(caliop_height)
    MAXHEIGHT_AVHRR = numpy.nanmax(avhrr_ctth_ok)
    MAXHEIGHT = numpy.nanmax([MAXHEIGHT_CALIPSO,MAXHEIGHT_CLOUDSAT,MAXHEIGHT_AVHRR])
    
    if file_type == 'eps' or file_type == 'ps':
        device = rpy.r.postscript
    elif file_type == 'jpeg' or file_type == 'jpg':
       device = rpy.r.jpeg
    elif file_type == 'npg':
        device = rpy.r.npg
    elif file_type == 'bmp':
        device = rpy.r.bmp
    elif file_type == 'pdf':
        device = rpy.r.pdf
    elif file_type == 'tiff':
        device = rpy.r.tiff   
    else:
        print('only formats eps, ps, jpeg, jpg, npg, bmp, pdf and tiff are suported')
        sys.exit()
    device("%s/%skm_%s_cloudsat_calipso_avhrr_clouds.%s"%(plotpath,RESOLUTION,basename,file_type),width=900,height=400) 

    # Let's use the CloudSat track as our reference
    dummy=clsatObj_cloudsat.latitude.shape[0]
    #print "Length of cloudsat array: ", dummy
    geometric_range_CloudSat = int(dummy*CLOUDSAT_TRACK_RESOLUTION)
    #print "Geometric length of cloudsat array: ", geometric_range_CloudSat
    pixel_position_plain=numpy.arange(dummy)
    
    pixel_position_geometric=numpy.arange(geometric_range_CloudSat)
    dummy_elevation_adjusted=numpy.zeros((geometric_range_CloudSat),'d')
    #time_diff_cloudsat=numpy.zeros((geometric_range_CloudSat),'d')
    for i in range(geometric_range_CloudSat): #Just extend original elevation array
        elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)
        dummy_elevation_adjusted[i] = elevation[elevation_index] # Scales elevation
        #time_diff_cloudsat[i] = clsatObj_diff_sec_1970[elevation_index] #For time ploting
    #dummy_elevation_adjusted = numpy.where(dummy_elevation_adjusted==-9,numpy.nan,dummy_elevation_adjusted)
    dummy = rpy.r.plot(pixel_position_geometric[::],dummy_elevation_adjusted,
                        ylim=[0,MAXHEIGHT],
                        xlab="Track Position",ylab="Cloud Height",
                        main=format_title("AVHRR-CloudSat-Caliop Cloud Top Heights"),
                        col="black",type='l',lwd=1)

    
    # Makes the surface black
    for i in range(geometric_range_CloudSat): #Let's also plot elevation in the vertical
        if dummy_elevation_adjusted[i] > 0.0:
            rpy.r.lines([i,i],[0.0,dummy_elevation_adjusted[i]],
                            col="black",lty=1,lwd=1)    
    
    # Colors:
    colors=[]
    for i in range(10):
        colors.append(rpy.r.rgb([(40-i)*6],[(40-i)*6],[(40-i)*6],maxColorValue=255))
    for i in range(10,30):
        colors.append(rpy.r.rgb([100+(i-10)*5],[100+(i-10)*5],[100-(i-10)*5],maxColorValue=255))
    for i in range(30,100):
        colors.append(rpy.r.rgb([255],[50],[50],maxColorValue=255))
    
    # Plot CloudSat
    clsat_max_height = numpy.repeat(clsatObj_cloudsat.Height[124,::],data_ok)
    
    for i in range(125):
        height = clsatObj_cloudsat.Height[i,::]
        cmask_ok = clsatObj_cloudsat.cloud_mask[i,::]
        
        for idx in range(len(pixel_position_plain)):
            nidx = int(cmask_ok[idx]+0.5)
            if nidx == 0:
                continue
            if height[idx] < 240*4: #or height[idx] > MAXHEIGHT:
                continue
            base_height = height[idx]-120
            top_height = height[idx]+120
            if nidx >= int(CLOUDSAT_CLOUDY_THR):
                rpy.r.lines([int(pixel_position_plain[idx]*CLOUDSAT_TRACK_RESOLUTION),\
                                int(pixel_position_plain[idx]*CLOUDSAT_TRACK_RESOLUTION)],[base_height,top_height],\
                            col=colors[nidx],lty=1,lwd=1)    
                clsat_max_height[idx] = max(clsat_max_height[idx],top_height)
    
    # Plot Caliop 
    calipso_displacement = 0
    if CALIPSO_DISPLACED: # Search for startpoint along CloudSat array
        for i in range(clsatObj_cloudsat.latitude.shape[0]):
            if (abs(clsatObj_cloudsat.latitude[i] - caObj_calipso.latitude[0]) < 0.005) and\
                (abs(clsatObj_cloudsat.longitude[i] - caObj_calipso.longitude[0]) < 1.0):
                calipso_displacement=int(i*CLOUDSAT_TRACK_RESOLUTION)
                break
    
    dummy=caObj_calipso.latitude.shape[0]
    pixel_position=numpy.arange(dummy)
    for i in range(10):
        base_ok = caliop_base[i]
        top_ok = caliop_height[i]
        for idx in range(len(pixel_position)):
            if cal_data_ok[idx]:
                if caObj_calipso.number_of_layers_found[idx] > i:
                    rpy.r.lines([pixel_position[idx]+calipso_displacement,pixel_position[idx]+calipso_displacement],[base_ok[idx],top_ok[idx]],
                                col="green",lty=1,lwd=0.5)

    
    
    #Plot Avhrr   
    rpy.r.points(pixel_position_ok+calipso_displacement,avhrr_ctth_ok,col="blue",pch="+",cex=0.7)
    # Close Figure
    rpy.r.dev_off()

# -----------------------------------------------------
def drawCalClsatCWCAvhrr1kmPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase):#, caObj, CALIPSO_DISPLACED, caliop_base, caliop_height, cal_data_ok, avhrr_ctth_cal_ok):
    import pdb
    import os, string
    import sys
    import rpy
    import numpy


    
    # -------------------------------------------------------------------------------------------------------------------------------
    # Create environment for plotting LWP and IWP        pixel_position_plain=numpy.arange(dummy)
	
    if phase=='IW':
        dataP=clsatObj.cloudsatcwc.RVOD_ice_water_path
        dataC=clsatObj.cloudsatcwc.RVOD_ice_water_content
    else:
        dataP=clsatObj.cloudsatcwc.RVOD_liq_water_path
        dataC=clsatObj.cloudsatcwc.RVOD_liq_water_content

    rpy.r.postscript("%s/1km_%s_calipso_%sP.eps"%(plotpath,basename,phase),width=900,height=400)	
    storlek_x=clsatObj.cloudsatcwc.latitude.shape[0] 				# Decides the size of the array latitude
    andrad_storlek_x = int(storlek_x*CLOUDSAT_TRACK_RESOLUTION) 		# Changes the resolution with a factor of 1.076
    #geometric_range_CloudSatcwc = dummycwc
    pixel_position=numpy.arange(storlek_x) 				# Creates an array with numbers 0 -> storlek_x
    andrad_pixel_position=numpy.arange(andrad_storlek_x)	# Creates an array with size geometric_range_CloudSatcwc numbers from 0->size
    y_zeros=numpy.zeros((storlek_x),'d')		# Creates an array with size andrad_storlek_x Filled with zeros of type double
    andrad_y_zeros=numpy.zeros(andrad_storlek_x,'d')		# Creates an array with size andrad_storlek_x Filled with zeros of type double
    andrad_dataP=numpy.ones(andrad_storlek_x)*-1
    for iP in range(storlek_x):
        niP=int(iP*CLOUDSAT_TRACK_RESOLUTION)
        andrad_dataP[niP]=dataP[iP]

    y=andrad_y_zeros.copy()

    for i in range(len(andrad_dataP)):
        if (andrad_dataP[i]<=0):
            y[i]=float('nan')
        else:
            y[i]=andrad_dataP[i]
	
    MAXVALUE=int(max(dataP)+100)							# Findes max value and add 100
	
    rpy.r.plot(andrad_pixel_position[::],y,ylim=[0,MAXVALUE],xlab="Track Position",ylab="%sP [g/m^2]"%(phase),main="CloudSat %sP"%(phase),col="black",type='p',lwd=1)
    rpy.r.dev_off()
    # ----------------------------------------------------------------------------------------------------------------------------------
    
    rpy.r.postscript("%s/1km_%s_calipso_%sP.eps"%(plotpath,'test',phase),width=900,height=400)
    yy=numpy.arange(1)
    ztemp=numpy.asmatrix(y.copy())
    z=numpy.reshape(ztemp,(-1,1))
    rpy.r.image(andrad_pixel_position,yy,z,xlab="Track Position",ylab="%sP [g/m^2]"%(phase),main="CloudSat %sP"%(phase))    
    rpy.r.dev_off()

    # ----------------------------------------------------------------------------------------------------------------------------------
    # Create environment for plotting LWC and IWC
    rpy.r.postscript("%s/1km_%s_calipso_%sC.eps"%(plotpath,basename,phase),width=900,height=400)
	#max(clsatObj.cloudsatcwc.RVOD_ice_water_content)
    clsatcwc_max_height=numpy.zeros([1,storlek_x],'d')
    y_adjust=andrad_y_zeros
    test=andrad_y_zeros
    for it in range(storlek_x):
        test[it]=int(pixel_position[it]*CLOUDSAT_TRACK_RESOLUTION)

    for i in range(andrad_storlek_x): #Just extend original elevation array
        elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)
        y_adjust[i] = elevationcwc[elevation_index] # Scales elevation

    rpy.r.plot(andrad_pixel_position[::],y_adjust, ylim=[0,MAXHEIGHT], xlab="Track Position",ylab="%sC Height"%(phase), main="CloudSat-%sC"%(phase), col="black",type='l',lwd=1)


	# Makes the surface black
    for i in range(andrad_storlek_x): #Let's also plot elevation in the vertical
        if y_adjust[i] > 0.0:
            rpy.r.lines([i,i],[0.0,y_adjust[i]],col="black",lty=1,lwd=1)  

#	for i in range(storlek_x):
#   clsatcwc_max_height[0,i]=int(clsatObj.cloudsatcwc.Height[124,i])
    clsatcwc_max_height[0,:]=numpy.repeat(clsatObj.cloudsatcwc.Height[124,::],data_okcwc)#Takes lowest Hight
    for i in range(125):
        heightcwc = clsatObj.cloudsatcwc.Height[i,::]
        RVOD_IWC = dataC[i,::]
        for idx in range(storlek_x):
            nidx = int(RVOD_IWC[idx]+0.5)
            if nidx == 0:	# If there are no value dont do anything
                continue
            if heightcwc[idx] < 240*4 or heightcwc[idx] > MAXHEIGHT:		# Treschold when measurements are valid
                continue
            base_heightcwc = heightcwc[idx]-120
            top_heightcwc = heightcwc[idx]+120
            if nidx >= int(CLOUDSAT_CLOUDY_THR):
                rpy.r.lines([int(pixel_position[idx]*CLOUDSAT_TRACK_RESOLUTION),\
                int(pixel_position[idx]*CLOUDSAT_TRACK_RESOLUTION)],[base_heightcwc,top_heightcwc],\
                col='red',lty=1,lwd=1)   
                clsatcwc_max_height[0,idx] = max(clsatcwc_max_height[0,idx],top_heightcwc) 

    rpy.r.dev_off()


# -----------------------------------------------------

def drawCalClsatAvhrrPlotTimeDiff(cllat, clsat_diff_sec_1970, cal_diff_sec_1970, plotpath, basename, res):
    import pdb
    import os, string
    import sys
    import rpy
    import numpy       
    # Plot time diff
        
    device = rpy.r.postscript
    device("%s/%skm_%s_time_diff.eps"%(plotpath,res,basename),width=900,height=400)
    # Let's use the CloudSat track as our reference
    dummy=cllat.shape[0]
    geometric_range_CloudSat = int(dummy*CLOUDSAT5KM_TRACK_RESOLUTION)
    
    x=numpy.arange(geometric_range_CloudSat)    
    y=numpy.ones(geometric_range_CloudSat)*-1000
    
    time_diff_cloudsat=numpy.zeros((geometric_range_CloudSat),'d')
    
    for i in range(geometric_range_CloudSat): #Just extend original time_diff array
        elevation_index = int(i/CLOUDSAT5KM_TRACK_RESOLUTION + 0.5)        
        time_diff_cloudsat[i] = clsat_diff_sec_1970[elevation_index] 
    
    MAXVALUE=(max(clsat_diff_sec_1970)/60+10)
    MINVALUE=(min(clsat_diff_sec_1970)/60-10)
    rpy.r.plot(x,y,ylim=[MINVALUE,MAXVALUE],
                        xlab="Track Position",ylab="Time diff [min]",
                        main=format_title("Time Difference Between Avhrr and CloudSat/Calispo"),
                        col="black",type='l',lwd=0)
    biggest_Cloudsat_diff_place = numpy.argmax(numpy.abs(time_diff_cloudsat))
    biggest_Cloudsat_diff = time_diff_cloudsat[biggest_Cloudsat_diff_place]
    biggest_Calipso_diff_place = numpy.argmax(numpy.abs(cal_diff_sec_1970))
    biggest_Calipso_diff = cal_diff_sec_1970[biggest_Calipso_diff_place]

    if time_diff_cloudsat.shape[0] >= cal_diff_sec_1970.shape[0]:
        length_time_diff = time_diff_cloudsat.shape[0]
    else:
        length_time_diff = cal_diff_sec_1970.shape[0]
        
    rpy.r.points(time_diff_cloudsat/60,col="red",pch="+",cex=0.7)
    rpy.r.points(cal_diff_sec_1970/60,col="green",pch="+",cex=0.7)
    rpy.r.legend((length_time_diff-1500),MAXVALUE-5,["CloudSat (max time diff %.2f min)" % round(biggest_Cloudsat_diff/60,2),"Calipso (max time diff %.2f min)" % round(biggest_Calipso_diff/60,2)],pch=[4,4],col=["red","green"])
    rpy.r.dev_off()
    
    
# -----------------------------------------------------   
def drawCalClsatAvhrrPlotSATZ(cllat, AvhrrClsatSatz, AvhrrCalSatz, plotpath, basename, res):
    import pdb
    import os, string
    import sys
    import rpy
    import numpy       
    # Plot time diff
        
    device = rpy.r.postscript
    device("%s/%skm_%s_satz.eps"%(plotpath,res,basename),width=900,height=400)
    # Let's use the CloudSat track as our reference
    dummy=cllat.shape[0]    
    geometric_range_CloudSat = int(dummy*CLOUDSAT5KM_TRACK_RESOLUTION)   
      
    x=numpy.arange(geometric_range_CloudSat)    
    y=numpy.ones(geometric_range_CloudSat)*-1000
    
    AvhrrClsatSatz_adjust=numpy.zeros((geometric_range_CloudSat),'d')
    
    for i in range(geometric_range_CloudSat): #Just extend original x array
        elevation_index = int(i/CLOUDSAT5KM_TRACK_RESOLUTION + 0.5)        
        AvhrrClsatSatz_adjust[i] = AvhrrClsatSatz[elevation_index] #For satz ploting
    
    MAXVALUE=(numpy.nanmax(AvhrrClsatSatz)+10)
    MINVALUE=(numpy.nanmin(AvhrrClsatSatz)-10)
    rpy.r.plot(x,y,ylim=[MINVALUE,MAXVALUE],                        
                        xlab="Track Position",ylab="satellite zenith angle [deg]",
                        main=format_title("Avhrr compared to CloudSat/Calipso"),
                        col="black",type='l',lwd=0)
        
    rpy.r.points(AvhrrClsatSatz_adjust,col="red",pch="+",cex=0.7)
    rpy.r.points(AvhrrCalSatz,col="green",pch="+",cex=0.7)
    rpy.r.legend(0,MAXVALUE,["CloudSat","Calipso"],pch=[4,4],col=["red","green"])
    rpy.r.dev_off()   
    
# -----------------------------------------------------
def drawCalClsatCWCAvhrr5kmPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase):#, caObj, CALIPSO_DISPLACED, caliop_base, caliop_height, cal_data_ok, avhrr_ctth_cal_ok):
    import pdb
    import os, string
    import sys
    import rpy
    import numpy


    
    # -------------------------------------------------------------------------------------------------------------------------------
    # Create environment for plotting LWP and IWP        pixel_position_plain=numpy.arange(dummy)
	
    if phase=='IW':
        dataP=clsatObj.cloudsat5kmcwc.RVOD_ice_water_path
        dataC=clsatObj.cloudsat5kmcwc.RVOD_ice_water_content
    else:
        dataP=clsatObj.cloudsat5kmcwc.RVOD_liq_water_path
        dataC=clsatObj.cloudsat5kmcwc.RVOD_liq_water_content

    rpy.r.postscript("%s/5km_%s_calipso_%sP.eps"%(plotpath,basename,phase),width=900,height=400)	
    storlek_x=clsatObj.cloudsat5kmcwc.latitude.shape[0] 				# Decides the size of the array latitude
    andrad_storlek_x = int(storlek_x*CLOUDSAT5KM_TRACK_RESOLUTION) 		# Changes the resolution with a factor of 1.076
    #geometric_range_CloudSatcwc = dummycwc
    pixel_position=numpy.arange(storlek_x) 				# Creates an array with numbers 0 -> storlek_x
    andrad_pixel_position=numpy.arange(andrad_storlek_x)	# Creates an array with size geometric_range_CloudSatcwc numbers from 0->size
    y_zeros=numpy.zeros((storlek_x),'d')		# Creates an array with size andrad_storlek_x Filled with zeros of type double
    andrad_y_zeros=numpy.zeros(andrad_storlek_x,'d')		# Creates an array with size andrad_storlek_x Filled with zeros of type double
    andrad_dataP=numpy.ones(andrad_storlek_x)*-1
    for iP in range(storlek_x):
        niP=int(iP*CLOUDSAT_TRACK_RESOLUTION)
        andrad_dataP[niP]=dataP[iP]

    y=andrad_y_zeros.copy()

    for i in range(len(andrad_dataP)):
        if (andrad_dataP[i]<=0):
            y[i]=float('nan')
        else:
            y[i]=andrad_dataP[i]
	
    MAXVALUE=int(max(dataP)+100)							# Findes max value and add 100
	
    rpy.r.plot(andrad_pixel_position[::],y,ylim=[0,MAXVALUE],xlab="Track Position",ylab="%sP [g/m^2]"%(phase),main="CloudSat %sP"%(phase),col="black",type='p',lwd=1)
    rpy.r.dev_off()
    # ----------------------------------------------------------------------------------------------------------------------------------
    
    rpy.r.postscript("%s/5km_%s_calipso_%sP.eps"%(plotpath,'test',phase),width=900,height=400)
    yy=numpy.arange(1)
    ztemp=numpy.asmatrix(y.copy())
    z=numpy.reshape(ztemp,(-1,1))
    rpy.r.image(andrad_pixel_position,yy,z,xlab="Track Position",ylab="%sP [g/m^2]"%(phase),main="CloudSat %sP"%(phase))    
    rpy.r.dev_off()

    # ----------------------------------------------------------------------------------------------------------------------------------
    # Create environment for plotting LWC and IWC
    rpy.r.postscript("%s/5km_%s_calipso_%sC.eps"%(plotpath,basename,phase),width=900,height=400)
	#max(clsatObj.cloudsatcwc.RVOD_ice_water_content)
    clsatcwc_max_height=numpy.zeros([1,storlek_x],'d')
    y_adjust=andrad_y_zeros
    test=andrad_y_zeros
    for it in range(storlek_x):
        test[it]=int(pixel_position[it]*CLOUDSAT_TRACK_RESOLUTION)

    for i in range(andrad_storlek_x): #Just extend original elevation array
        elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)
        y_adjust[i] = elevationcwc[elevation_index] # Scales elevation

    rpy.r.plot(andrad_pixel_position[::],y_adjust, ylim=[0,MAXHEIGHT], xlab="Track Position",ylab="%sC Height"%(phase), main="CloudSat-%sC"%(phase), col="black",type='l',lwd=1)


	# Makes the surface black
    for i in range(andrad_storlek_x): #Let's also plot elevation in the vertical
        if y_adjust[i] > 0.0:
            rpy.r.lines([i,i],[0.0,y_adjust[i]],col="black",lty=1,lwd=1)  

#	for i in range(storlek_x):
#   clsatcwc_max_height[0,i]=int(clsatObj.cloudsatcwc.Height[124,i])
    clsatcwc_max_height[0,:]=numpy.repeat(clsatObj.cloudsat5kmcwc.Height[124,::],data_okcwc)#Takes lowest Hight
    for i in range(125):
        heightcwc = clsatObj.cloudsat5kmcwc.Height[i,::]
        RVOD_IWC = dataC[i,::]
        for idx in range(storlek_x):
            nidx = int(RVOD_IWC[idx]+0.5)
            if nidx == 0:	# If there are no value dont do anything
                continue
            if heightcwc[idx] < 240*4 or heightcwc[idx] > MAXHEIGHT:		# Treschold when measurements are valid
                continue
            base_heightcwc = heightcwc[idx]-120
            top_heightcwc = heightcwc[idx]+120
            if nidx >= int(CLOUDSAT_CLOUDY_THR):
                rpy.r.lines([int(pixel_position[idx]*CLOUDSAT5KM_TRACK_RESOLUTION),\
                int(pixel_position[idx]*CLOUDSAT5KM_TRACK_RESOLUTION)],[base_heightcwc,top_heightcwc],\
                col='red',lty=1,lwd=1)   
                clsatcwc_max_height[0,idx] = max(clsatcwc_max_height[0,idx],top_heightcwc) 

    rpy.r.dev_off()
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
