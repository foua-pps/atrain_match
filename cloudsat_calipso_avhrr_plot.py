# -*- coding: utf-8 -*-
#Program cloudsat_calipso_avhrr_plot.py

import numpy
from config import AZIMUTH_RANGE, MAXHEIGHT, CLOUDSAT_CLOUDY_THR,\
    RESOLUTION, CLOUDSAT_TRACK_RESOLUTION

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

    from matplotlib import pyplot as plt
#    import rpy #@UnresolvedImport
#    from 
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
    
    
#    # Calculates Hihest Cloud Top   
    if MAXHEIGHT == None:
        maxheight_calipso = numpy.nanmax(caliop_height)
        maxheight_avhrr = numpy.nanmax(avhrr_ctth_ok)
        max_height_sat = [maxheight_calipso, maxheight_avhrr]
        if clsatObj_cloudsat != None:
            maxheight_cloudsat = numpy.nanmax(numpy.where(clsatObj_cloudsat.cloud_mask >= \
                                                          CLOUDSAT_CLOUDY_THR, \
                                                          clsatObj_cloudsat.Height, \
                                                          numpy.nan))+120
            max_height_sat.append(maxheight_cloudsat)
        maxheight = numpy.nanmax(max_height_sat)
        if maxheight < 12000.:
            maxheight = 12000
        elif maxheight < 25000.:
            maxheight = 25000.
        else:
            maxheight = maxheight + 1000
    else:
        maxheight = MAXHEIGHT

#    device("%s/%skm_%s_cloudsat_calipso_avhrr_clouds.%s"%(plotpath,RESOLUTION,basename,file_type),width=900,height=400) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if clsatObj_cloudsat != None:
        #: Let's use the CloudSat track as our reference
        dummy=clsatObj_cloudsat.latitude.shape[0]
        #print "Length of cloudsat array: ", dummy
        geometric_range_CloudSat = int(dummy*CLOUDSAT_TRACK_RESOLUTION)
        #print "Geometric length of cloudsat array: ", geometric_range_CloudSat
        pixel_position_plain=numpy.arange(dummy)
    
        pixel_position_geometric=numpy.arange(geometric_range_CloudSat)
        dummy_elevation_adjusted=numpy.zeros((geometric_range_CloudSat),'d')
        for i in range(geometric_range_CloudSat): #Just extend original elevation array
            elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)
            dummy_elevation_adjusted[i] = elevation[elevation_index] # Scales elevation
        #: Plots the Surface
        ax.plot(pixel_position_geometric[::],dummy_elevation_adjusted, 'k')

        #: Makes the surface black
        for i in range(geometric_range_CloudSat): #Let's also plot elevation in the vertical
            if dummy_elevation_adjusted[i] > 0.0:
                ax.vlines(i, 0.0, dummy_elevation_adjusted[i],
                                color="black", linestyle='solid', linewidth=1)    
        
        #: Colors
        colors=[]
        for i in range(10):
            colors.append(numpy.divide([(40-i)*6, (40-i)*6, (40-i)*6], 255.))
        for i in range(10,30):
            colors.append(numpy.divide([100+(i-10)*5, 100+(i-10)*5, 100-(i-10)*5], 255.))
        for i in range(30,100):
            colors.append(numpy.divide([255, 50, 50], 255.))
        
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
                    ax.vlines(int(pixel_position_plain[idx]*CLOUDSAT_TRACK_RESOLUTION),\
                                    base_height, top_height, color = colors[nidx], \
                                    linestyle = 'solid', linewidth = 1)    
                    clsat_max_height[idx] = max(clsat_max_height[idx],top_height)
        calipso_displacement = 0
        #: Search for startpoint along CloudSat array
        if CALIPSO_DISPLACED:
            for i in range(clsatObj_cloudsat.latitude.shape[0]):
                if (abs(clsatObj_cloudsat.latitude[i] - caObj_calipso.latitude[0]) < 0.005) and\
                    (abs(clsatObj_cloudsat.longitude[i] - caObj_calipso.longitude[0]) < 1.0):
                    calipso_displacement = int(i*CLOUDSAT_TRACK_RESOLUTION)
                    break
        title = "AVHRR-CloudSat-Caliop Cloud Top Heights"
    else:
        calipso_displacement = 0

        pixel_position = numpy.arange(len(caObj_calipso.elevation))
        ax.plot(pixel_position, caObj_calipso.elevation, 'k')
    
        for i in range(len(caObj_calipso.elevation)):
            if caObj_calipso.elevation[i] > 0:
                ax.vlines(i, 0 ,caObj_calipso.elevation[i], 'k')
        title = "AVHRR-Caliop Cloud Top Heights"
    #: Plot Caliop 
    
    dummy = caObj_calipso.latitude.shape[0]
    pixel_position = numpy.arange(dummy)
    for i in range(10):
        base_ok = caliop_base[i]
        top_ok = caliop_height[i]
        for idx in range(len(pixel_position)):
            if cal_data_ok[idx]:
#                if caObj_calipso.number_of_layers_found[idx] > i:
                if base_ok[idx] != -9.0:
                    ax.vlines(pixel_position[idx] + calipso_displacement ,base_ok[idx], top_ok[idx],
                                color = "green", linestyle = 'solid', linewidth = 0.5)
                
    
    
    #: Plot Avhrr   
    ax.plot(pixel_position_ok + calipso_displacement, avhrr_ctth_ok, 'b+')#\
            #color = "blue", linestyle = "+")#,cex=0.7)
    filename = "%s/%skm_%s_cloudsat_calipso_avhrr_clouds.%s" \
        %(plotpath, RESOLUTION, basename, file_type)
    ax.set_ylim(0, maxheight)
    ax.set_title(title)
    ax.set_xlabel("Track Position")
    ax.set_ylabel("Cloud Height")
    fig.savefig(filename, format = file_type)

# -----------------------------------------------------
def drawCalClsatCWCAvhrrPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase):#, caObj, CALIPSO_DISPLACED, caliop_base, caliop_height, cal_data_ok, avhrr_ctth_cal_ok):
    from matplotlib import pyplot as plt
    
    # -------------------------------------------------------------------------------------------------------------------------------
    # Create environment for plotting LWP and IWP        pixel_position_plain=numpy.arange(dummy)
    
    if phase=='IW':
        dataP=clsatObj.cloudsatcwc.RVOD_ice_water_path
        dataC=clsatObj.cloudsatcwc.RVOD_ice_water_content
    else:
        dataP=clsatObj.cloudsatcwc.RVOD_liq_water_path
        dataC=clsatObj.cloudsatcwc.RVOD_liq_water_content

    storlek_x=clsatObj.cloudsatcwc.latitude.shape[0] 				# Decides the size of the array latitude
    andrad_storlek_x = int(storlek_x*CLOUDSAT_TRACK_RESOLUTION) 		# Changes the resolution with a factor of 1.076
    #geometric_range_CloudSatcwc = dummycwc
    pixel_position=numpy.arange(storlek_x) 				# Creates an array with numbers 0 -> storlek_x
    andrad_pixel_position=numpy.arange(andrad_storlek_x)	# Creates an array with size geometric_range_CloudSatcwc numbers from 0->size
#    y_zeros=numpy.zeros((storlek_x),'d')		# Creates an array with size andrad_storlek_x Filled with zeros of type double
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
    
    maxvalue=int(max(dataP)+100)							# Findes max value and add 100
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(andrad_pixel_position[::], y, "ko")
    ax.set_ylim(0, maxvalue)
    ax.set_xlabel("Track Position")
    ax.set_ylabel("%sP [g/m^2]" %(phase))
    ax.set_title("CloudSat %sP"%(phase))
    filename = "%s/%ikm_%s_calipso_%sP."%(plotpath, RESOLUTION, basename, phase)
    fig.savefig(filename + 'eps')
    fig.savefig(filename + 'png')
    # ----------------------------------------------------------------------------------------------------------------------------------
    
    #------------------------------------------------------------------------------ 
    # rpy.r.postscript("%s/1km_%s_calipso_%sP.eps"%(plotpath,'test',phase),width=900,height=400)
    #-------------------------------------------------------- yy=numpy.arange(1)
    #-------------------------------------------- ztemp=numpy.asmatrix(y.copy())
    #--------------------------------------------- z=numpy.reshape(ztemp,(-1,1))
    # rpy.r.image(andrad_pixel_position,yy,z,xlab="Track Position",ylab="%sP [g/m^2]"%(phase),main="CloudSat %sP"%(phase))
    #----------------------------------------------------------- rpy.r.dev_off()

    # ----------------------------------------------------------------------------------------------------------------------------------
    # Create environment for plotting LWC and IWC
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    filename = "%s/%ikm_%s_calipso_%sC."%(plotpath, RESOLUTION, basename, phase)
    #max(clsatObj.cloudsatcwc.RVOD_ice_water_content)
    clsatcwc_max_height=numpy.zeros([1,storlek_x],'d')
    y_adjust=andrad_y_zeros
    test=andrad_y_zeros
    for it in range(storlek_x):
        test[it]=int(pixel_position[it]*CLOUDSAT_TRACK_RESOLUTION)

    for i in range(andrad_storlek_x): #Just extend original elevation array
        elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)
        y_adjust[i] = elevationcwc[elevation_index] # Scales elevation
    
    ax1.plot(andrad_pixel_position[::], y_adjust, 'k')

    # Makes the surface black
    for i in range(andrad_storlek_x): #Let's also plot elevation in the vertical
        if y_adjust[i] > 0.0:
            ax1.vlines(i, 0.0, y_adjust[i], color="black", linestyle='solid', linewidth=1)
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
                ax1.vlines(int(pixel_position[idx] * CLOUDSAT_TRACK_RESOLUTION), \
                               base_heightcwc, top_heightcwc, color="red", \
                               linestyle='solid', linewidth=1)
                clsatcwc_max_height[0,idx] = max(clsatcwc_max_height[0,idx],top_heightcwc) 
    
    ax1.set_xlabel("Track Position")
    ax1.set_ylabel("%sC Height"%(phase))
    ax1.set_title("CloudSat-%sC"%(phase))
    if MAXHEIGHT == None:
        maxheight = numpy.max(top_heightcwc) + 1000
    else:
        maxheight = MAXHEIGHT
    ax1.set_ylim(0, maxheight)
    fig1.savefig(filename + 'eps')
    fig1.savefig(filename + 'png')
# -----------------------------------------------------

def drawCalClsatAvhrrPlotTimeDiff(latitude, clsat_diff_sec_1970, cal_diff_sec_1970, plotpath, basename, res, file_type='eps'):
    from matplotlib import pyplot as plt

    # Plot time diff
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Let's use the CloudSat track as our reference
    dummy = latitude.shape[0]
    geometric_range_CloudSat = int(dummy * CLOUDSAT_TRACK_RESOLUTION)
    
#    x=numpy.arange(geometric_range_CloudSat)    
#    y=numpy.ones(geometric_range_CloudSat)*-1000
    
    time_diff_cloudsat=numpy.zeros((geometric_range_CloudSat),'d')
    maxvalue = numpy.nanmax(cal_diff_sec_1970)/60
    minvalue = numpy.nanmin(cal_diff_sec_1970)/60
    title = "Time Difference Between AVHRR and CALIPSO"
    if clsat_diff_sec_1970 != None:
        for i in range(geometric_range_CloudSat): #Just extend original time_diff array
            elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)  
            time_diff_cloudsat[i] = clsat_diff_sec_1970[elevation_index] 
        
        maxvalue = numpy.max([numpy.nanmax(clsat_diff_sec_1970)/60, maxvalue])
        minvalue = numpy.min([numpy.nanmin(clsat_diff_sec_1970)/60, minvalue])
        title = "Time Difference Between AVHRR and CloudSat/CALIPSO"
#        ax.plot(x, y, color="black", linestyle='solid', linewidth=0)

        biggest_Cloudsat_diff = numpy.nanmax(numpy.abs(time_diff_cloudsat))
#        if time_diff_cloudsat.shape[0] >= cal_diff_sec_1970.shape[0]:
#            length_time_diff = time_diff_cloudsat.shape[0]
#        else:
#            length_time_diff = cal_diff_sec_1970.shape[0]
#        rpy.r.legend((length_time_diff-1500),MAXVALUE-5,["CloudSat (max time diff %.2f min)" % round(biggest_Cloudsat_diff/60,2),"Calipso (max time diff %.2f min)" % round(biggest_Calipso_diff/60,2)],pch=[4,4],col=["red","green"])
        ax.plot(time_diff_cloudsat/60,"r+", label = "CloudSat (max time diff %.2f min)" % round(biggest_Cloudsat_diff/60,2))

    biggest_Calipso_diff = numpy.nanmax(numpy.abs(cal_diff_sec_1970))
    ax.set_title(title)    
    ax.set_xlabel("Track Position")
    ax.set_ylabel("Time diff [min]")
    ax.set_ylim(minvalue-10, maxvalue+10)
    ax.plot(cal_diff_sec_1970/60,"g+", label = "CALIPSO (max time diff %.2f min)" % round(biggest_Calipso_diff/60,2))
    ax.legend()
    fig.savefig("%s/%skm_%s_time_diff.%s"%(plotpath,RESOLUTION,basename, file_type))
    
# -----------------------------------------------------   
def drawCalClsatAvhrrPlotSATZ(latitude, AvhrrClsatSatz, AvhrrCalSatz, plotpath, basename, res, file_type='eps'):
    from matplotlib import pyplot as plt

    # Plot Satellite Zenith Angle
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    
    # Let's use the CloudSat track as our reference
    dummy = latitude.shape[0]
    geometric_range_CloudSat = int(dummy*CLOUDSAT_TRACK_RESOLUTION)   
      
#    x=numpy.arange(geometric_range_CloudSat)    
#    y=numpy.ones(geometric_range_CloudSat)*-1000
    maxvalue = numpy.nanmax(AvhrrCalSatz)
    minvalue = numpy.nanmin(AvhrrCalSatz)
    title = format_title("AVHRR compared to CALIPSO")
    if AvhrrClsatSatz != None:
        AvhrrClsatSatz_adjust=numpy.zeros((geometric_range_CloudSat),'d')
        
        for i in range(geometric_range_CloudSat): #Just extend original x array
            elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)        
            AvhrrClsatSatz_adjust[i] = AvhrrClsatSatz[elevation_index] #For satz ploting
        
        maxvalue = numpy.max([numpy.nanmax(AvhrrClsatSatz), maxvalue])
        minvalue = numpy.min([numpy.nanmin(AvhrrClsatSatz), minvalue])
        title = format_title("AVHRR compared to CloudSat/CALIPSO")
#        ax.plot(x, y, color="black", linestyle='solid', linewidth=0)
        ax.plot(AvhrrClsatSatz_adjust, "r+", label = "CloudSat")   

    ax.set_xlabel("Track Position")
    ax.set_ylabel("satellite zenith angle [deg]")
    ax.set_title(title)
    ax.set_ylim(minvalue-10, maxvalue+10)
    ax.plot(AvhrrCalSatz,"g+", label = "CALIPSO")
    ax.legend()
    fig.savefig("%s/%skm_%s_satz.%s"%(plotpath,RESOLUTION,basename, file_type))
    ax.clear()
    fig.clear()
