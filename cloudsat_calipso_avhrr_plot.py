# -*- coding: utf-8 -*-
#Program cloudsat_calipso_avhrr_plot.py

import numpy
from config import MAXHEIGHT, CLOUDSAT_CLOUDY_THR,\
    RESOLUTION, CLOUDSAT_TRACK_RESOLUTION
AZIMUTH_RANGE = [0, 360]
def format_title(title):
    """
    Format *title*, possibly adding satellite azimuth range (depending on range).
    
    Returns formatted title.
    """
    if AZIMUTH_RANGE[0] > 0 or AZIMUTH_RANGE[1] < 90:
        title = "%s (satz range [%.1f, %.1f] deg)" % (title, AZIMUTH_RANGE[0], AZIMUTH_RANGE[1])
    
    return title

# -----------------------------------------------------
def drawCalClsatGEOPROFAvhrrPlot(clsatObj_cloudsat, 
                                 caObj_calipso, 
                                 elevation, data_ok,
                                 CALIPSO_DISPLACED, caliop_base,
                                 caliop_height, cal_data_ok,
                                 avhrr_ctth_cal_ok, plotpath, basename,
                                 mode, emissfilt_calipso_ok=None, file_type='eps',
                                 **options):

    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'

    from matplotlib import pyplot as plt

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
        # Why this!? AD, 2012-Oct
        #if maxheight < 12000.:
        #    maxheight = 12000
        #elif maxheight < 25000.:
        #    maxheight = 25000.
        #else:
        #    maxheight = maxheight + 1000
        maxheight = maxheight + 1000
    else:
        maxheight = MAXHEIGHT

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
        title = "%s-CloudSat-CALIOP Cloud Top Heights" % instrument.upper()
    else:
        calipso_displacement = 0

        pixel_position = numpy.arange(len(caObj_calipso.elevation))
        ax.plot(pixel_position, caObj_calipso.elevation, 'k')
    
        for i in range(len(caObj_calipso.elevation)):
            if caObj_calipso.elevation[i] > 0:
                ax.vlines(i, 0 ,caObj_calipso.elevation[i], 'k')
        title = "%s-CALIOP Cloud Top Heights" % instrument.upper()
    #: Plot Caliop 
    
    dummy = caObj_calipso.latitude.shape[0]
    pixel_position = numpy.arange(dummy)
    caliop_label_set = False
    for i in range(10):
        base_ok = caliop_base[i]
        top_ok = caliop_height[i]
        for idx in range(len(pixel_position)):
            if cal_data_ok[idx]:
#                if caObj_calipso.number_of_layers_found[idx] > i:
                if base_ok[idx] != -9.0:
                    if not caliop_label_set:
                        ax.vlines(pixel_position[idx] + calipso_displacement,
                                  base_ok[idx], top_ok[idx],
                                  color = "green", linestyle = 'solid', linewidth = 0.5,
                                  label='caliop')
                        caliop_label_set = True
                    else:
                        ax.vlines(pixel_position[idx] + calipso_displacement ,base_ok[idx], top_ok[idx],
                                  color = "green", linestyle = 'solid', linewidth = 0.5)
    
    #: Plot Avhrr   
    ax.plot(pixel_position_ok + calipso_displacement, avhrr_ctth_ok, 'b+', 
            label=instrument.upper())
    ax.set_ylim(0, maxheight)
    ax.set_title(title)
    ax.set_xlabel("Track Position")
    ax.set_ylabel("Cloud Height (meter)")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    if isinstance(file_type, str) == True:
        filename = "%s/%skm_%s_cloudsat_calipso_%s_clouds.%s" \
            % (plotpath, RESOLUTION, basename, instrument, file_type)
        fig.savefig(filename, format = file_type)
    else:
        for filetype in file_type:
            filename = "%s/%skm_%s_cloudsat_calipso_%s_clouds.%s" \
                %(plotpath, RESOLUTION, basename, instrument, filetype)
            fig.savefig(filename, format = filetype)
      
# -----------------------------------------------------
def drawCalClsatCWCAvhrrPlot(clsatObj, elevationcwc, data_okcwc, 
                             plotpath, basename, phase, **options):

    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'

    from matplotlib import pyplot as plt
    
    # -----------------------------------------------------------------
    # Create environment for plotting LWP and IWP
    # pixel_position_plain=numpy.arange(dummy)
    
    if phase=='IW':
        dataP=clsatObj.cloudsatcwc.RVOD_ice_water_path
        dataC=clsatObj.cloudsatcwc.RVOD_ice_water_content
    else:
        dataP=clsatObj.cloudsatcwc.RVOD_liq_water_path
        dataC=clsatObj.cloudsatcwc.RVOD_liq_water_content

    # Decides the size of the array latitude:
    storlek_x=clsatObj.cloudsatcwc.latitude.shape[0]
    # Changes the resolution with a factor of 1.076:			
    andrad_storlek_x = int(storlek_x*CLOUDSAT_TRACK_RESOLUTION)
    # Creates an array with numbers 0 -> storlek_x
    pixel_position=numpy.arange(storlek_x)
    # Creates an array with size geometric_range_CloudSatcwc numbers from 0->size
    andrad_pixel_position=numpy.arange(andrad_storlek_x)
    # Creates an array with size andrad_storlek_x Filled with zeros of type double
    andrad_y_zeros=numpy.zeros(andrad_storlek_x,'d')
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
    
    # Findes max value and add 100
    maxvalue=int(max(dataP)+100)							
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
    # -------------------------------------------------------
    
    # Create environment for plotting LWC and IWC
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    filename = "%s/%ikm_%s_calipso_%sC."%(plotpath, RESOLUTION, basename, phase)
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

    # Takes lowest Height
    clsatcwc_max_height[0,:]=numpy.repeat(clsatObj.cloudsatcwc.Height[124,::],data_okcwc)
    for i in range(125):
        heightcwc = clsatObj.cloudsatcwc.Height[i,::]
        RVOD_IWC = dataC[i,::]
        for idx in range(storlek_x):
            nidx = int(RVOD_IWC[idx]+0.5)
            if nidx == 0:	# If there are no value dont do anything
                continue
            # Treschold when measurements are valid
            if heightcwc[idx] < 240*4 or heightcwc[idx] > MAXHEIGHT:
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

    return

# -----------------------------------------------------
def drawCalAvhrrTime(cal_sec_1970, avhrr_sec_1970, **options):
    from matplotlib import pyplot as plt

    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    maxvalue = numpy.nanmax([numpy.nanmax(cal_sec_1970), numpy.nanmax(avhrr_sec_1970)]) / 60.
    minvalue = numpy.nanmin([numpy.nanmin(cal_sec_1970), numpy.nanmin(avhrr_sec_1970)]) / 60.
    title = "Time since 1970"
    
    ax.set_title(title)
    ax.set_xlabel("Track Position")
    ax.set_ylabel("Time [min]")
    ax.set_ylim(minvalue-10, maxvalue+10)
    ax.plot(cal_sec_1970 / 60.,"g+", label = "CALIPSO")
    ax.plot(avhrr_sec_1970 / 60.,"r+", label = instrument.upper())
    ax.legend()
    fig.show()

# -----------------------------------------------------
def drawCalClsatAvhrrPlotTimeDiff(latitude, 
                                  clsat_diff_sec_1970, 
                                  cal_diff_sec_1970, 
                                  plotpath, basename, 
                                  resolution, file_type='eps',
                                  **options):
    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'

    from matplotlib import pyplot as plt

    # Plot time diff
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Let's use the CloudSat track as our reference
    dummy = latitude.shape[0]
    geometric_range_CloudSat = int(dummy * CLOUDSAT_TRACK_RESOLUTION)
    
    time_diff_cloudsat=numpy.zeros((geometric_range_CloudSat),'d')
    maxvalue = numpy.nanmax(cal_diff_sec_1970)/60
    minvalue = numpy.nanmin(cal_diff_sec_1970)/60
    title = "Time Difference CALIPSO - %s" % instrument.upper()
    if clsat_diff_sec_1970 != None:
        for i in range(geometric_range_CloudSat): #Just extend original time_diff array
            elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)  
            time_diff_cloudsat[i] = clsat_diff_sec_1970[elevation_index] 
        
        maxvalue = numpy.max([numpy.nanmax(clsat_diff_sec_1970)/60, maxvalue])
        minvalue = numpy.min([numpy.nanmin(clsat_diff_sec_1970)/60, minvalue])
        title = "Time Difference Between %s and CloudSat/CALIPSO" % instrument.upper()

        biggest_Cloudsat_diff = numpy.nanmax(numpy.abs(time_diff_cloudsat))
        ax.plot(time_diff_cloudsat/60,"r+", 
                label = "CloudSat (max time diff %.2f min)" % round(biggest_Cloudsat_diff/60,2))

    biggest_Calipso_diff = numpy.nanmax(numpy.abs(cal_diff_sec_1970))
    ax.set_title(title)    
    ax.set_xlabel("Track Position")
    ax.set_ylabel("Time diff (CALIPSO - %s)[min]" % instrument.upper())
    ax.set_ylim(minvalue-10, maxvalue+10)
    ax.plot(cal_diff_sec_1970/60,"g+", label = "CALIPSO (max time diff %.2f min)" % round(biggest_Calipso_diff/60,2))
    ax.legend()
    
    if isinstance(file_type, str) == True:
        fig.savefig("%s/%skm_%s_time_diff.%s" % (plotpath, 
                                                 resolution, basename, file_type))
    else:
        for filetype in file_type:
            fig.savefig("%s/%skm_%s_time_diff.%s" % (plotpath,
                                                     resolution, basename, filetype))

    return
            
# -----------------------------------------------------   
def drawCalClsatAvhrrPlotSATZ(latitude, 
                              AvhrrClsatSatz, 
                              AvhrrCalSatz, 
                              plotpath, basename, 
                              resolution, 
                              file_type='eps', 
                              **options):

    from matplotlib import pyplot as plt

    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'

    # Plot Satellite Zenith Angle
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    
    # Let's use the CloudSat track as our reference
    dummy = latitude.shape[0]
    geometric_range_CloudSat = int(dummy*CLOUDSAT_TRACK_RESOLUTION)   
      
    maxvalue = numpy.nanmax(AvhrrCalSatz)
    minvalue = numpy.nanmin(AvhrrCalSatz)
    title = format_title("%s compared to CALIPSO" % instrument.upper())
    if AvhrrClsatSatz != None:
        AvhrrClsatSatz_adjust=numpy.zeros((geometric_range_CloudSat),'d')
        
        for i in range(geometric_range_CloudSat): #Just extend original x array
            elevation_index = int(i/CLOUDSAT_TRACK_RESOLUTION + 0.5)        
            AvhrrClsatSatz_adjust[i] = AvhrrClsatSatz[elevation_index] #For satz ploting
        
        maxvalue = numpy.max([numpy.nanmax(AvhrrClsatSatz), maxvalue])
        minvalue = numpy.min([numpy.nanmin(AvhrrClsatSatz), minvalue])
        title = format_title("%s compared to CloudSat/CALIPSO" % instrument.upper())
        ax.plot(AvhrrClsatSatz_adjust, "r+", label = "%s - CloudSat" % instrument.upper())

    ax.set_xlabel("Track Position")
    ax.set_ylabel("satellite zenith angle [deg]")
    ax.set_title(title)
    ax.set_ylim(minvalue-10, maxvalue+10)
    ax.plot(AvhrrCalSatz,"g+", label = "%s - CALIPSO" % instrument.upper())
    ax.legend()
    if isinstance(file_type, str) == True:
        fig.savefig("%s/%skm_%s_satz.%s" % (plotpath, resolution, basename, file_type))
    else:
        for filetype in file_type:
            fig.savefig("%s/%skm_%s_satz.%s" % (plotpath, resolution, basename, filetype))

    del ax
    del fig
    del plt

    return


def map_avhrr_track(avhrr_lonlat, track_lonlat):
    """
    Plot *avhrr_lonlat* and *track_lonlat* on global map and return the figure.
    
    """
    from matplotlib.pyplot import figure
    from mpl_toolkits.basemap import Basemap
    
    fig = figure()
    ax = fig.add_subplot(111)
    m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-90,
                urcrnrlon=180, urcrnrlat=90, resolution='l', ax=ax)
    
    m.drawcoastlines(linewidth=.5, color='grey')
    
    # Don't draw each pixel, or the machine will choke!
    npixels = avhrr_lonlat[0].size
    from math import sqrt
    step = int(round(sqrt(npixels / 1e5))) # Will give a total of about 1e5 pixels
    _slice_2d = (slice(None, None, step),) * 2
    m.pcolormesh(avhrr_lonlat[0][_slice_2d], avhrr_lonlat[1][_slice_2d],
                 avhrr_lonlat[1][_slice_2d], alpha=.5)
    m.plot(track_lonlat[0], track_lonlat[1], 'o', markersize=1, alpha=.1,
           label='track')
    
    return fig
