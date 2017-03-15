# -*- coding: utf-8 -*-
#Program cloudsat_calipso_avhrr_plot.py

import numpy as np
from config import MAXHEIGHT, CLOUDSAT_CLOUDY_THR,\
    RESOLUTION

# -----------------------------------------------------
def drawCalClsatGEOPROFAvhrrPlot(clsatObj, 
                                 caObj,  
                                 avhrr_ctth_cal_ok, 
                                 plotpath, 
                                 basename,
                                 mode, 
                                 file_type='png',
                                 **options):

    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'

    from matplotlib import pyplot as plt

    caliop_height = caObj.calipso.layer_top_altitude*1000
    caliop_base = caObj.calipso.layer_base_altitude*1000
    caliop_base[caliop_base<0]=-9
    caliop_height[caliop_height<0]=-9
    pixel_position = np.arange(caObj.calipso.latitude.shape[0])                            
    # Calculates Hihest Cloud Top   
    if MAXHEIGHT == None:
        maxheight_calipso = np.nanmax(caliop_height)
        maxheight_avhrr = np.nanmax(avhrr_ctth_cal_ok)
        max_height_sat = [maxheight_calipso, maxheight_avhrr]
        maxheight = maxheight + 1000
    else:
        maxheight = MAXHEIGHT
    #PLOT    
    fig = plt.figure()
    #Plot ground
    ax = fig.add_subplot(111)
    ax.vlines(pixel_position, 0, caObj.calipso.elevation, 
              color='k', alpha=1.0)
    #plot cloudsat if we have it    
    if clsatObj != None:
        #: Colors
        colors=[]
        for i in range(10):
            colors.append(np.divide([(40-i)*6, (40-i)*6, (40-i)*6], 255.))
        for i in range(10,30):
            colors.append(np.divide([100+(i-10)*5, 100+(i-10)*5, 100-(i-10)*5], 255.))
        for i in range(30,100):
            colors.append(np.divide([255, 50, 50], 255.))
        
        # Plot CloudSat
        clsat_max_height = -9 + np.zeros(clsatObj.cloudsat.latitude.shape) 
        for i in range(125):
            height = clsatObj.cloudsat.Height[:,i]
            cmask_ok = clsatObj.cloudsat.cloud_mask[:,i] > CLOUDSAT_CLOUDY_THR
            base_height = height-120
            top_height = height+120
            plot_these =np.logical_and(cmask_ok, height>240*4)
            #TODO Fix colors!
            ax.vlines(pixel_position[clsatObj.cloudsat.calipso_index[plot_these]],
                      base_height[plot_these], top_height[plot_these], 'm',#color = colors[nidx], \
                      linestyle = 'solid', linewidth = 1)    
            clsat_max_height[clsat_max_height<top_height] = top_height[clsat_max_height<top_height]
        title = "%s-CloudSat-CALIOP Cloud Top Heights" % instrument.upper()
    else:
        title = "%s-CALIOP Cloud Top Heights" % instrument.upper()
    #: Plot Caliop 
    caliop_label_set = False
    #Plot all 10 calipso layers
    for i in range(10):
        base_ok = caliop_base[:,0]
        top_ok  = caliop_height[:,0]
        if np.min(top_ok<0):
            #no more clouds, quit plotting calipso
            break 
        if caliop_label_set:
            ax.vlines(pixel_position, 
                      base_ok, top_ok, linewidth=0.5,  
                      colors="g", linestyle='solid', 
                      alpha=1.0 )  
        else:
            ax.vlines(pixel_position, 
                      base_ok, top_ok,
                      colors="g", linewidth=0.5, linestyle='solid', 
                      alpha=1.0, label='caliop')
            caliop_label_set = True
    
    #: Plot Imager   
    got_height = avhrr_ctth_cal_ok>=0
    ax.plot(pixel_position[got_height], avhrr_ctth_cal_ok[got_height], 'b+', 
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

#added plot with two pps cloud-heights and no cloudsat
def drawCalPPSHeightPlot_PrototypePPSHeight(caObj_calipso, 
                                            data_ok,
                                            ctth_height1,
                                            ctth_height2,
                                            plotpath, 
                                            basename,
                                            file_type='png',
                                            xmin=0,
                                            xmax=-1,
                                            **options):
    if xmax<0:
        xmax = len(data_ok)
    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'
    from matplotlib import pyplot as plt
    # Prepare for Avhrr
    caliop_height = caObj_calipso.layer_top_altitude*1000
    caliop_base = caObj_calipso.layer_base_altitude*1000
    caliop_base[caliop_base<0]=-9
    caliop_height[caliop_height<0]=-9
    pixel_position=np.arange(caObj_calipso.latitude.shape[0])
    pixel_position_ok = pixel_position[data_ok]
    avhrr_ctth_ok1 = ctth_height1[data_ok]
    avhrr_ctth_ok2 = ctth_height2[data_ok]  
#    # Calculates Hihest Cloud Top  

    if MAXHEIGHT == None:
        maxheight_calipso = np.nanmax(caliop_height)
        maxheight_avhrr = np.nanmax(ctth_height1)
        max_height_sat = [maxheight_calipso, maxheight_avhrr]
        maxheight = maxheight + 1000
    else:
        maxheight = MAXHEIGHT
    maxheight = MAXHEIGHT
    fig = plt.figure(figsize = (20,15))
    title = "%s-CALIOP Cloud Top Heights" % instrument.upper()
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    for subplot_nr in [211, 212]:
        ax = fig.add_subplot(subplot_nr)

        ax.vlines(pixel_position, 0, caObj_calipso.elevation, 
                  color='k', alpha=1.0)
        #: Plot Caliop 
        caliop_label_set = False
        #Plot all 10 calipso layers
        for i in range(10):
            base_ok = caliop_base[:,0]
            top_ok  = caliop_height[:,0]
            if np.min(top_ok<0):
                #no more clouds, quit plotting calipso
                break 
            if caliop_label_set:
                ax.vlines(pixel_position, 
                          base_ok, top_ok, linewidth=0.5,  
                          colors="g", linestyle='solid', 
                          alpha=1.0 )  
            else:
                ax.vlines(pixel_position, 
                          base_ok, top_ok,
                          colors="g", linewidth=0.5, linestyle='solid', 
                          alpha=1.0, label='caliop')
                caliop_label_set = True
        ax.set_ylabel("Cloud Height (meter)", fontsize=22)
    #: Plot Avhrr   
    ax = fig.add_subplot(211)
    ax.plot(pixel_position_ok, avhrr_ctth_ok1, 'b+', linewidth=0.5, 
            label=instrument.upper() + " old-CTTH")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(xmin, xmax)
    ax = fig.add_subplot(212)
    ax.plot(pixel_position_ok, avhrr_ctth_ok2, 'b+', linewidth=0.5,  
            label=instrument.upper() + " nn-CTTH")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(xmin,xmax)
    #ax.set_ylim(0, maxheight)
    plt.suptitle(title, fontsize=24)
    ax.set_xlabel("Track Position", fontsize=22)

    #plt.show() 
    if isinstance(file_type, str) == True:
        filename = "%s/%skm_%s_calipso_%s_clouds.%s" \
            % (plotpath, RESOLUTION, basename, instrument, file_type)
        fig.savefig(filename, format = file_type)
    else:
        for filetype in file_type:
            filename = "%s/%skm_%s_calipso_%s_clouds.%s" \
                %(plotpath, RESOLUTION, basename, instrument, filetype)
            fig.savefig(filename, format = filetype)
      
# -----------------------------------------------------
def drawCalClsatCWCAvhrrPlot(clsatObj, elevationcwc, data_okcwc, 
                             plotpath, basename, phase, **options):

    CLOUDSAT_TRACK_RESOLUTION = 1.076
    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'

    from matplotlib import pyplot as plt
    
    # -----------------------------------------------------------------
    # Create environment for plotting LWP and IWP
    # pixel_position_plain=np.arange(dummy)
    
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
    pixel_position=np.arange(storlek_x)
    # Creates an array with size geometric_range_CloudSatcwc numbers from 0->size
    andrad_pixel_position=np.arange(andrad_storlek_x)
    # Creates an array with size andrad_storlek_x Filled with zeros of type double
    andrad_y_zeros=np.zeros(andrad_storlek_x,'d')
    andrad_dataP=np.ones(andrad_storlek_x)*-1
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
    clsatcwc_max_height=np.zeros([1,storlek_x],'d')
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
    clsatcwc_max_height[0,:]=np.repeat(clsatObj.cloudsatcwc.Height[124,::],data_okcwc)
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
        maxheight = np.max(top_heightcwc) + 1000
    else:
        maxheight = MAXHEIGHT
    ax1.set_ylim(0, maxheight)
    fig1.savefig(filename + 'eps')
    fig1.savefig(filename + 'png')

    return

# -----------------------------------------------------
def drawCalClsatAvhrrPlotTimeDiff(clsatObj, 
                                  caObj, 
                                  plotpath, basename, 
                                  resolution, file_type='png',
                                  **options):
    from matplotlib import pyplot as plt
    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'
    pixel_position=np.arange(caObj.calipso.latitude.shape[0])
    cal_diff_sec_1970 = caObj.diff_sec_1970/60.0
    # Plot time diff
    fig = plt.figure()
    ax = fig.add_subplot(111)
    maxvalue = np.nanmax(cal_diff_sec_1970)/60.0
    minvalue = np.nanmin(cal_diff_sec_1970)/60.0
    title = "Time Difference CALIPSO - %s" % instrument.upper()
    ylabel_str = "Time diff (CALIPSO - %s)[min]" % instrument.upper()
    if clsatObj != None:
        clsat_diff_sec_1970 = clsatObj.diff_sec_1970/60.0
        title = "Time Difference Between %s and CloudSat/CALIPSO" % instrument.upper()
        ylabel_str = "Time diff (CALIPSO/CloudSat - %s)[min]" % instrument.upper()
        maxvalue = np.max([np.nanmax(clsat_diff_sec_1970)/60.0, maxvalue])
        minvalue = np.min([np.nanmin(clsat_diff_sec_1970)/60.0, minvalue])

        biggest_Cloudsat_diff = np.nanmax(np.abs(clsat_diff_sec_1970/60.0))
        ax.plot(pixel_position[clsatObj.cloudsat.calipso_index],clsat_diff_sec_1970/60.0, 'r+', 
                label = "CloudSat (max time diff %.2f min)"%(biggest_Cloudsat_diff)) 
        
    biggest_Calipso_diff = np.nanmax(np.abs(cal_diff_sec_1970/60.0))
    ax.set_title(title)    
    ax.set_xlabel("Track Position")
    ax.set_ylabel(ylabel_str)
    ax.set_ylim(minvalue-5, maxvalue+5)
    ax.plot(pixel_position, cal_diff_sec_1970/60.0,"g", 
            label = "CALIPSO (max time diff %.2f min)" %(biggest_Calipso_diff))
    ax.legend(numpoints=4)
    
    if isinstance(file_type, str) == True:
        fig.savefig("%s/%skm_%s_time_diff.%s" % (plotpath, 
                                                 resolution, basename, file_type))
    else:
        for filetype in file_type:
            fig.savefig("%s/%skm_%s_time_diff.%s" % (plotpath,
                                                     resolution, basename, filetype))
    return
            
# -----------------------------------------------------   
def drawCalClsatAvhrrPlotSATZ(clsatObj, 
                              caObj, 
                              plotpath, basename, 
                              resolution, 
                              file_type='eps', 
                              **options):

    from matplotlib import pyplot as plt
    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'avhrr'
    pixel_position = np.arange(caObj.calipso.latitude.shape[0])
    # Plot Satellite Zenith Angle
    fig = plt.figure() 
    ax = fig.add_subplot(111) 
    if clsatObj != None:
        ax.plot(pixel_position[clsatObj.cloudsat.calipso_index],
                clsatObj.avhrr.satz,
                'r+', 
                label = "%s satz - CloudSat" % instrument.upper())
    ax.set_xlabel("Track Position")
    ax.set_ylabel("satellite zenith angle [deg]")
    ax.set_title("%s SATZ" % instrument.upper())
    ax.plot(caObj.avhrr.satz,"g", label = "%s - CALIPSO" % instrument.upper())
    ax.legend(numpoints=4)
    if isinstance(file_type, str) == True:
        fig.savefig("%s/%skm_%s_satz.%s" % (plotpath, resolution, basename, file_type))
    else:
        for filetype in file_type:
            fig.savefig("%s/%skm_%s_satz.%s" % (plotpath, resolution, basename, filetype))
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
