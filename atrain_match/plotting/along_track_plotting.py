# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
# -*- coding: utf-8 -*-
# Program truth_imager_plot.py

import numpy as np
from atrain_match.config import RESOLUTION

from matplotlib import pyplot as plt

# -----------------------------------------------------
def plot_cal_clsat_geoprof_imager(match_clsat,
                                 match_calipso,
                                 imager_ctth_m_above_seasurface,
                                 plotpath,
                                 basename,
                                 mode,
                                 file_type='png',
                                 **options):
    instrument = 'imager'
    if 'instrument' in options:
        instrument = options['instrument']
    MAXHEIGHT = None
    if 'MAXHEIGHT' in options:
        MAXHEIGHT = options["MAXHEIGHT"]
    caliop_height = match_calipso.calipso.layer_top_altitude*1000
    caliop_base = match_calipso.calipso.layer_base_altitude*1000
    calipso_val_h = match_calipso.calipso.validation_height
    caliop_base[caliop_base<0]=-9
    caliop_height[caliop_height<0]=-9
    pixel_position = np.arange(match_calipso.calipso.latitude.shape[0])
    # Calculates Hihest Cloud Top
    if MAXHEIGHT is None:
        maxheight_calipso = np.nanmax(caliop_height)
        maxheight_imager = np.nanmax(imager_ctth_m_above_seasurface)
        max_height_sat = np.max([maxheight_calipso, maxheight_imager])
        maxheight = max_height_sat + 1000
    else:
        maxheight =  MAXHEIGHT
    # PLOT
    fig = plt.figure()
    # Plot ground
    ax = fig.add_subplot(111)
    ax.vlines(pixel_position, 0, match_calipso.calipso.elevation,
              color='k', alpha=1.0)
    # plot cloudsat if we have it
    if match_clsat is not None:
        #  Colors
        colors=[]
        for i in range(10):
            colors.append(np.divide([(40-i)*6, (40-i)*6, (40-i)*6], 255.))
        for i in range(10,30):
            colors.append(np.divide([100+(i-10)*5, 100+(i-10)*5, 100-(i-10)*5], 255.))
        for i in range(30,100):
            colors.append(np.divide([255, 50, 50], 255.))

        # Plot CloudSat
        base_height = match_clsat.cloudsat.validation_height_base
        top_height = match_clsat.cloudsat.validation_height
        plot_these = top_height>0
        ax.vlines(pixel_position[match_clsat.cloudsat.calipso_index[plot_these]],
                  base_height[plot_these], top_height[plot_these], 'm',# color = colors[nidx], \
                  linestyle = 'solid', linewidth = 1,  rasterized=True, label='CPR (CloudSat)')
        title = "%s-CloudSat-CALIOP Cloud Top Heights" % instrument.upper()
    else:
        title = "%s-CALIOP Cloud Top Heights" % instrument.upper()
    #  Plot Caliop
    caliop_label_set = False
    # Plot all 10 calipso layers
    for i in range(10):
        base_ok = caliop_base[:,i]
        top_ok  = caliop_height[:,i]
        all_thin = calipso_val_h<=base_ok
        all_thick = np.logical_or(calipso_val_h>=top_ok, calipso_val_h<=0)
        half_thin = np.logical_and(calipso_val_h>base_ok,
                                  calipso_val_h<top_ok)
        all_thin = np.logical_and(all_thin,top_ok>0 )
        all_thick = np.logical_and(all_thick,top_ok>0 )
        half_thin = np.logical_and(half_thin,top_ok>0 )

        if np.min(top_ok<0):
            # no more clouds, quit plotting calipso
            break
        if caliop_label_set:
            ax.vlines(pixel_position[all_thick],
                      base_ok[all_thick], top_ok[all_thick], linewidth=0.5,
                      colors="g", linestyle='solid',
                      alpha=0.5 ,  rasterized=True)
        else:
            ax.vlines(pixel_position[all_thick],
                      base_ok[all_thick], top_ok[all_thick],
                      colors="g", linewidth=0.5, linestyle='solid',
                      alpha=0.5, label='caliop',  rasterized=True)
            caliop_label_set = True
        ax.vlines(pixel_position[all_thin],
                  base_ok[all_thin], top_ok[all_thin], linewidth=0.5,
                  colors="y", linestyle='solid',
                  alpha=0.5,  rasterized=True)
        ax.vlines(pixel_position[half_thin],
                  calipso_val_h[half_thin], top_ok[half_thin], linewidth=0.5,
                  colors="y", linestyle='solid',
                  alpha=0.5,  rasterized=True)
        # ax.vlines(pixel_position[half_thin],
        #          base_ok[half_thin], calipso_val_h[half_thin], linewidth=0.5,
        #          colors="g", linestyle='solid',
        #          alpha=0.5,  rasterized=True)

    #  Plot Imager
    got_height = imager_ctth_m_above_seasurface>=0
    ax.plot(pixel_position[got_height], imager_ctth_m_above_seasurface[got_height], 'b+',
            label=instrument.upper(),  rasterized=True)
    ax.set_ylim(0, maxheight)
    # plt.show()
    ax.set_title(title)
    ax.set_xlabel("Track Position")
    ax.set_ylabel("Cloud Height (meter)")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    if isinstance(file_type, str) == True:
        filename = "%s/%skm_%s_cloudsat_calipso_%s_clouds_%s.%s" \
            % (plotpath, RESOLUTION, basename, instrument, mode.lower(), file_type)
        fig.savefig(filename, format = file_type)
    else:
        for filetype in file_type:
            filename = "%s/%skm_%s_cloudsat_calipso_%s_clouds_%s.%s" \
                %(plotpath, RESOLUTION, basename, instrument,mode.lower(), filetype)
            fig.savefig(filename, format = filetype)

# added plot with two pps cloud-heights and no cloudsat
def drawCalPPSHeightPlot_PrototypePPSHeight(match_calipso_calipso,
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
    instrument = 'imager'
    if 'instrument' in options:
        instrument = options['instrument']
    MAXHEIGHT = None
    if 'MAXHEIGHT' in options:
        MAXHEIGHT = options["MAXHEIGHT"]


    # Prepare for Imager
    caliop_height = match_calipso_calipso.layer_top_altitude*1000
    caliop_base = match_calipso_calipso.layer_base_altitude*1000
    caliop_base[caliop_base<0]=-9
    caliop_height[caliop_height<0]=-9
    pixel_position=np.arange(match_calipso_calipso.latitude.shape[0])
    pixel_position_ok = pixel_position[data_ok]
    imager_ctth_ok1 = ctth_height1[data_ok]
    imager_ctth_ok2 = ctth_height2[data_ok]
#    # Calculates Hihest Cloud Top

    if MAXHEIGHT is None:
        maxheight_calipso = np.nanmax(caliop_height)
        maxheight_imager = np.nanmax(ctth_height1)
        max_height_sat = np.max([maxheight_calipso, maxheight_imager])
        maxheight = max_height_sat + 1000
    else:
        maxheight =  MAXHEIGHT
    fig = plt.figure(figsize = (20,15))
    title = "%s-CALIOP Cloud Top Heights" % instrument.upper()
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    for subplot_nr in [211, 212]:
        ax = fig.add_subplot(subplot_nr)

        ax.vlines(pixel_position, 0, match_calipso_calipso.elevation,
                  color='k', alpha=1.0)
        #  Plot Caliop
        caliop_label_set = False
        # Plot all 10 calipso layers
        for i in range(10):
            base_ok = caliop_base[:,0]
            top_ok  = caliop_height[:,0]
            if np.min(top_ok<0):
                # no more clouds, quit plotting calipso
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
    #  Plot Imager
    ax = fig.add_subplot(211)
    ax.plot(pixel_position_ok, imager_ctth_ok1, 'b+', linewidth=0.5,
            label=instrument.upper() + " old-CTTH")
    false = (caliop_height[:,0]<0)[data_ok]
    ax.plot(pixel_position_ok[false], imager_ctth_ok1[false], 'r+', linewidth=0.5,
            label=instrument.upper() + " old-CTTH (false")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, maxheight)
    ax = fig.add_subplot(212)
    ax.plot(pixel_position_ok, imager_ctth_ok2, 'c+', linewidth=0.5,
            label=instrument.upper() + " nn-CTTH")
    ax.plot(pixel_position_ok[false], imager_ctth_ok2[false], 'm+', linewidth=0.5,
            label=instrument.upper() + " nn-CTTH (false)")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(0, maxheight)
    plt.suptitle(title, fontsize=24)
    ax.set_xlabel("Track Position", fontsize=22)

    # plt.show()
    if isinstance(file_type, str) == True:
        filename = "%s/%skm_%s_calipso_%s_clouds.%s" \
            % (plotpath, RESOLUTION, basename, instrument, file_type)
        fig.savefig(filename, format = file_type)
    else:
        for filetype in file_type:
            filename = "%s/%skm_%s_calipso_%s_clouds.%s" \
                %(plotpath, RESOLUTION, basename, instrument, filetype)
            fig.savefig(filename, format = filetype)
    # plt.show()
    plt.close("all")

# -----------------------------------------------------
def plot_cal_clsat_cwc_imager(match_clsat, elevationcwc, data_okcwc,
                             plotpath, basename, phase, **options):

    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'imager'
    if phase=='IW':
        y1=match_clsat.cloudsat.RVOD_ice_water_path
        # dataC=match_clsat.cloudsat.RVOD_ice_water_content
        y2 = match_clsat.imager.cpp_iwp
    else:
        y1=match_clsat.cloudsat.RVOD_liq_water_path
        # dataC=match_clsat.cloudsat.RVOD_liq_water_content
        y2 = match_clsat.imager.cpp_lwp
    pixel_position=np.arange(match_clsat.cloudsat.latitude.shape[0])
    use = np.logical_and(y1>=0, y2>=0)
    # Findes max value and add 100
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(pixel_position[use], y1[use], "k." ,label='RVOD')
    ax.plot(pixel_position[use], y2[use], "b.", label='PPS')
    # ax.set_ylim(0, maxvalue)
    ax.set_xlabel("Track Position")
    ax.set_ylabel("%sP [g/m^2]" %(phase))
    ax.set_title("CloudSat %sP"%(phase))
    filename = "%s/%ikm_%s_cloudsat_rvod_%sP."%(plotpath, RESOLUTION, basename, phase)
    fig.savefig(filename + 'png')
    # -------------------------------------------------------
    return

# -----------------------------------------------------
def plot_cal_clsat_imager_time_diff(match_clsat,
                                  match_calipso,
                                  plotpath, basename,
                                  resolution, file_type='png',
                                  **options):
    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'imager'
    pixel_position=np.arange(match_calipso.calipso.latitude.shape[0])
    cal_diff_sec_1970 = match_calipso.diff_sec_1970/60.0
    # Plot time diff
    fig = plt.figure()
    ax = fig.add_subplot(111)
    maxvalue = np.nanmax(cal_diff_sec_1970)/60.0
    minvalue = np.nanmin(cal_diff_sec_1970)/60.0
    title = "Time Difference CALIPSO - %s" % instrument.upper()
    ylabel_str = "Time diff (CALIPSO - %s)[min]" % instrument.upper()
    if match_clsat is not None:
        clsat_diff_sec_1970 = match_clsat.diff_sec_1970/60.0
        title = "Time Difference Between %s and CloudSat/CALIPSO" % instrument.upper()
        ylabel_str = "Time diff (CALIPSO/CloudSat - %s)[min]" % instrument.upper()
        maxvalue = np.max([np.nanmax(clsat_diff_sec_1970)/60.0, maxvalue])
        minvalue = np.min([np.nanmin(clsat_diff_sec_1970)/60.0, minvalue])

        biggest_Cloudsat_diff = np.nanmax(np.abs(clsat_diff_sec_1970/60.0))
        ax.plot(pixel_position[match_clsat.cloudsat.calipso_index],clsat_diff_sec_1970/60.0, 'r+',
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
def plot_cal_clsat_imager_satz(match_clsat,
                              match_calipso,
                              plotpath, basename,
                              resolution,
                              file_type='eps',
                              **options):

    if 'instrument' in options:
        instrument = options['instrument']
    else:
        instrument = 'imager'
    pixel_position = np.arange(match_calipso.calipso.latitude.shape[0])
    # Plot Satellite Zenith Angle
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if match_clsat is not None:
        ax.plot(pixel_position[match_clsat.cloudsat.calipso_index],
                match_clsat.imager.satz,
                'r+',
                label = "%s satz - CloudSat" % instrument.upper())
    ax.set_xlabel("Track Position")
    ax.set_ylabel("satellite zenith angle [deg]")
    ax.set_title("%s SATZ" % instrument.upper())
    ax.plot(match_calipso.imager.satz,"g", label = "%s - CALIPSO" % instrument.upper())
    ax.legend(numpoints=4)
    if isinstance(file_type, str) == True:
        fig.savefig("%s/%skm_%s_satz.%s" % (plotpath, resolution, basename, file_type))
    else:
        for filetype in file_type:
            fig.savefig("%s/%skm_%s_satz.%s" % (plotpath, resolution, basename, filetype))
    return


def map_imager_track(imager_lonlat, track_lonlat):
    """
    Plot *imager_lonlat* and *track_lonlat* on global map and return the figure.

    """
    from mpl_toolkits.basemap import Basemap

    fig = figure()
    ax = fig.add_subplot(111)
    m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-90,
                urcrnrlon=180, urcrnrlat=90, resolution='l', ax=ax)

    m.drawcoastlines(linewidth=.5, color='grey')

    # Don't draw each pixel, or the machine will choke!
    npixels = imager_lonlat[0].size
    from math import sqrt
    step = int(round(sqrt(npixels / 1e5))) # Will give a total of about 1e5 pixels
    _slice_2d = (slice(None, None, step),) * 2
    m.pcolormesh(imager_lonlat[0][_slice_2d], imager_lonlat[1][_slice_2d],
                 imager_lonlat[1][_slice_2d], alpha=.5)
    m.plot(track_lonlat[0], track_lonlat[1], 'o', markersize=1, alpha=.1,
           label='track')

    return fig
