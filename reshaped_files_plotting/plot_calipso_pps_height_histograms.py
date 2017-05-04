#!/usr/bin/env python
# -*- coding: utf-8 -*-



"""Read all matched data and make some plotting of thresholds
"""

from glob import glob
#import os.path
import os
import numpy as np
from scipy import histogram
import re
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
import matplotlib.pyplot as plt

isACPGv2012=False
isGAC_v2014_morning_sat = True
isGAC_v2014 = True
if isGAC_v2014_morning_sat:
    num_files_to_read = 30*3
    isGAC=True
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa17/5km/20??/06/*/*h5")
    files = files + glob(ROOT_DIR + "metop*/5km/20??/06/*h5")
    figure_name = "figure_morning_sat_"
elif isGAC_v2014:
    num_files_to_read = 30
    isGAC=True
    figure_name = "figure_"
    ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa18/5km/20??/??/*/*noaa*h5")
    files = files + glob(ROOT_DIR + "noaa19/5km/20??/??/*/*noaa*h5")

def  my_plot_histogram_hexbin(x, y, 
                              bins='log', 
                              gridsize=None,
                              ttitle="", txlabel='x', tylabel='y', ymin=None,
                              file_name=None):
    """Plot a 2d-histogram using hexbins
    """
    fig = plt.figure(figsize = (9,8))
    ax = fig.add_subplot(111)
    xdata = x.ravel()
    ydata = y.ravel()
    if gridsize is not None:
        plt.hexbin(xdata,ydata,cmap=plt.cm.gray_r,bins='log',gridsize=gridsize)
    else:
        plt.hexbin(xdata,ydata,cmap=plt.cm.gray_r,bins='log')    
    plt.colorbar()
    plt.title(ttitle, color='k')
    if ymin is not None:
        ax.set_ylim(ymin,np.max(ydata))
    plt.xlabel(txlabel)
    plt.ylabel(tylabel)
    if file_name is not None:
        plt.savefig(file_name)
    plt.show()
    return 1

caObj = CalipsoAvhrrTrackObject()
for filename in files:
    print os.path.basename(filename)
    caObj = caObj + readCaliopAvhrrMatchObj(filename)
    
isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                             caObj.avhrr.all_arrays['cloudtype']<21)
isCloudyCaliop = caObj.calipso.all_arrays['number_layers_found']>0
isCloudyCaliop = np.logical_and(isCloudyCaliop, 
                                 caObj.calipso.all_arrays['total_optical_depth_5km']>0.5)

cal_temp = caObj.calipso.all_arrays['midlayer_temperature'][:,0]+273.15 
cal_height = 1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] 
pps_temp = caObj.avhrr.all_arrays['ctth_temperature']  
pps_height = caObj.avhrr.all_arrays['ctth_height']  
bt11 = caObj.avhrr.all_arrays['bt11micron']
tsur = caObj.avhrr.all_arrays['surftemp']
ctype = caObj.avhrr.all_arrays['cloudtype']

isCloudy =  np.logical_and(isCloudyPPS,isCloudyCaliop)
isCloudy =  np.logical_and(pps_temp>-9,isCloudy)
isCloudy = np.logical_and(pps_temp>150,isCloudy)
isCloudy = np.logical_and(cal_temp>150,isCloudy)


"""
for ind in range(5,15):
    print ind
    isCloudytypeT = np.logical_and(isCloudy,ctype==ind)
    my_plot_histogram_hexbin(bt11[isCloudytypeT], np.abs(pps_temp[isCloudytypeT]-cal_temp[isCloudytypeT]),gridsize=30,
                               txlabel='bt11', tylabel='diff cal pps temp',file_name='misc_png/bt11_cal_minus_pps_temp_klass_%d.png'%(ind))
"""
my_plot_histogram_hexbin(cal_temp[isCloudy], pps_temp[isCloudy],gridsize=30, 
                         txlabel='calipso temp', tylabel='pps temp',file_name='misc_png/cal_pps_temp.png')
my_plot_histogram_hexbin(bt11[isCloudy], cal_temp[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='cal temp',file_name='misc_png/bt11_cal_temp.png')
my_plot_histogram_hexbin(bt11[isCloudy], pps_temp[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='pps temp',file_name='misc_png/bt11_pps_temp.png')

my_plot_histogram_hexbin(cal_height[isCloudy], pps_height[isCloudy],gridsize=30, 
                         txlabel='calipso height', tylabel='pps height',file_name='misc_png/cal_pps_height.png')
my_plot_histogram_hexbin(bt11[isCloudy], cal_height[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='cal_height',file_name='misc_png/bt11_cal_height.png')
my_plot_histogram_hexbin(bt11[isCloudy], pps_height[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='pps height',file_name='misc_png/bt11_pps_height.png')
