#!/usr/bin/env python
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



"""Read all matched data and make some plotting of thresholds
"""

from glob import glob
#import os.path
import os
import numpy as np
from scipy import histogram
import re
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)
import matplotlib.pyplot as plt
from my_dir import ADIR


isACPGv2012=False
isGAC_v2014_morning_sat = True
isGAC_v2014 = True
if isGAC_v2014_morning_sat:
    num_files_to_read = 30*3
    isGAC=True
    ROOT_DIR = ADIR + "/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
    files = glob(ROOT_DIR + "noaa17/5km/20??/06/*/*h5")
    files = files + glob(ROOT_DIR + "metop*/5km/20??/06/*h5")
    figure_name = "figure_morning_sat_"
elif isGAC_v2014:
    num_files_to_read = 30
    isGAC=True
    figure_name = "figure_"
    ROOT_DIR = ADIR + "/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
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

caObj = CalipsoImagerTrackObject()
for filename in files:
    print os.path.basename(filename)
    caObj = caObj + readCaliopImagerMatchObj(filename)

isCloudyPPS = np.logical_and(caObj.imager.all_arrays['cloudtype']>4,
                             caObj.imager.all_arrays['cloudtype']<21)
isCloudyCaliop = caObj.calipso.all_arrays['number_layers_found']>0
isCloudyCaliop = np.logical_and(isCloudyCaliop,
                                 caObj.calipso.all_arrays['total_optical_depth_5km']>0.5)

cal_temp = caObj.calipso.all_arrays['midlayer_temperature'][:,0]+273.15
cal_height = 1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0]
pps_temp = caObj.imager.all_arrays['ctth_temperature']
pps_height = caObj.imager.all_arrays['ctth_height']
bt11 = caObj.imager.all_arrays['bt11micron']
tsur = caObj.imager.all_arrays['surftemp']
ctype = caObj.imager.all_arrays['cloudtype']

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
