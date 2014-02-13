#!/usr/bin/env python
# -*- coding: utf-8 -*-



"""Read all matched data and make some plotting of thresholds
"""

from glob import glob
import os.path
import os
import numpy as np
from scipy import histogram
import re
from pps_threshold_functions import (print_stats, 
                                 get_clear_and_cloudy_vectors)
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
import matplotlib.pyplot as plt

isACPGv2012=False
isGAC=False
ROOT_DIR = "/local_disk/nina_pps/data_validation_ctth_patch_nov2012/VALIDATION_20140116/Reshaped_Files/*/1km/"
files = glob(ROOT_DIR + "/????/??/arc*/*h5")

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
    
cloudobj = get_clear_and_cloudy_vectors(caObj, isACPGv2012, isGAC)




cal_temp = caObj.calipso.all_arrays['cloud_mid_temperature'][0]+273.15 
cal_height = 1000*caObj.calipso.all_arrays['cloud_top_profile'][0] 
pps_temp = caObj.avhrr.all_arrays['ctth_temperature']  
pps_height = caObj.avhrr.all_arrays['ctth_height']  
bt11 = caObj.avhrr.all_arrays['bt11micron']
tsur = caObj.avhrr.all_arrays['surftemp']
ctype = caObj.avhrr.all_arrays['cloudtype']
isCloudyPPS = cloudobj.isPpsCloudy
isCloudyCaliop = cloudobj.isCloudy.all_arrays['All']
isCloudy =  np.logical_and(isCloudyPPS,isCloudyCaliop)
isCloudy =  np.logical_and(pps_temp>-9,isCloudy)
#WantedCloudtypes = n
#isBad = pps_temp[np.logical_and(pps_temp<200,isCloudy)]
#print len(pps_temp[isBad])
#print pps_height[isBad]
#print bt11[isBad]
isCloudy = np.logical_and(pps_temp>150,isCloudy)
#bins=np.array(range(200,300,10))
import numpy as np

for ind in range(5,15):
    print ind
    isCloudytypeT = np.logical_and(isCloudy,ctype==ind)
    my_plot_histogram_hexbin(bt11[isCloudytypeT], np.abs(pps_temp[isCloudytypeT]-cal_temp[isCloudytypeT]),gridsize=30,
                               txlabel='bt11', tylabel='diff cal pps temp',file_name='bt11_cal_minus_pps_temp_klass_%d.png'%(ind))
my_plot_histogram_hexbin(cal_temp[isCloudy], pps_temp[isCloudy],gridsize=30, 
                         txlabel='calipso temp', tylabel='pps temp',file_name='cal_pps_temp.png')
my_plot_histogram_hexbin(bt11[isCloudy], cal_temp[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='cal temp',file_name='bt11_cal_temp.png')
my_plot_histogram_hexbin(bt11[isCloudy], pps_temp[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='pps temp',file_name='bt11_pps_temp.png')

my_plot_histogram_hexbin(cal_height[isCloudy], pps_height[isCloudy],gridsize=30, 
                         txlabel='calipso height', tylabel='pps height',file_name='cal_pps_height.png')
my_plot_histogram_hexbin(bt11[isCloudy], cal_height[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='cal_height',file_name='bt11_cal_height.png')
my_plot_histogram_hexbin(bt11[isCloudy], pps_height[isCloudy],gridsize=30, 
                         txlabel='bt11', tylabel='pps height',file_name='bt11_pps_height.png')
