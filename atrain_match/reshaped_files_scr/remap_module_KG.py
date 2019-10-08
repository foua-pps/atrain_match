#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read all COT profiles of Probability of detecting a cloud layer with a fixed COT: POD(cloudy_layer,cot). Then find the COT value for each grid point where POD first exceeds 50 % and plot the resulting map.
"""

import os
import pdb
import h5py
from glob import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import copy
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util_KG import (read_POD_CLOUDY_LAYER_afternoon, read_POD_CLOUDY_LAYER_morning)

layer_cots = [0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.55,0.65,0.75,0.85,0.95,1.5,2.5,3.5,4.5]
layer_cot_names = ["0_025","0_075","0_125","0_175","0_225","0_275","0_325","0_375","0_425","0_475","0_55","0_65","0_75","0_85","0_95","1_50","2_50","3_50","4_50"]

N_layers = 19

#BASE_PLOT_DIR = "/nobackup/smhid13/sm_kgkar/atrain_match_CCI_V3/plot_stats/figs_stats_global_cotranges_TCC_POD"
BASE_PLOT_DIR = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD"
BASE_PLOT_DIR_MORNING = "/nobackup/smhid12/sm_kgkar/atrain_match_CALIPSOv4_CPP_Simulator_results_CLARA_A2/plot_stats/figs_stats_global_cotranges_TCC_POD_morning"
DNT = "all"



def remap_and_plot_score_on_several_areas(satellites,detection_limit, lats, lons, radius_km, vmin=0.0, vmax=1.0, 
                                          score='Minimum_COT', screen_out_valid=False):
    print score
    PLOT_DIR_SCORE = BASE_PLOT_DIR + "/%s/%s/"%(score, satellites)
    PLOT_FILENAME_START = "fig_%s_filter_dnt_%s_%s_r%skm_"%(satellites,DNT,score,radius_km)
    
    if not os.path.exists(PLOT_DIR_SCORE):
        os.makedirs(PLOT_DIR_SCORE)
##     for plot_area_name in [
##         #'cea5km_test'
##         #'euro_arctic',
##         #'ease_world_test'
##         'euro_arctic',
##         'antarctica',
##         'npole',
##         'ease_nh_test',
##         'ease_sh_test' ]:
##         self._remap_a_score_on_an_area(plot_area_name=plot_area_name, 
##                                            vmin=vmin, vmax=vmax, score=score)
    #the real robinson projection
    if "metop" not in satellites:
        remap_a_score_on_an_robinson_projection(satellites, detection_limit, lats, lons, radius_km, vmin=vmin, vmax=vmax, 
                                                          score=score, screen_out_valid=False)




def remap_a_score_on_an_robinson_projection(satellites, detection_limit, lats, lons, radius_km, vmin=0.0, vmax=5.0, 
                                                 score='Minimum_COT', screen_out_valid=False):
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from scipy.interpolate import griddata
    PLOT_DIR_SCORE = BASE_PLOT_DIR + "/%s/%s/"%(score, satellites)
    PLOT_FILENAME_START = "fig_%s_filter_dnt_%s_%s_r%skm_"%(satellites,DNT,score,radius_km)
    plt.close('all')
    data = detection_limit
##     ma_data = getattr(self, score)
##         the_mask = ma_data.mask
##         data=ma_data.data
##         #data[np.logical_and(data>vmax,~the_mask)] = vmax
##         #data[np.logical_and(data<vmin,~the_mask)] = vmin
##         #reshape data a bit
##         ind = np.argsort(lats)
##         lons = lons[ind]
##         lats = lats[ind]
##         data = data[ind]
##         the_mask = the_mask[ind]
##         ind = np.argsort(lons)
##         lons = lons[ind]
##         lats = lats[ind]
##         data =data[ind]
##         the_mask = the_mask[ind]
##         lons =         lons.reshape(len(data),1)#*3.14/180
##         lats =         lats.reshape(len(data),1)#*3.14/180
##         data =         data.reshape(len(data),1)
##         the_mask =     the_mask.reshape(len(data),1)

    my_proj1 = Basemap(projection='robin',lon_0=0,resolution='c')
    numcols=1000
    numrows=500
    lat_min = -83.0
    lon_min = -179.9
    lat_max = 83.0
    lon_max = 179.9
            
    fig = plt.figure(figsize = (16,9))
    ax = fig.add_subplot(111)
    import copy; 
    #my_cmap=copy.copy(matplotlib.cm.Reds)
    #my_cmap=copy.copy(matplotlib.cm.bwr)
    my_cmap=make_my_own_color_map()
##         if score in "Bias" and screen_out_valid:
##             #This screens out values between -5 and +5% 
##             vmax=25
##             vmin=-25            
##             my_cmap=copy.copy(matplotlib.cm.get_cmap("BrBG", lut=100))
##             cmap_vals = my_cmap(np.arange(100)) #extractvalues as an array
##             cmap_vals[39:61] = [0.9, 0.9, 0.9, 1] #change the first value
##             my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
##                 "newBrBG", cmap_vals) 
##             print my_cmap
##         if score in "RMS" and screen_out_valid:
##             # This screens out values beteen 0 and 20%. 41/100=20%
##             vmax=50
##             vmin=0
##             my_cmap=copy.copy(matplotlib.cm.get_cmap("BrBG", lut=100))
##             cmap_vals = my_cmap(np.arange(100)) #extract values as an array
##             cmap_vals[0:41] = [0.9, 0.9, 0.9, 1] #change the first value
##             my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
##                 "newBrBG", cmap_vals) 
##             print my_cmap

      

        #to mask out where we lack data
##         data[np.logical_and(data>vmax,~the_mask)]=vmax
##         data[np.logical_and(data<vmin,~the_mask)]=vmin
##         data[the_mask]=2*vmax #give no data value that will be masked white
    xi = np.linspace(lon_min, lon_max, numcols)
    yi = np.linspace(lat_min, lat_max, numrows)
    xi, yi = np.meshgrid(xi, yi)
    # interpolate
    x, y, z = (np.array(lons.ravel()), 
               np.array(lats.ravel()), 
               np.array(data.ravel()))
    my_cmap.set_over(color='0.9',alpha=1)
    zi = griddata((x, y), z, (xi, yi), method='nearest')
    im1 = my_proj1.pcolormesh(xi, yi, zi, cmap=my_cmap,
                              vmin=vmin, vmax=vmax, latlon=True)
    #im1.set_clim([vmin,vmax]) #to get nice ticks in the colorbar
    #im1.set_clim([-4.775,5.0]) #to center around value 0.225!
    #draw som lon/lat lines
    my_proj1.drawparallels(np.arange(-90.,90.,30.))
    my_proj1.drawmeridians(np.arange(-180.,180.,60.))
    my_proj1.drawcoastlines()
    my_proj1.drawmapboundary(fill_color='1.0') #0.9 light grey"
    cb = my_proj1.colorbar(im1,"right", size="5%", pad="2%")
    tick_locator = matplotlib.ticker.MaxNLocator(nbins=10)
    cb.locator = tick_locator
    cb.ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    cb.update_ticks()

    ax.set_title(score)
    plt.savefig(PLOT_DIR_SCORE + PLOT_FILENAME_START+
                '_robinson_newLUT' +'.png')
##     plt.savefig(PLOT_DIR_SCORE + PLOT_FILENAME_START+
##                 '_robinson_newLUT' +'.eps')
    plt.close('all')


"""
A custom color map to be used for plot of scores at projected areas.
The colormap uses the bwr colormap for low values.
And for higher values the hot colormap in inverted directio is used
this adds darker red and brown colors for the really high values.
Just comment out the plotting if you do not want it.
"""

def make_my_own_color_map(minv=0.0, maxv=5.0, colorchange1=0.22, colorchange2=1.0):

    
    n_resolution = 1000
    step = (maxv-minv)*1.0/n_resolution
    n_res_darkred = int((maxv-colorchange2)/(step*0.4)) #will use ~40% of colormap
    #copy existing colormap, important to copy!
    my_cmap = copy.copy(matplotlib.cm.get_cmap("bwr",
                                             lut=n_resolution)) #blue to red
    cmap_vals = my_cmap(np.arange(n_resolution)) #extract values as an array
    new_cmap_vals = cmap_vals.copy()
    my_cmap_darkred=copy.copy(
        matplotlib.cm.get_cmap("hot",
                               lut=n_res_darkred)) #get red and brown colors
    cmap_vals_darkred = my_cmap_darkred(np.arange(n_res_darkred)) #values
    cmap_vals_darkred = cmap_vals_darkred[::-1]

    #first_color
    steps_first_color = int((colorchange1-minv)/step)
    every_xth_needed = int(n_resolution*0.5/steps_first_color)
    for i in  xrange(steps_first_color):
        new_cmap_vals[i] = cmap_vals[i*every_xth_needed]
    #second_color
    steps_second_color = int((colorchange2-colorchange1)/step)
    every_xth_needed = int(n_resolution*0.5/steps_second_color)
    for i in  xrange(steps_second_color):
        index_br = int(n_resolution*0.5) + i*every_xth_needed
        new_cmap_vals[steps_first_color+i] = cmap_vals[index_br ]
    #third_color
    steps_second_color = int((colorchange2-colorchange1)/step)
    every_xth_needed = int(n_resolution*0.5/steps_second_color)
    for i in  xrange(n_resolution-steps_second_color-steps_first_color):
        index = steps_first_color + steps_second_color+i
        print index, int(0.6*n_res_darkred), i, n_res_darkred, int(0.6*n_res_darkred)+i
        new_cmap_vals[index] = cmap_vals_darkred[int(0.6*n_res_darkred)+i]

    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "newBlueRedMoreRed", new_cmap_vals)
    a=np.array([np.linspace(0,5,num=1000)])
    a=a.reshape(20,50)
    #ax = plt.subplot(111)
    #ax.imshow(a, cmap=my_cmap, interpolation='nearest')
    #plt.show()
    return my_cmap
