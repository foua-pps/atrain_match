
"""Read all matched data and make some plotting
"""
import os
from glob import glob
import numpy as np
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       ppsMatch_Imager_CalipsoObject)
import matplotlib.pyplot as plt
from get_flag_info import get_calipso_clouds_of_type_i

cc_type_name={
   0: 'low overcast, transparent',
   1: 'low overcast, opaque',
   2: 'transition stratocumulus',
   3: 'low, broken cumulus',
   4: 'altocumulus (transparent)',
   5: 'altostratus (opaque)',
   6: 'cirrus (transparent)',
   7: 'deep convective (opaque)',
}

def plot_height_bias_histograms(caObj):
   height_c = 1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] 
   height_pps = caObj.avhrr.all_arrays['ctth_height']
   height_pps[height_pps>=0] = height_pps[height_pps>=0] + caObj.calipso.all_arrays['elevation'][height_pps>=0] 
   bias = 0.001*(height_pps - height_c)
   use = np.logical_and(caObj.avhrr.all_arrays['ctth_height']>-8,
                        #caObj.calipso.all_arrays['number_layers_found'][:]>0)
                        caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-8)

   bins = 0.5*np.array(xrange(-20,10, 1))
   fig = plt.figure(figsize = (16,10))
   ax = fig.add_subplot(331)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=0))], bins, alpha=0.5, label='0 transparent overcast', facecolor ='r')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(332)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=1))], bins, alpha=0.5, label='1 opaque overcast')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(333)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=2))], bins, alpha=0.5, label='2 transition cumulus')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(334)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=3))], bins, alpha=0.5, label='3 low broken cumuls')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(335)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=4))], bins, alpha=0.5, label='4 altocumuls')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(336)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=5))], bins, alpha=0.5, label='5 altostratus')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(337)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=6))], bins, alpha=0.5, label='6 cirrus', facecolor ='r')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(338)
   plt.hist(bias[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=7))],bins, alpha=0.5, label='7 deep convective')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(339)
   bins = np.array(xrange(0,20000, 1000))
   plt.hist(height_pps[use], bins, alpha=0.5,  label='pps heights')
   plt.legend(loc='upper right')
   plt.suptitle('Ctth height bias')
   #plt.legend(loc='upper left')
   #ax.set_ylabel('number of observations')
   #ax.set_xlabel("bias ")
   #plt.title("ctth height bias histogram (day)")
   #ax.set_xlim(-20,30)
   plt.savefig("ctth_bias_temp_hist.png")
   plt.show()

def calculate_lapse_rate_h(caObj):
    ttro = caObj.avhrr.all_arrays['ttro']
    tsur = caObj.avhrr.all_arrays['surftemp']
    tsur = caObj.avhrr.all_arrays['segment_nwp_temp'][:,0] 
    d45 = (caObj.avhrr.all_arrays['bt11micron'] - 
           caObj.avhrr.all_arrays['bt12micron'])
    temperature_pps = caObj.avhrr.all_arrays['bt12micron'].copy()#
    #temperature_pps = caObj.avhrr.all_arrays['ctth_temperature']    
    #temperature_pps[d45>1.0] = caObj.avhrr.all_arrays['ctth_temperature'][d45>1.0]
    temp_diff = temperature_pps - tsur
    rate = 6.50 #-5K per kilometer
    ##rate = 4.1329 - 0.0225*lat + 0.0005*lat*lat
    rate_neg = -1.0/6.0
    rate_pos = +1.0/3.0
    new_pps_h = rate_neg*temp_diff*1000 
    new_pps_h[temp_diff>0] = rate_pos*temp_diff[temp_diff>0]*1000 
    new_pps_h[new_pps_h<100] = 100
    #For cirrus we need to keep the old height somehow!
    #Question is how to find the high cirrus clouds?
    height_pps = caObj.avhrr.all_arrays['ctth_height']
    #new_pps_h[new_pps_h>4000] = height_pps[new_pps_h>4000]
    #new_pps_h[np.logical_and( d45>1.0,new_pps_h>2000)] = height_pps[np.logical_and( d45>1.0,new_pps_h>2000)]
    #OK for avhrr!!keep = ((tsur - caObj.avhrr.all_arrays['ctth_temperature'])/(tsur-ttro))>0.33
    #new_pps_h[keep] = height_pps[keep]
    #keep = (tsur - temperature_pps)>20
    #new_pps_h[keep] = height_pps[keep]
    #keep = new_pps_h>(1000*caObj.calipso.all_arrays['tropopause_height']*0.25)
    #keep = (caObj.avhrr.all_arrays['warmest_t11'] - caObj.avhrr.all_arrays['coldest_t12'])>10.0
    #test3 = (caObj.avhrr.all_arrays['modis_34'] - 
    #         caObj.avhrr.all_arrays['modis_33'])
    #keep = np.logical_or(test3>-4, test3<-10)
    #new_pps_h[keep] = height_pps[keep]
    #keep =  new_pps_h>4000
    #keep = caObj.avhrr.all_arrays['r13micron']>5.0
    #new_pps_h[keep] = height_pps[keep]
    keep = (caObj.avhrr.all_arrays['bt12micron'] - 
            caObj.avhrr.all_arrays['bt86micron'])<0.0
    new_pps_h[keep] = height_pps[keep]
    return new_pps_h

def plot_sub_scatter_plot_bias_diff_pps(caObj,cc_type, use, bias, bias_new, ax,
                                    xmin=-15000, xmax=15000, title=''):
   import time
   plt.title("%d"%(cc_type))
   my_use = np.logical_and(use,caObj.avhrr.all_arrays['cloudtype']==cc_type)
   plot_sub_scatter_plot_bias_diff_inner(caObj, my_use, bias, bias_new, ax,
                                         xmin=xmin, xmax=xmax)
def plot_sub_scatter_plot_bias_diff(caObj,cc_type, use, bias, bias_new, ax,
                                    xmin=-5000, xmax=15000, title=''):
   plt.title(cc_type_name[cc_type])
   my_use = np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=cc_type))
   plot_sub_scatter_plot_bias_diff_inner(caObj, my_use, bias, bias_new, ax,
                                         xmin=xmin, xmax=xmax)

def plot_sub_scatter_plot_bias_diff_inner(caObj, use, bias, bias_new, ax,
                                          xmin=-5000, xmax=15000):
   my_use = np.logical_and(use, bias<xmax)
   my_use = np.logical_and(my_use, bias_new<xmax)
   my_use = np.logical_and(my_use, bias>xmin)
   my_use = np.logical_and(my_use, bias_new>xmin)
   ymin=xmin
   ymax=xmax
   bins = (xmax-xmin)/100
   from scipy.stats import gaussian_kde
   x=bias_new[my_use].ravel()
   if len(x)==0:
      return
   y=bias[my_use].ravel()
   binsize=100
   edges=np.array(xrange(xmin,xmax+binsize,binsize))
   H, xe, ye = np.histogram2d(x,y,bins=edges)
   xi = np.floor((x - edges[0])/binsize).astype(np.int)               
   yi = np.floor((y - edges[0])/binsize).astype(np.int)  #-1?
   z=H[xi,yi] 
   #z=H[yi,xi] 
   
   perc90 = np.percentile(z,90)
   z[z>perc90]=perc90
   #nicer but too slow
   #points = np.vstack([x,y])   
   #kde= gaussian_kde(points)
   #z=kde(points)
   idx=z.argsort()
   xc=np.array(xrange(xmin,xmax,10))
   y1 =xc
   y2= -xc
   ax.fill_betweenx(xc, y2, y1, facecolor='green', alpha=0.3)
   #plt.scatter(x[idx], y[idx], c=z[idx], edgecolor='', cmap='OrRd', label=label)
   plt.scatter(x[idx], y[idx], c=z[idx], edgecolor='', cmap='viridis', alpha=0.2)
   ax.text(xmax,ymin, 
           "n %d\nbias old: %d\nbias new: %d\nstd old: %d\nstd new %d"%(
              len(bias), np.mean(bias[use]),np.mean(bias_new[use]), 
              np.std(bias[use]), np.std(bias_new[use])),
           verticalalignment='bottom', horizontalalignment='right',
           bbox={'facecolor':'white', 'alpha':1.0, 'pad':5}
        )
   #plt.scatter(x,y, alpha=0.1, label=label)
   #plt.legend(loc='upper left')
   min_val=np.min([np.abs(xmax),abs(xmin)])
   plt.plot(np.mean(x), np.mean(y),'r*')
   plt.plot([-min_val, min_val],[min_val,-min_val],'k')
   plt.plot([xmin, xmax],[ymin,ymax],'k')
   plt.plot([-1000,1000],[1000,1000],'k')
   plt.plot([-1000,1000],[-1000,-1000],'k')
   plt.plot([-1000,-1000],[-1000,1000],'k')
   plt.plot([1000,1000],[-1000,1000],'k')
   ax.set_ylim(xmin,xmax)
   ax.set_xlim(ymin,ymax)

def plot_feature_histograms(caObj,filename):
   new_pps_h = calculate_lapse_rate_h(caObj)
   height_c = 1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] - caObj.calipso.all_arrays['elevation']
   height_pps = caObj.avhrr.all_arrays['ctth_height']
   bias = height_pps -height_c
   bias_new = new_pps_h - height_c

   d45 = caObj.avhrr.all_arrays['bt11micron']-caObj.avhrr.all_arrays['bt12micron']
   feature = d45
   feature = caObj.avhrr.all_arrays['r09micron']/caObj.avhrr.all_arrays['r06micron']
   feature = caObj.avhrr.all_arrays['text_t37t12']
   feature = caObj.avhrr.all_arrays['r06micron']
   dt11ts =caObj.avhrr.all_arrays['bt11micron'] - caObj.avhrr.all_arrays['surftemp']
   #feature =caObj.avhrr.all_arrays['warmest_t37'] - caObj.avhrr.all_arrays['coldest_t37']
   test2 = caObj.avhrr.all_arrays['warmest_t37'] - caObj.avhrr.all_arrays['coldest_t37']
   test2[caObj.avhrr.all_arrays['warmest_t37']<0]=0
   test2[caObj.avhrr.all_arrays['coldest_t37']<0]=0
   test = caObj.avhrr.all_arrays['warmest_t11'] - caObj.avhrr.all_arrays['coldest_t12']
   test3 = caObj.avhrr.all_arrays['modis_34'] - caObj.avhrr.all_arrays['modis_33']
   feature = caObj.avhrr.all_arrays['bt37micron']-caObj.avhrr.all_arrays['bt12micron']

   feature = ((caObj.avhrr.all_arrays['surftemp'] - caObj.avhrr.all_arrays['ctth_temperature'])/(caObj.avhrr.all_arrays['surftemp'] - caObj.avhrr.all_arrays['ttro']))
   #feature = test - d45
   feature = caObj.avhrr.all_arrays['bt12micron']-caObj.avhrr.all_arrays['bt86micron']
   use = np.logical_and(caObj.avhrr.all_arrays['ctth_height']>-8,
                        #caObj.calipso.all_arrays['number_layers_found'][:]>0)
                        caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-8)
   use_temp = np.logical_and(
      use,
      np.logical_or(caObj.calipso.all_arrays['total_optical_depth_5km']>1.0,
                    caObj.calipso.all_arrays['total_optical_depth_5km']<0.0))
   #use= np.logical_and(use,#BRA
   #                    np.logical_and(test3<-4, test3>-10))
   #use= np.logical_and(use,#BRA
   #                    new_pps_h<4000)
   #use= np.logical_and(use,
   #                    caObj.avhrr.all_arrays['ctth_height']<5000)
   #use = np.logical_and(use, np.abs(bias_new)>np.abs(bias))
   #use = np.logical_and(use, test - d45<10)
   #use = np.logical_and(use, caObj.avhrr.all_arrays['bt11micron']>245)
   #use = np.logical_and(use,d45<2.0)
   #use = np.logical_and(use,test>0)
   #bins = 273.15 +4.0*np.array(xrange(-20,20, 1))
   bins = np.percentile(feature,[1,5,10,15,20,25,50,75,80,85,90,95,99])#2*0.5*np.array(xrange(-30,30, 1))
   #bins=0.1*np.array(xrange(-10,20,1))
   fig = plt.figure(figsize = (16,10))
   ax = fig.add_subplot(331)
   #font = {'family' : 'normal',
   #        'weight' : 'bold',
   #        'size'   : 12}
   #plt.rc('font', **font)
#0 = low overcast, transparent
#1 = low overcast, opaque
#2 = transition stratocumulus
#3 = low, broken cumulus
#4 = altocumulus (transparent)
#5 = altostratus (opaque)
#6 = cirrus (transparent)
#7 = deep convective (opaque)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=0))], bins, alpha=0.5, label='0 transparent overcast', facecolor ='r')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(332)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=1))], bins, alpha=0.5, label='1 opaque overcast')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(333)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=2))], bins, alpha=0.5, label='2 transition cumulus')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(334)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=3))], bins, alpha=0.5, label='3 low broken cumuls')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(335)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=4))], bins, alpha=0.5, label='4 altocumuls')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(336)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=5))], bins, alpha=0.5, label='5 altostratus')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(337)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=6))], bins, alpha=0.5, label='6 cirrus', facecolor ='r')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(338)
   plt.hist(feature[np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=7))],bins, alpha=0.5, label='7 deep convective')
   plt.legend(loc='upper left')
   #ax.set_ylim(0,600000)
   ax = fig.add_subplot(339)
   bins = np.array(xrange(0,20000, 1000))
   plt.hist(height_pps[use], bins, alpha=0.5,  label='pps heights')
   plt.legend(loc='upper right')
   plt.suptitle('Ctth height feature')
   plt.savefig("ctth_feature_temp_hist.png")


   fig = plt.figure(figsize = (16,10))
   fig.text(0.5, 0.04, 'new lapse rate height bias', ha='center')
   fig.text(0.04, 0.5, 'old height bias', va='center', rotation='vertical')
   fig.suptitle(filename)
   ax = fig.add_subplot(331) 
   plot_sub_scatter_plot_bias_diff(caObj,0, use, bias, bias_new, ax)
   ax = fig.add_subplot(332)
   plot_sub_scatter_plot_bias_diff(caObj,1, use, bias, bias_new, ax)
   ax = fig.add_subplot(333)
   plot_sub_scatter_plot_bias_diff(caObj,2, use, bias, bias_new, ax)
   ax = fig.add_subplot(334)
   plot_sub_scatter_plot_bias_diff(caObj,3, use, bias, bias_new, ax)
   ax = fig.add_subplot(335)
   plot_sub_scatter_plot_bias_diff(caObj,4, use, bias, bias_new, ax)
   ax = fig.add_subplot(336)
   plot_sub_scatter_plot_bias_diff(caObj,5, use, bias, bias_new, ax)
   ax = fig.add_subplot(337)
   plot_sub_scatter_plot_bias_diff(caObj,6, use, bias, bias_new, ax,
                                   xmin=-15000, xmax=5000)
   ax = fig.add_subplot(338)
   plot_sub_scatter_plot_bias_diff(caObj,7, use, bias, bias_new, ax) 
   ax = fig.add_subplot(339)
   use_temp = np.logical_and(use,caObj.calipso.all_arrays['total_optical_depth_5km']>1.0)
   use_temp = np.logical_and(use_temp,
                             np.logical_or(bias>-6000, bias_new>-6000))
   plot_sub_scatter_plot_bias_diff(caObj,6, use_temp, bias, bias_new, ax, 
                                   xmin=-15000, xmax=5000)
   plt.savefig("ctth_bias_bias_new_scatter_hist.png")


   fig = plt.figure(figsize = (16,10))
   fig.text(0.5, 0.04, 'new lapse rate height bias', ha='center')
   fig.text(0.04, 0.5, 'old height bias', va='center', rotation='vertical')
   fig.suptitle(filename)
   ax = fig.add_subplot(331) 
   plot_sub_scatter_plot_bias_diff_pps(caObj,5, use, bias, bias_new, ax)
   ax = fig.add_subplot(332)
   plot_sub_scatter_plot_bias_diff_pps(caObj,6, use, bias, bias_new, ax)
   ax = fig.add_subplot(333)
   plot_sub_scatter_plot_bias_diff_pps(caObj,7, use, bias, bias_new, ax)
   ax = fig.add_subplot(334)
   plot_sub_scatter_plot_bias_diff_pps(caObj,8, use, bias, bias_new, ax)
   ax = fig.add_subplot(335)
   plot_sub_scatter_plot_bias_diff_pps(caObj,9, use, bias, bias_new, ax)
   ax = fig.add_subplot(336)
   plot_sub_scatter_plot_bias_diff_pps(caObj,10, use, bias, bias_new, ax)
   ax = fig.add_subplot(337)
   plot_sub_scatter_plot_bias_diff_pps(caObj,11, use, bias, bias_new, ax, 
                                   xmin=-15000, xmax=5000)
   ax = fig.add_subplot(338)
   plot_sub_scatter_plot_bias_diff_pps(caObj,12, use, bias, bias_new, ax)
   ax = fig.add_subplot(339)
   plot_sub_scatter_plot_bias_diff_pps(caObj,13, use, bias, bias_new, ax)
   plt.savefig("ctth_bias_bias_new_scatter_hist.png")
   plt.show()
def find_suspicious_files(caObj, filename):
   height_c = (1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] - 
               caObj.calipso.all_arrays['elevation'])
   height_pps = caObj.avhrr.all_arrays['ctth_height']
   bias = height_pps - height_c
   new_pps_h = calculate_lapse_rate_h(caObj)
   bias_l_rate = new_pps_h - height_c

   d45 = caObj.avhrr.all_arrays['bt11micron'] - caObj.avhrr.all_arrays['bt12micron']
   use = np.logical_and(caObj.avhrr.all_arrays['ctth_height']>-8,
                        caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-999)
                        #caObj.calipso.all_arrays['number_layers_found'][:]>0)
   for cc_type in xrange(8):
       use_this_type = np.logical_and(use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=cc_type))
       if len(bias[use_this_type])>1:
          print cc_type, len(bias[use_this_type]), np.mean(bias[use_this_type]), 
          print np.mean(bias_l_rate[use_this_type]), np.mean(d45[use_this_type]), 
          print "%d %d %d %d"%(np.percentile(new_pps_h[use_this_type],10), 
                               np.percentile(new_pps_h[use_this_type],90),
                               np.percentile(height_pps[use_this_type],10), 
                               np.percentile(height_pps[use_this_type],90))
       if len(bias[use_this_type])>2 and np.mean(bias[use_this_type])>20000 :
           print filename, cc_type, np.max(bias[use_this_type]), np.percentile(bias[use_this_type],90)


def prototype_ctth_alg_lapse_rate(caObj, filename):
   from get_flag_info import (get_semi_opaque_info_pps2014,
                              get_calipso_high_clouds,
                              get_calipso_low_clouds)
   low_clouds = get_calipso_low_clouds(caObj)
   high_clouds = get_calipso_high_clouds(caObj)
   semi, opaq = get_semi_opaque_info_pps2014(caObj.avhrr.ctth_status)

   print "prototype new ctth height"
   height_c = (1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] - 
               caObj.calipso.all_arrays['elevation'])
   height_pps = caObj.avhrr.all_arrays['ctth_height']
   bias = height_pps - height_c
   tsur = caObj.avhrr.all_arrays['surftemp']
   tsur = caObj.avhrr.all_arrays['segment_nwp_temp'][:,0] 
   d45 = (caObj.avhrr.all_arrays['bt11micron'] - 
          caObj.avhrr.all_arrays['bt12micron'])
   temperature_pps = caObj.avhrr.all_arrays['bt12micron'].copy()#
   temp_diff = temperature_pps - tsur
   new_pps_h = calculate_lapse_rate_h(caObj)
   bias_l_rate = new_pps_h - height_c

   use = np.logical_and(caObj.avhrr.all_arrays['ctth_height']>-8,
                        caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-999)
   use = np.logical_and(use,caObj.calipso.all_arrays['elevation']<10)
   #use = np.logical_and(use,semi)
   #use = np.logical_and(use,low_clouds)
   #print "bias", np.mean(bias_nona[use]), "bias", np.mean(bias_l_rate_nona[use])
                        #caObj.calipso.all_arrays['number_layers_found'][:]>0)
   for clouds in [low_clouds, high_clouds]:
       use_this_type = np.logical_and(use,clouds)
       print "bias old", np.mean(bias[use_this_type]), "bias", np.mean(bias_l_rate[use_this_type])
       print "std old", np.std(bias[use_this_type]), "std", np.std(bias_l_rate[use_this_type])
       if len(bias[use_this_type])>1:
           print "n", len(bias[use_this_type]), "%d/"%(np.mean(bias[use_this_type])), 
           print "%d"%(np.mean(bias_l_rate[use_this_type])), 
           print "%d/"%(np.mean(bias[np.logical_and(temp_diff>0, use_this_type)])), 
           print "%d"%(np.mean(bias_l_rate[np.logical_and(temp_diff>0, use_this_type)])), 
           print "%0.2f"%(np.mean(d45[use_this_type])), 
           print "%d %d %d %d"%(np.percentile(new_pps_h[use_this_type],10), 
                                np.percentile(new_pps_h[use_this_type],90),
                                np.percentile(height_pps[use_this_type],10), 
                                np.percentile(height_pps[use_this_type],90))

def calculate_lapse_rate(caObj, filename):
   c_height = (1000*caObj.calipso.all_arrays['layer_top_altitude'][:,0] - 
               caObj.calipso.all_arrays['elevation'])
   new_pps_h = calculate_lapse_rate_h(caObj)
   pps_height = caObj.avhrr.all_arrays['ctth_height']
   #bias = 1000*(pps_height -c_height)
   tsur = caObj.avhrr.all_arrays['surftemp']
   tsur = caObj.avhrr.all_arrays['segment_nwp_temp'][:,0]
   temp_diff = 273.15 + caObj.calipso.all_arrays['layer_top_temperature'][:,0] - tsur
   temp_diff_pps = caObj.avhrr.all_arrays['ctth_temperature'][:] - tsur
   rate = 1000*temp_diff/c_height
   t_test = caObj.avhrr.all_arrays['ttro']-caObj.avhrr.all_arrays['surftemp']
   bias_l_rate = new_pps_h - 1000*c_height

   use =  np.logical_and(caObj.calipso.all_arrays['layer_top_temperature'][:,0]>-999,
                         caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-999)
   use =  np.logical_and(use,
                          caObj.avhrr.all_arrays['ctth_temperature'][:]>0)
   #use =  np.logical_and(use,
   #                      get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=6))
   fig = plt.figure(figsize = (16,10))
   ax = fig.add_subplot(321)
   plt.plot(temp_diff[use],c_height[use],'.b',alpha=0.01)
   plt.plot([0,-120],[0,1000*120/6.5],'k')
   ax.set_ylabel('true height')
   ax = fig.add_subplot(323)
   plt.plot(temp_diff[use],pps_height[use],'.r',alpha=0.01)
   plt.plot([0,-120],[0,1000*120/6.5],'k')
   ax.set_ylabel('pps height')
   ax = fig.add_subplot(325)
   plt.plot(temp_diff[use],new_pps_h[use],'.b',alpha=0.01)
   plt.plot([0,-120],[0,1000*120/6.5],'k')
   ax.set_ylabel('new lapse-rate height')
   ax.set_xlabel('true tempdiff')
   ax = fig.add_subplot(322)
   plt.plot(temp_diff_pps[use],c_height[use],'.c',alpha=0.01)
   plt.plot([0,-120],[0,1000*120/6.5],'k')
   ax = fig.add_subplot(324)
   plt.plot(temp_diff_pps[use],pps_height[use],'.g',alpha=0.1)
   ax = fig.add_subplot(326)
   plt.plot(temp_diff_pps[use],new_pps_h[use],'.m',alpha=0.01)
   plt.plot([0,-120],[0,1000*120/6.5],'k')
   ax.set_xlabel('ctth tempdiff')
   #plt.plot(caObj.avhrr.all_arrays['ctth_temperature'][use]-tsur[use],1000*c_height[use],'.b',alpha=0.1)


   fig = plt.figure(figsize = (16,10))
   ax = fig.add_subplot(211)
   plt.plot(t_test[use],c_height[use],'.b',alpha=0.1)
   ax = fig.add_subplot(212)
   plt.plot(t_test[use],pps_height[use],'.r',alpha=0.1)

   plt.show()

   
   my_use = np.logical_and(use, temp_diff>0)
   my_use = np.logical_and(use, temp_diff<0)
   d45 = caObj.avhrr.all_arrays['bt11micron'] - caObj.avhrr.all_arrays['bt12micron']
   for cc_type in xrange(8):
      my_use = np.logical_and(caObj.avhrr.all_arrays['sunz']>95,my_use)
      these = np.logical_and(my_use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=cc_type))
      if len(rate[these])>1:
         perc= np.percentile(rate[these],[10,50,90])
         print cc_type, "%2.1f"%(np.mean(d45[these])), "night", len(rate[these]),
         print np.percentile(pps_height[these],95), "%3.1f"%(np.max(c_height[these])),
         print "%2.1f %2.1f %2.1f"%(perc[0], perc[1], perc[2])
      my_use = np.logical_and(caObj.avhrr.all_arrays['sunz']<60,my_use)
      these = np.logical_and(my_use,get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=cc_type))
      if len(rate[these])>1:
         perc= np.percentile(rate[these],[10,50,90])
         print cc_type, "%2.1f"%(np.mean(d45[these])), "day", len(rate[these]),
         print np.percentile(pps_height[these],95), "%3.1f"%(np.max(c_height[these])),
         print "%2.1f %2.1f %2.1f"%(perc[0], perc[1], perc[2])
          



# ----------------------------------------
if __name__ == "__main__":
    isModis1km = True
    isNPP_v2014 = False
    isGAC_v2014_morning_sat = False
    isGAC_v2014 = True

    if isModis1km:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/global_modis_14th_created20160615/"
        #files = glob(ROOT_DIR + "Reshaped_Files/eos?/1km/????/12/2010??14_*/*h5")
        files = glob(ROOT_DIR + "Reshaped_Files/merged/*night*.h5")
    elif isNPP_v2014:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/sh_reshaped_patch_2014/"
        files = glob(ROOT_DIR + "Reshaped_Files/npp/1km/????/06/arc*/*h5")
    elif isGAC_v2014_morning_sat:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
        files = glob(ROOT_DIR + "noaa17/5km/20??/??/*/*h5")
        files = files + glob(ROOT_DIR + "metop*/5km/20??/??/*/*h5")
        #files = glob(ROOT_DIR + "noaa17/5km/20??/1*/*/*noaa*h5")
    elif isGAC_v2014:
        ROOT_DIR = "/home/a001865/DATA_MISC/reshaped_files/clara_a2_rerun/Reshaped_Files_CLARA_A2_final/"
        files = glob(ROOT_DIR + "merged/noaa18*h5")
        files = files + glob(ROOT_DIR + "merged/noaa19*h5")
        #files = glob(ROOT_DIR + "noaa18/5km/20??/??/*/*noaa*h5")
        #files = files +glob(ROOT_DIR + "noaa19/5km/20??/??/*/*noaa*h5")
        #files = glob(ROOT_DIR + "noaa18/5km/2014/*/*/*noaa*h5")
        #files = glob(ROOT_DIR + "noaa18/5km/2006/10/*/*noaa*20061004_1335*h5")
        #files = glob(ROOT_DIR + "noaa18/5km/2009/12/*/*noaa*20091220_0727*h5")
        #files = files + glob(ROOT_DIR + "noaa18/5km/2010/12/*/*noaa*20101221_0000*h5")
        #files = files + glob(ROOT_DIR + "noaa18/5km/2011/10/*/*noaa*20111011_0146*h5")
        #files = files + glob(ROOT_DIR + "noaa18/5km/2012/04/*/*noaa*20120415_0319*h5")

    caObj = CalipsoAvhrrTrackObject()
    num=0
    for filename in files:
        print  os.path.basename(filename)
        num+=1
        caObj_new=readCaliopAvhrrMatchObj(filename)#, var_to_skip='segment')
        #plot_feature_histograms(caObj_new, os.path.basename(filename))
        prototype_ctth_alg_lapse_rate(caObj_new, filename)
        #calculate_lapse_rate(caObj_new, filename)
        #find_suspicious_files(caObj_new, os.path.basename(filename))
        #plot_height_bias_histograms(caObj_new)
        if num>50:
            find_suspicious_files(caObj, os.path.basename(filename))
            caObj=caObj_new
            num=0
        #caObj = caObj + caObj_new
    #find_suspicious_files(caObj, os.path.basename(filename))
    #plot_height_bias_histograms(caObj)
