#!/usr/bin/env python

"""
Use this script to validate LWP against LWP from Cloudsat
"""

import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma

USE_LAND=False
USE_SEA=True

#For masking out data see code see:    #Selection
#one of the masks is at latitude
NORTH_LIMIT=75.0

#Make plots of LWP:s for each file (independent of selections)
PER_SCENE_PLOT=False

#do not change
CLOUDSAT_FLAG_PRECIP_BIT = range(2,3)  #i.e. only 2


def get_height_from_cloudwatercontent(lwcm,height,cloud_top):

    output=np.zeros(lwcm.shape[1],'f')
    if (cloud_top==1):
        #Find height of cloud top
        for i in range(lwcm.shape[1]):
            for j in range(lwcm.shape[0]):
                if (lwcm[j,i]>0.0):
                    output[i]=height[j,i]
                    break
    if (cloud_top==0):
        #Find height of cloud base
        for i in range(lwcm.shape[1]):
            j=lwcm.shape[0]
            while j>0:
                j=j-1
                if (lwcm[j,i]>0.0):
                    output[i]=height[j,i]
                    break
    return output

    

def get_match_data(filenames):
    """
    Read data from diff files *filenames*, return contented lwp, from the
    respective sources
    Returns lwp_avhrr, lwp_cloudsat, lon, lat
    """
    import h5py
    import matplotlib
    lwp_avhrr = []
    lwp_cloudsat = []
    lwp_cloudsat_uc = []  #uc=uncertainty
    height_cloudsat = []
    lon = []
    lat = []
    elevation = []
    flag = []
    dt = []

    for filename in filenames:
        with h5py.File(filename, 'r') as f:

            #Read avhrr data
            node=f['avhrr']
            try:
                lwp_avhrr.append(node['lwp'][:])
                if (PER_SCENE_PLOT):
                    tmp_lwpa=node['lwp'][:]
            except:
                print "lwp avhrr is missing for %s"%(filename)
                continue
            lon.append(node['longitude'][:])
            lat.append(node['latitude'][:])

            if (PER_SCENE_PLOT):
                tmp_lata=node['latitude'][:]
                tmp_lona=node['longitude'][:]
                tmp_timea=node['sec_1970'][:]

            #Read cloudsat data
            node=f['cloudsat']
            #  Is should be RVOD... not LO_RVOD...
            #  RVOD is the real value; LO_RVOD is if-assuming-only-water
            lwp_cloudsat.append(node['RVOD_liq_water_path'][:])
            lwp_cloudsat_uc.append(node['RVOD_liq_water_path_uncertainty'][:])
            #cloud height is calculated from liquid water content (i.e. where
            #it differs from 0.0), combined with Height
            cloud_top_height=get_height_from_cloudwatercontent(\
                node['RVOD_liq_water_content'][:],node['Height'][:],1)
            cloud_base_height=get_height_from_cloudwatercontent(\
                node['RVOD_liq_water_content'][:],node['Height'][:],0)
            height_cloudsat.append(cloud_top_height)
            #SHq: could possibly take cloud base height instead
            #height_cloudsat.append(cloud_base_height)
            elevation.append(node['elevation'][:])
            flag.append(node['RVOD_CWC_status'][:].astype(np.int))

            #This is already combined from avhrr and cloudsat
            diff_time=f['diff_sec_1970'][:]
            dt.append(diff_time)


            #Make plots for LWP:s of this scene
            if (PER_SCENE_PLOT):
                tmp_lwpc=node['RVOD_liq_water_path'][:]
                tmp_latc=node['latitude'][:]
                tmp_lonc=node['longitude'][:]
                tmp_timec=node['sec_1970'][:]

                #tmp_lwpc, setting nodata (neg. values) to 0
                tmp_lwpc=np.where(np.less(tmp_lwpc,0.0),0.0,tmp_lwpc)

                #Plot those data
                fig=matplotlib.pyplot.figure()
                ax=fig.add_subplot(111)
                ax.plot(tmp_lwpa,'r',label="AVHRR-LWP")
                ax.plot(tmp_lwpc,'b',label="Cloudsat-LWP")
                #figname="LWP_AVHRR_Cloudsat_%s.pdf"%(filename[102:128])
                figname="LWP_AVHRR_Cloudsat_%s.pdf"%(filename[112:138])
                fig.savefig(figname)

                #Also plot lat/lon and time
                fig=matplotlib.pyplot.figure()
                ax=fig.add_subplot(111)
                ax.plot(tmp_lata,'r')
                ax.plot(tmp_latc,'b')
                #figname="latitude_%s.pdf"%(filename[102:128])
                figname="latitude_%s.pdf"%(filename[112:138])
                fig.savefig(figname)
                fig=matplotlib.pyplot.figure()
                ax=fig.add_subplot(111)
                ax.plot(tmp_lona,'r')
                ax.plot(tmp_lonc,'b')
                #figname="longitude_%s.pdf"%(filename[102:128])
                figname="longitude_%s.pdf"%(filename[112:138])
                fig.savefig(figname)
                fig=matplotlib.pyplot.figure()
                ax=fig.add_subplot(111)
                ax.plot(tmp_timea,'r')
                ax.plot(tmp_timec,'b')
                #figname="sec1970_%s.pdf"%(filename[102:128])
                figname="sec1970_%s.pdf"%(filename[112:138])
                fig.savefig(figname)

            
    #Make data from different scenes into one array
    lwp_avhrr_array = np.concatenate(lwp_avhrr)
    lwp_cloudsat_array = np.concatenate(lwp_cloudsat)
    lwp_cloudsat_uc_array = np.concatenate(lwp_cloudsat_uc)
    height_cloudsat_array = np.concatenate(height_cloudsat)
                                         
    lon_array = np.concatenate(lon)
    lat_array = np.concatenate(lat)
    elevation_array = np.concatenate(elevation)
    flag_array = np.concatenate(flag)
    dt_array = np.concatenate(dt)


    #Selection
    #Make a selection of which data to keep/remove!
    #Any criterias can be combined, but remember first selection has got
    # a different call!
    

    #Mask for lwp_avhrr no-data (-9.0)
    selection=np.where(np.less(lwp_avhrr_array,-0.1),False,True)
    selection_str="LWPavhrr>=-0.1 "

    #Maks for lwp_cloudsat no-data (like -3333 or -4444) and no-clouds (0.0)
    selection=np.where(np.less(lwp_cloudsat_array,0.1),False,selection)
    selection_str=selection_str+"; LWPcloudsat>=0.1 "
    
    #mask out cloudsat lwp>180   (SHq: do no use!)
    #selection=np.where(np.greater(lwp_cloudsat_array,180.0),False,selection)
    #selection_str=selection_str+"; LWPcloudsat<=180 "

    #Mask out px with uncertain cloudsat
    selection=np.where(np.greater(lwp_cloudsat_uc_array,100.0),False,selection)
    selection_str=selection_str+"; LWPcloudsat_uncertain<=100 "

    #Mask out px with low cloudtop (as they are bad in cloudsat)
    selection=np.where(np.less(height_cloudsat_array,1000.0),False,selection)
    selection_str=selection_str+"; cloud top height(cloudsat)>=1000 "

    #Mask out px with flag=precip
    from validate_cph import get_bits
    selection=np.where(np.equal(get_bits(flag_array,CLOUDSAT_FLAG_PRECIP_BIT,
                                         shift=True),1),
                       False,selection)
    selection_str=selection_str+"; cloudsat flag!=precip "

    #Mask out lat>60
    selection=np.where(np.greater(lat_array,NORTH_LIMIT),False,selection)
    #selection=np.where(np.less_equal(lat_array,NORTH_LIMIT),False,selection)
    selection_str=selection_str+"; south of %d "%(NORTH_LIMIT)
     

        
    #select land and/or sea (setting done in the top of this file)
    if (USE_LAND and not(USE_SEA)):
        selection=np.where(np.less(elevation_array,0.0),False,selection)
        selection_str=selection_str+"; land "

    if (not(USE_LAND) and USE_SEA):
        selection=np.where(np.greater_equal(elevation_array,0.0),False,selection)
        selection_str=selection_str+"; sea "
                
    if (USE_LAND and USE_SEA):
        selection_str=selection_str+"; land and sea"


    print " "
    print "Selection made is:"
    print selection_str

    sz=np.sum(selection)


    #Masking is problematic, do it elementwise...
    lwp_avhrr_array_new=np.zeros(sz)
    lwp_cloudsat_array_new=np.zeros(sz)
    lon_array_new=np.zeros(sz)
    lat_array_new=np.zeros(sz)
    diff_new=np.zeros(sz)
    dt_new=np.zeros(sz)
    j=-1
    for i in range(len(selection)):
        if (selection[i]):
            j=j+1
            lwp_avhrr_array_new[j]=lwp_avhrr_array[i]
            lwp_cloudsat_array_new[j]=lwp_cloudsat_array[i]
            lon_array_new[j]=lon_array[i]
            lat_array_new[j]=lat_array[i]
            diff_new[j]=lwp_avhrr_array[i]-lwp_cloudsat_array[i]
            dt_new[j]=dt_array[i]
    

    return lwp_avhrr_array_new, lwp_cloudsat_array_new, \
           lon_array_new, lat_array_new,diff_new


def validate_all(filenames):
    """
    Use lwp_diff in all files in *filenames* for validation.
    
    """

    #Read data from match file, and make a selection of pixels
    lwp_avhrr, lwp_cloudsat, lon, lat, lwp_diff = \
        get_match_data(filenames)

    #Calculate statistics
    mean = np.mean(lwp_diff)
    median = np.median(lwp_diff)
    std = np.std(lwp_diff)

    print " "
    print "Statistics for LWP difference"
    print "================================"
    print("Number of pixels: %d" % len(lwp_diff))
    print("Mean: %.2f" % mean)
    print("Median: %.2f" % median)
    print("Standard deviation: %.2f" % std)

    #Calculate statistics per data set
    print " "
    print "Statistics for each LWP data set"
    print "================================"
    print "lwp_avhrr"
    print ("  mean: %.2f" % np.mean(lwp_avhrr))
    print ("  median: %.2f" % np.median(lwp_avhrr))
    print ("  std: %.2f" % np.std(lwp_avhrr))
    print ("  min: %.2f" % np.min(lwp_avhrr))
    print ("  max: %.2f" % np.max(lwp_avhrr))
    print "lwp_cloudsat"
    print ("  mean: %.2f" % np.mean(lwp_cloudsat))
    print ("  median: %.2f" % np.median(lwp_cloudsat))
    print ("  std: %.2f" % np.std(lwp_cloudsat))
    print ("  min: %.2f" % np.min(lwp_cloudsat))
    print ("  max: %.2f" % np.max(lwp_cloudsat))


    pixels_no=ma.count(lwp_diff)


    #Make diagrams and plots
    from amsr_avhrr.plotting import plot_hist, density, distribution_map
    hist_range = (np.percentile(lwp_diff, 1),
                  np.percentile(lwp_diff, 99))
    fig = plot_hist(lwp_diff, bins=500, range=hist_range)
    fig.axes[0].set_xlabel('lwp difference (g m**-2)')
    fig.suptitle("CPP lwp - cloudsat lwp\nPixels left: %d" %
                 (pixels_no))
    fig.savefig('validate_all.pdf')
    
    # Density plot
    fig2 = density(lwp_avhrr, lwp_cloudsat,
                   bins=xrange(0, 255))
    fig2.axes[0].set_xlabel('CPP lwp (g m**-2)')
    fig2.axes[0].set_ylabel('Cloudsat lwp (g m**-2)')
    fig2.suptitle("Number of pixels: %d" %
                 (pixels_no))
    fig2.savefig('density_all.pdf')
    fig2 = density(lwp_avhrr, lwp_cloudsat,
                   bins=xrange(0, 505))
    fig2.axes[0].set_xlabel('CPP lwp (g m**-2)')
    fig2.axes[0].set_ylabel('Cloudsat lwp (g m**-2)')
    fig2.suptitle("Number of pixels: %d" %
                 (pixels_no))
    fig2.savefig('density_all_big.pdf')

    # Map of pixel distribution
    #2012-06-26: This does not work any more! Why? /Sara Hornquist
    #fig3 = distribution_map(lon, lat)
    #fig3.suptitle("Distribution of valid pixels\n" +
    #              "Number of Pixels: %d" % pixels_no)
    #fig3.savefig('distribution_all.pdf')
    
    return mean, median, std, lwp_diff



if __name__ == '__main__':
    import sys
    filenames = sys.argv[1:]
    validate_all(filenames)
    #plt.show()
