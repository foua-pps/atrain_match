
"""Read all matched data and make some plotting
"""
import os
from glob import glob
import numpy as np
from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)

import matplotlib.pyplot as plt
from get_flag_info import get_calipso_clouds_of_type_i
from get_flag_info import (get_semi_opaque_info_pps2014,
                           get_day_night_twilight_info_pps2014,
                           get_land_coast_sea_info_pps2014,
                           get_mountin_info_pps2014,
                           get_inversion_info_pps2014,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_low_clouds)

cc_type_name={
   0: 'low overcast (tp) ',
    1: 'low overcast (oq)',
   2: 'transition stratocumulus',
   3: 'low broken cumulus',
   4: 'altocumulus (tp)',
   5: 'altostratus (oq)',
   6: 'cirrus (tp)',
   7: 'deep convective (op)'
}
ct_min_v =[1,2,3,4,5,10,7,8,11]
ct_max_v= [1,2,3,4,6,10,7,9,20]

def plot_ct_table(caObj):

    from get_flag_info import get_calipso_clouds_of_type_i

    print "N,          0     clear clear clear clear very-low, low, medium, high very-high, fractional, cirrus,  cirrus, cirrus,  cirrus-above-low"

    for type_i in xrange(0,8):
       if type_i ==1:
           continue
       is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
       pps_ctype = caObj.avhrr.all_arrays['cloudtype']
       use = np.logical_and(pps_ctype>0,pps_ctype<23)
       is_type_i =np.logical_and(is_type_i, use)
       N = np.sum(is_type_i)
       print "N", ("%d"%(N)).rjust(9,' '),
       for ind_ct in xrange(15):
           these = np.logical_and(is_type_i,pps_ctype== ind_ct)
           print ("%3.1f"%(np.sum(these)*1.0/N*100)).rjust(5,' '),
           if type_i==2 and ind_ct==2:
               plt.plot(caObj.avhrr.all_arrays['longitude'][these], caObj.avhrr.all_arrays['latitude'][these],'b.',alpha=0.4, label="cal-%d,pps-%d"%(type_i,ind_ct))
               plt.legend()
               plt.show()
       print cc_type_name[type_i]  
       
def plot_ct_table2(caObj):

    from get_flag_info import get_calipso_clouds_of_type_i
    fig = plt.figure(figsize = (30,64))
    print "N,           low, frac, medium, high, cirrus"
    for type_i in xrange(0,8):
        if type_i ==1:
            continue
        is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
        pps_ctype = caObj.avhrr.all_arrays['cloudtype']
        use = np.logical_and(pps_ctype>4,pps_ctype<23)
        is_type_i =np.logical_and(is_type_i, use)
        N = np.sum(is_type_i)
        print "N", ("%d"%(N)).rjust(9,' '),
        for ind_ct in xrange(4,9):
            pps_ok = np.logical_and(pps_ctype >= ct_min_v[ind_ct], pps_ctype <= ct_max_v[ind_ct])
            these = np.logical_and(is_type_i, pps_ok)
            print ("%3.1f"%(np.sum(these)*1.0/N*100)).rjust(5,' '),
        print cc_type_name[type_i], np.percentile(caObj.calipso.all_arrays['layer_top_pressure'][:,0][
            np.logical_and(these,
                           caObj.calipso.all_arrays['layer_top_pressure'][:,0]>0)],0.1) 

def plot_ct_table3(caObj):

    from get_flag_info import get_calipso_clouds_of_type_i
    fig = plt.figure(figsize = (30,64))

    print "N,          0     clear clear clear clear very-low, low, medium, high very-high, fractional, cirrus,  cirrus, cirrus,  cirrus-above-low"
    for type_i in xrange(0,8):
       if type_i ==1:
           continue
       is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
       pps_ctype = caObj.avhrr.all_arrays['cloudtype']
       use = np.logical_and(pps_ctype>0,pps_ctype<23)
       is_type_i =np.logical_and(is_type_i, use)
       N = np.sum(is_type_i)
       print "N", ("%d"%(N)).rjust(9,' '),
       for ind_ct in xrange(0,9):
           pps_ok = np.logical_and(pps_ctype >= ct_min_v[ind_ct], pps_ctype <= ct_max_v[ind_ct])
           these = np.logical_and(is_type_i, pps_ok)
           #if type_i==2 and ind_ct==2:
           sub_ind = type_i 
           if sub_ind > 0:
               sub_ind -=1
           ax = fig.add_subplot(7,9,sub_ind*9+ind_ct+1)
           frame1 = plt.gca()
           frame1.axes.get_xaxis().set_ticks([])
           frame1.axes.get_yaxis().set_ticks([])
           plt.plot(caObj.avhrr.all_arrays['longitude'][these], caObj.avhrr.all_arrays['latitude'][these],'b.',alpha=0.4, label="cal-%d,pps-%d-%d"%(type_i, ct_min_v[ind_ct], ct_max_v[ind_ct]))
           #plt.legend()
    plt.show()  

def plot_ct_table4(caObj, use_in=None):

    from get_flag_info import get_calipso_clouds_of_type_i
    fig = plt.figure(figsize = (30,64))
    print "N,           land, sea, snow, ice, low, frac, medium, high, cirrus"
    for type_i in xrange(0,8):
        if type_i ==1:
            continue
        is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
        pps_ctype = caObj.avhrr.all_arrays['cloudtype']
        use = np.logical_and(pps_ctype>0,pps_ctype<23)
        if use_in is not None:
            use = np.logical_and(use, use_in)
        is_type_i =np.logical_and(is_type_i, use)
        N = np.sum(is_type_i)
        print "N", ("%d"%(N)).rjust(9,' '),
        for ind_ct in xrange(0,9):
            pps_ok = np.logical_and(pps_ctype >= ct_min_v[ind_ct], pps_ctype <= ct_max_v[ind_ct])
            these = np.logical_and(is_type_i, pps_ok)
            print ("%3.1f"%(np.sum(these)*1.0/N*100)).rjust(5,' '),
        print cc_type_name[type_i]

def plot_ct_table5(caObj, use_in=None):

    from get_flag_info import get_calipso_clouds_of_type_i
    fig = plt.figure(figsize = (30,64))
    pps_ctype = caObj.avhrr.all_arrays['cloudtype']
    use = np.logical_and(pps_ctype>0,pps_ctype<23)
    use = np.logical_and(use, use_in)
    for type_i in [0,2,3,4,5,6,7]:
        is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
        is_type_i =np.logical_and(is_type_i, use)
        N = np.sum(is_type_i)
        print "N", ("%d"%(N)).rjust(9,' '),
        if N<50:
            print "------------------", cc_type_name[type_i]
            continue 
        these = np.logical_and(is_type_i, pps_ctype <5)
        print ("%3.1f"%(np.sum(these)*1.0/N*100)).rjust(5,' '),
        these = np.logical_and(is_type_i, pps_ctype >4)
        print ("%3.1f"%(np.sum(these)*1.0/N*100)).rjust(5,' '),
        print cc_type_name[type_i]

def plot_ct_table6(caObj, use_in=None):

    from get_flag_info import get_calipso_clouds_of_type_i
    fig = plt.figure(figsize = (30,64))
    pps_ctype = caObj.avhrr.all_arrays['cloudtype']
    use = np.logical_and(pps_ctype>0,pps_ctype<23)
    use = np.logical_and(use, use_in)
    use = np.logical_and(use,caObj.calipso.all_arrays['cloud_fraction']>0.9)
    #for type_i in [0,2,3,4,5,6,7]:
    #    print cc_type_name[type_i].replace(' ','_'),
    print "N      ", 
    for type_i in [0,2,3,4,5,6,7]:
        is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
        is_type_i =np.logical_and(is_type_i, use)
        N = np.sum(is_type_i)
        print ("%d"%(N)).rjust(5,' '),
    print ""
    print "clear: ",
    for type_i in [0,2,3,4,5,6,7]:
        is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
        is_type_i =np.logical_and(is_type_i, use)
        N = np.sum(is_type_i)
        if N<50:
            print ("---").rjust(5,' '),
        else:    
            these = np.logical_and(is_type_i, pps_ctype <5)
            print ("%3.1f"%(np.sum(these)*1.0/N*100)).rjust(5,' '),
    print ""
    print "cloudy:",
    for type_i in [0,2,3,4,5,6,7]:
        is_type_i = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=type_i)
        is_type_i =np.logical_and(is_type_i, use)
        N = np.sum(is_type_i)
        if N<50:
            print ("---").rjust(5,' '),
        else:    
            these = np.logical_and(is_type_i, pps_ctype >4)
            print ("%3.1f"%(np.sum(these)*1.0/N*100)).rjust(5,' '),
    print ""

def table_5_per_illumination(caObj):
    cloudtype_conditions = caObj.avhrr.all_arrays[ 'cloudtype_conditions']
    cloudtype_status = caObj.avhrr.all_arrays[ 'cloudtype_status']
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag) =     get_day_night_twilight_info_pps2014(cloudtype_conditions)
    (no_qflag, land_flag, sea_flag, coast_flag, all_lsc_flag) = get_land_coast_sea_info_pps2014(cloudtype_conditions)
    flag_all = caObj.avhrr.all_arrays[ 'sunz']
    rough_terrain = get_mountin_info_pps2014(cloudtype_conditions)
    inversion = get_inversion_info_pps2014(cloudtype_status)
    light_flag = caObj.avhrr.all_arrays[ 'sunz']<90
    dark_flag = caObj.avhrr.all_arrays[ 'sunz']>=90  
    snow = caObj.calipso.all_arrays['nsidc_surface_type']>50
    snow_free = caObj.calipso.all_arrays['nsidc_surface_type']<50
    
    #for use, name in zip( [ night_flag, twilight_flag, day_flag],
    #                      ["night", "twilight","day"]):
    for use3, name3 in zip( [ snow, snow_free,  all_dnt_flag],
                            ["snow", "snow_free","all"]):
        for use, name in zip( [ light_flag, dark_flag],
                              ["light", "dark"]):
            for use2, name2 in zip( [  land_flag, sea_flag, coast_flag],
                                    ["land", "sea","coast"]):
                print name, name2, name3
                plot_ct_table6(caObj, np.logical_and(np.logical_and(use,use2),use3))

def table_21_do_for_atbd(caObj):
    cloudtype_conditions = caObj.avhrr.all_arrays[ 'cloudtype_conditions']
    cloudtype_status = caObj.avhrr.all_arrays[ 'cloudtype_status']
    (no_qflag, night_flag, twilight_flag, day_flag, all_dnt_flag) =     get_day_night_twilight_info_pps2014(cloudtype_conditions)
    low_clouds = get_calipso_low_clouds(caObj)
    high_clouds = get_calipso_high_clouds(caObj)
    medium_clouds = get_calipso_medium_clouds(caObj)
    caliop_ok = np.logical_or(low_clouds,np.logical_or(high_clouds,medium_clouds))
    pps_ctype = caObj.avhrr.all_arrays['cloudtype']
    use = np.logical_and(pps_ctype>4,pps_ctype<23)
    use = np.logical_and(use, np.logical_and(caliop_ok,caObj.calipso.all_arrays['cloud_fraction']>0.9))
    pps_low = [pps_ctype_i in [5,6,10] for pps_ctype_i in pps_ctype]
    pps_medium = [pps_ctype_i in [7] for pps_ctype_i in pps_ctype] 
    pps_high = [pps_ctype_i in [8,9,11,12,13,14,15,16,17,18] for pps_ctype_i in pps_ctype] 
    print "POD low medium high FAR low medium high"                   
    for use_i, name in zip( [ all_dnt_flag, day_flag, night_flag, twilight_flag ],
                          ["all", "day", "night", "twilight"]):
        


        use_this = np.logical_and(use, use_i)
        
        pps_low_i = np.logical_and(pps_low,use_this)
        low_clouds_i =  np.logical_and(low_clouds,use_this)
        n_low_ok = 1.0*np.sum(np.logical_and(pps_low_i,low_clouds_i))
        POD_low =  100*n_low_ok/np.sum(low_clouds_i)
        FAR_low =  100*(np.sum(pps_low_i)-n_low_ok)*1.0/np.sum(pps_low_i) 
        print n_low_ok, np.sum(pps_low_i), np.sum(low_clouds_i)  
        pps_medium_i = np.logical_and(pps_medium,use_this)
        medium_clouds_i =  np.logical_and(medium_clouds,use_this)
        n_medium_ok = 1.0*np.sum(np.logical_and(pps_medium_i,medium_clouds_i))
        POD_medium =  100*n_medium_ok/np.sum(medium_clouds_i)
        FAR_medium =  100*(np.sum(pps_medium_i)-n_medium_ok)/np.sum(pps_medium_i)
        pps_high_i = np.logical_and(pps_high,use_this)
        high_clouds_i =  np.logical_and(high_clouds,use_this)
        n_high_ok = 1.0*np.sum(np.logical_and(pps_high_i,high_clouds_i))
        POD_high =  100*n_high_ok/np.sum(high_clouds_i)
        FAR_high = 100*(np.sum(pps_high_i)-n_high_ok)/np.sum(pps_high_i) 
        print "%3.1f %3.1f %3.1f %3.1f %3.1f %3.1f"%(POD_low, POD_medium,POD_high, FAR_low, FAR_medium,FAR_high)




# ----------------------------------------
if __name__ == "__main__":
    isModis1km = True
    isNPP_v2014 = False
    isGAC_v2014_morning_sat = False
    isGAC_v2014 = True


    ROOT_DIR_GAC_nn = ("/home/a001865/DATA_MISC/reshaped_files/"
                       "ATRAIN_RESULTS_GAC_nnavhrr_20161202/Reshaped_Files/noaa18/")
    ROOT_DIR_GAC_old = ("/home/a001865/DATA_MISC/reshaped_files/"
                        "ATRAIN_RESULTS_GAC_v2014/Reshaped_Files/noaa18/")

    files = glob(ROOT_DIR_GAC_nn + "5km/2009/*/*/*h5")
   

    caObj = CalipsoAvhrrTrackObject()
    for filename in files:
        #print  os.path.basename(filename)
        caObj +=readCaliopAvhrrMatchObj(filename)#, var_to_skip='segment')


    #table_5_per_illumination(caObj)
    plot_ct_table4(caObj)
    plot_ct_table2(caObj)
    table_21_do_for_atbd(caObj)
    #plot_ct_table3(caObj)
