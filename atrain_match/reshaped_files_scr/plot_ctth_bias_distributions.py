"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readCaliopImagerMatchObj,
                            DataObject,
                            CloudsatImagerTrackObject,
                            readCloudsatImagerMatchObj,
                            CalipsoImagerTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       ppsMatch_Imager_CalipsoObject)
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 18})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc('serif', "Times New Roman")
#plt.rcParams["font.serif"] = "times new roman"
#from matplotlib import rc
print matplotlib.rcParams
#rc('text', usetex=True)
#rc('mathtext', usetex=True)
#rc('font', size=18)
#rc('text.latex', preamble=r'\usepackage{times}')
delta_h = 100.0


from utils.get_flag_info import get_calipso_clouds_of_type_i
from utils.get_flag_info import (get_semi_opaque_info_pps2014,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_low_clouds)

class PlotAndDataObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'height_c': None,
            'height_c2': None,
            'low_clouds': None,
            'medium_clouds': None,
            'high_clouds': None,
            'height_mlvl2': None,
            'satz': None,
            'height_old': None,
            'height_pps': None,
            'cflag': None,
            #'bias_nnant': None,
            'pressure_c': None,
            'height_nna1nt': None,
            'height_nnvnt': None,
            'height_nnm2nt': None,
            'height_nnmintnco2': None, 
            'height_nnmint': None,
            'pressure_mlvl2': None,
            'pressure_old': None,
            'pressure_pps': None,
            'pressure_nnant': None,
            'pressure_nna1nt': None,
            'pressure_nnvnt': None,
            'pressure_nnm2nt': None,
            'pressure_nnmintnco2': None, 
            'pressure_nnmint': None,
            'use': None,
            'use_all': None,
            'old_bias': None,
            'pps_bias': None,
            'mlvl2_bias': None,
            'bias_nnant': None,
            'bias_nna1nt': None,
            'bias_nnvnt': None,
            'bias_nnm2nt': None,
            'bias_nnmintnco2': None, 
            'bias_nnmint': None,
            'ok_old': None,
            'ok_nnant': None,
            'ok_nna1nt': None,
            'ok_nnvnt': None,
            'ok_nnm2nt': None,
            'ok_nnmintnco2': None, 
            'ok_nnmint': None,
            'ok_cma_pps': None, 
            'ok_cma_mlvl2': None ,
            'ok_modis': None, 
            'ok_cma':None
        }

def extract_data(cObj, sat='cloudsat'):
    pltObj = PlotAndDataObject()
    print("max_pressure %d"%(0.01*np.max(cObj.imager.all_arrays['psur'])))
    if sat.lower() in 'cloudsat':
        pltObj.height_c2 = cObj.cloudsat.all_arrays['clsat_max_height']
        clsat_max_height = -9 + 0*np.zeros(cObj.cloudsat.latitude.shape)
        #for i in range(125):
        #    height = cObj.cloudsat.Height[:,i]
        #    cmask_ok = cObj.cloudsat.CPR_Cloud_mask[:,i]
        #    top_height = height+120
        #    #top_height[height<240*4] = -9999 #Do not use not sure why these are not used Nina 20170317
        #    is_cloudy = cmask_ok > 30
        #    top_height[~is_cloudy] = -9999
        #    clsat_max_height[clsat_max_height<top_height] =  top_height[clsat_max_height<top_height]
        #height_c =  clsat_max_height  
        height_c = cObj.cloudsat.validation_height

        pltObj.low_clouds = np.logical_and(height_c<cObj.imager.all_arrays['segment_nwp_h680'], height_c>-9)
        pltObj.medium_clouds = np.logical_and(height_c>=cObj.imager.all_arrays['segment_nwp_h680'], 
                                              height_c<=cObj.imager.all_arrays['segment_nwp_h440'])
        pltObj.high_clouds = np.logical_and(height_c>cObj.imager.all_arrays['segment_nwp_h440'], height_c>-9)
        elevation = cObj.cloudsat.all_arrays['elevation']
        elevation[elevation<0] = 0
    elif  sat.lower() in 'calipso': 
        height_c = 1000*cObj.calipso.all_arrays['layer_top_altitude'][:,0]
        pltObj.pressure_c = cObj.calipso.all_arrays['layer_top_pressure'][:,0]
        pltObj.low_clouds = get_calipso_low_clouds(cObj)
        pltObj.high_clouds = get_calipso_high_clouds(cObj)
        pltObj.medium_clouds = get_calipso_medium_clouds(cObj)
        elevation = cObj.calipso.all_arrays['elevation']
        elevation[elevation<0] = 0
        ninas_od = cObj.calipso.all_arrays['total_optical_depth_5km']
        print "min optical depth", np.min(ninas_od[np.logical_and(height_c>0,ninas_od>0)])
        use_part = False
        part = "all"
        if use_part:
            part = "single_layer_sea_below_60_in5km_as_art"
            use_part = np.logical_and(
                use_part,cObj.calipso.all_arrays['number_layers_found']==1)
            use_part = np.logical_and(
                use_part,np.abs(cObj.calipso.all_arrays['latitude'])<60)
            use = np.logical_and(
                use_part,np.equal(cObj.calipso.igbp_surface_type,17))
            #use = np.logical_and(use_part,height_c<3000+cObj.calipso.all_arrays['elevation'])
            use_part = np.logical_and(use_part,cObj.calipso.all_arrays[
                'feature_optical_depth_532_top_layer_5km']>0)
            use_part = np.logical_and(use_part,
                                 cObj.calipso.all_arrays[
                                     'feature_optical_depth_532_top_layer_5km']==
                                 cObj.calipso.all_arrays['total_optical_depth_5km'])
            part = "single_layer_sea_below_60_in5km_all_od_top"
            low = np.logical_and(low_clouds,use_part)
            #low = np.logical_and(height_c<low_clouds,use_part)
            low = np.logical_and(
                use_part,height_c<(3000+cObj.calipso.all_arrays['elevation']))
            medium = np.logical_and(medium_clouds,use_part)
            high = np.logical_and(height_c>8000,use_part)

        
        pltObj.cflag =  cObj.calipso.feature_classification_flags[::,0]

        elevation = cObj.calipso.all_arrays['elevation']
        elevation[elevation<0] = 0
    pltObj.height_c = height_c
    pltObj.height_mlvl2 = cObj.modis.all_arrays['height']#+elevation #??
    pltObj.pressure_mlvl2 = cObj.modis.all_arrays['pressure']
    #height_pps = cObj.imager.all_arrays['imager_ctth_m_above_seasurface']
    pltObj.satz = cObj.imager.all_arrays['satz']
    pltObj.height_pps = cObj.imager.all_arrays['ctthnnant_height']+elevation
    pltObj.height_old = cObj.imager.all_arrays['ctthold_height']+elevation
    pltObj.pps_bias = pltObj.height_pps - height_c
    pltObj.old_bias = pltObj.height_old - height_c
    pltObj.mlvl2_bias = pltObj.height_mlvl2 - height_c
    use =  height_c>=0
    use = np.logical_and(use, pltObj.height_mlvl2>-1)
    use = np.logical_and(use, cObj.modis.all_arrays['pressure']>=70) #no_effekt
    use = np.logical_and(use,cObj.imager.all_arrays['ctthnnant_height']>-1)
    use = np.logical_and(use,cObj.imager.all_arrays['ctthold_height']>-1)
    #use = np.logical_and(use,cObj.imager.all_arrays['ctthnnant_height']<55000)no_effekt
    use = np.logical_and(use,cObj.imager.all_arrays['ctthold_pressure']>=70*100)
    use = np.logical_and(use,cObj.modis.all_arrays['cloud_emissivity']<=100)
    pltObj.use = use
    use_all = use.copy()
    for var in ['ctthnnant_height','ctthnna1nt_height','ctthnnvnt_height','ctthnnm2nt_height','ctthnnmintnco2_height', 'ctthnnmint_height', 'ctthold_height']:
        name = var.replace('_height', '').replace('ctth', 'pressure_')
        #print name
        pltObj.all_arrays[name] = 0.01*cObj.imager.all_arrays[var.replace('height','pressure')]
        update = pltObj.all_arrays[name] > 0.01*cObj.imager.all_arrays['psur']
        if 'old' not in name:
            update = pltObj.all_arrays[name] > 0.01*cObj.imager.all_arrays['psur']
            pltObj.all_arrays[name][update] = 0.01*cObj.imager.all_arrays['psur'][update]
            #ok_pressure = pltObj.all_arrays[name]>=70
            pressure_name = name
        name = var.replace('_height', '').replace('ctth', 'height_')
        pltObj.all_arrays[name] = cObj.imager.all_arrays[var] +elevation 
        name = var.replace('_height', '').replace('ctth', 'bias_')
        pltObj.all_arrays[name] = cObj.imager.all_arrays[var] +elevation - height_c
        #ok_ctth = np.logical_and(cObj.imager.all_arrays[var]<550000, cObj.imager.all_arrays[var]>-1)
        ok_ctth = cObj.imager.all_arrays[var]>-1
        ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['text_t11']>-1)
        use_all = np.logical_and(use_all, ok_ctth)
        if 'nnmi' in var:                  
            ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['modis_27']>-1)
        if 'nna1' in var:                  
            ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['bt37micron']>-1)
        if 'nnv' in var or 'nnm' in var:                  
            ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['bt86micron']>-1)
        if 'nna1' not in var:                  
            ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['bt12micron']>-1)
        if 'nnm' in var:                  
            ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['modis_28']>-1)
        if 'ctthnnmint_height' == var:                  
            ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['modis_33']>-1)
        ok_ctth = np.logical_and(ok_ctth, cObj.imager.all_arrays['bt11micron']>-1)
        use_all = np.logical_and(use_all, ok_ctth)
        # for pressure scatter plot:
        #use_all = np.logical_and(use_all, ok_pressure)

        
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['warmest_t12']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['coldest_t12']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['psur']>-9) 
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['surftemp']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['t950']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['t850']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['t700']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['t500']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['t250']>-9)
        
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['text_t11']>-9) #without this 1793146 pixels
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['text_t11t12']>-9) 
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['warmest_t11']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['coldest_t11']>-9)                                             
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['modis_27']>-1)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['modis_28']>-1)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['modis_33']>-1)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['text_t37']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['warmest_t37']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['coldest_t37']>-9)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['bt11micron']>-1)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['bt12micron']>-1)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['bt37micron']>-1)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['bt86micron']>-1)
        use_all = np.logical_and(use_all, cObj.imager.all_arrays['ciwv']>-9)
        
        pltObj.all_arrays[name.replace('bias_','ok_')] =  ok_ctth
        print "excluded", np.sum(pltObj.all_arrays[pressure_name][use_all]<70)

    pltObj.use_all = use_all
    pltObj.ok_cma_mlvl2 = cObj.modis.all_arrays['cloud_emissivity']<=100
    pltObj.ok_cma_pps = np.logical_or(cObj.imager.all_arrays['cloudmask']==1,
                                      cObj.imager.all_arrays['cloudmask']==2)
    pltObj.ok_cma_pps = np.logical_and(pltObj.ok_cma_pps, height_c>=0)
    pltObj.ok_cma_mlvl2 = np.logical_and(pltObj.ok_cma_mlvl2, height_c>=0)                   
    pltObj.ok_modis = np.logical_and(pltObj.height_mlvl2>-1, pltObj.height_mlvl2<45000)
    pltObj.ok_cma = np.logical_and( pltObj.ok_modis,pltObj.ok_cma_pps)

    return pltObj

def print_stats(pltObj, month, day_strm, sat='CLOUDSAT_OR_CALIPSO'):
    pps_bias = pltObj.pps_bias
    old_bias = pltObj.old_bias
    mlvl2_bias = pltObj.mlvl2_bias
    print sat.upper()
    use_all =  pltObj.use_all
    low_print = np.logical_and(use_all, pltObj.low_clouds)
    medium_print = np.logical_and(use_all, pltObj.medium_clouds)
    high_print = np.logical_and(use_all, pltObj.high_clouds)
    for bias_v, name in zip([old_bias, mlvl2_bias, pps_bias],
                            ["PPS-v2014", "MODIS-C6", "NN-AVHRR"]):
        print "%s & %3.0f & %3.0f & %3.0f & %3.0f \\\\"%(name, 
                                                         np.mean(np.abs(bias_v[use_all])), 
                                                         np.mean(np.abs(bias_v[low_print])),
                                                         np.mean(np.abs(bias_v[medium_print])),
                                                         np.mean(np.abs(bias_v[high_print])))
    for var_v, name in zip(['bias_nnant','bias_nna1nt','bias_nnvnt','bias_nnm2nt','bias_nnmintnco2', 'bias_nnmint'],
                            ["NN-AHRR" ,"NN-AHRR1", "NN-VIIRS", "NN-MERSI2",  "NN-MetImage-NoCO2", "NN-Metimage" ]):
        bias_v = pltObj.all_arrays[var_v]
        print "%s & %3.0f & %3.0f & %3.0f & %3.0f \\\\"%(name, 
                                                         np.mean(np.abs(bias_v[use_all])), 
                                                         np.mean(np.abs(bias_v[low_print])),
                                                         np.mean(np.abs(bias_v[medium_print])),
                                                         np.mean(np.abs(bias_v[high_print])))
    for var_v, name in zip(['ok_old', 'ok_nnant','ok_nna1nt','ok_nnvnt','ok_nnm2nt','ok_nnmintnco2', 'ok_nnmint'],
                            ["PPS-v2014", "NN-AVHRR" ,"NN-AHRR1", "NN-VIIRS", "NN-MERSI2",  "NN-MetImage-NoCO2", "NN-Metimage" ]):
        ok_ctth = pltObj.all_arrays[var_v] 
        ok_cma = pltObj.ok_cma
        n_cma = np.sum(ok_cma)
        n_ctth = np.sum(np.logical_and(ok_cma,ok_ctth) )                 
        print "%s & %3.0d & %3.0d & %3.2f \\\\"%(name, 
                                                 n_ctth,
                                                 n_cma,
                                                 100*n_ctth*1.0/n_cma)
    ok_ctth = pltObj.ok_modis
    ok_cma = pltObj.ok_cma
    n_cma = np.sum(ok_cma)
    n_ctth = np.sum(np.logical_and(ok_cma,ok_ctth) )      
    print "MODIS-C6 & %3.0d & %3.0d & %3.2f \\\\"%( n_ctth,
                                                   n_cma,
                                                   100*n_ctth*1.0/n_cma)                    
    #Number of
    n_all = len(bias_v[use_all])*1.0
    n_all_3 = len(bias_v[ pltObj.use])*1.0
    print "%s & %d %d & & %d & %d & %d \\\\"%(name, 
                                         n_all,    n_all_3, 
                                         len(bias_v[low_print]),
                                         len(bias_v[medium_print]),
                                         len(bias_v[high_print]))
    print "%s & %d & %3.1f & %3.1f & %3.1f \\\\"%(name, 
                                         n_all,
                                         len(bias_v[low_print])*100/n_all,
                                         len(bias_v[medium_print])*100/n_all,
                                         len(bias_v[high_print])*100/n_all)

def print_data_for_figure_2(pltObj, month, day_str, sat='calipso'):
    out_filename = "%s/data_fig1_%s.txt"%(MYPATH,month)
    if sat.lower() in ["cloudsat"]:
        out_file_h = open(out_filename,'w')
    else:
        out_file_h = open(out_filename,'a') 
    use =  pltObj.use
    low = np.logical_and(pltObj.low_clouds,use)
    medium = np.logical_and(pltObj.medium_clouds,use)
    high = np.logical_and(pltObj.high_clouds,use)
    pps_bias = pltObj.pps_bias
    old_bias = pltObj.old_bias
    mlvl2_bias = pltObj.mlvl2_bias   
    out_text = "" 
    def print_one(out_text, bias_v, selection,  label, cloudc):
        bmin = -20
        bmax = 20
        bins = np.arange(bmin*1000,bmax*1000,delta_h)
        hist_heights,bins = np.histogram(bias_v[selection],bins=bins)
        n_pix = np.sum(selection)
        hist_heights = hist_heights*100.0/n_pix
        out_text += "%s_%s_%s_x "%(label, cloudc, sat)
        formated_x = ["%.2f"%(0.001*(item+delta_h*0.5)) for item in bins[0:-1]]
        formated_x = " ".join(formated_x)
        out_text +=  formated_x + "\n"
        out_text += "%s_%s_%s_y "%(label, cloudc, sat)
        formated_y = ["%.2f"%(item) for item in  hist_heights]
        formated_y = " ".join(formated_y)
        out_text +=  formated_y + "\n"
        out_text +=  "%s_%s_%s_bias %3.4f \n"%(label, cloudc, sat, 0.001*np.mean(bias_v[selection]))
        out_text +=  "%s_%s_%s_median %3.4f \n"%(label, cloudc, sat, 0.001*np.median(bias_v[selection]))
        out_text +=  "%s_%s_%s_std %3.4f \n"%(label, cloudc, sat, 0.001*np.std(bias_v[selection]))
        return out_text  
   
    out_text = print_one(out_text, pps_bias, high,   "NN-AVHRR", "High")
    out_text = print_one(out_text, mlvl2_bias, high,   "MODIS-C6", "High")
    out_text = print_one(out_text, old_bias, high,   "PPS-v2014", "High")
    out_text = print_one(out_text, pps_bias, medium, "NN-AVHRR", "Medium")
    out_text = print_one(out_text, mlvl2_bias, medium, "MODIS-C6", "Medium")
    out_text = print_one(out_text, old_bias, medium,   "PPS-v2014", "Medium")
    out_text = print_one(out_text, pps_bias, low,  "NN-AVHRR", "Low",)
    out_text = print_one(out_text, mlvl2_bias, low,   "MODIS-C6", "Low")
    out_text = print_one(out_text, old_bias, low,   "PPS-v2014", "Low") 
    out_text = print_one(out_text, pps_bias, use,  "NN-AVHRR", "All",)
    out_text = print_one(out_text, mlvl2_bias, use,   "MODIS-C6", "All")
    out_text = print_one(out_text, old_bias, use,   "PPS-v2014", "All") 
    out_file_h.write(out_text)
    out_file_h.write("\n") 

def replot_figure2_from_saved_data(month):
    out_filename = "%s/data_fig1_%s.txt"%(MYPATH,month)
    in_file_h = open(filename,'r')
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    text_d = {'cloudsat_High':   "(a) High clouds",
              'cloudsat_Medium': "(c) Medium clouds",
              'cloudsat_Low':    "(e) Low clouds",
              'calipso_High':    "(b) High clouds",
              'calipso_Medium':  "(d) Medium clouds",
              'calipso_Low':     "(f) Low clouds"}

    def plot_one(ax, x_data, y_data, sat, cloudc, legend=False):
        text_i = text_d[sat+"_"+cloudc]
        print_sat = sat.upper()
        if print_sat in ["CLOUDSAT"]:
            print_sat = "CPR (CloudSat)"
        if print_sat in ["CALIPSO"]:
            print_sat = "CALIOP"
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.plot(x_data[sat][cloudc]["PPS-v2014"], y_data[sat][cloudc]["PPS-v2014"],
                 "-b",label="PPS-v2014") 
        plt.plot(x_data[sat][cloudc]["MODIS-C6"], y_data[sat][cloudc]["MODIS-C6"],
                 "-k",label="MODIS-C6")
        plt.plot(x_data[sat][cloudc]["NN-AVHRR"], y_data[sat][cloudc]["NN-AVHRR"],
                 "-r",label="NN-AVHRR")
        
        if cloudc in ["Low"]:
            plt.xlabel(" Retrieved height - %s (km) "%(print_sat))
        if sat in ["cloudsat"]:    
            plt.ylabel("Percent of data")
        if legend:
            pass
            plt.legend(fancybox=True, loc=1,  numpoints=3, bbox_to_anchor=(1.25, 1.1), framealpha=1.0,edgecolor='w', fontsize=18)
        ax.grid(True)
        plt.text(0.024, 0.90, text_i, fontsize=18,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))

    x_data = {}
    y_data = {}
    bias_data = {}
    median_data = {}
    for sat in ["calipso", "cloudsat"]:
        x_data[sat] = {}
        y_data[sat] = {}
        bias_data[sat] = {}
        median_data[sat] = {}
        for cloudc in ["All", "Low", "Medium",  "High"]:
            y_data[sat][cloudc] = {}
            x_data[sat][cloudc] = {}
            bias_data[sat][cloudc] = {}
            median_data[sat][cloudc] = {}

    fig = plt.figure(figsize=(9, 11))
    plt.subplots_adjust(wspace=0, hspace=0)
    for line in in_file_h:
        line = line.rstrip()
        if line == "":
            continue
        data = line.split(" ")
        info = data.pop(0)
        info_list = info.split("_")
        #print info_list
        algorithm = info_list[0]
        cloudc = info_list[1]        
        sat = info_list[2]
        if '_x ' in line:
            x_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])
        if '_y' in line:
            y_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])
        if '_bias' in line:
            bias_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])[0]
        if '_median' in line:
            median_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])[0]


  
    ax = fig.add_subplot(321)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xticks(np.arange(-12,12,2.0))
    for label in ax.xaxis.get_ticklabels()[::]:
        label.set_visible(False)
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(True)
    plt.yticks(np.arange(0,4,0.5))                      
    ax.set_xlim(-11.5,11.5)
    ax.set_ylim(0,4.0)
    plot_one(ax, x_data, y_data,  "cloudsat", "High")
    ax = fig.add_subplot(322)
    plt.xticks(np.arange(-12,12,2.0))
    for label in ax.xaxis.get_ticklabels()[::]:
        label.set_visible(False)
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(True)

    plt.yticks(np.arange(0,4,0.5))                      
    ax.set_xlim(-11.5,11.5)
    ax.set_ylim(0,4.0)
    plot_one(ax, x_data, y_data,  "calipso", "High", True)
    ax = fig.add_subplot(323)
    plt.xticks(np.arange(-8,8,2.0))
    plt.yticks(np.arange(0,12,1.0))
    ax.set_xlim(-7,7)
    ax.set_ylim(0,6.5)
    plot_one(ax, x_data, y_data,  "cloudsat", "Medium")
    ax = fig.add_subplot(324)
    plt.xticks(np.arange(-8,8,2.0))
    plt.yticks(np.arange(0,12,1.0))
    ax.set_xlim(-7,7)
    ax.set_ylim(0,6.5)
    plot_one(ax, x_data, y_data, "calipso", "Medium")
    ax = fig.add_subplot(325)
    plt.xticks(np.arange(-5,5,1.0))
    plt.yticks(np.arange(0,14,2.0))
    ax.set_xlim(-4.5,4.5)
    ax.set_ylim(0,12.0)
    plot_one(ax, x_data, y_data, "cloudsat", "Low")
    ax = fig.add_subplot(326)
    plt.xticks(np.arange(-5,5,1.0))
    plt.yticks(np.arange(0,14,2.0))
    ax.set_xlim(-4.5,4.5)
    ax.set_ylim(0,12.0)
    plot_one(ax, x_data, y_data,  "calipso", "Low")

    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig02_%s.pdf"%(month), bbox_inches='tight')
    plt.close("all")


def plot_medium_bias_plot_from_saved_data(month):
    out_filename = "%s/data_fig1_%s.txt"%(MYPATH,month)
    filename = "data_fig1_%s.txt"%(month)
    in_file_h = open(filename,'r')
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    text_d = { "(a) High clouds",
              "(c) Medium clouds",
              "(e) Low clouds",
              "(b) High clouds",
                "(d) Medium clouds",
               "(f) Low clouds"}

    def plot_one(ax, x_data, y_data, median_data, bias_data, std_data, sat, cloudc, imager, text_i, maxx, maxy, legend=False):
        #text_i = text_d[sat+"_"+imager]
        print_sat = sat.upper()
        if print_sat in ["CLOUDSAT"]:
            print_sat = "CPR (CloudSat)"
        if print_sat in ["CALIPSO"]:
            print_sat = "CALIOP"
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        color1 = "-k"
        color2 = ":k"
        color3 = "--k"
        if imager in ["NN-AVHRR"]:
            color1 = "-r"
            color2 = ":r"
            color3 = "--r"
        if imager in ["PPS-v2014"]:
            color1 = "-b"
            color2 = ":b"
            color3 = "--b"

        n_pix =  10000000   
        temp_data = np.random.normal(1000*bias_data[sat][cloudc][imager],  1000*std_data[sat][cloudc][imager], n_pix)
        bmin = -20
        bmax = 20
        bins = np.arange(bmin*1000,bmax*1000,delta_h)
        hist_heights,bins = np.histogram(temp_data,bins=bins)
        hist_heights = hist_heights*100.0/n_pix

        ax.fill(x_data[sat][cloudc][imager], hist_heights, color='silver',  label='Gaussian')
        #plt.plot(x_data[sat][cloudc][imager], hist_heights, 0.1,  label='Gaussian')
        plt.plot(x_data[sat][cloudc][imager], y_data[sat][cloudc][imager],
                 color1,label=imager)
        yind_bias = np.argmin(np.abs(x_data[sat][cloudc][imager] - bias_data[sat][cloudc][imager]))
        ylim_bias = y_data[sat][cloudc][imager][yind_bias]
        yind_median = np.argmin(np.abs(x_data[sat][cloudc][imager] - median_data[sat][cloudc][imager]))
        ylim_median = y_data[sat][cloudc][imager][yind_median]
        plt.plot([median_data[sat][cloudc][imager],median_data[sat][cloudc][imager]], [0, ylim_median], color3, label="median")
        plt.plot([bias_data[sat][cloudc][imager],bias_data[sat][cloudc][imager]],[0, ylim_bias], color2, label="bias")
        #print bias_data[sat][cloudc]["NN-AVHRR"]
        #plt.show()
        #plt.title("Non Gaussian error distributions ")
        if imager in ["NN-AVHRR"]:
            plt.xlabel(" Retrieved height - %s (km) "%(print_sat))
        if sat in ["cloudsat"]:    
            plt.ylabel("Percent of data")
        if legend:
            pass
            plt.legend(fancybox=True, loc=1,  numpoints=3, bbox_to_anchor=(1.30, 1.1), framealpha=1.0,edgecolor='w', fontsize=18)
        ax.grid(True)
        plt.text(0.025, 0.90, text_i, fontsize=18,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))        
        plt.xticks(np.arange(-12,12,1.0))
        if cloudc in ["High", "Medium"]:
            plt.xticks(np.arange(-12,12,2.0))
        for label in ax.xaxis.get_ticklabels()[::]:
            label.set_visible(False)
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(True)
        ax.set_xlim(-maxx,maxx)
        ax.set_ylim(0, maxy)    

    x_data = {}
    y_data = {}
    bias_data = {}
    median_data = {}
    std_data = {}
    for sat in ["calipso", "cloudsat"]:
        x_data[sat] = {}
        y_data[sat] = {}
        bias_data[sat] = {}
        median_data[sat] = {}
        std_data[sat] = {}
        for cloudc in ["All", "Low", "Medium",  "High"]:
            y_data[sat][cloudc] = {}
            x_data[sat][cloudc] = {}
            bias_data[sat][cloudc] = {}
            median_data[sat][cloudc] = {}
            std_data[sat][cloudc] = {}

    for line in in_file_h:
        line = line.rstrip()
        if line == "":
            continue
        data = line.split(" ")
        info = data.pop(0)
        info_list = info.split("_")
        #print info_list
        algorithm = info_list[0]
        cloudc = info_list[1]        
        print cloudc
        sat = info_list[2]
        if '_x ' in line:
            x_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])
        if '_y' in line:
            y_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])
        if '_bias' in line:
            bias_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])[0]
            print bias_data[sat][cloudc][algorithm]
        if '_median' in line:
            median_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])[0]
        if '_std' in line:
            std_data[sat][cloudc][algorithm] = np.array([np.float(s) for s in data])[0]
    for cloudc, maxx, maxy, in zip(["All", "Low", "Medium", "High"],
                                   [4.5,4.5,7.0,8.5],
                                   [6.5,12.0,6.5,4.0]):        
        fig = plt.figure(figsize=(9, 11))
        if cloudc not in ["All"]:
            plt.suptitle(cloudc)
        ax = fig.add_subplot(321)
        plot_one(ax, x_data, y_data, median_data, bias_data, std_data, "cloudsat", cloudc, "PPS-v2014", " (a) ", maxx, maxy)#, True)
        ax = fig.add_subplot(323)
        plot_one(ax, x_data, y_data, median_data, bias_data, std_data,  "cloudsat", cloudc, "MODIS-C6", " (c) ", maxx, maxy)#,True)
        ax = fig.add_subplot(325)
        plot_one(ax, x_data, y_data, median_data, bias_data, std_data,  "cloudsat", cloudc, "NN-AVHRR", " (e) ", maxx, maxy)#,True)

        #plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_PLOTS/figxx_cloudsat_bias_median_%s.pdf"%(month), bbox_inches='tight')
        #plt.close("all")
        #fig = plt.figure(figsize=(9, 11))    
        ax = fig.add_subplot(322)
        plot_one(ax, x_data, y_data, median_data, bias_data, std_data,  "calipso", cloudc, "PPS-v2014", " (b) ", maxx, maxy,True)
        ax = fig.add_subplot(324)
        plot_one(ax, x_data, y_data, median_data, bias_data, std_data,  "calipso", cloudc, "MODIS-C6", " (d) ", maxx, maxy,True)
        ax = fig.add_subplot(326)
        plot_one(ax, x_data, y_data, median_data, bias_data,  std_data, "calipso", cloudc, "NN-AVHRR", " (f) ", maxx, maxy,True)

        plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig08_bias_median_%s_%s.pdf"%(cloudc, month), bbox_inches='tight')
        plt.close("all")
        fig = plt.figure(figsize=(9, 11))
        ax = fig.add_subplot(111)
        plot_one(ax, x_data, y_data, median_data, bias_data, std_data,  "cloudsat", cloudc, "NN-AVHRR", " (e) ", maxx, maxy,True)
        plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig08_bias_median_%s_%s_e.pdf"%(cloudc, month), bbox_inches='tight')

def make_profileplot(pltObj, month, day_str, sat='calipso'):
    print_stats(pltObj, month, day_str, sat=sat)
    print_data_for_figure_2(pltObj, month, day_str, sat=sat)
    use =  pltObj.use
    low = np.logical_and(pltObj.low_clouds,use)
    medium = np.logical_and(pltObj.medium_clouds,use)
    high = np.logical_and(pltObj.high_clouds,use)
    pps_bias = pltObj.pps_bias
    old_bias = pltObj.old_bias
    mlvl2_bias = pltObj.mlvl2_bias
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    delta_h = 100.0   
    text1 = "(b) High clouds"
    text2 = "(d) Medium clouds"
    text3 = "(f) Low clouds"
    print_sat = sat.upper()
    if print_sat in ["CLOUDSAT"]:
        print_sat = "CPR (CloudSat)"
        text1 = text1.replace("(b)", "(a)")
        text2 = text2.replace("(d)", "(c)")
        text3 = text3.replace("(f)", "(e)")
    def plot_one(bias_v, selection, bmin, bmax, color, label):
       bins = np.arange(bmin*1000,bmax*1000,delta_h)
       hist_heights,bins = np.histogram(bias_v[selection],bins=bins)
       n_pix = np.sum(selection)
       hist_heights = hist_heights*100.0/n_pix
       plt.plot(0.001*(bins[0:-1]+delta_h*0.5), hist_heights,
                color,label=label) 
       plt.xlabel(" Retrieved height - %s (km) "%(print_sat))
       plt.ylabel("Percent of data")
       ax.set_xlim(bmin,bmax)
    fig = plt.figure(figsize = (5.5,12))        
    ax = fig.add_subplot(313)
    #plt.title("Low")
    plt.text(0.02, 0.90, text3, fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-5,5,1.0))
    plt.yticks(np.arange(0,14,1.0))    
    plot_one(pps_bias, low, -14,14, '-r', "NN-AVHRR")
    plot_one(mlvl2_bias, low, -14,14, '-k', "MODIS-C6")
    plot_one(old_bias, low, -14,14, '-b', "PPS-v2014")
    #plot_one(nna1_bias, low, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-5,5)
    ax.set_ylim(0,12.0)
    ax.grid(True)
    ax = fig.add_subplot(312)
    #plt.title("Medium")
    plt.text(0.02, 0.90, text2, fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-8,8,2.0))
    plt.yticks(np.arange(0,12,1.0))
    plot_one(pps_bias, medium, -12,12, '-r', "NN-AVHRR")
    plot_one(mlvl2_bias, medium, -12,12, '-k', "MODIS-C6")
    plot_one(old_bias, medium, -12,12, '-b', "PPS-v2014")
    #plot_one(nna1_bias, medium, -12,12, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)
    ax.set_xlim(-7,7)
    ax.set_ylim(0,6.0)
    ax.grid(True)
    ax = fig.add_subplot(311)
    #plt.title("High")
    plt.text(0.02, 0.90, text1, fontsize=14,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
    plt.xticks(np.arange(-12,12,2.0))
    plt.yticks(np.arange(0,4,0.5))
    plot_one(pps_bias, high, -14,14, '-r', "NN-AVHRR")
    plot_one(mlvl2_bias, high, -14,14, '-k', "MODIS-C6")
    plot_one(old_bias, high, -14,14, '-b', "PPS-v2014")
    #plot_one(nna1_bias, high, -14,14, '0.3', "NN-AVHRR1")
    plt.legend(fancybox=True, loc=1,  numpoints=4)                        
    ax.set_xlim(-12,12)
    ax.set_ylim(0,4.0)
    ax.grid(True)
    #plt.show()   
    #plt.title("%s MAE = %3.0f"%(name,MAE))
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX_%s/ctth_bias_profile_%s_%s_%s.png"%(sat, day_str, month,sat))
    plt.savefig("/home/a001865/PICTURES_FROM_PYTHON/CTTH_BOX_%s/ctth_bias_profile_%s_%s_%s.pdf"%(sat, day_str, month,sat))
    plt.close("all")

def investigate_nn_ctth_modis_lvl2_cloudsat():
    day_str="01st"
    ROOT_DIR = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        "global_modis_%s_created20180316/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
        #"global_modis_%s_created20170330/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
    #clsatObj = CloudsatImagerTrackObject()
    plt_obj = PlotAndDataObject()
    name=""
    for month in [ "02", "04","06", "08","10", "12"]:#[ "03","05", "07","09", "11"]: #[ "02", "04","06", "08","10", "12"]:#, "03","05", "07","09", "11","01" ]:  #[ "06", "09"]: 
        print ROOT_DIR%(day_str,month)
        files = glob(ROOT_DIR%(day_str,month))
        name+=month 
        plt_obj_new = PlotAndDataObject()
        for filename in files:
            print filename
            clsatObj_new =  readCloudsatImagerMatchObj(filename)
            plt_obj_new += extract_data(clsatObj_new, sat='cloudsat')
        plt_obj +=  plt_obj_new
        make_profileplot(plt_obj_new, month=month,day_str=day_str, sat='cloudsat')
    make_profileplot(plt_obj, month=name,day_str=day_str, sat='cloudsat')

def investigate_nn_ctth_modis_lvl2():
    day_str="01st"
    ROOT_DIR = (
        "/home/a001865/DATA_MISC/reshaped_files/"
        #"global_modis_%s_created20170504/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")
        "global_modis_%s_created20180316/Reshaped_Files_merged_calipso_cbase/eos2/1km/2010/%s/*h5")
    plt_obj = PlotAndDataObject()
    name=""
    for month in [  "02", "04","06", "08","10","12"]:# ["01", "03","05", "07","09", "11"]: # #[ "06", "09"]:    
        print ROOT_DIR%(day_str,month)
        files = glob(ROOT_DIR%(day_str,month))
        name+=month 
        plt_obj_new = PlotAndDataObject()
        for filename in files:
            #print filename
            caObj_new = readCaliopImagerMatchObj(filename)  
            plt_obj_new += extract_data(caObj_new, sat='calipso')
        plt_obj +=  plt_obj_new
        make_profileplot(plt_obj_new, month=month,day_str=day_str, sat='calipso')
    make_profileplot(plt_obj, month=name,day_str=day_str, sat='calipso')

        
if __name__ == "__main__":
    for month in ["02", "04","06", "08","10","12", "020406081012"]:
        replot_figure2_from_saved_data(month)
        plot_medium_bias_plot_from_saved_data(month)
        print "hej"

   
    #investigate_nn_ctth_modis_lvl2_cloudsat()
    #investigate_nn_ctth_modis_lvl2()


