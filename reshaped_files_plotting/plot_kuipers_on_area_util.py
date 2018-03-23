import numpy as np
from pyresample import utils
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import get_neighbour_info
from pyresample.kd_tree import get_sample_from_neighbour_info
import pyresample as pr
import os
from scipy import ndimage
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use("TkAgg")
from matplotlib import rc
print matplotlib.rcParams
#matplotlib.rcParams.update({'font.size': 30})
#matplotlib.rcParams.update({'font.family': 'sans-serif'})
#matplotlib.rcParams.update({'font.sans-serif': ['Verdana']})
#rc('font', size=30)
#rc('text.latex', preamble=r'\usepackage{times}')
#rc('text', usetex=True)

matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 30})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

print matplotlib.rcParams
#plt.rc('font', family='sans serif')
#plt.rcParams["font.family"] = "Verdana"

from matchobject_io import DataObject        
class ppsMatch_Imager_CalipsoObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'detected_clouds': None,
            'undetected_clouds': None,
            'false_clouds': None,
            'detected_clear': None,
            'new_detected_clouds': None,
            'new_false_clouds': None,
            'lats': None,
            'lons': None}
    def get_some_info_from_caobj(self, caObj):
        self.set_false_and_missed_cloudy_and_clear(caObj=caObj)
        self.set_r13_extratest(caObj=caObj)
        self.get_thr_offset(caObj=caObj)
        self.get_lapse_rate(caObj=caObj)
        self.get_ctth_bias(caObj=caObj)
        self.get_ctth_bias_low(caObj=caObj)
        self.get_ctth_bias_high(caObj=caObj)
        self.get_ctth_bias_low_temperature(caObj=caObj)
        self.height_bias_type={}
        self.detected_height_type={}
        for cc_type in xrange(8):
            self.get_ctth_bias_type(caObj=caObj, calipso_cloudtype=cc_type)

    def set_false_and_missed_cloudy_and_clear(self, caObj):
        lat = caObj.avhrr.all_arrays['latitude']
        lon = caObj.avhrr.all_arrays['longitude']
        if caObj.avhrr.all_arrays['cloudmask'] is not None:
            isCloudyPPS = np.logical_or(caObj.avhrr.all_arrays['cloudmask']==1,
                                        caObj.avhrr.all_arrays['cloudmask']==2) 
            isClearPPS = np.logical_or(caObj.avhrr.all_arrays['cloudmask']==0,
                                       caObj.avhrr.all_arrays['cloudmask']==3)
        else:
            isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                                         caObj.avhrr.all_arrays['cloudtype']<21) 
            isClearPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                        caObj.avhrr.all_arrays['cloudtype']<5)


        isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                                     caObj.avhrr.all_arrays['cloudtype']<21) 
        isClearPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                    caObj.avhrr.all_arrays['cloudtype']<5)
        nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
        meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
        if self.cc_method == 'BASIC' and self.isGAC:
            isCalipsoCloudy =  caObj.calipso.all_arrays['cloud_fraction']>0.5
            isCalipsoClear = np.not_equal(isCalipsoCloudy, True)
        elif self.cc_method == 'KG' and self.isGAC:
            isCalipsoCloudy = np.logical_and(
                caObj.calipso.all_arrays['cloud_fraction']>0.5,
                caObj.calipso.all_arrays['total_optical_depth_5km']>0.15)
            isCalipsoClear = np.not_equal(isCalipsoCloudy, True)
        elif self.cc_method == 'Nina' and self.isGAC:    
            isCalipsoCloudy = np.logical_and(
                nlay > 0, 
                caObj.calipso.all_arrays['cloud_fraction']>0.5)
            isCalipsoCloudy = np.logical_and(
                isCalipsoCloudy, 
                caObj.calipso.all_arrays['total_optical_depth_5km']>0.15)
            isCalipsoClear = np.logical_and(nlay == 0, meancl<0.01)
            isCalipsoClear = np.logical_and(
                isCalipsoClear, 
                caObj.calipso.all_arrays['total_optical_depth_5km']<0)
        elif self.cc_method == 'KG':
            isCalipsoCloudy = nlay>0
            isCalipsoClear = np.not_equal(isCalipsoCloudy, True)   
        elif self.cc_method == 'Nina':
            isCalipsoCloudy = np.logical_or(
                caObj.calipso.all_arrays['total_optical_depth_5km']>0.15, 
                np.logical_and(caObj.calipso.all_arrays['total_optical_depth_5km']<0,
                               nlay>0))
            isCalipsoClear = np.logical_and(nlay == 0, meancl<0.01)
            isCalipsoClear = np.logical_and(
                isCalipsoClear,
                caObj.calipso.all_arrays['total_optical_depth_5km']<0)
        use_ok_lat_lon = np.logical_and(np.logical_and(lat>=-90, lat<=90),
                                        np.logical_and(lon>=-180, lat<=180))
        isCloudyPPS = np.logical_and(isCloudyPPS,  use_ok_lat_lon)
        isClearPPS =  np.logical_and(isClearPPS,  use_ok_lat_lon)
        sunz = caObj.avhrr.all_arrays['sunz']
        if self.filter_method == "satz":
            satz = caObj.avhrr.all_arrays['satz']
            isCloudyPPS = np.logical_and(isCloudyPPS,  satz<=30)
            isClearPPS =  np.logical_and(isClearPPS,  satz<=30) 
        if self.DNT in ["day"]:
            isCloudyPPS = np.logical_and(isCloudyPPS,  sunz<=80)
            isClearPPS =  np.logical_and(isClearPPS,  sunz<=80)
        if self.DNT in ["night"]:
            isCloudyPPS = np.logical_and(isCloudyPPS,  sunz>=95)
            isClearPPS =  np.logical_and(isClearPPS,  sunz>=95)
        if self.DNT in ["twilight"]:
            isCloudyPPS = np.logical_and(isCloudyPPS,  sunz>80)
            isClearPPS =  np.logical_and(isClearPPS,  sunz>80)
            isCloudyPPS = np.logical_and(isCloudyPPS,  sunz<95)
            isClearPPS =  np.logical_and(isClearPPS,  sunz<95)                  
        if self.DNT in ["all"]:
            pass
        undetected_clouds = np.logical_and(isCalipsoCloudy, isClearPPS)
        false_clouds = np.logical_and(isCalipsoClear, isCloudyPPS)
        detected_clouds = np.logical_and(isCalipsoCloudy, isCloudyPPS)
        detected_clear = np.logical_and(isCalipsoClear, isClearPPS)
        use = np.logical_or(np.logical_or(detected_clouds, detected_clear),
                            np.logical_or(false_clouds, undetected_clouds))
        detected_height = np.logical_and(detected_clouds,
                                         caObj.avhrr.all_arrays['ctth_height']>-9)
        detected_height = np.logical_and(detected_height,
                                         caObj.avhrr.all_arrays['ctth_height']<45000)

        #settings for article eos_modis:
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['warmest_t12']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['coldest_t12']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['psur']>-9) 
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['surftemp']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['t950']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['t850']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['t700']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['t500']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['t250']>-9)        
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['text_t11']>-9) #without this 1793146 pixels
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['text_t11t12']>-9) 
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['warmest_t11']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['coldest_t11']>-9)                                             
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['modis_27']>-1)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['modis_28']>-1)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['modis_33']>-1)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['text_t37']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['warmest_t37']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['coldest_t37']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['bt11micron']>-1)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['bt12micron']>-1)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['bt37micron']>-1)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['bt86micron']>-1)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays['ciwv']>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays["ctthold_height"]>-9)
        detected_height = np.logical_and(detected_height, caObj.avhrr.all_arrays["ctthnnant_height"]>-9)
        detected_height = np.logical_and(detected_height, caObj.modis.all_arrays["height"]>-9)


        detected_height = np.logical_and(detected_height,
                                         caObj.calipso.all_arrays['layer_top_altitude'][:,0]>-1)
        detected_temperature = np.logical_and(detected_clouds,
                                       caObj.avhrr.all_arrays['ctth_temperature']>-9)

        self.false_clouds = false_clouds[use]
        self.detected_clouds = detected_clouds[use]
        self.undetected_clouds = undetected_clouds[use]
        self.detected_clear = detected_clear[use]
        self.latitude = caObj.avhrr.latitude[use]
        self.longitude = caObj.avhrr.longitude[use] 
        self.use = use
        self.detected_height = detected_height[use]
        self.detected_temperature = detected_temperature[use]

    def set_r13_extratest(self,caObj):
        if np.size(caObj.avhrr.all_arrays['r13micron'])==1 and caObj.avhrr.all_arrays['r13micron'] is None:
            self.new_false_clouds = np.zeros(self.false_clouds.shape)
            self.new_detected_clouds = np.zeros(self.false_clouds.shape)  
            return
        r13 = caObj.avhrr.all_arrays['r13micron']
        sunz =  caObj.avhrr.all_arrays['sunz']
        sunz_cos = sunz.copy()
        sunz_cos[sunz>87] =87    
        r13[sunz<90] = r13[sunz<90]/np.cos(np.radians(sunz_cos[sunz<90]))
        isCloud_r13 = np.logical_and(r13>2.0, caObj.avhrr.all_arrays['ciwv']>3)
        new_detected_clouds = np.logical_and(self.detected_clouds,
                                             isCloud_r13[self.use])
        new_false_clouds = np.logical_and(self.detected_clear,
                                          isCloud_r13[self.use])
        self.new_false_clouds = new_false_clouds
        self.new_detected_clouds = new_detected_clouds

    def get_thr_offset(self, caObj):
        if np.size(caObj.avhrr.all_arrays['surftemp'])==1 and caObj.avhrr.all_arrays['surftemp'] is None:
            self.t11ts_offset = np.zeros(self.false_clouds.shape)
            return
        t11ts_offset = (caObj.avhrr.all_arrays['bt11micron'] - 
                        caObj.avhrr.all_arrays['surftemp'] - 
                        caObj.avhrr.all_arrays['thr_t11ts_inv'])
        t11t12_offset = (caObj.avhrr.all_arrays['bt11micron'] - 
                        caObj.avhrr.all_arrays['bt12micron'] - 
                        caObj.avhrr.all_arrays['thr_t11t12'])
        t37t12_offset = (caObj.avhrr.all_arrays['bt37micron'] - 
                        caObj.avhrr.all_arrays['bt12micron'] - 
                        caObj.avhrr.all_arrays['thr_t37t12'])
        t11t37_offset = (caObj.avhrr.all_arrays['bt11micron'] - 
                        caObj.avhrr.all_arrays['bt37micron'] - 
                        caObj.avhrr.all_arrays['thr_t11t37'])
        t11ts_offset = t11ts_offset[self.use]
        t11ts_offset[self.detected_clouds] = 99999
        t11ts_offset[self.undetected_clouds] = 99999
        self.t11ts_offset = t11ts_offset
        t11t12_offset = t11t12_offset[self.use]
        t11t12_offset[self.detected_clouds] = -99999
        t11t12_offset[self.undetected_clouds] = -99999
        self.t11t12_offset = t11t12_offset
        t37t12_offset = t37t12_offset[self.use]
        t37t12_offset[self.detected_clouds] = -99999
        t37t12_offset[self.undetected_clouds] = -99999
        self.t37t12_offset = t37t12_offset
        t11t37_offset = t11t37_offset[self.use]
        t11t37_offset[self.detected_clouds] = -99999
        t11t37_offset[self.undetected_clouds] = -99999
        self.t11t37_offset = t11t37_offset
    def get_lapse_rate(self, caObj):
        if np.size(caObj.avhrr.all_arrays['surftemp'])==1 and caObj.avhrr.all_arrays['surftemp'] is None:
            self.lapse_rate = np.zeros(self.false_clouds.shape)
            return
        from get_flag_info import get_calipso_low_clouds
        low_clouds = get_calipso_low_clouds(caObj)
        delta_h = caObj.calipso.all_arrays['layer_top_altitude'][:,0] - 0.001*caObj.calipso.all_arrays['elevation'][:]
        delta_t = (273.15 + caObj.calipso.all_arrays['layer_top_temperature'][:,0] - caObj.avhrr.all_arrays['surftemp'][:])
        lapse_rate = delta_t/delta_h
        lapse_rate[caObj.calipso.all_arrays['layer_top_temperature'][:,0]<-500] = 0
        lapse_rate[caObj.calipso.all_arrays['layer_top_altitude'][:,0]>35.0] = 0
        lapse_rate[~low_clouds] = 0.0
        self.lapse_rate = lapse_rate[self.use]

    def get_ctth_bias(self, caObj):
        height_c = (1000*(
            caObj.calipso.all_arrays['layer_top_altitude'][self.use,0]) - 
            caObj.calipso.all_arrays['elevation'][self.use])
        if "cci" in self.satellites:
            height_c = 1000*caObj.calipso.all_arrays[
                'layer_top_altitude'][self.use,0]   
        if "modis_lvl2" in self.satellites:
            height_c = 1000*caObj.calipso.all_arrays[
                'layer_top_altitude'][self.use,0]  

        height_pps = caObj.avhrr.all_arrays['ctth_height'][self.use]
        delta_h = height_pps - height_c
        self.height_bias = delta_h
        self.height_bias[~self.detected_height]=0
        
        self.detected_height_both = np.where(self.detected_height, False,self.detected_height)
        self.height_mae_diff = 0*self.height_bias
        if caObj.modis.all_arrays["height"] is not None:
            height_modis = caObj.modis.all_arrays["height"][self.use]   
            delta_h_modis = height_modis - 1000*caObj.calipso.all_arrays[
                'layer_top_altitude'][self.use,0] 
            self.detected_height_both = np.logical_and(self.detected_height,height_modis>0) 
            mae_pps = np.abs(self.height_bias.copy())
            mae_modis = np.abs(delta_h_modis)
            diff_mae = mae_modis - mae_pps 
            diff_mae[~self.detected_height_both] = 0
        
            self.height_mae_diff = diff_mae

        try:
            #tsur = caObj.avhrr.all_arrays['surftemp']
            tsur = caObj.avhrr.all_arrays['segment_nwp_temp'][:,0] 
            #temperature_pps = caObj.avhrr.all_arrays['ctth_temperature']
            temperature_pps = caObj.avhrr.all_arrays['bt12micron']
            temp_diff = temperature_pps - tsur
            rate_neg = -1.0/6.5
            rate_pos = +1.0/3.0
            new_pps_h = rate_neg*temp_diff*1000 
            new_pps_h[temp_diff>0] = rate_pos*temp_diff[temp_diff>0]*1000 
            new_pps_h[new_pps_h<100] = 100
            keep = (
                (tsur - caObj.avhrr.all_arrays['ctth_temperature'])/
                (tsur - caObj.avhrr.all_arrays['ttro']))>0.33
            keep = new_pps_h>3000
            new_pps_h[keep] = caObj.avhrr.all_arrays['ctth_height'][keep]
            self.lapse_bias = new_pps_h[self.use] - height_c
            self.lapse_bias[~self.detected_height]=0
            self.lapse_bias[temperature_pps<0]=0
        except:
            #raise
            self.lapse_bias = 0* height_c

    def get_ctth_bias_low(self, caObj):
        from get_flag_info import get_calipso_low_clouds
        low_clouds = get_calipso_low_clouds(caObj)
        detected_low = np.logical_and(self.detected_height, 
                                      low_clouds[self.use])
        delta_h = self.height_bias.copy()
        delta_h[~detected_low]=0
        self.height_bias_low = delta_h
        self.detected_height_low = detected_low
        delta_h = self.lapse_bias.copy()
        delta_h[~detected_low]=0
        self.lapse_bias_low = delta_h
 

    def get_ctth_bias_high(self, caObj):
        from get_flag_info import get_calipso_high_clouds
        high_clouds = get_calipso_high_clouds(caObj)
        detected_high = np.logical_and(self.detected_height, high_clouds[self.use])
        delta_h = self.height_bias.copy()
        delta_h[~detected_high]=0
        self.height_bias_high = delta_h.copy()
        self.detected_height_high = detected_high
        delta_h = self.lapse_bias.copy()
        delta_h[~detected_high]=0
        self.lapse_bias_high = delta_h


    def get_ctth_bias_low_temperature(self, caObj):
        from get_flag_info import get_calipso_low_clouds
        low_clouds = get_calipso_low_clouds(caObj)
        detected_low = np.logical_and(self.detected_height, low_clouds[self.use])
        temperature_pps = caObj.avhrr.all_arrays['ctth_temperature'][self.use]
        try:
            temperature_c = 273.15 + caObj.calipso.all_arrays['midlayer_temperature'][self.use,0]
            
            delta_t = temperature_pps - temperature_c
            delta_t[~detected_low]=0
            delta_t[temperature_c<0]=0
            delta_t[temperature_pps<0]=0
        except:
            delta_t = 0*temperature_pps
        try:
            temperature_pps = caObj.avhrr.all_arrays['bt11micron'][self.use]
            delta_t_t11 = temperature_pps - temperature_c
            delta_t_t11[~detected_low]=0
            delta_t_t11[temperature_c<0]=0
            delta_t_t11[temperature_pps<0]=0
        except:
            delta_t_t11 = 0*delta_t
        self.temperature_bias_low =delta_t
        #self.temperature_bias_low[~detected_low]=0
        self.temperature_bias_low_t11 = delta_t_t11

    def get_ctth_bias_type(self, caObj, calipso_cloudtype=0):
        from get_flag_info import get_calipso_clouds_of_type_i
        wanted_clouds = get_calipso_clouds_of_type_i(caObj, calipso_cloudtype=calipso_cloudtype)
        detected_typei = np.logical_and(self.detected_height, wanted_clouds[self.use])
        delta_h = self.height_bias.copy()
        delta_h[~detected_typei]=0
        self.height_bias_type[calipso_cloudtype] = delta_h
        self.detected_height_type[calipso_cloudtype] = detected_typei




class ppsStatsOnFibLatticeObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'definition': None,
            'lats': None,
            'lons': None,
            'N_false_clouds': None,
            'Min_lapse_rate': None,
            'N_detected_clouds': None,
            'N_new_false_clouds': None,
            'N_new_detected_clouds': None,
            'N_undetected_clouds': None,
            'N_detected_clear': None,
            'N_clear':None,
            'N_clouds':None,
            'N': None,
            'Kuipers': None}  
    def set_flattice(self, radius_km=200):
        self.radius_km = radius_km
        self.lons, self.lats =get_fibonacci_spread_points_on_earth(radius_km = radius_km)
        fig = plt.figure(figsize = (36,18))
        ax = fig.add_subplot(111)
        plt.plot(self.lons,self.lats,'b*')
        #plt.show()
        self.Sum_ctth_bias_low = np.zeros(self.lats.shape)
        self.Sum_lapse_bias_low = np.zeros(self.lats.shape)
        self.Sum_ctth_bias_high = np.zeros(self.lats.shape)
        self.Sum_ctth_mae_low = np.zeros(self.lats.shape)
        self.Sum_ctth_mae_high = np.zeros(self.lats.shape)
        self.Sum_ctth_mae = np.zeros(self.lats.shape)
        self.Sum_ctth_mae_diff = np.zeros(self.lats.shape)
        self.Sum_lapse_bias_high = np.zeros(self.lats.shape)
        self.Sum_ctth_bias_temperature_low = np.zeros(self.lats.shape)
        self.Sum_ctth_bias_temperature_low_t11 = np.zeros(self.lats.shape)
        self.Min_lapse_rate = np.zeros(self.lats.shape)
        self.N_new_false_clouds = np.zeros(self.lats.shape)
        self.N_new_detected_clouds = np.zeros(self.lats.shape)
        self.N_false_clouds = np.zeros(self.lats.shape)
        self.N_detected_clouds = np.zeros(self.lats.shape)
        self.N_undetected_clouds = np.zeros(self.lats.shape)
        self.N_detected_clear = np.zeros(self.lats.shape)
        self.N_detected_height_low = np.zeros(self.lats.shape)
        self.N_detected_height_high = np.zeros(self.lats.shape)
        self.N_detected_height = np.zeros(self.lats.shape)
        self.N_detected_height_both = np.zeros(self.lats.shape)
        self.Min_t11ts_offset = np.zeros(self.lats.shape)
        self.Max_t11t12_offset = np.zeros(self.lats.shape)
        self.Max_t37t12_offset = np.zeros(self.lats.shape)
        self.Max_t11t37_offset = np.zeros(self.lats.shape)
        self.Sum_height_bias_type={}
        self.N_detected_height_type={}
        for cc_type in xrange(8):
            self.Sum_height_bias_type[cc_type] = 1.0*np.zeros(self.lats.shape)
            self.N_detected_height_type[cc_type] = 1.0*np.zeros(self.lats.shape)
    def np_float_array(self):

        self.Sum_ctth_mae_low = 1.0*np.array(self.Sum_ctth_mae_low)
        self.Sum_ctth_mae_high = 1.0*np.array(self.Sum_ctth_mae_high)
        self.Sum_ctth_mae = 1.0*np.array(self.Sum_ctth_mae)
        self.Sum_ctth_mae_diff = 1.0*np.array(self.Sum_ctth_mae_diff)

        self.Sum_ctth_bias_low = 1.0*np.array(self.Sum_ctth_bias_low)
        self.Sum_lapse_bias_low = 1.0*np.array(self.Sum_lapse_bias_low)
        self.Sum_ctth_bias_temperature_low = 1.0*np.array(self.Sum_ctth_bias_temperature_low)
        self.Sum_ctth_bias_temperature_low_t11 = 1.0*np.array(self.Sum_ctth_bias_temperature_low_t11)
        self.Min_lapse_rate = 1.0*np.array(self.Min_lapse_rate)
        self.N_detected_clouds = 1.0*np.array(self.N_detected_clouds)
        self.N_undetected_clouds = 1.0*np.array(self.N_undetected_clouds)  
        self.N_false_clouds = 1.0*np.array(self.N_false_clouds)
        self.N_detected_clear = 1.0*np.array(self.N_detected_clear)
        self.N_new_detected_clouds = 1.0*np.array(self.N_new_detected_clouds)
        self.N_new_false_clouds = 1.0*np.array(self.N_new_false_clouds)
        self.N_detected_height_low = 1.0*np.array(self.N_detected_height_low)
    def find_number_of_clouds_clear(self):
        self.np_float_array()
        self.N_clear = self.N_detected_clear+self.N_false_clouds
        self.N_clouds = self.N_detected_clouds+self.N_undetected_clouds 
        self.N = self.N_clear + self.N_clouds
        self.Number_of = np.ma.masked_array(self.N, mask=self.N<0)
    
    def _remap_a_score_on_an_area(self, plot_area_name='npole', vmin=0.0, vmax=1.0, 
                                  score='Kuipers'):
        from pyresample import image, geometry
        area_def = utils.parse_area_file(
            'reshaped_files_plotting/region_config_test.cfg',  
            plot_area_name)[0]
        data = getattr(self, score)
        data = data.copy()
        if np.ma.is_masked(data):
            data[np.logical_and(np.equal(data.mask,False),data>vmax)]=vmax
            data[np.logical_and(np.equal(data.mask,False),data<vmin)]=vmin #do not wan't low ex hitrates set to nodata!
        else:
            data[data>vmax]=vmax
            data[data<vmin]=vmin   
        #lons = np.ma.masked_array(self.lons, mask=data.mask)
        #lats = np.ma.masked_array(self.lats, mask=data.mask)
        lons = self.lons
        lats = self.lats
        swath_def = geometry.SwathDefinition(lons=lons, lats=lats)
        swath_con = image.ImageContainerNearest(
            data, swath_def, 
            radius_of_influence=self.radius_km*1000*2.5,
            epsilon=1.0)
        area_con = swath_con.resample(area_def)
        result = area_con.image_data
        #pr.plot.show_quicklook(area_def, result,
        #                      vmin=vmin, vmax=vmax, label=score)
        if "FAR" in score:
            matplotlib.rcParams['image.cmap']= "BrBG"
        elif "diff" in score:  
            matplotlib.rcParams['image.cmap']= "BrBG"
        elif "mae" in score:
            matplotlib.rcParams['image.cmap']= "Reds"
        else:
            matplotlib.rcParams['image.cmap']= "BrBG"
        plot_label =score.replace('_','-')    
        if "mae" in score:
            plot_label =""
        pr.plot.save_quicklook(self.PLOT_DIR_SCORE + self.PLOT_FILENAME_START+
                               plot_area_name +'.png',
                               area_def, result, 
                               vmin=vmin, vmax=vmax, label=plot_label)

    def _remap_a_score_on_an_robinson_projection(self, vmin=0.0, vmax=1.0, 
                                                 score='Kuipers', screen_out_valid=False):
        from mpl_toolkits.basemap import Basemap
        from scipy.interpolate import griddata
        lons = self.lons
        lats = self.lats
        plt.close('all')
        ma_data = getattr(self, score)
        the_mask = ma_data.mask
        data=ma_data.data
        #data[np.logical_and(data>vmax,~the_mask)] = vmax
        #data[np.logical_and(data<vmin,~the_mask)] = vmin
        #reshape data a bit
        ind = np.argsort(lats)
        lons = lons[ind]
        lats = lats[ind]
        data = data[ind]
        the_mask = the_mask[ind]
        ind = np.argsort(lons)
        lons = lons[ind]
        lats = lats[ind]
        data =data[ind]
        the_mask = the_mask[ind]
        lons =         lons.reshape(len(data),1)#*3.14/180
        lats =         lats.reshape(len(data),1)#*3.14/180
        data =         data.reshape(len(data),1)
        the_mask =     the_mask.reshape(len(data),1)

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
        if "FAR" in score:
            my_cmap=copy.copy(matplotlib.cm.BrBG)
        elif "diff" in score:  
            my_cmap=copy.copy(matplotlib.cm.BrBG)
        elif "mae" in score:
            my_cmap=copy.copy(matplotlib.cm.Reds)
        else:    
            my_cmap=copy.copy(matplotlib.cm.BrBG)
        if score in "Bias" and screen_out_valid:
            #This screens out values between -5 and +5% 
            vmax=25
            vmin=-25            
            my_cmap=copy.copy(matplotlib.cm.get_cmap("BrBG", lut=100))
            cmap_vals = my_cmap(np.arange(100)) #extractvalues as an array
            cmap_vals[39:61] = [0.9, 0.9, 0.9, 1] #change the first value
            my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "newBrBG", cmap_vals) 
            print my_cmap
        if score in "RMS" and screen_out_valid:
            # This screens out values beteen 0 and 20%. 41/100=20%
            vmax=50
            vmin=0
            my_cmap=copy.copy(matplotlib.cm.get_cmap("BrBG", lut=100))
            cmap_vals = my_cmap(np.arange(100)) #extract values as an array
            cmap_vals[0:41] = [0.9, 0.9, 0.9, 1] #change the first value
            my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "newBrBG", cmap_vals) 
            print my_cmap

      

        #to mask out where we lack data
        data[np.logical_and(data>vmax,~the_mask)]=vmax
        data[np.logical_and(data<vmin,~the_mask)]=vmin
        data[the_mask]=2*vmax #give no data value that will be masked white
        xi = np.linspace(lon_min, lon_max, numcols)
        yi = np.linspace(lat_min, lat_max, numrows)
        xi, yi = np.meshgrid(xi, yi)
        # interpolate
        x, y, z = (np.array(lons.ravel()), 
                   np.array(lats.ravel()), 
                   np.array(data.ravel()))
        my_cmap.set_over(color='0.5',alpha=1)
        zi = griddata((x, y), z, (xi, yi), method='nearest')
        im1 = my_proj1.pcolormesh(xi, yi, zi, cmap=my_cmap,
                           vmin=vmin, vmax=vmax, latlon=True,  rasterized=True)
        im1.set_clim([vmin,vmax]) #to get nice ticks in the colorbar
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
        if not "mae" in score:
            ax.set_title(score.replace('_','-'), usetex=True)
        if score in ["ctth_mae"]:
            text_i = "(b)"
            if "v2014" in self.PLOT_FILENAME_START:
                text_i = "(a)"
            if "v2018" in self.PLOT_FILENAME_START:
                text_i = "(c)"    
            plt.text(0.01, 0.95, text_i, fontsize=36,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
                
        plt.savefig(self.PLOT_DIR_SCORE + self.PLOT_FILENAME_START+
                    '_robinson_' +'.pdf', bbox_inches='tight')
        #plt.savefig(self.PLOT_DIR_SCORE + self.PLOT_FILENAME_START+
        #            '_robinson_' +'.pdf', bbox_inches='tight')
        plt.close('all')

    def remap_and_plot_score_on_several_areas(self, vmin=0.0, vmax=1.0, 
                                              score='Kuipers', screen_out_valid=False):
        print score
        self.PLOT_DIR_SCORE = self.PLOT_DIR + "/%s/%s/"%(score, self.satellites)
        self.PLOT_FILENAME_START = "fig_%s_ccm_%s_%sfilter_dnt_%s_%s_r%skm_"%(
            self.satellites,self.cc_method,self.filter_method,
            self.DNT, score, self.radius_km)
         
        if not os.path.exists(self.PLOT_DIR_SCORE):
            os.makedirs(self.PLOT_DIR_SCORE)
        for plot_area_name in [
                #'cea5km_test'
                #'euro_arctic',
                #'ease_world_test'
                #'euro_arctic', #good
                #'antarctica', #good
                #'npole', #good
                'ease_nh_test',
                'ease_sh_test' ]:
            self._remap_a_score_on_an_area(plot_area_name=plot_area_name, 
                                           vmin=vmin, vmax=vmax, score=score)
        #the real robinson projection
        if "metop" not in self.satellites:
            self._remap_a_score_on_an_robinson_projection(vmin=vmin, vmax=vmax, 
                                                          score=score, screen_out_valid=False)
    def calculate_kuipers(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        N_clear = self.N_clear
        N_clouds = self.N_clouds
        N_detected_clouds = self.N_detected_clouds
        N_detected_clear = self.N_detected_clear
        #Typically we have N_clear/N_clouds = 30/70 
        #In areas with only clouds or only clears the Kuipers will be ==0
        #Even if all clouds/clears are classified correctly!
        #Do something for these set to none or update   
        Kuipers_devider = (N_clouds)*(N_clear)
        Kuipers_devider[Kuipers_devider==0] = 1.0
        Kuipers = (N_detected_clouds*N_detected_clear - 
                   self.N_false_clouds*self.N_undetected_clouds)/Kuipers_devider
        the_mask = np.logical_or(self.N_clear<20, self.N_clouds<20)
        the_mask = np.logical_or(the_mask, self.N_clouds < 0.01*self.N_clear)
        the_mask = np.logical_or(the_mask, self.N_clear < 0.01*self.N_clouds)
        use = np.logical_or(self.lats>=70,~the_mask)
        self.Kuipers_total_mean_polar = (
            (
                np.sum(N_detected_clouds[use])*np.sum(N_detected_clear[use]) - 
                np.sum(self.N_false_clouds[use])*np.sum(self.N_undetected_clouds[use]))/((np.sum(N_clouds[use]))*(np.sum(N_clear[use]))))
        Kuipers = np.ma.masked_array(Kuipers, mask=the_mask)
        self.Kuipers = Kuipers

    def calculate_lapse_rate(self):
        self.np_float_array()
        the_mask = self.Min_lapse_rate>-0.001
        lapse_rate = np.ma.masked_array(self.Min_lapse_rate, mask=the_mask)
        self.lapse_rate = lapse_rate
    def calculate_t11ts_offset(self):
        self.np_float_array()
        the_mask =  self.N_clear<10
        t11ts_offset = np.ma.masked_array(self.Min_t11ts_offset, mask=the_mask)
        self.t11ts_offset = t11ts_offset
    def calculate_t11t12_offset(self):
        self.np_float_array()
        the_mask =  self.N_clear<10
        t11t12_offset = np.ma.masked_array(self.Max_t11t12_offset, mask=the_mask)
        self.t11t12_offset = t11t12_offset
    def calculate_t37t12_offset(self):
        self.np_float_array()
        the_mask =  self.N_clear<10
        t37t12_offset = np.ma.masked_array(self.Max_t37t12_offset, mask=the_mask)
        self.t37t12_offset = t37t12_offset
    def calculate_t11t37_offset(self):
        self.np_float_array()
        the_mask =  self.N_clear<10
        t11t37_offset = np.ma.masked_array(self.Max_t11t37_offset, mask=the_mask)
        self.t11t37_offset = t11t37_offset
    def calculate_temperature_bias(self):
        self.np_float_array()
        the_mask = self.N_detected_height_low<10
        ctth_bias_temperature_low = self.Sum_ctth_bias_temperature_low*1.0/self.N_detected_height_low
        ctth_bias_temperature_low = np.ma.masked_array(ctth_bias_temperature_low, mask=the_mask)
        self.ctth_bias_temperature_low = ctth_bias_temperature_low
    def calculate_temperature_bias_t11(self):
        self.np_float_array()
        the_mask = self.N_detected_height_low<10
        ctth_bias_temperature_low_t11 = self.Sum_ctth_bias_temperature_low_t11*1.0/self.N_detected_height_low
        ctth_bias_temperature_low_t11 = np.ma.masked_array(ctth_bias_temperature_low_t11, mask=the_mask)
        self.ctth_bias_temperature_low_t11 = ctth_bias_temperature_low_t11
    def calculate_height_bias(self):
        self.np_float_array()
        the_mask = self.N_detected_height_low<10
        ctth_bias_low = self.Sum_ctth_bias_low*1.0/self.N_detected_height_low
        ctth_bias_low = np.ma.masked_array(ctth_bias_low, mask=the_mask)
        self.ctth_bias_low = ctth_bias_low
        the_mask = self.N_detected_height_high<10
        ctth_bias_high = self.Sum_ctth_bias_high*1.0/self.N_detected_height_high
        ctth_bias_high = np.ma.masked_array(ctth_bias_high, mask=the_mask)
        self.ctth_bias_high = ctth_bias_high
    def calculate_height_mae(self):
        self.np_float_array()
        the_mask = self.N_detected_height_low<10
        ctth_mae_low = self.Sum_ctth_mae_low*1.0/self.N_detected_height_low
        ctth_mae_low = np.ma.masked_array(ctth_mae_low, mask=the_mask)
        self.ctth_mae_low = ctth_mae_low
        the_mask = self.N_detected_height_high<10
        ctth_mae_high = self.Sum_ctth_mae_high*1.0/self.N_detected_height_high
        ctth_mae_high = np.ma.masked_array(ctth_mae_high, mask=the_mask)
        self.ctth_mae_high = ctth_mae_high
        the_mask = self.N_detected_height<10
        ctth_mae = self.Sum_ctth_mae*1.0/self.N_detected_height
        ctth_mae = np.ma.masked_array(ctth_mae, mask=the_mask)
        self.ctth_mae = ctth_mae
        the_mask = self.N_detected_height_both<10
        ctth_mae_diff = self.Sum_ctth_mae_diff*1.0/self.N_detected_height_both
        ctth_mae_diff = np.ma.masked_array(ctth_mae_diff, mask=the_mask)
        self.ctth_mae_diff = ctth_mae_diff
        the_mask = self.N_detected_height_both<0
        self.N_detected_height_both = np.ma.masked_array(self.N_detected_height_both, mask=the_mask)

    def calculate_height_bias_lapse(self):
        self.np_float_array()
        the_mask = self.N_detected_height_low<10
        lapse_bias_low = self.Sum_lapse_bias_low*1.0/self.N_detected_height_low
        lapse_bias_low = np.ma.masked_array(lapse_bias_low, mask=the_mask)
        self.lapse_bias_low = lapse_bias_low
        the_mask = self.N_detected_height_high<10
        lapse_bias_high = self.Sum_lapse_bias_high*1.0/self.N_detected_height_high
        lapse_bias_high = np.ma.masked_array(lapse_bias_high, mask=the_mask)
        self.lapse_bias_high = lapse_bias_high
    def calculate_height_bias_type(self):
        self.np_float_array()
        for cc_type in xrange(8):
            the_mask = self.N_detected_height_type[cc_type]<10
            ctth_bias_type_i = self.Sum_height_bias_type[cc_type]*1.0/self.N_detected_height_type[cc_type]
            ctth_bias_type_i = np.ma.masked_array(ctth_bias_type_i, mask=the_mask)
            setattr(self, "ctth_bias_type_%d"%(cc_type), ctth_bias_type_i)   

    def calculate_hitrate(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        Hitrate = (self.N_detected_clouds + self.N_detected_clear)*1.0/(
            self.N_clear + self.N_clouds)
        the_mask = self.N<20
        Hitrate = np.ma.masked_array(Hitrate, mask=the_mask)
        use = np.logical_or(self.lats>=70,~the_mask)
        self.Hitrate_total_mean_polar = (
            np.sum(self.N_detected_clouds[use]) + 
            np.sum(self.N_detected_clear[use]))*1.0/(
                np.sum(self.N[use]))
        self.Hitrate = Hitrate

    def calculate_increased_hitrate(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        Hitrate = (
            self.N_detected_clouds + self.N_detected_clear)*1.0/(
                self.N_clear + self.N_clouds)
        new_Hitrate = (
            self.N_detected_clouds + self.N_new_detected_clouds - 
            self.N_new_false_clouds+ self.N_detected_clear)*1.0/(
                self.N_clear + self.N_clouds)
        the_mask = self.N<20
        increased_Hitrate = np.ma.masked_array(new_Hitrate-Hitrate, mask=the_mask)
        self.increased_Hitrate = increased_Hitrate

    def calculate_threat_score(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        ThreatScore = (
            self.N_detected_clouds)*1.0/( self.N_clouds + self.N_false_clouds)
        the_mask = self.N_clouds<20
        ThreatScore = np.ma.masked_array(ThreatScore, mask=the_mask)
        self.Threat_Score = ThreatScore
    def calculate_threat_score_clear(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        ThreatScoreClear = (
            self.N_detected_clear)*1.0/( self.N_clear + 
                                         self.N_undetected_clouds)
        the_mask = self.N_clear<20
        ThreatScoreClear = np.ma.masked_array(ThreatScoreClear, 
                                              mask=the_mask)
        self.Threat_Score_Clear = ThreatScoreClear
    def calculate_pod_clear(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        PODclear = (
            self.N_detected_clear)*1.0/(self.N_clear)
        the_mask = self.N_clear<20
        PODclear = np.ma.masked_array(PODclear, mask=the_mask)
        self.PODclear = PODclear 
        use = np.logical_or(self.lats>=70,~the_mask)
        self.PODclear_total_mean_polar = (
            np.sum(self.N_detected_clear[use]))*1.0/(
                np.sum(self.N_clear[use]))

    def calculate_pod_cloudy(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        PODcloudy = (
            self.N_detected_clouds)*1.0/(self.N_clouds)
        the_mask = self.N_clouds<20
        PODcloudy = np.ma.masked_array(PODcloudy, mask=the_mask)
        self.PODcloudy = PODcloudy
        use = np.logical_or(self.lats>=70,~the_mask)
        self.PODcloudy_total_mean_polar = (
            np.sum(self.N_detected_clouds[use]))*1.0/(
                np.sum(self.N_clouds[use]))

    def calculate_far_clear(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        FARclear = (
            self.N_undetected_clouds)*1.0/(self.N_detected_clear +
                                           self.N_undetected_clouds)
        the_mask = self.N_clear<20
        FARclear = np.ma.masked_array(FARclear, mask=the_mask)
        self.FARclear = FARclear 
    def calculate_far_cloudy(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        FARcloudy = (
            self.N_false_clouds)*1.0/(self.N_detected_clouds+
                                      self.N_false_clouds)     
        the_mask = self.N_clouds<20
        FARcloudy = np.ma.masked_array(FARcloudy, mask=the_mask)
        self.FARcloudy = FARcloudy
        use = np.logical_or(self.lats>=70,~the_mask)
        self.FARcloudy_total_mean_polar = (
            np.sum(self.N_false_clouds[use]))*1.0/(
                np.sum(self.N_detected_clouds[use])+ 
                np.sum(self.N_false_clouds[use]))
    def calculate_calipso_cfc(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        calipso_cfc = 100*(
            self.N_detected_clouds + self.N_undetected_clouds)*1.0/(self.N)
        the_mask = self.N<20
        calipso_cfc = np.ma.masked_array(calipso_cfc, mask=the_mask)
        self.calipso_cfc = calipso_cfc
    def calculate_pps_cfc(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        pps_cfc = 100*(
            self.N_detected_clouds + self.N_false_clouds)*1.0/(self.N)
        the_mask = self.N<20
        pps_cfc = np.ma.masked_array(pps_cfc, mask=the_mask)
        self.pps_cfc = pps_cfc
    def calculate_bias(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        self.calculate_calipso_cfc()
        self.calculate_pps_cfc()
        Bias = self.pps_cfc - self.calipso_cfc
        the_mask = self.N<20 
        Bias = np.ma.masked_array(Bias, mask=the_mask)
        use = np.logical_or(self.lats>=70,~the_mask)
        self.Bias_total_mean_polar = 100*(
            -np.sum(self.N_undetected_clouds[use]) + 
            np.sum(self.N_false_clouds[use]))*1.0/(np.sum(self.N[use]))
        print "hej",self.Bias_total_mean_polar
        self.Bias = Bias
    def calculate_rms(self):
        self.np_float_array()
        self.calculate_calipso_cfc()
        self.find_number_of_clouds_clear()
        self.calculate_bias()
        RMS = np.sqrt((self.N_false_clouds*(100.0 - 0.0 - self.Bias)**2 + 
                       self.N_undetected_clouds*(0.0 - 100.0 -self.Bias)**2 + 
                       self.N_detected_clear*self.Bias**2 + 
                       self.N_detected_clouds*self.Bias**2)/(
                   self.N))
        the_mask = self.N<20 
        RMS = np.ma.masked_array(RMS, mask=the_mask)
        self.RMS = RMS


class PerformancePlottingObject:
    def __init__(self):
        self.flattice = ppsStatsOnFibLatticeObject()
    def add_detection_stats_on_fib_lattice(self, my_obj):
        #Start with the area and get lat and lon to calculate the stats:
        if  len(my_obj.longitude) == 0:
            print "Skipping file, no matches !"
            return
        lats = self.flattice.lats[:]
        max_distance=self.flattice.radius_km*1000*2.5
        area_def = SwathDefinition(*(self.flattice.lons,
                                     self.flattice.lats))
        target_def = SwathDefinition(*(my_obj.longitude, 
                                       my_obj.latitude)) 
        valid_in, valid_out, indices, distances = get_neighbour_info(
            area_def, target_def, radius_of_influence=max_distance, 
            epsilon=100, neighbours=1)
        cols = get_sample_from_neighbour_info('nn', target_def.shape,
                                              np.array(xrange(0,len(lats))),
                                              valid_in, valid_out,
                                              indices)
        cols = cols[valid_out]
        detected_clouds = my_obj.detected_clouds[valid_out]
        detected_clear = my_obj.detected_clear[valid_out]
        detected_height_low = my_obj.detected_height_low[valid_out]
        detected_height_high = my_obj.detected_height_high[valid_out]
        detected_height = my_obj.detected_height[valid_out]
        detected_height_both = my_obj.detected_height_both[valid_out]
        false_clouds = my_obj.false_clouds[valid_out]
        undetected_clouds = my_obj.undetected_clouds[valid_out]
        new_detected_clouds = my_obj.new_detected_clouds[valid_out]
        new_false_clouds = my_obj.new_false_clouds[valid_out]
        lapse_rate = my_obj.lapse_rate[valid_out]  
        t11ts_offset = my_obj.t11ts_offset[valid_out] 
        t11t12_offset = my_obj.t11t12_offset[valid_out] 
        t37t12_offset = my_obj.t37t12_offset[valid_out] 
        t11t37_offset = my_obj.t11t37_offset[valid_out] 
        height_bias_low = my_obj.height_bias_low[valid_out]
        height_bias = my_obj.height_bias[valid_out]
        height_mae_diff = my_obj.height_mae_diff[valid_out]
        temperature_bias_low = my_obj.temperature_bias_low[valid_out]
        temperature_bias_low_t11 = my_obj.temperature_bias_low_t11[valid_out]
        lapse_bias_low = my_obj.lapse_bias_low[valid_out]
        height_bias_high = my_obj.height_bias_high[valid_out]
        lapse_bias_high = my_obj.lapse_bias_high[valid_out]
        is_clear = np.logical_or(detected_clear,false_clouds)
        #lets make things faster, I'm tired of waiting!
        cols[distances>max_distance]=-9 #don't use pixles matched too far away!
        import time        
        tic = time.time()      
        arr, counts = np.unique(cols, return_index=False, return_counts=True)        
        for d in arr[arr>0]:
            use = cols==d
            ind = np.where(use)[0]
            #if ind.any():
            self.flattice.N_false_clouds[d] += np.sum(false_clouds[ind])
            self.flattice.N_detected_clouds[d] += np.sum(detected_clouds[ind])
            self.flattice.N_detected_clear[d] += np.sum(detected_clear[ind])
            self.flattice.N_undetected_clouds[d] += np.sum(undetected_clouds[ind])
            self.flattice.N_new_false_clouds[d] += np.sum(new_false_clouds[ind])
            self.flattice.N_new_detected_clouds[d] += np.sum(new_detected_clouds[ind])
            self.flattice.N_detected_height_low[d] += np.sum(detected_height_low[ind])
            self.flattice.N_detected_height_high[d] += np.sum(detected_height_high[ind])
            self.flattice.N_detected_height[d] += np.sum(detected_height[ind])
            self.flattice.N_detected_height_both[d] += np.sum(detected_height_both[ind])
            self.flattice.Sum_ctth_bias_low[d] += np.sum(height_bias_low[ind])
            self.flattice.Sum_ctth_mae_low[d] += np.sum(np.abs(height_bias_low[ind]))
            self.flattice.Sum_ctth_mae[d] += np.sum(np.abs(height_bias[ind]))
            self.flattice.Sum_ctth_mae_diff[d] += np.sum(height_mae_diff[ind])
            self.flattice.Sum_lapse_bias_low[d] += np.sum(lapse_bias_low[ind])
            self.flattice.Sum_ctth_bias_high[d] += np.sum(height_bias_high[ind])
            self.flattice.Sum_ctth_mae_high[d] += np.sum(np.abs(height_bias_high[ind]))
            self.flattice.Sum_lapse_bias_high[d] += np.sum(lapse_bias_high[ind])
            self.flattice.Sum_ctth_bias_temperature_low[d] += np.sum(temperature_bias_low[ind])
            self.flattice.Sum_ctth_bias_temperature_low_t11[d] += np.sum(temperature_bias_low_t11[ind])
            self.flattice.Min_lapse_rate[d] = np.min([self.flattice.Min_lapse_rate[d],
                                                      np.min(lapse_rate[ind])])  
            if np.sum(is_clear[ind])>0:
                self.flattice.Min_t11ts_offset[d] = np.min([self.flattice.Min_t11ts_offset[d],
                                                            np.percentile(t11ts_offset[ind][is_clear[ind]], 5)])
                self.flattice.Max_t11t12_offset[d] = np.max([self.flattice.Max_t11t12_offset[d],
                                                             np.percentile(t11t12_offset[ind][is_clear[ind]], 95)])
                self.flattice.Max_t37t12_offset[d] = np.max([self.flattice.Max_t37t12_offset[d],
                                                             np.percentile(t37t12_offset[ind][is_clear[ind]], 95)])
                self.flattice.Max_t11t37_offset[d] = np.max([self.flattice.Max_t11t37_offset[d],
                                                             np.percentile(t11t37_offset[ind][is_clear[ind]], 95)])
            cc_type = 0
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 1
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 2
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 3
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 4
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 5
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 6
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])
            cc_type = 7
            self.flattice.Sum_height_bias_type[cc_type][d] += np.sum(my_obj.height_bias_type[cc_type][ind])
            self.flattice.N_detected_height_type[cc_type][d] += np.sum(my_obj.detected_height_type[cc_type][ind])

        print "mapping took %1.4f seconds"%(time.time()-tic)   


def get_fibonacci_spread_points_on_earth(radius_km):
    #Earth area = 510072000km2
    #4000 point with radius~200km
    #1000 point with radium~100km
    #25000 radius 80km
    #64000 radius 5km
    EARTH_AREA = 510072000
    POINT_AREA = radius_km * radius_km * 3.14
    n = int(EARTH_AREA /POINT_AREA)
    #http://arxiv.org/pdf/0912.4540.pdf
    #Alvaro Gonzalez: Measurement of areas on sphere usig Fibonacci and latitude-longitude grid.
    #import math
    lin_space = np.array(xrange(-n/2,n/2))
    pi = 3.14
    theta = (1+np.sqrt(5))*0.5
    longitude = (lin_space % theta)*360/theta
    temp = (2.0 * lin_space) / (n)
    temp[temp>1.0]=0.999
    temp[temp<-1.0]=-0.999
    latitude = np.arcsin(temp)*180/pi
    longitude[longitude>180] = longitude[longitude>180] -360
    longitude[longitude<-180] = longitude[longitude<-180] +360
    #latitude[latitude>90]=180 - latitude[latitude>90]
    #latitude[latitude<-90]=-180 -latitude[latitude<-90]
    longitude =longitude[latitude<90]
    latitude =latitude[latitude<90]
    longitude =longitude[latitude>-90]
    latitude =latitude[latitude>-90]

    if np.isnan(np.max(latitude)):
        raise ValueError
    return longitude, latitude


    




