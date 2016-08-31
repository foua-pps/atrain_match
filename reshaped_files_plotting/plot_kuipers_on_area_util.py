import numpy as np
from pyresample import utils
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import get_neighbour_info
from pyresample.kd_tree import get_sample_from_neighbour_info
import pyresample as pr
from scipy import ndimage
import matplotlib
#matplotlib.use("TkAgg")
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
class ppsRemappedObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'definition': None,
            'lats': None,
            'lons': None,
            'N_false_clouds': None,
            'N_detected_clouds': None,
            'N_new_false_clouds': None,
            'N_new_detected_clouds': None,
            'N_undetected_clouds': None,
            'N_detected_clear': None,
            'Kuipers': None}  
    def set_area(self, radius_km=200):
        self.radius_km = radius_km
        self.lons, self.lats =get_fibonacci_spread_points_on_earth(radius_km = radius_km)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize = (36,18))
        ax = fig.add_subplot(111)
        plt.plot(self.lons,self.lats,'b*')
        #plt.show()
        self.N_new_false_clouds = np.zeros(self.lats.shape)
        self.N_new_detected_clouds = np.zeros(self.lats.shape)
        self.N_false_clouds = np.zeros(self.lats.shape)
        self.N_detected_clouds = np.zeros(self.lats.shape)
        self.N_undetected_clouds = np.zeros(self.lats.shape)
        self.N_detected_clear = np.zeros(self.lats.shape)
    def np_float_array(self):    
        self.N_detected_clouds = 1.0*np.array(self.N_detected_clouds)
        self.N_undetected_clouds = 1.0*np.array(self.N_undetected_clouds)  
        self.N_false_clouds = 1.0*np.array(self.N_false_clouds)
        self.N_detected_clear = 1.0*np.array(self.N_detected_clear)
        self.N_new_detected_clouds = 1.0*np.array(self.N_new_detected_clouds)
        self.N_new_false_clouds = 1.0*np.array(self.N_new_false_clouds)
    def find_number_of_clouds_clear(self):
        self.np_float_array()
        self.N_clear = self.N_detected_clear+self.N_false_clouds
        self.N_clouds = self.N_detected_clouds+self.N_undetected_clouds 
        self.N = self.N_clear + self.N_clouds

    def _remap_score(self, vmin=0.0, vmax=1.0, 
                     score='Kuipers', screen_out_valid=False):
        print score
        from pyresample import image, geometry
        for plot_area_name in [
                #'cea5km_test'
                #'euro_arctic',
                #'ease_world_test'
                'euro_arctic',
                'antarctica',
                'npole',
                'ease_nh_test',
                'ease_sh_test' ]:
            area_def = utils.parse_area_file(
                'reshaped_files_plotting/region_config_test.cfg',  
                plot_area_name)[0]
            data = getattr(self, score)
            data = data.copy()
            data[np.logical_and(np.equal(data.mask,False),data>vmax)]=vmax
            data[np.logical_and(np.equal(data.mask,False),data<vmin)]=vmin #do not wan't low ex hitrates set to nodata!
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
        
            pr.plot.save_quicklook(self.PLOT_DIR + self.figure_name + 
                                   score +'_' + plot_area_name +'.png',
                                   area_def, result, 
                                   vmin=vmin, vmax=vmax, label=score)
        #the real robinson projection
        if "morning" not in self.figure_name:
            import matplotlib.pyplot as plt
            from mpl_toolkits.basemap import Basemap
            from scipy.interpolate import griddata
            plt.close('all')
            data = getattr(self, score)
            the_mask = data.mask
            data=np.ma.masked_invalid(data)
            data[np.logical_and(data>vmax,~the_mask)] = vmax
            data[np.logical_and(data<vmin,~the_mask)] = vmin
            #lons = lons[np.not_equal(the_mask, True)]
            #lats = lats[np.not_equal(the_mask, True)]
            #data = data[np.not_equal(the_mask, True)]
            ind = np.argsort(lats)
            lons = lons[ind]
            lats = lats[ind]
            data =data[ind]
            ind = np.argsort(lons)
            lons = lons[ind]
            lats = lats[ind]
            data =data[ind]
            lons = lons.reshape(len(data),1)#*3.14/180
            lats = lats.reshape(len(data),1)#*3.14/180
            data =data.reshape(len(data),1)
            proj1 = Basemap(projection='robin',lon_0=0,resolution='c')
            m = proj1
            numcols=1000
            numrows=500
            lat_min = -83.0
            lon_min = -179.9
            lat_max = 83.0
            lon_max = 179.9
            
            fig = plt.figure(figsize = (16,9))
            ax = fig.add_subplot(111)
            # transform lon / lat coordinates to map projection
            plons, plats = m(*(lons, lats))        
            # grid the data to lat/lo grid, approximation but needed fo plotting
            # numcols, numrows = 1000, 500
            xi = np.linspace(lon_min, lon_max, numcols)
            yi = np.linspace(lat_min, lat_max, numrows)
            xi, yi = np.meshgrid(xi, yi)
            # interpolate
            x, y, z = (np.array(lons.ravel()), 
                       np.array(lats.ravel()), 
                       np.array(data.ravel()))
            import copy; 
            my_cmap=copy.copy(matplotlib.cm.coolwarm)
            my_cmap_for_masked=copy.copy(matplotlib.cm.coolwarm)
            if score in "Bias" and screen_out_valid:
                #This screens out values between -5 and +5% 
                vmax=25
                vmin=-25
                my_cmap=copy.copy(matplotlib.cm.get_cmap("coolwarm", lut=100))
                cmap_vals = my_cmap(np.arange(100)) #extractvalues as an array
                cmap_vals[39:61] = [0.9, 0.9, 0.9, 1] #change the first value
                my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                    "newwarmcool", cmap_vals) 
                print my_cmap
            if score in "RMS" and screen_out_valid:
                # This screens out values beteen 0 and 20%. 41/100=20%
                vmax=50
                vmin=0
                my_cmap=copy.copy(matplotlib.cm.get_cmap("coolwarm", lut=100))
                cmap_vals = my_cmap(np.arange(100)) #extract values as an array
                cmap_vals[0:41] = [0.9, 0.9, 0.9, 1] #change the first value
                my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                    "newwarmcool", cmap_vals) 
                print my_cmap
            my_cmap_for_masked.set_over('w') #masked should be white
            my_cmap_for_masked.set_under('k', alpha=0) #do not cover the data, set to transparent


            zi = griddata((x, y), z, (xi, yi), method='nearest')
            remapped_mask = griddata((x, y),  
                                     np.where(data.mask.ravel(),1.0,0.0), 
                                     (xi, yi), method='nearest')

            im1 = m.pcolormesh(xi, yi, zi, cmap=my_cmap,
                               vmin=vmin, vmax=vmax, latlon=True)
            im2 = m.pcolormesh(xi, yi, remapped_mask, cmap=my_cmap_for_masked,
                               vmin=0.4, vmax=0.5, latlon=True)
            #draw som lon/lat lines
            m.drawparallels(np.arange(-90.,90.,30.))
            m.drawmeridians(np.arange(-180.,180.,60.))
            m.drawcoastlines()
            m.drawmapboundary(fill_color='0.9')
            cb = m.colorbar(im1,"right", size="5%", pad="2%")
            ax.set_title(score)
            plt.savefig(self.PLOT_DIR + self.figure_name + 
                        'basemap_' + 
                        score +'_robinson_' +'.png')
            plt.close('all')
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
        Kuipers = np.ma.masked_array(Kuipers, mask=the_mask)
        self.Kuipers = Kuipers

    def calculate_hitrate(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        Hitrate = (
            self.N_detected_clouds + self.N_detected_clear)*1.0/(
                self.N_clear + self.N_clouds)
        the_mask = self.N<20
        Hitrate = np.ma.masked_array(Hitrate, mask=the_mask)
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
    def calculate_pod_cloudy(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        PODcloudy = (
            self.N_detected_clouds)*1.0/(self.N_clouds)
        the_mask = self.N_clouds<20
        PODcloudy = np.ma.masked_array(PODcloudy, mask=the_mask)
        self.PODcloudy = PODcloudy
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
        self.Bias = Bias
    def calculate_RMS(self):
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
        self.area = ppsRemappedObject()
    def add_detection_stats_on_area_map_grid(self, my_obj):
        #Start with the area and get lat and lon to calculate the stats:
        #area_def = self.area.definition 
        lats = self.area.lats[:]
        max_distance=self.area.radius_km*1000*2.5
        #print "lons shape", lons.shape
        if lats.ndim ==1:
            area_def = SwathDefinition(*(self.area.lons,
                                         self.area.lats))
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
            false_clouds = my_obj.false_clouds[valid_out]
            undetected_clouds = my_obj.undetected_clouds[valid_out]
            new_detected_clouds = my_obj.new_detected_clouds[valid_out]
            new_false_clouds = my_obj.new_false_clouds[valid_out]
            for c, ind in zip(cols.ravel(), xrange(len(cols.ravel()))):
                if distances[ind]<max_distance:
                    self.area.N_false_clouds[c] += false_clouds[ind]
                    self.area.N_detected_clouds[c] += detected_clouds[ind]
                    self.area.N_detected_clear[c] += detected_clear[ind]
                    self.area.N_undetected_clouds[c] += undetected_clouds[ind]
                    self.area.N_new_false_clouds[c] += new_false_clouds[ind]
                    self.area.N_new_detected_clouds[c] += new_detected_clouds[ind]
        else:
            target_def = SwathDefinition(*(my_obj.longitude, 
                                           my_obj.latitude)) 
            valid_in, valid_out, indices, distances = get_neighbour_info(
                area_def, target_def, radius_of_influence=max_distance, neighbours=1)
            #Use pyresampe code to find colmun and row numbers for each pixel
            cols_matrix, rows_matrix = np.meshgrid(np.array(xrange(0,lats.shape[1])),
                                                   np.array(xrange(0,lats.shape[0])))
            cols = get_sample_from_neighbour_info('nn', target_def.shape,
                                                  cols_matrix,
                                                  valid_in, valid_out,
                                                  indices)
            rows = get_sample_from_neighbour_info('nn', target_def.shape,
                                                  rows_matrix,
                                                  valid_in, valid_out,
                                                  indices)
            rows = rows[valid_out]
            cols = cols[valid_out]
            detected_clouds = my_obj.detected_clouds[valid_out]
            detected_clear = my_obj.detected_clear[valid_out]
            false_clouds = my_obj.false_clouds[valid_out]
            undetected_clouds = my_obj.undetected_clouds[valid_out]
            for r, c, ind in zip(rows.ravel(), cols.ravel(), xrange(len(cols.ravel()))):
                if distances[ind]<max_distance:
                    self.area.N_false_clouds[r,c] += false_clouds[ind]
                    self.area.N_detected_clouds[r,c] += detected_clouds[ind]
                    self.area.N_detected_clear[r,c] += detected_clear[ind]
                    self.area.N_undetected_clouds[r,c] += undetected_clouds[ind]  


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
    
def get_some_info_from_caobj(caObj, isGAC=True, isACPGv2012=False, 
                             method='KG', DNT='All'):
    my_obj = ppsMatch_Imager_CalipsoObject()
    #cloudObj = get_clear_and_cloudy_vectors(caObj, isACPGv2012, isGAC)  
    isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                                 caObj.avhrr.all_arrays['cloudtype']<21) 
    isClearPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                caObj.avhrr.all_arrays['cloudtype']<5)
    nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
    if method == 'KG' and isGAC:
        isCalipsoCloudy = np.logical_and(
            caObj.calipso.all_arrays['cloud_fraction']>0.5,
            caObj.calipso.all_arrays['total_optical_depth_5km']>0.15)
        isCalipsoClear = np.not_equal(isCalipsoCloudy,True)
    elif method == 'Nina' and isGAC:    
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
    elif method == 'KG':
        isCalipsoCloudy = nlay>0
        isCalipsoClear = np.not_equal(isCalipsoCloudy, True)   
    elif method == 'Nina':
        isCalipsoCloudy = np.logical_or(
            caObj.calipso.all_arrays['total_optical_depth_5km']>0.15, 
            np.logical_and(caObj.calipso.all_arrays['total_optical_depth_5km']<0,
                           nlay>0))
        isCalipsoClear = np.logical_and(nlay == 0, meancl<0.01)
    elif method =='KG_r13_extratest':
        isCalipsoCloudy = nlay>0  
        isCalipsoClear = np.not_equal(isCalipsoCloudy,True)
        r13 = caObj.avhrr.all_arrays['r13micron']
        sunz =  caObj.avhrr.all_arrays['sunz']
        sunz_cos = sunz.copy()
        sunz_cos[sunz>87] =87    
        r13[sunz<90] = r13[sunz<90]/np.cos(np.radians(sunz_cos[sunz<90]))
        isCloud_r13 = np.logical_and(r13>2.0, caObj.avhrr.all_arrays['ciwv']>3)
    elif method =='r13_extratest':
        isCalipsoCloudy = np.logical_or(
            caObj.calipso.all_arrays['total_optical_depth_5km']>0.15, 
            np.logical_and(caObj.calipso.all_arrays['total_optical_depth_5km']<0,
                           nlay>0))
        isCalipsoClear = np.logical_and(nlay == 0, meancl<0.01)
        r13 = caObj.avhrr.all_arrays['r13micron']
        sunz =  caObj.avhrr.all_arrays['sunz']
        sunz_cos = sunz.copy()
        sunz_cos[sunz>87] =87    
        r13[sunz<90] = r13[sunz<90]/np.cos(np.radians(sunz_cos[sunz<90]))
        isCloud_r13 = np.logical_and(r13>2.0, caObj.avhrr.all_arrays['ciwv']>3)
        

    if DNT in ["day"]:
        isCloudyPPS = np.logical_and(isCloudyPPS, 
                                     caObj.avhrr.all_arrays['sunz']<=80)
        isClearPPS =  np.logical_and(isClearPPS, 
                                     caObj.avhrr.all_arrays['sunz']<=80)
    if DNT in ["night"]:
        isCloudyPPS = np.logical_and(isCloudyPPS, 
                                     caObj.avhrr.all_arrays['sunz']>=95)
        isClearPPS =  np.logical_and(isClearPPS, 
                                     caObj.avhrr.all_arrays['sunz']>=95)
    if DNT in ["twilight"]:
        isCloudyPPS = np.logical_and(isCloudyPPS, 
                                     caObj.avhrr.all_arrays['sunz']>80)
        isClearPPS =  np.logical_and(isClearPPS, 
                                     caObj.avhrr.all_arrays['sunz']>80)
        isCloudyPPS = np.logical_and(isCloudyPPS, 
                                     caObj.avhrr.all_arrays['sunz']<95)
        isClearPPS =  np.logical_and(isClearPPS, 
                                     caObj.avhrr.all_arrays['sunz']<95)                  
               
    undetected_clouds = np.logical_and(isCalipsoCloudy, isClearPPS)
    false_clouds = np.logical_and(isCalipsoClear, isCloudyPPS)
    detected_clouds = np.logical_and(isCalipsoCloudy, isCloudyPPS)
    detected_clear = np.logical_and(isCalipsoClear, isClearPPS)
    use = np.logical_or(np.logical_or(detected_clouds, detected_clear),
                        np.logical_or(false_clouds, undetected_clouds))
    my_obj.false_clouds = false_clouds[use]
    my_obj.detected_clouds = detected_clouds[use]
    my_obj.undetected_clouds = undetected_clouds[use]
    my_obj.detected_clear = detected_clear[use]
    my_obj.latitude = caObj.avhrr.latitude[use]
    my_obj.longitude = caObj.avhrr.longitude[use]  
    if "r13" in method:
        new_detected_clouds = np.logical_and(
            isCalipsoCloudy,
            np.logical_and(isClearPPS, isCloud_r13))
        new_false_clouds = np.logical_and(
            isCalipsoClear,
            np.logical_and(isClearPPS, isCloud_r13))
        my_obj.new_false_clouds = new_false_clouds[use]
        my_obj.new_detected_clouds = new_detected_clouds[use]
    else:
       my_obj.new_false_clouds = np.zeros(
           my_obj.false_clouds.shape)
       my_obj.new_detected_clouds = np.zeros(
           my_obj.false_clouds.shape)
    return my_obj
