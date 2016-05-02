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
            'N_undetected_clouds': None,
            'N_detected_clear': None,
            'Kuipers': None}  
    def set_area(self, area_name = 'robinson_test', radius_km=200):
        #area_name = 'robinson_test'
        area_def = utils.parse_area_file(
            'reshaped_files_plotting/region_config_test.cfg',  
            area_name)[0]
        self.radius_km = radius_km
        self.area_name=area_name
        self.definition = area_def
        #self.lons = area_def.lons[:]
        #self.lats = area_def.lats[:]
        self.lons, self.lats =get_fibonacci_spread_points_on_earth(radius_km = radius_km)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize = (36,18))
        ax = fig.add_subplot(111)
        plt.plot(self.lons,self.lats,'b*')
        plt.show()
        self.N_false_clouds = np.zeros(self.lats.shape)
        self.N_detected_clouds = np.zeros(self.lats.shape)
        self.N_undetected_clouds = np.zeros(self.lats.shape)
        self.N_detected_clear = np.zeros(self.lats.shape)
    def np_float_array(self):    
        self.N_detected_clouds = 1.0*np.array(self.N_detected_clouds)
        self.N_undetected_clouds = 1.0*np.array(self.N_undetected_clouds)  
        self.N_false_clouds = 1.0*np.array(self.N_false_clouds)
        self.N_detected_clear = 1.0*np.array(self.N_detected_clear)
    def find_number_of_clouds_clear(self):
        self.np_float_array()
        self.N_clear = self.N_detected_clear+self.N_false_clouds
        self.N_clouds = self.N_detected_clouds+self.N_undetected_clouds 
        self.N = self.N_clear + self.N_clouds

    def _remap_score(self, vmin=0, vmax=1.0, score='Kuipers'):
        from pyresample import image, geometry
        for area_name in ['antarctica',
                          'npole',
                          'ease_nh_test',
                          'ease_sh_test',
                          #'cea5km_test'
                          #'euro_arctic',
                          #'ease_world_test'
        ]:
            area_def = utils.parse_area_file(
                'reshaped_files_plotting/region_config_test.cfg',  
                area_name)[0]
            data = getattr(self, score)
            #lons = np.ma.masked_array(self.lons, mask=data.mask)
            #lats = np.ma.masked_array(self.lats, mask=data.mask)
            lons = self.lons
            lats = self.lats
            swath_def = geometry.SwathDefinition(lons=lons, lats=lats)
            swath_con = image.ImageContainerNearest(
                data, swath_def, 
                radius_of_influence=self.radius_km*1000*1.5,
                epsilon=1.0)
            area_con = swath_con.resample(area_def)
            result = area_con.image_data
            #pr.plot.show_quicklook(area_def, result,
            #                      vmin=vmin, vmax=vmax, label=score)
        
            pr.plot.save_quicklook(self.PLOT_DIR + self.figure_name + 
                                   score +'_' + area_name +'.png',
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
            lons = lons[np.not_equal(the_mask, True)]
            lats = lats[np.not_equal(the_mask, True)]
            data = data[np.not_equal(the_mask, True)]
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
            #proj2 = Basemap(projection='splaea', boundinglat=-55,
            #                lat_0=-90, lon_0=0,resolution='c')  
            #proj3 = Basemap(projection='nplaea', boundinglat=55,
            #                lat_0=90, lon_0=0,resolution='c')  
            m = proj1
            numcols=10000
            numrows=5000
            lat_min = -83.0
            lon_min = -179.9
            lat_max = 83.0
            lon_max = 179.9
            
            fig = plt.figure(figsize = (11,9))
            ax = fig.add_subplot(111)
            #m = Basemap(projection='robin',lon_0=0,resolution='c')
            m.drawparallels(np.arange(-90.,90.,30.))
            m.drawmeridians(np.arange(-180.,180.,60.))
            # transform lon / lat coordinates to map projection
            plons, plats = m(*(lons, lats))        
            # grid data
            #numcols, numrows = 1000, 500
            xi = np.linspace(lon_min, lon_max, numcols)
            yi = np.linspace(lat_min, lat_max, numrows)
            xi, yi = np.meshgrid(xi, yi)
            # interpolate
            x, y, z = (np.array(lons.ravel()), 
                       np.array(lats.ravel()), 
                       np.array(data.ravel()))
            zi = griddata((x, y), z, (xi, yi), method='nearest')
            im1 = m.pcolormesh(xi, yi, zi, cmap='coolwarm', 
                               vmin=vmin, vmax=vmax, latlon=True)
            m.drawparallels(np.arange(-90.,90.,30.))
            m.drawmeridians(np.arange(-180.,180.,60.))
            m.drawcoastlines()
            m.drawmapboundary(fill_color='0.9')
            cb = m.colorbar(im1,"right", size="5%", pad="2%")
            ax.set_title(score)
            plt.savefig(self.PLOT_DIR + self.figure_name + 
                        'basemap_' + 
                        score +'_robinson_' +area_name  +'.png')
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
        print np.sum(np.logical_and(np.not_equal(the_mask,True), Hitrate==0))
        Hitrate = np.ma.masked_array(Hitrate, mask=the_mask)
        self.Hitrate = Hitrate

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
        calipso_cfc = 100*(
            self.N_detected_clouds + self.N_undetected_clouds)*1.0/(self.N)
        the_mask = self.N<20
        calipso_cfc = np.ma.masked_array(calipso_cfc, mask=the_mask)
        self.calipso_cfc = calipso_cfc
    def calculate_bias(self):
        self.np_float_array()
        self.calculate_calipso_cfc()
        self.find_number_of_clouds_clear()
        Bias = 100*(
            self.N_detected_clouds + 
            self.N_false_clouds)*1.0/(self.N) - self.calipso_cfc
        the_mask = self.N<20
        Bias = np.ma.masked_array(Bias, mask=the_mask)
        self.Bias = Bias

    def plot_kuipers(self, PLOT_DIR="/local_disk/temp/", figure_name="figure_"):
        
        pr.plot.show_quicklook(self.definition, self.Kuipers, 
                               vmin=0, vmax=1, label="kuipers")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'kuipers_' + 
                               self.area_name +'.png',
                               self.definition, self.Kuipers, 
                               vmin=0, vmax=1, label="kuipers")
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        fig = plt.figure(figsize = (9,16))
        ax = fig.add_subplot(111)
        
        m = Basemap(projection='kav7',lon_0=0,resolution=None)
        m.drawmapboundary(fill_color='0.0')
        im1 = m.pcolormesh(self.lons.ravel(), self.lats.ravel(), self.Kuipers.ravel(),shading='flat',cmap=plt.cm.coolwarm, latlon=True)
        m.drawparallels(np.arange(-90.,99.,30.))
        m.drawmeridians(np.arange(-180.,180.,60.))
        cb = m.colorbar(im1,"right", size="5%", pad="2%")
        ax.set_title('Kuipers')
        plt.show()        
        """
        """
        import matplotlib.pyplot as plt
        bmap = pr.plot.area_def2basemap(self.definition)
        bmng = bmap.bluemarble()
        col = bmap.imshow(self.Kuipers, origin='upper', cmap="hot", label="kuipers")
        plt.colorbar(label="kuipers")
        #plt.colorbar.set_label("kuipers")
        plt.show()
    
        from trollimage.colormap import rdbu
        from trollimage.image import Image

        img = Image(self.Kuipers, mode="L")        
        rdbu.set_range(-1, 1.0)
        img.colorize(rdbu)        
        img.show()
        plt.show()
        """
    def plot_hitrate(self, PLOT_DIR = "/local_disk/", figure_name="figure_"):
        pr.plot.show_quicklook( self.definition, self.Hitrate, 
                                vmin=0, vmax=1, label="hitrate")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'hitrate_' + 
                               self.area_name +'.png',
                               self.definition, self.Hitrate, 
                               vmin=0, vmax=1, label="hitrate") 
    def plot_threat_score(self, PLOT_DIR = "/local_disk/", 
                          figure_name="figure_"):
        pr.plot.show_quicklook( self.definition, self.Threat_Score, 
                                vmin=0, vmax=1, label="threat_score")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'threat_score_' + 
                               self.area_name +'.png',
                               self.definition, self.Threat_Score, 
                               vmin=0, vmax=1, label="threat_score")
    def plot_threat_score_clear(self, PLOT_DIR = "/local_disk/", 
                                figure_name="figure_"):
        pr.plot.show_quicklook( self.definition, self.Threat_Score_Clear, 
                                vmin=0, vmax=1, label="threat_score_clear")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'threat_score_clear_' 
                               + self.area_name +'.png',
                               self.definition, self.Threat_Score_Clear, 
                               vmin=0, vmax=1, label="threat_score_clear")

class PerformancePlottingObject:
    def __init__(self):
        self.match_imager_calipso = ppsMatch_Imager_CalipsoObject()
        self.area = ppsRemappedObject()
def get_detection_stats_on_area_map_grid(my_obj, max_distance=500*1000):
    #Start with the area and get lat and lon to calculate the stats:
    area_def = my_obj.area.definition 
    lats = my_obj.area.lats[:]
    #print "lons shape", lons.shape
    if lats.ndim ==1:
        area_def = SwathDefinition(*(my_obj.area.lons,
                                     my_obj.area.lats))
        target_def = SwathDefinition(*(my_obj.match_imager_calipso.longitude, 
                                       my_obj.match_imager_calipso.latitude)) 
        valid_in, valid_out, indices, distances = get_neighbour_info(
            area_def, target_def, radius_of_influence=max_distance, 
            epsilon=1000, neighbours=1)
        cols = get_sample_from_neighbour_info('nn', target_def.shape,
                                              np.array(xrange(0,len(lats))),
                                              valid_in, valid_out,
                                              indices)
        cols = cols[valid_out]
        detected_clouds = my_obj.match_imager_calipso.detected_clouds[valid_out]
        detected_clear = my_obj.match_imager_calipso.detected_clear[valid_out]
        false_clouds = my_obj.match_imager_calipso.false_clouds[valid_out]
        undetected_clouds = my_obj.match_imager_calipso.undetected_clouds[valid_out]
        for c, ind in zip(cols.ravel(), xrange(len(cols.ravel()))):
            if distances[ind]<max_distance:
                my_obj.area.N_false_clouds[c] += false_clouds[ind]
                my_obj.area.N_detected_clouds[c] += detected_clouds[ind]
                my_obj.area.N_detected_clear[c] += detected_clear[ind]
                my_obj.area.N_undetected_clouds[c] += undetected_clouds[ind]
    else:
        target_def = SwathDefinition(*(my_obj.match_imager_calipso.longitude, 
                                      my_obj.match_imager_calipso.latitude)) 
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
        detected_clouds = my_obj.match_imager_calipso.detected_clouds[valid_out]
        detected_clear = my_obj.match_imager_calipso.detected_clear[valid_out]
        false_clouds = my_obj.match_imager_calipso.false_clouds[valid_out]
        undetected_clouds = my_obj.match_imager_calipso.undetected_clouds[valid_out]
        for r, c, ind in zip(rows.ravel(), cols.ravel(), xrange(len(cols.ravel()))):
            if distances[ind]<max_distance:
                my_obj.area.N_false_clouds[r,c] += false_clouds[ind]
                my_obj.area.N_detected_clouds[r,c] += detected_clouds[ind]
                my_obj.area.N_detected_clear[r,c] += detected_clear[ind]
                my_obj.area.N_undetected_clouds[r,c] += undetected_clouds[ind]  
    return my_obj

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
    if np.isnan(np.max(latitude)):
        raise ValueError
    return longitude, latitude
    
def get_some_info_from_caobj(my_obj, caObj, isGAC=True, isACPGv2012=False, method='KG'):
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
            caObj.calipso.all_arrays['total_optical_depth_5km']>0.2)
        isCalipsoClear = np.not_equal(isCalipsoCloudy,True)
    elif method == 'Nina' and isGAc:    
        isCalipsoCloudy = np.logical_and(
            nlay > 0, 
            caObj.calipso.all_arrays['cloud_fraction']>0.5)
        isCalipsoCloudy = np.logical_and(
            isCalipsoCloudy, 
            caObj.calipso.all_arrays['total_optical_depth_5km']>0.2)
        isCalipsoClear = np.logical_and(nlay == 0, meancl<0.01)
        isCalipsoClear = np.logical_and(
            isCalipsoClear, 
            caObj.calipso.all_arrays['total_optical_depth_5km']<0)
    elif method == 'KG':
        isCalipsoCloudy = nlay>0
        isCalipsoClear = np.not_equal(isCalipsoCloudy,True)   
    elif method == 'Nina':
        isCalipsoCloudy = nlay>0  
        isCalipsoClear = np.logical_and(nlay == 0, meancl<0.01)


    undetected_clouds = np.logical_and(isCalipsoCloudy, isClearPPS)
    false_clouds = np.logical_and(isCalipsoClear, isCloudyPPS)
    detected_clouds = np.logical_and(isCalipsoCloudy, isCloudyPPS)
    detected_clear = np.logical_and(isCalipsoClear, isClearPPS)
    use = np.logical_or(np.logical_or(detected_clouds, detected_clear),
                        np.logical_or(false_clouds, undetected_clouds))
    my_obj.match_imager_calipso.false_clouds = false_clouds[use]
    my_obj.match_imager_calipso.detected_clouds = detected_clouds[use]
    my_obj.match_imager_calipso.undetected_clouds = undetected_clouds[use]
    my_obj.match_imager_calipso.detected_clear = detected_clear[use]
    my_obj.match_imager_calipso.latitude = caObj.avhrr.latitude[use]
    my_obj.match_imager_calipso.longitude = caObj.avhrr.longitude[use]  
    return my_obj
