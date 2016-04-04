import numpy as np
from pyresample import utils
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import get_neighbour_info
from pyresample.kd_tree import get_sample_from_neighbour_info
import pyresample as pr
from scipy import ndimage

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
    def set_area(self, area_name = 'euro_arctic_test'):     
        area_def = utils.parse_area_file(
            'reshaped_files_plotting/region_config_test.cfg',  
            area_name)[0]
        self.area_name=area_name
        self.definition = area_def
        self.lons = area_def.lons[:]
        self.lats = area_def.lats[:]
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
    def calculate_kuipers(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        N_clear = self.N_clear
        N_clouds = self.N_clouds
        N_detected_clouds = self.N_detected_clouds
        N_detected_clear = self.N_detected_clear
        """
        upd = self.N_clouds < 0.01* self.N_clear
        N_detected_clouds[upd] =  N_detected_clouds[upd] + 0.01* self.N_clear[upd]
        N_clouds[upd] = N_clouds[upd] + 0.01* self.N_clear[upd]
        upd =  self.N_clear < 0.01* self.N_clouds
        N_detected_clear[upd] =  N_detected_clear[upd] + 0.01* self.N_clouds[upd]
        N_clear[upd] = N_clear[upd] + 0.01* self.N_clouds[upd]
        """
        #Typically we have N_clear/N_clouds = 30/70 
        #In areas with only clouds or only clears the Kuipers will be ==0
        #Even if all clouds/clears are classified correctly!
        #Do something for these set to none or update   
        Kuipers_devider = (N_clouds)*(N_clear)
        Kuipers_devider[Kuipers_devider==0] = 1.0
        Kuipers = (N_detected_clouds*N_detected_clear - 
                   self.N_false_clouds*self.N_undetected_clouds)/Kuipers_devider
        Kuipers[self.N_clouds < 0.01* self.N_clear]=0.0
        Kuipers[self.N_clear < 0.01* self.N_clouds]=0.0
        Kuipers = np.ma.masked_array(Kuipers, mask=Kuipers==0)
        self.Kuipers = Kuipers

    def calculate_hitrate(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        Hitrate = (
            self.N_detected_clouds + self.N_detected_clear)*1.0/(
                self.N_clear + self.N_clouds)
        Hitrate = np.ma.masked_array(Hitrate, mask=Hitrate==0)
        self.Hitrate = Hitrate

    def calculate_threat_score(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        ThreatScore = (
            self.N_detected_clouds)*1.0/( self.N_clouds + self.N_false_clouds)
        ThreatScore = np.ma.masked_array(ThreatScore, mask=ThreatScore==0)
        self.Threat_Score = ThreatScore
    def calculate_threat_score_clear(self):
        self.np_float_array()
        self.find_number_of_clouds_clear()
        ThreatScoreClear = (
            self.N_detected_clear)*1.0/( self.N_clear + self.N_undetected_clouds)
        ThreatScoreClea = np.ma.masked_array(ThreatScoreClear, 
                                             mask=ThreatScoreClear==0)
        self.Threat_Score_Clear = ThreatScoreClear
    def plot_kuipers(self, PLOT_DIR="/local_disk/temp/", figure_name="figure_"):
        pr.plot.show_quicklook(self.definition, self.Kuipers, 
                               vmin=0, vmax=1, label="kuipers")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'kuipers_' + self.area_name +'.png',
                               self.definition, self.Kuipers, 
                               vmin=0, vmax=1, label="kuipers")
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
        """

    def plot_hitrate(self, PLOT_DIR = "/local_disk/temp/", figure_name="figure_"):
        pr.plot.show_quicklook( self.definition, self.Hitrate, 
                                vmin=0, vmax=1, label="hitrate")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'hitrate_' + self.area_name +'.png',
                               self.definition, self.Hitrate, 
                               vmin=0, vmax=1, label="hitrate") 
    def plot_threat_score(self, PLOT_DIR = "/local_disk/temp/", figure_name="figure_"):
        pr.plot.show_quicklook( self.definition, self.Threat_Score, 
                                vmin=0, vmax=1, label="threat_score")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'threat_score_' + self.area_name +'.png',
                               self.definition, self.Threat_Score, 
                               vmin=0, vmax=1, label="threat_score")
    def plot_threat_score_clear(self, PLOT_DIR = "/local_disk/temp/", figure_name="figure_"):
        pr.plot.show_quicklook( self.definition, self.Threat_Score_Clear, 
                                vmin=0, vmax=1, label="threat_score_clear")
        pr.plot.save_quicklook(PLOT_DIR + figure_name + 'threat_score_clear_' + self.area_name +'.png',
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
    """
    print len(distances[distances<max_distance])
    print caObj.avhrr.longitude.shape #avhrr
    print indices.shape
    print cols.shape, np.max(cols), cols, 
    print np.min(distances) #avhrr
    print rows.shape, np.max(rows) #avhrr
    print valid_in.shape, sum(valid_in)
    print valid_out.shape, sum(valid_out) #avhrr
    print 512*512
    """
    detected_clouds = my_obj.match_imager_calipso.detected_clouds[valid_out]
    detected_clear = my_obj.match_imager_calipso.detected_clear[valid_out]
    false_clouds = my_obj.match_imager_calipso.false_clouds[valid_out]
    undetected_clouds = my_obj.match_imager_calipso.undetected_clouds[valid_out]
    for r, c, ind in zip(rows.ravel(), cols.ravel(), xrange(len(cols.ravel()))):
        #print r, c, ind
        #print lats[r,c]
        #print lons[r,c]
        #print detected_clouds[ind]
        #print distances[ind]
        if distances[ind]<max_distance:
            my_obj.area.N_false_clouds[r,c] += false_clouds[ind]
            my_obj.area.N_detected_clouds[r,c] += detected_clouds[ind]
            my_obj.area.N_detected_clear[r,c] += detected_clear[ind]
            my_obj.area.N_undetected_clouds[r,c] += undetected_clouds[ind]
    return my_obj

def get_some_info_from_caobj(my_obj, caObj, isGAC=True, isACPGv2012=False):
    #cloudObj = get_clear_and_cloudy_vectors(caObj, isACPGv2012, isGAC)  
    isCloudyPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>4,
                                 caObj.avhrr.all_arrays['cloudtype']<21) 
    isClearPPS = np.logical_and(caObj.avhrr.all_arrays['cloudtype']>0,
                                caObj.avhrr.all_arrays['cloudtype']<5)
   
    print "should be zero", np.sum(np.logical_and(caObj.calipso.all_arrays['number_layers_found']==0,
                                                  caObj.calipso.all_arrays['layer_top_altitude'][::,0]>0))
    nlay =np.where(caObj.calipso.all_arrays['number_layers_found']>0,1,0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=10)
    isCalipsoCloudy = np.logical_and(nlay > 0, meancl>0.0)
    isCalipsoCloudy = np.logical_and(isCalipsoCloudy, 
                                     caObj.calipso.all_arrays['total_optical_depth_5km']>2.0)
    #isCalipsoCloudy = np.logical_and(nlay > 0, caObj.calipso.all_arrays['cloud_fraction']>0.2)
    #isCalipsoCloudy = np.logical_and(isCalipsoCloudy, 
    #                                 caObj.calipso.all_arrays['total_optical_depth_5km']>0.5)


    isCalipsoClear = np.logical_and(nlay == 0, meancl<0.1)
    isCalipsoClear = np.logical_and(isCalipsoClear, 
                                    caObj.calipso.all_arrays['total_optical_depth_5km']<0)


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
