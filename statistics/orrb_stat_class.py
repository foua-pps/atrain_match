'''
Created on Oct 18, 2010

@author: a001696
'''
import numpy as np

class OrrbStats():
    """Abstract class for accumulating statistics from atrain_match. (What does 
    orrb stand for?)"""
    
    def __init__(self, results_files=None, ac_data=None, truth_sat='calipso'):
        """Create an OrrbStats object with results from *results_files*."""
        self.results_files = results_files
        self.truth_sat = truth_sat
        self.ac_data = ac_data
        if self.results_files is not None:
            self.accumulate_data(results_files)        
        if self.results_files is not None or ac_data is not None:
             self.do_stats()


    def read_one_file(self, datafile): 
        print datafile
        data_dict = {}
        current_datafile = open(datafile, "r")
        for line in current_datafile:
            if ":" not in line:
                continue
            line = line.replace('CALIOP','CALIPSO')    
            #old before sept 2017 files have both cloudsat and calipso data in the same files   
            #Note both CALIPSO and CALIOP where used
            if self.truth_sat.upper()  not in line:
                continue
            what, data = line.rstrip().split(':')
            data = np.array([np.float(item) for item in data.lstrip().split(" ")])
            if what in data_dict.keys():
                print what
                raise KeyError("Key should not be already in list")
            data_dict[what] = data
        current_datafile.close()
        return data_dict

    def accumulate_data(self, results_files):
        print "reading data"
        acu = {}
        acu["scenes"] = len(results_files)
        #CFC DATA
        acu["n_clear_clear_cal"] = 0
        acu["n_clear_cloudy_cal"]  = 0
        acu["n_cloudy_clear_cal"]  = 0
        acu["n_cloudy_cloudy_cal"]  = 0
        acu["n_clear_clear_cal_MODIS"]  = 0
        acu["n_clear_cloudy_cal_MODIS"]  = 0
        acu["n_cloudy_clear_cal_MODIS"]  = 0
        acu["n_cloudy_cloudy_cal_MODIS"]  = 0
        got_cloudsat_modis_flag = False
        #CTY
        acu["n_low_low"] = 0
        acu["n_low_medium"] = 0
        acu["n_low_high"] = 0
        acu["n_medium_low"] = 0
        acu["n_medium_medium"] = 0
        acu["n_medium_high"] = 0
        acu["n_high_low"] = 0
        acu["n_high_medium"] = 0
        acu["n_high_high"] = 0
        acu["n_frac_low"] = 0
        acu["n_frac_medium"] = 0
        acu["n_frac_high"] = 0
        acu["n_clear_low"] = 0
        acu["n_clear_medium"] = 0
        acu["n_clear_high"] = 0
        acu["n_low_clear"] = 0
        acu["n_medium_clear"] = 0
        acu["n_high_clear"] = 0
        acu["n_frac_clear"] = 0

        #CTH
        acu["cal_all_samples"] = {}
        acu["cal_low_samples"] = {}
        acu["cal_medium_samples"] = {}
        acu["cal_high_samples"] = {}
        acu["mean_error_csa_sum"] = {}
        acu["mean_error_cal_all_sum"] = {}
        acu["mean_error_cal_low_sum"] = {}
        acu["mean_error_cal_medium_sum"] = {}
        acu["mean_error_cal_high_sum"] = {}
        acu["rms_error_csa_sum"] = {}
        acu["rms_error_cal_all_sum"] = {}
        acu["rms_error_cal_low_sum"] = {}
        acu["rms_error_cal_medium_sum"] = {}
        acu["rms_error_cal_high_sum"] =  {}
        acu["n_missed_ctth_all"] = {}
        acu["n_missed_ctth_low"] = {}
        acu["n_missed_ctth_medium"] = {}
        acu["n_missed_ctth_high"] = {}
        acu["n_missed_cma_all"] = {}
        acu["n_missed_cma_low"] = {}
        acu["n_missed_cma_medium"] = {}
        acu["n_missed_cma_high"] = {}
        cfc_stats_labels = ["CLOUD MASK %s-IMAGER TABLE"%(self.truth_sat.upper()),
                           "CLOUD MASK %s-PPS TABLE"%(self.truth_sat.upper())]
        cfc_stats_labels_modis = ["CLOUD MASK %s-MODIS TABLE"%(self.truth_sat.upper())]
        cty_stats_labels = ["CLOUD TYPE %s-IMAGER TABLE"%(self.truth_sat.upper()),
                           "CLOUD TYPE %s-PPS TABLE"%(self.truth_sat.upper())]
        cty_stats_labels_missed = ["CLOUD TYPE %s-IMAGER TABLE MISSED"%(self.truth_sat.upper()),
                                  "CLOUD TYPE %s-PPS TABLE MISSED"%(self.truth_sat.upper())] 

        for datafile in self.results_files:
            data_dict = self.read_one_file(datafile)
            # Accumulate CALIOP/ISS/CLOUDSAT statistics CFC
            for key in data_dict.keys():
                #If reprocessing old results files make sure to extract the right lines!
                #Ie do not use CLOUDSAT info when compiling stats for CALIPSO
                if  key in cfc_stats_labels:
                    cal_data = data_dict[key]
                    cal_data[cal_data<0] = 0
                    acu["n_clear_clear_cal"] += cal_data[0]
                    acu["n_clear_cloudy_cal"] += cal_data[1]
                    acu["n_cloudy_clear_cal"] += cal_data[2]
                    acu["n_cloudy_cloudy_cal"] += cal_data[3]
                if  key in cfc_stats_labels_modis:
                    modis_data = data_dict[key]
                    modis_data[modis_data<0] = 0
                    got_cloudsat_modis_flag = True
                    acu["n_clear_clear_cal_MODIS"] += modis_data[0]
                    acu["n_clear_cloudy_cal_MODIS"] += modis_data[1]
                    acu["n_cloudy_clear_cal_MODIS"] += modis_data[2]
                    acu["n_cloudy_cloudy_cal_MODIS"] += modis_data[3]

            # Accumulate CALIOP/ISS/CLOUDSAT statistics CTY    
            for key in data_dict.keys():
                if  key in cty_stats_labels:
                    cal_data = data_dict[key]
                    cal_data[cal_data<0] = 0
                    acu["n_low_low"] += cal_data[0]
                    acu["n_low_medium"] += cal_data[1]
                    acu["n_low_high"] += cal_data[2]
                    acu["n_medium_low"] += cal_data[3]
                    acu["n_medium_medium"] += cal_data[4]
                    acu["n_medium_high"] += cal_data[5]
                    acu["n_high_low"] += cal_data[6]
                    acu["n_high_medium"] += cal_data[7]
                    acu["n_high_high"] += cal_data[8]
                    acu["n_frac_low"] += cal_data[9]
                    acu["n_frac_medium"] += cal_data[10]
                    acu["n_frac_high"] += cal_data[11]

                if  key in cty_stats_labels_missed:
                    cal_data_missed = data_dict[key] 
                    cal_data_missed[cal_data_missed<0] = 0
                    acu["n_clear_low"] += cal_data_missed[0]
                    acu["n_clear_medium"] += cal_data_missed[1]
                    acu["n_clear_high"] += cal_data_missed[2]
                    acu["n_low_clear"] += cal_data_missed[3]
                    acu["n_medium_clear"] += cal_data_missed[4]
                    acu["n_high_clear"] += cal_data_missed[5]
                    acu["n_frac_clear"] += cal_data_missed[6]

            # Accumulate CALIOP/ISS/CLOUDSAT statistics CTH
            for key in data_dict.keys(): 
                print key
                if "HEIGHT" not in key:
                    continue 
                print key
                type_of_clouds = key.split(" ")[-2]      
                cloud_level  = key.split(" ")[-1] 
                data_one_cat = data_dict[key]   
                if data_one_cat[3]<0:
                    print "no pixels!"
                    continue
                
                if cloud_level == "MEDIUM":
                    if type_of_clouds not in acu["cal_medium_samples"].keys():                        
                        acu["cal_medium_samples"][type_of_clouds] = 0
                        acu["mean_error_cal_medium_sum"][type_of_clouds] = 0
                        acu["rms_error_cal_medium_sum"][type_of_clouds] = 0
                        acu["n_missed_ctth_medium"][type_of_clouds] = 0
                        acu["n_missed_cma_medium"][type_of_clouds] = 0
                    acu["n_missed_ctth_medium"][type_of_clouds] += data_one_cat[5]
                    acu["n_missed_cma_medium"][type_of_clouds] += data_one_cat[4]                    
                    acu["cal_medium_samples"][type_of_clouds] += data_one_cat[3]
                    acu["mean_error_cal_medium_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[1]
                    acu["rms_error_cal_medium_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[2]*data_one_cat[2] 
                elif cloud_level == "HIGH":
                    if type_of_clouds not in acu["cal_high_samples"].keys():
                        acu["cal_high_samples"][type_of_clouds] = 0
                        acu["mean_error_cal_high_sum"][type_of_clouds] = 0
                        acu["rms_error_cal_high_sum"][type_of_clouds] = 0 
                        acu["n_missed_ctth_high"][type_of_clouds] = 0
                        acu["n_missed_cma_high"][type_of_clouds] = 0
                    acu["n_missed_ctth_high"][type_of_clouds] += data_one_cat[5]
                    acu["n_missed_cma_high"][type_of_clouds] += data_one_cat[4]  
                    acu["cal_high_samples"][type_of_clouds] += data_one_cat[3]
                    acu["mean_error_cal_high_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[1]
                    acu["rms_error_cal_high_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[2]*data_one_cat[2]
                elif cloud_level == "LOW":
                    if type_of_clouds not in acu["cal_low_samples"].keys():
                        acu["cal_low_samples"][type_of_clouds] = 0
                        acu["mean_error_cal_low_sum"][type_of_clouds] = 0
                        acu["rms_error_cal_low_sum"][type_of_clouds] = 0 
                        acu["n_missed_ctth_low"][type_of_clouds] = 0
                        acu["n_missed_cma_low"][type_of_clouds] = 0
                    acu["n_missed_ctth_low"][type_of_clouds] += data_one_cat[5]
                    acu["n_missed_cma_low"][type_of_clouds] += data_one_cat[4]  
                    acu["cal_low_samples"][type_of_clouds] += data_one_cat[3]
                    acu["mean_error_cal_low_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[1]
                    acu["rms_error_cal_low_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[2]*data_one_cat[2]
                else:
                    if type_of_clouds not in acu["cal_all_samples"].keys():
                        acu["cal_all_samples"][type_of_clouds] = 0
                        acu["mean_error_cal_all_sum"][type_of_clouds] = 0
                        acu["rms_error_cal_all_sum"][type_of_clouds] = 0
                        acu["n_missed_ctth_all"][type_of_clouds] = 0
                        acu["n_missed_cma_all"][type_of_clouds] = 0
                    acu["n_missed_ctth_all"][type_of_clouds] += data_one_cat[5]
                    acu["n_missed_cma_all"][type_of_clouds] += data_one_cat[4]  
                    acu["cal_all_samples"][type_of_clouds] += data_one_cat[3]
                    acu["mean_error_cal_all_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[1]
                    acu["rms_error_cal_all_sum"][type_of_clouds] += data_one_cat[3]*data_one_cat[2]*data_one_cat[2] 
        acu['got_cloudsat_modis_flag'] = got_cloudsat_modis_flag
        acu["Num"] = (acu["n_cloudy_cloudy_cal"] + acu["n_cloudy_clear_cal"] +
             acu["n_clear_cloudy_cal"] + acu["n_clear_clear_cal"])
        self.ac_data = acu
    def do_stats(self):
        """
        Calculate all statistics and put results in instance attributes (?).
        
        """
        if not self.results_files and not self.ac_data:
            raise RuntimeError("No results files and no already loaded data")
    
    def printout(self):
        """Generate nice printout of the results."""
        raise NotImplementedError("The printout method should be implemented in"
                                  " a subclass of OrrbStat.")
    
    def write(self, filename, mode='w'):
        """Write printout to a file."""
        lines_to_write = self.printout()
        if len(lines_to_write) == 0:
            print "Not writing file %s"%(filename)
            print "No compiled results!"
        else:    
            print "Writing file %s"%(filename)
            f = open(filename, mode)
            for l in lines_to_write:
                f.write(l + '\n')
            f.close()
    

