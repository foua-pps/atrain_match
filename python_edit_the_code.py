import re, glob, os
ROOT_DIR = "/home/a001865/git/rename_atrainmatch/atrain_match/"
files = glob.glob(ROOT_DIR + "/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*/*.py")


var_name_dict ={
 "time":                      "profile_time",
 "utc_time":                  "profile_utc_time",
 #"cloud_top_profile":         "layer_top_altitude",
 #"cloud_top_profile_pressure":"layer_top_pressure",
 #"cloud_base_profile":        "layer_base_altitude",
 #"number_of_layers_found":    "number_layers_found",
 "elevation":                 "dem_surface_elevation",
 #"igbp":                      "igbp_surface_type",
 #"nsidc":                     "nsidc_surface_type",
 "optical_depth":             "feature_optical_depth_532"}

for filename in files:
    if os.path.basename(filename) in "python_edit_the_code.py":
        continue
        print "do not edit %s"%(os.path.basename(filename))
    print filename 
    all_file=""
    python_file = open(filename,'r') 

    for line in python_file:
        #line = line.replace("avhrr", "imager")
        #line = line.replace("AVHRR", "IMAGER")
        #line = line.replace("Avhrr", "Imager")
        line = line.replace("nnImager", "nnAvhrr")
        line = line.replace("nnavhrr", "nnimager")
        line = line.replace("NN-IMAGER", "NN-AVHRR")
        #line = re.sub(r"\.elevation", '.DEM_surface_elevation',line)
        #if re.search("alipso\.elevation",line) and 1==2:
        #    line = line.rstrip()
        #    line = re.sub(r"alipso\.elevation", 
        #                  'alipso.dem_surface_elevation',line)
        #    line = re.sub(r"alipsoObj\.elevation", 
        #                  'alipsoObj.dem_surface_elevation',line)
        #
        #    line = line + "\n"
        #
        #line = re.sub(r"nsidc", 'nsidc_surface_type',line)
        #line = re.sub(r"igbp", 'igbp_surface_type',line)
        #line = re.sub(r"number_of_layers_found", 'number_layers_found',line)
        #line = re.sub(r"cloud_top_profile_pressure", 
        #              'layer_top_pressure',line)
        #line = re.sub(r"cloud_base_profile", 
        #              'layer_base_altitude',line)
        #line = re.sub(r"cloud_top_profile", 
        #              'layer_top_altitude',line)
        #line = re.sub(r"\.optical_depth", 
        #              '.feature_optical_depth_532',line)
        #line = re.sub(r"\"optical_depth", 
        #              '"feature_optical_depth_532',line)
        #line = re.sub(r"\'optical_depth", 
        #              '\'feature_optical_depth_532',line)
        #line = re.sub(r"utc_time", 
        #              'profile_utc_time',line)
        #line = re.sub(r"time_tai", 
        #              'profile_time_tai',line)
        line = re.sub(r"feature_optical_depth_532_top_layer5km", 
                      'feature_optical_depth_532_top_layer_5km',line)

        """Maybe not do this!!
        line = re.sub(r"alipso\.time", 
                      'alipso.profile_time',line)
        line = re.sub(r"cal\.time", 
                      'cal.profile_time',line)
        line = re.sub(r"alipsoObj\.time", 
                      'alipsoObj.profile_time',line)
        """
        all_file += line 


    
    python_file.close()
    python_file = open(filename,'w')
    python_file.write(all_file)
