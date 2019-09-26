import re, glob, os
ROOT_DIR = "/home/a001865/git/atrain_match/atrain_match/"
files = glob.glob(ROOT_DIR + "/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*/*.py")

var_name_dict ={
    "imagerAngObj)": "angle_obj)",
    "imagerAngObj()": "ImagerAngObj()",
    "ppsImagerObject":"ExtractedImagerObject",
    "imagerAngObj =": "angle_obj =",
    "imagerAngObj.": "angle_obj.",
     " AngObj": " angle_obj",
    "(AngObj": "(angle_obj",
    "imagerObj": "imager_obj",
    "dataObj": "imager_obj",
    "caObj": "match_calipso",
    "clsatObj": "match_clsat",
    "Obj1": "calipso1km",
    "Obj5": "calipso5km",
    "CaObj":  "match_calipso",
    "moraObj": "match_mora",
    "synopObj": "match_synop",
    "cObj": "match_obj",
    "caObjAerosol": "calipso_aerosol",
    "CalipsoCloudOpticalDepth_new":"optical_depth_height_filtering",
    "CalipsoOpticalDepthHeightFiltering":"detection_height_filtering",
    "CalipsoOpticalDepthSetThinToClearFiltering1km":"set_thin_to_clear_filtering_1km",
    "drawCalClsatImagerPlotTimeDiff":"plot_cal_clsat_imager_time_diff",
    "drawCalClsatGEOPROFImagerPlot":"plot_cal_clsat_geoprof_imager",
    "drawCalClsatImagerPlotSATZ":"plot_cal_clsat_imager_satz",
    "drawCalClsatCWCImagerPlot":"plot_cal_clsat_cwc_imager",
    "plotSatelliteTrajectory":"plot_satellite_trajectory",
    "CalculateStatistics":"calculate_statistics",
    "readTruthImagerMatchObj":"read_truth_imager_match_obj",
    "writeTruthImagerMatchObj":"write_truth_imager_match_obj",
    "reshapeCloudsat ":"reshape_cloudsat ",
    "mergeCloudsat":"merge_cloudsat",
    "reshapeAmsr":"reshape_amsr",
    "reshapeMora":"reshape_mora",
    "reshapeSynop":"reshape_synop",
    "reshapeIss":"reshape_iss",
    "reshapeCalipso":"reshape_calipso",
    "discardCalipsoFilesOutsideTimeRange":"discard_calipso_files_outside_time_range",
    "add1kmTo5km":"add_1km_to_5km",
    "addSingleShotTo5km":"add_singleshot_to5km",
    "add5kmVariablesTo1kmresolution":"add_5km_variables_to_1km",
    "adjust5kmTo1kmresolution":"adjust_5km_to_1km_resolution",
}


for filename in files:
    CONCAT = 0
    linec = ""
    #if "reshaped_files_scr" in filename:
    #    continue
    if os.path.basename(filename) in "python_edit_the_code.py":
        continue
        print("do not edit %s"%(os.path.basename(filename)))
    print(filename )
    all_file=""
    python_file = open(filename,'r') 
    line_list = [line for line in python_file]  
    for line in line_list:
        for key in var_name_dict.keys():
            line = line.replace(key,var_name_dict[key])

        all_file += line

    
    python_file.close()
    python_file = open(filename,'w')
    python_file.write(all_file)
