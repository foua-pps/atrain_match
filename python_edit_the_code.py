import re, glob, os
ROOT_DIR = "/home/a001865/git/atrain_match/atrain_match/"
files = glob.glob(ROOT_DIR + "/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*/*.py")


var_name_dict ={
    "ppsMatch_Imager_CalipsoObject": "MatchImagerCalipso",
    "ppsStatsOnFibLatticeObject": "StatsOnFibonacciLattice",
    "N_": "num_",
    "Sum": "sum",
    "isCalipsoClear": "_calipso_clear", 
    "isCalipsoCloudy": "_calipso_cloudy",
    "isCloudyPPS": "cloudy_imager",
    "isClearPPS": "clear_imager",
}



for filename in files:
    CONCAT = 0
    linec = ""
    if "reshaped_files_scr" not in filename:
        pass
    if os.path.basename(filename) in "python_edit_the_code.py":
        continue
        print("do not edit %s"%(os.path.basename(filename)))
    print(filename )
    all_file=""
    python_file = open(filename,'r') 
    line_list = [line for line in python_file]  
    com = False
    inside_function_class_dif = False
    for line in line_list:
        for key in var_name_dict.keys():
            line = line.replace(key,var_name_dict[key])
        
        #if '"""' in line and com:
        #    com = False
        #if '"""' in line and not com:
        #    com = True
        #if line==line.lstrip() and len(line.lstrip())>=1 and  inside_function_class_dif:
        #    line = "\n\n" + line
        #    inside_function_class_dif = False
        #
        #
        #if len(line)>3 and (line[0:3]=="def" or line[0:5]=="class"):
        #    inside_function_class_dif = True



        #line = line.replace( '%d"%(', '{:d}".format(')
        # it was a one line comment
         
        all_file += line.rstrip() + "\n"


    all_file = all_file.rstrip() + "\n"
    python_file.close()
    python_file = open(filename,'w')
    python_file.write(all_file)
