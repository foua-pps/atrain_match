import re, glob, os
ROOT_DIR = "/home/a001865/git/atrain_match/atrain_match/"
files = glob.glob(ROOT_DIR + "/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*.py")
files = files + glob.glob(ROOT_DIR + "/*/*/*.py")


var_name_dict ={
    "readImagerData_nc": "read_imager_data_nc",
    "readImagerData_h5": "read_imager_data_h5",
    "createImagerTime": "create_imager_time",
    "Obt":"obt",
    "smallDataObject": "SmallDataObject"
}


for filename in files:
    CONCAT = 0
    linec = ""
    if "reshaped_files_scr" in filename:
        continue
    if os.path.basename(filename) in "python_edit_the_code.py":
        continue
        print("do not edit %s"%(os.path.basename(filename)))
    print(filename )
    all_file=""
    python_file = open(filename,'r') 
    line_list = [line for line in python_file]  
    for line in line_list:
        all_file += line.rstrip()+ "\n"

    
    python_file.close()
    python_file = open(filename,'w')
    python_file.write(all_file)
