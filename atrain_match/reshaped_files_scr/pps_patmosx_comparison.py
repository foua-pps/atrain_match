
"""Read all matched data and make tables
"""
import os
from glob import glob
import numpy as np
from scipy import ndimage
from matchobject_io import (read_files)


from utils.stat_util import (my_hist,
                       my_iqr, 
                       my_rms,
                       my_mae,
                       half_sample_mode,
                       half_sample_mode,
                       my_pe250m,
                       my_pe500m,
                       my_pe1000m,
                       my_pe2000m,
                       my_pe2500m,
                       my_pe5000m)

def print_stats(cObjPPS, cObjPATMOSX):

    patmosx_profile_id = cObjPATMOSX.calipso.sec_1970#profile_id[:,0]
    pps_profile_id = cObjPPS.calipso.sec_1970#profile_id[:,0]
    
    x = cObjPPS.calipso.all_arrays['validation_height']
    x_patmosx = cObjPATMOSX.calipso.all_arrays['validation_height'] 
    y_pps = cObjPPS.imager.all_arrays['imager_ctth_m_above_seasurface']
    y_patmosx =cObjPATMOSX.imager.all_arrays['ctth_height']

    use_pps = x>-9
    use_patmosx = x_patmosx>-9
    pps_profile_id[~use_pps]=-111 #remove nodata pps
    patmosx_profile_id[~use_patmosx]=-333 #remove nodata patmosx
    pps_profile_id[y_pps<0]=-111 #remove nodata pps
    patmosx_profile_id[y_patmosx<0]=-333 #remove nodata patmosx
    
    use_patmosx_same_profile = np.array([p_id in pps_profile_id[use_pps] for p_id in patmosx_profile_id])
    use_pps_same_profile =  np.array([p_id in patmosx_profile_id[use_patmosx] for p_id in pps_profile_id])
    use_patmosx =  np.logical_and(use_patmosx, use_patmosx_same_profile)
    use_pps =   np.logical_and(use_pps, use_pps_same_profile)
    pps_profile_id[~use_pps]=-111
    #remove doubles
    unique, index = np.unique(pps_profile_id, return_index=True)
    pps_profile_id_new = -111+0*pps_profile_id.copy()
    pps_profile_id_new[index] = unique
    pps_profile_id = pps_profile_id_new 

    unique, index = np.unique(patmosx_profile_id, return_index=True)
    patmosx_profile_id_new = -333+0*patmosx_profile_id.copy()
    patmosx_profile_id_new[index] = unique
    patmosx_profile_id = patmosx_profile_id_new 

    use_patmosx_same_profile = np.array([p_id in pps_profile_id[use_pps] for p_id in patmosx_profile_id])
    use_pps_same_profile =  np.array([p_id in patmosx_profile_id[use_patmosx] for p_id in pps_profile_id])
    use_patmosx =  np.logical_and(use_patmosx, use_patmosx_same_profile)
    use_pps =   np.logical_and(use_pps, use_pps_same_profile)

    
    print(sorted(patmosx_profile_id[use_patmosx])[-10:])
    print(sorted(pps_profile_id[use_pps])[-10:])

    print(np.sum(use_pps), np.sum(use_patmosx))

    use_patmosx =  np.logical_and(use_patmosx, x_patmosx>0)
    use_pps =   np.logical_and(use_pps, x>0)

    bias_pps = y_pps[use_pps]-x[use_pps]
    bias_patmosx = y_patmosx[use_patmosx]-x_patmosx[use_patmosx]
    abias_pps = np.abs(bias_pps)
    abias_patmosx = np.abs(bias_patmosx)

    print(" {:3.1f}, ({:3.1f}), {:3.1f}, ({:3.1f}), {:3.1f}, ({:3.1f}) {:3.1f}, ({:3.1f}) {:d}, ({:d})".format(np.mean(abias_pps),
                                                                                       np.mean(abias_patmosx),
                                                                                       my_pe500m(abias_pps),
                                                                                       my_pe500m(abias_patmosx),
                                                                                       np.median(bias_pps),
                                                                                       np.median(bias_patmosx),
                                                                                       my_iqr(bias_pps),
                                                                                       my_iqr(bias_patmosx),
                                                                                       len(bias_pps),
                                                                                       len(bias_patmosx),

                                             ))
    
if __name__ == "__main__":

    PATMOSX_ROOT_DIR = ("/home/a001865/VALIDATION_PATMOSX/Reshaped_Files/noaa18/5km/2009/*/*h5")
                 
    PPS_ROOT_DIR = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/2009/5km_noaa18_2009*cali*h5")
    PPS_ROOT_DIR = ("/home/a001865/DATA_MISC/reshaped_files_validation_2018/global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/2009/5km_noaa18_2009*cali*h5")
    PPS_ROOT_DIR = ("/home/a001865/DATA_MISC/reshaped_files_cci_noaa18_2009/V2/*h5")


    patmosx_files = glob(PATMOSX_ROOT_DIR)
    pps_files = glob(PPS_ROOT_DIR)
    cObjPATMOSX =  read_files(patmosx_files)
    cObjPPS =  read_files(pps_files)
    print_stats(cObjPPS, cObjPATMOSX)
