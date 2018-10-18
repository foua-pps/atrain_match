import h5py
import numpy as np
from datetime import datetime
from calendar import timegm
TAI93 = datetime(1993, 1, 1)
import config
from imager_cloud_products.read_cloudproducts_and_nwp_pps import NWPObj
from matchobject_io import (AmsrImagerTrackObject, 
                            AmsrObject)
from truths.calipso import (find_break_points, calipso_track_from_matched,
                     time_reshape_calipso, do_some_logging)

from utils.common import (ProcessingError, MatchupError, elements_within_range)
from libs.extract_imager_along_track import imager_track_from_matched
from utils.validate_lwp_util import LWP_THRESHOLD
import logging
logger = logging.getLogger(__name__)
AMSR_RADIUS = 5.4e3 #3.7e3 to include 5km pixels parly overlapping amsr-e footprint
def get_amsr(filename):

    if ".h5" in filename:
        retv = read_amsr_h5(filename)
    else:
        #hdf4 file:
        retv = read_amsr_hdf4(filename)

    density = 1e3 # Density of water [kg m**-3]
    n_lat_scans = len(retv.latitude)*1.0/(len(retv.sec1993)) #= 242!
    #print n_lat_scans
    epoch_diff = timegm(TAI93.utctimetuple())    
    nadir_sec_1970 = retv.sec1993 + epoch_diff
    retv.sec_1970 = np.repeat(nadir_sec_1970.ravel(), n_lat_scans)
    retv.sec1993 = None
    retv.lwp = retv.lwp_mm.ravel() * density # [mm * kg m**-3 = g m**-2]

    logger.info("Extract AMSR-E lwp between 0 and %d g/m-2", LWP_THRESHOLD)
    use_amsr = np.logical_and(retv.lwp >=0 ,
                              retv.lwp < LWP_THRESHOLD*100)
    retv = calipso_track_from_matched(retv, retv, use_amsr)
    #import matplotlib.pyplot as plt
    #plt.plot(retv.longitude, retv.latitude, '.')
    #plt.savefig('map_test.png')
    return retv 

def read_amsr_h5(filename):
    retv = AmsrObject()

    with h5py.File(filename, 'r') as f:
        #ravel AMSR-E data to 1 dimension
        retv.longitude = f['Swath1/Geolocation Fields/Longitude'][:].ravel() 
        retv.latitude = f['Swath1/Geolocation Fields/Latitude'][:].ravel()                        
        retv.sec1993 = f['Swath1/Geolocation Fields/Time']['Time'][:]
        #description='lwp (mm)',
        lwp_gain = f['Swath1/Data Fields/High_res_cloud'].attrs['Scale']#.ravel() 
        retv.lwp_mm = f['Swath1/Data Fields/High_res_cloud'][:].ravel() * lwp_gain
    if f:
        f.close() 
    return retv

def read_amsr_hdf4(filename):
    from pyhdf.SD import SD, SDC
    from pyhdf.HDF import HDF, HC
    import pyhdf.VS 

    retv = AmsrObject()
    h4file = SD(filename, SDC.READ)
    datasets = h4file.datasets()
    attributes = h4file.attributes()
    #for idx,attr in enumerate(attributes.keys()):
    #    print idx, attr
    for sds in ["Longitude", "Latitude", "High_res_cloud"]:
        data = h4file.select(sds).get()
        if sds in ["Longitude", "Latitude"]:
            retv.all_arrays[sds.lower()] = data.ravel()
        elif sds in ["High_res_cloud"]:
            lwp_gain = h4file.select(sds).attributes()['Scale']
            retv.all_arrays["lwp_mm"] = data.ravel() * lwp_gain

        #print h4file.select(sds).info()
    h4file = HDF(filename, SDC.READ)
    vs = h4file.vstart()
    data_info_list = vs.vdatainfo()
    #print "1D data compound/Vdata"
    for item in data_info_list:
        #1D data compound/Vdata
        name = item[0]
        #print name
        if name in ["Time"]:
            data_handle = vs.attach(name)
            data = np.array(data_handle[:])
            retv.all_arrays["sec1993"] = data 
            data_handle.detach()
        else:
            pass
            #print name
        #data = np.array(data_handle[:])
        #attrinfo_dic = data_handle.attrinfo()
        #factor = data_handle.findattr('factor')
        #offset = data_handle.findattr('offset')
        #print data_handle.factor
        #data_handle.detach()
    #print data_handle.attrinfo()
    h4file.close()
    #for key in retv.all_arrays.keys():
    #    print key, retv.all_arrays[key]
    return retv


def reshapeAmsr(amsrfiles, imager, SETTINGS):
    imager_end = imager.sec1970_end
    imager_start = imager.sec1970_start
    amsr = get_amsr(amsrfiles[0])
    for i in range(len(amsrfiles)-1):
        newAmsr = get_amsr(amsrfiles[i+1])
        amsr_start_all = amsr.sec_1970.ravel()
        amsr_new_all = newAmsr.sec_1970.ravel()
        if not amsr_start_all[0]<amsr_new_all[0]:
            raise ProcessingError("AMSR files are in the wrong order")
        amsr_break = np.argmin(np.abs(amsr_start_all - amsr_new_all[0]))+1
        # Concatenate the feature values
        #arname = array name from truths.amsrObj
        for arname, value in amsr.all_arrays.items(): 
            if value is not None:
                if value.size != 1:
                    amsr.all_arrays[arname] = np.concatenate(
                        (amsr.all_arrays[arname],
                         newAmsr.all_arrays[arname]),axis=0)          
    # Finds Break point
    #import pdb; pdb.set_trace()
    startBreak, endBreak = find_break_points(amsr, imager, SETTINGS)
    amsr = time_reshape_calipso(amsr, startBreak, endBreak)
    return amsr



def match_amsr_imager(amsrObj, imagerGeoObj, imagerObj, ctype, cma, ctth, nwp,
                      imagerAngObj, cpp, nwp_segments, SETTINGS):
    retv = AmsrImagerTrackObject()
    retv.imager_instrument = imagerObj.instrument.lower()

    if (getattr(cpp, "cpp_lwp")<0).all():
        logger.warning("Not matching AMSR-E with scene with no lwp.")
        return None
        #return MatchupError("No imager Lwp.") # if only LWP matching?

    from utils.common import map_imager_distances
    n_neighbours = 8
    if config.RESOLUTION == 5:
        n_neighbours = 5
    mapper_and_dist = map_imager_distances(imagerGeoObj, 
                                          amsrObj.longitude.ravel(), 
                                          amsrObj.latitude.ravel(),
                                          radius_of_influence=AMSR_RADIUS,
                                          n_neighbours=n_neighbours)
    cal, cap = mapper_and_dist["mapper"]
    distances = mapper_and_dist["distances"]
    cal_1 = cal[:,0]
    cap_1 = cap[:,0]

    calnan = np.where(cal_1 == config.NODATA, np.nan, cal_1)
    if (~np.isnan(calnan)).sum() == 0:
        logger.warning("No matches within region.")
        return None   
    #check if it is within time limits:
    if len(imagerGeoObj.time.shape)>1:
        imager_time_vector = [imagerGeoObj.time[line,pixel] for line, pixel in zip(cal_1,cap_1)]
        imager_lines_sec_1970 = np.where(cal_1 != config.NODATA, imager_time_vector, np.nan)
    else:
        imager_lines_sec_1970 = np.where(cal_1 != config.NODATA, imagerGeoObj.time[cal_1], np.nan)
    # Find all matching Amsr pixels within +/- sec_timeThr from the IMAGER data
    imager_sunz_vector = np.array([imagerAngObj.sunz.data[line,pixel] for line, pixel in zip(cal_1,cap_1)])
    idx_match = np.logical_and(
        elements_within_range(amsrObj.sec_1970, imager_lines_sec_1970, SETTINGS["sec_timeThr"]),
        imager_sunz_vector<=84) #something larger than 84 (max for lwp)

    if idx_match.sum() == 0:
        logger.warning("No light (sunz<84)  matches in region within time threshold %d s.", SETTINGS["sec_timeThr"])
        return None
    retv.amsr = calipso_track_from_matched(retv.amsr, amsrObj, idx_match)
 
    # Amsr line,pixel inside IMAGER swath (one neighbour):
    retv.amsr.imager_linnum = np.repeat(cal_1, idx_match).astype('i')
    retv.amsr.imager_pixnum = np.repeat(cap_1, idx_match).astype('i')
    retv.amsr.imager_linnum_nneigh = np.repeat(cal, idx_match, axis=0)
    retv.amsr.imager_pixnum_nneigh = np.repeat(cap, idx_match, axis=0)
    retv.amsr.imager_amsr_dist = np.repeat(distances, idx_match, axis=0)

    # Imager time
    retv.imager.sec_1970 = np.repeat(imager_lines_sec_1970, idx_match)
    retv.diff_sec_1970 = retv.amsr.sec_1970 - retv.imager.sec_1970

    do_some_logging(retv, amsrObj)
    logger.debug("Extract imager lwp along track!")
    
    nwp_small = NWPObj({'fractionofland': getattr(nwp,'fractionofland') ,
                        'landuse': getattr(nwp, 'landuse')})

    
    retv = imager_track_from_matched(retv, SETTINGS,
                                     imagerGeoObj, None, imagerAngObj, 
                                    #nwp, ctth, ctype, cma,  
                                    nwp_small, None, None, cma,
                                    cpp=cpp, nwp_segments=None,
                                    extract_some_data_for_x_neighbours=True)
    return retv
