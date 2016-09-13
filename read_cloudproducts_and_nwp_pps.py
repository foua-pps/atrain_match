import numpy as np
import logging
logger = logging.getLogger(__name__)
#logger.debug('Just so you know: this module has a logger...')
class NWPObj(object):
    def __init__(self, array_dict):
        self.surftemp = None
        self.t500 = None
        self.t700 = None
        self.t850 = None
        self.t950 = None
        self.ttro = None
        self.ciwv = None
        self.text_r06 = None
        self.text_t11 = None
        self.text_t37t12 = None
        self.text_t37 = None
        self.thr_t11ts_inv = None
        self.thr_t11t37_inv = None
        self.thr_t37t12_inv = None
        self.thr_t11t12_inv = None
        self.thr_t85t11_inv = None
        self.thr_t11ts = None
        self.thr_t11t37 = None
        self.thr_t37t12 = None
        self.thr_t11t12 = None
        self.thr_t85t11 = None
        self.thr_r06 = None
        self.thr_r09 = None
        self.emis1 = None
        self.emis6 = None
        self.emis8 = None
        self.emis9 = None
        self.__dict__.update(array_dict) 

from config import (VAL_CPP)
def readCpp(filename, cpp_type):
    import h5py 
    h5file = h5py.File(filename, 'r')
    if cpp_type in h5file.keys():
        value = h5file[cpp_type].value
        gain = h5file[cpp_type].attrs['gain']
        intercept = h5file[cpp_type].attrs['intercept']
        nodat = h5file[cpp_type].attrs['no_data_value']
        product = np.where(value != nodat,value * gain + intercept, value)   
    h5file.close()
    return product

def read_nwp_h5(filename, type_of_nwp,nwp_key):
    import h5py 
    h5file = h5py.File(filename, 'r')
    if nwp_key in h5file.keys():
        logger.info("Read NWP %s"%(type_of_nwp))
        value = h5file[nwp_key].value
        gain = h5file[nwp_key].attrs['gain']
        intercept = h5file[nwp_key].attrs['intercept']
        nodat = h5file[nwp_key].attrs['nodata']
        return  np.where(value != nodat,value * gain + intercept, value)
    else:
        logger.info("NO NWP %s File, Continue"%(type_of_nwp))
        return None

def read_segment_data(filename):
    import h5py
    product = {}
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        for attribute in list(h5file.attrs):
            product[attribute] = h5file.attrs[attribute]
        for attribute in list(h5file['satdef'].attrs):
            product[attribute] = h5file['satdef'].attrs[attribute]
        product['colidx'] = h5file['satdef']['colidx'].value
        product['rowidx'] = h5file['satdef']['rowidx'].value
        logger.info("Read segment info moisture")
        for moist_data in ['moist', 'surfaceMoist']:
            data = h5file['segment_area'][moist_data]
            gain = h5file.attrs['m_gain']
            intercept = h5file.attrs['m_intercept']
            nodata = h5file.attrs['m_nodata']
            #data[data!=nodata] = data[data!=nodata] * (gain) + intercept
            product[moist_data] = data
        logger.info("Read segment info pressure")
        for pressure_data in ['pressure', 'surfacePressure', 'ptro']:
            #pressure is in Pa in segments file
            data = h5file['segment_area'][pressure_data]
            gain = h5file.attrs['p_gain']
            intercept = h5file.attrs['p_intercept']
            nodata = h5file.attrs['p_nodata']
            data[data!=nodata] = data[data!=nodata] * (gain/100) + intercept/100 #Pa => hPa
            product[pressure_data] = data
        logger.info("Read segment info height")
        for geoheight_data in ['geoheight', 'surfaceGeoHeight']:
            #geo height is in meters in segment file
            data = h5file['segment_area'][geoheight_data]
            gain = h5file.attrs['h_gain']
            intercept = h5file.attrs['h_intercept']
            nodata = h5file.attrs['h_nodata']
            data[data!=nodata] = data[data!=nodata] * gain + intercept
            product[geoheight_data] = data
        logger.info("Read segment info temperature")
        for temperature_data in ['temp']:
            # temperature are measured in Kelvin in segmentfile
            data = h5file['segment_area'][temperature_data]
            gain = h5file.attrs['t_gain']
            intercept = h5file.attrs['t_intercept']
            nodata = h5file.attrs['t_nodata']
            data[data!=nodata] = data[data!=nodata] * gain + intercept
            product[temperature_data] = data
        for temperature_data in ['t850', 'ttro', 'surfaceLandTemp', 'surfaceSeaTemp']:
            data = h5file['segment_area'][temperature_data]
            gain = h5file.attrs['t_gain']
            intercept = h5file.attrs['t_intercept']
            nodata = h5file.attrs['t_nodata']
            data_float = np.array(data, dtype=np.float)
            data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
            product[temperature_data] = data_float
        for misc_data in ['meanElevation', 'fractionOfLand']:
            product[misc_data] = h5file['segment_area'][misc_data]
        logger.info("Read segment info brightness temperature")
        try:
            for tb_data in ['tb11clfree_sea',
                            'tb12clfree_sea',
                            'tb11clfree_land',
                            'tb12clfree_land']:
                data = h5file['segment_area'][tb_data]
                gain = h5file.attrs['t_gain']
                intercept = h5file.attrs['tb_intercept']
                nodata = h5file.attrs['t_nodata']
                data_float = np.array(data, dtype=np.float)
                data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
                product[tb_data] = data_float
        except ValueError:
            for tb_data, h5name in zip(['tb11clfree_sea',
                                        'tb12clfree_sea',
                                        'tb11clfree_land',
                                        'tb12clfree_land'],
                                       ['tb4clfree_sea',
                                        'tb5clfree_sea',
                                        'tb4clfree_land',
                                       'tb5clfree_land']):
                data = h5file['segment_area'][h5name]
                gain = h5file.attrs['t_gain']
                intercept = h5file.attrs['tb_intercept']
                nodata = h5file.attrs['t_nodata']
                data_float = np.array(data, dtype=np.float)
                data_float[data!=nodata] = data_float[data!=nodata] * gain + intercept
                product[tb_data] = data_float
        h5file.close()
        return product
    else:
        logger.info("NO segment %s File, Continue"%(filename))
        return None


def read_thr_h5(filename, h5_obj_type, thr_type):
    import h5py 
    product = None
    if thr_type in ["emis1","emis6", "emis8", "emis9"]:
        if filename is not None: 
            h5file = h5py.File(filename, 'r')
            if 1==1:#h5_obj_type in h5file.keys():
                value = h5file[h5_obj_type].value
                gain = h5file.attrs['gain']
                intercept = h5file.attrs['intercept']
                product = value * gain + intercept
                logger.info("Read EMIS: %s"%(thr_type))
            else:
                logger.info("ERROR","Could not read %s File, Continue"%(thr_type))
            h5file.close()   
        else:
            logger.info("NO EMIS %s File, Continue"%(thr_type))
        return product  
    if filename is not None: 
        h5file = h5py.File(filename, 'r')
        if h5_obj_type in h5file.keys():
            value = h5file[h5_obj_type].value
            gain = h5file[h5_obj_type].attrs['gain']
            intercept = h5file[h5_obj_type].attrs['intercept']
            product = value * gain + intercept
            logger.info("Read THR: %s"%(thr_type))
        else:
            logger.info("ERROR","Could not read %s File, Continue"%(thr_type))
        h5file.close()   
    else:
        logger.info("NO THR %s File, Continue"%(thr_type))
    return product

def readViirsData_h5(filename):
    import h5py
    from pps_io import NewAvhrrData, AvhrrChannelData
    avhrrdata = NewAvhrrData()
    ich=0
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        nof_images = h5file['what'].attrs['sets']
        for num in xrange(1, nof_images+1,1):
            image = "image%d"%(num)
            channel = h5file[image].attrs['channel']
            if channel in ["M05",
                           "M07",
                           "M10",
                           "M12",
                           "M14",
                           "M15",
                           "M16",
                           "M11",
                           "M09"]:
                      
                avhrrdata.channel.append(AvhrrChannelData())
                avhrrdata.channel[ich].data = h5file[image]['data'].value
                avhrrdata.channel[ich].des = "VIIRS %s"%(channel)
                avhrrdata.channel[ich].gain = h5file[image]['what'].attrs['gain']
                avhrrdata.channel[ich].intercept = h5file[image]['what'].attrs['offset']
                avhrrdata.nodata = h5file[image]['what'].attrs['nodata']
                avhrrdata.missingdata = h5file[image]['what'].attrs['missingdata']
                avhrrdata.no_data = h5file[image]['what'].attrs['nodata']
                avhrrdata.missing_data = h5file[image]['what'].attrs['missingdata']
                ich = ich + 1
    h5file.close()
    return avhrrdata

def readModisData_h5(filename):
    import h5py
    from pps_io import NewImagerData, ImagerChannelData
    avhrrdata = NewImagerData()
    ich=0
    if filename is not None:
        h5file = h5py.File(filename, 'r')
        nof_images = h5file['what'].attrs['sets']
        for num in xrange(1, nof_images+1,1):
            image = "image%d"%(num)
            channel = h5file[image].attrs['channel']
            if channel in [
                    "1", "2", "3", "4", 
                    "5", "6", "7", "8", 
                    "9", "10", "11", "12", 
                    "13lo", "13hi", "14lo", "14hi", 
                    "15", "16", "17", "18",
                    "19", "20", "21", "22",
                    "23", "24", "25", "26",
                    "27", "28", "29", "30",
                    "31", "32", "33", "34",
                    "35", "36"]:
                avhrrdata.channel.append(ImagerChannelData())
                avhrrdata.channel[ich].data = h5file[image]['data'].value
                avhrrdata.channel[ich].des = "MODIS %s"%(channel)
                avhrrdata.channel[ich].gain = h5file[image]['what'].attrs['gain']
                avhrrdata.channel[ich].intercept = h5file[image]['what'].attrs['offset']
                avhrrdata.nodata = h5file[image]['what'].attrs['nodata']
                avhrrdata.missingdata = h5file[image]['what'].attrs['missingdata']
                avhrrdata.no_data = h5file[image]['what'].attrs['nodata']
                avhrrdata.missing_data = h5file[image]['what'].attrs['missingdata']
                ich = ich + 1
    h5file.close()
    return avhrrdata

def pps_read_all(pps_files, avhrr_file, cross):
    import pps_io #@UnresolvedImport
    import epshdf #@UnresolvedImport
    logger.info("Read AVHRR geolocation data")
    avhrrGeoObj = pps_io.readAvhrrGeoData(avhrr_file)    
    logger.info("Read AVHRR Sun -and Satellites Angles data")
    avhrrAngObj = pps_io.readSunSatAngles(pps_files.sunsatangles) #, withAbsoluteAzimuthAngles=True)
    if 'npp' in [cross.satellite1, cross.satellite2]:
        logger.info("Read VIIRS data")
        avhrrObj = pps_io.readViirsData(avhrr_file)
        #avhrrObj = readViirsData_h5(avhrr_file)
    elif (cross.satellite1 in ['eos1','eos2'] or 
        cross.satellite2 in ['eos1', 'eos2']):
        #avhrrObj = pps_io.readModisData(avhrr_file)
        avhrrObj = readModisData_h5(avhrr_file)
    else:
        logger.info("Read AVHRR data")
        avhrrObj = pps_io.readAvhrrData(avhrr_file)    
    cppLwp = None
    cppCph = None
    if VAL_CPP:    
        logger.info("Read CPP data")
        from ppshdf_cloudproducts import CppProducts #@UnresolvedImport
        if PPS_FORMAT_2012_OR_EARLIER:
            try:
                cpp = CppProducts.from_h5(pps_files.cpp,
                                          product_names=['cph','lwp'])
                cppLwp = cpp.products['lwp'].array
                cppCph = cpp.products['cph'].array
                logger.info("CPP chp and lwp data read")
            except KeyError:
                #import traceback
                #traceback.print_exc()
                cppLwp = readCpp(pps_files.cpp, 'lwp')
                cppCph = readCpp(pps_files.cpp, 'cph')
                logger.info("CPP lwp, cph data read")
        else:
            cpp = CppProducts.from_h5(pps_files.cpp,
                                      product_names=['cpp_phase','cpp_lwp'],
                                      scale_up=True)
            # LWP convert from kg/m2 to g/m2
            cppLwp = 1000. * cpp.products['cpp_lwp'].array
            cppCph = cpp.products['cpp_phase'].array
            
    logger.info("Read PPS Cloud Type")
    try:
        ctype = epshdf.read_cloudtype(pps_files.cloudtype, 1, 1, 1)  
    except:
        logger.info("Could not use pps program to read, use mpop instead")    
        #read with mpop instead
        from mpop.satin.nwcsaf_pps import NwcSafPpsChannel
        ctype_mpop = NwcSafPpsChannel()
        print ctype_mpop
        ctype_mpop.read(pps_files.cloudtype)
        logger.info("Done reading cloudtype") 
        print vars(ctype_mpop).keys()
        for varname in vars(ctype_mpop).keys():
            logger.info(varname) 
        logger.info("Done reading cloudtype") 
        #need to make empty ctype object here
        ctype_mpop.ct_quality.data
        #ctype=ctype_mpop
        ctype.cloudtype = ctype_mpop.cloudtype.data
        ctype.quality_flag = ctype_mpop.ct_quality.data
        ctype.conditions_flag = ctype_mpop.ct_conditions.data
        ctype.status_flag = ctype_mpop.ct_statusflag.data
        #print ctype_mpop.ct_quality

    logger.info("Read PPS CTTH")
    try:
        ctth = epshdf.read_cloudtop(pps_files.ctth, 1, 1, 1, 0, 1)
    except:
        logger.info("No CTTH")
        ctth = None

    logger.info("Read PPS NWP data")   
    nwp_dict={}
    nwp_dict['surftemp'] = read_nwp_h5(pps_files.nwp_tsur, "surface temperature" ,"tsur")
    nwp_dict['t500'] = read_nwp_h5(pps_files.nwp_t500, "temperature 500hPa","t500")
    nwp_dict['t700'] = read_nwp_h5(pps_files.nwp_t700, "temperature 700HPa","t700")
    nwp_dict['t850'] = read_nwp_h5(pps_files.nwp_t850, "temperature 850hPa","t850")
    nwp_dict['t950'] = read_nwp_h5(pps_files.nwp_t950, "temperature 950hPa","t950")
    nwp_dict['ttro'] = read_nwp_h5(pps_files.nwp_ttro, "tropopause temperature","ttro")
    nwp_dict['ciwv'] = read_nwp_h5(pps_files.nwp_ciwv, 
                                "atmosphere_water_vapor_content","ciwv")

    logger.info("Read PPS threshold data")  
    for ttype in ['r06', 't11', 't37t12', 't37']:
        h5_obj_type = ttype +'_text'
        text_type = 'text_' + ttype
        nwp_dict[text_type] = read_thr_h5(getattr(pps_files,text_type), 
                                       h5_obj_type,text_type)
    for h5_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv', 
                        't11ts', 't11t37', 't37t12', 't11t12',
                        'r09', 'r06', 't85t11_inv', 't85t11']:
        thr_type = 'thr_' + h5_obj_type
        nwp_dict[thr_type] = read_thr_h5(getattr(pps_files,thr_type), 
                                      h5_obj_type, thr_type)
    logger.info("Read PPS Emissivity data")  
    for h5_obj_type in ['emis1',"emis6", 'emis8','emis9']:
        emis_type = h5_obj_type
        nwp_dict[emis_type] = read_thr_h5(getattr(pps_files,"emis"), 
                                       h5_obj_type, emis_type)
    nwp_obj = NWPObj(nwp_dict)
    logger.info("Read PPS NWP segment resolution data") 
    segment_data_object = read_segment_data(getattr(pps_files,'nwp_segments'))
    return avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, nwp_obj, cppLwp, cppCph, segment_data_object 
