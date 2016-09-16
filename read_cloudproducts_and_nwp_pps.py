import numpy as np
import netCDF4
import h5py
import logging
import time
import calendar
logger = logging.getLogger(__name__)
from config import (VAL_CPP)
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

class smallDataObject(object):
    def __init__(self):
        self.data=None
        self.gain = 1.0
        self.intercept = 0.0
        self.no_data = 0.0
        self.missing_data = 0.0
class imagerAngObj(object):
    def __init__(self):        
        self.satz = smallDataObject()
        self.sunz = smallDataObject()
        self.azidiff = smallDataObject()

class imagerGeoObj(object):
    def __init__(self):
        self.latitude = None
        self.longitude = None
        self.azidiff = None
        self.num_of_lines = None

class NewImagerData:
    def __init__(self):
        self.des = ""
        self.no_data = -999
        self.missing_data = -999
        self.nodata = -999
        self.missingdata = -999
        self.channel = []

class ImagerChannelData:
    def __init__(self):
        self.channel = ""
        self.des = ""
        self.gain = 1.0
        self.intercept = 0.0
        self.data = None

class CtypeObj:
    #skeleton container for v2014 cloudtype
    def __init__(self):
        self.cloudtype = None
        self.ct_statusflag = None
        self.ct_quality = None
        self.ct_conditions = None
        self.phaseflag = None #v2012
class CtthObj:
    #skeleton container for v2014 cloudtype
    def __init__(self):
        self.height = None
        self.temperature = None
        self.pressure = None
        self.ctth_statusflag = None

def read_ctth_h5(filename):
    h5file = h5py.File(filename, 'r')
    ctth = CtthObj()
    ctth.height = h5file['ctth_alti'].value
    ctth.temperature = h5file['ctth_tempe'].value
    ctth.pressure = h5file['ctth_pres'].value
    ctth.ctth_statusflag = h5file['ctth_status_flag'].value
    #Currently unpacked arrays later in calipso.py
    #TODO: move this here also for h5!
    ctth.h_gain = h5file['ctth_alti'].attrs['scale_factor']
    ctth.h_intercept = h5file['ctth_alti'].attrs['add_offset']
    ctth.t_gain = h5file['ctth_tempe'].attrs['scale_factor']
    ctth.t_intercept = h5file['ctth_tempe'].attrs['add_offset']
    ctth.p_gain = h5file['ctth_pres'].attrs['scale_factor']
    ctth.p_intercept = h5file['ctth_pres'].attrs['add_offset']
    ctth.h_nodata = h5file['ctth_alti'].attrs['_FillValue']
    ctth.t_nodata = h5file['ctth_tempe'].attrs['_FillValue']
    ctth.p_nodata = h5file['ctth_pres'].attrs['_FillValue']
    return ctth

def read_ctth_nc(filename):
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    ctth = CtthObj()
    ctth.height = pps_nc.variables['ctth_alti'][0,:,:]
    ctth.temperature = pps_nc.variables['ctth_tempe'][0,:,:]
    ctth.pressure = pps_nc.variables['ctth_pres'][0,:,:]
    ctth.ctth_statusflag = pps_nc.variables['ctth_status_flag'][0,:,:]
    #Currently unpacked arrays later in calipso.py
    #TODO: move this here also for h5!
    ctth.h_gain=1.0
    ctth.h_intercept=0.0
    ctth.p_gain=1.0
    ctth.p_intercept=0.0
    ctth.t_gain=1.0
    ctth.t_intercept=0.0
    ctth.h_nodata = pps_nc.variables['ctth_alti']._FillValue
    ctth.t_nodata = pps_nc.variables['ctth_tempe']._FillValue
    ctth.p_nodata = pps_nc.variables['ctth_pres']._FillValue
    return ctth

def read_cloudtype_h5(filename):
    h5file = h5py.File(filename, 'r')
    ctype = CtypeObj()
    ctype.cloudtype = h5file['ct'].value
    ctype.ct_conditions = h5file['ct_conditions'].value
    ctype.ct_statusflag = h5file['ct_status_flag'].value
    ctype.ct_quality = h5file['ct_quality'].value
    return ctype

def read_cloudtype_nc(filename):
    pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    ctype = CtypeObj()
    ctype.cloudtype = pps_nc.variables['ct'][0,:,:]
    ctype.ct_conditions = pps_nc.variables['ct_conditions'][0,:,:]
    ctype.ct_statusflag = pps_nc.variables['ct_status_flag'][0,:,:]
    ctype.ct_quality = pps_nc.variables['ct_quality'][0,:,:]
    return ctype


def readImagerData_nc(pps_nc):
    imager_data = NewImagerData()
    for var in pps_nc.variables.keys():
        if 'image' in var:
            image = pps_nc.variables[var]
            logger.info("reading channel %s", image.description)
            one_channel = ImagerChannelData()
            #channel = image.channel                     
            one_channel.data = image[0,:,:]
            one_channel.des = image.description
            #Currently unpacked arrays later in calipso.py:
            #TODO: move this herealso for h5!
            one_channel.gain = 1.0 #data already scaled
            one_channel.intercept = 0.0 #data already scaled
            imager_data.channel.append(one_channel) 
            fill_value = pps_nc.variables[var]._FillValue
            imager_data.nodata = fill_value
            imager_data.missingdata = fill_value
            imager_data.no_data = fill_value
            imager_data.missing_data = fill_value
    return imager_data

def read_pps_angobj_nc(pps_nc):
    """Read angles info from filename
    """
    AngObj = imagerAngObj()
    AngObj.satz.data = pps_nc.variables['satzenith'][0,:,:]
    AngObj.sunz.data = pps_nc.variables['sunzenith'][0,:,:]
    AngObj.azidiff.data = pps_nc.variables['azimuthdiff'][0,:,:]
    AngObj.satz.no_data = pps_nc.variables['satzenith']._FillValue
    AngObj.sunz.no_data = pps_nc.variables['sunzenith']._FillValue
    AngObj.azidiff.no_data = pps_nc.variables['azimuthdiff']._FillValue
    AngObj.satz.missing_data = pps_nc.variables['satzenith']._FillValue
    AngObj.sunz.missing_data = pps_nc.variables['sunzenith']._FillValue
    AngObj.azidiff.missing_data = pps_nc.variables['azimuthdiff']._FillValue
    return AngObj

def read_pps_angobj_h5(filename):
    """Read angles info from filename
    """
    h5file = h5py.File(filename, 'r')
    AngObj = imagerAngObj()    
    for var in h5file.keys():
        if 'image' in var:
            image = h5file[var]            
            if image.attrs['description'] == "Solar zenith angle":
                AngObj.sunz.data = image['data'].value
                AngObj.sunz.gain = image['what'].attrs['gain']
                AngObj.sunz.intercept = image['what'].attrs['offset']
                AngObj.sunz.no_data = image['what'].attrs['nodata']
                AngObj.sunz.missing_data = image['what'].attrs['missingdata']
            elif image.attrs['description'] == "Satellite zenith angle":
                AngObj.satz.data = image['data'].value
                AngObj.satz.gain = image['what'].attrs['gain']
                AngObj.satz.intercept = image['what'].attrs['offset']
                AngObj.satz.no_data = image['what'].attrs['nodata']
                AngObj.satz.missing_data = image['what'].attrs['missingdata']
            elif (image.attrs['description'] == 
                  "Relative satellite-sun azimuth angle"):
                AngObj.azidiff.data = image['data'].value
                AngObj.azidiff.gain = image['what'].attrs['gain']
                AngObj.azidiff.intercept = image['what'].attrs['offset']
                AngObj.azidiff.no_data = image['what'].attrs['nodata']
                AngObj.azidiff.missing_data = image['what'].attrs['missingdata']

    return AngObj

def read_pps_geoobj_nc(pps_nc):
    """Read geolocation and time info from filename
    """
    GeoObj = imagerGeoObj()
    GeoObj.longitude = pps_nc.variables['lon'][::]
    GeoObj.nodata = pps_nc.variables['lon']._FillValue
    GeoObj.latitude = pps_nc.variables['lat'][::]
    GeoObj.num_of_lines = GeoObj.latitude.shape[0]
    time_temp = pps_nc.variables['time'].units #to 1970 s
    time_obj = time.strptime(time_temp,'seconds since %Y-%m-%d %H:%M:%S.%f +00:00')
    sec_since_1970 = calendar.timegm(time_obj)
    GeoObj.sec1970_start = sec_since_1970 + np.min(pps_nc.variables['time_bnds'][::])
    GeoObj.sec1970_end = sec_since_1970 + np.max(pps_nc.variables['time_bnds'][::])
    GeoObj.sec1970_start = int(GeoObj.sec1970_start) 
    GeoObj.sec1970_end = int(GeoObj.sec1970_end)
    tim1 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(GeoObj.sec1970_start))
    tim2 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(GeoObj.sec1970_end))
    logger.info("Starttime: %s, end time: %s"%(tim1, tim2))
    logger.info("Min lon: %s, max lon: %d"%(
        np.min(GeoObj.longitude),
        np.max(GeoObj.longitude)))
    logger.info("Min lat: %d, max lat: %d"%(
        np.min(GeoObj.latitude),np.max(GeoObj.latitude)))
    return  GeoObj

def read_pps_geoobj_h5(filename):
    """Read geolocation and time info from filename
    """
    h5file = h5py.File(filename, 'r')
    GeoObj = imagerGeoObj()
    in_fillvalue1 = h5file['where/lon/what'].attrs['nodata']
    in_fillvalue2 = h5file['where/lon/what'].attrs['missingdata']
    GeoObj.nodata = -999.0
    gain = h5file['where/lon/what'].attrs['gain']
    intercept = h5file['where/lon/what'].attrs['offset']
    GeoObj.longitude = h5file['where/lon']['data'].value*gain + intercept
    GeoObj.latitude = h5file['where/lat']['data'].value*gain + intercept

    GeoObj.longitude[h5file['where/lon']['data'].value==in_fillvalue1] = GeoObj.nodata
    GeoObj.latitude[h5file['where/lat']['data'].value==in_fillvalue1] = GeoObj.nodata
    GeoObj.longitude[h5file['where/lon']['data'].value==in_fillvalue2] = GeoObj.nodata
    GeoObj.latitude[h5file['where/lat']['data'].value==in_fillvalue2] = GeoObj.nodata

    GeoObj.num_of_lines = GeoObj.latitude.shape[0]
    GeoObj.sec1970_start = h5file['how'].attrs['startepochs']
    GeoObj.sec1970_end =  h5file['how'].attrs['endepochs']
    tim1 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(GeoObj.sec1970_start))
    tim2 = time.strftime("%Y%m%d %H:%M", 
                         time.gmtime(GeoObj.sec1970_end))
    logger.info("Starttime: %s, end time: %s"%(tim1, tim2))
    logger.info("Min lon: %s, max lon: %d"%(
        np.min(GeoObj.longitude),
        np.max(GeoObj.longitude)))
    logger.info("Min lat: %d, max lat: %d"%(
        np.min(GeoObj.latitude),np.max(GeoObj.latitude)))
    return  GeoObj

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

def read_nwp_h5(filename, nwp_key):
    import h5py 
    h5file = h5py.File(filename, 'r')
    if nwp_key in h5file.keys():
        logger.info("Read NWP %s"%(nwp_key))
        value = h5file[nwp_key].value
        gain = h5file[nwp_key].attrs['gain']
        intercept = h5file[nwp_key].attrs['intercept']
        nodat = h5file[nwp_key].attrs['nodata']
        return  np.where(value != nodat,value * gain + intercept, value)
    else:
        logger.info("NO NWP %s File, Continue"%(nwp_key))
        return None

def read_etc_nc(ncFile, etc_key):
    if etc_key in ncFile.variables.keys():
        logger.info("Read %s"%(etc_key))
        return  ncFile.variables[etc_key][0,:,:]
    else:
        logger.info("NO %s File, Continue"%(etc_key))
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
            logger.error("Could not read %s File, Continue"%(thr_type))
        h5file.close()   
    else:
        logger.info("NO THR %s File, Continue"%(thr_type))
    return product


def readImagerData_h5(filename):
    h5file = h5py.File(filename, 'r')
    imager_data = NewImagerData()
    for var in h5file.keys():
        if 'image' in var:
            image = h5file[var]
            logger.info("reading channel %s", image.attrs['description'])
            one_channel = ImagerChannelData()                   
            one_channel.data = image['data'].value
            one_channel.des = image.attrs['description']
            #Currently unpacked arrays later in calipso.py:
            #TODO: move this herealso for h5!
            one_channel.gain = image['what'].attrs['gain']
            one_channel.intercept = image['what'].attrs['offset']
            imager_data.channel.append(one_channel) 
            imager_data.nodata = image['what'].attrs['nodata']
            imager_data.missingdata = image['what'].attrs['missingdata']
            imager_data.no_data = imager_data.nodata
            imager_data.missing_data = imager_data.missingdata
    h5file.close()
    return imager_data


def pps_read_all(pps_files, avhrr_file, cross):
    logger.info("Read AVHRR geolocation data")
    if '.nc' in avhrr_file:
        pps_nc = netCDF4.Dataset(avhrr_file, 'r', format='NETCDF4')
        avhrrGeoObj = read_pps_geoobj_nc(pps_nc)
    else:    
        #use mpop?
        avhrrGeoObj = read_pps_geoobj_h5(avhrr_file)    
    logger.info("Read AVHRR Sun -and Satellites Angles data")
    if '.nc' in pps_files.sunsatangles:
        pps_nc_ang = netCDF4.Dataset(pps_files.sunsatangles, 'r', format='NETCDF4')
        avhrrAngObj = read_pps_angobj_nc(pps_nc_ang)
    else:
        #use mpop?
        avhrrAngObj = read_pps_angobj_h5(pps_files.sunsatangles)
    logger.info("Read Imagerdata data")
    if '.nc' in avhrr_file:
        avhrrObj = readImagerData_nc(pps_nc)
    else:
        avhrrObj = readImagerData_h5(avhrr_file)

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
    if '.nc' in pps_files.cloudtype:
        ctype = read_cloudtype_nc(pps_files.cloudtype)
    else:
        ctype = read_cloudtype_h5(pps_files.cloudtype)
        #Removed not working read with mpop 2016-09-14
    logger.info("Read PPS CTTH")
    if '.nc' in pps_files.ctth:
        ctth = read_ctth_nc(pps_files.ctth)
    else:
        ctth = read_ctth_h5(pps_files.ctth)

    logger.info("Read PPS NWP data")
    nwp_dict={}
    if '.nc' in pps_files.nwp_tsur:
        pps_nc_nwp = netCDF4.Dataset(pps_files.nwp_tsur, 'r', format='NETCDF4')
        nwp_dict['surftemp'] = read_etc_nc(pps_nc_nwp, "tsur")
        nwp_dict['t500'] = read_etc_nc(pps_nc_nwp, "t500")
        nwp_dict['t700'] = read_etc_nc(pps_nc_nwp, "t700")
        nwp_dict['t850'] = read_etc_nc(pps_nc_nwp, "t850")
        nwp_dict['t950'] = read_etc_nc(pps_nc_nwp, "t950")
        nwp_dict['ttro'] = read_etc_nc(pps_nc_nwp, "ttro")
        nwp_dict['ciwv'] = read_etc_nc(pps_nc_nwp, "ciwv")
    else:   
         nwp_dict['surftemp'] = read_nwp_h5(pps_files.nwp_tsur,"tsur")
         nwp_dict['t500'] = read_nwp_h5(pps_files.nwp_t500, "t500")
         nwp_dict['t700'] = read_nwp_h5(pps_files.nwp_t700, "t700")
         nwp_dict['t850'] = read_nwp_h5(pps_files.nwp_t850, "t850")
         nwp_dict['t950'] = read_nwp_h5(pps_files.nwp_t950, "t950")
         nwp_dict['ttro'] = read_nwp_h5(pps_files.nwp_ttro, "ttro")
         nwp_dict['ciwv'] = read_nwp_h5(pps_files.nwp_ciwv, "ciwv")

    logger.info("Read PPS texture data")  
    if '.nc' in pps_files.text_t11:
        pps_nc_txt = netCDF4.Dataset(pps_files.text_t11, 'r', format='NETCDF4')
        for ttype in ['r06', 't11', 't37t12', 't37']:
            text_type = 'text_' + ttype
            nwp_dict[text_type] = read_etc_nc(pps_nc_txt, ttype)
    else:    
        for ttype in ['r06', 't11', 't37t12', 't37']:
            h5_obj_type = ttype +'_text'
            text_type = 'text_' + ttype
            nwp_dict[text_type] = read_thr_h5(getattr(pps_files,text_type), 
                                              h5_obj_type,text_type)
    logger.info("Read PPS threshold data") 
    if '.nc' in pps_files.thr_t11ts:
        pps_nc_thr = netCDF4.Dataset(pps_files.thr_t11ts, 'r', format='NETCDF4')
        for nc_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv', 
                            't11ts', 't11t37', 't37t12', 't11t12',
                            'r09', 'r06', 't85t11_inv', 't85t11']:
            thr_type = 'thr_' + nc_obj_type
            nwp_dict[thr_type] = read_etc_nc(pps_nc_thr,nc_obj_type)
    else:    
        for h5_obj_type in ['t11ts_inv', 't11t37_inv', 't37t12_inv', 't11t12_inv', 
                            't11ts', 't11t37', 't37t12', 't11t12',
                            'r09', 'r06', 't85t11_inv', 't85t11']:
            thr_type = 'thr_' + h5_obj_type
            nwp_dict[thr_type] = read_thr_h5(getattr(pps_files,thr_type), 
                                             h5_obj_type, thr_type)
    logger.info("Read PPS Emissivity data") 
    if '.nc' in pps_files.emis:
        pps_nc_thr = netCDF4.Dataset(pps_files.emis, 'r', format='NETCDF4')
        for emis_type in ['emis1',"emis6", 'emis8','emis9']:
            nwp_dict[emis_type] = read_etc_nc(pps_nc_thr, emis_type)
    else:
        for h5_obj_type in ['emis1',"emis6", 'emis8','emis9']:
            emis_type = h5_obj_type
            nwp_dict[emis_type] = read_thr_h5(getattr(pps_files,"emis"), 
                                              h5_obj_type, emis_type)
    nwp_obj = NWPObj(nwp_dict)
    logger.info("Read PPS NWP segment resolution data") 
    segment_data_object = read_segment_data(getattr(pps_files,'nwp_segments'))
    return avhrrAngObj, ctth, avhrrGeoObj, ctype, avhrrObj, nwp_obj, cppLwp, cppCph, segment_data_object 
