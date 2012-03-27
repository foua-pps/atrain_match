"""
Program cloudsat_calipso_avhrr_match.py

This program is used to process and output statistics for the inter-comparison of AVHRR PPS
results and CloudSat/CALIPSO observations. It may be run
repeatedly and supervised by program cloudsat_calipso_process_master.py.

This particular version of Adam's original CloudSat/CALIPSO matchup and analysis program has been complemented
with the following:

 * A method to calculate cloud emissivities for the uppermost CALIPSO cloud layer. With the
   use of parameters EMISS_FILTERING, EMISS_MIN_HEIGHT and EMISS_LIMIT the thinnest uppermost
   CALIPSO cloud layers can be analysed and the entire column can be disregarded if the
   emissivity falls below the EMISS_LIMIT value.
   Cloud emissivities Ec are calculated as follows:

                    Ec = (I-Iclear)/(B(Tc)-Iclear)
   where
       I = Measured radiance in AVHRR channel 4 (11 micron)
           To be calculated as the Planck radiance for the associated brightness temperature
       Iclear = Estimated radiance in cloud free situations
                To be calculate as the Planck radiance for the NWP-analysed surface temperature
                (i.e., neglecting further atmospheric contributions)
       B(Tc) = Planck radiance for the uppermost cloud layer using CALIPSO mid-layer temperatures

 * Adjusted scales between CloudSat and CALIPSO datasets. The previous assumption that both datasets
   had 1 km resolution resulted in that datasets went out of phase for distances longer than about
   1000 km. An empirical scale factor (CLOUDSAT_TRACK_RESOLUTION) of 1.076 is used to get the most
   optimal match.

 * AVHRR cloud top height datasets have been recalculated to heights above mean sea level using
   CloudSat and CALIPSO elevation data

 * The MODIS cloud flag has been added to the extracted CALIPSO dataset. This enables direct
   comparisons to the MODIS cloud mask! Consequently, corresponding MODIS Cloud Mask statistics
   are calculated and printed.

 * The Vertical Feature Mask parameter in the CALIPSO dataset has been used to subdivide results
   into three cloud groups: Low, Medium and High. This has enabled an evaluation of PPS Cloud Type
   results and a further sub-division of Cloud Top Height results

 * The National Snow and Ice Data Center (NSIDC) ice and snow mapping results have been added
   to the extracted Calipso parameters. Together with the IGBP land use classification it is then
   possible to isolate the study to focus on one of the following categories:

       ICE_COVER_SEA, ICE_FREE_SEA, SNOW_COVER_LAND, SNOW_FREE_LAND or COASTAL_ZONE

RUNNING INSTRUCTIONS
--------------------

The program is capable of running in a wide range of modes (according to description above). These
various modes are selected by enabling (disabling) the following parameters:

PLOT_OPTION, EMISS_FILTERING, ICE_COVER_SEA, ICE_FREE_SEA, SNOW_COVER_LAND, SNOW_FREE_LAND, COASTAL_ZONE

However, notice that only one mode can be chosen for each run. The only exception is PLOT_OPTION
(i.e., the generation of a PNG plot) which can be combined with EMISS_FILTERING. This also means
that processing of individual surface categories only generates statistics and not any plots.

Input data has to be supplied at directories defined by MAIN_DIR and SUB_DIR parameters below.
 
Every exection of the program prints out statistics for PPS Cloud Mask, Cloud Type and Cloud Top Height
directly on the screen.

Don't forget to set all parameters needed for standard ACPG/AHAMAP execution since the matchup
software uses parts of the ACPG/AHAMAP software.

Dependencies: For a successful run of the program the following supporting python modules must be
              available in the default run directory:

              cloudsat.py
              calipso.py
              calipso_avhrr_matchup.py
              cloudsat_avhrr_matchup.py
              radiance_tb_tables_kgtest.py

              For full consistency make sure that MAIN_DIR and SUB_DIR parameters are the same also
              in modules cloudsat.py and calipso.py.

Output data files: Main results are generally written in the directory MAIN_DIR/SUB_DIR but plotting
                   results are stored at ./Plot and temporary results at ./Data directories. Thus,
                   make sure that these directories exist as subdirectories at the default run directory.
                  This is now made automatic /Erik

Finally, notice that the matching of the PPS, CloudSat and CALIPSO datasets have been calculated using the
fix area arctic_super_5010 defined over the Arctic region in Lambert Azimuthal Equal Area projection. Thus,
for matching data to other regions please modify modules calipso.py and cloudsat.py and replace area
arctic_super_5010 with the desired area.
This is now made below /Erik

/KG March 2010

"""
import os
import sys
from numpy import NaN
#import inspect
from radiance_tb_tables_kgtest import * #Just use the brightness temperature to radiance conversion/KG @UnusedWildImport

from pps_basic_configure import *
from pps_error_messages import * #@UnusedWildImport

#import config
from common import attach_subdir_from_config, MatchupError

from cloudsat_calipso_avhrr_statistics import *
from trajectory_plot import * #@UnusedWildImport
from cloudsat_calipso_avhrr_prepare import *

from cloudsat import reshapeCloudsat, match_cloudsat_avhrr,\
    writeCloudsatAvhrrMatchObj, readCloudsatAvhrrMatchObj
from calipso import reshapeCalipso, match_calipso_avhrr,\
    writeCaliopAvhrrMatchObj, readCaliopAvhrrMatchObj, add1kmTo5km
import inspect
import numpy
from cloudsat_calipso_avhrr_plot import drawCalClsatAvhrrPlotTimeDiff,\
    drawCalClsatGEOPROFAvhrrPlot, drawCalClsatAvhrrPlotSATZ,\
    drawCalClsatCWCAvhrrPlot
import pdb #@UnusedImport
from config import CALIPSO_CLOUD_FRACTION
#from cloudsat_avhrr_matchup5km import *
#from calipso_avhrr_matchup5km import *
test = 0
# -----------------------------------------------------

def find_avhrr_file(cross):
    """
    Find the *satellite* avhrr file closest to *datetime*.
    """
    if cross.satellite1.lower() in ['cloudsat', 'calipso']:
        cross_satellite = cross.satellite2.lower()
        cross_time = cross.time2
    else:
        cross_satellite = cross.satellite1.lower()
        cross_time = cross.time1
    
    from file_finders import PpsFileFinder #@UnresolvedImport
    avhrr_finder = PpsFileFinder(config.PPS_DATA_DIR, 'avhrr.h5')
    attach_subdir_from_config(avhrr_finder)
    
    # Set time window
    if cross.time_window is not None:
        avhrr_finder.set_time_window(-(config.SAT_ORBIT_DURATION + cross.time_window), 0)
    else:
        avhrr_finder.set_time_window(-(config.SAT_ORBIT_DURATION + config.sec_timeThr), 0)
    try:
        try:            
            if config.DEBUG is True:
                write_log('DEBUG', "Looking for avhrr file corresponding to %s: %s" % \
                          (str(cross), avhrr_finder.pattern(cross_time, cross_satellite))) #@UndefinedVariable
        except:
            pass

        avhrr_file = avhrr_finder.find(cross_time, cross_satellite)[0]
    except IndexError:
        def dpp_subdir(time, *args, **kwargs):
            return os.path.join('%d' % time.year, '%02d' % time.month,
                                '%d%02d%02d' %(time.year, time.month, time.day))
        avhrr_finder.subdir = dpp_subdir
        try:
            avhrr_file = avhrr_finder.find(cross_time, cross_satellite)[0]
        except IndexError:
            raise MatchupError("No avhrr file found for %s." % cross)

    return avhrr_file


def find_files_from_avhrr(avhrr_file):
    """
    Find all files needed to process matchup from source data files.
    """
    import file_finders #@UnresolvedImport
    
    pps_finder = file_finders.PpsFileFinder(config.PPS_DATA_DIR, time_window=5*60)
    attach_subdir_from_config(pps_finder)
    parsed = pps_finder.parse(avhrr_file)
    satname = parsed['satellite']
    datetime = parsed['datetime']
    
    cloudsat_finder = file_finders.CloudsatFileFinder(config.CLOUDSAT_DIR,
                                                      config.RESOLUTION,
                                                      config.CLOUDSAT_TYPE)
    attach_subdir_from_config(cloudsat_finder)
    cloudsat_finder.set_time_window(config.SAT_ORBIT_DURATION + config.sec_timeThr)
    cloudsat_files = sorted(cloudsat_finder.find(datetime))
    if len(cloudsat_files) == 0:
        #TODO: Ordna det pa ett snyggare satt =)
        pass
#        raise MatchupError("No cloudsat-%s files found corresponding to %s." % \
#                           (config.CLOUDSAT_TYPE, avhrr_file))
    
    calipso_finder = file_finders.CalipsoFileFinder(config.CALIPSO_DIR,
                                                    config.RESOLUTION)
    attach_subdir_from_config(calipso_finder)
    calipso_finder.set_time_window(config.SAT_ORBIT_DURATION + config.sec_timeThr)
    calipso_files = sorted(calipso_finder.find(datetime))
    if len(calipso_files) == 0:
        raise MatchupError("No calipso files found corresponding to %s." % avhrr_file)
    # nedan kollar om filerna man hittar ar i hdf format eller h5 format. I fall
    # de inte finns i h5 men i hdf sa converteras de till h5.
    calipso_h5_med_filendelse = []
    calipso_hdf_utan_filendelse = []
    calipso_h5_utan_filendelse = []
    for caFiles in calipso_files:
        if caFiles.split('.')[-1] == 'hdf':
            calipso_hdf_utan_filendelse.append('.'.join(caFiles.split('.')[:-1]))
        elif caFiles.split('.')[-1] == 'h5':
            calipso_h5_med_filendelse.append(caFiles)
            calipso_h5_utan_filendelse.append('.'.join(caFiles.split('.')[:-1]))
    
    for hdfU in calipso_hdf_utan_filendelse:
        if hdfU not in calipso_h5_utan_filendelse:
            cmdConv = '/home/pps/opt/H4H5/2.1.1/bin/h4toh5' + ' ' + hdfU + '.hdf'
            print('.hdf file exist but not .h5. Therefore convert')
            os.system(cmdConv)            
            calipso_h5_med_filendelse.append(hdfU + '.h5')
    calipso_files = sorted(calipso_h5_med_filendelse)
    try:
        cloudtype_file = pps_finder.find(datetime, satname, ending='cloudtype.h5')[0]
    except IndexError:
        def dpp_subdir(time, *args, **kwargs):
            return os.path.join('%d' % time.year, '%02d' % time.month,
                                '%d%02d%02d' %(time.year, time.month, time.day))
        pps_finder.subdir = dpp_subdir
        try:
            cloudtype_file = pps_finder.find(datetime, satname, ending='cloudtype.h5')[0]
        except IndexError:
            raise MatchupError("No cloudtype file found corresponding to %s." % avhrr_file)
    
    try:
        ctth_file = pps_finder.find(datetime, satname,
                                    ending='%s.h5' % config.CTTH_FILE)[0]
    except IndexError:
        def dpp_subdir(time, *args, **kwargs):
            return os.path.join('%d' % time.year, '%02d' % time.month,
                                '%d%02d%02d' %(time.year, time.month, time.day))
        pps_finder.subdir = dpp_subdir
        try:
            ctth_file = pps_finder.find(datetime, satname,
                                    ending='%s.h5' % config.CTTH_FILE)[0] 
        except IndexError:
            raise MatchupError("No %s file found corresponding to %s." % \
                               (config.CTTH_FILE, avhrr_file))
    
    try:
        nwp_tsur_file = pps_finder.find(datetime, satname, ending='nwp_tsur.h5')[0]
    except IndexError:
        def dpp_subdir(time, *args, **kwargs):
            return os.path.join('%d' % time.year, '%02d' % time.month,
                                '%d%02d%02d' %(time.year, time.month, time.day))
        pps_finder.subdir = dpp_subdir
        try:
            nwp_tsur_file = pps_finder.find(datetime, satname, ending='nwp_tsur.h5')[0]
        except IndexError:
            #TODO: Fixa pa ett snyggare satt
            nwp_tsur_file = []
#            raise MatchupError("No nwp_tsur file found corresponding to %s." % avhrr_file)
    
    try:
        sunsatangles_file = pps_finder.find(datetime, satname, ending='sunsatangles.h5')[0]
    except IndexError:
        def dpp_subdir(time, *args, **kwargs):
            return os.path.join('%d' % time.year, '%02d' % time.month,
                                '%d%02d%02d' %(time.year, time.month, time.day))
        pps_finder.subdir = dpp_subdir
        try:
            sunsatangles_file = pps_finder.find(datetime, satname, ending='sunsatangles.h5')[0]
        except IndexError:
            raise MatchupError("No sunsatangles file found corresponding to %s." % avhrr_file)

    return (cloudsat_files, calipso_files, cloudtype_file, ctth_file, 
            nwp_tsur_file, sunsatangles_file)


def get_cloudsat_matchups(cloudsat_files, cloudtype_file, avhrrGeoObj, avhrrObj, ctype, ctth, surft, avhrrAngObj):
    """
    Read Cloudsat data and match with the given PPS data.
    """
    if config.CLOUDSAT_TYPE == 'CWC-RVOD':
        if config.RESOLUTION == 1:      
            if test == 1:
                from cloudsat_cwc import match_cloudsatCwc_avhrr
                from cloudsat_cwc import reshapeCloudsat1kmCwc
                reshape_fun = reshapeCloudsat1kmCwc
                match_fun = match_cloudsatCwc_avhrr
            else:
                reshape_fun = reshapeCloudsat
                match_fun = match_cloudsat_avhrr 
        elif config.RESOLUTION == 5:
            from cloudsat5km_cwc import match_cloudsat5kmCwc_avhrr5km
            from cloudsat5km_cwc import reshapeCloudsat5kmCwc
            reshape_fun = reshapeCloudsat5kmCwc
            match_fun = match_cloudsat5kmCwc_avhrr5km
    else:
        reshape_fun = reshapeCloudsat
        match_fun = match_cloudsat_avhrr
    cloudsat = reshape_fun(cloudsat_files, avhrrGeoObj, cloudtype_file)
    cl_matchup, cl_min_diff, cl_max_diff = match_fun(cloudtype_file, cloudsat,
                                                     avhrrGeoObj, avhrrObj, ctype,
                                                     ctth, surft, avhrrAngObj)
    
    return cl_matchup, (cl_min_diff, cl_max_diff)


def get_calipso_matchups(calipso_files, cloudtype_file, avhrrGeoObj, avhrrObj, ctype, ctth, surft, avhrrAngObj, cafiles1km=None):
    """
    Read Calipso data and match with the given PPS data.
    """
    if cafiles1km != None:
        calipso1km = reshapeCalipso(cafiles1km, avhrrGeoObj, cloudtype_file, False, 1)[0]
        calipso5km, startBreak, endBreak = reshapeCalipso(calipso_files,avhrrGeoObj, cloudtype_file, False)
        calipso = add1kmTo5km(calipso1km, calipso5km, startBreak, endBreak)
        calipso1km = None
        calipso5km = None
    else:
        calipso = reshapeCalipso(calipso_files,avhrrGeoObj, cloudtype_file, True)[0]

    ca_matchup, ca_min_diff, ca_max_diff = match_calipso_avhrr(cloudtype_file, calipso,
                                                     avhrrGeoObj, avhrrObj, ctype,
                                                     ctth, surft, avhrrAngObj)
    
    return ca_matchup, (ca_min_diff, ca_max_diff)


def get_matchups_from_data(cross):
    """
    Retrieve Cloudsat- and Calipso-AVHRR matchups from Cloudsat, Calipso, and
    PPS files.
    """
    import pps_io #@UnresolvedImport
    import epshdf #@UnresolvedImport
    import os #@Reimport
    
    avhrr_file = find_avhrr_file(cross)
    
    cloudsat_files, calipso_files, cloudtype_file, ctth_file, nwp_tsur_file, \
        sunsatangles_file = find_files_from_avhrr(avhrr_file)

    write_log("INFO","Read AVHRR geolocation data") #@UndefinedVariable
    avhrrGeoObj = pps_io.readAvhrrGeoData(avhrr_file)
    
    write_log("INFO","Read AVHRR Sun -and Satellites Angles data") #@UndefinedVariable
    avhrrAngObj = pps_io.readSunSatAngles(sunsatangles_file) #, withAbsoluteAzimuthAngles=True)
    
    write_log("INFO","Read AVHRR data") #@UndefinedVariable
    avhrrObj = pps_io.readAvhrrData(avhrr_file)

    write_log("INFO","Read PPS Cloud Type") #@UndefinedVariable
    ctype = epshdf.read_cloudtype(cloudtype_file,1,1,0)
    try:
        ctth = epshdf.read_cloudtop(ctth_file,1,1,1,0,1)
    except:
        ctth = None

    if isinstance(nwp_tsur_file, str) == True:
        write_log("INFO","Read NWP surface temperature") #@UndefinedVariable
        nwpinst = epshdf.read_nwpdata(nwp_tsur_file)
        surft = nwpinst.gain*nwpinst.data.astype('d') + nwpinst.intercept
    else:
        write_log("INFO","NO NWP surface temperature File, Continue") #@UndefinedVariable
        surft = None
    
    if isinstance(cloudsat_files, str) == True:
        write_log("INFO","Read CLOUDSAT %s data" % config.CLOUDSAT_TYPE) #@UndefinedVariable
        cl_matchup, cl_time_diff = get_cloudsat_matchups(cloudsat_files, cloudtype_file,
                                                         avhrrGeoObj, avhrrObj, ctype,
                                                         ctth, surft, avhrrAngObj)
    else:
        write_log("INFO","NO CLOUDSAT File, Continue") #@UndefinedVariable
    
    write_log("INFO","Read CALIPSO data") #@UndefinedVariable
    if config.RESOLUTION == 5:
        import glob
        calipso1km = []
        for file5km in calipso_files:
            file1km = file5km.replace('/5km/', '/1km/').\
                                        replace('05kmCLay', '01kmCLay').\
                                        replace('-Prov-V3-01.', '*')
            calipso1km.append(glob.glob(file1km)[0])
        if len(calipso_files) != len(calipso1km):
            pdb.set_trace()
    ca_matchup, ca_time_diff = get_calipso_matchups(calipso_files, cloudtype_file,
                                                    avhrrGeoObj, avhrrObj, ctype,
                                                    ctth, surft, avhrrAngObj, calipso1km)
    pdb.set_trace()
    
    # Get satellite name, time, and orbit number from avhrr_file
    from file_finders import PpsFileFinder #@UnresolvedImport
    pps_finder = PpsFileFinder()
    parsed = pps_finder.parse(avhrr_file)
    satellite = parsed['satellite']
    datetime = parsed['datetime']
    basename = '_'.join(os.path.basename(avhrr_file).split('_')[:4])
    
    # Build base file name
    from file_finders import CloudsatCalipsoAvhrrMatchFileFinder #@UnresolvedImport
    match_finder = CloudsatCalipsoAvhrrMatchFileFinder(resolution=config.RESOLUTION,
                                                       region=config.AREA)
    attach_subdir_from_config(match_finder)
    rematched_file_base = os.path.join(config.RESHAPE_DIR,
                                       match_finder.subdir(datetime, satname=satellite),
                                       "%dkm_%s_atrain_datatype_avhrr_match.h5" % \
                                       (config.RESOLUTION, basename))
    
    # Create directories if they don't exist yet
    if not os.path.exists(os.path.dirname(rematched_file_base)):
        os.makedirs(os.path.dirname(rematched_file_base))
    
    # Write cloudsat matchup
    try:
        cl_match_file = rematched_file_base.replace('atrain_datatype', 'cloudsat-%s' % config.CLOUDSAT_TYPE)
        writeCloudsatAvhrrMatchObj(cl_match_file, cl_matchup)
    except NameError:
        cl_matchup = None
        cl_time_diff = None
        print('CloudSat is not defined. No CloudSat Match File created')
    # Write calipso matchup
    ca_match_file = rematched_file_base.replace('atrain_datatype', 'caliop')
    writeCaliopAvhrrMatchObj(ca_match_file,ca_matchup)

    return {'cloudsat': cl_matchup, 'cloudsat_time_diff': cl_time_diff,
            'calipso': ca_matchup, 'calipso_time_diff': ca_time_diff,
            'basename': basename}


def get_matchups(cross, reprocess=False):
    """
    Retrieve Cloudsat- and Calipso-AVHRR matchups. If *reprocess* is False, and if
    matchup files exist, get matchups directly from the processed files.
    Otherwise process Cloudsat, Calipso, and PPS files first.
    """
    caObj = None
    clObj = None
    
    try:
        satellite = cross.satellite1.lower()
        datetime = cross.time1
    except AttributeError:
        raise ValueError('cross is not a valid SNO cross. (cross: %s)' % cross)
    
    if satellite in ['calipso', 'cloudsat']:
        satellite = cross.satellite2.lower()
        datetime = cross.time2
        
    if reprocess is False:
        from file_finders import CloudsatCalipsoAvhrrMatchFileFinder #@UnresolvedImport

        match_finder = CloudsatCalipsoAvhrrMatchFileFinder(config.RESHAPE_DIR,
                                                           config.RESOLUTION,
                                                           region=config.AREA)
        # Set time window
        if cross.time_window is not None:
            match_finder.set_time_window(-(config.SAT_ORBIT_DURATION + cross.time_window), 0)
        else:
            match_finder.set_time_window(-(config.SAT_ORBIT_DURATION + config.sec_timeThr), 0)
        
        attach_subdir_from_config(match_finder)
        try:
            cl_match_file = match_finder.find(datetime, satellite, atrain_datatype='cloudsat-%s' % config.CLOUDSAT_TYPE)[0]
            clObj = readCloudsatAvhrrMatchObj(cl_match_file) 
            basename = '_'.join(os.path.basename(cl_match_file).split('_')[1:5])
            write_log('INFO', "CloudSat Matchups read from previously processed data.") #@UndefinedVariable
        except IndexError:
            write_log('INFO', "No processed CloudSat match files found. Generating from source data.") #@UndefinedVariable
        try:
            ca_match_file = match_finder.find(datetime, satellite, atrain_datatype='caliop')[0]
            caObj = readCaliopAvhrrMatchObj(ca_match_file)
            basename = '_'.join(os.path.basename(ca_match_file).split('_')[1:5])    
            write_log('INFO', "CALIPSO Matchups read from previously processed data.") #@UndefinedVariable
        except IndexError:
            write_log('INFO', "No processed CALIPSO match files found. Generating from source data.") #@UndefinedVariable
    #TODO: Fix a better solution for below so it can handle missing cloudsat better.
    if None in [caObj] and None in [clObj]:
        return get_matchups_from_data(cross)
    elif None in [caObj]:
        return get_matchups_from_data(cross)
    elif None in [clObj]:
        return {'calipso': caObj, 'cloudsat': clObj, 'basename': basename}
    else:
        return {'calipso': caObj, 'cloudsat': clObj, 'basename': basename}


def run(cross, process_mode_dnt, reprocess=False):
    """
    The main work horse.
    
    """
    
    write_log('INFO', "Case: %s" % str(cross)) #@UndefinedVariable
    write_log('INFO', "Process mode: %s" % process_mode_dnt) #@UndefinedVariable
    cloudsat_type = config.CLOUDSAT_TYPE
    # split process_mode_dnt into two parts. One with process_mode and one dnt_flag
    mode_dnt = process_mode_dnt.split('_')
    if len(mode_dnt) == 1:
        process_mode = process_mode_dnt
        dnt_flag = None
    elif mode_dnt[-1] in ['DAY', 'NIGHT', 'TWILIGHT']:
        process_mode = '_'.join(mode_dnt[0:-1])
        dnt_flag = mode_dnt[-1]
    else:
        process_mode = process_mode_dnt
        dnt_flag = None

    sno_satname = cross.satellite1.lower()
    if sno_satname in ['calipso', 'cloudsat']:
        sno_satname = cross.satellite2.lower()
        
    if sno_satname == "noaa17":
        noaa_number=17
    elif sno_satname == "noaa18":
        noaa_number=18
    elif sno_satname == "metop02":
        noaa_number=2 # Poor man's solution!
    elif sno_satname == 'noaa19':
        noaa_number = 19
    else:
        raise NotImplementedError("Support for satellite %s is not yet implemented." % sno_satname)

    # Now fetch all the datasets for the section of the AREA where all
    # three datasets match. Also get maximum and minimum time differences to AVHRR (in seconds)
    matchup_results = get_matchups(cross, reprocess)
    caObj = matchup_results['calipso']
    clsatObj = matchup_results['cloudsat']
    clsat_min_diff, clsat_max_diff = matchup_results.get('cl_time_diff', (NaN, NaN))
    ca_min_diff, ca_max_diff = matchup_results.get('ca_time_diff', (NaN, NaN))
    
    basename = matchup_results['basename']
    base_sat = basename.split('_')[0]
    base_year = basename.split('_')[1][:4]
    base_month = basename.split('_')[1][4:6]
    
        #clsatObj,clsat_min_diff,clsat_max_diff = getCloudsat5kmAvhrr5kmMatch(avhrrfile,cloudsatfile,ctypefile,ctthfile,surftfile,sunanglefile,cloudsat_type)
        #caObj,ca_min_diff,ca_max_diff = getCaliop5kmAvhrr5kmMatch(avhrrfile,calipsofile,ctypefile,ctthfile,surftfile,sunanglefile)        
        #clsatObj,clsat_min_diff,clsat_max_diff = getCloudsatAvhrrMatch(avhrrfile,cloudsatfile,ctypefile,ctthfile,surftfile,sunanglefile,cloudsat_type)
        #caObj,ca_min_diff,ca_max_diff = getCaliopAvhrrMatch(avhrrfile,calipsofile,ctypefile,ctthfile,surftfile,sunanglefile)  
        #caObj,ca_min_diff,ca_max_diff = getCaliop5kmAvhrr5kmMatch(avhrrfile,calipsofile1,calipsofile2,calipsofile3,ctypefile,ctthfile,surftfile,int(Resolution[-1]))
    if config.RESOLUTION not in [1, 5]:
        write_log("INFO","Define resolution") #@UndefinedVariable
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit(-9) 
    ## Cloudsat ##
    if clsatObj != None:
        if cloudsat_type == 'GEOPROF':
            clsatObj = CloudsatAvhrrSatz(clsatObj)
            cllon = clsatObj.cloudsat.longitude.copy()
            cllat = clsatObj.cloudsat.latitude.copy()
        elif cloudsat_type == 'CWC-RVOD':
            if config.RESOLUTION == 1:
                cllon = clsatObj.cloudsatcwc.longitude.copy()
                cllat = clsatObj.cloudsatcwc.latitude.copy()
            elif config.RESOLUTION == 5:
                cllon = clsatObj.cloudsat5kmcwc.longitude.copy()
                cllat = clsatObj.cloudsat5kmcwc.latitude.copy()
        # Issue a warning if startpoint or endpoint latitude of CloudSat and CALIPSO differ by more than 0.1 degrees
        # Furthermore, if startpoint differs it means that CALIPSO data is not available for the first part of the matchup
        # cross section. This means that we must find the first corresponding CloudSat point to this CALIPSO start point
        # before we can do the plotting (statistics calculations are not affected). Consequently, set the CALIPSO_DISPLACED
        # flag and find correct startpoint just before starting the plotting of CALIPSO data!
        CALIPSO_DISPLACED = 0
        latdiff = abs(clsatObj.cloudsat.latitude[0] - caObj.calipso.latitude[0])
        print "latdiff: ", latdiff
        if latdiff > 0.1:
            write_log('INFO', "CloudSat/CALIPSO startpoint differ by %f degrees." % latdiff) #@UndefinedVariable
            write_log('INFO', "Cloudsat start lon, lat: %f, %f" % \
                      (clsatObj.cloudsat.longitude[0], clsatObj.cloudsat.latitude[0])) #@UndefinedVariable
            write_log('INFO', "CALIPSO start lon, lat: %f, %f" % \
                      (caObj.calipso.longitude[0], caObj.calipso.latitude[0])) #@UndefinedVariable
            CALIPSO_DISPLACED = 1
            for j in range(clsatObj.cloudsat.latitude.shape[0]):
                if (abs(clsatObj.cloudsat.latitude[j] - caObj.calipso.latitude[0]) < 0.05) and\
                       (abs(clsatObj.cloudsat.longitude[j] - caObj.calipso.longitude[0]) < 0.1):
                    calipso_displacement=int(j*config.CLOUDSAT_TRACK_RESOLUTION)
                    write_log('INFO', "CALIPSO_DISPLACEMENT: %d" % calipso_displacement) #@UndefinedVariable
                    break
        # First make sure that PPS cloud top heights are converted to height above sea level
        # just as CloudSat and CALIPSO heights are defined. Use corresponding DEM data.            
        elevation = numpy.where(numpy.less_equal(clsatObj.cloudsat.elevation,0),\
                            -9,clsatObj.cloudsat.elevation)			# If clsatObj.cloudsat.elevation is <= 0 elevation(i,j)=-9, else the value = clsatObj.cloudsat.elevation(i,j)
        data_ok = numpy.ones(clsatObj.cloudsat.elevation.shape,'b')
        print "Length of CLOUDSAT array: ", len(data_ok)
        lat_ok = numpy.repeat(clsatObj.cloudsat.latitude[::],data_ok)
        avhrr_ctth_csat_ok = numpy.repeat(clsatObj.avhrr.ctth_height[::],data_ok)
        avhrr_ctth_csat_ok = numpy.where(numpy.greater(avhrr_ctth_csat_ok,0.0),avhrr_ctth_csat_ok[::]+elevation*1.0,avhrr_ctth_csat_ok)
        if len(data_ok) == 0:
            print "Processing stopped: Zero lenght of matching arrays!"
            print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
            sys.exit()
    else:
        data_ok = None
        avhrr_ctth_csat_ok = None
    ## Calipso ##        
    caObj = CalipsoAvhrrSatz(caObj)
    calon = caObj.calipso.longitude.copy()
    calat = caObj.calipso.latitude.copy()
    avhrlon = caObj.avhrr.longitude.copy()
    avhrlat = caObj.avhrr.latitude.copy()

    # First make sure that PPS cloud top heights are converted to height above sea level
    # just as CloudSat and CALIPSO heights are defined. Use corresponding DEM data.
    cal_elevation = numpy.where(numpy.less_equal(caObj.calipso.elevation,0),
                                -9,caObj.calipso.elevation)
    cal_data_ok = numpy.ones(caObj.calipso.elevation.shape,'b')
    print "Length of CALIOP array: ", len(cal_data_ok)
    avhrr_ctth_cal_ok = numpy.repeat(caObj.avhrr.ctth_height[::],cal_data_ok)
    avhrr_ctth_cal_ok = numpy.where(numpy.greater(avhrr_ctth_cal_ok,0.0),avhrr_ctth_cal_ok[::]+cal_elevation,avhrr_ctth_cal_ok)                    

    if (len(cal_data_ok) == 0):
        print "Processing stopped: Zero lenght of matching arrays!"
        print("Program cloudsat_calipso_avhrr_match.py at line %i" %(inspect.currentframe().f_lineno+1))
        sys.exit()
        
    # If everything is OK, now create filename for statistics output file and open it for writing.
    # Notice that more than one file
    # (but maximum 2) can be created for one particular noaa orbit.
    resultpath = "%s/%s/%ikm/%s/%s/%s/%s" % (config.RESULT_DIR, base_sat,int(config.RESOLUTION), base_year, base_month, config.AREA, process_mode_dnt)
    if process_mode == 'OPTICAL_DEPTH':
        resultpath = "%s/%s/%ikm/%s/%s/%s/%s-%.1f" % (config.RESULT_DIR, base_sat,int(config.RESOLUTION), base_year, base_month,config.AREA, process_mode_dnt, config.MIN_OPTICAL_DEPTH)        

    if not os.path.exists(resultpath):
        os.makedirs(resultpath)
    statfilename = "%s/%ikm_%s_cloudsat_calipso_avhrr_stat.dat" % (resultpath,int(config.RESOLUTION),basename)
    if CALIPSO_CLOUD_FRACTION == True:
        statfilename = "%s/CCF_%ikm_%s_cloudsat_calipso_avhrr_stat.dat" % (resultpath,int(config.RESOLUTION),basename)        
    statfile = open(statfilename,"w")
    if process_mode == "BASIC":
        if clsatObj != None:
            statfile.write("CloudSat min and max time diff: %f %f \n" %(clsat_min_diff,clsat_max_diff))
        else:
            statfile.write('No CloudSat \n')
        statfile.write("CALIPSO min and max time diff: %f %f \n" %(ca_min_diff,ca_max_diff))
    else:
        if clsatObj != None:
            statfile.write("CloudSat min and max time diff: See results for BASIC! \n")
        else:
            statfile.write('No CloudSat \n')
        statfile.write("CALIPSO min and max time diff: See results for BASIC! \n")
    if clsatObj != None:
        statfile.write("Start-Stop-Length Cloudsat: %f %f %f %f %s \n" %(clsatObj.cloudsat.latitude[0],clsatObj.cloudsat.longitude[0],clsatObj.cloudsat.latitude[len(clsatObj.cloudsat.latitude)-1],clsatObj.cloudsat.longitude[len(clsatObj.cloudsat.latitude)-1],len(data_ok)))
    else:
        statfile.write('No CloudSat \n')
    statfile.write("Start-Stop-Length CALIPSO: %f %f %f %f %s \n" %(caObj.calipso.latitude[0],caObj.calipso.longitude[0],caObj.calipso.latitude[len(caObj.calipso.latitude)-1],caObj.calipso.longitude[len(caObj.calipso.latitude)-1],len(cal_data_ok)))
    
    if CALIPSO_CLOUD_FRACTION == True:
        (new_cloud_top, new_cloud_base, new_optical_depth, new_cloud_fraction, new_fcf, new_ssccf) = \
            CalipsoCloudFraction(caObj.calipso.cloud_top_profile, \
                                 caObj.calipso.cloud_base_profile, \
                                 caObj.calipso.optical_depth, \
                                 caObj.calipso.cloud_fraction, \
                                 caObj.calipso.feature_classification_flags, \
                                 caObj.calipso.single_shot_cloud_cleared_fraction)
        caObj.calipso.cloud_top_profile = new_cloud_top
        caObj.calipso.cloud_base_profile = new_cloud_base
        caObj.calipso.optical_depth = new_optical_depth
        caObj.calipso.cloud_fraction = new_cloud_fraction
        caObj.calipso.feature_classification_flags = new_fcf
        caObj.calipso.single_shot_cloud_cleared_fraction = new_ssccf
    "If mode = OPTICAL_DEPTH -> Change cloud -top and -base profile"
    if process_mode == 'OPTICAL_DEPTH':
        (new_cloud_top, new_cloud_base, new_cloud_fraction, new_fcf) = \
            CloudsatCloudOpticalDepth(caObj.calipso.cloud_top_profile, caObj.calipso.cloud_base_profile, \
            caObj.calipso.optical_depth, caObj.calipso.cloud_fraction, caObj.calipso.feature_classification_flags)
        caObj.calipso.cloud_top_profile = new_cloud_top
        caObj.calipso.cloud_base_profile = new_cloud_base
        caObj.calipso.cloud_fraction = new_cloud_fraction
        caObj.calipso.feature_classification_flags = new_fcf
    # Extract CALIOP Vertical Feature Classification (bits 10-12) from 16 bit representation
    # for topmost cloud layer
    cal_vert_feature = numpy.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9
    feature_array = 4*numpy.bitwise_and(numpy.right_shift(caObj.calipso.feature_classification_flags[0,::],11),1) + 2*numpy.bitwise_and(numpy.right_shift(caObj.calipso.feature_classification_flags[0,::],10),1) + numpy.bitwise_and(numpy.right_shift(caObj.calipso.feature_classification_flags[0,::],9),1)

    cal_vert_feature = numpy.where(numpy.not_equal(caObj.calipso.feature_classification_flags[0,::],1),feature_array[::],cal_vert_feature[::])   
    # Prepare for plotting, cloud emissivity and statistics calculations
    
    midlayer_temp_kelvin = caObj.calipso.cloud_mid_temperature + 273.15
    caliop_toplay_thickness = numpy.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9

    caliop_max_height_midlaytemp = numpy.where(numpy.greater(caObj.calipso.cloud_top_profile[0,::],-9),midlayer_temp_kelvin[0,::],-9)
    thickness = (caObj.calipso.cloud_top_profile[0,::]-caObj.calipso.cloud_base_profile[0,::])*1000.
    caliop_toplay_thickness = numpy.where(numpy.greater(caObj.calipso.cloud_top_profile[0,::],-9),thickness,caliop_toplay_thickness)

    
    caliop_height = []
    caliop_base = []
    caliop_max_height = numpy.ones(caObj.calipso.cloud_top_profile[0,::].shape)*-9
    
    for i in range(10):
        hh = numpy.where(numpy.greater(caObj.calipso.cloud_top_profile[i,::],-9),
                            caObj.calipso.cloud_top_profile[i,::] * 1000.,-9)
                                        
        caliop_max_height = numpy.maximum(caliop_max_height,
                                            caObj.calipso.cloud_top_profile[i,::] * 1000.)
        # This is actually unnecessary - we know that layer 1 is always the highest layer!!
        # However, arrays caliop_height and caliop_base are needed later for plotting/ KG

        caliop_height.append(hh)
        bb = numpy.where(numpy.greater(caObj.calipso.cloud_base_profile[i,::],-9),
                            caObj.calipso.cloud_base_profile[i,::] * 1000.,-9)
        caliop_base.append(bb)
        thickness = hh - bb
        
    x = numpy.repeat(caObj.calipso.number_of_layers_found.ravel(),
                        numpy.greater(caObj.calipso.number_of_layers_found.ravel(),0))
    #print "Number of points with more than 0 layers: ",x.shape[0]
    cal_data_ok = numpy.greater(caliop_max_height,-9.)
    cal_lat_ok = numpy.repeat(caObj.calipso.latitude[::],cal_data_ok)
    cal_avhrr_ctth_ok = numpy.repeat(avhrr_ctth_cal_ok[::],cal_data_ok)
    cal_topmidlay_temp_ok = numpy.repeat(caliop_max_height_midlaytemp[::],cal_data_ok)
    cal_toplay_thickness_ok = numpy.repeat(caliop_toplay_thickness[::],cal_data_ok)
    cal_maxheight_ok = numpy.repeat(caliop_max_height[::],cal_data_ok)
    cal_bt11temp_ok = numpy.repeat(caObj.avhrr.bt11micron[::],cal_data_ok)

    if caObj.avhrr.surftemp != None:
        cal_surftemp_ok = numpy.repeat(caObj.avhrr.surftemp[::],cal_data_ok)

    if clsatObj != None:
        # Transfer CloudSat MODIS cloud flag to CALIPSO representation
        cal_MODIS_cflag = numpy.zeros(len(cal_data_ok),'b')
        for i in range(len(cal_data_ok)):
            if not CALIPSO_DISPLACED: 
                cloudsat_index = int(i/config.CLOUDSAT_TRACK_RESOLUTION + 0.5)
                #print len(clsatObj.cloudsat.MODIS_cloud_flag),cloudsat_index
                if len(clsatObj.cloudsat.MODIS_cloud_flag) > (cloudsat_index):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
            else:
                cloudsat_index = int((i+calipso_displacement)/config.CLOUDSAT_TRACK_RESOLUTION + 0.5)
                if len(clsatObj.cloudsat.MODIS_cloud_flag) > (cloudsat_index-1):
                    cal_MODIS_cflag[i] = clsatObj.cloudsat.MODIS_cloud_flag[cloudsat_index]
                else:
                    cal_MODIS_cflag[i] = -9
    else:
        cal_MODIS_cflag = None
                
    # Now calculate cloud emissivity Ec for topmost CALIOP and CloudSat layer

    #    Ec = numpy.zeros(cal_lat_ok.shape[0],'d')
    Ec = numpy.zeros(cal_data_ok.shape[0],'d')
    dum1,cwnum,dum2=get_central_wavenumber(noaa_number,273.15)
    if caObj.avhrr.surftemp != None:
        #    for i in range(cal_lat_ok.shape[0]):
        for i in range(cal_data_ok.shape[0]):
            if cal_data_ok[i]:
                I=tb2radiance_using_central_wavenumber_klm(cwnum,caObj.avhrr.bt11micron[i],noaa_number,'4')
                Iclear=tb2radiance_using_central_wavenumber_klm(cwnum,caObj.avhrr.surftemp[i],noaa_number,'4')
                if I > Iclear:
                    Iclear = I   # Just avoiding too much mismatch between forecasted and real surface temps
                BTc=tb2radiance_using_central_wavenumber_klm(cwnum,caliop_max_height_midlaytemp[i],noaa_number,'4')
                if BTc < Iclear:
                    Ec[i]=(I-Iclear)/(BTc-Iclear)
                else:
                    Ec[i]=-9.0 # Give up on all temperature inversion cases!
            else:
                Ec[i]=-9.0                   
    if process_mode == 'EMISSFILT': #Apply only when cloud heights are above EMISS_MIN_HEIGHT
        #print "Emissivity filtering applied!"
        caliop_min_height_ok = numpy.greater(caliop_max_height, config.EMISS_MIN_HEIGHT)
        emissfilt_calipso_ok = numpy.logical_or(numpy.logical_and(numpy.greater(Ec, config.EMISS_LIMIT),caliop_min_height_ok),numpy.logical_or(numpy.equal(caliop_max_height,-9.),numpy.less_equal(caliop_max_height, config.EMISS_MIN_HEIGHT)))
    ##########################################################################################################################################
    ### 1 KM DATA CWC-RVOD ###                       
    elif config.RESOLUTION == 1 and cloudsat_type == 'CWC-RVOD':
        elevationcwc = numpy.where(numpy.less_equal(clsatObj.cloudsatcwc.elevation,0),
                            -9, clsatObj.cloudsatcwc.elevation)

        data_okcwc = numpy.ones(clsatObj.cloudsatcwc.elevation.shape,'b')
                
    ### 5 KM DATA CWC-RVOD ###                       
    elif config.RESOLUTION == 5 and cloudsat_type == 'CWC-RVOD':   
        elevationcwc = numpy.where(numpy.less_equal(clsatObj.cloudsat5kmcwc.elevation,0),
                            -9, clsatObj.cloudsat5kmcwc.elevation,-9)

        data_okcwc = numpy.ones(clsatObj.cloudsat5kmcwc.elevation.shape,'b')   
    #==============================================================================================
    #Draw plot
    pltttt = 0
    if pltttt == 1:
#    if process_mode_dnt in config.PLOT_MODES:
        plotpath = "%s/%s/%ikm/%s/%s/%s" %(config.PLOT_DIR, base_sat, config.RESOLUTION, base_year, base_month, config.AREA)
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
            
        #trajectorypath = "%s/trajectory_plot/%ikm/%s" %(MAIN_RUNDIR,int(config.RESOLUTION),AREA)         
        trajectorypath = "%s/trajectory_plot" %(plotpath)
        if not os.path.exists(trajectorypath):
                os.makedirs(trajectorypath)
        trajectoryname = "%s/%skm_%s_trajectory" %(trajectorypath,int(config.RESOLUTION),basename)
        # To make it possible to use the same function call to drawCalClsatGEOPROFAvhrr*kmPlot
        # in any processing mode:
        if 'emissfilt_calipso_ok' not in locals():
            emissfilt_calipso_ok = None
        if clsatObj == None:
            file_type = ['eps', 'png']
            drawCalClsatAvhrrPlotTimeDiff(calat, None, caObj.diff_sec_1970, plotpath, basename, config.RESOLUTION, file_type)
            drawCalClsatAvhrrPlotSATZ(calat, None, caObj.avhrr.satz, plotpath, basename, config.RESOLUTION, file_type)
            plotSatelliteTrajectory(calon,calat,trajectoryname, file_type)
            drawCalClsatGEOPROFAvhrrPlot(None, None, caObj.calipso, caObj.avhrr, None, data_ok,
                                            None, caliop_base,
                                            caliop_height, cal_data_ok,
                                            avhrr_ctth_cal_ok, plotpath,
                                            basename, process_mode, emissfilt_calipso_ok, file_type)
        else:                    
            if cloudsat_type=='GEOPROF':
                file_type = ['eps', 'png']
                drawCalClsatAvhrrPlotSATZ(cllat, clsatObj.avhrr.satz, caObj.avhrr.satz, plotpath, basename, config.RESOLUTION, file_type)
                drawCalClsatGEOPROFAvhrrPlot(clsatObj.cloudsat, clsatObj.avhrr, caObj.calipso, caObj.avhrr, elevation, data_ok,
                                                CALIPSO_DISPLACED, caliop_base,
                                                caliop_height, cal_data_ok,
                                                avhrr_ctth_cal_ok, plotpath,
                                                basename, process_mode, emissfilt_calipso_ok, file_type)
                drawCalClsatAvhrrPlotTimeDiff(cllat, clsatObj.diff_sec_1970, caObj.diff_sec_1970, plotpath, basename, config.RESOLUTION, file_type)
                plotSatelliteTrajectory(cllon,cllat,trajectoryname, file_type)
                
            elif cloudsat_type=='CWC-RVOD':
                drawCalClsatAvhrrPlotTimeDiff(cllat, clsatObj.diff_sec_1970, caObj.diff_sec_1970, plotpath, basename, config.RESOLUTION)
                phase='LW'  
                drawCalClsatCWCAvhrrPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase) #caObj, CALIPSO_DISPLACED, caliop_base, caliop_height, cal_data_ok, avhrr_ctth_cal_ok)
                phase='IW'  
                drawCalClsatCWCAvhrrPlot(clsatObj, elevationcwc, data_okcwc, plotpath, basename, phase) #caObj, CALIPSO_DISPLACED, caliop_base, caliop_height, cal_data_ok, avhrr_ctth_cal_ok)
                
    #================================================================================================
    #Calculate Statistics
    if cloudsat_type=='GEOPROF':
        if process_mode == 'EMISSFILT':
            process_calipso_ok = emissfilt_calipso_ok
        else:
            process_calipso_ok = 0

        CalculateStatistics(process_mode, clsatObj, statfile, caObj, cal_MODIS_cflag,
                            cal_vert_feature, avhrr_ctth_csat_ok, data_ok,
                            cal_data_ok, avhrr_ctth_cal_ok, caliop_max_height,
                            process_calipso_ok, dnt_flag)
