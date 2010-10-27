# program process_master.py

# source /data/proj/saf/ejohansson/runscripts/source_erik_cmsaf
# python process_master.py

# Writer Erik Johansson
import os, string
import time
from pps_error_messages import write_log
import setup

#--------------------------------------------------------------------------------------------------------
def MonthLessTen(month_org):
    if int(month_org) < 10:
        month_new = "0%s" %month_org
    else:
        month_new = "%s" %month_org
    
    return month_new

#--------------------------------------------------------------------------------------------------------
def AvhrrSec(avhrrfile_name):
    avhrr_split = os.path.basename(avhrrfile_name).split("_")
    avhrr_year = avhrr_split[1][0:4]
    avhrr_month = avhrr_split[1][4:6]
    avhrr_day = avhrr_split[1][6:8]
    avhrr_hour = avhrr_split[2][0:2]
    avhrr_min = avhrr_split[2][2:4]
    
    avhrr_date = "%s %s %s %s %s" %(avhrr_year,avhrr_month,avhrr_day,avhrr_hour,avhrr_min)
    avhrr_date = time.strptime(avhrr_date, "%Y %m %d %H %M")
    nsec = time.mktime(avhrr_date)
    
    return(nsec)

#-------------------------------------------------------------------------------------------------------- 

def FindAvhrrFile(match_avhrr):
    fixed_avhrr = string.capwords(match_avhrr)
    fixed_avhrr_split = fixed_avhrr.split(" ")
    fixed_avhrr_split[5] = str(int(float(fixed_avhrr_split[5])))
    fixed_avhrr = string.join(fixed_avhrr_split," ")
    match_avhrr_date = time.strptime(fixed_avhrr,"%Y %m %d %H %M %S")
    match_avhrr_sec = time.mktime(match_avhrr_date)
    avhrr_dir_date = [] # to return to save dubbel work
    
    # if it is the first day of the month then add the files from the month before
    if fixed_avhrr_split[2]=='1':
        if fixed_avhrr_split[1]=='1':
            year_avhrr_temp = str(int(fixed_avhrr_split[0])-1)
            month_avhrr_temp = '12'
        else:
            year_avhrr_temp = fixed_avhrr_split[0]
            month_avhrr_temp = str(int(fixed_avhrr_split[1])-1)
        
        month_avhrr_temp = MonthLessTen(month_avhrr_temp) # Month before
        avhrr_dir_temp = "%s/%s/%skm/%s/%s/export" %(setup.SAT_DIR,avhrr_sat,setup.RESOLUTION,year_avhrr_temp,month_avhrr_temp) #Dir month before
        all_avhrr = os.listdir(avhrr_dir_temp)  #Files month before
        month_avhrr = MonthLessTen(fixed_avhrr_split[1]) # month
        avhrr_dir = "%s/%s/%skm/%s/%s/export" %(setup.SAT_DIR,avhrr_sat,setup.RESOLUTION,fixed_avhrr_split[0],month_avhrr) # Dir month
        all_avhrr_temp = os.listdir(avhrr_dir) #file month
        all_avhrr.extend(all_avhrr_temp) #All files
        avhrr_dir_date.append(month_avhrr_temp)
        avhrr_dir_date.append(year_avhrr_temp)
        avhrr_dir_date.append(month_avhrr)
        avhrr_dir_date.append(fixed_avhrr_split[0])

    # if it is the "last" days of the month then add the files from the month after. Last day > 28. This is the smalest last day. Easy solution
    elif int(fixed_avhrr_split[2]) >= 28:
        if fixed_avhrr_split[1]=='12': 
            year_avhrr_temp = str(int(fixed_avhrr_split[0])+1)
            month_avhrr_temp = '1'
        else:
            year_avhrr_temp = fixed_avhrr_split[0]
            month_avhrr_temp = str(int(fixed_avhrr_split[1])+1)
        
        month_avhrr = MonthLessTen(fixed_avhrr_split[1]) # month
        avhrr_dir = "%s/%s/%skm/%s/%s/export" %(setup.SAT_DIR,avhrr_sat,setup.RESOLUTION,fixed_avhrr_split[0],month_avhrr) # dir month
        all_avhrr = os.listdir(avhrr_dir) #file month
        month_avhrr_temp = MonthLessTen(month_avhrr_temp) # month after
        avhrr_dir_temp = "%s/%s/%skm/%s/%s/export" %(setup.SAT_DIR,avhrr_sat,setup.RESOLUTION,year_avhrr_temp,month_avhrr_temp) #dir month after
        all_avhrr_temp = os.listdir(avhrr_dir_temp) # files month after
        all_avhrr.extend(all_avhrr_temp) # All files
        avhrr_dir_date.append(month_avhrr)
        avhrr_dir_date.append(fixed_avhrr_split[0])
        avhrr_dir_date.append(month_avhrr_temp)
        avhrr_dir_date.append(year_avhrr_temp)
    else:
        month_avhrr = MonthLessTen(fixed_avhrr_split[1])
        avhrr_dir = "%s/%s/%skm/%s/%s/export" %(setup.SAT_DIR,avhrr_sat,setup.RESOLUTION,fixed_avhrr_split[0],month_avhrr)
        all_avhrr = os.listdir(avhrr_dir)
        avhrr_dir_date.append(month_avhrr)
        avhrr_dir_date.append(fixed_avhrr_split[0])
        avhrr_dir_date.append(month_avhrr)
        avhrr_dir_date.append(fixed_avhrr_split[0])
    all_avhrr.sort()
    
    
    k=[]
    for i in range(len(all_avhrr)):
        if all_avhrr[i].split(".h5")[0].split("_")[-1]=='ctth':
            filesec = AvhrrSec(all_avhrr[i])
            k.append(i)
            #print((match_avhrr_sec-filesec)/60)
            if (filesec+(0*60)) > match_avhrr_sec:
                avhrr_file = all_avhrr[k[-2]]
                break
    #pdb.set_trace()
    return avhrr_file, avhrr_dir_date         
        
#--------------------------------------------------------------------------------------------------------

def FindFiles(avhrrfile, avhrr_dir_date=None):
    if avhrr_dir_date is None:
        # Use the file_finders package for finding files
        import file_finders
        
        pps_finder = file_finders.PpsFileFinder(setup.SAT_DIR, time_window=5*60)
        try:
            pps_finder.set_subdir_method(setup.subdir)
        except AttributeError:
            pass
        parsed = pps_finder.parse(avhrrfile)
        satname = parsed['satellite']
        datetime = parsed['datetime']
        cloudsat_finder = file_finders.CloudsatFileFinder(setup.CLOUDSAT_DIR,
                                                          setup.RESOLUTION,
                                                          setup.CLOUDSAT_TYPE)
        cloudsat_finder.set_time_window(setup.SAT_ORBIT_DURATION + setup.sec_timeThr)
        cloudsat_files = sorted(cloudsat_finder.find(datetime))
        calipso_finder = file_finders.CalipsoFileFinder(setup.CALIPSO_DIR,
                                                        setup.RESOLUTION)
        calipso_finder.set_time_window(setup.SAT_ORBIT_DURATION + setup.sec_timeThr)
        calipso_files = sorted(calipso_finder.find(datetime))
        cloudtype_file = pps_finder.find(satname, datetime, ending='cloudtype.h5')[0]
        ctth_file = pps_finder.find(satname, datetime,
                                    ending='%s.h5' % setup.CTTH_FILE)[0]
        avhrr_file = avhrrfile
        nwp_tsur_file = pps_finder.find(satname, datetime, ending='nwp_tsur.h5')[0]
        sunsatangles_file = pps_finder.find(satname, datetime, ending='sunsatangles.h5')[0]
        
        return (cloudsat_files, calipso_files, cloudtype_file, ctth_file, 
                avhrr_file, nwp_tsur_file, sunsatangles_file)
    
    avhrr_split = os.path.basename(avhrrfile).split("_")
    avhrr_year = avhrr_split[1][0:4]
    avhrr_month = avhrr_split[1][4:6]
    #avhrr_day = avhrr_split[1][6:8]
    #avhrr_hour = avhrr_split[2][0:2]
    #avhrr_min = avhrr_split[2][2:4]
    
    #avhrr_date = "%s %s %s %s %s" %(avhrr_year,avhrr_month,avhrr_day,avhrr_hour,avhrr_min)
    #avhrr_date = time.strptime(avhrr_date, "%Y %m %d %H %M")
    #avhrr_sec = time.mktime(avhrr_date)-sec_timeThr
    avhrr_sec = AvhrrSec(avhrrfile)
    avhrr_sec = avhrr_sec-setup.sec_timeThr
    new_date = time.localtime(avhrr_sec)
    # Controls if more than one cloudsat/calipso dir is necessary
    all_cl = []
    all_cal = []

    if avhrr_dir_date[0]!=avhrr_dir_date[2]:
        CL_DIR = "%s/CloudSat/%skm/%s/%s" %(setup.SAT_DIR,setup.RESOLUTION,avhrr_dir_date[1],avhrr_dir_date[0])
        CAL_DIR = "%s/Calipso/%skm/%s/%s" %(setup.SAT_DIR,setup.RESOLUTION,avhrr_dir_date[1],avhrr_dir_date[0])
        CL_DIR_1 = "%s/CloudSat/%skm/%s/%s" %(setup.SAT_DIR,setup.RESOLUTION,avhrr_dir_date[3],avhrr_dir_date[2])
        CAL_DIR_1 = "%s/Calipso/%skm/%s/%s" %(setup.SAT_DIR,setup.RESOLUTION,avhrr_dir_date[3],avhrr_dir_date[2])
        all_cl = os.listdir(CL_DIR)
        all_cl_1 = os.listdir(CL_DIR_1)
        all_cl.extend(all_cl_1)
        
        all_cal = os.listdir(CAL_DIR)
        all_cal_1 = os.listdir(CAL_DIR_1)
        all_cal.extend(all_cal_1)
    else:
        CL_DIR = "%s/CloudSat/%skm/%s/%s" %(setup.SAT_DIR,setup.RESOLUTION,avhrr_dir_date[1],avhrr_dir_date[0])
        CAL_DIR = "%s/Calipso/%skm/%s/%s" %(setup.SAT_DIR,setup.RESOLUTION,avhrr_dir_date[1],avhrr_dir_date[0])
        all_cl = os.listdir(CL_DIR)
        all_cal = os.listdir(CAL_DIR)
  
    all_cl_geo = []
    all_cl_cwc = []
    
    #This just beqause the file in one and 5 km data set have different names
    if setup.RESOLUTION==1:
        clsat_type_place=3
        clsat_time_place=0
    elif setup.RESOLUTION==5:
        clsat_type_place=4
        clsat_time_place=1
        
    for i in range(len(all_cl)):
                
        cloudsat_type = all_cl[i].split(".h5")[0]
        cloudsat_type = string.join(cloudsat_type.split("_")[clsat_type_place].split("-")[1:],"-")
        if cloudsat_type == 'GEOPROF':
            all_cl_geo.append(all_cl[i])
        elif cloudsat_type == 'CWC-RVOD':
            all_cl_cwc.append(all_cl[i])
        
    
    if setup.clsat_type==1:
        all_clsat=all_cl_geo
    elif setup.clsat_type==2:
        all_clsat=all_cl_cwc
        
    CLOUDSAT_DIR=[]
    ncl=0  

    for j in range(len(all_clsat)):
        all_clsat.sort()
        cloudsat_split = all_clsat[j].split("_")
        cloudsat_date = "%s %s %s %s %s" %(cloudsat_split[clsat_time_place][0:4], cloudsat_split[clsat_time_place][4:7],cloudsat_split[clsat_time_place][7:9],cloudsat_split[clsat_time_place][9:11],cloudsat_split[clsat_time_place][11:13])
        cloudsat_date = time.strptime(cloudsat_date, "%Y %j %H %M %S")
        cloudsat_sec = time.mktime(cloudsat_date)
        if cloudsat_date[1] == int(avhrr_dir_date[0]):
            CLOUDSAT_DIR.append(CL_DIR)
        else:
            CLOUDSAT_DIR.append(CL_DIR_1)
            
        if cloudsat_sec >= avhrr_sec and ncl==0:
            #clsat_file = "'%s/%s'   '%s/%s'    '%s/%s'" %(CLOUDSAT_DIR[j-1],all_clsat[j-1],CLOUDSAT_DIR[j],all_clsat[j],CLOUDSAT_DIR[j+1],all_clsat[j+1])
            clsat_file = ['%s/%s' % (CLOUDSAT_DIR[j-1], all_clsat[j-1])]
            clsat_file.append('%s/%s' % (CLOUDSAT_DIR[j], all_clsat[j]))
            ncl=ncl+1
        elif cloudsat_sec >= avhrr_sec and ncl==1:
            clsat_file.append('%s/%s' % (CLOUDSAT_DIR[j], all_clsat[j]))
            ncl=ncl+1
        if ncl>1:
            break
    
            
    
    CALIPSO_DIR=[] 
    ncal=0 
    for a in range(len(all_cal)):
        all_cal.sort()
        cal_split = all_cal[a].split("-")
        cal_year = cal_split[3].split(".")[1]
        cal_month = cal_split[4]
        cal_day = cal_split[5].split("T")[0]
        cal_hour = cal_split[5].split("T")[1]
        cal_min = cal_split[6]
        cal_sec = cal_split[7][0:2]
        cal_date = "%s %s %s %s %s %s" %(cal_year,cal_month,cal_day,cal_hour,cal_min,cal_sec)
        cal_date = time.strptime(cal_date, "%Y %m %d %H %M %S")
        cal_tot_sec = time.mktime(cal_date)
        if cal_date[1] == int(avhrr_dir_date[0]):
            CALIPSO_DIR.append(CAL_DIR)
        else:
            CALIPSO_DIR.append(CAL_DIR_1)
        
        if cal_tot_sec >= avhrr_sec and ncal==0:         
            cal_file = ['%s/%s' % (CALIPSO_DIR[a-2], all_cal[a-2])]
            cal_file.append('%s/%s' % (CALIPSO_DIR[a-1], all_cal[a-1]))
            cal_file.append('%s/%s' % (CALIPSO_DIR[a], all_cal[a]))
            ncal=ncal+1
        elif cal_tot_sec >= avhrr_sec and ncal>=1 and ncal<4:
            cal_file.append('%s/%s' % (CALIPSO_DIR[a], all_cal[a]))
            ncal=ncal+1
        if ncl >3:
            break
            
    AVHRR_DIR = "%s/%s/%skm/%s/%s" %(setup.SAT_DIR,avhrr_sat,setup.RESOLUTION,avhrr_year,avhrr_month)
    #test=all_cal.sort()
    avhrr_join = string.join(avhrr_split[0:-1],"_")
    cloudtype_file = "%s/export/%s_cloudtype.h5" %(AVHRR_DIR, avhrr_join)
    ctth_file = "%s/export/%s_ctth.h5" %(AVHRR_DIR, avhrr_join)
    avhrr_file = "%s/import/%s_avhrr.h5" %(AVHRR_DIR, avhrr_join)
    nwp_tsur_file = "%s/import/%s_nwp_tsur.h5" %(AVHRR_DIR, avhrr_join)
    sunsatangles_file = "%s/import/%s_sunsatangles.h5" %(AVHRR_DIR, avhrr_join)
    
    raise Warning
    #   avhrr_date = (avhrr_year,avhrr_month,avhrr_day,avhrr_hour,avhrr_min)
    return clsat_file, cal_file, cloudtype_file, ctth_file, avhrr_file, nwp_tsur_file, sunsatangles_file

#------------------------------------------------------------------------------------------------------------------  
if __name__=='__main__':
    #pdb.set_trace()
    from optparse import OptionParser
    import find_crosses
    import file_finders
    import cloudsat_calipso_avhrr_match
    from common import MatchupError
    
    parser = OptionParser()
    parser.set_usage("usage: %prog [options]\n"
                     "Some influential environment variables:\n"
                     "SAT_DIR                   Base directory where satellite data files are stored.\n"
                     "VALIDATION_RESULTS_DIR    Base directory where results will be stored.\n"
                     "                          Used indirectly by cloudsat_calipso_avhrr_match.py.")
    parser.add_option('-m', '--matches', type='string', metavar='FILE',
                      default=setup.match_file, help="Use FILE for matchups (SNO output)")
    parser.add_option('-M', '--mode', type='string', action='append',
                      help="Run validation software in MODE (valid modes are %s)" % \
                      ', '.join(setup.ALLOWED_MODES))
    parser.add_option('-d', '--debug', action='store_true')
    (options, args) = parser.parse_args()
    
    if options.mode is not None:
        run_modes = options.mode
    else:
        run_modes = setup.ALLOWED_MODES
    
    matchups = find_crosses.parse_crosses_file(options.matches)
    if matchups[0].satellite1 is not None:
        # Assume AVHRR satellite name in satellite1 (first argument to SNO executable)
        avhrr_sat = matchups[0].satellite1.lower()
    #match_times_file = open(options.matches, "r")

    #match_times_list = match_times_file.readlines()
    #match_times_file.close()
    resolution = "%i.%i" %(setup.RESOLUTION,setup.clsat_type)
    
    avhrr_finder = file_finders.PpsFileFinder(setup.SAT_DIR, 'avhrr.h5')
    if matchups[0].time_window is not None:
        avhrr_finder.set_time_window(-(setup.SAT_ORBIT_DURATION + matchups[0].time_window), matchups[0].time_window)
    else:
        avhrr_finder.set_time_window(-(setup.SAT_ORBIT_DURATION + setup.sec_timeThr), setup.sec_timeThr)
    avhrr_finder.set_subdir_method(setup.subdir)
    
    problem_files = set()
    no_matchup_files = []
    no_matchup_cross = []
    processed_avhrr_files = set()
    for match in matchups:
        
        try:
            # A bit backwards, but makes it possible to use find_crosses for parsing
            # the SNO output file. We don't have to prune the file in advance...
            #(avhrr_file, avhrr_dir_date) = FindAvhrrFile(match.as_sno_line()[:24])
            
            avhrr_file = avhrr_finder.find(match.satellite1.lower(), match.time1)[0]
            parsed = avhrr_finder.parse(avhrr_file)
            
            #avhrr_dir_date = ("%02d" % parsed['datetime'].month, "%d" % parsed['datetime'].year,
            #                  "%02d" % parsed['datetime'].month, "%d" % parsed['datetime'].year)
            avhrr_dir_date = None
        except IndexError:
            write_log('WARNING', "No matching AVHRR file found for %s %s." % \
                      (match.satellite1, match.time1.strftime("%F %T")))
            no_matchup_cross.append(match)
            continue
        if avhrr_file in processed_avhrr_files:
            write_log('WARNING', "Matching AVHRR file has already been processed: %s" % avhrr_file)
            continue
        processed_avhrr_files.add(avhrr_file)
        
        # This is what I would like to do instead...
        # import file_finders
        # pps_finder = file_finders.PpsFileFinder(basedir=???, ending='avhrr.h5')
        #avhrr_file = pps_finder.find(match)
        try:
            cloudsat_file, calipso_file, cloudtype_file, ctth_file, avhrr_file, nwp_tsur_file, sunsatangles_file = FindFiles(avhrr_file, avhrr_dir_date)
        except IndexError:
            problem_files.add(avhrr_file)
            write_log('WARNING', "Couldn't find all needed files for %s." % avhrr_file)
            if options.debug is True:
                raise
            break
        
        for mode in run_modes:
            try:
                cloudsat_calipso_avhrr_match.run(cloudsat_file, calipso_file,
                                                 cloudtype_file, ctth_file,
                                                 avhrr_file, nwp_tsur_file,
                                                 sunsatangles_file, mode, 
                                                 resolution)
            except MatchupError, err:
                write_log('WARNING', "Matchup problem for %s %s: %s" % \
                          (match.satellite1, match.time1.strftime("%F %T"),
                           err.message))
                no_matchup_files.append(avhrr_file)
                break
            except:
                problem_files.add(avhrr_file)
                write_log('WARNING', "Couldn't run cloudsat_calipso_avhrr_match.")
                if options.debug is True:
                    raise
                break
    if len(no_matchup_cross) > 0:
        write_log('WARNING', "%d of %d cases had no matching AVHRR files:\n%s" % \
                  (len(no_matchup_cross), len(matchups),
                   '\n'.join([str(c) for c in no_matchup_cross])))
    if len(no_matchup_files) > 0:
        write_log('WARNING', "%d of %d cases had no matchups in region, within the time window:\n%s" % \
                  (len(no_matchup_files), len(matchups), '\n'.join(no_matchup_files)))
    if len(problem_files) > 0:
        write_log('WARNING', "%d of %d cases had unknown problems:\n%s" % \
                  (len(problem_files), len(matchups), '\n'.join(problem_files)))