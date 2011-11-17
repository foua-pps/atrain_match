# Program trajectory_plot





# python trajectory_plot.py
def makesTheActuallyTrajectoryPlot(m,lon,lat):
    N = len(lon)/500 + 2 # Number of segments (make sure we have at least two)
    upplosning = len(lon)/N
    #print("Draw CloudSat trajectory")
    for i in range(1,len(lon)-upplosning,upplosning):#len(cllon)-1,500):
        if lon[i] < 0 and lon[i+upplosning] > 0:
            pass
        else:
            m.drawgreatcircle(lon[i],lat[i],lon[i+upplosning],lat[i+upplosning],color='red')
    
    m.drawgreatcircle(lon[i],lat[i],lon[-1],lat[-1],color='red')
    return(m)

def plotSatelliteTrajectory(longitude,latitude,trajectoryname,fig_type='eps'):
    # TODO: This plotting function needs to be looked over...
    
    from mpl_toolkits.basemap import Basemap #@UnresolvedImport
    from matplotlib import pyplot as plt
    import numpy
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                llcrnrlon=-180,urcrnrlon=180,resolution='c', ax=ax)

    m.drawmapboundary()
    m.drawcoastlines()
    
    lon_nanix = numpy.where(~numpy.isnan(longitude))[0]    
    lon_ix_split = numpy.where(numpy.diff(lon_nanix) != 1)[0]
    
    if lon_ix_split.shape[0]==0 and  lon_nanix.shape[0]!=0:
        lon = longitude[lon_nanix]
        lat = latitude[lon_nanix]
        m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color='red')
        m1 = m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color='red')
        m = makesTheActuallyTrajectoryPlot(m,lon,lat)
    elif lon_ix_split.shape[0]==1:
        lon = longitude[lon_nanix[0:lon_ix_split[0]+1]]
        lat = latitude[lon_nanix[0:lon_ix_split[0]+1]]
        m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color='red')
        m1 = m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color='red')
        m = makesTheActuallyTrajectoryPlot(m,lon,lat)
        
        lon = longitude[lon_nanix[lon_ix_split[-1]+1:]]
        lat = latitude[lon_nanix[lon_ix_split[-1]+1:]]
        m = makesTheActuallyTrajectoryPlot(m,lon,lat)
    elif lon_ix_split.shape[0]>1:
        lon = longitude[lon_nanix[0:lon_ix_split[0]+1]]
        lat = latitude[lon_nanix[0:lon_ix_split[0]+1]]
        m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color='red')
        m1 = m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color='red')
        m = makesTheActuallyTrajectoryPlot(m,lon,lat)
        
        for c in range(len(lon_ix_split)-1):
            lon = longitude[lon_nanix[lon_ix_split[c]+1:lon_ix_split[c+1]+1]]
            lat = latitude[lon_nanix[lon_ix_split[c]+1:lon_ix_split[c+1]+1]]
            m = makesTheActuallyTrajectoryPlot(m,lon,lat)
        
        lon = longitude[lon_nanix[lon_ix_split[-1]+1:]]
        lat = latitude[lon_nanix[lon_ix_split[-1]+1:]]
        m = makesTheActuallyTrajectoryPlot(m,lon,lat)    

    ax.legend((m1),['CloudSat/Calipso'],loc=0)
    
    if isinstance(fig_type, str) == True:
        figname = '%s.%s' %(trajectoryname, fig_type)
        fig.savefig(figname)
    else:
        for figtype in fig_type:
            figname = '%s.%s' %(trajectoryname, figtype)
            fig.savefig(figname)
    
def drawTrajectoryOfAllSNO(sno_output_file, file_type='eps'):
    """
    This function plots all trajectoies in one plot
    Limitations:
    It can only handle match files. i.e. the files that are saved when running atrain_match
    Only SNO in year 2006-2009 can be handled. Other year might be ploted or 
    there might be an error message. How knows?
    The filename is fixed.
    """
    import find_crosses
    from file_finders import CloudsatCalipsoAvhrrMatchFileFinder #@UnresolvedImport
    from mpl_toolkits.basemap import Basemap #@UnresolvedImport
    import pylab
    import config
    from calipso import readCaliopAvhrrMatchObj
        
    found_matchups = find_crosses.parse_crosses_file(sno_output_file)
    
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.drawmapboundary()
    m.drawcoastlines()
    k6 = 0
    k7 = 0
    k8 = 0
    k9 = 0
    k = 0
    
    for cross in found_matchups:
        print('file %i out of %i' %(k, len(found_matchups)))
        k = k + 1
        satellite = cross.satellite1.lower()
        datetime = cross.time1
        if satellite in ['calipso', 'cloudsat']:
            satellite = cross.satellite2.lower()
            datetime = cross.time2
        match_finder = CloudsatCalipsoAvhrrMatchFileFinder(config.RESHAPE_DIR,config.RESOLUTION,region=config.AREA)
        match_finder.set_time_window(-(config.SAT_ORBIT_DURATION + cross.time_window), 0)
        try:
            ca_match_file = match_finder.find(datetime, satellite, atrain_datatype='caliop')[0]
        except IndexError:
            #print(datetime.year)
            continue
        caObj = readCaliopAvhrrMatchObj(ca_match_file)
        
        lon = caObj.avhrr.longitude.astype('float64')
        lat = caObj.avhrr.latitude.astype('float64')
        if datetime.year == 2006:
            colour = 'red'
            if k6 == 0:
                m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color=colour, label='2006')
                k6 = 1
        elif datetime.year == 2007:
            colour = 'blue'
            if k7 == 0:
                m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color=colour, label='2007')
                k7 = 1
        elif datetime.year == 2008:
            colour = 'green'
            if k8 == 0:
                m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color=colour, label='2008')
                k8 = 1
        elif datetime.year == 2009:
            colour = 'black'
            if k9 == 0:
                m.drawgreatcircle(lon[0],lat[0],lon[1],lat[1],color=colour, label='2009')
                k9 = 1
      
        upplosning = 10
        if len(lon)<=upplosning:
            upplosning = 1
        for i in range(1,len(lon)-upplosning,upplosning):
            if lon[i] < 0 and lon[i+upplosning] > 0:
                pass
            else:
                m.drawgreatcircle(lon[i],lat[i],lon[i+upplosning],lat[i+upplosning],color=colour)
       
    
        m.drawgreatcircle(lon[i],lat[i],lon[-1],lat[-1],color=colour)
    pylab.legend()
    
    if isinstance(file_type, str) == True:
        filename = "trajektory.%s" %(file_type)
        pylab.savefig(filename, format = file_type)
        print("saves the file %s" %filename)
    else:
        for filetype in file_type:
            filename = "trajektory.%s" %(filetype)
            pylab.savefig(filename, format = filetype)
            print("saves the file %s" %filename)
            
if __name__=='__main__':
    # TODO: Fix so that the snofile is given in comand line instead of fixed
    snofile = '../SNO_tools/Snotimes/calipso_noaa18_20061001_20091231_0.2min.dat'
    file_type = ['png', 'eps']
    drawTrajectoryOfAllSNO(snofile, file_type)



       

