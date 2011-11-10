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
    import pylab
    import numpy
    
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                llcrnrlon=-180,urcrnrlon=180,resolution='c')

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

    pylab.legend((m1),['CloudSat/Calipso'],loc=0)
    
    figname = '%s.%s' %(trajectoryname, fig_type)
    pylab.savefig(figname)
    



