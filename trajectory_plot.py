# Program trajectory_plot





# python trajectory_plot.py

def plotSatelliteTrajectory(cllon,cllat,calon,calat,avhrlon,avhrlat,trajectoryname):
    from mpl_toolkits.basemap import Basemap
    import pylab
    #import pdb
    #pdb.set_trace()
    #from setup import AREA
    
    # TODO: This plotting function needs to be looked over...
    N = len(cllon)/500 + 2 # Number of segments (make sure we have at least two)
    upplosning = len(cllon)/N
    print('Resolution = %i' %(upplosning))
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                llcrnrlon=-180,urcrnrlon=180,resolution='c')

    m.drawmapboundary()
    m.drawcoastlines()
    
    m.drawgreatcircle(cllon[0],cllat[0],cllon[1],cllat[1],color='red')
    m1 = m.drawgreatcircle(cllon[0],cllat[0],cllon[1],cllat[1],color='red')
    #m.drawgreatcircle(calon[0],calat[0],calon[1],calat[1],color='green')
    #m.drawgreatcircle(avhrlon[0],avhrlat[0],avhrlon[1],avhrlat[1],color='blue')
    
    print("Draw CloudSat trajectory")
    for i in range(1,len(cllon)-upplosning,upplosning):#len(cllon)-1,500):
        if cllon[i] < 0 and cllon[i+upplosning] > 0:
            pass
        else:
            m.drawgreatcircle(cllon[i],cllat[i],cllon[i+upplosning],cllat[i+upplosning],color='red')
    print("Draw Calipso trajectory")
    #for j in range(1,len(calon)-upplosning,upplosning):
    #    if calon[j] < 0 and calon[j+upplosning] > 0:
    #        pass
    #    else:
    #        m.drawgreatcircle(calon[j],calat[j],calon[j+upplosning],calat[j+upplosning],color='green')
    #print("Draw Avhrr trajectory")
    #for k in range(1,len(avhrlon)-upplosning,upplosning):
    #    if avhrlon[k] < 0 and avhrlon[k+upplosning] > 0:
    #        pass
    #    else:
    #        m.drawgreatcircle(avhrlon[k],avhrlat[k],avhrlon[k+upplosning],avhrlat[k+upplosning],color='blue')    
    
    
    
    m.drawgreatcircle(cllon[i],cllat[i],cllon[-1],cllat[-1],color='red')
    #m.drawgreatcircle(calon[j],calat[j],calon[-1],calat[-1],color='green')
    #m.drawgreatcircle(avhrlon[k],avhrlat[k],avhrlon[-1],avhrlat[-1],color='blue')
    
    #pdb.set_trace()
    
    #pylab.legend(('CloudSat','Calipso','Avhrr'),loc=0)
    pylab.legend((m1),['CloudSat/Calipso'],loc=0)
    #pylab.legend(("CloudSat"),loc=0)
    
    epsfig = '%s.eps' %trajectoryname
    #pngfig = '%s.png' %trajectoryname
    pylab.savefig(epsfig)
    #pylab.savefig(pngfig)


if __name__ == "__main__":

    lon11 = -128.11390686035156
    lon12 = -162.1396484375
    lat11 = 75.187126159667969
    lat12 = 18.735273361206055

    lon21 = 4.0472955703735352
    lon22 = -50.375637054443359
    lat21 = 64.620552062988281
    lat22 = 78.766883850097656

     
    #m = Basemap(projection='ortho', lon_0=-145, lat_0=0)
    #m = Basemap(width=10000000,height=8000000,
    #            resolution='l',projection='laea',\
    #            lat_ts=50,lat_0=50,lon_0=-120.)
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                llcrnrlon=-180,urcrnrlon=180,resolution='c')

    m.drawmapboundary()
    m.drawcoastlines()

    m.drawgreatcircle(lon11,lat11,lon12,lat12,color='red')
    m.drawgreatcircle(lon21,lat21,lon22,lat22,color='blue')
    
    pylab.savefig('./trajectory_plot/test.eps')
    pylab.savefig('./trajectory_plot/test.png')



