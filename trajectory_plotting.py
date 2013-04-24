"""New plotting routines using pyresample functionality"""

def plotSatelliteTrajectory(longitude, latitude,
                            trajectoryname, 
                            area_config_file,
                            fig_type='eps',
                            **options):
    """Plot a trajectory of geolocations (lon,lat) on a map"""
    import numpy as np
    import os
    
    if "area_id" in options:
        area_id = options['area_id']
    else:
        area_id = "mill10km_test"
        
    track = np.ones(latitude.shape)

    import pyresample as pr
    swath_def = pr.geometry.SwathDefinition(longitude, latitude)

    area_def = pr.utils.load_area(area_config_file, area_id)

    result = pr.kd_tree.resample_nearest(swath_def, track, area_def,
                                         radius_of_influence=20000, fill_value=None)

    import matplotlib.pyplot as plt
    import pylab

    bmap = pr.plot.area_def2basemap(area_def)
    bmng = bmap.bluemarble()
    col = bmap.imshow(result, cmap=pylab.get_cmap('rainbow'), origin='upper')
    plt.title('Calipso matchup track')

    for figt in fig_type:
        figname = trajectoryname + '.' + figt
        plt.savefig(figname, bbox_inches='tight')

    del plt
    return
