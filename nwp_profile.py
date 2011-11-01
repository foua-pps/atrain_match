"""
Extract and plot NWP profiles

"""

import logging
logging.basicConfig(level=logging.INFO)

def profile(field, y_coords, selection):
    """
    Produce a profile plot of *field* vs *y_coords* in *selection* pixels.
    
    """
    from matplotlib import pyplot as plt
    import numpy as np
    
    from pps_nwp.fields import NWPField
    assert isinstance(field, NWPField)
    assert isinstance(y_coords, NWPField)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    lons, lats = field.lonlat() #@UnusedVariable
    lats = lats[selection[1:]] # First item in *selection* selects vertical coords
    LATS = np.array([lats for i in field[selection]]) #@UnusedVariable
    im = ax.pcolormesh(LATS, y_coords[selection], field[selection])
    if 'Pa' in y_coords.units:
        ax.invert_yaxis()
    ax.set_xlabel("latitude [deg]")
    ax.set_ylabel("%s [%s]" % (y_coords.cfName, y_coords.units))
    cbar = fig.colorbar(im, orientation='horizontal')
    cbar.set_label("%s [%s]" % (field.cfName, field.units))
    
    return fig


def pressure_to_height(pressure):
    """
    Convert pressure field to height according to the following relationship,
    found at http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html:
    
    p = 101325 * (1 - 2.25577e-5 * h)**5.25588
    <==>
    h = ((p / 101325)**(1 / 5.25588) - 1) / (-2.25577e-5)
    
    """
    from pps_nwp.fields import from_array
    
    origin = dict(method="((p / 101325)**(1 / 5.25588) - 1) / (-2.25577e-5)",
                  p=pressure.origin)
    if pressure.units == 'Pa':
        p = pressure[:]
    elif pressure.units == 'hPa':
        p = pressure[:] * 100
        origin['p scaled'] = 100
    else:
        raise ValueError("Pressure units %r not recognized" % pressure.units)
    
    h_array = ((p / 101325)**(1/5.25588) - 1) / (-2.25577e-5)
    return from_array(h_array, pressure, origin, 'height_above_sea', 'm')


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser("Usage: %prog [options] GRIBFILE lons.txt lats.txt")
    parser.add_option('-p', '--pressure', action='store_true',
                      help="Plot pressure on y-axis (default geopotential height)")
    parser.add_option('-l', '--level', type='int', help="top NWP level")
    parser.add_option('-q', '--humidity', action='store_true',
                      help="Plot specific humidity, in addition to temperature")
    (options, arguments) = parser.parse_args()
    
    grib_name, lon_name, lat_name = arguments
    
    from numpy import loadtxt
    lons = loadtxt(lon_name).reshape((-1, 1)) # pps_nwp requires 2-d lon/lat
    lats = loadtxt(lat_name).reshape((-1, 1))
    
    from pps_nwp import GRIBFile
    gribfile = GRIBFile(grib_name, (lons, lats))
    
    t = gribfile.get_t_vertical()
    if options.pressure:
        y_coords = gribfile.get_p_vertical()
    else:
        y_coords = gribfile.get_gh_vertical()
    
    selection = (slice(options.level), slice(None), 0)
    fig = profile(t, y_coords, selection)
    fig.suptitle("Temperature profile from %r" % grib_name)
    
    if options.humidity:
        q = gribfile.get_q_vertical()
        fig_q = profile(q, y_coords, selection)
        fig_q.suptitle("Specific humidity profile from %r" % grib_name)
    
    from pylab import show
    show()