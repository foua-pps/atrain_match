# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
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
    lons, lats = field.lonlat()  # @UnusedVariable
    lats = lats[selection[1:]]  # First item in *selection* selects vertical coords
    LATS = np.array([lats for i in field[selection]])  # @UnusedVariable
    im = ax.pcolormesh(LATS, y_coords[selection], field[selection])
    if 'Pa' in y_coords.units:
        ax.invert_yaxis()
    ax.set_xlabel("latitude [deg]")
    ax.set_ylabel("%s [%s]" % (y_coords.cfName, y_coords.units))
    cbar = fig.colorbar(im, orientation='horizontal')
    cbar.set_label("%s [%s]" % (field.cfName, field.units))

    return fig


def write_profile(filename, fields, selection, y_coords=None, lonlat=None):
    """
    Write profiles *fields*, y coordinates *y_coords*, and (lon, lat) to HDF5
    file *filename*. *selection* must be a 3-tuple, used to select (levels,
    rows, cols) for writing.

    """
    from pps_nwp.fields import from_array
    for field in fields:
        a = field[selection]
        origin = dict(field=field.origin, n_selected=a.shape[1])
        from_array(a, field, origin).write(filename, field.cfName)

    if y_coords is not None:
        y = y_coords[selection]
        origin = dict(field=y_coords.origin, n_selected=a.shape[1])
        from_array(y, y_coords, origin).write(filename, 'y_coords')

    if lonlat is not None:
        lon, lat = lonlat
        import h5py
        with h5py.File(filename, 'a') as f:
            f.create_dataset('lon', data=lon[selection[1:]], compression=True)
            f.create_dataset('lat', data=lat[selection[1:]], compression=True)


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


def generate_profile(grib_name, lon_name, lat_name, pressure=False,
                     level=None, humidity=False, out_filename=None):
    from numpy import loadtxt
    lons = loadtxt(lon_name).reshape((-1, 1))  # pps_nwp requires 2-d lon/lat
    lats = loadtxt(lat_name).reshape((-1, 1))

    from pps_nwp import GRIBFile
    gribfile = GRIBFile(grib_name, (lons, lats))

    t = gribfile.get_t_vertical()
    if t.cfName == 'unknown':
        t.cfName = 'air_temperature'
    if humidity:
        q = gribfile.get_q_vertical()
        if q.cfName == 'unknown':
            q.cfName = 'air_specific_humidity'

    if pressure:
        y_coords = gribfile.get_p_vertical()
    else:
        y_coords = gribfile.get_gh_vertical()

    selection = (slice(level), slice(None), 0)

    if out_filename:
        # Just write to file
        fields = [t]
        if humidity:
            fields.append(q)
        write_profile(out_filename, fields, selection, y_coords, (lons, lats))
        return

    fig = profile(t, y_coords, selection)
    fig.suptitle("Temperature profile from %r" % grib_name)

    if humidity:
        fig_q = profile(q, y_coords, selection)
        fig_q.suptitle("Specific humidity profile from %r" % grib_name)

    from pylab import show
    show()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser("Usage: %prog [options] GRIBFILE lons.txt lats.txt")
    parser.add_option('-p', '--pressure', action='store_true',
                      help="Plot pressure on y-axis (default geopotential height)")
    parser.add_option('-l', '--level', type='int', help="top NWP level")
    parser.add_option('-q', '--humidity', action='store_true',
                      help="Plot specific humidity, in addition to temperature")
    parser.add_option('-o', '--out', type='string', metavar='FILE',
                      help="Write profile(s) to FILE")
    (opts, args) = parser.parse_args()

    try:
        grib_name, lon_name, lat_name = args
    except ValueError:
        parser.error("Wrong number of input arguments")
    generate_profile(grib_name, lon_name, lat_name, opts.pressure, opts.level,
                     opts.humidity, opts.out)
