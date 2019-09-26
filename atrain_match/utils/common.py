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
'''
Created on Oct 19, 2010

'''
import numpy as np
import logging
logger = logging.getLogger(__name__)


class Cross:
    """A cross where two satellites (almost) meet."""

    def __init__(self, satellite, time):
        self.satellite1 = satellite
        self.satellite2 = 'xxx'
        self.time = time

    def __repr__(self):
        return "Cross: %s at time %s" % (self.satellite1, self.time)

    def __str__(self):
        return self.__repr__()


class MatchupError(Exception):
    """This exception is used when a problem matching IMAGER data with
    Cloudsat / CALIPSO data has occured."""
    pass


class TimeMatchError(Exception):
    """This exception is used when the time in a file is not
    the same as the time in the filename."""
    pass


class InputError(Exception):
    """This exception is used when the input does
    not match what is expected."""
    pass


class ProcessingError(Exception):
    """This exception is used when the processing fails."""
    pass


def elements_within_range(compare, base, _range):
    """Compare arrays *compare* and *base*, elementwise. Returns an array with
    elements set to True if compare[i] is within (base[i]-_range, base[i]+_range),
    otherwise false."""
    c = np.array(compare)
    b = np.array(base)
    return np.logical_and(c > b - _range, c < b + _range)


def map_imager_distances(imager, lon, lat, radius_of_influence, n_neighbours=1):
    """
    Map IMAGER object *imager* to (lon, lat).

    A better use of this function would be to return *mapper*! But the calling
    functions would need some adjustment...

    """
    from atrain_match.config import NODATA
    from atrain_match.utils.match import match_lonlat
    source = (imager.longitude, imager.latitude)
    target = (lon, lat)
    # if imager.longitude.dtype != lon.dtype or  imager.latitude.dtype != lat.dtype:
    source = (imager.longitude.astype(np.float64),
              imager.latitude.astype(np.float64))
    target = (lon.astype(np.float64), lat.astype(np.float64))
    # print imager.longitude.dtype, lon.dtype, imager.latitude.dtype, lat.dtype
    mapper, distances = match_lonlat(source, target, radius_of_influence,
                                     n_neighbours=n_neighbours)
    # Return the nearest (and the only calculated) neighbour
    # return mapper.rows.filled(NODATA)[:, 0], mapper.cols.filled(NODATA)[:, 0]
    # Nina 2016-01-19 changed mapper.rows to be 1D arrays not 2D-arrays with
    # Nina 2019-08-22 changed back to (n, 1) arrays needed for amsr
    # one column as that is mostly "needed for array magic."
    # Note that ravel() transform array (n, 1) array to (n, )
    # Array2D[:, 0] gives (n, )
    # np.logical_and(array_of_size(n, 1), array_of_size(n, )) => (n, n)
    out = {"mapper": (mapper.rows.filled(NODATA)[:], mapper.cols.filled(NODATA)[:]),
           "distances": distances}

    return out


def map_imager(imager, lon, lat, radius_of_influence, n_neighbours=1):
    """
    Map IMAGER object *imager* to (lon, lat).

    A better use of this function would be to return *mapper*! But the calling
    functions would need some adjustment...

    """
    retv = map_imager_distances(imager, lon, lat, radius_of_influence, n_neighbours=1)
    return retv["mapper"]


def write_match_objects(filename, datasets, groups, group_attrs_dict, SETTINGS=None):
    """
    Write match objects to HDF5 file *filename*.

    Arguments:

        *diff_sec_1970*: `numpy.ndarray`
            time diff between matched satellites
        *groups*: dict
            each key/value pair should hold a list of `numpy.ndarray` instances
            to be written under HDF5 group with the same name as the key

    E.g. to write a calipso match:

    >>> groups = {'calipso': ca_obj.calipso.all_arrays,
    ...           'imager': ca_obj.imager.all_arrays}
    >>> write_match_objects('match.h5', ca_obj.diff_sec_1970, groups)

    The match object data can then be read using `read_match_objects`:

    >>> diff_sec_1970, groups = read_match_objects('match.h5')

    """
    from atrain_match.config import COMPRESS_LVL
    import h5py
    from atrain_match.matchobject_io import the_used_variables
    with h5py.File(filename, 'w') as f:
        for name in datasets.keys():
            f.create_dataset(name, data=datasets[name],
                             compression=COMPRESS_LVL)

        for group_name, group_object in groups.items():
            skip_some = False
            if SETTINGS is not None and group_name in ['calipso',
                                                       'calipso_aerosol'
                                                       'cloudsat',
                                                       'iss',
                                                       'mora',
                                                       'synop',
                                                       'modis']:
                skip_some = SETTINGS["WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE"]

            try:
                attrs_dict = group_attrs_dict[group_name]
            except KeyError:
                attrs_dict = {}
            g = f.create_group(group_name)
            for key in attrs_dict:
                g.attrs[key] = attrs_dict[key]
            for array_name, array in group_object.items():
                if array is None:
                    continue
                if (skip_some and
                        array_name not in the_used_variables):
                    logger.debug("Not writing unimportant %s to file",
                                 array_name)
                    continue
                try:
                    if len(array) == 0:
                        continue
                except:
                    # Scalar data can't be compressed
                    # TODO: Write it as and attribute instead?
                    g.create_dataset(array_name, data=array)
                else:
                    # print "writing", array_name
                    g.create_dataset(array_name, data=array,
                                     compression=COMPRESS_LVL)


if __name__ == "__main__":
    pass
