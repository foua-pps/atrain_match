#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Read Calipso/VIIRS/AVHRR matchup data object from hdf5 file
"""

TESTDIR = "/local_disk/laptop/NowcastingSaf/FA/cloud_week_2013may/atrain_matchdata/2012/10/arctic_europe_1km"
import os.path
TESTFILE = os.path.join(TESTDIR, 
                        "1km_npp_20121012_1246_04968_caliop_viirs_match.h5")


class DataObject(object):
    """
    Class to handle data objects with several arrays.
    
    """
    def __getattr__(self, name):
        try:
            return self.all_arrays[name]
        except KeyError:
            raise AttributeError("%s instance has no attribute '%s'" % (self.__class__.__name__, name))
    
    def __setattr__(self, name, value):
        if name == 'all_arrays':
            object.__setattr__(self, name, value)
        else:
            self.all_arrays[name] = value
            
class ppsAvhrrObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'sec_1970': None,
            'ctth_height': None,
            'ctth_pressure': None,
            'ctth_temperature': None,
            'ctth_opaque': None,  # True if opaque retrieval was applied
            'cloudtype': None,
            'cloudtype_qflag': None,
            'surftemp': None,
            'r06micron': None,
            'r09micron': None,
            'bt37micron': None,
            'bt11micron': None,
            'bt12micron': None,
            'satz': None,
            'lwp': None
            }
        
class CalipsoObject(DataObject):
    def __init__(self):
        DataObject.__init__(self)                            
        self.all_arrays = {
            'longitude': None,
            'latitude': None,
            'avhrr_linnum': None,
            'avhrr_pixnum': None,
            'cloud_fraction': None,
            'cloud_top_profile': None,
            'cloud_base_profile': None,
            'cloud_mid_temperature': None,
            'number_of_layers_found': None,
            'igbp': None,
            'nsidc': None,
            'elevation': None,
            'time': None,
            'utc_time': None, 
            'sec_1970': None,
            'feature_classification_flags': None,
            'day_night_flag': None,
            'optical_depth': None,
            'optical_depth_uncertainty': None,
            'single_shot_cloud_cleared_fraction': None,
            'Horizontal_Averaging': None
            }
class CalipsoAvhrrTrackObject:
    def __init__(self):
        self.avhrr=ppsAvhrrObject()
        self.calipso=CalipsoObject()
        self.diff_sec_1970=None

# ----------------------------------------
def readCaliopAvhrrMatchObj(filename):
    import h5py #@UnresolvedImport
    
    retv = CalipsoAvhrrTrackObject()
    
    h5file = h5py.File(filename, 'r')
    for group, data_obj in [(h5file['/calipso'], retv.calipso),
                            (h5file['/avhrr'], retv.avhrr)]:
        for dataset in group.keys():        
            if dataset in data_obj.all_arrays.keys():
                data_obj.all_arrays[dataset] = group[dataset].value

    retv.diff_sec_1970 = h5file['diff_sec_1970'].value

    h5file.close()

    return retv

# ----------------------------------------
def writeCaliopAvhrrMatchObj(filename, ca_obj):
    """
    Write *ca_obj* to *filename*.
    
    """
    from common import write_match_objects
    groups = {'calipso': ca_obj.calipso.all_arrays,
              'avhrr': ca_obj.avhrr.all_arrays}
    write_match_objects(filename, ca_obj.diff_sec_1970, groups)    

    status = 1
    return status

# ----------------------------------------
if __name__ == "__main__":

    caObj = readCaliopAvhrrMatchObj(TESTFILE)

