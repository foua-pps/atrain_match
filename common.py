'''
Created on Oct 19, 2010

@author: a001696
'''
import numpy as np

class MatchupError(Exception):
    """This exception is used when a problem matching AVHRR data with 
    Cloudsat / CALIPSO data has occured."""
    pass

class TimeMatchError(Exception):
    """This exception is used when the time in a file is not 
    the same as the time in the filename."""
    pass


def elements_within_range(compare, base, _range):
    """Compare arrays *compare* and *base*, elementwise. Returns an array with
    elements set to True if compare[i] is within (base[i]-_range, base[i]+_range),
    otherwise false."""
    import numpy
    
    c = numpy.array(compare)
    b = numpy.array(base)
    
    return numpy.logical_and(c > b - _range, c < b + _range)


def attach_subdir_from_config(finder):
    """Attach ``config.subdir`` to file finder *finder*."""
    import config
    
    try:
        subdir = config.subdir
        finder.set_subdir_method(subdir)
    except AttributeError:
        pass


def map_avhrr(avhrr, lon, lat, radius_of_influence):
    """
    Map AVHRR object *avhrr* to (lon, lat).
    
    A better use of this function would be to return *mapper*! But the calling
    functions would need some adjustment...
    
    """
    from config import NODATA
    from amsr_avhrr.match import match_lonlat
    source = (avhrr.longitude, avhrr.latitude)
    target = (lon, lat)
    if avhrr.longitude.dtype != lon.dtype:
        source = (avhrr.longitude.astype(np.float64), avhrr.latitude.astype(np.float64))
        target = (lon.astype(np.float64), lat.astype(np.float64))      
    mapper = match_lonlat(source, target, radius_of_influence, n_neighbours=1)
    
    # Return the nearest (and the only calculated) neighbour
    #return mapper.rows.filled(NODATA)[:, 0], mapper.cols.filled(NODATA)[:, 0]
    # Nina 2016-01-19 changed mapper.rows to be 1D arrays not 2D-arrays with 
    # one column as that is mostly "needed for array magic."
    return mapper.rows.filled(NODATA)[:], mapper.cols.filled(NODATA)[:]

def write_match_objects(filename, diff_sec_1970, groups):
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
    ...           'avhrr': ca_obj.avhrr.all_arrays}
    >>> write_match_objects('match.h5', ca_obj.diff_sec_1970, groups)
    
    The match object data can then be read using `read_match_objects`:
    
    >>> diff_sec_1970, groups = read_match_objects('match.h5')
    
    """
    from config import COMPRESS_LVL
    import h5py
    with h5py.File(filename, 'w') as f:
        f.create_dataset('diff_sec_1970', data=diff_sec_1970,
                         compression=COMPRESS_LVL)
        
        for group_name, group_object in groups.items():
            g = f.create_group(group_name)
            for array_name, array in group_object.items():
                if array is None:
                    continue
                try:
                    if len(array) == 0:
                        continue
                except:
                    # Scalar data can't be compressed
                    # TODO: Write it as and attribute instead?
                    g.create_dataset(array_name, data=array)
                else:
                    g.create_dataset(array_name, data=array,
                                     compression=COMPRESS_LVL)


#This function "read_match_objects()" is not used. 
#Not sure when it was last used. 20160314 a001865
#Function in matchobject_io is used for reading match objects. 
#Removed the code, look for it in git if needed.
