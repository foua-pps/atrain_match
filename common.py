'''
Created on Oct 19, 2010

@author: a001696
'''
import numpy as np
import logging
logger = logging.getLogger(__name__)

class Cross:
    """A cross where two satellites (almost) meet."""
    def __init__(self, satellite,  time):
        self.satellite1 = satellite
        self.satellite2 = 'xxx'
        self.time = time
    def __repr__(self):
        return "Cross: %s at time %s" %(self.satellite1, self.time)
    def __str__(self):
        return self.__repr__()

class MatchupError(Exception):
    """This exception is used when a problem matching AVHRR data with 
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

def map_avhrr(avhrr, lon, lat, radius_of_influence, n_neighbours=1):
    """
    Map AVHRR object *avhrr* to (lon, lat).
    
    A better use of this function would be to return *mapper*! But the calling
    functions would need some adjustment...
    
    """
    from config import NODATA
    from amsr_avhrr.match import match_lonlat
    source = (avhrr.longitude, avhrr.latitude)
    target = (lon, lat)
    #if avhrr.longitude.dtype != lon.dtype or  avhrr.latitude.dtype != lat.dtype:
    source = (avhrr.longitude.astype(np.float64), 
              avhrr.latitude.astype(np.float64))
    target = (lon.astype(np.float64), lat.astype(np.float64))   
    #print avhrr.longitude.dtype, lon.dtype, avhrr.latitude.dtype,  lat.dtype    
    mapper = match_lonlat(source, target, radius_of_influence, 
                          n_neighbours=n_neighbours)    
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
    from config import COMPRESS_LVL, WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE
    import h5py
    from matchobject_io import the_used_variables 
    with h5py.File(filename, 'w') as f:
        f.create_dataset('diff_sec_1970', data=diff_sec_1970,
                         compression=COMPRESS_LVL)
        
        for group_name, group_object in groups.items():
            g = f.create_group(group_name)
            for array_name, array in group_object.items():
                if array is None:
                    continue
                if (WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE and 
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
                    #print "writing", array_name
                    g.create_dataset(array_name, data=array,
                                     compression=COMPRESS_LVL)



if __name__ == "__main__":
    pass
