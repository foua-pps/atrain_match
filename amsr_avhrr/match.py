"""
Match AMSR and AVHRR data

"""

from __future__ import with_statement
import numpy as np
import logging
import h5py
from config import RESOLUTION, NODATA
logger = logging.getLogger(__name__)


class MatchMapper(object):
    """
    Map arrays from one swath to another.
    
    Note that MatchMapper always works with an extra dimension: neighbour
    
    """
    def __init__(self, rows, cols, pixel_mask, time_diff=None,
                 time_threshold=None):
        self._rows = rows
        self._cols = cols
        self._pixel_mask = pixel_mask
        self._time_diff = time_diff
        self.time_threshold = time_threshold
    
    def __call__(self, array):
        """
        Maps *array* to target swath.
        
        """
        return np.ma.array(array[self.rows, self.cols], mask=self.mask)
    
    @property
    def rows(self):
        #self._rows = np.array(self._rows, dtype=np.int64)
        return np.ma.array(self._rows, mask=self.mask, fill_value=NODATA,
                           hard_mask=True)
    
    @property
    def cols(self):

        #self._cols = np.array(self._cols, dtype=np.int64)
        return  np.ma.array(self._cols, mask=self.mask, fill_value=-NODATA,
                            hard_mask=True)
    
    @property
    def time_diff(self):
        """Time difference in seconds"""
        if self._time_diff is None:
            return None
        # Only use pixel mask
        return np.ma.array(self._time_diff, mask=self._pixel_mask,
                           fill_value=np.inf, hard_mask=True)
    
    @time_diff.setter
    def time_diff(self, value):
        self._time_diff = value
    
    @property
    def mask(self):
        if not None in (self.time_diff, self.time_threshold):
            return (self._pixel_mask +
                    (abs(self.time_diff) > self.time_threshold))
        return self._pixel_mask
    
    def write(self, filename, compression=True):
        """
        Write mapper to hdf5 file *filename*.
        
        """
        with h5py.File(filename, 'w') as f:
            f.create_dataset('rows', data=self.rows.filled(),
                             compression=compression)
            f.create_dataset('cols', data=self.cols.filled(),
                             compression=compression)
            f.create_dataset('pixel_mask', data=self._pixel_mask,
                             compression=compression)
            if self.time_diff is not None:
                f.create_dataset('time_diff', data=self.time_diff.filled(),
                                 compression=compression)
            if self.time_threshold is not None:
                f.attrs['time_threshold'] = self.time_threshold
        if f:
            f.close()
    
    @classmethod
    def from_file(cls, filename):
        """
        Create a mapper from contents of *filename*.
        
        """
        with h5py.File(filename, 'r') as f:
            rows = f['rows'][:]
            cols = f['cols'][:]
            pixel_mask = f['pixel_mask'][:]
            if 'time_diff' in f.keys():
                time_diff = f['time_diff'][:]
                if 'time_threshold' in f['time_diff'].attrs.keys():
                    time_threshold = f['time_diff'].attrs['time_threshold']
                elif  'time_threshold' in f.attrs.keys():
                    time_threshold = f.attrs['time_threshold']
                else:
                    print "Did not find time_threshold!"
                    time_threshold = None
            else:
                time_diff = None
                time_threshold = None
        if f:
            f.close()
        
        return cls(rows=rows, cols=cols, pixel_mask=pixel_mask,
                   time_diff=time_diff, time_threshold=time_threshold)


def match_lonlat(source, target,
                 radius_of_influence=0.7*RESOLUTION*1000.0,
                 n_neighbours=1):
    """
    Produce a masked array of the same shape as the arrays in *target*, with
    indices of nearest neighbours in *source*. *source* and *target* should be
    tuples (lon, lat) of the source and target swaths, respectively.
    
    Note::
    
        * Fastest matching is obtained when *target* has lower resolution than
        *source*.
        
        * *source* should have 2-dimensional lon and lat arrays.
    
    """
    from pyresample.geometry import SwathDefinition
    from pyresample.kd_tree import get_neighbour_info
    from pyresample.kd_tree import get_sample_from_neighbour_info

    lon, lat = source
    mask_out_lat = np.logical_or(lat<-90, lat>90)
    mask_out_lon = np.logical_or(lon>180, lat<-180)
    mask_out = np.logical_or(mask_out_lat, mask_out_lon)
    lat = np.ma.masked_array(lat, mask=mask_out)
    lon = np.ma.masked_array(lon, mask=mask_out)
    #lat = np.around(lat, decimals=4)
    #lon = np.around(lon, decimals=4)
    source_def = SwathDefinition(*(lon,lat))
    target_def = SwathDefinition(*target)
    logger.debug("Matching %d nearest neighbours", n_neighbours)
    valid_in, valid_out, indices, distances = get_neighbour_info(
        source_def, target_def, radius_of_influence, neighbours=n_neighbours)
    #Use pyresampe code to find colmun and row numbers for each pixel
    #This is works also with no-data in imager lat/lon.
    cols_matrix, rows_matrix = np.meshgrid(np.array(xrange(0,lat.shape[1])),
                                           np.array(xrange(0,lat.shape[0])))
    cols = get_sample_from_neighbour_info('nn', target_def.shape,
                                          cols_matrix,
                                          valid_in, valid_out,
                                          indices)
    rows = get_sample_from_neighbour_info('nn', target_def.shape,
                                          rows_matrix,
                                          valid_in, valid_out,
                                          indices)
    rows = np.array(rows)
    cols = np.array(cols)

                                          
    """ Code used during debugging, leaving it here for now
    #Hopfully not needed anymore as indices is not used directly
    if indices.dtype in ['uint32']:
        #With pykdtree installed get_neighbour_info returns indices
        # as type uint32
        #This does not combine well with a nodata value of -9.
        indices = np.array(indices,dtype=np.int64)
    #get_expected_output even for nodata in lat/lon!
    #print "indices", indices
    if 1==1:
        print distances, indices
        print max(indices)
        print min(indices)
        print len(valid_in)
        print len(valid_in[valid_in])
        # But why is +1 item needed??
        from_one_to_many = np.array(xrange(0,len(valid_in)+1))
        print from_one_to_many
        valid_in_new = np.append(valid_in,np.array([True]), axis=0)
        print valid_in_new
        use_these = indices[valid_out]
        print use_these
        new_numbers = from_one_to_many[valid_in_new]
        print new_numbers
        indices[valid_out] = new_numbers[use_these]
    #print "indices", indices
    shape = list(target_def.shape)
    shape.append(n_neighbours)
    indices.shape = shape
    distances.shape = shape
    rows = indices // source_def.shape[1]
    cols = indices % source_def.shape[1]
    print "c", cols, "r", rows
    print rows.shape, cols.shape
    """

    # Make sure all indices are valid
    #import ipdb; ipdb.set_trace()
    rows[rows >= source_def.shape[0]] = NODATA
    cols[cols >= source_def.shape[1]] = NODATA
    mask = distances > radius_of_influence
    return MatchMapper(rows, cols, mask)

def match(amsr_filename, avhrr_filename, sunsat_filename, radius_of_influence=1e3,
          time_threshold=None, n_neighbours=8):
    """
    Find matching indices in AVHRR array for each element in AMSR swath.
    
    Arguments:
    
        amsr_filename: string
            full path of AMSR-E HDF5 file
        avhrr_filename: string
            full path of AVHRR PPS HDF5 file
        radius_of_influence: float
            radius of influence in meters in pixel-pixel matching (default:
            1000 m)
        time_threshold: float
            largest absolute time difference to include in match
        n_neighbours: int
            number of nearest AVHRR neighbours to use
    
    Returns:
    
        mapper: `MatchMapper` instance.
    
    """
    from .util import get_amsr_lonlat, get_avhrr_lonlat
    from .util import get_amsr_time, get_avhrr_time
    
    avhrr_lonlat = get_avhrr_lonlat(sunsat_filename)
    amsr_lonlat = get_amsr_lonlat(amsr_filename)
    
    mapper = match_lonlat(avhrr_lonlat, amsr_lonlat, radius_of_influence,
                          n_neighbours=n_neighbours)
    
    avhrr_time = get_avhrr_time(sunsat_filename)
    amsr_time = get_amsr_time(amsr_filename)
    
    time_diff = np.abs(avhrr_time[mapper.rows] -
                       amsr_time.reshape((amsr_time.size, 1, 1))).astype(np.float32)
    
    mapper.time_diff = time_diff
    mapper.time_threshold = time_threshold
    
    logger.debug("Time diff (min, max): %3.1f, %3.1f", 
                 time_diff.min(),time_diff.max())
    
    return mapper

