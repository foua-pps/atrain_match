'''
Created on Oct 19, 2010

@author: a001696
'''
class MatchupError(Exception):
    """This exception is used when a problem with matching AVHRR data with 
    Cloudsat / CALIPSO data has occured."""
    pass


def elements_within_range(compare, base, range):
    """Compare arrays *compare* and *base*, elementwise. Returns an array with
    elements set to True if compare[i] is within (base[i]-range, base[i]+range),
    otherwise false."""
    import numpy
    
    c = numpy.array(compare)
    b = numpy.array(base)
    
    return numpy.logical_and(c > b - range, c < b + range)