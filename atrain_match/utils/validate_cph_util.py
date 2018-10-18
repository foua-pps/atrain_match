"""
Use this script to validate the CPP cloud phase (cph) product

"""
import numpy as np

#: Calipso cloud phase bits
CALIPSO_PHASE_BITS = range(5, 7)

#: Calipso cloud phase values
CALIPSO_PHASE_VALUES = dict(unknown=0,
                            ice=1,
                            water=2,
                            horizontal_oriented_ice=3)

#: Water (no mixed) value
CALIPSO_WATER_VALUE = 2

#: Calipso quality bits
CALIPSO_QUAL_BITS = range(7, 9)

#: Calipso quality values
CALIPSO_QUAL_VALUES = dict(none=0,
                           low=1,
                           medium=2,
                           high=3)


#: CPP cph value meanings, from v2014
CPP_PHASE_VALUES = dict(liquid=1,
                        ice=2,
                        no_data=255)
#: CPP cph value meanings, in v2012
CPP_PHASE_VALUES_v2012 = dict(no_cloud=0,
                        liquid=1,
                        ice=2,
                        mixed=3,
                        non_opice=4,
                        non_opwater=5,
                        uncertain=6,
                        no_observation=-1)


#: Cloud type phase value equivalents
CTYPE_PHASE_BITS = {'Not processed or undefined': 1,
                    'Water': 2,
                    'Ice': 4,
                    'Tb11 below 260K': 8}




def get_bits(value, bits, shift=False):
    """
    Returns value for bits *bits* in *value*.
    
    Examples
    
    >>> get_bits(6, [0, 1])
    2
    >>> get_bits(6, [1, 2])
    6
    
    If *shift* is True, shift the obtained value by min(bits) bits:
    
    >>> get_bits(6, [1, 2], shift=True)
    3
    
    """
    selected = value & sum([2**i for i in bits])
    if shift:
        return selected >> min(bits)
    return selected

def get_calipso_phase_inner(features, qual_min=CALIPSO_QUAL_VALUES['medium'],
                            max_layers=1, same_phase_in_top_three_lay=True):
    """
    Returns Calipso cloud phase.    
    Pixels with quality lower than *qual_min* are masked out.    
    Screen out pixels with more than *max_layers* layers.    
    """
    if same_phase_in_top_three_lay:
        phase1 = get_bits(features[:,0], CALIPSO_PHASE_BITS, shift=True)
        phase2 = get_bits(features[:,1], CALIPSO_PHASE_BITS, shift=True)
        phase3 = get_bits(features[:,2], CALIPSO_PHASE_BITS, shift=True)
        two_layer_pixels = features[:, 2] >1
        three_layer_pixels = features[:, 3] >1
        lay1_lay2_differ = np.logical_and(two_layer_pixels,
                                          np.not_equal(phase1, phase2))
        lay2_lay3_differ = np.logical_and(three_layer_pixels,
                                          np.not_equal(phase2, phase3))
        varying_phases_in_top_3lay = np.logical_or(lay1_lay2_differ,
                                                      lay2_lay3_differ)
    # Reduce to single layer, masking any multilayer pixels
    features = np.ma.array(features[:, 0],
                           mask=(features[:, max_layers:] > 1).any(axis=-1))
    if same_phase_in_top_three_lay:
        features = np.ma.array(features,                               
                                mask = varying_phases_in_top_3lay)
    phase = get_bits(features, CALIPSO_PHASE_BITS, shift=True)
    qual = get_bits(features, CALIPSO_QUAL_BITS, shift=True)    
    # Don't care about pixels with lower than *qual_min* quality
    return np.ma.array(phase, mask=qual < qual_min)



    #from utils.plotting import distribution_map
    #fig = distribution_map(lon, lat)
    #fig.suptitle("Distribution of valid pixels in cloud phase validation\n" +
    #              "Number of Pixels: %d" % lon.size)
    #fig.savefig('cph_distribution_all.pdf')


if __name__ == '__main__':
    pass
