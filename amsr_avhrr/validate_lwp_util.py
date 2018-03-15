"""
Validation functions

"""

import numpy as np
import logging
logger = logging.getLogger(__name__)


#: Default threshold for lwp screening [g m**-2]
LWP_THRESHOLD = 170
LWP_THRESHOLD_CPP = 3000
DO_PLOT = False

def get_lwp_diff(aObj, val_subset, threshold=LWP_THRESHOLD):
    """
    Screen lwp pixels based on *sea* mask, amsr < threshold, and cpp_lwp < 0.
    
    Returns array with the differences for used selected pixels.    
    """
    use_sea = np.logical_or(aObj.avhrr.fractionofland <=0,
                            aObj.amsr.imager_linnum_nneigh <=0) # might have less than 8 neighbours                          
    use_phase = np.logical_or(aObj.avhrr.cpp_phase == 1,
                              aObj.amsr.imager_linnum_nneigh <=0) # might have less than 8 neighbours   
    use_lwp = np.logical_or(aObj.avhrr.cpp_lwp>0,
                            aObj.amsr.imager_linnum_nneigh <=0)  # might have less than 8 neighbours 
    use_lwp_upper = np.logical_or(aObj.avhrr.cpp_lwp<LWP_THRESHOLD_CPP,
                            aObj.amsr.imager_linnum_nneigh <=0)
    use = np.logical_and(use_sea, use_phase)
    use = np.logical_and(use, use_lwp)
    use = np.logical_and(use, use_lwp_upper)
    selection = use.all(axis=-1)
    selection = np.logical_and(val_subset, selection)
    #import pdb; pdb.set_trace()
    cpp_lwp = aObj.avhrr.cpp_lwp
    cpp_lwp[cpp_lwp<0] = 0
    n_cpp = np.sum(cpp_lwp>0, axis=-1)
    sum_cpp = np.sum(cpp_lwp, axis=-1)
    cpp_mean = sum_cpp * 1.0/ n_cpp 

    lwp_diff = cpp_mean - aObj.amsr.lwp 

    use_amsr = np.logical_and(aObj.amsr.lwp >0 ,
                              aObj.amsr.lwp < threshold)
    selection = np.logical_and(use_amsr,  selection)
    selection = np.logical_and(cpp_mean<LWP_THRESHOLD_CPP,  selection)
    selection = np.logical_and(aObj.avhrr.sunz<72,  selection)
    
    return lwp_diff[selection]




"""
def validate_all(filenames):
    from .plotting import plot_hist, density, distribution_map
    mean = lwp_diff.mean()
    median = np.median(lwp_diff)
    std = lwp_diff.std()
    
    print("Restrictions: %s" % '; '.join(restrictions))
    print("Number of pixels: %d" % lwp_diff.size)
    print("Mean: %.2f" % mean)
    print("Median: %.2f" % median)
    print("Standard deviation: %.2f" % std)
    

    
    # Density plot
    fig2 = density(cwp, lwp,
                   bins=xrange(0, 171))
    fig2.axes[0].set_xlabel('CPP cwp (g m**-2)')
    fig2.axes[0].set_ylabel('AMSR-E lwp (g m**-2)')
    fig2.suptitle("Restrictions: %s\nNumber of pixels: %d" %
                 ('; '.join(restrictions), cwp.size))
    fig2.savefig('density_all.pdf')
    
    # Map of pixel distribution
    fig3 = distribution_map(lon, lat)
    fig3.suptitle("Distribution of valid pixels\n" +
                  #("Restrictions: %s\n" % '; '.join(restrictions)) +
                  "Number of Pixels: %d" % lon.size)
    fig3.savefig('distribution_all.pdf')
    
    return mean, median, std, lwp_diff
"""
