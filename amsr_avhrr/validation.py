"""
Validation functions

"""

import numpy as np
import logging
logger = logging.getLogger(__name__)


#: Default threshold for lwp screening [kg m**-2]
LWP_THRESHOLD = 170
DO_PLOT = False

def screen_lwp(amsr_lwp, cpp_lwp, mask, screened, threshold=LWP_THRESHOLD):
    """
    Screen lwp pixels based on *sea* mask, amsr < threshold, and cpp_lwp < 0.
    
    Returns an bool array with screened out pixels set to True, suitable for
    creating masked arrays.
    
    """
    def show_mask(mask, screened):
        if DO_PLOT:
            from matplotlib import pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)
            im = ax.imshow(mask[..., 0])
            fig.colorbar(im)
            ax.set_title(', '.join(screened))
            fig.savefig('mask%d.png' % len(screened))
    
    screened = list(screened)
    
    mask |= (amsr_lwp > threshold)
    screened.append('AMSR-E lwp > %r g m**-2' % threshold)
    show_mask(mask, screened)
    
    mask |= (amsr_lwp < 0)
    screened.append('AMSR-E lwp < 0')
    show_mask(mask, screened)
    
    mask |= (cpp_lwp < 0)
    screened.append('CPP lwp < 0')
    show_mask(mask, screened)
    
    return mask, screened


def validate_lwp(amsr_lwp, cpp_lwp, selection):
    """
    Compare liquid water path, lwp, in *amsr_filename* and *cpp_filename* files.
    Use only `True` elements in selection.
    
    """
    # Use only AMSR-E pixels for which all corresponding AVHRR pixel are valid
    selection_amsr = selection.all(axis=-1)
    if not selection_amsr.any():
        return None
    
    amsr_masked = np.ma.array(amsr_lwp, mask=~selection_amsr)
    cpp_masked = np.ma.array(cpp_lwp, mask=~selection)
    
    # Use average of all AVHRR pixels in AMSR footprint
    assert len(cpp_masked.shape) == 3
    cpp_masked = cpp_masked.mean(axis=-1)
    
    lwp_diff = cpp_masked - amsr_masked
    
    print('=' * 40)
    print("CPP cwp - AMSR-E lwp")
    print("Number of pixels in comparison: %d" % lwp_diff.compressed().size)
    print("bias:    %.4g" % lwp_diff.mean())
    print("std:     %.4g" % lwp_diff.std())
    print("rel std: %.4g %%" % abs(100. * lwp_diff.std() / lwp_diff.mean()))
    
    return lwp_diff


def validate_all(filenames):
    """
    Use lwp_diff in all files in *filenames* for validation.
    
    """
    from .util import get_diff_data
    lwp_diff, restrictions, cwp, lwp, lon, lat = \
        get_diff_data(filenames, ('cpp_cwp', 'amsr_lwp', 'longitudes',
                                  'latitudes'))
    
    mean = lwp_diff.mean()
    median = np.median(lwp_diff)
    std = lwp_diff.std()
    
    print("Restrictions: %s" % '; '.join(restrictions))
    print("Number of pixels: %d" % lwp_diff.size)
    print("Mean: %.2f" % mean)
    print("Median: %.2f" % median)
    print("Standard deviation: %.2f" % std)
    
    from .plotting import plot_hist, density, distribution_map
    hist_range = (np.percentile(lwp_diff, 1),
                  np.percentile(lwp_diff, 99))
    fig = plot_hist(lwp_diff, bins=500, range=hist_range)
    fig.axes[0].set_xlabel('lwp difference (g m**-2)')
    fig.suptitle("CPP cwp - AMSR-E lwp\nRestrictions: %s\nPixels left: %d" %
                 ('; '.join(restrictions), lwp_diff.size))
    fig.savefig('validate_all.pdf')
    
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
