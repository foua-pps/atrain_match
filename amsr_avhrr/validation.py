"""
Validation functions

"""

import numpy as np
import logging
logger = logging.getLogger(__name__)


#: Default threshold for lwp screening [kg m**-2]
LWP_THRESHOLD = 170


def screen_lwp(amsr_lwp, cpp_lwp, mask, screened, threshold=LWP_THRESHOLD):
    """
    Screen lwp pixels based on *sea* mask, amsr < threshold, and cpp_lwp < 0.
    
    Returns an bool array with screened out pixels set to True, suitable for
    creating masked arrays.
    
    """
    def show_mask(mask, screened):
        return # Don't use
        
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
    # Use pixels with at least one selected AVHRR pixel
    selection_amsr = selection.any(axis=2)
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
    import h5py
    
    lwp_diffs = []
    for filename in filenames:
        with h5py.File(filename, 'r') as f:
            data = f['lwp_diff'][:]
            restrictions = f['lwp_diff'].attrs['restrictions']
            if len(lwp_diffs) > 0:
                if not (restrictions == lwp_diffs[-1][-1]).all():
                    raise RuntimeError("Inconsistent restrictions: %r != %r" %
                                       (restrictions, lwp_diffs[-1][-1]))
            lwp_diffs.append((filename, data, restrictions))
    
    lwp_diff_array = np.hstack(data for filename, data, restrictions in lwp_diffs)
    mean = lwp_diff_array.mean()
    median = np.median(lwp_diff_array)
    std = lwp_diff_array.std()
    
    print("Restrictions: %s" % '; '.join(restrictions))
    print("Number of pixels: %d" % lwp_diff_array.size)
    print("Mean: %.2f" % mean)
    print("Median: %.2f" % median)
    print("Standard deviation: %.2f" % std)
    
    from .plotting import plot_hist
    fig = plot_hist(lwp_diff_array, bins=np.linspace(-200, 200, 500))
    fig.axes[0].set_xlabel('lwp difference (g m**-2)')
    fig.suptitle("CPP cwp - AMSR-E lwp\nRestrictions: %s\nPixels left: %d" %
                 ('; '.join(restrictions), lwp_diff_array.size))
    fig.savefig('validate_all.pdf')
    
    return mean, median, std, lwp_diff_array
