"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readAmsrAvhrrMatchObj,
                            DataObject,
                            AmsrAvhrrTrackObject)

from amsr_avhrr.validate_lwp_util import ( get_lwp_diff)
def plot_hist_lwp(lwp_diff, filename):
    from histogram_plotting import plot_hist
    hist_range = (np.percentile(lwp_diff, 1),
                  np.percentile(lwp_diff, 99))
    fig = plot_hist(lwp_diff, bins=100, range=hist_range)
    fig.axes[0].set_xlabel('lwp difference (g m**-2)')
    fig.suptitle("CPP lwp - AMSR-E lwp npix: %d"% lwp_diff.size)
    my_path, my_file = os.path.split(filename)
    fig.savefig(my_path + "/fig2_" + my_file.replace('h5','pdf'))

filename = "/home/a001865/FromCollegues/forJanFokke/for_JanFokke/before_sza/1km_eos2_20100414_1040_00000_amsr_modis_match.h5"
#filename = "/home/a001865/FromCollegues/forJanFokke/for_JanFokke/before_sza/1km_meteosat9_20100414_1045_99999_amsr_seviri_match.h5"

if __name__ == "__main__":
    aObj = readAmsrAvhrrMatchObj(filename)
    val_subset = np.bool_(np.ones(aObj.amsr.latitude.shape))
    lwp_diff = get_lwp_diff(aObj, val_subset)
    plot_hist_lwp(lwp_diff,filename)
