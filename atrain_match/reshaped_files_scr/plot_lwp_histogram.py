# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readAmsrImagerMatchObj,
                            DataObject,
                            AmsrImagerTrackObject)
from my_dir import ADIR
from utils.validate_lwp_util import (get_lwp_diff)


def plot_hist_lwp(lwp_diff, filename):
    from histogram_plotting import plot_hist
    hist_range = (np.percentile(lwp_diff, 1),
                  np.percentile(lwp_diff, 99))
    fig = plot_hist(lwp_diff, bins=100, range=hist_range)
    fig.axes[0].set_xlabel('lwp difference (g m**-2)')
    fig.suptitle("CPP lwp - AMSR-E lwp npix: %d" % lwp_diff.size)
    my_path, my_file = os.path.split(filename)
    fig.savefig(my_path + "/fig2_" + my_file.replace('h5', 'pdf'))


filename = ADIR + "/FromCollegues/forJanFokke/for_JanFokke/before_sza/1km_eos2_20100414_1040_00000_amsr_modis_match.h5"
# filename = ADIR + "/FromCollegues/forJanFokke/for_JanFokke/before_sza/1km_meteosat9_20100414_1045_99999_amsr_seviri_match.h5"

if __name__ == "__main__":
    aObj = readAmsrImagerMatchObj(filename)
    val_subset = np.bool_(np.ones(aObj.amsr.latitude.shape))
    lwp_diff = get_lwp_diff(aObj, val_subset)
    plot_hist_lwp(lwp_diff, filename)
