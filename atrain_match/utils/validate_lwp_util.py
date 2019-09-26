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
"""
Validation functions

"""

import numpy as np
import logging
logger = logging.getLogger(__name__)

#  Default threshold for lwp screening [g m**-2]
LWP_THRESHOLD = 170
LWP_THRESHOLD_CPP = 3000
DO_PLOT = False


def get_lwp_diff(aObj, val_subset, threshold=LWP_THRESHOLD):
    retv = get_lwp_diff_inner(aObj, val_subset, threshold=LWP_THRESHOLD)
    return retv[0]


def get_lwp_diff_inner(aObj, val_subset, threshold=LWP_THRESHOLD):
    """
     Screen lwp pixels based on *sea* mask, amsr < threshold, and cpp_lwp < 0.

     Returns array with the differences for used selected pixels.
    """
    use_sea = np.logical_or(aObj.imager.fractionofland <= 0,
                            aObj.amsr.imager_linnum_nneigh <= 0)  # might have less than 8 neighbours
    use_phase = np.logical_or(aObj.imager.cpp_phase == 1,
                              aObj.amsr.imager_linnum_nneigh <= 0)  # might have less than 8 neighbours
    # exclude very high values
    aObj.imager.cpp_lwp[aObj.imager.cpp_lwp > LWP_THRESHOLD_CPP] = -9

    use_lwp = np.logical_or(aObj.imager.cpp_lwp >= 0,
                            aObj.amsr.imager_linnum_nneigh <= 0)  # might have less than 8 neighbours
    use_lwp_upper = np.logical_or(aObj.imager.cpp_lwp < LWP_THRESHOLD_CPP,
                                  aObj.amsr.imager_linnum_nneigh <= 0)

    # use = use_sea
    use = np.logical_and(use_sea, use_phase)
    use = np.logical_and(use, use_lwp)
    # use = np.logical_and(use, use_lwp_upper)
    selection = use.all(axis=-1)
    selection = np.logical_and(val_subset, selection)
    # import pdb; pdb.set_trace()
    cpp_lwp = aObj.imager.cpp_lwp
    n_cpp = np.sum(cpp_lwp >= 0, axis=-1)  # before
    cpp_lwp[cpp_lwp < 0] = 0
    cpp_lwp[np.isnan(cpp_lwp)] = 0
    sum_cpp = np.sum(cpp_lwp, axis=-1)
    n_cpp[n_cpp == 0] = 1.0
    cpp_mean = sum_cpp * 1.0 / n_cpp
    lwp_diff = cpp_mean - aObj.amsr.lwp
    use_amsr = np.logical_and(aObj.amsr.lwp >= 0,
                              aObj.amsr.lwp < threshold)
    selection = np.logical_and(use_amsr, selection)
    selection = np.logical_and(cpp_mean >= 0, selection)
    selection = np.logical_and(cpp_mean < LWP_THRESHOLD_CPP, selection)
    selection = np.logical_and(aObj.imager.sunz < 72, selection)

    return [lwp_diff[selection], cpp_mean[selection], aObj.amsr.lwp[selection], selection]


"""


def validate_all(filenames):
    from atrain_match..plotting import plot_hist, density, distribution_map
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
                   bins=range(0, 171))
    fig2.axes[0].set_xlabel('CPP cwp (g m**-2)')
    fig2.axes[0].set_ylabel('AMSR-E lwp (g m**-2)')
    fig2.suptitle("Restrictions: %s\nNumber of pixels: %d" %
                 ('; '.join(restrictions), cwp.size))
    fig2.savefig('density_all.pdf')

    # Map of pixel distribution
    fig3 = distribution_map(lon, lat)
    fig3.suptitle("Distribution of valid pixels\n" +
                  # ("Restrictions: %s\n" % '; '.join(restrictions)) +
                  "Number of Pixels: %d" % lon.size)
    fig3.savefig('distribution_all.pdf')

    return mean, median, std, lwp_diff


"""
