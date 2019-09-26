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
import numpy as np
import logging
logger = logging.getLogger(__name__)


def my_hist(data, use, bmin=-3000, bmax=3000, delta_h=5,
            return_also_corresponding_gaussian=False):
    if use is None:
        use = np.ones(data.shape, dtype=bool)
    n_pix = len(data[use])
    bins = np.arange(bmin,bmax,delta_h)
    x_ = np.array([item+delta_h*0.5 for item in bins[0:-1]])
    hist_heights, bins = np.histogram(data[use],bins=bins)
    hist_heights = hist_heights*100.0/n_pix

    hist_heights_gaussian = None
    if return_also_corresponding_gaussian:
        n_pix_g = 1000000
        temp_data = np.random.normal(np.mean(data[use]),np.std(data[use]), n_pix_g)
        hist_heights_gaussian,bins = np.histogram(temp_data,bins=bins)
        hist_heights_gaussian =  hist_heights_gaussian *100.0/n_pix_g
    return hist_heights, x_, hist_heights_gaussian


def my_iqr(data):
    return np.percentile(data,75)- np.percentile(data,25)

def my_rms(data):
    return np.sqrt(np.mean(data*data))
def my_mae(data):
    return np.mean(np.abs(data))

#https:/
#https://stats.stackexchange.com/questions/278237/half-sample-mode-estimate-of-sample-of-weighted-data
def half_sample_mode(x, already_sorted=False):
    if len(x) < 3:
        return np.mean(x)
    if already_sorted:
        sorted_x = x # No need to sort
    else:
        sorted_x = np.sort(x)
    half_idx = int((len(x) + 1) / 2) # Round up to include the middle value, in the case of an odd-length array

    # Calculate all interesting ranges that span half of all data points
    ranges = sorted_x[-half_idx:] - sorted_x[:half_idx]
    smallest_range_idx = np.argmin(ranges)

    # Now repeat the procedure on the half that spans the smallest range
    x_subset = sorted_x[smallest_range_idx : (smallest_range_idx+half_idx)]
    return half_sample_mode(x_subset, already_sorted=True)

def my_pex(data, x):
    return len(data[np.abs(data)>x])*100.0/len(data)

def my_pe250m(data):
    return my_pex(data, 250)
def my_pe500m(data):
    return my_pex(data, 500)
def my_pe1000m(data):
    return my_pex(data, 1000)
def my_pe2000m(data):
    return my_pex(data, 2000)
def my_pe2500m(data):
    return my_pex(data, 2500)
def my_pe5000m(data):
    return my_pex(data, 5000)

def my_mode(bias):
    bmin = -40
    bmax = 40
    delta_h = 100.0
    bins = np.arange(bmin*1000,bmax*1000,delta_h)
    hist_heights,bins = np.histogram(bias,bins=bins)
    n_pix = len(bias)
    hist_heights = hist_heights*100.0/n_pix
    maxind = np.argmax(hist_heights)
    maxind2 = len(hist_heights)-1 - np.argmax(hist_heights[::-1])
    if maxind != maxind2:
        print( maxind, maxind2)
        raise ValueError
    mode =bins[maxind] + delta_h*0.5
    return mode


def HR_cma(indict):
    det_clear = indict["det_clear"]
    det_cloudy = indict["det_cloudy"]
    N = indict["N"]
    return (det_clear + det_cloudy)*1.0/N
def K_cma(indict):
    det_clear = indict["det_clear"]
    det_cloudy = indict["det_cloudy"]
    undet_cloudy = indict["undet_cloudy"]
    false_cloudy = indict["false_cloudy"]
    return (det_clear * det_cloudy - false_cloudy * undet_cloudy)*1.0/(
        (det_clear + false_cloudy)* (det_cloudy + undet_cloudy))
def PODcy(indict):
    det_cloudy = indict["det_cloudy"]
    undet_cloudy = indict["undet_cloudy"]
    return det_cloudy*100.0/(det_cloudy + undet_cloudy)
def FARcy(indict):
    det_cloudy = indict["det_cloudy"]
    false_cloudy = indict["false_cloudy"]
    return  false_cloudy*100.0/(false_cloudy + det_cloudy)
def PODcl(indict):
    det_clear = indict["det_clear"]
    false_cloudy = indict["false_cloudy"]
    return det_clear*100.0/(det_clear + false_cloudy)
def FARcl(indict):
    det_clear = indict["det_clear"]
    undet_cloudy = indict["undet_cloudy"]
    return undet_cloudy*100.0/(det_clear + undet_cloudy)


if __name__ == "__main__":
    pass
