#!/usr/bin/env python
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
# -*- coding: utf-8 -*-

# Author(s):

#  Adam.Dybbroe

"""Read all matched data and make some plotting
"""

from glob import glob
import os.path
import numpy as np
from scipy import histogram

ROOT_DIR = "/local_disk/laptop/NowcastingSaf/FA/cloud_week_2013may/atrain_matchdata"


files = glob(ROOT_DIR + "/2012/??/arctic_europe_1km/*h5")

from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)

match_calipso = CalipsoImagerTrackObject()
for filename in files:
    print os.path.basename(filename)
    match_calipso = match_calipso + readCaliopImagerMatchObj(filename)

isCloud = match_calipso.calipso.all_arrays['number_layers_found'] > 0
r16 = match_calipso.imager.all_arrays['r16micron']
t37 = match_calipso.imager.all_arrays['bt37micron']
t86 = match_calipso.imager.all_arrays['bt86micron']
t11 = match_calipso.imager.all_arrays['bt11micron']
t12 = match_calipso.imager.all_arrays['bt12micron']

# nsidc_surface_type
# 255 = Coast
# 1-100 = Sea ice concentration %
# 101 = Permamnent ice
# 0 = Ice free
isIce = np.logical_and(match_calipso.calipso.all_arrays['nsidc_surface_type'] >= 15,
                       match_calipso.calipso.all_arrays['nsidc_surface_type'] <= 100)
isWater = np.logical_and(np.equal(match_calipso.calipso.all_arrays['igbp_surface_type'], 17),
                         np.equal(isIce, False))
isClearWater = np.logical_and(isWater, np.equal(isCloud, False))
isClearIce = np.logical_and(isIce, np.equal(isCloud, False))

nodata = np.logical_or(t11 <= -9, t12 <= -9)
t11t12 = np.ma.array(t11 - t12, mask=nodata)

# Masking out what is not a cloud ==> gives all cloudy pixels
t11t12_cloud = np.ma.masked_where(isCloud == False, t11t12).compressed()
t11t12_clear = np.ma.masked_where(isCloud == True, t11t12).compressed()
t11t12_clear_ice = np.ma.masked_where(isClearIce == False, t11t12).compressed()
t11t12_clear_water = np.ma.masked_where(isClearWater == False, t11t12).compressed()



import matplotlib.pyplot as plt

fig = plt.figure(figsize=(9, 8))
ax = fig.add_subplot(311)

n, bins, patches = ax.hist(t11t12_cloud,
                           100, range=[-10, 20],
                           normed=1, facecolor='green', alpha=0.75,
                           label='cloudy')
ax.set_ylabel('Frequency')
ax.set_title('T11-T12: Cloudy and clear pixels according to Caliop')
ax.set_xlim(-8, 16)
ax.legend()

ax = fig.add_subplot(312)
n, bins, patches = ax.hist(t11t12_clear_water,
                           100, range=[-10, 20],
                           normed=1, facecolor='blue', alpha=0.75,
                           label='open water (ice conc < 15%')
ax.legend()
ax.set_ylabel('Frequency')
ax.set_xlim(-8, 16)

ax = fig.add_subplot(313)
n, bins, patches = ax.hist(t11t12_clear_ice,
                           100, range=[-10, 20],
                           normed=1, facecolor='red', alpha=0.75,
                           label='cloudfree sea ice')
ax.legend()
ax.set_ylabel('Frequency')
ax.set_xlim(-8, 16)
ax.set_xlabel('T11-T12 (K)')

plt.savefig('./t11t12.png')
