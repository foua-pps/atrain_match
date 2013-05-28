#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Read all matched data and make some plotting
"""

from glob import glob
import os.path
import numpy as np


DATADIR = "/local_disk/laptop/NowcastingSaf/FA/cloud_week_2013may/atrain_matchdata/2012/10/arctic_europe_1km"

files = glob(DATADIR + "/*h5")

from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)

caObj = CalipsoAvhrrTrackObject()
for filename in files:
    print os.path.basename(filename)

    caObj = caObj + readCaliopAvhrrMatchObj(filename)


isCloud = caObj.calipso.all_arrays['cloud_base_profile'][0] > 0.
t11 = caObj.avhrr.all_arrays['bt11micron']
t12 = caObj.avhrr.all_arrays['bt12micron']
nodata = np.logical_and(t11<=-9, t12<=-9)

t11t12_ok = np.ma.array(t11 - t12, mask=nodata)
cloud_ok =  np.logical_and(isCloud, np.equal(nodata, False))
clear_ok =  np.logical_and(np.equal(isCloud, False), 
                           np.equal(nodata, False))

t11t12_cloud = np.take(t11t12_ok, cloud_ok).data
t11t12_clear = np.take(t11t12_ok, clear_ok).data

#not_okay = np.logical_or(t11t12_cloud < -8, t11t12_cloud > 8)
#t11t12_cloud = np.ma.array(t11t12_cloud, mask=not_okay)
#t11t12_clear = np.ma.array(t11t12_clear, mask=not_okay)


import matplotlib.pyplot as plt

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(211)

n, bins, patches = ax.hist(t11t12_cloud, 
                           100, range=[-10,20],
                           normed=1, facecolor='green', alpha=0.75,
                           label='cloudy')
ax.set_ylabel('Frequency')
ax.set_title('All Cloudy - Caliop')
ax.set_xlim(-12, 22)
ax.legend()

ax = fig.add_subplot(212)
n, bins, patches = ax.hist(t11t12_clear, 
                           100, range=[-10,20],
                           normed=1, facecolor='blue', alpha=0.75,
                           label='clear')
ax.legend()
ax.set_ylabel('Frequency')
ax.set_title('All Clear - Caliop')
ax.set_xlim(-12, 22)
ax.set_xlabel('T11-T12 (K)')

plt.savefig('./t11t12.png')
