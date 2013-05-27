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

DATADIR = "/local_disk/laptop/NowcastingSaf/FA/cloud_week_2013may/atrain_matchdata/2012/10/arctic_europe_1km"

files = glob(DATADIR + "/*h5")

from matchobject_io import (readCaliopAvhrrMatchObj,
                            CalipsoAvhrrTrackObject)

caObj = CalipsoAvhrrTrackObject()
for filename in files:
    print os.path.basename(filename)

    caObj = caObj + readCaliopAvhrrMatchObj(filename)

    
