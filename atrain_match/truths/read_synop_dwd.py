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

# Copyright (c) 2014, 2015, 2016 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe 

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

"""Reader and data handler for the Synop reports from DWD
"""

#TESTFILE = "./DataFromDwd/201403/sy_SYNOP_20140306.qc"
TESTFILE = "/data/proj6/saf/adybbroe/satellite_synop_matchup/DataFromDwd/201403/sy_SYNOP_20140306.qc"
filename = TESTFILE

#from astropy.io import ascii
import pandas as pd
import numpy as np
from datetime import datetime


def get_synop_data(filename):
    """Get all Synop data from one file"""

    convert_datefunc = lambda x: datetime.strptime(x, '%Y%m%d%H%M')

    dtype = [('date', object), ('station', '|S5'),
             ('lat', 'f8'), ('lon', 'f8'),
             ('msl', 'f8'), ('nix', 'i4'),
             ('pressure', 'f8'), ('temp', 'f8'),
             ('dtemp', 'f8'),
             ('total_cloud_cover', 'i4'),
             ('nh', 'i4'),
             ('cl', 'i4'),
             ('cm', 'i4'),
             ('ch', 'i4'),
             ('vvvv', 'i4'),
             ('ww', 'i4'),
             ]

    data = np.genfromtxt(filename,
                         skip_header=1,
                         skip_footer=35,
                         usecols=(
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14, 15, 28, 29),
                         dtype=dtype,
                         unpack=True,
                         converters={0: convert_datefunc,
                                     2: lambda x: float(x) / 100.,
                                     3: lambda x: float(x) / 100.,
                                     6: lambda x: float(x) / 10.,
                                     7: lambda x: float(x) / 10.,
                                     8: lambda x: float(x) / 10., })
    return pd.DataFrame(data)

if __name__ == "__main__":

    synop = get_data(TESTFILE)

    mystation = synop[synop['station'] == '01023']
    print mystation['total_cloud_cover']
