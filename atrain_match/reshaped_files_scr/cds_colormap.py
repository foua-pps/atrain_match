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
# Author: Nina Hakansson
import copy
import matplotlib
import numpy as np
"""A custom color map to be used for plot of scores at projected areas.
The colormap uses the bwr colormap for low values.
And for higher values the hot colormap in inverted directio is used
this adds darker red and brown colors for the really high values.
Just comment out the plotting if you do not want it.
"""

def get_cds_colormap(minv=0.0, maxv=5.0, colorchange1=0.22, colorchange2=1.0):

    
    n_resolution = 1000
    step = (maxv-minv) * 1.0 / n_resolution
    n_res_darkred = int((maxv-colorchange2) / (step*0.4)) #will use ~40% of colormap
    #copy existing colormap, important to copy!
    my_cmap = copy.copy(matplotlib.cm.get_cmap("bwr",
                                               lut=n_resolution)) #blue to red
    cmap_vals = my_cmap(np.arange(n_resolution)) #extract values as an array
    new_cmap_vals = cmap_vals.copy()
    my_cmap_darkred=copy.copy(
        matplotlib.cm.get_cmap("hot",
                               lut=n_res_darkred)) #get red and brown colors
    cmap_vals_darkred = my_cmap_darkred(np.arange(n_res_darkred)) #values
    cmap_vals_darkred = cmap_vals_darkred[::-1]

    #first_color
    steps_first_color = int((colorchange1-minv) / step)
    every_xth_needed = int(n_resolution * 0.5 / steps_first_color)
    for i in  range(steps_first_color):
        new_cmap_vals[i] = cmap_vals[i * every_xth_needed]
    #second_color
    steps_second_color = int((colorchange2-colorchange1) / step)
    every_xth_needed = int(n_resolution * 0.5 / steps_second_color)
    for i in  range(steps_second_color):
        index_br = int(n_resolution*0.5) + i*every_xth_needed
        new_cmap_vals[steps_first_color+i] = cmap_vals[index_br ]
    #third_color
    steps_second_color = int((colorchange2-colorchange1) / step)
    every_xth_needed = int(n_resolution * 0.5 / steps_second_color)
    for i in  range(n_resolution-steps_second_color-steps_first_color):
        index = steps_first_color + steps_second_color+i
        #  print(index, int(0.6 * n_res_darkred), i, n_res_darkred, int(0.6 * n_res_darkred)+i)
        new_cmap_vals[index] = cmap_vals_darkred[int(0.6*n_res_darkred)+i]

    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "newBlueRedMoreRed", new_cmap_vals)
    a = np.array([np.linspace(0, 5, num=1000)])
    a = a.reshape(20, 50)
    #ax = plt.subplot(111)
    #ax.imshow(a, cmap=my_cmap, interpolation='nearest')
    #plt.show()
    return my_cmap
