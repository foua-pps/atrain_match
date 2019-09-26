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
import os
from glob import glob
import numpy as np
from scipy import ndimage
from matchobject_io import (readCaliopImagerMatchObj,
                            CalipsoImagerTrackObject)
from utils.get_flag_info import get_pixels_where_test_is_passed
from my_dir import ADIR


def print_common_stats(match_calipso, use, name_dict):
    nlay =np.where(match_calipso.calipso.all_arrays['number_layers_found']>0, 1, 0)
    meancl=ndimage.filters.uniform_filter1d(nlay*1.0, size=3)
    isCalipsoCloudy = np.logical_and(
        nlay > 0,
        match_calipso.calipso.all_arrays['cloud_fraction']>0.5)
    isCalipsoCloudy = np.logical_and(
        isCalipsoCloudy,
        match_calipso.calipso.all_arrays['total_optical_depth_5km']>0.15)
    # isCalipsoClear = match_calipso.calipso.all_arrays['cloud_fraction']<0.5
    isCalipsoClear = nlay == 0
    isCalipsoClear = np.logical_and(isCalipsoClear, meancl < 0.01)
    isCalipsoClear = np.logical_and(
        isCalipsoClear,
        match_calipso.calipso.all_arrays['total_optical_depth_5km']<0)
    isCalipsoSnowIce = np.logical_and(
        isCalipsoClear,
        match_calipso.calipso.all_arrays['nsidc_surface_type']>50)
    print np.sum(isCalipsoSnowIce)
    isCalipsoNotSnowIce = np.logical_and(
        isCalipsoClear,
        match_calipso.calipso.all_arrays['nsidc_surface_type']<=0)

    test_1(match_calipso, isCalipsoCloudy, isCalipsoClear)
    test_2(match_calipso, isCalipsoCloudy, isCalipsoClear)
    test_3(match_calipso, isCalipsoCloudy, isCalipsoClear)
    test_4(match_calipso, isCalipsoCloudy, isCalipsoClear)


def test_1(match_calipso, isCloudy, isClear):
    test_is_on = get_pixels_where_test_is_passed(
        match_calipso.imager.all_arrays['cma_testlist1'], bit_nr=5)
    import matplotlib.pyplot as plt
    feature1 = match_calipso.imager.all_arrays['text_t11t12']
    feature2 = match_calipso.imager.all_arrays['bt11micron']
    bad = np.logical_and(test_is_on, isClear)
    good = np.logical_and(test_is_on, isCloudy)
    fig = plt.figure()
    plt.plot( feature2[bad], feature1[bad], 'b.')
    plt.plot( feature2[good], feature1[good], 'g.')
    plt.show()
    fig.savefig("pps_investigation_test1.png", format = 'png')


def test_2(match_calipso, isCloudy, isClear):
    test_is_on = get_pixels_where_test_is_passed(
        match_calipso.imager.all_arrays['cma_testlist1'], bit_nr=11)
    import matplotlib.pyplot as plt
    feature1 = match_calipso.imager.all_arrays['bt11micron'] - match_calipso.imager.all_arrays['bt37micron']
    feature2 = match_calipso.imager.all_arrays['bt37micron']
    bad = np.logical_and(test_is_on, isClear)
    good = np.logical_and(test_is_on, isCloudy)
    fig = plt.figure()
    plt.plot( feature2[bad], feature1[bad], 'b.')
    plt.plot( feature2[good], feature1[good], 'g.')
    plt.show()
    fig.savefig("pps_investigation_test2.png", format = 'png')


def test_3(match_calipso, isCloudy, isClear):
    test_is_on = get_pixels_where_test_is_passed(
        match_calipso.imager.all_arrays['cma_testlist1'], bit_nr=12)
    import matplotlib.pyplot as plt
    feature1 = match_calipso.imager.all_arrays['text_t11t12']
    feature2 = match_calipso.imager.all_arrays['bt11micron']
    bad = np.logical_and(test_is_on, isClear)
    good = np.logical_and(test_is_on, isCloudy)
    fig = plt.figure()
    plt.plot( feature2[bad], feature1[bad], 'b.')
    plt.plot( feature2[good], feature1[good], 'g.')
    plt.show()
    fig.savefig("pps_investigation_test3.png", format = 'png')


def test_4(match_calipso, isCloudy, isClear):
    test_is_on = get_pixels_where_test_is_passed(
        match_calipso.imager.all_arrays['cma_testlist2'], bit_nr=4)
    import matplotlib.pyplot as plt
    feature1 = match_calipso.imager.all_arrays['bt11micron'] - match_calipso.imager.all_arrays['bt12micron']
    feature2 = match_calipso.imager.all_arrays['surftemp']
    bad = np.logical_and(test_is_on, isClear)
    good = np.logical_and(test_is_on, isCloudy)
    fig = plt.figure()
    plt.plot( feature2[bad], feature1[bad], 'b.')
    plt.plot( feature2[good], feature1[good], 'g.')
    plt.show()
    fig.savefig("pps_investigation_test4.png", format = 'png')


ROOT_DIR = (ADIR + "/DATA_MISC/reshaped_files/"
            "global_modis_14th_created20161108/")
ROOT_DIR_GAC = (ADIR + "/DATA_MISC/reshaped_files/"
            "ATRAIN_RESULTS_GAC/Reshaped_Files/noaa18/")

files = glob(ROOT_DIR + "Reshaped_Files/merged/modis*h5")
files = glob(ROOT_DIR_GAC + "5km/2009/*/*/*h5")

test_list_file = ADIR + "/git/acpg_develop/acpg/pges/cloudmask/pps_pge01_cmasktests.h"
TEST_NAMEFILE = open(test_list_file, 'r')
name_dict = {'cma_testlist0': {},
             'cma_testlist1': {},
             'cma_testlist2': {},
             'cma_testlist3': {},
             'cma_testlist4': {},
             'cma_testlist5': {}}
for line in TEST_NAMEFILE:
    if 'define SM_ACMG_' in line:
        list_of_line = line.split(' ')
        (name, list_nr, bit) = (list_of_line[1].replace('SM_ACMG_', ''), list_of_line[2].replace(', ', ''), list_of_line[3])
        var = 'cma_testlist' + list_nr
        name_dict[var][int(bit)] = name

match_calipsoPPS = CalipsoImagerTrackObject()
for filename in files:
    match_calipsoPPS += readCaliopImagerMatchObj(filename)
    print "Scene %s"%(os.path.basename(filename))
    use = match_calipsoPPS.imager.all_arrays['bt11micron']>-9
print_common_stats(match_calipsoPPS, use, name_dict)
