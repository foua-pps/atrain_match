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
import re
import numpy as np
from scipy import ndimage
from my_dir import ADIR

from matchobject_io import read_files
limits = [0.1*ind for ind in range(0, 10)]
limits = limits + [ind for ind in range(1, 6)]

def make_pod_vector(match_calipso):
    pod_d = []
    pod_n = []
    feature_n = []
    feature_d = []
    od = match_calipso.calipso.all_arrays["total_optical_depth_5km"]
    # try:
    #   pps_cloudy = match_calipso.imager.all_arrays['cma_prob']>50
    #
    # except:
    if match_calipso.imager.all_arrays['cloudmask'] is not None:
        pps_cloudy = np.logical_or(match_calipso.imager.all_arrays['cloudmask']==1,
                                   match_calipso.imager.all_arrays['cloudmask']==2)
    else:
        pps_cloudy = np.logical_and(np.greater(match_calipso.imager.cloudtype, 4), np.less(match_calipso.imager.cloudtype, 20))
    igbp_st = getattr(match_calipso.calipso, 'igbp_surface_type')
    alat = np.abs(match_calipso.imager.all_arrays['latitude'])
    use_all = np.logical_and(np.equal(igbp_st, 17), alat < 45)
    use_all = np.logical_and(use_all, match_calipso.calipso.all_arrays["cloud_fraction"]>=1.0)
    sunz = match_calipso.imager.all_arrays['sunz']
    if np.max(sunz) < 5:
        sunz = 100*sunz
    day = sunz < 90

    # feature = np.array(match_calipso.imager.all_arrays['bt11micron'])-match_calipso.imager.all_arrays['surftemp']
    # feature = match_calipso.imager.all_arrays['thr_t37t12'] - np.array(match_calipso.imager.all_arrays['bt37micron']) +match_calipso.imager.all_arrays['bt12micron']
    feature = np.array(match_calipso.imager.all_arrays['bt11micron'])# -match_calipso.imager.all_arrays['bt12micron'] - match_calipso.imager.all_arrays['thr_t11t12']
    feature = od
    for i, lower in enumerate(limits):
        try:
            upper = limits[i + 1]
        except:
            upper = 100000
        use = np.logical_and(use_all, np.logical_and(od >= lower, od < upper))
        use = np.logical_and(use, day)

        pod_d.append(np.sum(np.logical_and(use, pps_cloudy)) * 100.0/np.sum(use))
        feature_d.append(np.sum(feature[np.logical_and(use, np.not_equal(pps_cloudy, True))] > 297)*100.0/np.sum(use) )
        # feature_d.append(np.mean(feature[np.logical_and(use, np.equal(pps_cloudy, True))]))

        use = np.logical_and(use_all, np.logical_and(od >= lower, od < upper))
        use = np.logical_and(use, np.not_equal(day, True))
        pod_n.append(np.sum(np.logical_and(use, pps_cloudy)) * 100.0/np.sum(use))
        feature_n.append(np.sum(feature[np.logical_and(use, np.not_equal(pps_cloudy, True))] > 297)*100.0/np.sum(use) )
        # feature_n.append(np.mean(feature[np.logical_and(use, np.equal(pps_cloudy, True))]))
    use = np.logical_and(use_all, np.logical_and(od >= 0.2, od < 0.5))
    use = np.logical_and(use, np.not_equal(pps_cloudy, True))
    use_i = np.logical_and(use, day)
    from collections import Counter
    try :
        print(Counter(match_calipso.imager.all_arrays['cma_testlist0'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist1'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist2'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist3'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist4'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist5'][use_i]))
        use_i = np.logical_and(use, np.not_equal(day, True))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist0'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist1'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist2'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist3'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist4'][use_i]))
        print(Counter(match_calipso.imager.all_arrays['cma_testlist5'][use_i]))
    except:
        pass
    return np.array(pod_d), np.array(pod_n), np.array(feature_d), np.array(feature_n)


BASE_DIR = ADIR + "/DATA_MISC/reshaped_files_validation_2018/"
ROOT_DIR_v2014_GAC = (BASE_DIR + "global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/2009/*cali*h5")
ROOT_DIR_v2018_GAC = (BASE_DIR + "global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/2009/*cali*h5")
ROOT_DIR_v2014_NPP = (BASE_DIR + "global_viirs_v2014_created20180914/Reshaped_Files_merged_caliop/npp/1km/2015/*/*h5")
ROOT_DIR_v2018_NPP = (BASE_DIR + "global_viirs_v2018_created20180907/Reshaped_Files_merged_caliop/npp/1km/2015/*/*h5")
# ROOT_DIR_v2018_NPP = (BASE_DIR + "global_viirs_v2018_created20181002_new_cmaprobv5/Reshaped_Files_merged_caliop/npp/1km/2015/*/*h5")
ROOT_DIR_v2014 = (BASE_DIR + "global_modis_v2014_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/*01_*cali*h5")
ROOT_DIR_v2018 = (BASE_DIR + "global_modis_v2018_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/*01_*cali*h5")
ROOT_DIR_PATMOSX = ADIR + "/VALIDATION_PATMOSX/Reshaped_Files/noaa18/5km/2009/*/*h5"
ROOT_DIR_CCI = (ADIR + "/DATA_MISC/reshaped_files_cci_noaa18_2009/V2/*2009*h5")

re_name = re.compile("_global_(\w+_\w+_\w+)\/")

files = glob(ROOT_DIR_v2014_GAC)
match_obj2014 = read_files(files)
pod14_d, pod14_n, f14_d, f14_n = make_pod_vector(match_obj2014)
match_obj2014 = None
files = glob(ROOT_DIR_v2018_GAC)
match_obj2018 = read_files(files)
pod18_d, pod18_n, f18_d, f18_n = make_pod_vector(match_obj2018)
files = glob(ROOT_DIR_PATMOSX)
match_objP = read_files(files)
podP_d, podP_n, fP_d, fP_n = make_pod_vector(match_objP)
files = glob(ROOT_DIR_CCI)
match_objC = read_files(files)
podC_d, podC_n, fC_d, fC_n = make_pod_vector(match_objC)
name = "GAC"


from matplotlib import pyplot as plt
fig = plt.figure(figsize=(9, 11))
ax = fig.add_subplot(211)
plt.plot(limits, pod18_d-pod18_n, '-r.', label="PPS-v2018")
# plt.plot(limits, pod18_d, 'c*')
plt.plot(limits, pod14_d-pod14_n, '-k.', label="PPS-v2014")
# plt.plot(limits, pod14_n, 'c*')
plt.plot(limits, podP_d-podP_n, '-m.', label="PATMOSX")
plt.plot(limits, podC_d-podC_n, '-y.', label="CCI-V2!")
plt.plot(limits, 0*np.array(limits), 'k:')
plt.ylabel("POD day - POD night")
plt.xlabel("Calipso total optical depth")
plt.legend()
ax = fig.add_subplot(212)
plt.plot(limits, pod18_d, '-r.', label="PPS-v2018")
plt.plot(limits, pod14_d, '-k.', label="PPS-v2014")
plt.plot(limits, podP_d, '-m.', label="PATMOSX")
plt.plot(limits, podC_n, '-y.', label="CCI-V2!")
plt.legend()
plt.ylabel("POD day")
plt.xlabel("Calipso total optical depth")
plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/PODdayMinusPODnight_latitude45_%s_and_od_sea.png"%(name), bbox_inches='tight')
plt.show()
