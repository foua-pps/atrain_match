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

#   Adam.Dybbroe

"""Read all matched data and make some plotting
"""

from glob import glob
import os.path
import numpy as np
from scipy import histogram
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 18})

#ROOT_DIR= "/home/a001865/DATA_MISC/reflection/some_files/5km_meteosat*2010*"
#my_title= "SEVIR_ONE_DAY"
ROOT_DIR = "/home/a001865/NN_CTTH/data/training_data_seviri_c4_ok_angles/5km_meteosat*2015*"
my_title = "SEVIRI_2015"
ROOT_DIR = "/home/a001865/NN_CTTH/data/training_data_seviri_c4_ok_angles/5km_meteosat*2010*"
my_title = "SEVIRI"
ROOT_DIR = "/home/a001865/NN_CTTH/data/sunglint_data_noaa18_gac/*"
my_title="AVHRR"

blinn=True
if blinn:
    my_title= my_title+"blinn"

def get_specular_refl_phong(azidiff, sunz, satz, exponent=25):
    a=1.0
    B =np.where(azidiff>=90, np.cos(np.radians(azidiff)),0)#??
    B = np.cos(np.radians(azidiff))
    cosalpha =  np.where(sunz>=90,0,np.cos(np.radians(sunz)))*np.cos(np.radians(satz)) -a*np.where(sunz>=90,0,np.sin(np.radians(sunz)))*np.sin(np.radians(satz))*B
    #refl = cosalpha**exponent
    #refl[refl==0] = -1
    return np.arccos(cosalpha)*180/np.pi

def get_specular_refl_blinn_phong(azidiff, sunz, satz, exponent=4*25):

    B =np.where(azidiff>=90, np.cos(np.radians(azidiff)),0)#??
    B = np.cos(np.radians(azidiff))
    dev =  np.where(sunz>=90,0,np.cos(np.radians(sunz)))*np.cos(np.radians(satz)) + np.where(sunz>=90,0,np.sin(np.radians(sunz)))*np.sin(np.radians(satz))*B
    cosalpha = (np.cos(np.radians(sunz))+np.cos(np.radians(satz)))/np.sqrt(2*(1+dev))
    #cosalpha = (np.cos(np.radians(sunz))+np.cos(np.radians(satz)))/(2*(1+dev))
    return np.arccos(cosalpha)*180/np.pi

def get_sunglint_info_pps2014(cloudtype_conditions):

    temp_val = (cloudtype_conditions>>3 & 1)
    sunglint_flag = temp_val == 1
    return  sunglint_flag


files = glob(ROOT_DIR)

from matchobject_io import (read_truth_imager_match_obj,
                            TruthImagerTrackObject)

match_calipso = TruthImagerTrackObject(truth='calipso')
for filename in files:
    print(os.path.basename(filename))
    i_match_calipso = read_truth_imager_match_obj(filename, skip_var=['cal_MODIS_cflag'])
    """

    i_sunz = i_match_calipso.imager.all_arrays['sunz']
    i_satz = i_match_calipso.imager.all_arrays['satz']
    i_azidiff = i_match_calipso.imager.all_arrays['azidiff']
    i_degree_from_reflection = get_specular_refl_phong(i_azidiff, i_sunz, i_satz)
    i_cloudtype_conditions = getattr(i_match_calipso.imager, 'cloudtype_conditions')
    i_ppssg=get_sunglint_info_pps2014(i_cloudtype_conditions)
    i_phongsg = i_degree_from_reflection<20
    print(np.sum(i_phongsg),np.sum(i_ppssg))
    if np.sum(i_phongsg)<np.sum(i_ppssg):
        print("!!!!!!!!!!!!!!!!!!!!!")
    """
    match_calipso = match_calipso + i_match_calipso


isCloudfree = np.logical_and(
    match_calipso.calipso.all_arrays['number_layers_found'] == 0,
    match_calipso.calipso.all_arrays['cloud_fraction'] ==0)
#isCloudfree = np.logical_and(
#    isCloudfree,
#    match_calipso.calipso.all_arrays['column_optical_depth_tropospheric_aerosols_532'] ==0)
sunz = match_calipso.imager.all_arrays['sunz']
satz = match_calipso.imager.all_arrays['satz']
azidiff = match_calipso.imager.all_arrays['azidiff']
tsurf = match_calipso.imager.all_arrays['surftemp']

print(sunz)

print(match_calipso.imager.all_arrays.keys())

print(sunz)
mu0 = np.cos(np.radians(sunz))
scaler = 24.35 / (2 * mu0 + np.sqrt(498.5225 * mu0 * mu0 + 1))

#r16 = match_calipso.imager.all_arrays['r16micron']*scaler
r06 = match_calipso.imager.all_arrays['r06micron']*scaler
r09 = match_calipso.imager.all_arrays['r09micron']*scaler
t11 = match_calipso.imager.all_arrays['bt11micron']
isCloudfree = np.logical_and(isCloudfree, t11>275)
isCloudfree = np.logical_and(isCloudfree, np.abs(t11-tsurf)<5)
nsidc_st = getattr(match_calipso.calipso, 'nsidc_surface_type')
igbp_st = getattr(match_calipso.calipso, 'igbp_surface_type')
cloudtype_conditions = getattr(match_calipso.imager, 'cloudtype_conditions')

#r09 = r06/r09


isCloudfree = np.logical_and(isCloudfree, np.equal(nsidc_st,0))
isCloudfree = np.logical_and(isCloudfree, np.equal(igbp_st,17))
isPPSSG = np.logical_and(
    isCloudfree,
    get_sunglint_info_pps2014(cloudtype_conditions))
isPPSSGn = np.logical_and(
    isCloudfree,
    ~get_sunglint_info_pps2014(cloudtype_conditions))




import matplotlib.pyplot as plt

degree_from_reflection = get_specular_refl_phong(azidiff, sunz, satz)
if blinn:
    degree_from_reflection = get_specular_refl_blinn_phong(azidiff, sunz, satz)


ppssg=get_sunglint_info_pps2014(cloudtype_conditions)
phongsg = degree_from_reflection<15
both = np.logical_and(ppssg,phongsg)
neither = np.logical_and(~ppssg,~phongsg)
onlypps = np.logical_and(ppssg,~phongsg)
onlyphong = np.logical_and(~ppssg,phongsg)
clear = match_calipso.calipso.all_arrays['cloud_fraction'] ==0
cloudy =  match_calipso.calipso.all_arrays['cloud_fraction'] >=0.98
#pps_clear =  match_calipso.imager.all_arrays['cloudmask'] ==0
#pps_cloudy = match_calipso.imager.all_arrays['cloudmask'] ==1
pps_clear = np.logical_or(np.equal(match_calipso.imager.cloudmask,3),
                          np.equal(match_calipso.imager.cloudmask,0))
pps_cloudy = np.logical_or(np.equal(match_calipso.imager.cloudmask,1),
                           np.equal(match_calipso.imager.cloudmask,2))
use = np.logical_or(pps_cloudy,pps_clear)
use = np.logical_and(use, np.equal(nsidc_st,0))
use = np.logical_and(use, np.equal(igbp_st,17))
use = np.logical_and(use, sunz<90)
clear[~use] = False
cloudy[~use]= False
print(my_title)
for selection, name in zip([both, neither, onlypps, onlyphong],
                           ["both", "neither", "onlypps", "onlyphong"]):
    n_ok_clear = np.sum(np.logical_and(clear, pps_clear)[selection])
    n_clear = np.sum(clear[selection])
    n_ok_cloudy = np.sum(np.logical_and(cloudy, pps_cloudy)[selection])
    n_cloudy = np.sum(cloudy[selection])
    print("{:s} & {:3.0f} & {:3.0f} & {:d}\\\\".format(
        name,
        100*n_ok_cloudy/n_cloudy,
        100*n_ok_clear/n_clear,
        np.sum(selection))
    )


xdat = degree_from_reflection[isCloudfree]
ydat = r09[isCloudfree]
xdat2 = degree_from_reflection[isPPSSG]
ydat2 = r09[isPPSSG]
data = []
#for ang in range(0,160,5):
#   data.append(list(ydat[np.logical_and(xdat>ang,xdat<ang+5)].ravel()))

#for ang
fig = plt.figure(figsize=(11,8))
ax = fig.add_subplot(111)
plt.fill_between([0,6.5],0,100, color='g', alpha=0.3)
plt.fill_between([0,5.5],0,100, color='g', alpha=0.3)
plt.fill_between([0,4.5],0,100, color='g', alpha=0.3)
data=[np.array(ydat[np.logical_and(xdat>=ang,xdat<ang+5)]) for ang in range(0,150,5)]
data2=[np.array(ydat2[np.logical_and(xdat2>=ang,xdat2<ang+5)]) for ang in range(0,150,5)]
box = plt.boxplot(data,
                  labels=[str(ang+2)+"-"+str(ang+5) for ang in range(0,150,5)],
                  #positions=[str(ang+2) for ang in range(0,150,5)],
                  whis='range', patch_artist=True)
for patch in box['boxes']:
    patch.set_facecolor('c')
box = plt.boxplot(data2,
                  labels=[str(ang)+"-"+str(ang+5) for ang in range(0,150,5)],
                  #positions=[str(ang+2) for ang in range(0,150,5)],
                  whis='range', patch_artist=True)
for patch in box['boxes']:
    patch.set_facecolor('r')
plt.title(my_title)
plt.xlabel("Angle between mirror reflection and satellite")
plt.ylabel("Reflectance R09 (%)")
ax1 =plt.gca()
ax1.set_ylim([0,100])
#xticks = ax1.xaxis.get_major_ticks()
plt.xticks(rotation=90)
#plt.xticks([950, 550, 150])
#        plt.yticks([950, 550, 150])
#for ind in range(8,150/5,2):
#    xticks[ind].set_visible(False)
plt.savefig("/home/a001865/DATA_MISC/reflection/figbox_"+my_title+".pdf", bbox_inches='tight', dpi=500, pad_inches = 0)
plt.savefig("/home/a001865/DATA_MISC/reflection/figbox_"+my_title+".png", bbox_inches='tight', dpi=500, pad_inches = 0)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
#plt.fill_between([0,6.5],0,100, color='g', alpha=0.3)
#plt.fill_between([0,5.5],0,100, color='g', alpha=0.3)
#plt.fill_between([0,4.5],0,100, color='g', alpha=0.3)
data=[len(ydat[np.logical_and(xdat>=ang,xdat<ang+5)]) for ang in range(0,150,5)]
data2=[len(ydat2[np.logical_and(xdat2>=ang,xdat2<ang+5)]) for ang in range(0,150,5)]
plt.hist(xdat, bins=list(range(0,150,5)), color = 'b', alpha=0.5)
plt.hist(xdat[xdat<30], bins=list(range(0,150,5)), color = 'c', alpha=1.0)
plt.hist(xdat2, bins=list(range(0,150,5)), color = 'r', alpha=0.5)
plt.hist(xdat2[xdat2<30], bins=list(range(0,150,5)), color = 'r', alpha=1.0)
plt.title(my_title)
plt.xlabel("Angle between mirror reflection and satellite")
plt.ylabel("Number of observations")
ax1 =plt.gca()
ax1.set_ylim([0,np.max(data[0:7])])
plt.savefig("/home/a001865/DATA_MISC/reflection/fighist_"+my_title+".pdf", bbox_inches='tight', dpi=500, pad_inches = 0)
plt.savefig("/home/a001865/DATA_MISC/reflection/fighist_"+my_title+".png", bbox_inches='tight', dpi=500, pad_inches = 0)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.fill_between([0,30],0,100, color='g', alpha=0.3)
plt.fill_between([0,25],0,100, color='g', alpha=0.3)
plt.fill_between([0,20],0,100, color='g', alpha=0.3)
plt.plot(degree_from_reflection[isCloudfree] , r09[isCloudfree], 'b*', alpha=1.0, rasterized=True)
ax1 = plt.plot(degree_from_reflection[isPPSSG] , r09[isPPSSG], '.r', alpha=1.0, rasterized=True)
plt.title(my_title)
plt.xlabel("Angle between mirror reflection and satellite")
plt.ylabel("Reflectance R09 (%)")
ax1 = plt.gca()
#ax1.set_xlim([0,60])
ax1.set_ylim([0,100])
plt.savefig("/home/a001865/DATA_MISC/reflection/fig_"+my_title+".pdf", bbox_inches='tight', dpi=500, pad_inches = 0)
plt.savefig("/home/a001865/DATA_MISC/reflection/fig_"+my_title+".png", bbox_inches='tight', dpi=500, pad_inches = 0)
fig = plt.figure()
ax = fig.add_subplot(221)
degree_from_reflection = get_specular_refl_phong(azidiff, sunz, satz)
if blinn:
    degree_from_reflection = get_specular_refl_blinn_phong(azidiff, sunz, satz)
plt.plot(degree_from_reflection[isCloudfree] , r09[isCloudfree], '.', alpha=1)
ax1 = plt.plot(degree_from_reflection[isPPSSG] , r09[isPPSSG], '.r', alpha=1)
ax1 = plt.gca()
#ax1.set_xlim([0,60])
ax1.set_ylim([0,100])
plt.xlabel("Angle between mirror reflection and satellite")
plt.ylabel("Reflectance R09 (%)")
ax = fig.add_subplot(222)
ax1 = plt.plot(degree_from_reflection[isPPSSG] , r09[isPPSSG], '.r', alpha=1)

plt.plot(degree_from_reflection[isPPSSGn] , r09[isPPSSGn], '.', alpha=1)
ax1 = plt.gca()
#ax1.set_xlim([0,60])
ax1.set_ylim([0,100])
plt.xlabel("Angle between mirror reflection and satellite")
plt.ylabel("Reflectance R09 (%)")
ax = fig.add_subplot(223)
degree_from_reflection = get_specular_refl_phong(azidiff, sunz, satz)
plt.plot(degree_from_reflection[isCloudfree] , r06[isCloudfree], '.', alpha=1)
ax1 = plt.plot(degree_from_reflection[isPPSSG] , r06[isPPSSG], '.r', alpha=1)
ax1 = plt.gca()
#ax1.set_xlim([0,60])
ax1.set_ylim([0,100])
plt.xlabel("Angle between mirror reflection and satellite")
plt.ylabel("Reflectance R06 (%)")
ax = fig.add_subplot(224)
ax1 = plt.plot(degree_from_reflection[isPPSSG] , r06[isPPSSG], '.r', alpha=1)

plt.plot(degree_from_reflection[isPPSSGn] , r06[isPPSSGn], '.', alpha=1)
ax1 = plt.gca()
#ax1.set_xlim([0,60])
ax1.set_ylim([0,100])
plt.xlabel("Angle between mirror reflection and satellite")
plt.ylabel("Reflectance R06 (%)")
#plt.show()





#plot for polar!
