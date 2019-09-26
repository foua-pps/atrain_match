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
"""Read all matched data and make some plotting
"""
from plot_ctth_bias_distributions import PlotAndDataObject, extract_data
from my_dir import ADIR
import os
import re
from glob import glob
import numpy as np
from matchobject_io import (readCaliopImagerMatchObj,
                            DataObject,
                            CloudsatImagerTrackObject,
                            readCloudsatImagerMatchObj,
                            CalipsoImagerTrackObject)
from plot_kuipers_on_area_util import (PerformancePlottingObject,
                                       ppsMatch_Imager_CalipsoObject)
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 18})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

tag_dict = {"old": "(a) PPS-v2014",
            "mlvl2": "(b) MODIS-C6",
            "nnant": "(c) NN-AVHRR",
            "pps": "(c) NN-AVHRR",
            "nna1nt": "(d) NN-AVHRR1",
            "nnvnt": "(e) NN-VIIRS",
            "nnm2nt": "(f) NN-MERSI-2",
            "nnmint": "(h) NN-MetImage",
            "nnmintnco2": "(g) NN-MetImage-NoCO$_2$",

            }


def do_one_subplot(plt_obj, ax, fig, compare, truth='height_c', vmax=250, do_colorbar=False, ptype='scatter', height_calipso=False):
    x = getattr(plt_obj, truth)
    y = getattr(plt_obj, compare)
    use = np.logical_and(x >= 0, getattr(plt_obj, 'use_all'))
    x = x[use]
    y = y[use]
    xmin = 0
    # xmax = 22000 # m~22km
    # if height_calipso:
    xmax = 22000  # ~24km
    binsize = 250  # binsixe in meters
    #binsize = xmax//4
    # vmax=vmax*100
    # print binsize
    to_km = 0.001
    if 'pressure' in truth:
        xmax = 1100
        binsize = 10  # in hPa
        to_km = 1.0  # Keep it in hPa
    ymin = xmin
    ymax = xmax
    the_title = "xxxxx"
    for tag in sorted(tag_dict.keys()):
        print tag
        if tag in compare:
            the_title = tag_dict[tag]

    from scipy.stats import gaussian_kde
    n_edges = int(xmax*1.0/binsize) + 1
    # n_edges=3 for testing
    # print "bins", n_edges
    edgesx = np.linspace(xmin, xmax, n_edges)
    edgesy = np.linspace(ymin, ymax, n_edges)
    # print edgesx
    H, xe, ye = np.histogram2d(x, y, bins=[edgesx, edgesy])
    # print H
    xi = np.searchsorted(edgesx, x)  # - edgesx[0])/(n_edges+1)).astype(np.int)
    yi = np.searchsorted(edgesy, y)  # np.floor((y - edgesy[0])/(n_edges+1)).astype(np.int)  #-1?
    # print max(y), max(x)
    # print max(xi), max(yi)
    # print x[1:5], xi[1:5]
    xi = xi - 1
    yi = yi - 1
    # only needed if max to small !
    # n_edges =3,
    # edges: [0, 1, 2]
    # np.searchsorted [(0), 1, 2, (3)]
    # xi-t [(-1), 0, 1, (2)]
    # dim(H)=2x2
    # H(0) and H(1) is ok. H(2) is not!
    if 'pressure' not in truth:
        print compare
        print " n to large x values {:d}".format(np.sum(xi >= (n_edges-1)))
        print " n to large y values {:d}".format(np.sum(yi >= (n_edges-1)))
        # yi(yi==3)==2 no one
        yi[yi == n_edges] = n_edges - 1
        xi[xi == n_edges] = n_edges - 1
        # yi(yi==2)==1 This is what we need!
        yi[yi == n_edges - 1] = n_edges - 2
        xi[xi == n_edges - 1] = n_edges - 2
    z = H[xi, yi]
    # z=H[yi, xi]

    import copy
    my_cmap = copy.copy(matplotlib.cm.viridis)
    my_cmap.set_under(color='1.0', alpha=1)

    perc90 = np.percentile(z, 75)
    print "perc90", perc90
    # z[z>perc90]=perc90
    # nicer but too slow
    # points = np.vstack([x, y])
    # kde= gaussian_kde(points)
    # z=kde(points)
    idx = z.argsort()
    # xc = np.array(range(xmin, xmax, 10))
    # y1 = xc
    # y2 = -xc
    # ax.fill_betweenx(xc, y2, y1, facecolor='green', alpha=0.3)
    # plt.scatter(x[idx], y[idx], c=z[idx], edgecolor='', cmap='OrRd', label=label)
    # my_cmap=copy.copy(matplotlib.cm.BrBG)

    my_cmap = copy.copy(matplotlib.cm.get_cmap("inferno_r", lut=100))
    cmap_vals = my_cmap(np.arange(100))  # extractvalues as an array
    print cmap_vals[0]
    cmap_vals[0:5] = cmap_vals[5]  # change the first values
    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "new_inferno_r", cmap_vals)

    z_values_to_use = z[idx]
    few = z_values_to_use <= 10
    more = z_values_to_use > 10

    if ptype in ['hexbin']:
        b = plt.hexbin(to_km*x, to_km*y, gridsize=binsize, cmap=my_cmap,
                       vmin=10, vmax=vmax)
    if ptype in ['hist2d']:
        b = plt.hist2d(to_km*x, to_km*y, bins=binsize, cmap=my_cmap,
                       vmin=1, vmax=vmax)
    if ptype in ['scatter']:

        # plt.scatter(to_km*x[idx][few], to_km*y[idx][few], c=z[idx][few],
        #           edgecolor='', cmap='inferno_r', vmin=1, vmax=vmax,
        #           alpha=0.2, marker='.', edgecolors=None)
        # b = plt.scatter(to_km*x[idx][more], to_km*y[idx][more], c=z[idx][more],
        #               edgecolor='', cmap='inferno_r', vmin=1, vmax=vmax,
        #               alpha=1.0, marker='.', edgecolors=None)
        b = plt.scatter(to_km*x[idx], to_km*y[idx], c=z[idx],
                        edgecolor='', cmap=my_cmap, vmin=1, vmax=vmax,
                        alpha=1.0, marker='.', edgecolors=None, rasterized=True)
        # b = plt.plot(to_km*x[idx], to_km*y[idx], c=z[idx],
        #             cmap=my_cmap, vmin=1, vmax=vmax,
        #            alpha=1.0, marker='.', rasterized=True)

        if 'pressure' in truth:
            ax.text(to_km*xmax*0.95, to_km*xmax*0.05, the_title, bbox={'facecolor': 'blue', 'alpha': 0.0, 'pad': 1})
            # plt.plot([xmin, 0.9*to_km*xmax], [xmin+25, 0.9*to_km*xmax+25], 'w', alpha=0.5)
            # plt.plot([xmin, 0.9*to_km*xmax], [xmin-25, 0.9*to_km*xmax-25], 'w', alpha=0.5)
            plt.plot([0.1*to_km*xmax, to_km*xmax], [0.1*to_km*xmax, to_km*xmax], 'b:')
        else:
            ax.text(xmax*to_km*0.05, 0.90*to_km*xmax, the_title, bbox={'facecolor': 'blue', 'alpha': 0.0, 'pad': 1})
            # plt.plot([xmin, 0.9*to_km*xmax], [xmin+0.5, 0.9*to_km*xmax+0.5], color='grey', alpha=1.0)
            # plt.plot([xmin, 0.9*to_km*xmax], [xmin-0.5, 0.9*to_km*xmax-0.5], color='grey', alpha=10)
            plt.plot([xmin, 0.85*to_km*xmax], [xmin, 0.85*to_km*xmax], 'b:')
    # plt.scatter(x, y, alpha=0.1, label=label)

    ax.set_ylim(ymin, to_km*ymax)
    ax.set_xlim(xmin, to_km*xmax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if 'pressure' in truth:
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        plt.xticks([950, 550, 150])
        plt.yticks([950, 550, 150])
    else:
        plt.xticks([5, 10, 15])
        plt.yticks([0, 5, 10, 15])

    if do_colorbar and ptype in ['scatter']:
        cax = fig.add_axes([0.80, 0.65, 0.06, 0.22])
        cbar = fig.colorbar(b, cax=cax, )


def do_the_scatter_plot(plt_obj_cali_new, plt_obj_csat_new, month):
    vmax = 300
    if len(month) > 3:
        print("Much data")
        vmax = 1800*2

    fig = plt.figure(figsize=(11, 11))

    ax = fig.add_subplot(331, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_old', vmax=vmax, height_calipso=True)
    ax = fig.add_subplot(332, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_mlvl2', vmax=vmax, do_colorbar=True, height_calipso=True)
    ax = fig.add_subplot(334, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_pps', vmax=vmax, height_calipso=True)
    ax = fig.add_subplot(335, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_nna1nt', vmax=vmax, height_calipso=True)
    ax = fig.add_subplot(336, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_nnvnt', vmax=vmax, height_calipso=True)
    ax = fig.add_subplot(337, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_nnm2nt', vmax=vmax, height_calipso=True)
    ax = fig.add_subplot(338, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_nnmintnco2', vmax=vmax, height_calipso=True)
    ax = fig.add_subplot(339, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'height_nnmint', vmax=vmax, height_calipso=True)
    fig.text(0.5, 0.04, 'Cloud height CALIOP (km)', ha='center', fontsize=18)
    fig.text(0.04, 0.5, 'Retrieved height (km)', va='center', rotation='vertical', fontsize=18)
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig_calipso_scatter_inferno_r_modified_max_height_%s.png" %
                (month), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig_calipso_scatter_inferno_r_modified_max_height_%s.pdf" %
                (month), bbox_inches='tight')
    plt.close("all")

    fig = plt.figure(figsize=(11, 11))

    ax = fig.add_subplot(331, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_old', vmax=vmax, truth='pressure_c')
    ax = fig.add_subplot(332, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_mlvl2', vmax=vmax, truth='pressure_c', do_colorbar=True)
    ax = fig.add_subplot(334, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_nnant', vmax=vmax, truth='pressure_c')
    ax = fig.add_subplot(335, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_nna1nt', vmax=vmax, truth='pressure_c')
    ax = fig.add_subplot(336, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_nnvnt', vmax=vmax, truth='pressure_c')
    ax = fig.add_subplot(337, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_nnm2nt', vmax=vmax, truth='pressure_c')
    ax = fig.add_subplot(338, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_nnmintnco2', vmax=vmax, truth='pressure_c')
    ax = fig.add_subplot(339, aspect='equal')
    do_one_subplot(plt_obj_cali_new, ax, fig, 'pressure_nnmint', vmax=vmax, truth='pressure_c')
    fig.text(0.5, 0.04, 'Cloud top pressure CALIOP (hPa)', ha='center', fontsize=18)
    fig.text(0.04, 0.5, 'Retrieved cloud top pressure (hPa)', va='center', rotation='vertical', fontsize=18)
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig04_scatter_inferno_r_modified_max_pressure_%s.png" %
                (month), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig04_%s.pdf" % (month), bbox_inches='tight')
    vmax = 200
    if len(month) > 3:
        print("Much data")
        vmax = 1250*2

    fig = plt.figure(figsize=(11, 11))

    ax = fig.add_subplot(331, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_old', vmax=vmax)
    ax = fig.add_subplot(332, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_mlvl2', vmax=vmax, do_colorbar=True)
    ax = fig.add_subplot(334, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_pps', vmax=vmax)
    ax = fig.add_subplot(335, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_nna1nt', vmax=vmax)
    ax = fig.add_subplot(336, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_nnvnt', vmax=vmax)
    ax = fig.add_subplot(337, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_nnm2nt', vmax=vmax)
    ax = fig.add_subplot(338, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_nnmintnco2', vmax=vmax)
    ax = fig.add_subplot(339, aspect='equal')
    do_one_subplot(plt_obj_csat_new, ax, fig, 'height_nnmint', vmax=vmax)
    fig.text(0.5, 0.04, 'Cloud top height CPR (CloudSat) (km)', ha='center', fontsize=18)
    fig.text(0.04, 0.5, 'Retrieved cloud top height (km)', va='center', rotation='vertical', fontsize=18)
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig03_cloudsat_scatter_inferno_r_modified_max_height_%s.png" %
                (month), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/CTTH_PLOTS/fig03_%s.pdf" % (month), bbox_inches='tight')


def get_plot_object_nn_ctth_modis_lvl2_cloudsat(month):
    day_str = "01st"
    ROOT_DIR = (
        ADIR + "/DATA_MISC/reshaped_files/"
        "global_modis_%s_created20180316/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
    # "global_modis_%s_created20170330/Reshaped_Files_merged_cloudsat/eos2/1km/2010/%s/*h5")
    # match_clsat = CloudsatImagerTrackObject()
    plt_obj = PlotAndDataObject()
    print ROOT_DIR % (day_str, month)
    files = glob(ROOT_DIR % (day_str, month))
    for filename in files:
        print filename
        match_clsat_new = readCloudsatImagerMatchObj(filename)
        plt_obj += extract_data(match_clsat_new, sat='cloudsat')
    return plt_obj


def get_plot_object_nn_ctth_modis_lvl2(month):
    day_str = "01st"
    ROOT_DIR = (
        ADIR + "/DATA_MISC/reshaped_files/"
        # "global_modis_%s_created20170504/Reshaped_Files_merged/eos2/1km/2010/%s/*h5")
        "global_modis_%s_created20180316/Reshaped_Files_merged_calipso_cbase/eos2/1km/2010/%s/*h5")
    plt_obj = PlotAndDataObject()
    print ROOT_DIR % (day_str, month)
    files = glob(ROOT_DIR % (day_str, month))
    for filename in files:
        print filename
        match_calipso_new = readCaliopImagerMatchObj(filename)
        plt_obj += extract_data(match_calipso_new, sat='calipso')
    return plt_obj


def do_the_plots():

    plt_obj_cali = PlotAndDataObject()
    plt_obj_csat = PlotAndDataObject()

    merged_months = ""
    for month in ["02", "04", "06", "08", "10", "12"]:
        # for month in ["08"]:
        merged_months += month
        plt_obj_cali_new = get_plot_object_nn_ctth_modis_lvl2(month)
        plt_obj_csat_new = get_plot_object_nn_ctth_modis_lvl2_cloudsat(month)
        do_the_scatter_plot(plt_obj_cali_new, plt_obj_csat_new, month)
        plt_obj_cali += plt_obj_cali_new
        plt_obj_csat += plt_obj_csat_new

    do_the_scatter_plot(plt_obj_cali, plt_obj_csat, merged_months)


if __name__ == "__main__":
    do_the_plots()
